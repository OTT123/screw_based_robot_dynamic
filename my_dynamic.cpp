#include "./my_dynamic.h"


//SE3 求逆矩阵
Eigen::MatrixXd inv_SE3(Eigen::MatrixXd T)
{
    Eigen::MatrixXd R = Eigen::MatrixXd(3, 3);
    Eigen::MatrixXd p = Eigen::MatrixXd(3, 1);
    Eigen::MatrixXd inv_T = Eigen::MatrixXd(4, 4);
    R = T.block(0, 0, 3, 3);
    p = T.block(0, 3, 3, 1);

    //std::cout<<-R.transpose() * p<<std::endl;
    inv_T << R.transpose(), -R.transpose() * p,
        0, 0, 0, 1;
    return inv_T;
}
//@output :[p]
//3D向量的反对称矩阵
Eigen::MatrixXd vec2so3(Eigen::MatrixXd p)
{
    Eigen::MatrixXd R = Eigen::MatrixXd(3, 3);
    R(0, 0) = 0;
    R(0, 1) = -p(2);
    R(0, 2) = p(1);
    R(1, 0) = p(2);
    R(1, 1) = 0;
    R(1, 2) = -p(0);
    R(2, 0) = -p(1);
    R(2, 1) = p(0);
    R(2, 2) = 0;
    return R;
}


//SE3的伴随矩阵
//@param T:4x4 SE3
//@output :ADT 6x6
Eigen::MatrixXd AD_SE3(Eigen::MatrixXd T)
{
    Eigen::MatrixXd ADT = Eigen::MatrixXd(6, 6);
    Eigen::MatrixXd R = Eigen::MatrixXd(3, 3);
    Eigen::MatrixXd p = Eigen::MatrixXd(3, 1);
    Eigen::MatrixXd zero = Eigen::MatrixXd(3, 3);
    R = T.block(0, 0, 3, 3);
    p = T.block(0, 3, 3, 1);
    zero.setZero();
    ADT << R, zero,
        vec2so3(p)* R, R;
    return ADT;
}
//se3的伴随矩阵
//@param S:6x1 screw axis
//@output :[s] or ads,6x6
Eigen::MatrixXd ad_screw(Eigen::MatrixXd S)
{
    Eigen::MatrixXd w = Eigen::MatrixXd(3, 1);
    Eigen::MatrixXd v = Eigen::MatrixXd(3, 1);
    Eigen::MatrixXd ads = Eigen::MatrixXd(6, 6);
    Eigen::MatrixXd zero = Eigen::MatrixXd(3, 3); zero.setZero();
    w = S.block(0, 0, 3, 1);
    v = S.block(3, 0, 3, 1);
    ads << vec2so3(w), zero,
        vec2so3(v), vec2so3(w);

    return ads;
}

//旋转矩阵指数运算
//输入螺旋轴和旋转角度，输出对应的SE3
//@param w:3x1旋转轴
//@param theta:旋转角度
//@output :旋转矩阵R SO3
Eigen::MatrixXd exp_so3_2_SO3(Eigen::MatrixXd w, double theta)
{
    Eigen::MatrixXd I = Eigen::MatrixXd(3, 3);
    Eigen::MatrixXd R = Eigen::MatrixXd(3, 3);
    I << 1, 0, 0,
        0, 1, 0,
        0, 0, 1;
    R = I + std::sin(theta) * vec2so3(w) + (1 - std::cos(theta)) * vec2so3(w) * vec2so3(w);
    return R;
}


//螺旋轴指数运算
//输入螺旋轴和旋转角度，输出对应的SE3
//@param s:6x1螺旋轴[w v]
//@param theta:旋转角度
//@output :齐次变换矩阵 SE3
Eigen::MatrixXd exp_twist_2_SE3(Eigen::MatrixXd S, double theta)
{
    Eigen::MatrixXd R = Eigen::MatrixXd(3, 3);
    Eigen::MatrixXd G = Eigen::MatrixXd(3, 3);
    Eigen::MatrixXd I = Eigen::MatrixXd(3, 3);
    Eigen::MatrixXd T = Eigen::MatrixXd(4, 4);
    Eigen::MatrixXd w = Eigen::MatrixXd(3, 1);
    Eigen::MatrixXd v = Eigen::MatrixXd(3, 1);
    Eigen::MatrixXd zero = Eigen::MatrixXd(1, 3); zero.setZero();
    w = S.block(0, 0, 3, 1);
    v = S.block(3, 0, 3, 1);
    R = exp_so3_2_SO3(w, theta);
    I << 1, 0, 0,
        0, 1, 0,
        0, 0, 1;
    G = I * theta + (1 - std::cos(theta)) * vec2so3(w) + (theta - std::sin(theta)) * vec2so3(w) * vec2so3(w);
    T << R, G* v,
        zero, 1;
    return T;
}



//惯性矩阵转换为空间惯量矩阵
//@param inertial:3x1
//@output inerG:6x6
Eigen::MatrixXd Inertial2inerG(Eigen::MatrixXd inertial_i, double m)
{
    Eigen::MatrixXd inerG = Eigen::MatrixXd(6, 6);
    Eigen::MatrixXd zero = Eigen::MatrixXd(3, 3);
    Eigen::MatrixXd inertial = Eigen::MatrixXd(3, 3);
    Eigen::MatrixXd I = Eigen::MatrixXd(3, 3);
    I << 1, 0, 0,
        0, 1, 0,
        0, 0, 1;
    inertial << inertial_i(0), 0, 0,
        0, inertial_i(1), 0,
        0, 0, inertial_i(2);
    zero.setZero();

    inerG << inertial, zero,
        zero, m* I;
    return inerG;
}

/// <summary>
/// 质量分布在不同坐标系的转换
/// </summary>
/// <param name="T">
/// B坐标系 相对与 A坐标系 的位姿
/// </param>
/// <param name="iner">
/// 连杆在B坐标系的质量分布 3x3
/// </param>
/// <param name="m">
/// 连杆质量
/// </param>
/// <returns>
/// 连杆在A坐标系的质量分布
/// </returns>
Eigen::Matrix<double, 3, 3> iner_trans(Eigen::Matrix<double, 4, 4> T, Eigen::Matrix<double, 3, 3> iner, double m)
{
    Eigen::Matrix<double, 3, 3> R;
    Eigen::Matrix<double, 3, 1> p;
    Eigen::Matrix<double, 3, 3> one3x3;
    Eigen::Matrix<double, 3, 3> trans_iner;
    R = T.block(0, 0, 3, 3);;
    p = T.block(0, 2, 3, 1);
    one3x3.setIdentity();
    trans_iner = R * iner * R.transpose() + m * (p.transpose() * p * one3x3 - p * p.transpose());
    return trans_iner;
}


Eigen::Matrix<double, 6, 1> Jaccobi_s_i(Eigen::Matrix<double, 6, 1> theta, Eigen::Matrix<double, 6, 1> screw[], int col_i)
{
    Eigen::Matrix<double, 6, 1> Jaccobi_s_col_i;
    Eigen::Matrix<double, 4, 4> temp1;
    temp1 << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1;
    for (int i = 0; i < col_i - 1; ++i)
    {
        temp1 = temp1 * exp_twist_2_SE3(screw[i], theta(i));
    }
    Jaccobi_s_col_i = AD_SE3(temp1) * screw[col_i - 1];
    return Jaccobi_s_col_i;
}

Eigen::Matrix<double, 6, 6> Jaccobi_s(Eigen::Matrix<double, 6, 1> theta, Eigen::Matrix<double, 6, 1> screw[], int n)
{
    Eigen::Matrix<double, 6, 6> Jaccobi_s;
    //计算第一列
    for (int i = 0; i < n; ++i)
    {
        Jaccobi_s.col(i) = Jaccobi_s_i(theta, screw, i + 1);
    }
    return Jaccobi_s;
}

//RPY2rot
Eigen::Matrix<double, 3, 3> rpyToRotationMatrix(const Eigen::Vector3d& rpy)
{
    double roll = rpy(0);
    double pitch = rpy(1);
    double yaw = rpy(2);

    Eigen::Matrix3d rotX, rotY, rotZ;
    rotX << 1, 0, 0,
        0, std::cos(roll), -std::sin(roll),
        0, std::sin(roll), std::cos(roll);

    rotY << std::cos(pitch), 0, std::sin(pitch),
        0, 1, 0,
        -std::sin(pitch), 0, std::cos(pitch);

    rotZ << std::cos(yaw), -std::sin(yaw), 0,
        std::sin(yaw), std::cos(yaw), 0,
        0, 0, 1;

    Eigen::Matrix<double, 3, 3> rotationMatrix = rotZ * rotY * rotX;
    return rotationMatrix;
}


Eigen::VectorXd parse_str_vec(std::string scr, int n)
{
    Eigen::VectorXd des(n);
    std::istringstream iss(scr);
    for (int i = 0; i < n; ++i)
    {
        double value;
        if (iss >> value)
        {
            des(i) = value;
        }
        else {
            std::cerr << "Error parsing centroid frame transepose " << i << std::endl;
            exit(0);
        }
    }
    return des;
}

// [Ixx Iyy Izz Iyz Ixz Ixy] to 3x3 matrix
Eigen::Matrix<double, 3, 3> roatate_iner_line_2_matrix(Eigen::Matrix<double, 6, 1> line)
{
    Eigen::Matrix<double, 3, 3> iner;
    iner(0, 0) = line[0];
    iner(1, 1) = line[1];
    iner(2, 2) = line[2];
    iner(0, 1) = line[5];
    iner(0, 2) = line[4];
    iner(1, 2) = line[3];
    iner(1, 0) = line[0, 1];
    iner(2, 0) = line[0, 2];
    iner(2, 1) = line[1, 2];
    return iner;
}


my_Robot_dynamic::my_Robot_dynamic()
{
    this->stack_A.setZero();
    this->stack_inerG.setZero();
    this->stack_ad_V.setZero();
    this->stack_ad_A_theta_dot.setZero();
    this->stack_W.setZero();
    this->stack_V_base.setZero();
    this->stack_V_base_dot.setZero();
    this->stack_L.setZero();
    this->stack_F_tip.setZero();
}


/// <summary>
/// 和github相比我的corMat计算错误
/// 困扰了一个月的bug
/// 其实就是下面这段代码出错了
/// iter_twistV_dot[i + 1] = AD_SE3(iter_T[i]) * iter_twistV_dot[i] + ad_screw(iter_twistV[i]) * this->links_A[i] * iter_q_dot[i] + this->links_A[i] * iter_q_ddot[i];  
/// 注意看ad_screw(iter_twistV[i])
/// 其中标号应该是iter_twistV[i+1]
/// 回看书上算法，，这一步使用的iter_twistV就是该次迭代中的上一步算出来的，但是我这样写就变成了上一次迭代的上一步计算结果
/// 以后要多加注意角标问题

Eigen::Matrix<double, DOF, 1> my_Robot_dynamic::forword_dyn(
    Eigen::Matrix<double, DOF, 1> iter_q,
    Eigen::Matrix<double, DOF, 1> iter_q_dot,
    Eigen::Matrix<double, DOF, 1> iter_q_ddot,
    Eigen::Matrix<double, 3, 1> gravity)
{
    //迭代所需变量
    Eigen::Matrix<double, 4, 4> iter_T[DOF + 1];
    Eigen::Matrix<double, 6, 1> iter_twistV[DOF + 1];
    Eigen::Matrix<double, 6, 1> iter_twistV_dot[DOF + 1];
    Eigen::Matrix<double, 6, 1> iter_wrenchF[DOF + 1];
    Eigen::Matrix<double, DOF, 1> iter_tor;

    iter_twistV[0] << 0, 0, 0, 0, 0, 0;
    iter_twistV_dot[0] << 0, 0, 0, -gravity;
    iter_wrenchF[DOF] << 0, 0, 0, 0, 0, 0;    //末端执行器的力

    for (int i = 0; i < DOF; ++i)
    {
        iter_T[i] = exp_twist_2_SE3(-this->links_A[i], iter_q[i]) * inv_SE3(this->links_M[i]);
        iter_twistV[i + 1] = AD_SE3(iter_T[i]) * iter_twistV[i] + this->links_A[i] * iter_q_dot[i];
        iter_twistV_dot[i + 1] = AD_SE3(iter_T[i]) * iter_twistV_dot[i] + ad_screw(iter_twistV[i + 1]) * this->links_A[i] * iter_q_dot[i] + this->links_A[i] * iter_q_ddot[i];
        //std::cout << "AdTi" << i << ": " << std::endl;
        //std::cout << AD_SE3(iter_T[i]) << std::endl;
    
    }
    iter_T[DOF].setIdentity();
    for (int i = DOF - 1; i >= 0; --i)
    {
        iter_wrenchF[i] = AD_SE3(iter_T[i + 1]).transpose() * iter_wrenchF[i + 1] + this->inertial_G[i] * iter_twistV_dot[i + 1] - ad_screw(iter_twistV[i + 1]).transpose() * this->inertial_G[i] * iter_twistV[i + 1];
        iter_tor[i] = iter_wrenchF[i].transpose() * links_A[i];
    }
    return iter_tor;
}



Eigen::Matrix<double, DOF, DOF> my_Robot_dynamic::cal_massMat(
    Eigen::Matrix<double, DOF, 1> iter_q)
{
    Eigen::Matrix<double, 3, 1> gravity;
    gravity.setZero();
    Eigen::Matrix<double, DOF, 1> iter_q_dot; iter_q_dot.setZero();
    Eigen::Matrix<double, DOF, 1> iter_q_ddot; iter_q_ddot.setZero();
    Eigen::Matrix<double, DOF, DOF> massMat;
    for (int i = 0; i < DOF; ++i)
    {
        iter_q_ddot.setZero();
        iter_q_ddot[i] = 1;
        massMat.col(i) = forword_dyn(iter_q, iter_q_dot, iter_q_ddot, gravity);
    }
    return massMat;
}


Eigen::Matrix<double, DOF, 1>  my_Robot_dynamic::cal_corMat(
    Eigen::Matrix<double, DOF, 1> iter_q,
    Eigen::Matrix<double, DOF, 1> iter_q_dot)
{
    Eigen::Matrix<double, 3, 1> gravity;
    gravity.setZero();
    Eigen::Matrix<double, DOF, 1> corMat;
    Eigen::Matrix<double, DOF, 1> iter_q_ddot; iter_q_ddot.setZero();
    corMat = forword_dyn(iter_q, iter_q_dot, iter_q_ddot, gravity);
    return corMat;
}

Eigen::Matrix<double, DOF, 1>  my_Robot_dynamic::cal_graMat(
    Eigen::Matrix<double, DOF, 1> iter_q,
    Eigen::Matrix<double, 3, 1> gravity)
{
    Eigen::Matrix<double, DOF, 1> graMat;
    Eigen::Matrix<double, DOF, 1> iter_q_dot; iter_q_dot.setZero();
    Eigen::Matrix<double, DOF, 1> iter_q_ddot; iter_q_ddot.setZero();
    graMat = forword_dyn(iter_q, iter_q_dot, iter_q_ddot, gravity);
    return graMat;
}

//机器人参数
my_Robot::my_Robot()
{
    ;
};

void my_Robot::set_param_from_src_code()
{
    this->links_T[0] << 1, 0, 0, 2.4327E-05,   //Tmesh{1,0}
        0, 1, 0, -2.2695E-05,
        0, 0, 1, -0.105,
        0, 0, 0, 1;
    this->links_T[1] << 1, 0, 0, -0.080113,
        0, 1, 0, 2.2695E-05,
        0, 0, 1, 0,
        0, 0, 0, 1;
    this->links_T[2] << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, -0.17447,
        0, 0, 0, 1;
    this->links_T[3] << 1, 0, 0, 0.084516,
        0, 1, 0, 0,
        0, 0, 1, -0.17442,
        0, 0, 0, 1;
    this->links_T[4] << 1, 0, 0, -0.080113,
        0, 1, 0, -2.2695E-05,
        0, 0, 1, 0,
        0, 0, 0, 1;
    this->links_T[5] << 1, 0, 0, 0,
        0, 1, 0, 4.5391E-05,
        0, 0, 1, -0.080089,
        0, 0, 0, 1;
    this->links_T[6] << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1;
    this->T_C[0] << 1, 0, 0, 4.33680868994202E-18,
        0, 1, 0, -6.00214600243731E-11,
        0, 0, 1, 0.0271404483614374,
        0, 0, 0, 1;
    this->T_C[1] << 1, 0, 0, 0.00464247258603146,
        0, 1, 0, -3.55800695980568E-05,
        0, 0, 1, -0.000556506974114818,
        0, 0, 0, 1;
    this->T_C[2] << 1, 0, 0, -0.00354399449686532,
        0, 1, 0, 6.83167072432211E-06,
        0, 0, 1, 0.0535596075762168,
        0, 0, 0, 1;
    this->T_C[3] << 1, 0, 0, -0.046552753550141,
        0, 1, 0, -6.68259662081616E-06,
        0, 0, 1, 0.0523242797053055,
        0, 0, 0, 1;
    this->T_C[4] << 1, 0, 0, 0.0182269184920331,
        0, 1, 0, 1.24412678982888E-05,
        0, 0, 1, 0.00395855146065832,
        0, 0, 0, 1;
    this->T_C[5] << 1, 0, 0, 0.00395855202736296,
        0, 1, 0, -3.51379306736399E-05,
        0, 0, 1, 0.0182025903111476,
        0, 0, 0, 1;
    this->T_C[6] << 1, 0, 0, -0.00231632591813243,
        0, 1, 0, 4.03819823737578E-05,
        0, 0, 1, -5.08515816782795E-05,
        0, 0, 0, 1;
    for (int i = 0; i < DOF; ++i)
    {
        //links_M[i] = T_C[i] * links_T[i] * inv_SE3(T_C[i + 1]);
        this->links_M[i] = this->links_T[i];    //近似处理
        //std::cout << links_M[i] << std::endl;
        //std::cout << "********************************" << std::endl;
    }
    this->links_M[DOF] << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1;
    //旋量轴在mesh中的定义
    this->links_B[0] << 0, 0, -1, 0, 0, 0;
    this->links_B[1] << -1, 0, 0, 0, 0, 0;
    this->links_B[2] << 1, 0, 0, 0, 0, 0;
    this->links_B[3] << -1, 0, 0, 0, 0, 0;
    this->links_B[4] << 0, 0, -1, 0, 0, 0;
    this->links_B[5] << -1, 0, 0, 0, 0, 0;
    //std::cout << "A" << std::endl;
    for (int i = 0; i < DOF; ++i)
    {
        //links_A[i] = AD_SE3(inv_SE3(T_C[i])) * links_B[i];
        this->links_A[i] = this->links_B[i]; //近似处理
        //std::cout << links_A[i] << std::endl;
        //std::cout << "*" << std::endl;
    }

    this->links_mass << 0.166259789544396, 0.316094041078938, 0.323135746651779, 0.174947986817189, 0.174947985033356, 0.122091979083374;
    Eigen::Matrix<double, 3, 3> inertial_i[DOF];
    inertial_i[0] << 6.76518298922659E-05, -1.5556489245956E-08, -8.15879264327807E-06,
        -1.5556489245956E-08, 6.70385565720365E-05, 4.11369338798553E-08,
        -8.15879264327807E-06, 4.11369338798553E-08, 6.03698353342124E-05;
    inertial_i[1] << 0.000174638310500152, -4.1835416811351E-08, -2.08841744154906E-07,
        -4.1835416811351E-08, 0.000184936237600061, 1.60668786323061E-08,
        -2.08841744154906E-07, 1.60668786323061E-08, 0.000107780595065675;
    inertial_i[2] << 0.000176611314431077, 4.18325643591624E-08, 2.08841259569185E-07,
        4.18325643591624E-08, 0.000185932099122903, 1.40317731158265E-08,
        2.08841259569185E-07, 1.40317731158265E-08, 0.000108776454638344;
    inertial_i[3] << 6.61039276820243E-05, 4.2340099438174E-08, -8.15879264340351E-06,
        4.2340099438174E-08, 6.76311730651952E-05, -1.4353323683088E-08,
        -8.15879264340351E-06, -1.4353323683088E-08, 6.34296690009448E-05;
    inertial_i[4] << 6.34296685392068E-05, 1.43441999593359E-08, -8.15879467982627E-06,
        1.43441999593359E-08, 6.76311733496419E-05, -4.23418981680198E-08,
        -8.15879467982627E-06, -4.23418981680198E-08, 6.61039240611986E-05;
    inertial_i[5] << 3.76561532469651E-05, 4.18325643595058E-08, -2.08841259371089E-07,
        4.18325643595058E-08, 4.00422747874576E-05, -1.50489542382975E-08,
        -2.08841259371089E-07, -1.50489542382975E-08, 3.97966719504564E-05;

    Eigen::Matrix<double, 3, 3> zero3x3;
    Eigen::Matrix<double, 3, 3> eye3x3;
    eye3x3.setIdentity();
    zero3x3.setZero();
    for (int i = 0; i < DOF; ++i)
    {
        this->inertial_G[i] <<
            inertial_i[i], zero3x3,
            zero3x3, this->links_mass(i)* eye3x3;
    }

    //DH
    this->T_dh_link[0] <<
        0, 1, 0, 0,
        -1, 0, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1;
    this->T_dh_link[1] <<
        0, 1, 0, 0,
        -1, 0, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1;
    this->T_dh_link[2] <<
        0, 0, 1, 0,
        0, -1, 0, 0,
        1, 0, 0, 0,
        0, 0, 0, 1;
    this->T_dh_link[3] <<
        0, 0, 1, 0,
        0, -1, 0, 0,
        1, 0, 0, 0,
        0, 0, 0, 1;
    this->T_dh_link[4] <<
        0, 1, 0, 0,
        0, 0, 1, 0,
        1, 0, 0, -0.075,
        0, 0, 0, 1;
    this->T_dh_link[5] <<
        0, 1, 0, 0,
        -1, 0, 0, 0,
        0, 0, 1, -0.08,
        0, 0, 0, 1;
    this->T_dh_link[6] <<
        0, 1, 0, 0,
        0, 0, 1, 0,
        1, 0, 0, 0,
        0, 0, 0, 1;
    for (int i = 0; i < DOF + 1; ++i)
    {
        this->T_dh_c[i] = this->T_C[i] * this->T_dh_link[i];
    }
    for (int i = 0; i < DOF; ++i)
    {
        iner_DH[i] = iner_trans(this->T_dh_c[i + 1], inertial_G[i].block(0, 0, 3, 3), links_mass(i));
    }
}

//解析urdf接口
void my_Robot::urdf_parser(std::string urdf_name, std::string links_names[], std::string axies_names[])
{
    // 输入参数为urdf文件，连杆名列表，关节名列表。
    // 任意一个现存的爬虫用的urdf解析器都可以读取M，inertialG, A
    std::ifstream urdfFile(urdf_name);
    if (!urdfFile.is_open()) {
        std::cerr << "Failed to open URDF file." << std::endl;
        exit(1);
    }
    std::string urdfContent((std::istreambuf_iterator<char>(urdfFile)), std::istreambuf_iterator<char>());
    // 2. 解析URDF文件内容
    pugi::xml_document doc;
    pugi::xml_parse_result result = doc.load_string(urdfContent.c_str());

    if (!result) {
        std::cerr << "XML parsing error: " << result.description() << std::endl;
        exit(1);
    }

    // 3. 获取每个连杆的质量和惯性矩阵
    int foundIndex = 0;
    Eigen::Matrix<double, 3, 3> inertia_Mat_temp;
    Eigen::Matrix<double, 6, 6> inertia_G_temp;
    Eigen::Matrix<double, 3, 3> zero3x3; zero3x3.setZero();
    Eigen::Matrix<double, 3, 3> one3x3; one3x3.setIdentity();
    Eigen::Matrix<double, 3, 1> transpose_temp;
    Eigen::Matrix<double, 3, 3> rotate_temp;
    Eigen::Matrix<double, 3, 1> rpy_temp;
    Eigen::Matrix<double, 4, 4> T_temp;
    Eigen::Matrix<double, 3, 1> screw_temp;

    pugi::xml_node robotNode = doc.child("robot");

    for (pugi::xml_node linkNode = robotNode.child("link"); linkNode; linkNode = linkNode.next_sibling("link"))
    {
        //查找连杆属性，包括连杆质量，惯性矩阵，质心相对mesh的位姿
        std::string link_name = std::string(linkNode.attribute("name").as_string());
        auto it = std::find(links_names, links_names + DOF + 1, link_name);
        if (it != links_names + DOF + 1)
        {
            foundIndex = std::distance(links_names, it);
        }
        else {
            std::cout << "can't find link: \"" + link_name + "\" in your input links names!" << std::endl;
            continue;
        }


        pugi::xml_node inertialNode = linkNode.child("inertial");
        if (inertialNode)
        {
            //质心偏移
            std::string offsetStr_transe = inertialNode.child("origin").attribute("xyz").as_string();
            transpose_temp = parse_str_vec(offsetStr_transe, 3);
            std::string offsetStr_rotate = inertialNode.child("origin").attribute("rpy").as_string();
            rpy_temp = parse_str_vec(offsetStr_rotate, 3);

            this->T_C[foundIndex] << rpyToRotationMatrix(rpy_temp), transpose_temp,
                0, 0, 0, 1;

            std::cout << "init T_C at: \"" + link_name + "\"" << std::endl;
            std::cout << this->T_C[foundIndex] << std::endl;

            //要求baselink放在第零位，但是不需要baselink的质量与质量分布。
            if (foundIndex == 0)
            {
                continue;
            }
            //连杆质量
            double mass = inertialNode.child("mass").attribute("value").as_double();
            this->links_mass(foundIndex - 1) = mass;
            //连杆质量分布
            inertia_Mat_temp(0, 0) = inertialNode.child("inertia").attribute("ixx").as_double();
            inertia_Mat_temp(1, 1) = inertialNode.child("inertia").attribute("iyy").as_double();
            inertia_Mat_temp(2, 2) = inertialNode.child("inertia").attribute("izz").as_double();
            // urdf如果是sw导出的，那么除了对角线意外的元素都要添加负号?
            inertia_Mat_temp(0, 1) = inertialNode.child("inertia").attribute("ixy").as_double();
            inertia_Mat_temp(0, 2) = inertialNode.child("inertia").attribute("ixz").as_double();
            inertia_Mat_temp(1, 2) = inertialNode.child("inertia").attribute("iyz").as_double();
            inertia_Mat_temp(1, 0) = inertia_Mat_temp(0, 1);
            inertia_Mat_temp(2, 0) = inertia_Mat_temp(0, 2);
            inertia_Mat_temp(2, 1) = inertia_Mat_temp(1, 2);
            inertia_G_temp << inertia_Mat_temp, zero3x3,
                zero3x3, this->links_mass(foundIndex - 1)* one3x3;
            this->inertial_G[foundIndex - 1] = inertia_G_temp;
            std::cout << "init inerG at: \"" + link_name + "\"" << std::endl;
            std::cout << this->inertial_G[foundIndex - 1] << std::endl;
        }
    }
    std::cout << "init link mass:" << std::endl;
    std::cout << my_Robot::links_mass.transpose() << std::endl;


    //查找关节属性，包括旋量轴，初始位姿
    for (pugi::xml_node linkNode = robotNode.child("joint"); linkNode; linkNode = linkNode.next_sibling("joint"))
    {
        std::string joint_name = std::string(linkNode.attribute("name").as_string());
        auto it = std::find(axies_names, axies_names + DOF, joint_name);
        if (it != axies_names + DOF)
        {
            foundIndex = std::distance(axies_names, it);
        }
        else {
            std::cout << "can't find joint: \"" + joint_name + "\" in your input axies names!" << std::endl;
            continue;
        }
        //相对位姿
        std::string offsetStr_transe = linkNode.child("origin").attribute("xyz").as_string();
        transpose_temp = parse_str_vec(offsetStr_transe, 3);
        std::string offsetStr_rotate = linkNode.child("origin").attribute("rpy").as_string();
        rpy_temp = parse_str_vec(offsetStr_rotate, 3);

        T_temp << rpyToRotationMatrix(rpy_temp), transpose_temp,
            0, 0, 0, 1;
        this->links_T[foundIndex] = inv_SE3(T_temp);

        //旋量轴
        std::string offsetStr_screw = linkNode.child("axis").attribute("xyz").as_string();
        screw_temp = parse_str_vec(offsetStr_screw, 3);

        Eigen::Matrix <double, 3, 1> zero3x1;
        Eigen::Matrix <double, 6, 1> cat_temp;
        zero3x1.setZero();
        cat_temp << screw_temp,
            zero3x1;
        this->links_B[foundIndex] = cat_temp;

    }
    this->links_T[DOF] << 1, 0, 0, 0,
        0, 1, 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1;


    //处理质心与连杆坐标系的差别
    //temp_T[i]: 连杆坐标系i相对于基坐标系0的位姿
    Eigen::Matrix<double, 4, 4> temp_T[DOF];
    Eigen::Matrix<double, 3, 1> temp_p[DOF];
    Eigen::Matrix<double, 3, 1> temp_w[DOF];
    Eigen::Matrix<double, 4, 4> cor_M[DOF + 1];
    for (int i = 0; i < DOF; ++i)
    {
        int j = 0;
        temp_T[i].setIdentity();
        while (j <= i)
        {
            temp_T[i] *= inv_SE3(this->links_T[j]);
            j += 1;
        }

        temp_p[i] = temp_T[i].block(0, 3, 3, 1);
        temp_w[i] = temp_T[i].block(0, 0, 3, 3) * this->links_B[i].block(0, 0, 3, 1);

        this->links_S[i] << temp_w[i], -vec2so3(temp_w[i]) * temp_p[i];
        // temp_T[i].block(0, 0, 3, 3).setIdentity();
        cor_M[i] = temp_T[i] * this->T_C[i + 1]; //此处不算基座，因此从i+1开始
    }

    cor_M[DOF].setIdentity();
    for (int i = 0; i < DOF + 1; ++i)
    {
        cor_M[DOF] *= inv_SE3(this->links_T[i]);
    }

    Eigen::Matrix<double, 4, 4> one4x4;
    one4x4.setIdentity();
    for (int i = 0; i < DOF + 1; ++i)
    {
        if (i == 0)
        {
            this->links_M[i] = one4x4 * cor_M[i];
            continue;
        }
        this->links_M[i] = inv_SE3(cor_M[i - 1]) * cor_M[i];
    }

    //处理旋量轴S到A
    Eigen::Matrix<double, 4, 4> fk_temp_M;
    fk_temp_M.setIdentity();
    for (int i = 0; i < DOF; ++i)
    {
        fk_temp_M *= this->links_M[i];
        links_A[i] = AD_SE3(inv_SE3(fk_temp_M)) * this->links_S[i];
    }

    std::cout << "parser urdf success!" << std::endl;
}
