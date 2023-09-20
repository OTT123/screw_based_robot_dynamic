#pragma once
#pragma once
#include <iostream>
#include "Eigen/Dense"
#include "pugixml.hpp"
#include <string>
#include <stdlib.h>
#include <fstream>
#include <sstream>
#include <algorithm>
//自由度
#include "robot_basic_params.h"

//SE3 求逆矩阵
Eigen::MatrixXd inv_SE3(Eigen::MatrixXd T);

//@output :[p]
Eigen::MatrixXd vec2so3(Eigen::MatrixXd p);

//SE3的伴随矩阵
//@param T:4x4 SE3
//@output :ADT 6x6
Eigen::MatrixXd AD_SE3(Eigen::MatrixXd T);

//se3的伴随矩阵
//@param S:6x1 screw axis
//@output :[s] or ads,6x6
Eigen::MatrixXd ad_screw(Eigen::MatrixXd S);

//旋转矩阵指数运算
//输入螺旋轴和旋转角度，输出对应的SE3
//@param w:3x1旋转轴
//@param theta:旋转角度
//@output :旋转矩阵R SO3
Eigen::MatrixXd exp_so3_2_SO3(Eigen::MatrixXd w, double theta);

//螺旋轴指数运算
//输入螺旋轴和旋转角度，输出对应的SE3
//@param s:6x1螺旋轴[w v]
//@param theta:旋转角度
//@output :齐次变换矩阵 SE3
Eigen::MatrixXd exp_twist_2_SE3(Eigen::MatrixXd S, double theta);

//解析字符串向量，以空格分开
Eigen::VectorXd parse_str_vec(std::string scr, int num);

//惯性矩阵转换为空间惯量矩阵
//@param inertial:3x1
//@output inerG:6x6
Eigen::MatrixXd Inertial2inerG(Eigen::MatrixXd inertial_i, double m);

//雅可比矩阵，测试中，不保真
Eigen::Matrix<double, 6, 1> Jaccobi_s_i(Eigen::Matrix<double, 6, 1> theta, Eigen::Matrix<double, 6, 1> screw[], int col_i);
Eigen::Matrix<double, 6, 6> Jaccobi_s(Eigen::Matrix<double, 6, 1> theta, Eigen::Matrix<double, 6, 1> screw[], int n);

//RPY->RotateMat
Eigen::Matrix<double, 3, 3> rpyToRotationMatrix(const Eigen::Vector3d& rpy);

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
Eigen::Matrix<double, 3, 3> iner_trans(Eigen::Matrix<double, 4, 4> T, Eigen::Matrix<double, 3, 3> iner, double m);


/// class: 机器人参数
/// urdf中连杆坐标系（mesh）大多数情况下与质心坐标系不重合。urdf的初始相对位姿和关节旋转轴都是定义在mesh坐标系中的。
/// 旋量理论动力学中所有量都在质心坐标系中讨论。
/// 因此实现了方法urdf_parser，它可以将urdf中的信息转换为旋量动力学所需的量
/// 经过反复验证，现在算法已经和KDL，Pinocchio，matlab robotic tool box的结果相同。
/// modern robotic 的github仓库中给出的算法是正确的，但是它的算法api所需的参数和书上相比有两点不同
/// 1. 他需要质心i相对于基坐标系的位姿，而非质心i-1相对于质心i的位姿
/// 2. 他需要关节旋量轴在基坐标系的表示，而非关节i旋量轴在质心坐标系i的表示。
/// 但是，modern robotic github上推荐的仓库mr_urdf_loader （作者： MinchangSung0223）在解析urdf的代码中有一段意义不明的语句，我将其修改，发现
/// KDL，Pinocchio，matlab robotic tool box， 我的算法 以及 modern robotic， 总共五方算法结果一致。问题解决，是mr_urdf_loader源码有问题。
/// 计算一次惯性矩阵(最复杂)耗时150us左右(60kHz)，可以达到实时控制要求
class my_Robot
{
public:
    //初始状态下mesh坐标系的位姿
    Eigen::Matrix<double, 4, 4> links_T[DOF + 1];
    //质心相对于mesh坐标系的SE3
    Eigen::Matrix<double, 4, 4> T_C[DOF + 1];
    //质心坐标系的初始相对位姿
    // 质心坐标系i相对于质心坐标系i-1的相对位姿
    // 迭代算法时注意求逆
    Eigen::Matrix<double, 4, 4> links_M[DOF + 1];
    //关节i的旋量轴在i坐标系的定义
    Eigen::Matrix<double, 6, 1> links_A[DOF];
    Eigen::Matrix<double, 6, 1> links_B[DOF];
    //关节i的旋量轴在基坐标系的定义
    Eigen::Matrix<double, 6, 1> links_S[DOF];
    //连杆i质量(不含基座)
    Eigen::MatrixXd links_mass = Eigen::MatrixXd(DOF, 1);
    //连杆的惯性矩阵
    Eigen::Matrix<double, 6, 6> inertial_G[DOF];

    //DH参数中的i坐标系为基准,link i 坐标系的SE3(含基座)
    Eigen::Matrix<double, 4, 4> T_dh_link[DOF + 1];
    //DH参数中的i坐标系为基准,link i 质心坐标系的SE3(含基座)
    Eigen::Matrix<double, 4, 4> T_dh_c[DOF + 1];
    //DH坐标系下各连杆的质量分布(含基座)
    Eigen::Matrix<double, 3, 3> iner_DH[DOF + 1];

    my_Robot();
    //初始化了Gluon机械臂参数,这个函数允许用户按照模板修改源码，以实现自定义机械臂参数的设置
    void set_param_from_src_code();
    //urdf解析器
    void urdf_parser(std::string urdf_name, std::string links_names[], std::string axies_names[]);
};

//动力学算法
class my_Robot_dynamic :public my_Robot
{
private:
    Eigen::Matrix<double, 6 * DOF, DOF> stack_A;
    Eigen::Matrix<double, 6 * DOF, 6 * DOF> stack_inerG;
    Eigen::Matrix<double, 6 * DOF, 6 * DOF> stack_ad_V;
    Eigen::Matrix<double, 6 * DOF, 6 * DOF> stack_ad_A_theta_dot;
    Eigen::Matrix<double, 6 * DOF, 6 * DOF> stack_W;
    Eigen::Matrix<double, 6 * DOF, 6 * DOF> stack_L;
    Eigen::Matrix<double, 6 * DOF, 1> stack_V_base;
    Eigen::Matrix<double, 6 * DOF, 1> stack_V_base_dot;
    Eigen::Matrix<double, 6 * DOF, 1> stack_F_tip;
public:
    Eigen::Matrix<double, DOF, DOF> stack_massMat;
    Eigen::Matrix<double, DOF, 1> stack_corMat;
    Eigen::Matrix<double, DOF, 1> stack_graMat;

    my_Robot_dynamic();
    //逆向动力学
    Eigen::Matrix<double, DOF, 1> forword_dyn(
        Eigen::Matrix<double, DOF, 1> iter_q,
        Eigen::Matrix<double, DOF, 1> iter_q_dot,
        Eigen::Matrix<double, DOF, 1> iter_q_ddot,
        Eigen::Matrix<double, 3, 1> gravity);
    //惯性矩阵
    Eigen::Matrix<double, DOF, DOF> cal_massMat(
        Eigen::Matrix<double, DOF, 1> iter_q);
    //科里奥利力与向心力矩阵
    Eigen::Matrix<double, DOF, 1> cal_corMat(
        Eigen::Matrix<double, DOF, 1> iter_q,
        Eigen::Matrix<double, DOF, 1> iter_q_dot);
    //重力矩阵
    Eigen::Matrix<double, DOF, 1> cal_graMat(
        Eigen::Matrix<double, DOF, 1> iter_q,
        Eigen::Matrix<double, 3, 1> gravity);
};
