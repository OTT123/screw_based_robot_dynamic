#include "./my_dynamic.h"
#include <windows.h>


int main()
{
    //测试
    my_Robot_dynamic robot;
    //ur5机械臂 
    //std::string link_name[DOF + 1] = { "base_link" , "shoulder_link" ,"upper_arm_link" ,"forearm_link" ,"wrist_1_link" ,"wrist_2_link" ,"wrist_3_link" };
    //std::string axis_name[DOF] = { "shoulder_pan_joint","shoulder_lift_joint","elbow_joint","wrist_1_joint","wrist_2_joint","wrist_3_joint" };
    //std::string urdf_name = ".\\urdf\\ur5.urdf";
    //gluon机械臂 
    std::string link_name[DOF + 1] = { "base_link" , "1_Link" ,"2_Link", "3_Link" ,"4_Link" ,"5_Link" ,"6_Link" };
    std::string axis_name[DOF] = { "axis_joint_1","axis_joint_2","axis_joint_3","axis_joint_4","axis_joint_5","axis_joint_6" };
    std::string urdf_name = ".\\urdf\\gluon.urdf";

    //puma560机械臂 
    //std::string link_name[DOF + 1] = { "base_link" , "Link1" ,"Link2" ,"Link3" ,"Link4" ,"Link5" ,"Link6" };
    //std::string axis_name[DOF] = { "J1","J2","J3","J4","J5","J6" };
    //std::string urdf_name = ".\\urdf\\puma560_robot.urdf";

    //sw导出的珞石机械臂urdf
    //std::string link_name[DOF + 1] = { "base_link" , "1_Link" ,"2_Link" ,"3_Link" ,"4_Link" ,"5_Link" ,"6_Link" };
    //std::string axis_name[DOF] = { "1_joint","2_joint","3_joint","4_joint","5_joint","6_joint" };
    //std::string urdf_name = ".\\urdf\\test.SLDASM.urdf";

    std::cout << "parsing robot urdf @: " + urdf_name << std::endl;
    robot.urdf_parser(urdf_name, link_name, axis_name);

    //关节变量
    Eigen::Matrix<double, DOF, 1> iter_q;
    Eigen::Matrix<double, DOF, 1> iter_q_dot;
    Eigen::Matrix<double, DOF, 1> iter_q_ddot;
    iter_q << 1, 1, 1, 1, 1, 1;
    iter_q_dot << 0.5, 0.8, 1.1, 0.4, 0.2, 0.5;
    iter_q_ddot << 0, 0, 0, 0, 0, 0;

    //重力常数
    Eigen::Matrix<double, 3, 1> gravity;
    gravity << 0, 0, -9.81;

    //动力学参数
    Eigen::Matrix<double, DOF, DOF> massMat;
    Eigen::Matrix<double, DOF, 1> corMat;
    Eigen::Matrix<double, DOF, 1> graMat;

    int count = pow(2,10);   //计算count次用于测量时间性能
    LARGE_INTEGER t1, t2, tc;
    std::cout << "***********************" << std::endl;
    QueryPerformanceFrequency(&tc);
    QueryPerformanceCounter(&t1);
    for (int i = 0; i < count; ++i)
    {
        massMat = robot.cal_massMat(iter_q);
    }
    std::cout << "massMat:" << std::endl;
    std::cout << massMat << std::endl;
    QueryPerformanceCounter(&t2);
    double time = (double)(t2.QuadPart - t1.QuadPart) / (double)tc.QuadPart;
    std::cout << "Calculate massMat "<< count <<" times costs: "<< time <<"s, frenquency:" << count/time <<" Hz" << std::endl; 

    std::cout << "***********************" << std::endl;
    QueryPerformanceFrequency(&tc);
    QueryPerformanceCounter(&t1);
    for (int i = 0; i < count; ++i)
    {
        corMat = robot.cal_corMat(iter_q, iter_q_dot);
    }
    std::cout << "corMat:" << std::endl;
    std::cout << corMat << std::endl;
    QueryPerformanceCounter(&t2);
    time = (double)(t2.QuadPart - t1.QuadPart) / (double)tc.QuadPart;
    std::cout << "Calculate corMat " << count << " times costs: " << time << "s, frenquency:" << count / time << " Hz" << std::endl; 

    std::cout << "***********************" << std::endl;
    QueryPerformanceFrequency(&tc);
    QueryPerformanceCounter(&t1);
    for (int i = 0; i < count; ++i)
    {
        graMat = robot.cal_graMat(iter_q, gravity);
    }
    std::cout << "graMat:" << std::endl;
    std::cout << graMat << std::endl;
    QueryPerformanceCounter(&t2);
    time = (double)(t2.QuadPart - t1.QuadPart) / (double)tc.QuadPart;
    std::cout << "Calculate graMat " << count << " times costs: " << time << "s, frenquency:" << count / time << " Hz" << std::endl;  
    return 0;
}