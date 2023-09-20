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
//���ɶ�
#include "robot_basic_params.h"

//SE3 �������
Eigen::MatrixXd inv_SE3(Eigen::MatrixXd T);

//@output :[p]
Eigen::MatrixXd vec2so3(Eigen::MatrixXd p);

//SE3�İ������
//@param T:4x4 SE3
//@output :ADT 6x6
Eigen::MatrixXd AD_SE3(Eigen::MatrixXd T);

//se3�İ������
//@param S:6x1 screw axis
//@output :[s] or ads,6x6
Eigen::MatrixXd ad_screw(Eigen::MatrixXd S);

//��ת����ָ������
//�������������ת�Ƕȣ������Ӧ��SE3
//@param w:3x1��ת��
//@param theta:��ת�Ƕ�
//@output :��ת����R SO3
Eigen::MatrixXd exp_so3_2_SO3(Eigen::MatrixXd w, double theta);

//������ָ������
//�������������ת�Ƕȣ������Ӧ��SE3
//@param s:6x1������[w v]
//@param theta:��ת�Ƕ�
//@output :��α任���� SE3
Eigen::MatrixXd exp_twist_2_SE3(Eigen::MatrixXd S, double theta);

//�����ַ����������Կո�ֿ�
Eigen::VectorXd parse_str_vec(std::string scr, int num);

//���Ծ���ת��Ϊ�ռ��������
//@param inertial:3x1
//@output inerG:6x6
Eigen::MatrixXd Inertial2inerG(Eigen::MatrixXd inertial_i, double m);

//�ſɱȾ��󣬲����У�������
Eigen::Matrix<double, 6, 1> Jaccobi_s_i(Eigen::Matrix<double, 6, 1> theta, Eigen::Matrix<double, 6, 1> screw[], int col_i);
Eigen::Matrix<double, 6, 6> Jaccobi_s(Eigen::Matrix<double, 6, 1> theta, Eigen::Matrix<double, 6, 1> screw[], int n);

//RPY->RotateMat
Eigen::Matrix<double, 3, 3> rpyToRotationMatrix(const Eigen::Vector3d& rpy);

/// <summary>
/// �����ֲ��ڲ�ͬ����ϵ��ת��
/// </summary>
/// <param name="T">
/// B����ϵ ����� A����ϵ ��λ��
/// </param>
/// <param name="iner">
/// ������B����ϵ�������ֲ� 3x3
/// </param>
/// <param name="m">
/// ��������
/// </param>
/// <returns>
/// ������A����ϵ�������ֲ�
/// </returns>
Eigen::Matrix<double, 3, 3> iner_trans(Eigen::Matrix<double, 4, 4> T, Eigen::Matrix<double, 3, 3> iner, double m);


/// class: �����˲���
/// urdf����������ϵ��mesh����������������������ϵ���غϡ�urdf�ĳ�ʼ���λ�˺͹ؽ���ת�ᶼ�Ƕ�����mesh����ϵ�еġ�
/// �������۶���ѧ��������������������ϵ�����ۡ�
/// ���ʵ���˷���urdf_parser�������Խ�urdf�е���Ϣת��Ϊ��������ѧ�������
/// ����������֤�������㷨�Ѿ���KDL��Pinocchio��matlab robotic tool box�Ľ����ͬ��
/// modern robotic ��github�ֿ��и������㷨����ȷ�ģ����������㷨api����Ĳ�����������������㲻ͬ
/// 1. ����Ҫ����i����ڻ�����ϵ��λ�ˣ���������i-1���������i��λ��
/// 2. ����Ҫ�ؽ��������ڻ�����ϵ�ı�ʾ�����ǹؽ�i����������������ϵi�ı�ʾ��
/// ���ǣ�modern robotic github���Ƽ��Ĳֿ�mr_urdf_loader �����ߣ� MinchangSung0223���ڽ���urdf�Ĵ�������һ�����岻������䣬�ҽ����޸ģ�����
/// KDL��Pinocchio��matlab robotic tool box�� �ҵ��㷨 �Լ� modern robotic�� �ܹ��巽�㷨���һ�¡�����������mr_urdf_loaderԴ�������⡣
/// ����һ�ι��Ծ���(���)��ʱ150us����(60kHz)�����Դﵽʵʱ����Ҫ��
class my_Robot
{
public:
    //��ʼ״̬��mesh����ϵ��λ��
    Eigen::Matrix<double, 4, 4> links_T[DOF + 1];
    //���������mesh����ϵ��SE3
    Eigen::Matrix<double, 4, 4> T_C[DOF + 1];
    //��������ϵ�ĳ�ʼ���λ��
    // ��������ϵi�������������ϵi-1�����λ��
    // �����㷨ʱע������
    Eigen::Matrix<double, 4, 4> links_M[DOF + 1];
    //�ؽ�i����������i����ϵ�Ķ���
    Eigen::Matrix<double, 6, 1> links_A[DOF];
    Eigen::Matrix<double, 6, 1> links_B[DOF];
    //�ؽ�i���������ڻ�����ϵ�Ķ���
    Eigen::Matrix<double, 6, 1> links_S[DOF];
    //����i����(��������)
    Eigen::MatrixXd links_mass = Eigen::MatrixXd(DOF, 1);
    //���˵Ĺ��Ծ���
    Eigen::Matrix<double, 6, 6> inertial_G[DOF];

    //DH�����е�i����ϵΪ��׼,link i ����ϵ��SE3(������)
    Eigen::Matrix<double, 4, 4> T_dh_link[DOF + 1];
    //DH�����е�i����ϵΪ��׼,link i ��������ϵ��SE3(������)
    Eigen::Matrix<double, 4, 4> T_dh_c[DOF + 1];
    //DH����ϵ�¸����˵������ֲ�(������)
    Eigen::Matrix<double, 3, 3> iner_DH[DOF + 1];

    my_Robot();
    //��ʼ����Gluon��е�۲���,������������û�����ģ���޸�Դ�룬��ʵ���Զ����е�۲���������
    void set_param_from_src_code();
    //urdf������
    void urdf_parser(std::string urdf_name, std::string links_names[], std::string axies_names[]);
};

//����ѧ�㷨
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
    //������ѧ
    Eigen::Matrix<double, DOF, 1> forword_dyn(
        Eigen::Matrix<double, DOF, 1> iter_q,
        Eigen::Matrix<double, DOF, 1> iter_q_dot,
        Eigen::Matrix<double, DOF, 1> iter_q_ddot,
        Eigen::Matrix<double, 3, 1> gravity);
    //���Ծ���
    Eigen::Matrix<double, DOF, DOF> cal_massMat(
        Eigen::Matrix<double, DOF, 1> iter_q);
    //���������������������
    Eigen::Matrix<double, DOF, 1> cal_corMat(
        Eigen::Matrix<double, DOF, 1> iter_q,
        Eigen::Matrix<double, DOF, 1> iter_q_dot);
    //��������
    Eigen::Matrix<double, DOF, 1> cal_graMat(
        Eigen::Matrix<double, DOF, 1> iter_q,
        Eigen::Matrix<double, 3, 1> gravity);
};
