#pragma once
#include <cmath>
#include "Eigen/Dense"
#include "robot_basic_params.h"
#include <iostream>

Eigen::Matrix<double, DOF* DOF, 1> M_gluon(const double* parms, const double* q);
Eigen::Matrix<double, DOF, DOF> sym_M_gluon(const double* parms, const double* q);