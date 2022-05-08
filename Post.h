#include <vector>
#include <Eigen/Dense>

#include "Global.h"

#pragma once

//std::vector<double> assignNodes(Eigen::MatrixXd &sol);
std::vector<double> assignNodes(double* sol);

void write_solution(std::vector<double> &vals);
