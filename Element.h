#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <sstream>
#include <fstream>
#include <iostream>
#include <stdio.h>

#pragma once
typedef Eigen::Triplet<double> T;

class Element {
public:
  Element(int id, std::vector<int> nodes);
  // Id
  int id;
  // List of nodes in counterclockwise direction
  std::vector<int> nodes;

  void computeXY();
  void computeJ();
  void computeZE();
  void computePhi();
  void initQ();
  void precompute();

  void computeKF();
  //void assemble(std::vector<T> &t, Eigen::MatrixXd &F);
  void assemble(Eigen::SparseMatrix<double> &K, Eigen::MatrixXd &F);

  /* Gaussian points in global coordinates */
  double XGP[3][3];
  double YGP[3][3];

  /* Jacobian at each Gauss point */
  double J[3];

  /* Derivatives of ref. coordinates, e.g. dzeta/dx  at each Gauss point */
  double Zeta[2][3];
  double Eta[2][3];

  /* Derivatives of elemental shape functions in x and y */
  double PhiX[3][3];
  double PhiY[3][3];

  Eigen::MatrixXd Ke;
  Eigen::MatrixXd Fe;
  Eigen::Vector3d Qe;

};

void read_elems(std::vector<Element> &elements);
