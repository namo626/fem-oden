#include <vector>

#pragma once

extern int nodesCount;
extern std::vector<double> nodesX;
extern std::vector<double> nodesY;
extern int quadDegree;
extern double DPHI[3][3][3];
extern int* ID; // converts global node number to equation number in global matrix
extern int neq;  // number of global equations = nodesCount - dirichlet nodes
extern std::vector<int> dirichletNodes;  // list of nodes with specified essential BC
extern double* Q; // values at the essential BC, zero otherwise

/* Gaussian integration points and weights on the ref. triangle */
extern double GaussPoints[3][3];
extern double GaussWeights[3];

double ue(double x, double y);
double f(double x, double y);
void read_nodes();
void init_dirichlet();
void init_ID();
void init_dPhi();
void initialize();
