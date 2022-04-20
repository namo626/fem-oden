#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <string>

#include "Global.h"
#include "Element.h"
#include "Post.h"

using namespace std;
using namespace Eigen;
typedef Triplet<double> T;

int main() {
  initialize();

  vector<Element> elements;
  read_elems(elements);

  /* Global matrices */
  MatrixXd F;
  SparseMatrix<double> K(neq, neq);
  vector<T> t;
  F.resize(neq, 1);
  F.setZero();

  cout << "Number of nodes: " << nodesCount << endl;
  cout << "Number of equations: " << neq << endl;

  for (auto e: elements) {
    e.precompute();
    e.computeKF();
    e.assemble(t, F);
  }
  K.setFromTriplets(t.begin(), t.end());
  cout << "Finished assembly" << endl;

  /* Solve the global equation */
  SimplicialLDLT<SparseMatrix<double>> solver(K);
  MatrixXd sol = solver.solve(F);

  vector<double> vals = assignNodes(sol);
  write_solution(vals);

  return 0;
}
