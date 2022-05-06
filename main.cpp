#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <string>
#include <ctime>
#include "omp.h"
#include <mpi.h>

#include "Global.h"
#include "Element.h"
#include "Post.h"

using namespace std;
using namespace Eigen;
typedef Triplet<double> T;

double getTime(double i1, double i2) {
  double t1 = (i2 - i1);
  return t1;
}

/* Choosing between serial or parallel */
double timeit() {
#ifdef _OPENMP
  return omp_get_wtime();
#else
  return (double) clock() / CLOCKS_PER_SEC ;
#endif
}

int main() {

  MPI_Init(NULL, NULL);

  int rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  // timer vars
  double i1, i2;
  double t1, t2;

  i1 = MPI_Wtime();
  read_nodes_par();
  i2 = MPI_Wtime();

  if (rank == 0) {
    printf("%.3f\n", i2-i1);
  }
  // if (rank == 3) {
  //   for (int i = 0; i < nodesCount; i++) {
  //     printf("%d, %.1f\n", i, nodesX[i]);
  //   }
  // }
  cout << "FInished reading" << endl;
  // initialize();

  MPI_Finalize();
  // vector<Element> elements;
  // read_elems(elements);

  // /* Global matrices */
  // MatrixXd F;
  // SparseMatrix<double> K(neq, neq);
  // K.reserve(VectorXi::Constant(neq,7));
  // F.resize(neq, 1);
  // F.setZero();

  // cout << "Number of nodes: " << nodesCount << endl;
  // cout << "Number of elements: " << elements.size() << endl;
  // cout << "Number of equations: " << neq << endl;

  // i1 = timeit();
  // #pragma omp parallel for
  // for (auto e: elements) {
  //   e.precompute();
  //   e.computeKF();
  //   //e.assemble(t, F);
  //   e.assemble(K,F);
  // }
  // cout << "Finished assembly" << endl;
  // i2 = timeit();
  // t1 = getTime(i1, i2);

  // /* Solve the global equation */
  // i1 = timeit();
  // SimplicialLDLT<SparseMatrix<double> > solver(K);
  // //BiCGSTAB<SparseMatrix<double,RowMajor> > solver(K);
  // //LeastSquaresConjugateGradient<SparseMatrix<double>> solver;
  // //solver.compute(K);
  // //solver.setTolerance(1e-6);
  // MatrixXd sol = solver.solve(F);
  // i2 = timeit();
  // t2 = getTime(i1, i2);
  // // cout << "#iterations:     " << solver.iterations() << endl;
  // // cout << "estimated error: " << solver.error()      << endl;
  // vector<double> vals = assignNodes(sol);
  // write_solution(vals);

  // printf("Assembly: %.3f\n", t1);
  // printf("Solving : %.3f\n", t2);


  return 0;
}
