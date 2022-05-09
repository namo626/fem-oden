#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <string>
#include <ctime>
#include "omp.h"
//#include <mpi.h>
#include <petscksp.h>

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

int main(int argc,char **args) {

  char help[] = "help";

  PetscInitialize(&argc,&args,(char*)0,help);
  //MPI_Init(NULL, NULL);

  int rank;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);


  // timer vars
  double i1, i2;
  double t1, t2;

  i1 = MPI_Wtime();
  initialize();

  vector<Element> elements;
  elements = read_elems_par((char*)"mesh");
  // for (auto e: elements) {
  //   e.precompute();
  //   e.computeKF();
  // }
  i2 = MPI_Wtime();


  Element e = elements[2];
  // e.precompute();
  // if (rank == 1) {
  //   //cout << elements[2].Fe << endl;
  //   //cout << elements[2].J[1] << endl;
  //   //cout << DPHI[1][1][1] << endl;
  //   for (int i = 0; i < 3; i++) {
  //     for (int j = 0; j < 3; j++) {
  //       printf("%.2f ", e.YGP[i][j]);
  //     }
  //   }
  //   cout << endl;

  //   cout << e.Ke << endl;
  //   cout << e.Fe << endl;
  // }

  if (rank == 0) {
    cout << "Number of nodes: " << nodesCount << endl;
    cout << "Number of elements: " << elements.size() << endl;
    cout << "Number of equations: " << neq << endl;
  }

  /* Create and assemble the global matrices in PETSC */
  Vec F;
  VecCreateMPI(MPI_COMM_WORLD, PETSC_DECIDE, neq, &F);
  VecSet(F, 0.);

  Mat K;
  MatCreateAIJ(MPI_COMM_WORLD,PETSC_DECIDE,PETSC_DECIDE, neq, neq, 7, NULL, 7, NULL, &K);
  //MatCreate(MPI_COMM_WORLD, &K);
  // MatCreateDense(MPI_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE, neq,neq, NULL, &K);
  MatZeroEntries(K);

  // for (auto e: elements) {
  //   e.assemble(K, F);
  // }

  //MatView(K, PETSC_VIEWER_STDOUT_WORLD);
  //i1 = MPI_Wtime();
  for (int i = 0; i < elements.size(); i++) {
    elements[i].precompute();
    elements[i].computeKF();
    elements[i].assemble(K,F);
  }


  VecAssemblyBegin(F);
  MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY);

  VecAssemblyEnd(F);
  MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY);

  i2 = MPI_Wtime();
  if (rank == 0) {
    cout << "Total time: " << i2-i1 << endl;
  }
  // debug
  //VecView(F, PETSC_VIEWER_STDOUT_WORLD);
  //MatView(K, PETSC_VIEWER_STDOUT_WORLD);

  /* Solve the system */
  i1 = MPI_Wtime();
  Vec x;
  VecCreateMPI(MPI_COMM_WORLD, PETSC_DECIDE, neq, &x);

  KSP ksp;
  KSPCreate(MPI_COMM_WORLD, &ksp);
  KSPSetOperators(ksp, K, K);
  KSPSetFromOptions(ksp);
  KSPSolve(ksp,F,x);
  i2 = MPI_Wtime();

  if (rank == 0) {
    cout << "Solving time: " << i2 - i1 << endl;
  }

  /* Write solution locally, then gather and write to file */

  Vec vout;
  VecScatter ctx;
  double* sol;
  VecScatterCreateToZero(x,&ctx,&vout);
  VecScatterBegin(ctx,x,vout,INSERT_VALUES,SCATTER_FORWARD);
  VecScatterEnd(ctx,x,vout,INSERT_VALUES,SCATTER_FORWARD);

  if (rank == 0) {
    VecGetArray(vout, &sol);
    vector<double> vals = assignNodes(sol);
    write_solution(vals);
  }

  // PetscViewer viewer;
  // PetscViewerASCIIOpen(MPI_COMM_WORLD, (char*)"output.txt", &viewer);
  // for (int i = 0; i < nodesCount; i++) {
  //   PetscViewerASCIIPrintf(viewer,"%.3f %.3f %.3f", nodesX[i], nodesY[i], sol[i]);
  // }



  PetscFinalize();
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
