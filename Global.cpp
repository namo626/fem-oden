#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include "Global.h"
#include <math.h>
#include <algorithm>
#include <mpi.h>

using namespace std;
int nodesCount;
vector<double> nodesX;
vector<double> nodesY;
int quadDegree;
double DPHI[3][3][3];
int* ID; // converts global node number to equation number in global matrix
int neq;  // number of global equations = nodesCount - dirichlet nodes
vector<int> dirichletNodes;  // list of nodes with specified essential BC
double* Q; // values at the essential BC, zero otherwise

double GaussPoints[3][3] = {
  {0.5, 0},
  {0, 0.5},
  {0.5, 0.5}
};

double GaussWeights[3] = {1./3, 1./3, 1./3};

/* Test solution */
double ue(double x, double y) {
  //return 1 + pow(x,2) + 2*pow(y,2);
  return exp(-(x*x/9. + y*y/9.));
}

/* Load function */
double f(double x, double y) {
  //return -6.0;
  return -(4/81.)*exp(-(x*x/9.+y*y/9.))*(-9 + x*x + y*y);
}

/* Parallel read */
void read_nodes_par() {
  MPI_File fh;
  MPI_Status status;

  int rank, nprocs;
  int bufsize, nints, total, numrows;
  MPI_Offset filesize;

  double* buf;
  double* nodesX_local;
  double* nodesY_local;
  int* counts;
  int* disp;

  int rows_per_proc, nodes_per_proc, remainder, offset;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  counts = new int[nprocs];
  disp = new int[nprocs];

  MPI_File_open(MPI_COMM_WORLD, "nodes", MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  MPI_File_get_size(fh, &filesize);

  total = filesize / sizeof(double);
  numrows = total / 2;
  rows_per_proc = floor((double)numrows / nprocs);

  if (rank == 0) {
    printf("Total data = %d\n", total);
    printf("Numrows = %d\n", numrows);
  }
  nodesCount = numrows;
  nodesX.resize(nodesCount);
  nodesY.resize(nodesCount);

  // last rank gets the remainder
  remainder = numrows - (nprocs-1)*rows_per_proc;
  if (rank == nprocs-1) {
    nodes_per_proc = remainder;
  } else {
    nodes_per_proc = rows_per_proc;
  }

  if (rank == 0) {
    printf("rows_per_proc = %d\n", rows_per_proc);
    printf("remainder = %d\n", remainder);
  }
  bufsize = 2 * sizeof(double) * nodes_per_proc;
  nints   = bufsize/sizeof(double);
  buf = new double[nints];

  offset = 2 * sizeof(double) * rows_per_proc;

  MPI_File_seek(fh, rank*offset, MPI_SEEK_SET);
  MPI_File_read(fh, buf, nints, MPI_DOUBLE, &status);
  MPI_File_close(&fh);

  nodesX_local = new double[nodes_per_proc];
  nodesY_local = new double[nodes_per_proc];

  if (rank == 0) {
    printf("Finished initializing\n");
  }
  // convert the binary array to nodesX and nodesY
  for (int i = 0; i < nodes_per_proc; i++) {
    nodesX_local[i] = buf[2*i];
    nodesY_local[i] = buf[2*i+1];
  }

  // if (rank == 1) {
  //   //printf("Finished allocating local arrays\n");
  //   for (int i = 0; i < nodes_per_proc; i++) {
  //     printf("%.1f\n", nodesX_local[i]);
  //     //printf("%.1f\n", buf[i]);
  //   }
  // }
  // fill the count array for each rank
  for (int i = 0; i < nprocs; i++) {
    counts[i] = rows_per_proc;
    disp[i] = i * rows_per_proc;
  }
  counts[nprocs-1] = remainder;

  // if (rank == 0) {
  //   for (int i = 0; i < nprocs; i++) {
  //     cout << disp[i] << endl;
  //   }
  // }
  // gather all the local nodelists from other ranks
  MPI_Allgatherv(nodesX_local, nodes_per_proc, MPI_DOUBLE, &nodesX[0], counts, disp, MPI_DOUBLE, MPI_COMM_WORLD);
  MPI_Allgatherv(nodesY_local, nodes_per_proc, MPI_DOUBLE, &nodesY[0], counts, disp, MPI_DOUBLE, MPI_COMM_WORLD);

  delete [] nodesX_local;
  delete [] nodesY_local;
  delete [] buf;
  delete []  counts;
  delete [] disp;
}

/* Initializes nodesCount, nodesX, nodesY */
void read_nodes() {
  int nodeID = 0;
  vector<double> row;
  string line, word;
  fstream file("nodes.csv", ios::in);

  if (file.is_open()) {
    while(getline(file, line)) {
      row.clear();
      stringstream str(line);

      while(getline(str, word, ',')) {
        row.push_back(stof(word));
      }
      nodesX.push_back(row[0]);
      nodesY.push_back(row[1]);

      nodeID++;
    }
  } else {
    cout << "File doesn't exist" << endl;
  }

  nodesCount = nodeID;
}

/* Initialize list of boundary nodes */
void init_dirichlet() {
  int n = (int) sqrt(nodesCount);
  // bottom
  for (int i = 0; i <= n-1; i++) {
    dirichletNodes.push_back(i);
  }

  // top
  for (int i = nodesCount-n; i <= nodesCount-1; i++) {
    dirichletNodes.push_back(i);
  }

  // left
  for (int i = 1; i <= n-2; i++) {
    dirichletNodes.push_back(i*n);
  }

  // right
  for (int i = 1; i <= n-2; i++) {
    dirichletNodes.push_back((i+1)*n-1);
  }


  neq = nodesCount - dirichletNodes.size();
}


/* Initialize ID array */
void init_ID() {
  ID = new int[nodesCount];
  Q = new double[nodesCount];

  int count = 0;
  for (int i = 0; i < nodesCount; i++) {
    if(find(dirichletNodes.begin(), dirichletNodes.end(), i) != dirichletNodes.end()) {
      ID[i] = -1;
      Q[i] = ue(nodesX[i], nodesY[i]);
    } else {
      ID[i] = count;
      Q[i] = 0;
      count++;
    }
  }
}


/* Evaluate phi and its derivatives at ref. gaussian points */
void init_dPhi() {
  // Phi_1: 1 - zeta - eta
  for (int i = 0; i < quadDegree; i++) {
    double zeta = GaussPoints[i][0];
    double eta = GaussPoints[i][1];

    DPHI[0][0][i] = 1 - zeta - eta;
    DPHI[0][1][i] = -1;
    DPHI[0][2][i] = -1;

    /* Phi_2: zeta */
    DPHI[1][0][i] = zeta;
    DPHI[1][1][i] = 1;
    DPHI[1][2][i] = 0;

    /* Phi_3: eta */
    DPHI[2][0][i] = eta;
    DPHI[2][1][i] = 0;
    DPHI[2][2][i] = 1;
  }
}


/* Initialize everything */
void initialize() {
  quadDegree = 3;
  read_nodes();
  init_dirichlet();
  init_ID();
  init_dPhi();
}
