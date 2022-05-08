#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Element.h"
#include "Global.h"
//#include "omp.h"
//#include <mpi.h>
#include <petscksp.h>

using namespace std;
using namespace Eigen;
typedef Triplet<double> T;

/* Constructor for element class */
Element::Element(vector<int> nnodes) {
  //id = idd;
  nodes = nnodes;
  Ke.resize(3, 3);
  Fe.resize(3, 1);
  Ke.setZero();
  Fe.setZero();

}

void Element::initQ() {
  for (int i = 0; i < 3; i++) {
    Qe(i) = Q[nodes[i]];
  }
}

ostream &operator<<(ostream &os, const Element &elem)
{
  //os << "ID: " << elem.id << endl;
  os << "Nodes: " << endl;
  for (auto e: elem.nodes) {
    os << e << " ";
  }
  os << endl;
  os << "NodesX: " << endl;
  for (auto e: elem.nodes) {
    os << nodesX[e] << " ";
  }
  os << endl;
  return os;
}

/* Calculate K and F matrices for each element */
void Element::computeXY() {

  double x, y;
  /* For each of x, dx/dzeta, dx/deta */
  for (int i = 0; i < 3; i++) {
    /* For each quadrature point */
    for (int j = 0; j < quadDegree; j++) {
      /* Sum over the ref. shape functions */
      XGP[i][j] = 0;
      YGP[i][j] = 0;
      for (int k = 0; k < 3; k++) {
        x = nodesX[nodes[k]];
        y = nodesY[nodes[k]];
        XGP[i][j] += x * DPHI[k][i][j];
        YGP[i][j] += y * DPHI[k][i][j];
      }
    }
    //printf("x = %.2f\n", x);
  }

}

/* Compute the Jacobian at each Gauss point */
void Element::computeJ() {
  for (int i = 0; i < quadDegree; i++) {
    J[i] = XGP[1][i]*YGP[2][i] - XGP[2][i]*YGP[1][i];
  }
}

/* Compute derivatives of zeta and eta */
void Element::computeZE() {
  for (int i = 0; i < quadDegree; i++) {
    Zeta[0][i] = (1/J[i]) * YGP[2][i];
    Zeta[1][i] = (-1/J[i]) * XGP[2][i];

    Eta[0][i] = (-1/J[i]) * YGP[1][i];
    Eta[1][i] = (1/J[i]) * XGP[1][i];
  }
}

void Element::computePhi() {
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < quadDegree; j++) {
      PhiX[i][j] = DPHI[i][1][j]*Zeta[0][j] + DPHI[i][2][j]*Eta[0][j];
      PhiY[i][j] = DPHI[i][1][j]*Zeta[1][j] + DPHI[i][2][j]*Eta[1][j];
    }
  }
}

void Element::precompute() {
  computeXY();
  computeJ();
  computeZE();
  computePhi();
  initQ();
}

void Element::computeKF() {
  Ke.setZero();
  Fe.setZero();

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      /* Sum over gauss points */
      for (int k = 0; k < 3; k++) {
        Ke(i,j) += J[k]*GaussWeights[k]*(PhiX[i][k]*PhiX[j][k] + PhiY[i][k]*PhiY[j][k]);
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    for (int k = 0; k < quadDegree; k++) {
      Fe(i) += f(XGP[0][k], YGP[0][k])*DPHI[i][0][k]*J[k]*GaussWeights[k];
    }
  }

  /* Add the essential BC */
  Fe = Fe - Ke*Qe;
}

/* Assemble elemental matrices to global matrix via ID array */
void Element::assemble(Mat &K, Vec &F) {
  int ii, jj;

  vector<int> idxm;
  vector<int> idxn;
  vector<double> values;
  vector<int> inds;
  vector<double> fvals;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ii = ID[nodes[i]];
      jj = ID[nodes[j]];
      if (ii != -1 && jj != -1) {
        //K(ii,jj) += Ke(i,j);
        //t.push_back(Triplet<double>(ii, jj, Ke(i,j)));
        //#pragma omp atomic
        //K.coeffRef(ii,jj) += Ke(i,j);
        // idxm.push_back(ii);
        // idxn.push_back(jj);
        // values.push_back(Ke(i,j));
        MatSetValues(K, 1, &ii, 1, &jj, &Ke(i,j), ADD_VALUES);
      }
    }
  }
  // for (int i = 0; i < values.size(); i++) {
  //   cout << idxm[i] << endl;
  // }
  //cout << idxm.size() << " " << idxn.size() << " " << values.size() << endl;
  for (int i = 0; i < 3; i++) {
    ii = ID[nodes[i]];
    if (ii != -1) {
      //F(ii) += Fe(i);
      // inds.push_back(ii);
      // fvals.push_back(Fe(i));
      VecSetValues(F, 1, &ii, &Fe(i), ADD_VALUES);
    }
  }
  //printf("inds.size() = %d\n", inds.size());
}

/* Initialize the vector of elements */
void read_elems(vector<Element> &elements) {
  int elem = 0;
  vector<int> row;
  string line, word;
  fstream file("mesh.csv", ios::in);

  if (file.is_open()) {
    while(getline(file, line)) {
      row.clear();
      stringstream str(line);

      while(getline(str, word, ',')) {
        row.push_back(stoi(word));
      }
      elements.push_back(Element(row));
      elem++;
    }
  } else {
    cout << "File doesn't exist" << endl;
  }

}

/* MPI IO for parallel reading of the element list */
vector<Element> read_elems_par(char* filename) {
  MPI_File fh;
  MPI_Status status;

  int rank, nprocs;
  int bufsize, nints, total, numrows;
  MPI_Offset filesize;

  int* buf;
  int rows_per_proc, elems_per_proc, remainder, offset;

  // amount of entries per row
  int datum = 3;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &fh);
  MPI_File_get_size(fh, &filesize);

  total = filesize / sizeof(int);
  numrows = total / datum;
  rows_per_proc = numrows / nprocs; //floor division

  // if (rank == 0) {
  //   printf("Total data = %d\n", total);
  //   printf("Numrows = %d\n", numrows);
  // }

  // last rank gets the remainder
  remainder = numrows - (nprocs-1)*rows_per_proc;
  if (rank == nprocs-1) {
    elems_per_proc = remainder;
  } else {
    elems_per_proc = rows_per_proc;
  }

  bufsize = datum * sizeof(int) * elems_per_proc;
  nints   = bufsize/sizeof(int);
  buf = new int[nints];

  offset = datum * sizeof(int) * rows_per_proc;

  MPI_File_seek(fh, rank*offset, MPI_SEEK_SET);
  MPI_File_read(fh, buf, nints, MPI_INT, &status);
  MPI_File_close(&fh);

  // Now initialize the elements
  vector<int> nnodes(3);
  vector<Element> elems;

  for (int i = 0; i < elems_per_proc; i++) {
    nnodes = {buf[datum*i], buf[datum*i+1], buf[datum*i+2]};
    elems.push_back(Element(nnodes));
  }
  return elems;

}
