#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include "Element.h"
#include "Global.h"

using namespace std;
using namespace Eigen;
typedef Triplet<double> T;

/* Constructor for element class */
Element::Element(int idd, vector<int> nnodes) {
  id = idd;
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
  os << "ID: " << elem.id << endl;
  for (auto e: elem.nodes) {
    os << e << " ";
  }
  os << endl;
  return os;
}

/* Calculate K and F matrices for each element */
void Element::computeXY() {

  /* For each of x, dx/dzeta, dx/deta */
  for (int i = 0; i < 3; i++) {
    /* For each quadrature point */
    for (int j = 0; j < quadDegree; j++) {
      /* Sum over the ref. shape functions */
      XGP[i][j] = 0;
      YGP[i][j] = 0;
      for (int k = 0; k < 3; k++) {
        double x = nodesX[nodes[k]];
        double y = nodesY[nodes[k]];
        XGP[i][j] += x * DPHI[k][i][j];
        YGP[i][j] += y * DPHI[k][i][j];
      }
    }
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
void Element::assemble(vector<T> &t, MatrixXd &F) {
  int ii, jj;

  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      ii = ID[nodes[i]];
      jj = ID[nodes[j]];
      if (ii != -1 && jj != -1) {
        //K(ii,jj) += Ke(i,j);
        t.push_back(Triplet<double>(ii, jj, Ke(i,j)));
      }
    }
  }

  for (int i = 0; i < 3; i++) {
    ii = ID[nodes[i]];
    if (ii != -1) {
      F(ii) += Fe(i);
    }
  }
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
      elements.push_back(Element(elem, row));
      elem++;
    }
  } else {
    cout << "File doesn't exist" << endl;
  }

}
