#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

using namespace std;
using namespace Eigen;
typedef Triplet<double> T;

double* Q; // values at the essential BC, zero otherwise

/* Gaussian integration points and weights on the ref. triangle */
double GaussPoints[3][3] = {
  {0.5, 0},
  {0, 0.5},
  {0.5, 0.5}
};

double GaussWeights[3] = {1./3, 1./3, 1./3};

class Element {
public:
  Element(int id, vector<int> nodes);
  // Id
  int id;
  // List of nodes in counterclockwise direction
  vector<int> nodes;

  void computeXY();
  void computeJ();
  void computeZE();
  void computePhi();
  void initQ();
  void precompute();

  void computeKF();
  void assemble(vector<T> &t, MatrixXd &F);

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

  MatrixXd Ke;
  MatrixXd Fe;
  Vector3d Qe;

};

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

int nodesCount;
vector<double> nodesX;
vector<double> nodesY;
int quadDegree;
double DPHI[3][3][3];
vector<int> dirichletNodes;  // list of nodes with specified essential BC
int* ID; // converts global node number to equation number in global matrix
vector<Element> elements;
int neq;  // number of global equations = nodesCount - dirichlet nodes

/* Test solution */
double ue(double x, double y) {
  return 1 + pow(x,2) + 2*pow(y,2);
}

/* Load function */
double f(double x, double y) {
  return -6.0;
}

/* Initialize the vector of elements */
void read_elems() {
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

vector<double> assignNodes(MatrixXd &sol) {
  vector<double> vals(nodesCount, 0);
  for (int i = 0; i < nodesCount; i++) {
    if (ID[i] == -1) {
      vals[i] = Q[i];
    } else {
      vals[i] = sol(ID[i]);
    }
  }
  return vals;
}

void write_solution(vector<double> &vals) {
  ofstream outfile("output.csv");
  cout << "Writing solution to output.csv" << endl;

  for (int i = 0; i < nodesCount; i++) {
    outfile << nodesX[i] << "," << nodesY[i] << "," << vals[i] << "\n";
  }
}


int main() {
  read_elems();
  // for (auto e: elements) {
  //   cout << e << endl;
  // }

  read_nodes();
  init_dirichlet();
  init_ID();
  quadDegree = 3;
  init_dPhi();
  // for (int i = 0; i < 3; i++) {
  //   for (int j = 0; j < 3; j++) {
  //     printf("%.2f ", DPHI[i][j][2]);
  //   }
  //   printf("\n");
  // }
  /* Global matrices */
  MatrixXd F;
  SparseMatrix<double> K(neq, neq);
  vector<T> t;
  F.resize(neq, 1);
  F.setZero();


  cout << "Number of nodes: " << nodesCount << endl;
  cout << "Number of equations: " << neq << endl;

  // cout << "Dirichlet nodes: " << endl;
  // for (auto n: dirichletNodes) {
  //   cout << n << " ";
  // }
  // printf("\n");

  // cout << "ID array: " << endl;
  // for (int i = 0; i < nodesCount; i++) {
  //   printf("%d ", ID[i]);
  // }
  // printf("\n");


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
  //MatrixXd sol = K.fullPivLu().solve(F);

  vector<double> vals = assignNodes(sol);
  write_solution(vals);

  return 0;
}
