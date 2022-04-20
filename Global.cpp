#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include "Global.h"
#include <math.h>
#include <algorithm>

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
  return 1 + pow(x,2) + 2*pow(y,2);
}

/* Load function */
double f(double x, double y) {
  return -6.0;
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
