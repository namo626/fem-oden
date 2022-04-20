#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <Eigen/Dense>
#include <string>
#include <sstream>
#include <fstream>
#include <iostream>

#include "Global.h"
#include "Post.h"

using namespace std;
using namespace Eigen;

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
