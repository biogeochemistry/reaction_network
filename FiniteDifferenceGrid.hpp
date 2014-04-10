#ifndef FINITE_DIFFERENCE_GRID_HPP
#define FINITE_DIFFERENCE_GRID_HPP

#include "iostream"
#include "Eigen/Dense"
#include "Node.hpp"
#include "gtest/gtest_prod.h"
#include <vector>

using namespace Eigen;
using namespace std;

class FiniteDifferenceGrid {
  
private:
  FRIEND_TEST(FiniteDifferenceGrid, mesh_formation);
  std::vector<Node> VecNod;
  VectorXd mNodes;
public:
  FiniteDifferenceGrid(int numNodes, double xMin, double xMax){
    mNodes.resize(numNodes);
    for (int i = 0; i < numNodes; ++i)
    {
      mNodes[i] = xMin + (xMax - xMin) / (numNodes-1) * i;
    };
  };
};

#endif // FINITE_DIFFERENCE_GRID_HPP