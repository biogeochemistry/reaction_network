#ifndef FINITE_DIFFERENCE_GRID_HPP
#define FINITE_DIFFERENCE_GRID_HPP

#include "iostream"
#include "Eigen/Sparse"
#include "gtest/gtest_prod.h"
#include <vector>

using namespace Eigen;
using namespace std;

class FiniteDifferenceGrid {
        friend class BvpOde;
    private:
        FRIEND_TEST(FiniteDifferenceGrid, mesh_formation);
        VectorXd mNodes;
    public:
        FiniteDifferenceGrid(int numNodes, double xMin, double xMax){
            mNodes.resize(numNodes);
            for (int i = 0; i < numNodes; ++i) {
                mNodes[i] = xMin + (xMax - xMin) / (numNodes-1) * i;
            };
        };
};

#endif // FINITE_DIFFERENCE_GRID_HPP