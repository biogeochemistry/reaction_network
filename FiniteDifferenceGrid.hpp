#ifndef FINITE_DIFFERENCE_GRID_HPP
#define FINITE_DIFFERENCE_GRID_HPP

#include "iostream"
#include "Eigen/Sparse"
#include <vector>
#include "Node.hpp"
#include "gtest/gtest_prod.h"

using namespace Eigen;
using namespace std;

class FiniteDifferenceGrid {

    public:
        // 1D constructor
        FiniteDifferenceGrid(int xNumNodes, double xMin, double xMax);
        // 2D constructor
        FiniteDifferenceGrid(int xNumNodes, double xMin, double xMax, int yNumNodes, double yMin, double yMax);
        
        vector<Node> mNodes;
        VectorXd xGrid;
        VectorXd yGrid;
        friend class BvpOde; 
        void yGridFormation(int yNumNodes, double yMin, double yMax);
        void xGridForamtion(int xNumNodes, double xMin, double xMax);
};

#endif // FINITE_DIFFERENCE_GRID_HPP