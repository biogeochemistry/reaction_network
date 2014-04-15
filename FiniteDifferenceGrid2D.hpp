#ifndef FINITE_DIFFERENCE_GRID_2D_HPP
#define FINITE_DIFFERENCE_GRID_2D_HPP

#include "iostream"
#include "Eigen/Sparse"
#include <vector>
#include "Node2D.hpp"
#include "gtest/gtest_prod.h"
#include "FiniteDifferenceGrid.hpp"

using namespace Eigen;
using namespace std;

class FiniteDifferenceGrid2D : public FiniteDifferenceGrid {
    private:
        // test framework
        FRIEND_TEST(FiniteDifferenceGrid2D, mesh_formation);
        vector<Node2D> mNodes;
        friend class BvpOde1D; 
        friend class BvpOde; 
    public:

        FiniteDifferenceGrid2D(int xNumNodes, int yNumNodes, double xMin, double xMax, double yMin, double yMax);
};

#endif // FINITE_DIFFERENCE_GRID_HPP