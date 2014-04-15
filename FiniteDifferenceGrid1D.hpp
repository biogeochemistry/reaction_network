#ifndef FINITE_DIFFERENCE_GRID_1D_HPP
#define FINITE_DIFFERENCE_GRID_1D_HPP

#include "iostream"
#include "Eigen/Sparse"
#include <vector>
#include "Node1D.hpp"
#include "gtest/gtest_prod.h"
#include "FiniteDifferenceGrid.hpp"

using namespace Eigen;
using namespace std;

class FiniteDifferenceGrid1D : public FiniteDifferenceGrid {
    private:
        // test framework
        FRIEND_TEST(FiniteDifferenceGrid, mesh_formation);
        vector<Node1D> mNodes;
        friend class BvpOde1D; 
    public:

        FiniteDifferenceGrid1D(int xNumNodes, double xMin, double xMax);
};

#endif // FINITE_DIFFERENCE_GRID_HPP