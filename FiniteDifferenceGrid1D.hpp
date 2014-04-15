#ifndef FINITE_DIFFERENCE_GRID_1D_HPP
#define FINITE_DIFFERENCE_GRID_1D_HPP

#include "iostream"
#include "Eigen/Sparse"
#include <vector>
#include "Node.hpp"
#include "gtest/gtest_prod.h"
#include "FiniteDifferenceGrid.hpp"

using namespace Eigen;
using namespace std;

class FiniteDifferenceGrid1D : public FiniteDifferenceGrid {
    private:
        // test framework
        FRIEND_TEST(FiniteDifferenceGrid, mesh_formation);
        friend class BvpOde1D; 
        friend class BvpOde; 
    public:

        FiniteDifferenceGrid1D(int xNumNodes, double xMin, double xMax);
};

#endif // FINITE_DIFFERENCE_GRID_HPP