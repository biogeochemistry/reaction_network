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
    private:
        // test framework
        FRIEND_TEST(FiniteDifferenceGrid, mesh_formation);

    public:
        std::vector<Node> mNodes;
        VectorXd xGrid;
        friend class BvpOde;
        void xGridForamtion(int xNumNodes, double xMin, double xMax);
        FiniteDifferenceGrid(int xNumNodes, double xMin, double xMax);
};

#endif // FINITE_DIFFERENCE_GRID_HPP