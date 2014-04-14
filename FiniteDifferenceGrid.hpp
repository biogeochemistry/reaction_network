#ifndef FINITE_DIFFERENCE_GRID_HPP
#define FINITE_DIFFERENCE_GRID_HPP

#include "iostream"
#include "Eigen/Sparse"
#include <vector>
#include "Node.hpp"

using namespace Eigen;
using namespace std;

class FiniteDifferenceGrid {
    private:

    public:
        std::vector<Node> mNodes;
        VectorXd xGrid;
        friend class BvpOde;
        FiniteDifferenceGrid(int numNodes, double xMin, double xMax);
};

#endif // FINITE_DIFFERENCE_GRID_HPP