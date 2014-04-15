#include <cassert>
#include "FiniteDifferenceGrid.hpp"
#include "Node.hpp"

using namespace Eigen;

void FiniteDifferenceGrid::xGridForamtion(int xNumNodes, double xMin, double xMax){
    double stepsize = (xMax-xMin)/(xNumNodes-1);
    xGrid.resize(xNumNodes);
    for(unsigned i = 0; i < xNumNodes; ++i) {
        xGrid(i) =  xMin+i*stepsize;
    }
}