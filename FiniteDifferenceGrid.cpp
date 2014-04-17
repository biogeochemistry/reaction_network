#include <cassert>
#include "FiniteDifferenceGrid.hpp"

using namespace Eigen;

void FiniteDifferenceGrid::xGridForamtion(int xNumNodes, double xMin, double xMax){
    assert(xMin < xMax);
    double stepsize = (xMax-xMin)/(xNumNodes-1);
    xGrid.resize(xNumNodes);
    for(int i = 0; i < xNumNodes; ++i) {
        xGrid(i) =  xMin+i*stepsize;
    }
}

void FiniteDifferenceGrid::yGridFormation(int yNumNodes, double yMin, double yMax){
    assert(yMin < yMax);
    double stepsize = (yMax-yMin)/(yNumNodes-1);
    yGrid.resize(yNumNodes);
    for(int i = 0; i < yNumNodes; ++i) {
        yGrid(i) =  yMin+i*stepsize;
    }
}
