#include <cassert>
#include "FiniteDifferenceGrid.hpp"

using namespace Eigen;

void FiniteDifferenceGrid::xGridForamtion(int xNumNodes, double xMin, double xMax){
    double stepsize = (xMax-xMin)/(xNumNodes-1);
    xGrid.resize(xNumNodes);
    for(unsigned i = 0; i < xNumNodes; ++i) {
        xGrid(i) =  xMin+i*stepsize;
    }
}

void FiniteDifferenceGrid::yGridFormation(int yNumNodes, double yMin, double yMax){
    double stepsize = (yMax-yMin)/(yNumNodes-1);
    yGrid.resize(yNumNodes);
    for(unsigned i = 0; i < yNumNodes; ++i) {
        yGrid(i) =  yMin+i*stepsize;
    }
}
