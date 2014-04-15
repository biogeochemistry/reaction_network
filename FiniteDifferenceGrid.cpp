#include <cassert>
#include "FiniteDifferenceGrid.hpp"
#include "Node.hpp"

FiniteDifferenceGrid::FiniteDifferenceGrid(int xNumNodes, double xMin, double xMax) {
    xGridForamtion(xNumNodes, xMin, xMax);
    for (int i=0; i<xNumNodes; i++){
        Node node;
        node.coordinate = xGrid(i);
        mNodes.push_back(node);
        }
    assert(mNodes.size() == xNumNodes);
}

void FiniteDifferenceGrid::xGridForamtion(int xNumNodes, double xMin, double xMax){
    double stepsize = (xMax-xMin)/((double)(xNumNodes-1));
    xGrid.resize(xNumNodes);
    for(unsigned i = 0; i < xNumNodes; ++i) {
        xGrid(i) =  xMin+i*stepsize;
    }
}