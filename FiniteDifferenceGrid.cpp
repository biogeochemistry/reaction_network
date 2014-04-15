#include <cassert>
#include "FiniteDifferenceGrid.hpp"
#include "Node.hpp"

FiniteDifferenceGrid::FiniteDifferenceGrid(int numNodes, double xMin, double xMax) {
    xGridForamtion(numNodes, xMin, xMax);
    for (int i=0; i<numNodes; i++){
        Node node;
        node.coordinate = xGrid(i);
        mNodes.push_back(node);
        }
    assert(mNodes.size() == numNodes);
}

void FiniteDifferenceGrid::xGridForamtion(int numNodes, double xMin, double xMax){
    double stepsize = (xMax-xMin)/((double)(numNodes-1));
    xGrid.resize(numNodes);
    for(unsigned i = 0; i < numNodes; ++i) {
        xGrid(i) =  xMin+i*stepsize;
    }
}