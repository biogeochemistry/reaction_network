#include <cassert>
#include "FiniteDifferenceGrid.hpp"
#include "Node.hpp"

FiniteDifferenceGrid::FiniteDifferenceGrid(int numNodes, double xMin, double xMax) {

double stepsize = (xMax-xMin)/((double)(numNodes-1));
    xGrid.resize(numNodes);
    for (int i=0; i<numNodes; i++){
        Node node;
        node.coordinate = xMin+i*stepsize;
        mNodes.push_back(node);
        xGrid(i) = node.coordinate;
        }
    assert(mNodes.size() == numNodes);
}