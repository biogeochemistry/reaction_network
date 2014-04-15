#include <cassert>
#include "FiniteDifferenceGrid1D.hpp"
#include "Node1D.hpp"

using namespace Eigen;

FiniteDifferenceGrid1D::FiniteDifferenceGrid1D(int xNumNodes, double xMin, double xMax) {
    xGridForamtion(xNumNodes, xMin, xMax);
    for (int i=0; i<xNumNodes; i++){
        Node1D node;
        node.coordinate = xGrid(i);
        mNodes.push_back(node);
        }
    assert(mNodes.size() == xNumNodes);
}
