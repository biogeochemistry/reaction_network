#include <cassert>
#include "FiniteDifferenceGrid.hpp"

using namespace Eigen;

// 1D constructor
FiniteDifferenceGrid::FiniteDifferenceGrid(int xNumNodes, double xMin, double xMax) {
    xGridForamtion(xNumNodes, xMin, xMax);
    for (int i=0; i<xNumNodes; i++){
        Node node;
        node.C.x = xGrid(i);
        mNodes.push_back(node);
        }
    assert(mNodes.size() == xNumNodes);
}

// 2D constructor
FiniteDifferenceGrid::FiniteDifferenceGrid(int xNumNodes, double xMin, double xMax, int yNumNodes, double yMin, double yMax) {
    xGridForamtion(xNumNodes, xMin, xMax);
    yGridFormation(yNumNodes, yMin, yMax);
    for (int j=0; j < yNumNodes; ++j) {
        for (int i = 0; i < xNumNodes ; ++i) {
            Node node;

            if(j==0 && i!=0 && i!=xNumNodes-1) {
                // BC all except South
                // Central node
                node.C.x = xGrid(i);
                node.C.y = yGrid(j);

                // North node
                node.N.x = xGrid(i);
                node.N.y = yGrid(j+1);

                // East node
                node.E.x = xGrid(i+1);
                node.E.y = yGrid(j);

                // West node
                node.W.x = xGrid(i-1);
                node.W.y = yGrid(j);

                // South node
                node.S.x = NULL;
                node.S.y = NULL;
            }
            else if(j==yNumNodes-1 && i!=0 && i!=xNumNodes-1) {
                // BC all except North
                // Central node
                node.C.x = xGrid(i);
                node.C.y = yGrid(j);

                // South node 
                node.S.x = xGrid(i);
                node.S.y = yGrid(j-1);

                // East node
                node.E.x = xGrid(i+1);
                node.E.y = yGrid(j);

                // West node
                node.W.x = xGrid(i-1);
                node.W.y = yGrid(j);

               // North node
                node.N.x = NULL;
                node.N.y = NULL;
            }
            else if (i==0 && j!=0 && j!=yNumNodes-1)
            {
                // BC all except West
                // Central node
                node.C.x = xGrid(i);
                node.C.y = yGrid(j);

                // North node
                node.N.x = xGrid(i);
                node.N.y = yGrid(j+1);

                // South node 
                node.S.x = xGrid(i);
                node.S.y = yGrid(j-1);

                // East node
                node.E.x = xGrid(i+1);
                node.E.y = yGrid(j);               

               // West node
                node.W.x = NULL;
                node.W.y = NULL;
            }
            else if (i==xNumNodes-1 && j!=0 && j!=yNumNodes-1)
            {
                // BC all except East
                // Central node
                node.C.x = xGrid(i);
                node.C.y = yGrid(j);

                // North node
                node.N.x = xGrid(i);
                node.N.y = yGrid(j+1);

                // South node 
                node.S.x = xGrid(i);
                node.S.y = yGrid(j-1);

                // West node
                node.W.x = xGrid(i-1);
                node.W.y = yGrid(j);

                // East node
                node.E.x = NULL; 
                node.E.y = NULL;
            }
            else if (i==0 && j==0)
            {
                // Central node
                node.C.x = xGrid(i);
                node.C.y = yGrid(j);

                // North node
                node.N.x = xGrid(i);
                node.N.y = yGrid(j+1);

                // East node
                node.E.x = xGrid(i+1);
                node.E.y = yGrid(j);

                // West node
                node.W.x = NULL;
                node.W.y = NULL;

                // South node
                node.S.x = NULL;
                node.S.y = NULL;
            }
            else if (i==0 && j==yNumNodes-1)
            {
                // Central node
                node.C.x = xGrid(i);
                node.C.y = yGrid(j);

                // South node 
                node.S.x = xGrid(i);
                node.S.y = yGrid(j-1);

                // East node
                node.E.x = xGrid(i+1);
                node.E.y = yGrid(j);

                // West node
                node.W.x = NULL;
                node.W.y = NULL;

                // North node
                node.N.x = NULL;
                node.N.y = NULL;
            }
            else if (i==xNumNodes-1 && j==0)
            {
                // Central node
                node.C.x = xGrid(i);
                node.C.y = yGrid(j);

                // North node
                node.N.x = xGrid(i);
                node.N.y = yGrid(j+1);

                // West node
                node.W.x = xGrid(i-1);
                node.W.y = yGrid(j);

                // South node
                node.S.x = NULL; 
                node.S.y = NULL;

                // East node
                node.E.x = NULL; 
                node.E.y = NULL;
            }
            else if (i==xNumNodes-1 && j==yNumNodes-1)
            {
                // Central node
                node.C.x = xGrid(i);
                node.C.y = yGrid(j);

                // South node 
                node.S.x = xGrid(i);
                node.S.y = yGrid(j-1);

                // West node
                node.W.x = xGrid(i-1);
                node.W.y = yGrid(j);

                // East node
                node.E.x = NULL;
                node.E.y = NULL;

                // North node
                node.N.x = NULL;
                node.N.y = NULL;
            }
            else {
                // Formation of inner part of the stencil
                // Central node
                node.C.x = xGrid(i);
                node.C.y = yGrid(j);

                // North node
                node.N.x = xGrid(i);
                node.N.y = yGrid(j+1);

                // South node 
                node.S.x = xGrid(i);
                node.S.y = yGrid(j-1);

                // East node
                node.E.x = xGrid(i+1);
                node.E.y = yGrid(j);

                // West node
                node.W.x = xGrid(i-1);
                node.W.y = yGrid(j);
            }
            node.num = j*yNumNodes+i;
            mNodes.push_back(node);
        }
    }
    assert(mNodes.size() == yNumNodes*xNumNodes);
}

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
