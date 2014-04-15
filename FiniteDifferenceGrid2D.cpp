#include <cassert>
#include "FiniteDifferenceGrid2D.hpp"

using namespace Eigen;

FiniteDifferenceGrid2D::FiniteDifferenceGrid2D(int xNumNodes, int yNumNodes, double xMin, double xMax, double yMin, double yMax) {
    xGridForamtion(xNumNodes, xMin, xMax);
    yGridFormation(yNumNodes, yMin, yMax);
    for (int j=0; j < yNumNodes; ++j) {
        for (int i = 0; i < xNumNodes ; ++i) {
            Node2D node;

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

            mNodes.push_back(node);
        }
    }
    // assert(mNodes.size() == xNumNodes);
}
