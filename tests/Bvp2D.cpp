#include "FiniteDifferenceGrid2D.hpp"
#include "gtest/gtest.h"

using namespace std;
// using namespace Eigen;

double model_prob_1_rhs(double x){return 1.0;}
double model_prob_2_rhs(double x){return 34.0*sin(x);}

TEST(FiniteDifferenceGrid2D, mesh_formation) {
    FiniteDifferenceGrid2D grid = FiniteDifferenceGrid2D(3, 1, 2, 3, 1, 2);
    Vector3d v(1,1.5,2);  
    for (int i = 0; i < v.size(); ++i) {
        ASSERT_EQ(v[i], grid.xGrid[i]) << "xGrid formation";
        ASSERT_EQ(v[i], grid.yGrid[i]) << "yGrid formation";
    } 
    // testing mesh 
    // format: Node#[x,y]
    // 
    //   y^ 
    //    |   6[1,2]-----7[1.5,2]------8[2,2]
    //    |    |           |            |
    //    |    |           |            |
    //    |   3[1,1.5]--4[1.5,1.5]----5[2,1.5]
    //    |    |            |           |
    //    |    |            |            |
    //    |    0[1,1] ---1[1.5,1]------2[2,1]
    //    |
    //    0----------------------------------->
    //                                        x
    
    // Node # 0
    ASSERT_EQ(grid.mNodes[0].C.x, 1.0);
    ASSERT_EQ(grid.mNodes[0].C.y, 1.0);
    ASSERT_EQ(grid.mNodes[0].N.x, 1.0);
    ASSERT_EQ(grid.mNodes[0].N.y, 1.5);
    ASSERT_EQ(grid.mNodes[0].E.x, 1.5);
    ASSERT_EQ(grid.mNodes[0].E.y, 1);
    ASSERT_EQ(grid.mNodes[0].S.x, NULL);
    ASSERT_EQ(grid.mNodes[0].S.y, NULL);
    ASSERT_EQ(grid.mNodes[0].W.x, NULL);
    ASSERT_EQ(grid.mNodes[0].W.y, NULL);

    // Node # 1
    ASSERT_EQ(grid.mNodes[1].C.x, 1.5);
    ASSERT_EQ(grid.mNodes[1].C.y, 1.0);
    ASSERT_EQ(grid.mNodes[1].N.x, 1.5);
    ASSERT_EQ(grid.mNodes[1].N.y, 1.5);
    ASSERT_EQ(grid.mNodes[1].E.x, 2);
    ASSERT_EQ(grid.mNodes[1].E.y, 1);
    ASSERT_EQ(grid.mNodes[1].S.x, NULL);
    ASSERT_EQ(grid.mNodes[1].S.y, NULL);
    ASSERT_EQ(grid.mNodes[1].W.x, 1);
    ASSERT_EQ(grid.mNodes[1].W.y, 1);

    // Node # 2
    ASSERT_EQ(grid.mNodes[2].C.x, 2);
    ASSERT_EQ(grid.mNodes[2].C.y, 1);
    ASSERT_EQ(grid.mNodes[2].N.x, 2);
    ASSERT_EQ(grid.mNodes[2].N.y, 1.5);
    ASSERT_EQ(grid.mNodes[2].E.x, NULL);
    ASSERT_EQ(grid.mNodes[2].E.y, NULL);
    ASSERT_EQ(grid.mNodes[2].S.x, NULL);
    ASSERT_EQ(grid.mNodes[2].S.y, NULL);
    ASSERT_EQ(grid.mNodes[2].W.x, 1.5);
    ASSERT_EQ(grid.mNodes[2].W.y, 1);

    // Node # 3
    ASSERT_EQ(grid.mNodes[3].C.x, 1);
    ASSERT_EQ(grid.mNodes[3].C.y, 1.5);
    ASSERT_EQ(grid.mNodes[3].N.x, 1);
    ASSERT_EQ(grid.mNodes[3].N.y, 2);
    ASSERT_EQ(grid.mNodes[3].E.x, 1.5);
    ASSERT_EQ(grid.mNodes[3].E.y, 1.5);
    ASSERT_EQ(grid.mNodes[3].S.x, 1);
    ASSERT_EQ(grid.mNodes[3].S.y, 1);
    ASSERT_EQ(grid.mNodes[3].W.x, NULL);
    ASSERT_EQ(grid.mNodes[3].W.y, NULL);

    // Node # 4
    ASSERT_EQ(grid.mNodes[4].C.x, 1.5);
    ASSERT_EQ(grid.mNodes[4].C.y, 1.5);
    ASSERT_EQ(grid.mNodes[4].N.x, 1.5);
    ASSERT_EQ(grid.mNodes[4].N.y, 2);
    ASSERT_EQ(grid.mNodes[4].E.x, 2);
    ASSERT_EQ(grid.mNodes[4].E.y, 1.5);
    ASSERT_EQ(grid.mNodes[4].S.x, 1.5);
    ASSERT_EQ(grid.mNodes[4].S.y, 1);
    ASSERT_EQ(grid.mNodes[4].W.x, 1);
    ASSERT_EQ(grid.mNodes[4].W.y, 1.5);

    // Node # 5
    ASSERT_EQ(grid.mNodes[5].C.x, 2);
    ASSERT_EQ(grid.mNodes[5].C.y, 1.5);
    ASSERT_EQ(grid.mNodes[5].N.x, 2);
    ASSERT_EQ(grid.mNodes[5].N.y, 2);
    ASSERT_EQ(grid.mNodes[5].E.x, NULL);
    ASSERT_EQ(grid.mNodes[5].E.y, NULL);
    ASSERT_EQ(grid.mNodes[5].S.x, 2);
    ASSERT_EQ(grid.mNodes[5].S.y, 1);
    ASSERT_EQ(grid.mNodes[5].W.x, 1.5);
    ASSERT_EQ(grid.mNodes[5].W.y, 1.5);

    // Node # 6
    ASSERT_EQ(grid.mNodes[6].C.x, 1);
    ASSERT_EQ(grid.mNodes[6].C.y, 2);
    ASSERT_EQ(grid.mNodes[6].N.x, NULL);
    ASSERT_EQ(grid.mNodes[6].N.y, NULL);
    ASSERT_EQ(grid.mNodes[6].E.x, 1.5);
    ASSERT_EQ(grid.mNodes[6].E.y, 2);
    ASSERT_EQ(grid.mNodes[6].S.x, 1);
    ASSERT_EQ(grid.mNodes[6].S.y, 1.5);
    ASSERT_EQ(grid.mNodes[6].W.x, NULL);
    ASSERT_EQ(grid.mNodes[6].W.y, NULL);

    // Node # 7
    ASSERT_EQ(grid.mNodes[7].C.x, 1.5);
    ASSERT_EQ(grid.mNodes[7].C.y, 2);
    ASSERT_EQ(grid.mNodes[7].N.x, NULL);
    ASSERT_EQ(grid.mNodes[7].N.y, NULL);
    ASSERT_EQ(grid.mNodes[7].E.x, 2);
    ASSERT_EQ(grid.mNodes[7].E.y, 2);
    ASSERT_EQ(grid.mNodes[7].S.x, 1.5);
    ASSERT_EQ(grid.mNodes[7].S.y, 1.5);
    ASSERT_EQ(grid.mNodes[7].W.x, 1);
    ASSERT_EQ(grid.mNodes[7].W.y, 2);

    // Node # 8
    ASSERT_EQ(grid.mNodes[8].C.x, 2);
    ASSERT_EQ(grid.mNodes[8].C.y, 2);
    ASSERT_EQ(grid.mNodes[8].N.x, NULL);
    ASSERT_EQ(grid.mNodes[8].N.y, NULL);
    ASSERT_EQ(grid.mNodes[8].E.x, NULL);
    ASSERT_EQ(grid.mNodes[8].E.y, NULL);
    ASSERT_EQ(grid.mNodes[8].S.x, 2);
    ASSERT_EQ(grid.mNodes[8].S.y, 1.5);
    ASSERT_EQ(grid.mNodes[8].W.x, 1.5);
    ASSERT_EQ(grid.mNodes[8].W.y, 2);
}