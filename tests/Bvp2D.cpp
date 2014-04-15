#include "FiniteDifferenceGrid2D.hpp"
#include "gtest/gtest.h"

using namespace std;
// using namespace Eigen;

double model_prob_1_rhs(double x){return 1.0;}
double model_prob_2_rhs(double x){return 34.0*sin(x);}

TEST(FiniteDifferenceGrid2D, mesh_formation) {
    FiniteDifferenceGrid2D grid = FiniteDifferenceGrid2D(5, 5, 1, 2, 2, 3);
    cout << grid.mNodes[0].C.x;
  // Vector3d v(1,1.5,2);  
  // for (int i = 0; i < v.size(); ++i)
  // {
    // EXPECT_EQ(v[i], grid.mNodes[i].coordinate) << "Vectors x and y differ at index ";
  // } 
}

