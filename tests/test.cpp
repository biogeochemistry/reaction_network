#include "gtest/gtest.h"
// #include "/BvpOde.hpp"
#include "Eigen/dense"
// #include "FiniteDifferenceGrid.hpp"
// #include "BoundaryConditions.hpp"
// #include "SecondOrderOde.hpp"
#include "BvpOde.hpp"
#include "iostream"


// using namespace std;
// using namespace Eigen;

double model_prob_1_rhs(double x){return 1.0;}

TEST(FiniteDifferenceGrid, mesh_formation) {
  FiniteDifferenceGrid grid = FiniteDifferenceGrid(3,1,2);
  Vector3d v(1,1.5,2);  
  for (int i = 0; i < v.size(); ++i)
  {
    EXPECT_EQ(v[i], grid.mNodes[i]) << "Vectors x and y differ at index ";
  } 
}

TEST(Boundary_Conditions, def_constractor) {
  BoundaryConditions bc;
  ASSERT_EQ(bc.mLhsBcIsDirichlet, false);
  ASSERT_EQ(bc.mRhsBcIsDirichlet, false);
  ASSERT_EQ(bc.mLhsBcIsNeumann, false);
  ASSERT_EQ(bc.mRhsBcIsNeumann, false);
}

TEST(second_order_ode, assigning_var){
  SecondOrderOde ode_mp1(-1.0, 0.0, 0.0, model_prob_1_rhs, 0.0, 1.0);
  ASSERT_EQ(ode_mp1.mCoeffOfUxx, -1);
  ASSERT_EQ(ode_mp1.mCoeffOfUx, 0);
  ASSERT_EQ(ode_mp1.mCoeffOfU, 0);
  ASSERT_EQ(ode_mp1.mXmin, 0);
  ASSERT_EQ(ode_mp1.mXmax, 1);
}
