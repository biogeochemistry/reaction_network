#include "gtest/gtest.h"
#include "BvpOde1D.hpp"
#include "iostream"
#include "gnuplot_i.hpp"


// using namespace std;
// using namespace Eigen;

double model_prob_1_rhs(double x){return 1.0;}
double model_prob_2_rhs(double x){return 34.0*sin(x);}

TEST(FiniteDifferenceGrid, mesh_formation) {
  FiniteDifferenceGrid1D grid = FiniteDifferenceGrid1D(3,1,2);
  Vector3d v(1,1.5,2);  
  for (int i = 0; i < v.size(); ++i)
  {
    EXPECT_EQ(v[i], grid.mNodes[i].coordinate) << "Vectors x and y differ at index ";
  } 
}

TEST(boundary_conditions, constractor) {
  BoundaryConditions bc;
  ASSERT_EQ(bc.mLhsBcIsDirichlet, false);
  ASSERT_EQ(bc.mRhsBcIsDirichlet, false);
  ASSERT_EQ(bc.mLhsBcIsNeumann, false);
  ASSERT_EQ(bc.mRhsBcIsNeumann, false);
}

TEST(second_order_ode, assigning_var){
  SecondOrderOde1D ode_mp1(-1.0, 0.0, 0.0, model_prob_1_rhs, 0.0, 1.0);
  ASSERT_EQ(ode_mp1.mCoeffOfUxx, -1);
  ASSERT_EQ(ode_mp1.mCoeffOfUx, 0);
  ASSERT_EQ(ode_mp1.mCoeffOfU, 0);
  ASSERT_EQ(ode_mp1.mXmin, 0);
  ASSERT_EQ(ode_mp1.mXmax, 1);
}

TEST(bvpode, error_of_the_solution){
    Gnuplot g1;
    SecondOrderOde1D ode_mp1(-1.0, 0.0, 0.0, model_prob_1_rhs, 0.0, 1);
    BoundaryConditions bc_mp1;
    bc_mp1.SetLhsDirichletBc(0);
    bc_mp1.SetRhsDirichletBc(0);
    BvpOde1D bvpode_mp1(&ode_mp1, &bc_mp1, 128);
    bvpode_mp1.SetFilename("model_problem_results1.dat");
    bvpode_mp1.Solve();
    g1.set_style("points").plot_xy(bvpode_mp1.mpGrid->xGrid,bvpode_mp1.mpSolVec,"differentiation");
    g1.set_style("lines").plot_equation("0.5*x*(1-x)","exact solution");

    
    Gnuplot g2;
    SecondOrderOde1D ode_mp2(1.0, 3.0, -4.0, model_prob_2_rhs, 0.0, M_PI);
    BoundaryConditions bc_mp2;
    bc_mp2.SetLhsNeumannBc(-5.0);
    bc_mp2.SetRhsDirichletBc(4.0);
    BvpOde1D bvpode_mp2(&ode_mp2, &bc_mp2, 128);
    bvpode_mp2.SetFilename("model_problem_results2.dat");
    bvpode_mp2.Solve();
    g2.set_style("points").plot_xy(bvpode_mp2.mpGrid->xGrid,bvpode_mp2.mpSolVec);
    g2.set_style("lines").plot_equation("(4*exp(x) + exp(-4*x) ) / (4*exp(pi)+exp(-4*pi))-5*sin(x)-3*cos(x)","exact solution");
}
