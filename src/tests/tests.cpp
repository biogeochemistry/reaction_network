#include "gtest/gtest.h"
#include "BvpOde.hpp"
#include "iostream"
#include "gnuplot_i.hpp"

double model_prob_1_rhs(double x){return 1.0;}
double model_prob_2_rhs(double x){return 34.0*sin(x);}

double bc_x0(double y){return -pow(y,2);}
double bc_xN(double y){return pow(y,2);}
double bc_y0(double x){return -pow(x,4);}
double bc_yN(double x){return pow(x,4);}

TEST(boundary_conditions, assigning_var){
    BoundaryConditions bc_dir, bc_neu;
    // Dirichlet
    bc_dir.SetX0DirichletBc2D(bc_x0);
    bc_dir.SetXNDirichletBc2D(bc_xN);
    bc_dir.SetY0DirichletBc2D(bc_y0);
    bc_dir.SetYNDirichletBc2D(bc_yN);

    // Neumann
    bc_neu.SetX0NeumannBc2D(bc_x0);
    bc_neu.SetXNNeumannBc2D(bc_xN);
    bc_neu.SetY0NeumannBc2D(bc_y0);
    bc_neu.SetYNNeumannBc2D(bc_yN);

    Vector3d bc_v(1,2,3);
    for (int i = 0; i < bc_v.size(); ++i) {
        /*
        Testing of the assigning of the value of given function to the BC
         */
        // Dirichlet
        EXPECT_EQ( bc_dir.mpX0BcFunc2D(bc_v[i]), bc_x0(bc_v[i]) );
        EXPECT_EQ( bc_dir.mpXNBcFunc2D(bc_v[i]), bc_xN(bc_v[i]) );
        EXPECT_EQ( bc_dir.mpY0BcFunc2D(bc_v[i]), bc_y0(bc_v[i]) );
        EXPECT_EQ( bc_dir.mpYNBcFunc2D(bc_v[i]), bc_yN(bc_v[i]) );
        EXPECT_EQ( bc_dir.mX0BcIsDirichlet,  true);
        EXPECT_EQ( bc_dir.mXNBcIsDirichlet,  true);
        EXPECT_EQ( bc_dir.mYNBcIsDirichlet,  true);
        EXPECT_EQ( bc_dir.mY0BcIsDirichlet,  true);
        EXPECT_EQ( bc_dir.mX0BcIsNeumann, false);
        EXPECT_EQ( bc_dir.mXNBcIsNeumann, false);
        EXPECT_EQ( bc_dir.mYNBcIsNeumann, false);
        EXPECT_EQ( bc_dir.mY0BcIsNeumann, false);

        // Neumann
        EXPECT_EQ( bc_neu.mpX0BcFunc2D(bc_v[i]), bc_x0(bc_v[i]) );
        EXPECT_EQ( bc_neu.mpXNBcFunc2D(bc_v[i]), bc_xN(bc_v[i]) );
        EXPECT_EQ( bc_neu.mpY0BcFunc2D(bc_v[i]), bc_y0(bc_v[i]) );
        EXPECT_EQ( bc_neu.mpYNBcFunc2D(bc_v[i]), bc_yN(bc_v[i]) );
        EXPECT_EQ( bc_neu.mX0BcIsDirichlet,  false);
        EXPECT_EQ( bc_neu.mXNBcIsDirichlet,  false);
        EXPECT_EQ( bc_neu.mYNBcIsDirichlet,  false);
        EXPECT_EQ( bc_neu.mY0BcIsDirichlet,  false);
        EXPECT_EQ( bc_neu.mX0BcIsNeumann, true);
        EXPECT_EQ( bc_neu.mXNBcIsNeumann, true);
        EXPECT_EQ( bc_neu.mYNBcIsNeumann, true);
        EXPECT_EQ( bc_neu.mY0BcIsNeumann, true);
    }

}

TEST(FiniteDifferenceGrid, mesh_formation) {
    FiniteDifferenceGrid grid = FiniteDifferenceGrid(3,1,2);
    Vector3d v(1,1.5,2);  
    for (int i = 0; i < v.size(); ++i) {
        EXPECT_EQ(v[i], grid.mNodes[i].C.x) << "Vectors x and y differ at index ";
    } 
}

TEST(boundary_conditions, constructor) {
    BoundaryConditions bc;
    ASSERT_EQ(bc.mX0BcIsDirichlet, false);
    ASSERT_EQ(bc.mXNBcIsDirichlet, false);
    ASSERT_EQ(bc.mX0BcIsNeumann, false);
    ASSERT_EQ(bc.mXNBcIsNeumann, false);
}

TEST(second_order_ode, assigning_var){
    SecondOrderOde ode_mp1(-1.0, 0.0, 0.0, model_prob_1_rhs, 0.0, 1.0);
    ASSERT_EQ(ode_mp1.mCoeffOfUxx, -1);
    ASSERT_EQ(ode_mp1.mCoeffOfUx, 0);
    ASSERT_EQ(ode_mp1.mCoeffOfU, 0);
    ASSERT_EQ(ode_mp1.mXmin, 0);
    ASSERT_EQ(ode_mp1.mXmax, 1);
}

TEST(bvpode, error_of_the_solution){
    Gnuplot g1;
    SecondOrderOde ode_mp1(-1.0, 0.0, 0.0, model_prob_1_rhs, 0.0, 1);
    BoundaryConditions bc_mp1;
    bc_mp1.SetX0DirichletBc1D(0);
    bc_mp1.SetXNDirichletBc1D(0);
    BvpOde bvpode_mp1(&ode_mp1, &bc_mp1, 128);
    bvpode_mp1.SetFilename("model_problem_results1.dat");
    bvpode_mp1.Solve();
    // g1.set_style("points").plot_xy(bvpode_mp1.mpGrid->xGrid,bvpode_mp1.mpSolVec,"differentiation");
    // g1.set_style("lines").plot_equation("0.5*x*(1-x)","exact solution");

    
    Gnuplot g2;
    SecondOrderOde ode_mp2(1.0, 3.0, -4.0, model_prob_2_rhs, 0.0, M_PI);
    BoundaryConditions bc_mp2;
    bc_mp2.SetX0NeumannBc1D(-5.0);
    bc_mp2.SetXNDirichletBc1D(4.0);
    BvpOde bvpode_mp2(&ode_mp2, &bc_mp2, 128);
    bvpode_mp2.SetFilename("model_problem_results2.dat");
    bvpode_mp2.Solve();
    // g2.set_style("points").plot_xy(bvpode_mp2.mpGrid->xGrid,bvpode_mp2.mpSolVec);
    // g2.set_style("lines").plot_equation("(4*exp(x) + exp(-4*x) ) / (4*exp(pi)+exp(-4*pi))-5*sin(x)-3*cos(x)","exact solution");
}

TEST(FiniteDifferenceGrid2d, mesh_formation) {
    FiniteDifferenceGrid grid = FiniteDifferenceGrid(3, 1, 2, 3, 1, 2);
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
    //    |    |            |           |
    //    |    0[1,1] ---1[1.5,1]------2[2,1]
    //    |
    //    0----------------------------------->
    //                                        x
    
    // Node # 0
    ASSERT_EQ(grid.mNodes[0].num, 0);
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
    ASSERT_EQ(grid.mNodes[1].num, 1);
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
    ASSERT_EQ(grid.mNodes[2].num, 2);
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
    ASSERT_EQ(grid.mNodes[3].num, 3);
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
    ASSERT_EQ(grid.mNodes[4].num, 4);
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
    ASSERT_EQ(grid.mNodes[5].num, 5);
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
    ASSERT_EQ(grid.mNodes[6].num, 6);
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
    ASSERT_EQ(grid.mNodes[7].num, 7);
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
    ASSERT_EQ(grid.mNodes[8].num, 8);
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

TEST(boundary_conditions2d, constructor) {
    BoundaryConditions bc;
    ASSERT_EQ(bc.mX0BcIsDirichlet, false);
    ASSERT_EQ(bc.mXNBcIsDirichlet, false);
    ASSERT_EQ(bc.mX0BcIsNeumann, false);
    ASSERT_EQ(bc.mXNBcIsNeumann, false);
    ASSERT_EQ(bc.mYNBcIsDirichlet, false);
    ASSERT_EQ(bc.mYNBcIsDirichlet, false);
    ASSERT_EQ(bc.mYNBcIsNeumann, false);
    ASSERT_EQ(bc.mY0BcIsNeumann, false);
}
