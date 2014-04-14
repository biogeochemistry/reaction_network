#include "iostream"
#include <cmath>
#include <string>
#include "BvpOde.hpp"
#include "gnuplot_i.hpp"


double model_prob_1_rhs(double x){return 1.0;}
double model_prob_2_rhs(double x){return 34.0*sin(x);}
double model_prob_3_rhs(double x){return 0;}
int main(int argc, char* argv[]) {
    Gnuplot g1;

    SecondOrderOde ode_mp1(-1.0, 0.0, 0.0, model_prob_1_rhs, 0.0, 1);
    BoundaryConditions bc_mp1;
    bc_mp1.SetLhsDirichletBc(0);
    bc_mp1.SetRhsDirichletBc(0);
    BvpOde bvpode_mp1(&ode_mp1, &bc_mp1, 101);
    bvpode_mp1.SetFilename("model_problem_results1.dat");
    bvpode_mp1.Solve();
    g1.set_style("points").plot_xy(bvpode_mp1.mpGrid->xGrid,bvpode_mp1.mpSolVec);

    
    SecondOrderOde ode_mp2(1.0, 3.0, -4.0, model_prob_2_rhs, 0.0, M_PI);
    BoundaryConditions bc_mp2;
    bc_mp2.SetLhsNeumannBc(-5.0);
    bc_mp2.SetRhsDirichletBc(4.0);
    BvpOde bvpode_mp2(&ode_mp2, &bc_mp2, 101);
    bvpode_mp2.SetFilename("model_problem_results2.dat");
    bvpode_mp2.Solve();
    g1.set_style("points").plot_xy(bvpode_mp2.mpGrid->xGrid,bvpode_mp2.mpSolVec);

    SecondOrderOde ode_mp3(1.0, 6.0, 1.0, model_prob_3_rhs, 0.0, 1.0);
    BoundaryConditions bc_mp3;
    bc_mp3.SetLhsDirichletBc(10.0);
    bc_mp3.SetRhsDirichletBc(-4.0);
    BvpOde bvpode_mp3(&ode_mp3, &bc_mp3, 101);
    bvpode_mp3.SetFilename("model_problem_results3.dat");
    bvpode_mp3.Solve();
    g1.set_style("points").plot_xy(bvpode_mp3.mpGrid->xGrid,bvpode_mp3.mpSolVec);


   
    

    return 0;
}