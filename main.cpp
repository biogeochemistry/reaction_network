#include "iostream"
#include <cmath>
#include <string>
#include "BvpOde1D.hpp"
#include "gnuplot_i.hpp"


double model_prob_1_rhs(double x){return 1.0;}
double model_prob_2_rhs(double x){return 34.0*sin(x);}
double model_prob_3_rhs(double x){return 0;}
int main(int argc, char* argv[]) {
    Gnuplot g1;

    SecondOrderOde1D ode_mp1(-1.0, 0.0, 0.0, model_prob_1_rhs, 0.0, 1.0);
    BoundaryConditions1D bc_mp1;
    bc_mp1.SetLhsDirichletBc1D(0);
    bc_mp1.SetRhsDirichletBc1D(0);
    BvpOde1D bvpode_mp1(&ode_mp1, &bc_mp1, 1024);
    bvpode_mp1.SetFilename("model_problem_results1.dat");
    bvpode_mp1.Solve();
    g1.set_style("points").plot_xy(bvpode_mp1.mpGrid->xGrid,bvpode_mp1.mpSolVec,"differentiation");
    g1.set_style("lines").plot_equation("0.5*x*(1-x)","exact solution");

    
    // Gnuplot g2;
    // SecondOrderOde1D ode_mp2(1.0, 3.0, -4.0, model_prob_2_rhs, 0.0, M_PI);
    // BoundaryConditions bc_mp2;
    // bc_mp2.SetLhsNeumannBc(-5.0);
    // bc_mp2.SetRhsDirichletBc(4.0);
    // BvpOde1D bvpode_mp2(&ode_mp2, &bc_mp2, 1024);
    // bvpode_mp2.SetFilename("model_problem_results2.dat");
    // bvpode_mp2.Solve();
    // g2.set_style("points").plot_xy(bvpode_mp2.mpGrid->xGrid,bvpode_mp2.mpSolVec);
    // g2.set_style("lines").plot_equation("(4*exp(x) + exp(-4*x) ) / (4*exp(pi)+exp(-4*pi))-5*sin(x)-3*cos(x)","exact solution");

    return 0;
}