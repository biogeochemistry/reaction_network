#include "iostream"
#include <cmath>
#include <string>
#include "CoupledPde.hpp"
#include "gnuplot_i.hpp"
#include <omp.h>
// #include <boost/test/unit_test.hpp>


using namespace Eigen;
using namespace std;


double model_prob_1_rhs(double x){return 2.0;}
double model_prob_2_rhs(double x){return 34.0*sin(x);}
double model_prob_3_rhs(double x, double y){return exp(-(x*x+y*y));}
double uj0(double x){return 100*sin(x);}



int main(int argc, char* argv[]) {
    // int nthreads, tid;
    // #pragma omp parallel
    // {
    //     /* Obtain thread number */
    //     tid = omp_get_thread_num();
    //     printf("Hello World from thread = %d\n", tid);

    //     /* Only master thread does this */
    //     if (tid == 0) {
    //         nthreads = omp_get_num_threads();
    //         printf("Number of threads = %d\n", nthreads);
    //     }
    // }
    
    

    SecondOrderOde ode_mp1(500, -5.0, 0.0, model_prob_1_rhs, 0.0, 15.0);
    SecondOrderOde ode_mp2(5, -5.0, 0.0, model_prob_1_rhs, 0.0, 15.0);
    BoundaryConditions bc_mp1;
    BoundaryConditions bc_mp2;
    bc_mp1.SetX0DirichletBc1D(0.15);
    bc_mp1.SetXNNeumannBc1D(0);
    bc_mp2.SetX0NeumannBc1D(1.0);
    bc_mp2.SetXNNeumannBc1D(0);
    // bc_mp1.SetX0NeumannBc1D(0);
    // bc_mp1.SetXNDirichletBc1D(0);
    bc_mp1.SetXNNeumannBc1D(0);
    // bc_mp1.SetXNRobinBc1D(0);
    BvpPde oxygen(&ode_mp1, &bc_mp1,0.01, 1.0, 10.0, 64, uj0);
    BvpPde organic_mater_1(&ode_mp2, &bc_mp2,0.01, 1.0, 10.0, 64, uj0);
    // BvpOde oxygen(&ode_mp1, &bc_mp1, 10);
    oxygen.SetFilename("oxygen_results.dat");
    organic_mater_1.SetFilename("om1_results.dat");
    // oxygen.Solve();
    CoupledPde cpde(&oxygen, &organic_mater_1);
    cpde.SolveSystem();


    Gnuplot g1, g2;
    g1.set_style("lines").plot_xy(oxygen.mpGrid->xGrid,oxygen.solution,"differentiation");
    g2.set_style("lines").plot_xy(organic_mater_1.mpGrid->xGrid,organic_mater_1.solution,"differentiation");
    // g1.set_style("lines").plot_equation("0.5*x*(1-x)","exact solution");

    // Gnuplot g2;
    // SecondOrderOde ode_mp2(10.0, 5.0, 0.0, uj0, 0.0, 20.0);
    // BoundaryConditions bc_mp2;
    // bc_mp2.SetX0NeumannBc1D(10.0);
    // bc_mp2.SetXNDirichletBc1D(0.0);
    // BvpOde organic_mater_1(&ode_mp2, &bc_mp2, 64);
    // organic_mater_1.SetFilename("model_problem_results2.dat");
    // organic_mater_1.Solve();
    // g2.set_style("points").plot_xy(organic_mater_1.mpGrid->xGrid,organic_mater_1.solution);
    // g2.set_style("lines").plot_equation("(4*exp(x) + exp(-4*x) ) / (4*exp(pi)+exp(-4*pi))-5*sin(x)-3*cos(x)","exact solution");

    return 0;
}