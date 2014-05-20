#include "iostream"
#include <cmath>
#include <string>
#include "BvpPde.hpp"
#include "gnuplot_i.hpp"
#include <omp.h>

using namespace Eigen;
using namespace std;


double model_prob_1_rhs(double x){return 1.0;}
double model_prob_2_rhs(double x){return 34.0*sin(x);}
double model_prob_3_rhs(double x, double y){return exp(-(x*x+y*y));}
double uj0(double x){return 0.5;}


int main(int argc, char* argv[]) {
    int nthreads, tid;
    #pragma omp parallel
    {
        /* Obtain thread number */
        tid = omp_get_thread_num();
        printf("Hello World from thread = %d\n", tid);

        /* Only master thread does this */
        if (tid == 0) {
            nthreads = omp_get_num_threads();
            printf("Number of threads = %d\n", nthreads);
        }
    }

    Gnuplot g1;

    SecondOrderOde ode_mp1(1.0, 0.0, 0.0, model_prob_1_rhs, 0.0, 1.0);
    BoundaryConditions bc_mp1;
    bc_mp1.SetX0DirichletBc1D(1);
    bc_mp1.SetXNNeumannBc1D(0);
    BvpPde bvpode_mp1(&ode_mp1, &bc_mp1,0.01, 1, 1, 8, uj0);
    bvpode_mp1.SetFilename("model_problem_results1.dat");
    bvpode_mp1.SolvePde();

    // g1.set_style("points").plot_xy(bvpode_mp1.mpGrid->xGrid,bvpode_mp1.mSolVec,"differentiation");
    // g1.set_style("lines").plot_equation("0.5*x*(1-x)","exact solution");

    // Gnuplot g2;
    // SecondOrderOde ode_mp2(1.0, 3.0, -4.0, model_prob_2_rhs, 0.0, M_PI);
    // BoundaryConditions bc_mp2;
    // bc_mp2.SetX0NeumannBc1D(-5.0);
    // bc_mp2.SetXNDirichletBc1D(4.0);
    // BvpOde bvpode_mp2(&ode_mp2, &bc_mp2, 64);
    // bvpode_mp2.SetFilename("model_problem_results2.dat");
    // bvpode_mp2.Solve();
    // g2.set_style("points").plot_xy(bvpode_mp2.mpGrid->xGrid,bvpode_mp2.mSolVec);
    // g2.set_style("lines").plot_equation("(4*exp(x) + exp(-4*x) ) / (4*exp(pi)+exp(-4*pi))-5*sin(x)-3*cos(x)","exact solution");

    return 0;
}