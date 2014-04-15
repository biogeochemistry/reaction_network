#ifndef BVP_ODE_1D_HPP
#define BVP_ODE_1D_HPP

#include "BvpOde.hpp"
#include "gtest/gtest_prod.h"
#include <string>
#include "FiniteDifferenceGrid1D.hpp"
#include "SecondOrderOde1D.hpp"
#include "BoundaryConditions1D.hpp"
#include "Eigen/Dense"
#include <fstream>

using namespace Eigen;
using namespace std;

class BvpOde1D : public BvpOde {
    private:
        // test framework
        FRIEND_TEST(bvpode, error_of_the_solution);

        BvpOde1D(const BvpOde1D& otherBvpOde1D){};
        void PopulateMatrix();
        void PopulateVector();
        void ApplyBoundaryConditions();
        void WriteSolutionFile();
        SecondOrderOde1D* mpOde;
        BoundaryConditions1D* mpBconds;
    public:
        BvpOde1D(SecondOrderOde1D *pOde, BoundaryConditions1D *pBcs, int numNodes);
        ~BvpOde1D();
        FiniteDifferenceGrid1D* mpGrid;

};

#endif // BVP_ODE_HPP