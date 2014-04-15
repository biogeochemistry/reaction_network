#ifndef BVP_ODE_1D_HPP
#define BVP_ODE_1D_HPP

#include "BvpOde.hpp"
#include "gtest/gtest_prod.h"
#include <string>
#include "FiniteDifferenceGrid.hpp"
#include "SecondOrderOde.hpp"
#include "BoundaryConditions.hpp"
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
    public:
        BvpOde1D(SecondOrderOde *pOde, BoundaryConditions *pBcs, int numNodes);
        ~BvpOde1D();


};

#endif // BVP_ODE_HPP