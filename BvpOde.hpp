#ifndef BVP_ODE_HPP
#define BVP_ODE_HPP

#include "gtest/gtest_prod.h"
#include <string>
#include "FiniteDifferenceGrid.hpp"
#include "SecondOrderOde.hpp"
#include "BoundaryConditions.hpp"
#include "Eigen/Sparse"

using namespace Eigen;
using namespace std;

class BvpOde {
    private:
        int mNumNodes;
        FiniteDifferenceGrid* mpGrid;
        SecondOrderOde* mpOde;
        BoundaryConditions* mpBconds;
        VectorXd *mpSolVec, *mpRhsVec;
        SparseMatrix<double> *mpLhsMat;
        void PopulateMatrix();
        void PopulateVector();
        void ApplyBoundaryConditions();

    public:
        BvpOde(SecondOrderOde *pOde, BoundaryConditions *pBcs, int numNodes);
        ~BvpOde();
        void Solve();
        void Solve_cg();
        void Solve_sldlt();
        void Solve_sllt();

};

#endif // BVP_ODE_HPP