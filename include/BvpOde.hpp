#ifndef BVP_ODE_HPP
#define BVP_ODE_HPP

#include "gtest/gtest_prod.h"
#include <string>
#include "FiniteDifferenceGrid.hpp"
#include "SecondOrderOde.hpp"
#include "BoundaryConditions.hpp"
#include "Eigen/Dense"
#include <fstream>
#include "LinearSolver.hpp"

using namespace Eigen;
using namespace std;

class BvpOde{
    FRIEND_TEST(bvpode, error_of_the_solution);
    protected:
        // test framework

        int mNumNodes;
        VectorXd *mpRhsVec;
        SparseMatrix<double> *mpLhsMat;
        void PopulateMatrix();
        void PopulateMatrix6thOrder();
        void PopulateVector();
        void ApplyBoundaryConditions();
        void WriteSolutionFile();
        SecondOrderOde* mpOde;
        BoundaryConditions* mpBconds;
        LinearSolver *mpLinearSolver;
        string mFilename;
    public:
        BvpOde(SecondOrderOde *pOde, BoundaryConditions *pBcs, int numNodes);
        ~BvpOde();
        FiniteDifferenceGrid* mpGrid;

        void Solve();
        void SetFilename(const std::string& name);
        VectorXd solution;

};

#endif // BVP_ODE_HPP