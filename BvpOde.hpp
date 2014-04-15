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

class BvpOde {
    private:
        // test framework
        FRIEND_TEST(bvpode, error_of_the_solution);

        BvpOde(const BvpOde& otherBvpOde){};
        int mNumNodes;
        SecondOrderOde* mpOde;
        BoundaryConditions* mpBconds;
        VectorXd *mpRhsVec;
        SparseMatrix<double> *mpLhsMat;
        void PopulateMatrix();
        void PopulateVector();
        void ApplyBoundaryConditions();
        void WriteSolutionFile();
        LinearSolver *mpLinearSolver;
        string mFilename;
    public:
        BvpOde(SecondOrderOde *pOde, BoundaryConditions *pBcs, int numNodes);
        ~BvpOde();
        void Solve();
        void SetFilename(const std::string& name){
            mFilename = name;
        }
        VectorXd mpSolVec;
        FiniteDifferenceGrid* mpGrid;

};

#endif // BVP_ODE_HPP