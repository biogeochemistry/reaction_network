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
    protected:
        // test framework

        int mNumNodes;
        VectorXd *mpRhsVec;
        SparseMatrix<double> *mpLhsMat;
        virtual void PopulateMatrix() = 0;
        virtual void PopulateVector() = 0;
        virtual void ApplyBoundaryConditions() = 0;
        virtual void WriteSolutionFile() = 0;
        LinearSolver *mpLinearSolver;
        string mFilename;
        BoundaryConditions* mpBconds;
    public:

        void Solve();
        void SetFilename(const std::string& name){
            mFilename = name;
        }
        VectorXd mpSolVec;

};

#endif // BVP_ODE_HPP