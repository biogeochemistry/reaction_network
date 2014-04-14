#ifndef BVP_ODE_HPP
#define BVP_ODE_HPP

#include "gtest/gtest_prod.h"
#include <string>
#include "FiniteDifferenceGrid.hpp"
#include "SecondOrderOde.hpp"
#include "BoundaryConditions.hpp"
#include "Eigen/Dense"
#include <fstream>

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
        void Solve_dense_ldlt(); // for symmetric positive defined Matrices
        void Solve_dense_fullPivHouseholderQr(); // for any, slowest but accurate.
        void Solve_dense_colPivHouseholderQR(); //for any faster than fullPivHouseholderQr
        void Solve_sparse_cg(); // ConjugateGradiend. for symmetric positive defined Matrices 
        void Solve_sparse_ldlt(); // SparseCholesky: Direct LDLt factorization. for symmetric positive defined Matrices
        void Solve_sparse_llt(); // SparseCholesky: Direct LLt factorization. for symmetric positive defined Matrices
        void Solve_sparse_qr(); // SparseQR: Any, rectangular. for symmetric positive defined Matrices
        void Solve_sparse_BiCGSTAB(); // A bi conjugate gradient stabilized solver for sparse square problems
        void Solve_sparse_lu(); // LU factorization. ****
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