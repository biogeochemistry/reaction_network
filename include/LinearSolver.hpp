#ifndef LINEAR_SOLVER_HPP
#define LINEAR_SOLVER_HPP
#include "Eigen/Dense"
#include "Eigen/Sparse"

using namespace Eigen;

class LinearSolver  {
public:
    LinearSolver(const SparseMatrix<double>& A,const VectorXd& b);
    ~LinearSolver();
    SparseMatrix<double> *mpA;
    VectorXd *mpb, *mpx;
    VectorXd SolveLinearSystem();
    void Solve_dense_ldlt(); // for symmetric positive defined Matrices
    void Solve_dense_fullPivHouseholderQr(); // for any, slowest but accurate.
    void Solve_dense_colPivHouseholderQR(); //for any faster than fullPivHouseholderQr
    void Solve_sparse_cg(); // ConjugateGradiend. for symmetric positive defined Matrices 
    void Solve_sparse_ldlt(); // SparseCholesky: Direct LDLt factorization. for symmetric positive defined Matrices
    void Solve_sparse_llt(); // SparseCholesky: Direct LLt factorization. for symmetric positive defined Matrices
    void Solve_sparse_qr(); // SparseQR: Any, rectangular. for symmetric positive defined Matrices
    void Solve_sparse_BiCGSTAB(); // A bi conjugate gradient stabilized solver for sparse square problems
    void Solve_sparse_lu(); // LU factorization. ****
};

#endif // LINEAR_SOLVER_HPP