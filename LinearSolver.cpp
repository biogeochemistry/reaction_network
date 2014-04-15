#include "BvpOde.hpp"
#include "LinearSolver.hpp"

using namespace Eigen;

LinearSolver::LinearSolver(const SparseMatrix<double>& A,const VectorXd& b) {
    mpA = new SparseMatrix<double>(A);
    mpb = new VectorXd(b);
    mpx = new VectorXd(b);
}

LinearSolver::~LinearSolver(){
    delete mpA;
    delete mpb;
    delete mpx;
}

VectorXd LinearSolver::SolveLinearSystem(){
    Solve_sparse_lu();
    return *mpx;
}

void LinearSolver::Solve_sparse_cg() {
    // symmetric positive definite
    ConjugateGradient<SparseMatrix<double> > solver;
    solver.compute(*mpA);
    *mpx  = solver.solve(*mpb);
}

void LinearSolver::Solve_sparse_ldlt() {
    // symmetric positive definite
    SimplicialLDLT<SparseMatrix<double> > solver;
    solver.compute(*mpA);
    if(solver.info()==Eigen::Success) {
        std::cout << 'Success';
    }
    *mpx  = solver.solve(*mpb);
}


void LinearSolver::Solve_sparse_llt() {
    // symmetric positive definite
    SimplicialLLT<SparseMatrix<double> > solver;
    solver.compute(*mpA);
    *mpx = solver.solve(*mpb);
}

void LinearSolver::Solve_dense_ldlt() {
    // // symmetric positive definite
    *mpx = MatrixXd(*mpA).ldlt().solve(*mpb);
}

void LinearSolver::Solve_dense_fullPivHouseholderQr() {
    // for any matrix
    *mpx = MatrixXd(*mpA).fullPivHouseholderQr().solve(*mpb);
}

void LinearSolver::Solve_dense_colPivHouseholderQR() {
    // for any matrix
    MatrixXd A = MatrixXd(*mpA);
    ColPivHouseholderQR<MatrixXd> dec(A);
    *mpx = dec.solve(*mpb);
}

void LinearSolver::Solve_sparse_BiCGSTAB() {
    // for any sparse matrix
    BiCGSTAB<SparseMatrix<double> > solver;
    solver.compute(*mpA);
    *mpx  = solver.solve(*mpb);
}

void LinearSolver::Solve_sparse_lu() {
    // for any sparse matrix
    SparseLU<SparseMatrix<double,ColMajor>, AMDOrdering<int> > slu;
    slu.compute(*mpA);
    *mpx=slu.solve(*mpb);
}
