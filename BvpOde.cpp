#include <iostream>
#include <fstream>
#include <cassert>
#include "BvpOde.hpp"

BvpOde::BvpOde(SecondOrderOde* pOde,BoundaryConditions* pBcs, int numNodes){
    mpOde = pOde; mpBconds = pBcs;
    mNumNodes = numNodes;
    mpGrid = new FiniteDifferenceGrid(mNumNodes, pOde->mXmin, pOde->mXmax);
    mpRhsVec = new VectorXd(mNumNodes);
    mpLhsMat = new SparseMatrix<double> (mNumNodes, mNumNodes);
    mFilename = "ode_output.dat";
}

BvpOde::~BvpOde() {
    delete mpRhsVec;
    delete mpLhsMat;
    delete mpGrid;
}

void BvpOde::Solve(){
    PopulateMatrix();
    PopulateVector();
    ApplyBoundaryConditions();
    // cout << *mpLhsMat << endl;
    // cout << *mpRhsVec<<endl;
    Solve_dense_colPivHouseholderQR();
    WriteSolutionFile();
}

void BvpOde::PopulateMatrix() {
    for (int i=1; i<mNumNodes-1; i++) {
        double xm = mpGrid->mNodes[i-1].coordinate; 
        double x = mpGrid->mNodes[i].coordinate; 
        double xp = mpGrid->mNodes[i+1].coordinate; 
        double diffusion_alpha = 2.0/(xp-xm)/(x-xm);
        double diffusion_beta = -2.0/(xp-x)/(x-xm);
        double diffusion_gamma = 2.0/(xp-xm)/(xp-x);
        (*mpLhsMat).insert(i,i-1) = (mpOde->mCoeffOfUxx)*diffusion_alpha - (mpOde->mCoeffOfUx)/(xp-xm);
        (*mpLhsMat).insert(i,i) = (mpOde->mCoeffOfUxx)*diffusion_beta + mpOde->mCoeffOfU;
        (*mpLhsMat).insert(i,i+1) = (mpOde->mCoeffOfUxx)*diffusion_gamma +(mpOde->mCoeffOfUx)/(xp-xm);
    }
}

void BvpOde::PopulateVector() {
    for (int i=1; i<mNumNodes-1; i++) {
        double x = mpGrid->mNodes[i].coordinate;
        (*mpRhsVec)(i) = mpOde->mpRhsFunc(x);
    }
}

void BvpOde::ApplyBoundaryConditions() {
    bool left_bc_applied = false; bool right_bc_applied = false;
    if (mpBconds->mLhsBcIsDirichlet) {
        (*mpLhsMat).insert(0,0) = 1.0;
        (*mpRhsVec)(0) = mpBconds->mLhsBcValue; 
        left_bc_applied = true;
    }
    if (mpBconds->mRhsBcIsDirichlet) {
        (*mpLhsMat).insert(mNumNodes-1,mNumNodes-1) = 1.0;
        (*mpRhsVec)(mNumNodes-1) = mpBconds->mRhsBcValue; 
        right_bc_applied = true;
    }
    if (mpBconds->mLhsBcIsNeumann) {
        assert(left_bc_applied == false);
        double h = mpGrid->mNodes[1].coordinate - mpGrid->mNodes[0].coordinate;
        (*mpLhsMat).insert(0,0) = -1.0/h; 
        (*mpLhsMat).insert(0,1) = 1.0/h;
        (*mpRhsVec)(0) = mpBconds->mLhsBcValue; 
        left_bc_applied = true;
    }

    if (mpBconds->mRhsBcIsNeumann) {
        assert(right_bc_applied == false);
        double h = mpGrid->mNodes[mNumNodes-1].coordinate - mpGrid->mNodes[mNumNodes-2].coordinate; 
        (*mpLhsMat).insert(mNumNodes-1,mNumNodes-2) = -1.0/h; 
        (*mpLhsMat).insert(mNumNodes-1,mNumNodes-1) = 1.0/h; 
        (*mpRhsVec)(mNumNodes-1) = mpBconds->mRhsBcValue; 
        right_bc_applied = true;
    }
}

void BvpOde::WriteSolutionFile() {
    std::ofstream output_file(mFilename.c_str()); 
    assert(output_file.is_open());
    for (int i=0; i<mNumNodes; i++) {
        double x = mpGrid->mNodes[i].coordinate;
          output_file << x << "  " << mpSolVec(i) << "\n";
    }
   output_file.flush();
   output_file.close();
   std::cout<<"Solution written to "<<mFilename<<"\n";
}

void BvpOde::Solve_sparse_cg() {
    // symmetric positive definite
    ConjugateGradient<SparseMatrix<double> > solver;
    solver.compute(*mpLhsMat);
    mpSolVec  = solver.solve(*mpRhsVec);
}

void BvpOde::Solve_sparse_ldlt() {
    // symmetric positive definite
    SimplicialLDLT<SparseMatrix<double> > solver;
    solver.compute(*mpLhsMat);
    if(solver.info()==Eigen::Success) {
        cout << 'Success';
    }
    mpSolVec  = solver.solve(*mpRhsVec);
}


void BvpOde::Solve_sparse_llt() {
    // symmetric positive definite
    SimplicialLLT<SparseMatrix<double> > solver;
    solver.compute(*mpLhsMat);
    mpSolVec = solver.solve(*mpRhsVec);
}

void BvpOde::Solve_dense_ldlt() {
    // // symmetric positive definite
    mpSolVec = MatrixXd(*mpLhsMat).ldlt().solve(*mpRhsVec);
}

void BvpOde::Solve_dense_fullPivHouseholderQr() {
    // for any matrix
    mpSolVec = MatrixXd(*mpLhsMat).fullPivHouseholderQr().solve(*mpRhsVec);
}

void BvpOde::Solve_dense_colPivHouseholderQR() {
    // for any matrix
    MatrixXd A = MatrixXd(*mpLhsMat);
    ColPivHouseholderQR<MatrixXd> dec(A);
    mpSolVec = dec.solve(*mpRhsVec);
}
