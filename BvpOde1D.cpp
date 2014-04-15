#include <iostream>
#include <fstream>
#include <cassert>
#include "BvpOde1D.hpp"


BvpOde1D::BvpOde1D(SecondOrderOde1D* pOde,BoundaryConditions* pBcs, int numNodes){
    mpOde = pOde; 
    mpBconds = pBcs;
    mNumNodes = numNodes;
    mpGrid = new FiniteDifferenceGrid1D(mNumNodes, pOde->mXmin, pOde->mXmax);
    mpRhsVec = new VectorXd(mNumNodes);
    mpLhsMat = new SparseMatrix<double> (mNumNodes, mNumNodes);
    mFilename = "ode_output.dat";
}

BvpOde1D::~BvpOde1D() {
    delete mpRhsVec;
    delete mpLhsMat;
    delete mpGrid;
}

void BvpOde1D::PopulateMatrix() {
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

void BvpOde1D::PopulateVector() {
    for (int i=1; i<mNumNodes-1; i++) {
        double x = mpGrid->mNodes[i].coordinate;
        (*mpRhsVec)(i) = mpOde->mpRhsFunc(x);
    }
}

void BvpOde1D::ApplyBoundaryConditions() {
    bool left_bc_applied = false; 
    bool right_bc_applied = false;

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

void BvpOde1D::WriteSolutionFile() {
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

