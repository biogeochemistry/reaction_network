#include <iostream>
#include <fstream>
#include <cassert>
#include "BvpOde.hpp"


void BvpOde::Solve(){
    PopulateMatrix();
    PopulateVector();
    ApplyBoundaryConditions();
    mpLinearSolver = new LinearSolver(*mpLhsMat, *mpRhsVec);
    mpSolVec = mpLinearSolver->SolveLinearSystem();
    WriteSolutionFile();
}

BvpOde::BvpOde(SecondOrderOde* pOde,BoundaryConditions* pBcs, int numNodes){
    mpOde = pOde; 
    mpBconds = pBcs;
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

void BvpOde::PopulateMatrix() {
    for (int i=1; i<mNumNodes-1; i++) {
        double xm = mpGrid->mNodes[i-1].C.x; 
        double x = mpGrid->mNodes[i].C.x;
        double xp = mpGrid->mNodes[i+1].C.x;
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
        double x = mpGrid->mNodes[i].C.x;
        (*mpRhsVec)(i) = mpOde->mpRhsFunc(x);
    }
}

void BvpOde::ApplyBoundaryConditions() {
    bool left_bc_applied = false; 
    bool right_bc_applied = false;

    if (mpBconds->mX0BcIsDirichlet) {
        (*mpLhsMat).insert(0,0) = 1.0;
        (*mpRhsVec)(0) = mpBconds->mX0BcValue; 
        left_bc_applied = true;
    }

    if (mpBconds->mXNBcIsDirichlet) {
        (*mpLhsMat).insert(mNumNodes-1,mNumNodes-1) = 1.0;
        (*mpRhsVec)(mNumNodes-1) = mpBconds->mXNBcValue; 
        right_bc_applied = true;
    }

    if (mpBconds->mX0BcIsNeumann) {
        assert(left_bc_applied == false);
        double h = mpGrid->mNodes[1].C.x - mpGrid->mNodes[0].C.x;
        (*mpLhsMat).insert(0,0) = -1.0/h; 
        (*mpLhsMat).insert(0,1) = 1.0/h;
        (*mpRhsVec)(0) = mpBconds->mX0BcValue; 
        left_bc_applied = true;
    }

    if (mpBconds->mXNBcIsNeumann) {
        assert(right_bc_applied == false);
        double h = mpGrid->mNodes[mNumNodes-1].C.x - mpGrid->mNodes[mNumNodes-2].C.x; 
        (*mpLhsMat).insert(mNumNodes-1,mNumNodes-2) = -1.0/h; 
        (*mpLhsMat).insert(mNumNodes-1,mNumNodes-1) = 1.0/h; 
        (*mpRhsVec)(mNumNodes-1) = mpBconds->mXNBcValue; 
        right_bc_applied = true;
    }
}

void BvpOde::WriteSolutionFile() {
    std::ofstream output_file(mFilename.c_str()); 
    assert(output_file.is_open());
    for (int i=0; i<mNumNodes; i++) {
        double x = mpGrid->mNodes[i].C.x;
          output_file << x << "  " << mpSolVec(i) << "\n";
    }
   output_file.flush();
   output_file.close();
   std::cout<<"Solution written to "<<mFilename<<"\n";
}

void BvpOde::SetFilename(const std::string& name){
    mFilename = name;
}

