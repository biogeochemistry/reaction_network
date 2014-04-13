#include <iostream>
#include <fstream>
#include <cassert>
#include "BvpOde.hpp"

BvpOde::BvpOde(SecondOrderOde* pOde,BoundaryConditions* pBcs, int numNodes){
    mpOde = pOde; 
    mpBconds = pBcs;
    mNumNodes = numNodes;
    mpGrid = new FiniteDifferenceGrid(mNumNodes, pOde->mXmin, pOde->mXmax);
    VecCreate(PETSC_COMM_WORLD, &mpRhsVec);
    VecSetSizes(mpRhsVec,PETSC_DECIDE,mNumNodes);
    VecSetFromOptions(mpRhsVec);
    VecDuplicate(mpRhsVec,&mpSolVec);

    mpLhsMat = new Mat(mNumNodes, mNumNodes);
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
    // cout << (*mpRhsVec)<< endl;
    Solve_sldlt();
    WriteSolutionFile();
}

void BvpOde::PopulateMatrix() {
    for (int i=1; i<mNumNodes-1; i++) {
        double xm = mpGrid->mNodes[i-1].coordinate; 
        double x = mpGrid->mNodes[i].coordinate; 
        double xp = mpGrid->mNodes[i+1].coordinate; 
        double alpha = 2.0/(xp-xm)/(x-xm);
        double beta = -2.0/(xp-x)/(x-xm);
        double gamma = 2.0/(xp-xm)/(xp-x);
        (*mpLhsMat).insert(i,i-1) = (mpOde->mCoeffOfUxx)*alpha - (mpOde->mCoeffOfUx)/(xp-xm);
        (*mpLhsMat).insert(i,i) = (mpOde->mCoeffOfUxx)*beta + mpOde->mCoeffOfU;
        (*mpLhsMat).insert(i,i+1) = (mpOde->mCoeffOfUxx)*gamma +(mpOde->mCoeffOfUx)/(xp-xm);
    }
}

void BvpOde::PopulateVector() {
    for (int i=1; i<mNumNodes-1; i++) {
        double x = mpGrid->mNodes[i].coordinate;
        (*mpRhsVec)(i) = mpOde->mpRhsFunc(x);
    }
}

void BvpOde::ApplyBoundaryConditions() {
    bool left_bc_applied = false; 
    bool right_bc_applied = false;
    if (mpBconds->mLhsBcIsDirichlet) {
        assert(left_bc_applied == false);
        (*mpLhsMat).insert(0,0) = 1.0;
        (*mpRhsVec)(0) = mpBconds->mLhsBcValue; 
        left_bc_applied = true;
    }
    if (mpBconds->mRhsBcIsDirichlet) {
        assert(right_bc_applied == false);
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

void BvpOde::Solve_cg() {
    ConjugateGradient<SparseMatrix<double> > cg;
    cg.compute(*mpLhsMat);
    mpSolVec  = cg.solve(*mpRhsVec);
}

void BvpOde::Solve_sldlt() {
    SimplicialLDLT<SparseMatrix<double> > sldlt;
    sldlt.compute(*mpLhsMat);
    mpSolVec  = sldlt.solve(*mpRhsVec);
}


void BvpOde::Solve_sllt() {
    SimplicialLLT<SparseMatrix<double> > sllt;
    sllt.compute(*mpLhsMat);
    mpSolVec = sllt.solve(*mpRhsVec);
}
