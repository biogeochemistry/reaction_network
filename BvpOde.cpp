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
    MatCreate(PETSC_COMM_WORLD,&mpLhsMat);
    MatSetSizes(mpLhsMat,PETSC_DECIDE,PETSC_DECIDE,mNumNodes,mNumNodes);
    MatSetFromOptions(mpLhsMat);
    mFilename = "ode_output.dat";
}

BvpOde::~BvpOde() { 
    delete mpGrid;
}

void BvpOde::Solve(){

    PopulateMatrix();
    PopulateVector();
    ApplyBoundaryConditions();
    // MatView(mpLhsMat,PETSC_VIEWER_STDOUT_WORLD);
    // VecView(mpRhsVec,PETSC_VIEWER_STDOUT_WORLD);
    Solve_petsc();
    VecView(mpSolVec,PETSC_VIEWER_STDOUT_WORLD);
    WriteSolutionFile();
}

void BvpOde::PopulateMatrix() {
    MatSetUp(mpLhsMat);
    PetscInt col[3];
    PetscScalar temp[3];
    for (PetscInt i=1; i<mNumNodes-1; i++) {
        double xm = mpGrid->mNodes[i-1].coordinate; 
        double x = mpGrid->mNodes[i].coordinate; 
        double xp = mpGrid->mNodes[i+1].coordinate; 
        double alpha = 2.0/(xp-xm)/(x-xm);
        double beta = -2.0/(xp-x)/(x-xm);
        double gamma = 2.0/(xp-xm)/(xp-x);
        col[0]=i-1; col[1]=i; col[2]=i+1;
        temp[0] = (mpOde->mCoeffOfUxx)*alpha - (mpOde->mCoeffOfUx)/(xp-xm);
        temp[1] = (mpOde->mCoeffOfUxx)*beta + mpOde->mCoeffOfU;
        temp[2] = (mpOde->mCoeffOfUxx)*gamma +(mpOde->mCoeffOfUx)/(xp-xm);
        
        MatSetValues(mpLhsMat, 1, &i, 3, col, temp, INSERT_VALUES);
    }
}

void BvpOde::PopulateVector() {
    PetscScalar temp[mNumNodes];
    for (PetscInt i=1; i<mNumNodes-1; i++) {
        double x = mpGrid->mNodes[i].coordinate;
        temp[i] = mpOde->mpRhsFunc(x);
        VecSetValues(mpRhsVec,1,&i,&temp[i],INSERT_VALUES);
    }

}

void BvpOde::ApplyBoundaryConditions() {
    bool left_bc_applied = false; 
    bool right_bc_applied = false;
    if (mpBconds->mLhsBcIsDirichlet) {
        PetscScalar temp[1];PetscInt col[1], row[1];
        assert(left_bc_applied == false);
        temp[0] = 1.0; col[0] = 0; row[0] = 0;
        MatSetValues(mpLhsMat, 1, row, 1, col, temp, INSERT_VALUES);
        temp[0] = mpBconds->mLhsBcValue; 
        VecSetValues(mpRhsVec,1,row,temp,INSERT_VALUES);
        left_bc_applied = true;
    }
    if (mpBconds->mRhsBcIsDirichlet) {
        PetscScalar temp[1];PetscInt col[1], row[1];
        assert(right_bc_applied == false);
        temp[0] = 1.0; col[0] = mNumNodes-1; row[0] = mNumNodes-1;
        MatSetValues(mpLhsMat, 1, row, 1, col, temp, INSERT_VALUES);
        temp[0] = mpBconds->mRhsBcValue; 
        VecSetValues(mpRhsVec,1,row,temp,INSERT_VALUES);
        right_bc_applied = true;
    }
    if (mpBconds->mLhsBcIsNeumann) {
        PetscScalar temp[2], temp1[1];PetscInt col[2], row[1];
        assert(left_bc_applied == false);
        double h = mpGrid->mNodes[1].coordinate - mpGrid->mNodes[0].coordinate;
        temp[0] = -1.0/h; col[0] = 0; row[0] = 0;
        temp[1] = 1.0/h; col[1] = 1;
        MatSetValues(mpLhsMat, 1, row, 2, col, temp, INSERT_VALUES);
        temp1[0] = mpBconds->mLhsBcValue; 
        VecSetValues(mpRhsVec,1,row,temp1,INSERT_VALUES);
        left_bc_applied = true;
    }

    if (mpBconds->mRhsBcIsNeumann) {
        PetscScalar temp[2], temp1[1];PetscInt col[2], row[1];
        assert(right_bc_applied == false);
        double h = mpGrid->mNodes[mNumNodes-1].coordinate - mpGrid->mNodes[mNumNodes-2].coordinate; 
        temp[0] = -1.0/h; col[0] = mNumNodes-2; row[0] = mNumNodes-1;
        temp[1] = 1.0/h;  col[1] = mNumNodes-1;
        MatSetValues(mpLhsMat, 1, row, 2, col, temp, INSERT_VALUES);
        temp1[0] = mpBconds->mRhsBcValue;
        VecSetValues(mpRhsVec,1,row,temp1,INSERT_VALUES);
        right_bc_applied = true;
    }
    
    /* need to assemble after setting values! do necessary
    message passing etc to propagate matrix to all ranks */
    MatAssemblyBegin(mpLhsMat,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(mpLhsMat,MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(mpRhsVec);
    VecAssemblyEnd(mpRhsVec);
}

void BvpOde::WriteSolutionFile() {
    std::ofstream output_file(mFilename.c_str()); 
    assert(output_file.is_open());
    for (int i=0; i<mNumNodes; i++) {
        double x = mpGrid->mNodes[i].coordinate;
          // output_file << x << "  " << mpSolVec(i) << "\n";
    }
    output_file.flush();
    output_file.close();
    std::cout<<"Solution written to "<<mFilename<<"\n";
}

void BvpOde::Solve_petsc() { 
    KSP ksp;   /* linear solver context */
    PC pc;     /* preconditioner context */
    KSPCreate(PETSC_COMM_WORLD,&ksp);
    KSPCreate(PETSC_COMM_WORLD,&ksp);
    /* operator is A matrix, also set matrix for preconditioning here */
    KSPSetOperators(ksp,mpLhsMat,mpLhsMat,DIFFERENT_NONZERO_PATTERN);
    /* get pc context from ksp context */
    KSPGetPC(ksp,&pc);
    /* set preconditioner type */
    PCSetType(pc,PCJACOBI);
    KSPSetTolerances(ksp,1e-6,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
    KSPSetFromOptions(ksp);
    /* solve! */
    KSPSolve(ksp,mpRhsVec,mpSolVec);
    // KSPView(ksp,PETSC_VIEWER_STDOUT_WORLD);
}
