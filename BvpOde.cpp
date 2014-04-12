#include "BvpOde.hpp"
#include <Eigen/SparseCholesky>

BvpOde::BvpOde(SecondOrderOde *pOde, BoundaryConditions *pBcs, int numNodes){
    mpOde = pOde;
    mpBconds = pBcs;
    mNumNodes = numNodes;
    mpLhsMat = new SparseMatrix<double>(numNodes,numNodes);
    mpRhsVec = new VectorXd(numNodes);
    mpGrid = new FiniteDifferenceGrid(numNodes, mpOde->mXmin, mpOde->mXmax);
}

BvpOde::~BvpOde(){
    delete mpLhsMat;
    delete mpGrid;
    delete mpRhsVec;
};

void BvpOde::Solve(){
    PopulateMatrix();
    PopulateVector();
    ApplyBoundaryConditions();
    Solve_sldlt();
    WriteSolution();
    // cout << (*mpRhsVec)<< endl;
    // cout << mpSolVec << endl;
}

void BvpOde::WriteSolution(){
    std::ofstream output_file(mFilename.c_str()); 
    assert(output_file.is_open()); 
    for (int i=0; i<mNumNodes; i++) {
        double x = mpGrid->mNodes[i]; 
        output_file << x << " " << mpSolVec(i) << "\n"; 
    }
    output_file.flush(); 
    output_file.close();
    std::cout<<"Solution written to "<<mFilename<<"\n";
}
    

void BvpOde::PopulateMatrix(){
    // (*mpLhsMat).resize(mNumNodes, mNumNodes);
    double diffusion_coef = mpOde->mCoeffOfUxx;
    double advection_coef = mpOde->mCoeffOfUx;
    double function_coef = mpOde->mCoeffOfU;
    double diffusionLeft, diffusionCenter, diffusionRight, xLeft, xCenter, xRight, advection;
    // cout << diffusion_coef << advection_coef << function_coef << endl;

    // Formation of inner part of the matrix without BC
    for (int i = 1; i < mNumNodes-1; ++i)    {
        xLeft = mpGrid->mNodes(i-1);
        xCenter = mpGrid->mNodes(i);
        xRight = mpGrid->mNodes(i+1);
        diffusionLeft   = 2.0 * diffusion_coef / (xRight - xLeft) / (xCenter - xLeft);
        diffusionCenter = - 2.0 * diffusion_coef / (xRight - xCenter) / (xCenter - xLeft);
        diffusionRight  = 2.0 * diffusion_coef /  (xRight - xLeft) / (xRight - xCenter);
        advection = 1.0 * advection_coef / (xRight - xLeft);
        mpLhsMat->insert(i,i-1) = (diffusionLeft - advection) ;
        mpLhsMat->insert(i,i) = (diffusionCenter + function_coef);
        mpLhsMat->insert(i,i+1) = (diffusionRight + advection) ;
    }
}

void BvpOde::PopulateVector() {
    for (int i = 1; i < mNumNodes-1; ++i)   {
        double x = mpGrid->mNodes(i);
        (*mpRhsVec)(i) = mpOde->mpRhsFunc(x);
    }
    // cout << (*mpRhsVec) << endl << endl;
}

void BvpOde::ApplyBoundaryConditions(){
    bool left_bc_applied = false; 
    bool right_bc_applied = false;
    
    // Formation of left BC
    if (mpBconds->mLhsBcIsNeumann) {
        assert(left_bc_applied == false);
        mpLhsMat->insert(0,0) = - 1.f / (mpGrid->mNodes(1) - mpGrid->mNodes(0));
        mpLhsMat->insert(0,1) = 1.f / (mpGrid->mNodes(1) - mpGrid->mNodes(0));
        left_bc_applied = true;
    }
    else if (mpBconds->mLhsBcIsDirichlet) {
        assert(left_bc_applied == false);
        mpLhsMat->insert(0,0) = 1.f ;
        left_bc_applied = true;
    }
    else {cout << "Exception: Left boundary condition is not specified\n";
    }
    (*mpRhsVec)(0) = mpBconds->mLhsBcValue;

    // Formation of right BC
    if (mpBconds->mRhsBcIsNeumann) {
        assert(right_bc_applied == false);
        mpLhsMat->insert(mNumNodes-1,mNumNodes-1) = 1.f / (mpGrid->mNodes(mNumNodes-1) - mpGrid->mNodes(mNumNodes-2));
        mpLhsMat->insert(mNumNodes-1,mNumNodes-2) = - 1.f / (mpGrid->mNodes(mNumNodes-1) - mpGrid->mNodes(mNumNodes-2));
        right_bc_applied = false;
    }
    else if (mpBconds->mRhsBcIsDirichlet) {
        assert(right_bc_applied == false);
        mpLhsMat->insert(mNumNodes-1,mNumNodes-1) = 1.f ;
        right_bc_applied = false;
    }
    else {cout << "Exception: Right boundary condition is not specified\n";  }
    (*mpRhsVec)(mNumNodes-1) = mpBconds->mRhsBcValue;
    // cout << (*mpLhsMat) << endl << endl;
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
