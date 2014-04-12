#include "BvpOde.hpp"
#include <Eigen/SparseCholesky>

void BvpOde::PopulateMatrix(){
  // (*mpLhsMat).resize(mNumNodes, mNumNodes);
  double diffusion_coef = (*mpOde).mCoeffOfUxx;
  double advection_coef = (*mpOde).mCoeffOfUx;
  double function_coef = (*mpOde).mCoeffOfU;
  double diffusionLeft, diffusionCenter, diffusionRight, xLeft, xCenter, xRight, advection;
  // cout << diffusion_coef << advection_coef << function_coef << endl;

  // Formation of inner part of the matrix without BC
  for (int i = 1; i < mNumNodes-1; ++i)
  {
    xLeft = (*mpGrid).mNodes(i-1);
    xCenter = (*mpGrid).mNodes(i);
    xRight = (*mpGrid).mNodes(i+1);
    diffusionLeft   = 2.f * diffusion_coef / (xRight - xLeft) / (xCenter - xLeft);
    diffusionCenter = - 2.f * diffusion_coef / (xRight - xCenter) / (xCenter - xLeft);
    diffusionRight  = 2.f * diffusion_coef /  (xRight - xLeft) / (xRight - xCenter);
    advection = advection_coef / (xRight - xLeft);
    mpLhsMat->insert(i,i-1) = (diffusionLeft - advection) ;
    mpLhsMat->insert(i,i) = (diffusionCenter + function_coef);
    mpLhsMat->insert(i,i+1) = (diffusionRight + advection) ;
  }
  
  // Formation of left BC
  if ((*mpBconds).mLhsBcIsNeumann == 1) {
    mpLhsMat->insert(0,0) = - 1.f / ((*mpGrid).mNodes(1) - (*mpGrid).mNodes(0));
    mpLhsMat->insert(0,1) = 1.f / ((*mpGrid).mNodes(1) - (*mpGrid).mNodes(0));
  }
  else if ((*mpBconds).mLhsBcIsDirichlet == 1) {
    mpLhsMat->insert(0,0) = 1.f ;
  }
  else {
    cout << "Exception: Left boundary condition is not specified\n";
  }

  // Formation of right BC
  if ((*mpBconds).mRhsBcIsNeumann == 1) {
    mpLhsMat->insert(mNumNodes-1,mNumNodes-1) = 1.f / ((*mpGrid).mNodes(mNumNodes-1) - (*mpGrid).mNodes(mNumNodes-2));
    mpLhsMat->insert(mNumNodes-1,mNumNodes-2) = - 1.f / ((*mpGrid).mNodes(mNumNodes-1) - (*mpGrid).mNodes(mNumNodes-2));
  }
  else if ((*mpBconds).mRhsBcIsDirichlet == 1) {
    mpLhsMat->insert(mNumNodes-1,mNumNodes-1) = 1.f ;
  }
  else {
    cout << "Exception: Right boundary condition is not specified\n";
  }
  // cout << (*mpLhsMat) << endl << endl;
}

void BvpOde::PopulateVector() {
  for (int i = 1; i < mNumNodes-1; ++i)
  {
    (*mpRhsVec)(i) = mpOde->mpRhsFunc(mpGrid->mNodes(i));
  }
  (*mpRhsVec)(0) = (*mpBconds).mLhsBcValue;
  (*mpRhsVec)(mNumNodes-1) = (*mpBconds).mRhsBcValue;
  // cout << (*mpRhsVec) << endl << endl;
}

void BvpOde::Solve_cg() {
  VectorXd x(mNumNodes);
  ConjugateGradient<SparseMatrix<double> > cg;
  cg.compute(*mpLhsMat);
  x = cg.solve(*mpRhsVec);
  mpSolVec = &x;
  cout << x << endl << endl;
}

void BvpOde::Solve_sldlt() {
  VectorXd x(mNumNodes);
  SimplicialLDLT<SparseMatrix<double> > sldlt;
  sldlt.compute(*mpLhsMat);
  x = sldlt.solve(*mpRhsVec);
  mpSolVec = &x;
  cout << (*mpSolVec) << endl << endl;
}


void BvpOde::Solve_sllt() {
  VectorXd x(mNumNodes);
  SimplicialLLT<SparseMatrix<double> > sllt;
  sllt.compute(*mpLhsMat);
  x = sllt.solve(*mpRhsVec);
  mpSolVec = &x;
  cout << (*mpSolVec) << endl << endl;
}