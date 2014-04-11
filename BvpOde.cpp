#include "BvpOde.hpp"

void BvpOde::PopulateMatrix(){
  (*mpLhsMat).resize(mNumNodes, mNumNodes);
  double diffusion = (*mpOde).mCoeffOfUxx;
  double advection = (*mpOde).mCoeffOfUx;
  double funct_coef = (*mpOde).mCoeffOfU;
  double diffusionLeft, diffusionCenter, diffusionRight, xLeft, xCenter, xRight, advectionLeft, advectionRight;


  // Formation of inner part of the matrix without BC
  for (int i = 1; i < mNumNodes-1; ++i)
  {
    xLeft = (*mpGrid).mNodes(i-1);
    xCenter = (*mpGrid).mNodes(i);
    xRight = (*mpGrid).mNodes(i+1);
    diffusionLeft   = 2.f * diffusion / ( (xRight - xLeft) * (xCenter - xLeft) );
    diffusionCenter = 2.f * diffusion / ( (xRight - xCenter) * (xCenter - xLeft));
    diffusionRight  = 2.f * diffusion / ( (xRight - xLeft) * (xRight - xCenter) );
    advectionLeft = - advection / (xRight - xLeft);
    advectionRight = advection / (xRight - xLeft);
    mpLhsMat->insert(i,i-1) = (diffusionLeft + advectionLeft) ;
    mpLhsMat->insert(i,i) = (diffusionCenter + funct_coef);
    mpLhsMat->insert(i,i+1) = (diffusionRight + advectionRight) ;
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
    mpLhsMat->insert(mNumNodes-1,mNumNodes-1) = - 1.f / ((*mpGrid).mNodes(mNumNodes-1) - (*mpGrid).mNodes(mNumNodes-2));
    mpLhsMat->insert(mNumNodes-1,mNumNodes-2) = 1.f / ((*mpGrid).mNodes(mNumNodes-1) - (*mpGrid).mNodes(mNumNodes-2));
  }
  else if ((*mpBconds).mRhsBcIsDirichlet == 1) {
    mpLhsMat->insert(mNumNodes-1,mNumNodes-1) = 1.f ;
  }
  else {
    cout << "Exception: Right boundary condition is not specified\n";
  }
cout << (*mpLhsMat);
}

