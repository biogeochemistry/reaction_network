#include <cassert>
#include "BoundaryConditions1D.hpp"


BoundaryConditions1D::BoundaryConditions1D() {
    mLhsBcIsDirichlet = false;
    mRhsBcIsDirichlet = false;
    mLhsBcIsNeumann = false;
    mRhsBcIsNeumann = false;
}

void BoundaryConditions1D::SetLhsDirichletBc1D(double lhsValue)
{
    assert(!mLhsBcIsNeumann);
    mLhsBcIsDirichlet = true;
    mLhsBcValue1D = lhsValue;
}

void BoundaryConditions1D::SetRhsDirichletBc1D(double rhsValue){
    assert(!mRhsBcIsNeumann);
    mRhsBcIsDirichlet = true;
    mRhsBcValue1D = rhsValue;
}

void BoundaryConditions1D::SetLhsNeumannBc1D(double lhsDerivValue){
    assert(!mLhsBcIsDirichlet);
    mLhsBcIsNeumann = true;
    mLhsBcValue1D = lhsDerivValue;
}

void BoundaryConditions1D::SetRhsNeumannBc1D(double rhsDerivValue){
    assert(!mRhsBcIsDirichlet); 
    mRhsBcIsNeumann = true;
    mRhsBcValue1D = rhsDerivValue; 
}