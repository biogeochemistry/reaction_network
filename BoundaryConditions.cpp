#include "BoundaryConditions.hpp"

BoundaryConditions::BoundaryConditions()
{
    mLhsBcIsDirichlet = false;
    mRhsBcIsDirichlet = false;
    mLhsBcIsNeumann = false;
    mRhsBcIsNeumann = false;
}

void BoundaryConditions::SetLhsDirichletBc(double lhsValue){
    mLhsBcIsDirichlet = true;
    mLhsBcValue = lhsValue;
}
void BoundaryConditions::SetRhsDirichletBc(double rhsValue){
    mRhsBcIsDirichlet = true;
    mRhsBcValue = rhsValue;
}
void BoundaryConditions::SetLhsNeumannBc(double lhsDerivValue){
    mLhsBcIsNeumann = true;
    mLhsBcValue = lhsDerivValue;
}
void BoundaryConditions::SetRhsNeumannBc(double rhsDerivValue){
    mRhsBcIsNeumann = true;
    mRhsBcValue = rhsDerivValue;
}