#include "cassert"
#include "BoundaryConditions.hpp"

BoundaryConditions::BoundaryConditions() {
    mYNBcIsDirichlet = false;
    mYNBcIsDirichlet = false;
    mYNBcIsNeumann = false;
    mY0BcIsNeumann = false;
    mX0BcIsDirichlet = false;
    mXNBcIsDirichlet = false;
    mX0BcIsNeumann = false;
    mXNBcIsNeumann = false;
}

// X space
void BoundaryConditions::SetX0DirichletBc1D(double x0Value)
{
    assert(!mX0BcIsNeumann);
    mX0BcIsDirichlet = true;
    mX0BcValue = x0Value;
}

void BoundaryConditions::SetXNDirichletBc1D(double xNValue){
    assert(!mXNBcIsNeumann);
    mXNBcIsDirichlet = true;
    mXNBcValue = xNValue;
}

void BoundaryConditions::SetX0NeumannBc1D(double x0DerivValue){
    assert(!mX0BcIsDirichlet);
    mX0BcIsNeumann = true;
    mX0BcValue = x0DerivValue;
}

void BoundaryConditions::SetXNNeumannBc1D(double xNDerivValue){
    assert(!mXNBcIsDirichlet); 
    mXNBcIsNeumann = true;
    mXNBcValue = xNDerivValue; 
}

void BoundaryConditions::SetX0DirichletBc2D(double (*leftBCFunc)(double)){
    assert(!mX0BcIsNeumann);
    mX0BcIsDirichlet = true;
    mpX0BcFunc2D = leftBCFunc;
}

void BoundaryConditions::SetXNDirichletBc2D(double (*NBCFunc)(double)){
    assert(!mXNBcIsNeumann);
    mXNBcIsDirichlet = true;
    mpXNBcFunc2D = NBCFunc;
}


// Y space
void BoundaryConditions::SetY0DirichletBc(double y0Value)
{
    assert(!mY0BcIsNeumann);
    mY0BcIsDirichlet = true;
    mY0BcValue = y0Value;
}

void BoundaryConditions::SetYNDirichletBc(double yNValue){
    assert(!mYNBcIsNeumann);
    mYNBcIsDirichlet = true;
    mYNBcValue = yNValue;
}

void BoundaryConditions::SetY0NeumannBc(double y0DerivValue){
    assert(!mYNBcIsDirichlet);
    mY0BcIsNeumann = true;
    mY0BcValue = y0DerivValue;
}

void BoundaryConditions::SetYNNeumannBc(double yNDerivValue){
    assert(!mYNBcIsDirichlet); 
    mYNBcIsNeumann = true;
    mYNBcValue = yNDerivValue; 
}