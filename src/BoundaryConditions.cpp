#include "cassert"
#include "BoundaryConditions.hpp"
#include "iostream"

BoundaryConditions::BoundaryConditions() {
    mYNBcisDirichlet = false;
    mY0BcisDirichlet = false;
    mYNBcIsNeumann = false;
    mY0BcIsNeumann = false;
    mX0BcisDirichlet = false;
    mXNBcisDirichlet = false;
    mX0BcIsNeumann = false;
    mXNBcIsNeumann = false;
    mX0BcIsRobin = false;
    mXNBcIsRobin = false;
}


// =================1D======================
/*
in 1D case method "set_bc..." assign the value provided by user
to the variable mBCvalue.
*/

// X space 1D
void BoundaryConditions::SetX0DirichletBc1D(double x0Value)
{
    assert(!mX0BcIsNeumann  & !mX0BcIsRobin);
    mX0BcisDirichlet = true;
    mX0BcValue = x0Value;
}

void BoundaryConditions::SetXNDirichletBc1D(double xNValue){
    assert(!mXNBcIsNeumann  & !mXNBcIsRobin);
    mXNBcisDirichlet = true;
    mXNBcValue = xNValue;
}

void BoundaryConditions::SetX0NeumannBc1D(double x0DerivValue){
    assert(!mX0BcisDirichlet & !mX0BcIsRobin);
    mX0BcIsNeumann = true;
    mX0BcValue = x0DerivValue;
}

void BoundaryConditions::SetXNNeumannBc1D(double xNDerivValue){
    assert(!mXNBcisDirichlet & !mXNBcIsRobin); 
    mXNBcIsNeumann = true;
    mXNBcValue = xNDerivValue; 
}

void BoundaryConditions::SetX0RobinBc1D(double x0RobinValue){
    assert(!mX0BcisDirichlet & !mX0BcIsNeumann);
    mX0BcIsRobin = true;
    mX0BcValue = x0RobinValue;
}

void BoundaryConditions::SetXNRobinBc1D(double xNRobinValue){
    assert(!mXNBcisDirichlet & !mXNBcIsNeumann); 
    mXNBcIsRobin = true;
}


// Y space 1D
void BoundaryConditions::SetY0DirichletBc1D(double y0Value)
{
    assert(!mY0BcIsNeumann);
    mY0BcisDirichlet = true;
    mY0BcValue = y0Value;
}

void BoundaryConditions::SetYNDirichletBc1D(double yNValue){
    assert(!mYNBcIsNeumann);
    mYNBcisDirichlet = true;
    mYNBcValue = yNValue;
}

void BoundaryConditions::SetY0NeumannBc1D(double y0DerivValue){
    assert(!mYNBcisDirichlet);
    mY0BcIsNeumann = true;
    mY0BcValue = y0DerivValue;
}

void BoundaryConditions::SetYNNeumannBc1D(double yNDerivValue){
    assert(!mYNBcisDirichlet); 
    mYNBcIsNeumann = true;
    mYNBcValue = yNDerivValue; 
}


// =================2D======================
/*
In 2d case method "set_boundary_condition" returns pointer to the 
boundary condition which user should provide, it is not assigning of the 
variable as it is in 1D case.
*/
// X space 2D
void BoundaryConditions::SetX0DirichletBc2D(double (*X0BCFunc)(double)){
    assert(!mX0BcIsNeumann);
    mX0BcisDirichlet = true;
    mpX0BcFunc2D = X0BCFunc;
}

void BoundaryConditions::SetXNDirichletBc2D(double (*XNBCFunc)(double)){
    assert(!mXNBcIsNeumann);
    mXNBcisDirichlet = true;
    mpXNBcFunc2D = XNBCFunc;
}

void BoundaryConditions::SetX0NeumannBc2D(double (*X0DerivBCFunc)(double)){
    assert(!mX0BcisDirichlet);
    mX0BcIsNeumann = true;
    mpX0BcFunc2D = X0DerivBCFunc;
}

void BoundaryConditions::SetXNNeumannBc2D(double (*XNDerivBCFunc)(double)){
    assert(!mXNBcisDirichlet); 
    mXNBcIsNeumann = true;
    mpXNBcFunc2D = XNDerivBCFunc; 
}

// Y space 2D
void BoundaryConditions::SetY0DirichletBc2D(double (*Y0BCFunc)(double)){
    assert(!mY0BcIsNeumann);
    mY0BcisDirichlet = true;
    mpY0BcFunc2D = Y0BCFunc;
}

void BoundaryConditions::SetYNDirichletBc2D(double (*YNBCFunc)(double)){
    assert(!mYNBcIsNeumann);
    mYNBcisDirichlet = true;
    mpYNBcFunc2D = YNBCFunc;
}

void BoundaryConditions::SetY0NeumannBc2D(double (*X0DerivBCFunc)(double)){
    assert(!mY0BcisDirichlet);
    mY0BcIsNeumann = true;
    mpY0BcFunc2D = X0DerivBCFunc;
}

void BoundaryConditions::SetYNNeumannBc2D(double (*XNDerivBCFunc)(double)){
    assert(!mYNBcisDirichlet); 
    mYNBcIsNeumann = true;
    mpYNBcFunc2D = XNDerivBCFunc; 
}