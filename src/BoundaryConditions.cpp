#include "cassert"
#include "BoundaryConditions.hpp"
#include "iostream"

BoundaryConditions::BoundaryConditions() {
    mYNBcIsConst = false;
    mY0BcIsConst = false;
    mYNBcIsNoFlux = false;
    mY0BcIsNoFlux = false;
    mX0BcIsConst = false;
    mXNBcIsConst = false;
    mX0BcIsNoFlux = false;
    mXNBcIsNoFlux = false;
    mX0BcIsFlux = false;
    mXNBcIsFlux = false;
}


// =================1D======================
/*
in 1D case method "set_bc..." assign the value provided by user
to the variable mBCvalue.
*/

// X space 1D
void BoundaryConditions::SetX0ConstBc1D(double x0Value)
{
    assert(!mX0BcIsNoFlux  & !mX0BcIsFlux);
    mX0BcIsConst = true;
    mX0BcValue = x0Value;
}

void BoundaryConditions::SetXNConstBc1D(double xNValue){
    assert(!mXNBcIsNoFlux  & !mXNBcIsFlux);
    mXNBcIsConst = true;
    mXNBcValue = xNValue;
}

void BoundaryConditions::SetX0NoFluxBc1D(double x0DerivValue){
    assert(!mX0BcIsConst & !mX0BcIsFlux);
    mX0BcIsNoFlux = true;
    mX0BcValue = x0DerivValue;
}

void BoundaryConditions::SetXNNoFluxBc1D(double xNDerivValue){
    assert(!mXNBcIsConst & !mXNBcIsFlux); 
    mXNBcIsNoFlux = true;
    mXNBcValue = xNDerivValue; 
}

void BoundaryConditions::SetX0FluxBc1D(double x0FluxValue){
    assert(!mX0BcIsConst & !mX0BcIsNoFlux);
    mX0BcIsFlux = true;
    mX0BcValue = x0FluxValue;
}

void BoundaryConditions::SetXNFluxBc1D(double xNFluxValue){
    assert(!mXNBcIsConst & !mXNBcIsNoFlux); 
    mXNBcIsFlux = true;
}


// Y space 1D
void BoundaryConditions::SetY0ConstBc1D(double y0Value)
{
    assert(!mY0BcIsNoFlux);
    mY0BcIsConst = true;
    mY0BcValue = y0Value;
}

void BoundaryConditions::SetYNConstBc1D(double yNValue){
    assert(!mYNBcIsNoFlux);
    mYNBcIsConst = true;
    mYNBcValue = yNValue;
}

void BoundaryConditions::SetY0NoFluxBc1D(double y0DerivValue){
    assert(!mYNBcIsConst);
    mY0BcIsNoFlux = true;
    mY0BcValue = y0DerivValue;
}

void BoundaryConditions::SetYNNoFluxBc1D(double yNDerivValue){
    assert(!mYNBcIsConst); 
    mYNBcIsNoFlux = true;
    mYNBcValue = yNDerivValue; 
}


// =================2D======================
/*
In 2d case method "set_boundary_condition" returns pointer to the 
boundary condition which user should provide, it is not assigning of the 
variable as it is in 1D case.
*/
// X space 2D
void BoundaryConditions::SetX0ConstBc2D(double (*X0BCFunc)(double)){
    assert(!mX0BcIsNoFlux);
    mX0BcIsConst = true;
    mpX0BcFunc2D = X0BCFunc;
}

void BoundaryConditions::SetXNConstBc2D(double (*XNBCFunc)(double)){
    assert(!mXNBcIsNoFlux);
    mXNBcIsConst = true;
    mpXNBcFunc2D = XNBCFunc;
}

void BoundaryConditions::SetX0NoFluxBc2D(double (*X0DerivBCFunc)(double)){
    assert(!mX0BcIsConst);
    mX0BcIsNoFlux = true;
    mpX0BcFunc2D = X0DerivBCFunc;
}

void BoundaryConditions::SetXNNoFluxBc2D(double (*XNDerivBCFunc)(double)){
    assert(!mXNBcIsConst); 
    mXNBcIsNoFlux = true;
    mpXNBcFunc2D = XNDerivBCFunc; 
}

// Y space 2D
void BoundaryConditions::SetY0ConstBc2D(double (*Y0BCFunc)(double)){
    assert(!mY0BcIsNoFlux);
    mY0BcIsConst = true;
    mpY0BcFunc2D = Y0BCFunc;
}

void BoundaryConditions::SetYNConstBc2D(double (*YNBCFunc)(double)){
    assert(!mYNBcIsNoFlux);
    mYNBcIsConst = true;
    mpYNBcFunc2D = YNBCFunc;
}

void BoundaryConditions::SetY0NoFluxBc2D(double (*X0DerivBCFunc)(double)){
    assert(!mY0BcIsConst);
    mY0BcIsNoFlux = true;
    mpY0BcFunc2D = X0DerivBCFunc;
}

void BoundaryConditions::SetYNNoFluxBc2D(double (*XNDerivBCFunc)(double)){
    assert(!mYNBcIsConst); 
    mYNBcIsNoFlux = true;
    mpYNBcFunc2D = XNDerivBCFunc; 
}