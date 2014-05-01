#include "cassert"
#include "BoundaryConditions.hpp"

BoundaryConditions::BoundaryConditions() {
    mTopBcIsDirichlet = false;
    mBotBcIsDirichlet = false;
    mTopBcIsNeumann = false;
    mBotBcIsNeumann = false;
    mLhsBcIsDirichlet = false;
    mRhsBcIsDirichlet = false;
    mLhsBcIsNeumann = false;
    mRhsBcIsNeumann = false;
}

// X space
void BoundaryConditions::SetLhsDirichletBc(double lhsValue)
{
    assert(!mLhsBcIsNeumann);
    mLhsBcIsDirichlet = true;
    mLhsBcValue = lhsValue;
}

void BoundaryConditions::SetRhsDirichletBc(double rhsValue){
    assert(!mRhsBcIsNeumann);
    mRhsBcIsDirichlet = true;
    mRhsBcValue = rhsValue;
}

void BoundaryConditions::SetLhsNeumannBc(double lhsDerivValue){
    assert(!mLhsBcIsDirichlet);
    mLhsBcIsNeumann = true;
    mLhsBcValue = lhsDerivValue;
}

void BoundaryConditions::SetRhsNeumannBc(double rhsDerivValue){
    assert(!mRhsBcIsDirichlet); 
    mRhsBcIsNeumann = true;
    mRhsBcValue = rhsDerivValue; 
}


// Y space
void BoundaryConditions::SetBotDirichletBc(double BotValue)
{
    assert(!mBotBcIsNeumann);
    mBotBcIsDirichlet = true;
    mBotBcValue = BotValue;
}

void BoundaryConditions::SetTopDirichletBc(double TopValue){
    assert(!mTopBcIsNeumann);
    mTopBcIsDirichlet = true;
    mTopBcValue = TopValue;
}

void BoundaryConditions::SetBotNeumannBc(double BotDerivValue){
    assert(!mBotBcIsDirichlet);
    mBotBcIsNeumann = true;
    mBotBcValue = BotDerivValue;
}

void BoundaryConditions::SetTopNeumannBc(double TopDerivValue){
    assert(!mTopBcIsDirichlet); 
    mTopBcIsNeumann = true;
    mTopBcValue = TopDerivValue; 
}