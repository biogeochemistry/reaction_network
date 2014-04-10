#ifndef BOUNDARYCONDITIONSHEADERDEF
#define BOUNDARYCONDITIONSHEADERDEF

#include "gtest/gtest_prod.h"

class BoundaryConditions

{
  public:
    friend class BvpOde; 
  private:
    FRIEND_TEST(Boundary_Conditions, def_constractor);
    bool mLhsBcIsDirichlet, mRhsBcIsDirichlet, mLhsBcIsNeumann, mRhsBcIsNeumann; 
    double mLhsBcValue, mRhsBcValue;
  public:
    BoundaryConditions();
    void SetLhsDirichletBc(double lhsValue);
    void SetRhsDirichletBc(double rhsValue);
    void SetLhsNeumannBc(double lhsDerivValue);
    void SetRhsNeumannBc(double rhsDerivValue);
};

#endif