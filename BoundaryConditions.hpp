#ifndef BOUNDARYCONDITIONSHEADERDEF
#define BOUNDARYCONDITIONSHEADERDEF

#include "gtest/gtest_prod.h"


class BoundaryConditions {
    private:
         // test framework
        FRIEND_TEST(boundary_conditions, constractor);
        
        bool mLhsBcIsDirichlet, mRhsBcIsDirichlet, mLhsBcIsNeumann, mRhsBcIsNeumann; 
        double mLhsBcValue, mRhsBcValue;
    public:
        friend class BvpOde; 
        BoundaryConditions();
        void SetLhsDirichletBc(double lhsValue);
        void SetRhsDirichletBc(double rhsValue);
        void SetLhsNeumannBc(double lhsDerivValue);
        void SetRhsNeumannBc(double rhsDerivValue);
};

#endif