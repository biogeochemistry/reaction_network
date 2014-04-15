#ifndef BOUNDARYCONDITIONS1DHEADERDEF
#define BOUNDARYCONDITIONS1DHEADERDEF

#include "gtest/gtest_prod.h"
#include "BoundaryConditions.hpp"

class BoundaryConditions1D : public BoundaryConditions {
    private:
         // test framework
        FRIEND_TEST(boundary_conditions, constractor);
        
        double mLhsBcValue1D, mRhsBcValue1D;
    public:
        friend class BvpOde1D; 
        friend class BvpOde; 
        BoundaryConditions1D();
        void SetLhsDirichletBc1D(double lhsValue);
        void SetRhsDirichletBc1D(double rhsValue);
        void SetLhsNeumannBc1D(double lhsDerivValue);
        void SetRhsNeumannBc1D(double rhsDerivValue);
};

#endif