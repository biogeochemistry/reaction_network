#ifndef BOUNDARYCONDITIONSHEADERDEF
#define BOUNDARYCONDITIONSHEADERDEF

#include "gtest/gtest_prod.h"


class BoundaryConditions {
    FRIEND_TEST(boundary_conditions, constructor);
    FRIEND_TEST(boundary_conditions2d, constructor);    
    protected:
        
        bool mLhsBcIsDirichlet, mRhsBcIsDirichlet, mLhsBcIsNeumann, mRhsBcIsNeumann; 
        bool mTopBcIsDirichlet, mBotBcIsDirichlet, mTopBcIsNeumann, mBotBcIsNeumann; 
        double mLhsBcValue, mRhsBcValue, mTopBcValue, mBotBcValue;
    public:
        BoundaryConditions();
        friend class BvpOde;
        void SetLhsDirichletBc(double lhsValue);
        void SetRhsDirichletBc(double rhsValue);
        void SetLhsNeumannBc(double lhsDerivValue);
        void SetRhsNeumannBc(double rhsDerivValue);
        void SetTopDirichletBc(double TopValue);
        void SetBotDirichletBc(double BotValue);
        void SetTopNeumannBc(double TopDerivValue);
        void SetBotNeumannBc(double BotDerivValue);
};

#endif