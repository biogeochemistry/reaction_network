#ifndef BOUNDARYCONDITIONSHEADERDEF
#define BOUNDARYCONDITIONSHEADERDEF

#include "gtest/gtest_prod.h"


class BoundaryConditions {
    protected:
        
        bool mLhsBcIsDirichlet, mRhsBcIsDirichlet, mLhsBcIsNeumann, mRhsBcIsNeumann; 
    public:
};

#endif