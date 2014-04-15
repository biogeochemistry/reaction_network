#ifndef SECONDORDERODE1DHEADERDEF
#define SECONDORDERODE1DHEADERDEF

#include "gtest/gtest_prod.h"
#include "SecondOrderOde.hpp"


class SecondOrderOde1D : public SecondOrderOde {
     
    
    protected:
        // test framework
        FRIEND_TEST(second_order_ode, assigning_var);
        friend class BvpOde;
        friend class BvpOde1D;
        // friend class FiniteDifferenceGrid;
    public: 
    SecondOrderOde1D(double coeffUxx, double coeffUx, double coeffU, double (*righthandSide)(double), double xMinimum, double xMaximum)
    {
        mCoeffOfUxx = coeffUxx;
        mCoeffOfUx = coeffUx;
        mCoeffOfU = coeffU;
        mpRhsFunc = righthandSide;
        mXmin = xMinimum;
        mXmax = xMaximum;
    };
};

#endif