#ifndef SECONDORDERODEHEADERDEF
#define SECONDORDERODEHEADERDEF

#include "gtest/gtest_prod.h"


class SecondOrderOde {
    friend class BvpOde;
    friend class BvpOde1D;  
    FRIEND_TEST(second_order_ode, assigning_var);
    protected:

        //Coefficients on LHS of ODE
        double mCoeffOfUxx; // diffusion
        double mCoeffOfUx;  // advection
        double mCoeffOfU; // function
        //Function on RHS of ODE
        double (*mpRhsFunc)(double x);
        //Interval for domain
        double mXmin; 
        double mXmax;
    public: 
    SecondOrderOde(double coeffUxx, double coeffUx, double coeffU, double (*righthandSide)(double), double xMinimum, double xMaximum);
};

#endif