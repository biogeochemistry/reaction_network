#ifndef SECONDORDERODEHEADERDEF
#define SECONDORDERODEHEADERDEF

#include "gtest/gtest_prod.h"


class SecondOrderOde {
    friend class BvpOde;
    friend class BvpOde1D;  
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
};

#endif