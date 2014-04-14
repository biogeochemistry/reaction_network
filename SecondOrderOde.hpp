#ifndef SECONDORDERODEHEADERDEF
#define SECONDORDERODEHEADERDEF


class SecondOrderOde

{
  friend class BvpOde;
    private:
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
    SecondOrderOde(double coeffUxx, double coeffUx, double coeffU, double (*righthandSide)(double), double xMinimum, double xMaximum)
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