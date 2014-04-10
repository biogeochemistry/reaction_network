#ifndef SECONDORDERODEHEADERDEF
#define SECONDORDERODEHEADERDEF

#include "gtest/gtest_prod.h"

class SecondOrderOde

{
  friend class BvpOde;
private:
  FRIEND_TEST(second_order_ode, assigning_var);
  //Coefficients on LHS of ODE
  double mCoeffOfUxx;
  double mCoeffOfUx; 
  double mCoeffOfU;
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