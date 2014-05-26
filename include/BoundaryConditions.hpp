#ifndef BOUNDARYCONDITIONSHEADERDEF
#define BOUNDARYCONDITIONSHEADERDEF

#include "gtest/gtest_prod.h"


class BoundaryConditions {
    FRIEND_TEST(boundary_conditions, constructor);
    FRIEND_TEST(boundary_conditions2d, constructor);    
    FRIEND_TEST(boundary_conditions, assigning_var);
    protected:
        
        bool mX0BcIsConst, mXNBcIsConst, mX0BcIsNoFlux, mXNBcIsNoFlux, mX0BcIsFlux, mXNBcIsFlux; 
        bool mY0BcIsConst, mYNBcIsConst, mYNBcIsNoFlux, mY0BcIsNoFlux; 
        double mX0BcValue, mXNBcValue, mYNBcValue, mY0BcValue, x0FluxValue, xNFluxValue;
        double (*mpX0BcFunc2D)(double y);
        double (*mpXNBcFunc2D)(double y);
        double (*mpY0BcFunc2D)(double x);
        double (*mpYNBcFunc2D)(double x);

    public:
        BoundaryConditions();
        friend class BvpOde;
        // 1D
        /*
        NOTE:NoFlux BC can be only = 0. 
         */
        void SetX0ConstBc1D(double x0Value);
        void SetXNConstBc1D(double xNValue);
        void SetX0NoFluxBc1D(double x0DerivValue);
        void SetXNNoFluxBc1D(double xNDerivValue);
        void SetX0FluxBc1D(double x0FluxValue);
        void SetXNFluxBc1D(double xNFluxValue);
        void SetYNConstBc1D(double yNValue);
        void SetY0ConstBc1D(double y0Value);
        void SetYNNoFluxBc1D(double yNDerivValue);
        void SetY0NoFluxBc1D(double y0DerivValue);
        // 2D BC as a function 
        void SetX0ConstBc2D(double (*leftBCFunc)(double));
        void SetXNConstBc2D(double (*XNBCFunc)(double));
        void SetX0NoFluxBc2D(double (*x0DerivBCFunc)(double));
        void SetXNNoFluxBc2D(double (*xNDerivBCFunc)(double));
        // 2D BC as a function 
        void SetY0ConstBc2D(double (*X0BCFunc)(double));
        void SetYNConstBc2D(double (*YNBCFunc)(double));
        void SetY0NoFluxBc2D(double (*Y0DerivBCFunc)(double));
        void SetYNNoFluxBc2D(double (*YNDerivBCFunc)(double));
};

#endif