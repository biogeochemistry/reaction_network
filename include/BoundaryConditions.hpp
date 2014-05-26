#ifndef BOUNDARYCONDITIONSHEADERDEF
#define BOUNDARYCONDITIONSHEADERDEF

#include "gtest/gtest_prod.h"


class BoundaryConditions {
    FRIEND_TEST(boundary_conditions, constructor);
    FRIEND_TEST(boundary_conditions2d, constructor);    
    FRIEND_TEST(boundary_conditions, assigning_var);
    protected:
        
        bool mX0BcisDirichlet, mXNBcisDirichlet, mX0BcIsNeumann, mXNBcIsNeumann, mX0BcIsRobin, mXNBcIsRobin; 
        bool mY0BcisDirichlet, mYNBcisDirichlet, mYNBcIsNeumann, mY0BcIsNeumann; 
        double mX0BcValue, mXNBcValue, mYNBcValue, mY0BcValue, x0RobinValue, xNRobinValue;
        double (*mpX0BcFunc2D)(double y);
        double (*mpXNBcFunc2D)(double y);
        double (*mpY0BcFunc2D)(double x);
        double (*mpYNBcFunc2D)(double x);

    public:
        BoundaryConditions();
        friend class BvpOde;
        // 1D
        /*
        Robin BC doesnt work
         */
        void SetX0DirichletBc1D(double x0Value);
        void SetXNDirichletBc1D(double xNValue);
        void SetX0NeumannBc1D(double x0DerivValue);
        void SetXNNeumannBc1D(double xNDerivValue);
        void SetX0RobinBc1D(double x0RobinValue);
        void SetXNRobinBc1D(double xNRobinValue);
        void SetYNDirichletBc1D(double yNValue);
        void SetY0DirichletBc1D(double y0Value);
        void SetYNNeumannBc1D(double yNDerivValue);
        void SetY0NeumannBc1D(double y0DerivValue);
        // 2D BC as a function 
        void SetX0DirichletBc2D(double (*leftBCFunc)(double));
        void SetXNDirichletBc2D(double (*XNBCFunc)(double));
        void SetX0NeumannBc2D(double (*x0DerivBCFunc)(double));
        void SetXNNeumannBc2D(double (*xNDerivBCFunc)(double));
        // 2D BC as a function 
        void SetY0DirichletBc2D(double (*X0BCFunc)(double));
        void SetYNDirichletBc2D(double (*YNBCFunc)(double));
        void SetY0NeumannBc2D(double (*Y0DerivBCFunc)(double));
        void SetYNNeumannBc2D(double (*YNDerivBCFunc)(double));
};

#endif