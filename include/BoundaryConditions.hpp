#ifndef BOUNDARYCONDITIONSHEADERDEF
#define BOUNDARYCONDITIONSHEADERDEF

#include "gtest/gtest_prod.h"


class BoundaryConditions {
    FRIEND_TEST(boundary_conditions, constructor);
    FRIEND_TEST(boundary_conditions2d, constructor);    
    FRIEND_TEST(boundary_conditions, assigning_var);
    protected:
        
        bool mX0BcIsDirichlet, mXNBcIsDirichlet, mX0BcIsNeumann, mXNBcIsNeumann; 
        bool mY0BcIsDirichlet, mYNBcIsDirichlet, mYNBcIsNeumann, mY0BcIsNeumann; 
        double mX0BcValue, mXNBcValue, mYNBcValue, mY0BcValue;
        double (*mpX0BcFunc2D)(double y);
        double (*mpXNBcFunc2D)(double y);

    public:
        BoundaryConditions();
        friend class BvpOde;
        void SetX0DirichletBc1D(double x0Value);
        void SetXNDirichletBc1D(double xNValue);
        void SetX0NeumannBc1D(double x0DerivValue);
        void SetXNNeumannBc1D(double xNDerivValue);
        void SetYNDirichletBc(double yNValue);
        void SetY0DirichletBc(double y0Value);
        void SetYNNeumannBc(double yNDerivValue);
        void SetY0NeumannBc(double y0DerivValue);
        void SetX0DirichletBc2D(double (*leftBCFunction)(double));
        void SetXNDirichletBc2D(double (*NBCFunc)(double));
};

#endif