#ifndef BVPPDEHPP
#define BVPPDEHPP

#include "BvpOde.hpp"


class BvpPde: public BvpOde {
    FRIEND_TEST(PDE_init, constractor);
    private:
        SparseMatrix<double> mUj0Operator;
        SparseMatrix<double> mUj1Operator;
        void PopulateInitCondiotions();
        void PopulateOperators();
        VectorXd mb;
        VectorXd mp_temp;
        double (*mUj0Func)(double x);
        BvpOde* ode; //
        double mtau; // tau, shift in time
        LinearSolver *mpLinearSolver;
    public:
        BvpPde(SecondOrderOde *pOde, BoundaryConditions *pBcs, double dt, double tau, double timer, int numNodes, double (*uj0)(double) = 0);
        void Solve();
        MatrixXd solutionInTime;
        VectorXd solution;
        double mdt; // dt
        double mT; // last point in time
        void SolvePde1TS(int i);
};
#endif // BVPPDEHPP
