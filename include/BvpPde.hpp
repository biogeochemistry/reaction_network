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
        void SolvePdeInTime();
        VectorXd mb;
        VectorXd mp_temp;
        MatrixXd solutionInTime;
        VectorXd solution;
        double (*mUj0Func)(double x);
        BvpOde* ode; //
        double mtau; // tau, shift in time
        double mdt; // dt
        double T; // last point in time
        LinearSolver *mpLinearSolver;
    public:
        BvpPde(SecondOrderOde *pOde, BoundaryConditions *pBcs, double dt, double tau, double timer, int numNodes, double (*uj0)(double) = 0);
        void SolvePde();
};
#endif // BVPPDEHPP
