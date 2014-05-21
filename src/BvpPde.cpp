#include "BvpPde.hpp"
using namespace Eigen;
using namespace std;

BvpPde::BvpPde(SecondOrderOde *pOde, BoundaryConditions *pBcs, double dt, double tau, double timer, int numNodes, double (*uj0)(double)) : BvpOde(pOde, pBcs, numNodes){
    mtau = tau;
    mdt = dt;
    mUj0Func = uj0;
    T = timer;
}

void BvpPde::SolvePde(){
    PopulateOperators();
    PopulateInitCondiotions();
    SolvePdeInTime();
}

void BvpPde::SolvePdeInTime(){
    int columns = T/mdt;
    cout << mb;
    solutionInTime.resize(BvpOde::mNumNodes, columns+1);
    solutionInTime.col(0) = mb;
    // Applying BC
    mb(0) = (*BvpOde::mpRhsVec)(0);
    mb(mNumNodes-1) = (*BvpOde::mpRhsVec)(mNumNodes-1);
    // cout << mb;
    // cout << mb;
    for (int i = 1; i <2; ++i) {
        VectorXd temp = mUj0Operator*mb;
        // cout << temp;
        temp(0) = (*BvpOde::mpRhsVec)(0);
        temp(mNumNodes-1) = (*BvpOde::mpRhsVec)(mNumNodes-1);
        mpLinearSolver = new LinearSolver(mUj1Operator, temp);
        // cout << "mUj1Operator" << endl << mUj1Operator << endl;
        // cout << "temp"<<endl << temp << endl;
        // NOTE: solver works correctly (Checked with MATLAB)
        solution = mpLinearSolver->SolveLinearSystem();
        mb = solution;
        // solutionInTime.col(i) = solution;
        // cout << "mb"<<endl << mb << endl;
        mb(0) = (*BvpOde::mpRhsVec)(0);
        mb(mNumNodes-1) = (*BvpOde::mpRhsVec)(mNumNodes-1);
    }
    // cout << mb;
}

void BvpPde::PopulateOperators(){
    int x = BvpOde::mNumNodes;
    SparseMatrix<double> I = MatrixXd::Identity(x, x).sparseView();
    SparseMatrix<double> F = *BvpOde::mpLhsMat;
    mUj0Operator = (I +     mtau*mdt*F);
    mUj1Operator = (I - (1-mtau)*mdt*F);
    cout << "mUj1Operator" << endl << mUj1Operator << endl;
}

void BvpPde::PopulateInitCondiotions(){
    int numNodes = BvpOde::mNumNodes;
    mb.resize(mNumNodes);
    for (int i = 0; i < numNodes; ++i){
        double x = BvpOde::mpGrid->mNodes[i].C.x;
        if (mUj0Func!=0) {
            mb(i) = mUj0Func(x);
        } else {
            mb(i) = 0;
        }
    }
}


