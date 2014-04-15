#include <iostream>
#include <fstream>
#include <cassert>
#include "BvpOde.hpp"


void BvpOde::Solve(){
    PopulateMatrix();
    PopulateVector();
    ApplyBoundaryConditions();
    mpLinearSolver = new LinearSolver(*mpLhsMat, *mpRhsVec);
    mpSolVec = mpLinearSolver->SolveLinearSystem();
    WriteSolutionFile();
}

