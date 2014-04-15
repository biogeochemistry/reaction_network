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
}

void BvpOde::WriteSolutionFile() {
    std::ofstream output_file(mFilename.c_str()); 
    assert(output_file.is_open());
    for (int i=0; i<mNumNodes; i++) {
        double x = mpGrid->mNodes[i].coordinate;
          output_file << x << "  " << mpSolVec(i) << "\n";
    }
   output_file.flush();
   output_file.close();
   std::cout<<"Solution written to "<<mFilename<<"\n";
}
