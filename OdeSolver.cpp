#include "OdeSolver.hpp"

double OdeSolver::RightHandSide(double t, double y){
    return 1 + t;
}

void OdeSolver::WritingResultsIntoFile(std::ofstream &file, double t, double y) {
    file << t << "\t|\t" << y << "\n";
}

void OdeSolver::SetStepSize (double h) {
    stepSize = h;
}

void OdeSolver::SetTimeInterval (double t0, double t1) {
    initialTime = t0;
    finalTime = t1;
}

void OdeSolver::SetInitialValue (double y0) {
    initialValue = y0;
}