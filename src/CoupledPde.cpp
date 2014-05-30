#include "CoupledPde.hpp"
using namespace Eigen;
using namespace std;

CoupledPde::CoupledPde(BvpPde *FirstPde, BvpPde *SecondPde) {
    Pde1 = FirstPde;
    Pde2 = SecondPde;
    mT = Pde1->mT;
    mdt = Pde1->mdt;
    counter = mT/mdt;
}

void CoupledPde::SolveSystem(){
    int tid;
    for (int i = 0; i < counter-1; ++i) {
        #pragma omp parallel num_threads(2) private(tid)
        {
            tid = omp_get_thread_num();
            if (tid==0) {
                Pde1->SolvePde1TS(i);
            }
            
            if (tid==1)
            {
                Pde2->SolvePde1TS(i);
            }
        }
    }
    WriteSolution();
}

void CoupledPde::WriteSolution(){
    string name1 = Pde1->mFilename;
    string name2 = Pde2->mFilename;
    
    ofstream file(name1);
    if (file.is_open()) {
        file << Pde1->solutionInTime << '\n';
    }

    ofstream file2(name2);
    if (file2.is_open()) {
        file2 << Pde2->solutionInTime << '\n';
    }
}


