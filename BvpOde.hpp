#ifndef BVP_ODE_HPP
#define BVP_ODE_HPP

#include <string>
#include "FiniteDifferenceGrid.hpp"
#include "SecondOrderOde.hpp"
#include "BoundaryConditions.hpp"
#include "Eigen/Dense"
#include <petsc.h>
#include <fstream>

using namespace Eigen;
using namespace std;

class BvpOde {
    private:
        BvpOde(const BvpOde& otherBvpOde){};
        int mNumNodes;
        SecondOrderOde* mpOde;
        BoundaryConditions* mpBconds;
        Vec mpRhsVec;
        Mat mpLhsMat;
        void PopulateMatrix();
        void PopulateVector();
        void ApplyBoundaryConditions();
        void WriteSolutionFile();
        string mFilename;

    public:
        BvpOde(SecondOrderOde *pOde, BoundaryConditions *pBcs, int numNodes);
        ~BvpOde();
        void Solve();
        void Solve_cg();
        void Solve_sldlt();
        void Solve_sllt();
        void Solve_sqr();
        void SetFilename(const std::string& name){
            mFilename = name;
        }
        Vec mpSolVec;
        FiniteDifferenceGrid* mpGrid;

};

#endif // BVP_ODE_HPP