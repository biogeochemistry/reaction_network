#ifndef BVP_ODE_HPP
#define BVP_ODE_HPP

#include "gtest/gtest_prod.h"
#include <string>
#include "FiniteDifferenceGrid.hpp"
#include "SecondOrderOde.hpp"
#include "BoundaryConditions.hpp"
#include "Eigen/Sparse"

using namespace Eigen;
using namespace std;

class BvpOde {
private:
  BvpOde(const BvpOde& otherBvpOde) {}
  int mNumNodes;
  FiniteDifferenceGrid* mpGrid;
  SecondOrderOde* mpOde;
  BoundaryConditions* mpBconds;
  VectorXd *mpSolVec, *mpRhsVec;
  SparseMatrix<double> *mpLhsMat;
  void PopulateMatrix();
  void PopulateVector();
  void ApplyBoundaryConditions();

public:
    BvpOde(SecondOrderOde *pOde, BoundaryConditions *pBcs, int numNodes) {
      mpOde = pOde;
      mpBconds = pBcs;
      mNumNodes = numNodes;
      mpLhsMat = new SparseMatrix<double>(numNodes,numNodes);
      mpRhsVec = new VectorXd(numNodes);
      mpGrid = new FiniteDifferenceGrid(numNodes, mpOde->mXmin, mpOde->mXmax);
      PopulateMatrix();
      PopulateVector();
    };
    ~BvpOde(){
      delete mpLhsMat;
      delete mpGrid;
      delete mpRhsVec;
    };
    void Solve_cg();
    void Solve_sldlt();
    void Solve_sllt();

};

#endif // BVP_ODE_HPP