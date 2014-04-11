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
  VectorXd* mpSolVec;
  SparseMatrix<double> *mpLhsMat;
  string mFilename;

  void PopulateMatrix();
  void PopulateVector();
  void ApplyBoundaryConditions();

public:
    BvpOde(SecondOrderOde *pOde, BoundaryConditions *pBcs, int numNodes) {
      mpOde = pOde;
      mpBconds = pBcs;
      mNumNodes = numNodes;
      mpLhsMat = new SparseMatrix<double>;
      mpGrid = new FiniteDifferenceGrid(numNodes, mpOde->mXmin, mpOde->mXmax);
      PopulateMatrix();
    };
    ~BvpOde(){
      delete mpLhsMat;
      delete mpGrid;
    };

    void SetFileName(const string& name){
      mFilename = name;
    }
    void Solve();
    // void WriteSolutionFile();
};

#endif // BVP_ODE_HPP