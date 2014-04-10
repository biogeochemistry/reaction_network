#include "BvpOde.hpp"

void BvpOde::PopulateMatrix(){
  mpLhsMat->resize(mNumNodes, mNumNodes);
}