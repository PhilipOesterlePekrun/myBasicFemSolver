#pragma once
#include <Global.hpp>

#include <myUtils.hpp>

#include "Element_line2.hpp"

namespace Problem {
  
#include <LinAlg.hpp>

class Linear1D {
  
 public:
  Linear1D() {};
  
  void runNoInputExample();
  
 private:
  dynArray<Element::line2*> elements_;
  dynMatrixd K_;
  dynArrayd rhs_;
  
 private:
  void readMeshTxt(std::string inputFilePath);
  
  dynMatrixd assembleK();
  //dynArrayd assembleRhs();// later when we have Wext
  
  void applyDirichlet(int globalNodeId/*, dofIndex or localDofindex of the node, needed for higher dim*/, double val);
  
  // TODO we will eventually put all the solver methods into a different class. This class is just a manager which also does the assembly.
  dynArrayd solveSystem_gaussSeidel();
  dynArrayd solveSystem_Jacobi();
}; // class Linear1D

} // namespace Problem
