#pragma once
#include <Global.hpp>

#include <myUtils.hpp>

#include "Element_Line2.hpp"

namespace Problem {
  
#include <LinAlg.hpp>

using namespace LinAlg;

class Linear1D {
 public:
  Linear1D() {};
  
  void runNoInputExample();
  
 private:
  Array<Element::Line2*> elements_;
  Matrixd K_;
  Vectord rhs_;
  
 private:
  void readMeshTxt(std::string inputFilePath);
  
  Matrixd assembleK();
  //dynArrayd assembleRhs();// later when we have Wext
  
  void applyDirichlet(int globalNodeId/*, dofIndex or localDofindex of the node, needed for higher dim*/, double val);
  
  // TODO we will eventually put all the solver methods into a different class. This class is just a manager which also does the assembly.
  ///Vectord solveSystem_gaussSeidel();
  Vectord solveSystem_Jacobi();
}; // class Linear1D

} // namespace Problem
