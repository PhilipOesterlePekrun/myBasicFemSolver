#pragma once
#include <Global.hpp>

#include <myUtils.hpp>

#include "Element_Tri3.hpp"

namespace Problem {
  
#include <LinAlg.hpp>

using namespace LinAlg;

class Linear2D {
// ctor(s)
 public:
  Linear2D() {};
  
  void runNoInputExample();
  void runNoInputExample2();
  
 private:
  std::size_t nnode_;
  std::size_t ndofn_;// = 2;
  ///std::size_t nele_;
  
  Array<Element::Tri3*> elements_;
  Matrix2d K_;
  Vectord rhs_;
  
  Array<int> dirichletDofIds_;
  Vectord dirichletVect_;
  Array<int> solutionDofIds_; // with removed ids after applying dirichlet
  Vectord solutionVect_;
  
 public:
  Vectord X_0_;
  // Full solution vector including dirichlet
  Vectord fullSolution();
  
 private:
  void readMeshTxt(std::string inputFilePath);
  
  Matrix2d assembleK();
  //dynArrayd assembleRhs();// later when we have Wext
  
  void applyDirichlet(int globalNodeId/*, dofIndex or localDofindex of the node, needed for higher dim*/, double val);
  
  // TODO we will eventually put all the solver methods into a different class. This class is just a manager which also does the assembly.
  ///Vectord solveSystem_gaussSeidel();
  Vectord solveSystem_Jacobi(int maxiter, double maxResNorm = -1.0);
  
  
  // TODO: we move this suff outside later I suppose
  void visualize();
}; // class Linear1D

} // namespace Problem
