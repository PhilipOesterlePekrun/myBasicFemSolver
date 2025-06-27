#pragma once
#include <Global.hpp>

#include <myUtils.hpp>

#include "Element_Tri3.hpp"

namespace MyFem {

namespace Problem {
  
#include <LinAlg.hpp>

using namespace LinAlg;

class Linear2D {
// ctor(s)
 public:
  Linear2D() {};
  
  void runNoInputExample();
  void runNoInputExample1();
  void runNoInputExample2();
  void runNoInputExample_SingleEle();
  
 private:
  size_t nnode_;
  size_t ndofn_;// = 2;
  ///std::size_t nele_;
  
  Array<Element::Tri3*> elements_;
  Matrix2d K_;
  Vectord rhs_;
  
  Array<size_t> globalDofIds_; // The actual global dofs including dirichlet and solution. Shouldn't change after assembly gets that info from the elements. (in any case it should just be contiguous with size ndofn*nnode)
  Array<size_t> dirichletDofIds_;
  Vectord dirichletVect_;
  Vectord solutionVect_;
  
 private:
  Vectord X_0_;
  
 public:
  Vectord getX_0() {return X_0_;}
  Vectord getX_t();
  // Full solution vector including dirichlet
  Vectord fullSolution();
  
  size_t get_ndof() {return nnode_*ndofn_;}
  int get_nnode() {return nnode_;}
  int get_ndofn() {return ndofn_;}
  
  Array<Element::Tri3*> getElements() {
    return elements_;
  }
  
 private:
  void readMeshTxt(std::string inputFilePath);
  
  Matrix2d assembleK();
  //dynArrayd assembleRhs();// later when we have Wext
  
  void applyDirichlet(int globalNodeId/*, dofIndex or localDofindex of the node, needed for higher dim*/, double val);
  void applyDirichlet(const Array<size_t>& ids, const Vectord& vals);
  
  // TODO we will eventually put all the solver methods into a different class. This class is just a manager which also does the assembly.
  ///Vectord solveSystem_gaussSeidel();
  Vectord solveSystem_Jacobi(int maxiter, double maxResNorm = -1.0);
  Vectord solveSystem_GaussSeidel(int maxiter, double maxResNorm = -1.0);
  
  
  // TODO: we move this suff outside later I suppose
  void visualize();
}; // class Linear1D

} // namespace Problem

} // namespace MyFem
