#pragma once
#include "Global.hpp"

#include "mu.hpp"
#include "mu_core_LinAlg.hpp"

#include "Element_Tri3.hpp"

namespace MyFem {

namespace Problem {

using namespace MyUtils;
using namespace LinAlg;

class Linear2D {
// ctor(s)
 public:
  Linear2D() {};
  
  void example_beam(double lx, double ly, int nx, int ny, int maxIter = 200, double tol = 1e-4, double density=0);
  
  // nc = n circumferential; nr = n radial
  // n = number of nodes, not eles
  void example_torus(double x0, double y0, double ri, double ro, int nc, int nr);
  
  void runNoInputExample();
  void runNoInputExample1();
  void runNoInputExample2();
  void runNoInputExample_SingleEle();
  
 private:
  size_t nnode_;
  size_t ndofn_;// = 2;
  ///std::size_t nele_;
  
  void setNodeDofInfo();
  
  std::vector<Element::Tri3*> elements_;
  Matrix2d K_;
  Matrix2d M_;
  Vectord rhs_;
  
  std::vector<size_t> globalDofIds_; // The actual global dofs including dirichlet and solution. Shouldn't change after assembly gets that info from the elements. (in any case it should just be contiguous with size ndofn*nnode)
  
  std::vector<size_t> dirichletDofIds_;
  Vectord dirichletVect_;
  std::vector<size_t> neumannDofIds_;
  Vectord neumannVect_;
  
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
  
  std::vector<Element::Tri3*> getElements() {
    return elements_;
  }
  
 private:
  void readMeshTxt(std::string inputFilePath);
  
  void assembleK();
  void assembleM();
  void assembleFGravity(double gravityAccel = 9.81);
  void assembleFull(double gravityAccel = 9.81);
  //dynArrayd assembleRhs();// later when we have Wext?
  
  void applyNeumann(const std::vector<size_t>& ids, const Vectord& vals);
  void applyDirichlet(const std::vector<size_t>& ids, const Vectord& vals);
  void applyBCs(const std::vector<size_t>& neumannIds, const Vectord& neumannVals, const std::vector<size_t>& dirichletIds, const Vectord& dirichletVals);
  
 public:
  std::string infoString();
  void printInfo();
}; // class Linear1D

} // namespace Problem

} // namespace MyFem
