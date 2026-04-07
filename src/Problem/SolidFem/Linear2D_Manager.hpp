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
  
  void example_beam_static(double lx, double ly, int nx, int ny, int maxIter = 200, double tol = 1e-4);
  void example_beam_dyn(double lx, double ly, int nx, int ny, int maxIter = 200, double tol = 1e-4);
  
  // nc = n circumferential; nr = n radial
  // n = number of nodes, not eles
  ///void example_torus(double x0, double y0, double ri, double ro, int nc, int nr);
  
 private:
  size_t nnode_;
  size_t ndofn_;// = 2;
  ///std::size_t nele_;
  
  void setNodeDofInfo();
  
  vector<Element::Tri3*> elements_;
  Matrix2d KFull_; // without dirichlet applied
  Matrix2d KRed_; // with dirichlet applied
  Matrix2d MFull_;
  Matrix2d MRed_;
  Vectord rhsFull_;
  Vectord rhsDirich_; // with dirichlet applied (changes values)
  
  
  vector<size_t> globalDofIds_; // The actual global dofs including dirichlet and solution. Shouldn't change after assembly gets that info from the elements. (in any case it should just be contiguous with size ndofn*nnode)
  vector<size_t> dirichDofIds_; // CURRENT RESTRICTION: dirichletDofIds_ is the same for all time steps
  vector<size_t> freeDofIds_; // = {globalDofIds_} \ {dirichDofIds_}
  
  vector<Vectord> U_t_; // Displacement solution in time (not including dirichlet). If the problem is static, element 1 of the vector is the final solution to the static problem. Element 0 of the vector may be empty sometimes. TODOi: Fix the way the data and computation are handled so that it is guaranteed that all the data is computed
  
  vector<Vectord> UFull_t_; // Displacement solution in time (not including dirichlet). If the problem is static, element 1 of the vector is the final solution to the static problem. Element 0 of the vector may be empty sometimes. TODOi: Fix the way the data and computation are handled so that it is guaranteed that all the data is computed
  
 private:
  Vectord X_0_;
  Vectord V_0_;
  vector<Vectord> X_t_; // The actual position field, including dirichlet
  
 public:
  Vectord get_X_0() {return X_0_;}
  vector<Vectord> get_X_t() {return X_t_;}
  vector<Vectord> get_U_t() {return UFull_t_;} // Full solution including dirichlet
  
  // Computes and sets UFull_t_ and X_t_
  void setFinalSolution(const Vectord& URed_t, const vector<size_t>& dirichletDofIds, const Vectord& dirichletVals, size_t n);
  
  size_t get_ndof() const {return nnode_*ndofn_;}
  size_t get_nnode() const {return nnode_;}
  size_t get_ndofn() const {return ndofn_;}
  
  double finalT_;
  double deltaT_;
  int get_timeSteps() {return finalT_/deltaT_;}
  
  vector<Element::Tri3*> getElements() {
    return elements_;
  }
  
 private:
  void readMeshTxt(std::string inputFilePath);
  
  Matrix2d assembleKfull() const;
  Matrix2d assembleMfull() const;
  Vectord assembleFGravity(double gravityAccel = 9.81) const;
  void assembleAll(double gravityAccel = 9.81);
  //dynArrayd assembleRhs();// later when we have Wext?
  
  // Removes duplicates; if there are duplicates, the later one is ignored (no override)
  static void removeDuplicates(vector<size_t>& ids, Vectord& vals);
  
  void applyNeumannToRhsFull(Vectord& rhsFull, std::vector<size_t>& ids, Vectord& vals) const; // Constitutes applyNeumann entirely; duplicates allowed because we call removeDuplicates()
  inline void applyNeumann(std::vector<size_t>& ids, Vectord& vals) {
    applyNeumannToRhsFull(rhsFull_, ids, vals);
  }
    
  void setFreeAndDirichDofIds(std::vector<size_t>& ids, Vectord& vals);
  static Matrix2d reducedMatrix(const Matrix2d& A, const vector<size_t>& freeDofIds); // Duplicates allowed
  static Vectord reducedVector(const Vectord& x, const vector<size_t>& freeDofIds); // Duplicates allowed
  
  void applyDirichletToRhs(Vectord& rhsFull, const Matrix2d& KFull, vector<size_t>& dirichletDofIds, Vectord& dirichletVals) const; // Duplicates allowed because we call removeDuplicates()
  
  void applyDirichlet(vector<size_t>& ids, Vectord& vals);
  
  void applyBCs(vector<size_t>& neumannDofIds, Vectord& neumannVals, vector<size_t>& dirichletDofIds, Vectord& dirichletVals);
  
 public:
  std::string infoString();
  void printInfo();
}; // class Linear1D

} // namespace Problem

} // namespace MyFem
