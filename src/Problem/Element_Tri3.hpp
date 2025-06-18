#pragma once
#include <Global.hpp>

#include <myUtils.hpp>

#include <LinAlg.hpp>

namespace Element {
  
using namespace LinAlg; // mTODO: consider removing this and replacing it by specific usings. Maybe with a good macro I can reuse in other files too; In any case, except for structures like matrix vector array, I will try to be explicit about LinAlg::

class Tri3 {
// Characteristic of the element type
 public:
  static constexpr int nnode_ = 3;
  static constexpr int ndofn_ = 2;
  static constexpr int ndof_ = ndofn_*nnode_;
 private:
  static constexpr int dim_ = ndofn_;
  static constexpr double Xi_[nnode_][ndofn_] = {{0, 0}, {1, 0}, {0, 1}};
  /*
      |\
    ^ | \
xi1 | |  \
    | |___\
      --->
       xi0 (and first dof)
       
   N2    
    |\
    | \
    |  \
    |___\
   N0   N1
  
  */


// Node/dof mappings
 public:
  // node and node-local dof to element-global dof
  int nodeldof2dof(int n, int i) {
    return n*ndofn_ + i;
  }
  int nodeldof2dof(Array<int> ni) {
    return ni(0)*ndofn_ + ni(1);
  }
  // element-global dof to node and node-local dof
  Array<int> dof2nodeldof(int g) {
    return Array<int>({g/ndofn_, g%ndofn_});
  }

// Constitutive
 public:
  bool planeStressElsePlaneStrain_;
  
  // {E, \nu}
  Vectord YoungPoisson_x(double x0, double x1) const {
    return Vectord({200, 0});
  }
  
  // {\lambda, \mu}
  Vectord lameConsts_x(double x0, double x1) const {
    auto yp = YoungPoisson_x(x0, x1);
    if(planeStressElsePlaneStrain_) // plane stress
      return Vectord({yp(0)*yp(1) / (1.0 - yp(1)*yp(1)), yp(0) / (2*(1.0 + yp(1)))});
    else // plane strain
      return Vectord({yp(0)*yp(1) / ((1.0 + yp(1))*(1.0 - 2.0*yp(1))), yp(0) / (2.0*(1.0 + yp(1)))});
  }
 private:
  // St. Venant-Kirchoff
  Matrix2d stressFromStrain(const Matrix2d& strain, double x0, double x1) const {
    Matrix2d stress = Matrix2d(strain.nRows(), strain.nCols());
    const auto lambdaMu = lameConsts_x(x0, x1);
    double trE = strain.trace();
    for(int i=0; i<strain.nRows(); ++i)
      for(int j=0; j<strain.nCols(); ++j)
        stress(i, j) = lambdaMu(0)*diracDelta(i, j)*trE + 2.0*lambdaMu(1)*strain(i, j);
    return stress;
  }
  
  // St. Venant-Kirchoff tensor
  Matrix4d C_VK(double x0, double x1) {
    Matrix4d m(ndofn_, ndofn_, ndofn_, ndofn_);
    auto dD = LinAlg::diracDelta; // function alias
    
    // INEFFICIENT: C_VK inherently has a lot of symmetry and I'm not using that here, so a bit inefficient no doubt; but to use it, I really need symmetry-exploiting operators (if I can even use the operators rather than pure index notation)
    for(int i=0; i<ndofn_; ++i)
      for(int j=0; j<ndofn_; ++j)
        for(int k=0; k<ndofn_; ++k)
          for(int l=0; l<ndofn_; ++l)
            m(i,j,k,l) = lameConsts_x(x0, x1)(0) * dD(i, j) * dD(k, l) +
              lameConsts_x(x0, x1)(1) * (dD(i, k) * dD(j, l) + dD(i, l) * dD(j, k));
    return m;
  }
  
// Kinematic
 private:
  //strainMeasure
  
  enum strainMeasure {
    LINEAR,
    GREEN_LAGRANGE,
    EULER_ALMANSI,
    LOGARITHMIC
  };
  
  strainMeasure strainMeasure_ = LINEAR; // support linear first, then others
  
  // J_lj = \frac{\del x_j}{\del \xi_l}
  Matrix2d jacobian(double xi0, double xi1) {
    Matrix2d J(ndofn_, ndofn_);
    for(int l=0; l<ndofn_; ++l)
      for(int j=0; j<ndofn_; ++j) {
        double sumK = 0;
        for(int k=0; k<nnode_; ++k)
          sumK += gradL_shFct_xi(xi0, xi1)(k,l) * X_0_(nodeldof2dof(k, l));//#j or l in last l?
        J(l,j) = sumK;
      }
  }
  // (J^{-1})_lj = \frac{\del \xi_l}{\del x_j}
  Matrix2d invJacobian(double xi0, double xi1) {
    return LinAlg::invertMat2d(jacobian(xi0, xi1));
  }

// ctor(s)
 public:
  Tri3(Array<int> nodes, Vectord X_0, bool planeStressElsePlaneStrain)
  : globalNodeIds_(nodes), X_0_(X_0), planeStressElsePlaneStrain_(planeStressElsePlaneStrain) {}
  
// Shape functions
 private:
  ///static double lagrangePolynomialLin_xi(double xi, int i/*, int j for dim>1*/) {return (1.0/2) * (1.0 + Xi_[i]*xi);}
  
 private:
  static constexpr int shFct_xi_size_ = nnode_; // this and similar are ultimately unnecessary but I keep for now
  static inline Vectord shFct_xi(double xi0, double xi1) {
    return Vectord({1 - xi0 - xi1, xi0, xi1});
  }
  static constexpr int shFctMatrixForm_nRows_ = ndofn_;
  static constexpr int shFctMatrixForm_nCols_ = ndof_;
  static inline Matrix2d shFctMatrixForm_xi(double xi0, double xi1) {
    Matrix2d mat(ndofn_, ndof_);
    auto shFct = shFct_xi(xi0, xi1);
    for(int i=0; i<ndofn_; ++i)
      for(int j=0; j<ndof_; ++j) {
        int nodeBlock = j/ndofn_; // floor
        mat(i, j) = diracDelta(i, j%ndofn_) * shFct(nodeBlock);
      }
    return mat;
  }
  // const because line shFcts
  static inline Matrix2d gradL_shFct_xi(double xi0 = -2, double xi = -2) {
    // \partial N_i / (\partial xi_j) (each row is grad(N_i) where i const)
    Matrix2d mat(3/*nnode_*/, 2/*ndofn_*/,
    {
      -1, -1, // grad shFct(0)
      1, 0, // grad shFct(1)
      0, 1 // grad shFct(2)
    }
    );
    return mat;
  }
  
  // row = node, col = node-local dof
  inline Matrix2d X_0_asMat() const {
    Matrix2d mat(nnode_, ndofn_);
    for(int i=0; i<nnode_; ++i)
      for(int j=0; j<ndofn_; ++j)
        mat(i, j) = X_0_(i*ndofn_+j);
    return mat;
  }
  
  // Jacobian, with right hand gradient I guess
  inline Matrix2d gradR_x_xi(double xi0, double xi1) {
    
    Matrix2d mat(ndofn_, ndofn_);
    auto gradL_shFct = gradL_shFct_xi(xi0, xi1);
    mat = mat2dTimesMat2d(transposeMat2d(gradL_shFct), X_0_asMat());
    /*for(int i=0; i<ndofn_; ++i)
      for(int j=0; j<ndofn_; ++j)
        for(int k=0; k<ndofn_; ++k)
          mat(i, j) += gradL_shFct_xi(i, k) * deriv*/
  }
  
  inline Matrix2d gradL_shFct_x(double x0, double x1) {
    mat2dTimesMat2d(invertMat2d(gradR_x_xi(x0, x1)), gradL_shFct_xi(x0, x1));
  }
  
  // Q and q will be my in general "some quantity"
  ///inline double q_xi(const Vectord& Q_xi_, double xi) {return vectDotProduct(shFct_xi(xi), Q_xi_);}
  
  ///inline static double deriv(double (*f)(double), double xix, double h) {return (f(xix+h)-f(xix-h))/2.0;}
  
  
  static double integrateTrapz(double (*f)(double), double l, double r, double n) {
    double h = (r - l) / n;
    double result = 0;
    for(int i=0; i<n-1; ++i) {
      result += (f(i*h) + f((i+1)*h))/2.0 * h;
      //db::pr("result ="+std::to_string(result));
    }
    return result;
  }
  template<typename Fct>
  static double integrateTrapzTemplated(Fct f, double l, double r, double n) {
    double h = (r - l) / n;
    double result = 0;
    for(int i=0; i<n-1; ++i) {
      result += (f(i*h) + f((i+1)*h))/2.0 * h;
    }
    return result;
  }
  //\del E : S | For engineering/infinitesimal/linear strains: E = (\Grad u + (\Grad u)^T)/2 == makeMatSymmetric(\Grad u)
  //\del (1+F) : S(E)
  //(\del du/dx)*Cvk*(du/dx) | u = ND, N is a horizontally-long matrix, D is a vertical vector
  //\del D^T dN^T/dx*Cvk* dN/dx D | dN/dx = dN/dxi * dxi/dx = (+- 1/2) * 1/(dx/dxi); x = NX -> dx/dxi = dN/xi * X = (+- 1/2)*X --> dN/dx = (+- 1/2) * 1/((+- 1/2)*X)
  //           ^now a vertical vector                                                                               "  ^=:J   "                          "  ^=J    "
  // Integrating over xi: dx = J*dxi
  // --> \del D^T int_xi(dN^T/dx*Cvk* dN/dx * J dxi) * D
  // = \del D^T int_xi((deriv_shFct_xi^T*{J^-1}) * (deriv_shFct_xi*{J^-1}) * J dxi) * D * Cvk
  // = \del D^T int_xi(deriv_shFct_xi^T * deriv_shFct_xi dxi) * Jinv * D * Cvk // Jinv const so can be taken out of integral
  //                                    ^dyadic product because vertical vector times horizontal vector
  // Also, for Cvk(x) to Cvk(xi), we need x(xi). We have x(xi), it is just N(xi)*X
  Matrix2d integrandKmat(double xi0, double xi1) {
    Matrix2d Jinv = invJacobian(xi0, xi1);
    double detJ = detMat2d(Jinv);
  }
 public: //tmp
  Matrix2d Kmat() {
    Matrix2d result(ndof_, ndof_);
    for(int i=0;i<ndof_; ++i)
      for(int j=0;j<ndof_; ++j) {
        auto integrandKmatij = [this, i, j](double xi) {return integrandKmat(xi)(i, j);};
        result(i, j) = integrateTrapzTemplated(integrandKmatij, Xi_[0], Xi_[1], 500);
      }
    return result;
  };
  static void test() {
    std::cout<<"\n";
  }

 public:
  // This is the local (pos) to global (val) mapping
  Array<int> globalNodeIds_;
  Vectord X_0_;
}; // class Tri3

} // namespace Element