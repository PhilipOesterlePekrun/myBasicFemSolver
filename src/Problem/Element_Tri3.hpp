#pragma once
#include <Global.hpp>

#include <myUtils.hpp>

#include <LinAlg.hpp>

namespace Element {
  
using namespace LinAlg;
// for now, how we do it is we give element the info of both the reference config and the displacement/current config. This is probably necessary; if not, whatever
// most of this class is actually general enough to be the 1D element base class. But not for dim>1
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
  // St. Venant - Kirchoff
  Matrixd stressFromStrain(const Matrixd& strain, double x0, double x1) const {
    Matrixd stress = Matrixd(strain.nRows(), strain.nCols());
    const auto lambdaMu = lameConsts_x(x0, x1);
    double trE = strain.trace();
    for(int i=0; i<strain.nRows(); ++i)
      for(int j=0; j<strain.nCols(); ++j)
        stress(i, j) = lambdaMu(0)*diracDelta(i, j)*trE + 2.0*lambdaMu(1)*strain(i, j);
    return stress;
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
  static inline Matrixd shFctMatrixForm_xi(double xi0, double xi1) {
    Matrixd mat(ndofn_, ndof_);
    auto shFct = shFct_xi(xi0, xi1);
    for(int i=0; i<ndofn_; ++i)
      for(int j=0; j<ndof_; ++j) {
        int nodeBlock = j/ndofn_; // floor
        mat(i, j) = diracDelta(i, j%ndofn_) * shFct(nodeBlock);
      }
    return mat;
  }
  // const because line shFcts
  static inline Matrixd gradL_shFct_xi(double xi0 = -2, double xi = -2) {
    // \partial N_i / (\partial xi_j) (each row is grad(N_i) where i const)
    Matrixd mat(nnode_, ndofn_,
    {
      -1, -1, // grad shFct(0)
      1, 0, // grad shFct(1)
      0, 1 // grad shFct(2)
    }
    );
    return mat;
  }
  
  // row = node, col = node-local dof
  inline Matrixd X_0_asMat() const {
    Matrixd mat(nnode_, ndofn_);
    for(int i=0; i<nnode_; ++i)
      for(int j=0; j<ndofn_; ++j)
        mat(i, j) = X_0_(i*ndofn_+j);
    return mat;
  }
  
  // Jacobian, with right hand gradient I guess
  inline Matrixd gradR_x_xi(double xi0, double xi1) {
    
    Matrixd mat(ndofn_, ndofn_);
    auto gradL_shFct = gradL_shFct_xi(xi0, xi1);
    mat = matTimesMat(transposeMat(gradL_shFct), X_0_asMat()),
    /*for(int i=0; i<ndofn_; ++i)
      for(int j=0; j<ndofn_; ++j)
        for(int k=0; k<ndofn_; ++k)
          mat(i, j) += gradL_shFct_xi(i, k) * deriv*/
  }
  
  inline Matrixd gradL_shFct_x(double x0, double x1) {
    matTimesMat(invertMat2(gradR_x_xi(x0, x1)), gradL_shFct_xi(x0, x1));
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
  Matrixd integrandKmat(double xi0, double xi1) {
    double J = vectDotProduct(deriv_shFct_xi(xi), X_0_);
    double Jinv = 1.0/J; // Jinv const in this case, can be taken out of integral
    return scaleMat(Cvk_xi(xi) * Jinv,
      vectDyadicProduct(deriv_shFct_xi(xi), deriv_shFct_xi(xi))
    );/*
    return {
      {2, 2,
        3, 4}
    };*/
  }
 public: //tmp
  Matrixd Kmat() {
    Matrixd result(ndof_, ndof_);
    for(int i=0;i<ndof_; ++i)
      for(int j=0;j<ndof_; ++j) {
        auto integrandKmatij = [this, i, j](double xi) {return integrandKmat(xi)(i, j);};
        result(i, j) = integrateTrapzTemplated(integrandKmatij, Xi_[0], Xi_[1], 500);
      }
    return result;
  };
  static void test() {
    std::cout<<std::endl;
  }

 public:
  // This is the local (pos) to global (val) mapping
  Array<int> globalNodeIds_;
  Vectord X_0_;
}; // class Tri3

} // namespace Element