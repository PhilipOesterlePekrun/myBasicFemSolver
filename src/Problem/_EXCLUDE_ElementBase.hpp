#pragma once
#include <Global.hpp>

#include <myUtils.hpp>
#include <array>

#include <LinAlg.hpp>

namespace Element {

// for now, how we do it is we give element the info of both the reference config and the displacement/current config. This is probably necessary; if not, whatever
// most of this class is actually general enough to be the 1D element base class. But not for dim>1
class Base { // TODO
// Characteristic of the element type
 public:
  static constexpr int nnodee_ = 3;
  static constexpr int ndofn_ = 2;
  static const int ndof_ = ndofn_*nnodee_;
 private:
  // 0 and to omitted because it does not matter to this quantity. these are the discete parametric coordinates of the nodes
  static constexpr double Xi_[ndof_] = {-1, 1};
  
// Properties
 public:
  static double Cvk_x(double x) {
    return 5.7;
  }
 private:
  double Cvk_xi(double xi) {
    return Cvk_x(LinAlg::Static::vectDotProduct(shFct_xi(xi), X_0_));
  }

// ctor(s)
 public:
  tri3(arrayi<nnodee_> nodes, arrayd<ndof_> X_0)
  : globalNodeIds_(nodes), X_0_(X_0) {}
  
 private:
  static double lagrangePolynomialLin_xi(double xi, int i/*, int j for dim>1*/) {
    // mTODO: could add checks for valid input
    return (1.0/2) * (1.0 + Xi_[i]*xi);
  };
  
 private:
  static inline arrayd<ndof_> shFct_xi(double xi) {
    return {lagrangePolynomialLin_xi(xi, 0), lagrangePolynomialLin_xi(xi, 1)};
  }
  static inline arrayd<ndof_> deriv_shFct_xi(double xi = -2) {
    return {-0.5, 0.5};
  }
  
  
  inline static double deriv(double (*f)(double), double xix, double h) {
    return (f(xix+h)-f(xix-h))/2.0;
  }
  
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
  //\del E : S
  //\del (1+F) *Cvk*(1+F)
  //(\del du/dx)*Cvk*(du/dx) | u = ND, N is a horizontal vector, D is a vertical vector
  //\del D^T dN^T/dx*Cvk* dN/dx D | dN/dx = dN/dxi * dxi/dx = (+- 1/2) * 1/(dx/dxi); x = NX -> dx/dxi = dN/xi * X = (+- 1/2)*X --> dN/dx = (+- 1/2) * 1/((+- 1/2)*X)
  //           ^now a vertical vector                                                                               "  ^=:J   "                          "  ^=J    "
  // Integrating over xi: dx = J*dxi
  // --> \del D^T int_xi(dN^T/dx*Cvk* dN/dx * J dxi) * D
  // = \del D^T int_xi((deriv_shFct_xi^T*{J^-1}) * (deriv_shFct_xi*{J^-1}) * J dxi) * D * Cvk
  // = \del D^T int_xi(deriv_shFct_xi^T * deriv_shFct_xi dxi) * Jinv * D * Cvk // Jinv const so can be taken out of integral
  //                                    ^dyadic product because vertical vector times horizontal vector
  // Also, for Cvk(x) to Cvk(xi), we need x(xi). We have x(xi), it is just N(xi)*X
  matrixd<ndof_, ndof_> integrandKmat(double xi) {
    double J = LinAlg::Static::vectDotProduct(deriv_shFct_xi(xi), X_0_);
    double Jinv = 1.0/J; // Jinv const in this case, can be taken out of integral
    return LinAlg::Static::scaleMat(Cvk_xi(xi) * Jinv,
      LinAlg::Static::vectDyadicProduct(deriv_shFct_xi(xi), deriv_shFct_xi(xi))
    );/*
    return {
      {2, 2,
        3, 4}
    };*/
  }
 public: //tmp
  matrixd<ndof_, ndof_> Kmat() {
    matrixd<ndof_, ndof_> result;
    for(int i=0;i<ndof_; ++i)
      for(int j=0;j<ndof_; ++j) {
        auto integrandKmatij = [this, i, j](double xi) {return integrandKmat(xi)[i][j];};
        result[i][j] = integrateTrapzTemplated(integrandKmatij, Xi_[0], Xi_[1], 500);
      }
    return result;
  }

 public:
  // This is the local (pos) to global (val) mapping
  arrayi<nnodee_> globalNodeIds_;
  arrayd<ndof_> X_0_;
}; // class tri3

} // namespace Element