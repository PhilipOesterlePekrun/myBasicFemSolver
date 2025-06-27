#pragma once
#include <Global.hpp>

#include <myUtils.hpp>

#include <LinAlg.hpp>

namespace MyFem {

namespace Element {
  
using namespace LinAlg;
// for now, how we do it is we give element the info of both the reference config and the displacement/current config. This is probably necessary; if not, whatever
// most of this class is actually general enough to be the 1D element base class. But not for dim>1
class Line2 {
// Characteristic of the element type
 public:
  static constexpr int nnode_ = 2;
  static constexpr int ndofn_ = 1;
  static const int ndof_ = ndofn_*nnode_;
 private:
  // 0 and to omitted because it does not matter to this quantity. these are the discete parametric coordinates of the nodes
  static constexpr double Xi_[ndof_] = {-1, 1};
  
// Properties
 public:
  static double Cvk_x(double x) {
    return 2.0+1.5*x;
  }
 private:
  double Cvk_xi(double xi) {
    return Cvk_x(vect2dDotVect2d(shFct_xi(xi), X_0_));
  }

// ctor(s)
 public:
  Line2(Array<int> nodes, Vectord X_0)
  : globalNodeIds_(nodes), X_0_(X_0) {}
  
 private:
  static double lagrangePolynomialLin_xi(double xi, int i/*, int j for dim>1*/) {
    return (1.0/2) * (1.0 + Xi_[i]*xi);
  }
  
 private:
  static inline Vectord shFct_xi(double xi) {
    return Vectord({lagrangePolynomialLin_xi(xi, 0), lagrangePolynomialLin_xi(xi, 1)});
  }
  static inline Vectord deriv_shFct_xi(double xi = -2) {
    return Vectord({-0.5, 0.5});
  }
  
  
  // Q and q will be my in general "some quantity"
  ///inline double q_xi(const Vectord& Q_xi_, double xi) {return vectDotProduct(shFct_xi(xi), Q_xi_);}
  
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
  Matrix2d integrandKmat(double xi) {
    double J = vect2dDotVect2d(deriv_shFct_xi(xi), X_0_);
    double Jinv = 1.0/J; // Jinv const in this case, can be taken out of integral
    return scaleMat2d(Cvk_xi(xi) * Jinv,
      vect2dOuterVect2d(deriv_shFct_xi(xi), deriv_shFct_xi(xi))
    );/*
    return {
      {2, 2,
        3, 4}
    };*/
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
  }
  void test() {
    std::cout<<"\n";
  }

 public:
  // This is the local (pos) to global (val) mapping
  Array<int> globalNodeIds_;
  Vectord X_0_;
}; // class Line2

} // namespace Element

} // namespace MyFem
