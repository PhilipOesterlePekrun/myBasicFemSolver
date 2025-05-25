#pragma once
#include <Global.hpp>

#include <myUtils.hpp>

#include <array>

namespace Element {

// for now, how we do it is we give element the info of both the reference config and the displacement/current config. This is probably necessary; if not, whatever
// most of this class is actually general enough to be the 1D element base class. But not for dim>1
class Element_line2 {

protected:
template <typename T, std::size_t N>
using array = std::array<T, N>;

template <std::size_t N>
using arrayd = array<double, N>;

template <std::size_t N>
using arrayi = array<int, N>;


  
  
// Characteristic of the problem
public:
  static constexpr int nnodee_ = 2;
  static constexpr int ndofn_ = 1;
  static const int ndof_ = ndofn_*nnodee_;
private:
  // 0 and to omitted because it does not matter to this quantity. these are the discete parametric coordinates of the nodes
  static constexpr double Xi_[ndof_] = {-1, 1};

// ctor(s)
public:
  //Element_line2(int (&nodes)[nnodee_], int (&nodeX_0)[nnodee_]);
  //Element_line2(arrayi<nnodee_> nodes, arrayd<ndof_> X_0, double (&X_t)[ndof_]);
  //Element_line2(arrayi<nnodee_> nodes, arrayd<ndof_> X_0, arrayd<ndof_>* X_t);
  Element_line2(arrayi<nnodee_> nodes, arrayd<ndof_> X_0);
  
  double X_wrt_xi(double xi) {
    
  };
  
private:
  double lagrangePolynomialLin_xi(double xi, int i/*, int j for dim>1*/) {
    // mTODO: could add checks for valid input
    return (1/2) * (1 - Xi_[i]*xi);
  };
  
private:
  // shape function wrt xi for node0
  double shFct0_xi_(double xi) {return lagrangePolynomialLin_xi(xi, 0);};
  
  double shFct1_xi_(double xi) {return lagrangePolynomialLin_xi(xi, 1);};
  
  inline arrayd<ndof_> shFct_xi_(double xi) {
    return {shFct0_xi_(xi), shFct1_xi_(xi)};
  };
  
  double deriv_shFct0_xi(double xi = -2) {};
  
  // Q and q will be my in general "some quantity"
  inline double q_xi(const arrayd<ndof_>& Q_xi_, double xi) {return shFct0_xi_(xi)*Q_xi_[0] + shFct1_xi_(xi)*Q_xi_[1];};

public:
  arrayi<nnodee_> nodeGlobalIds_;
  
// Not all of these may be necessary
  // Nodal coordinates (capital X is nodal discrete, lowercase x is continuous function but discretized) in reference configuration (0) in parametric space (xi)
  ///arrayd<ndof_> X_0_xi_ = {-1, 1};
  // if parametric is not specified, they are not parametric
  arrayd<ndof_> X_0_;
  
  
  // !the displacement in reference configuration is obviously 0 so not defined
  ///! double (*D_0_xi)[ndof_]; // these are still the old raw arrays I was using first but whatever
  ///! double (*D_0_)[ndof_];
  
  arrayd<ndof_>* D_t_xi;
  arrayd<ndof_>* D_t_;
  
  // in deformed configuration
  ///arrayd<ndof_>* X_t_xi_;
  arrayd<ndof_>* X_t_;
  
  // again, lowercase x means continuous function but discretized
  inline double x_0_xi(double xi) {return q_xi(X_0_, xi);};
  double x_0_(double x) {};
  
  ///! double d_0_xi(double xi) {};
  ///! double d_0_(double xi) {};
  
  // could rename to u from d idk
  double d_t_xi(double xi) {return q_xi(*D_t_xi, xi);};
  double d_t_(double x) {};
  
  ///double x_t_xi(double xi) {return q_xi(*X_t_, xi);}; //equiv
  double x_t_xi(double xi) {return x_0_xi(xi) + d_t_xi(xi);}; //equiv
  double x_t_(double x) {return x_0_(x) + d_t_(x);};
  
  
  // // Inverse
  double xi_t_x(double x) {return 1/x_t_xi}
  
  
  
  // // Derived quantities
  
  double F_(double d) {return 1 + deriv0_(d)};
  
  
  
  
  

  
}; // class Element_line2

} // namespace Element