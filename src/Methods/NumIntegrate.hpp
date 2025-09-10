#pragma once
#include <Global.hpp>

// basic numerical methods
namespace MyFem::NumMethods {

template<typename Fct>
static double integrateTrapz(Fct f, double l, double r, double n) {
  double h = (r - l) / n;
  double result = 0;
  for(int i=0; i<n-1; ++i) {
    result += (f(i*h) + f((i+1)*h))/2.0 * h;
  }
  return result;
}

namespace { int w = 2; }
template<typename Fct>
static double integrateGaussianQuadrature(Fct f, double l, double r, double highestOrder) {
  double result = 0;
  // TODO
  return result;
}

} // namespace MyFem::NumMethods
