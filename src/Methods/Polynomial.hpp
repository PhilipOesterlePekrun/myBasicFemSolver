#pragma once
#include <Global.hpp>

#include <myUtils.hpp>

#include <LinAlg.hpp>

// basic numerical methods
namespace MyFem::NumMethods {

using namespace LinAlg;

// returns [smallerX, largerX], real solution of ax^2 + bx + c = 0; throws error if any imaginary solution
inline Vectord solveScalarQuadraticEq(double a, double b, double c) {
  if(a==0)
    return Vectord({-c/b, -c/b});
  double underSqrt = b*b - 4*a*c;
  if(underSqrt < 0) {
    db::throwAndExit("");
    return Vectord();
  }
  else if(underSqrt == 0) return Vectord({-b/(2*a), -b/(2*a)});
  //else if(underSqurt > 0)
    if(a>0) return Vectord({(-b + sqrt(underSqrt))/(2*a), (-b - sqrt(underSqrt))/(2*a)});
    if(a<0) return Vectord({(-b - sqrt(underSqrt))/(2*a), (-b + sqrt(underSqrt))/(2*a)});
    
  return Vectord(); //# shouldnt ever get reached I think
}

} // namespace MyFem::NumMethods