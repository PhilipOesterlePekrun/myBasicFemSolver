#pragma once
#include <Global.hpp>

#include <myUtils.hpp>
#include <LinAlg.hpp>

// geometrical object of type convex polygon
namespace MyFem::Geo2D {
  
using namespace LinAlg;

namespace Utils {
  
// Euclidian distance between two points
inline double dist(const Vectord& p1, const Vectord& p2) {
  const double dx = p2(0) - p1(0);
  const double dy = p2(1) - p1(1);
  return sqrt(dx*dx+dy*dy);
}

// Linear interpolation between two points; c \in [0, 1]
inline Vectord lInterp(const Vectord& p1, const Vectord& p2, double c) {
  // equiv to p1(0)+c*(p2(0) - p1(0))
  return Vectord({(1-c)*p1(0)+c*p2(0), (1-c)*p1(1)+c*p2(1)});
}
  
} // namespace Utils

} // namespace MyFem::Geo2D