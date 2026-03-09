#pragma once
#include "Global.hpp"

#include "mu.hpp"
#include "mu_core_LinAlg.hpp"

// geometrical object of type convex polygon
namespace MyFem::Geo2D {
  
using namespace MyUtils;
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

// Angle between two separate line segments, in radians
/*
e.g.
  /l2(1)
 /
/l2(0)

           l1(0)___l1(1)
returns ~ 2\pi / 3 radians (60deg)
(if you switch node l1 and l2, you would get -60deg)
*/
inline double lAngle(std::vector<Vectord/*2*/>/*2*/ l1, std::vector<Vectord/*2*/>/*2*/ l2) {
  double angle1 = ::MyUtils::Math::CommonFunctions::atan2(l1[1](1) - l1[0](1), l1[1](0) - l1[0](0));
  
  double angle2 = ::MyUtils::Math::CommonFunctions::atan2(l2[1](1) - l2[0](1), l2[1](0) - l2[0](0));
  
  return angle2 - angle1; //# correct order? should be
}

// Angle between three points (in order as two segments), in radians
/*
e.g.
   /2
  /
1/___0
returns ~ 2\pi / 3 radians (60deg)
(if you switch node 0 and 2 in order, you would get -60deg)
*/
inline double lAngle(std::vector<Vectord/*2*/>/*3*/ l) {
  double angle10 = ::MyUtils::Math::CommonFunctions::atan2(l[0](1) - l[1](1), l[0](0) - l[1](0));
  
  double angle12 = ::MyUtils::Math::CommonFunctions::atan2(l[2](1) - l[1](1), l[2](0) - l[1](0));
  
  return angle12 - angle10; //# correct order? should be
}
  
} // namespace Utils

} // namespace MyFem::Geo2D