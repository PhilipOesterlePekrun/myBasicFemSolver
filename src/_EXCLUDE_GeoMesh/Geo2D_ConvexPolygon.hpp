#pragma once
#include "Global.hpp"

#include "mu.hpp"
#include "mu_core_LinAlg.hpp"

// geometrical object of type convex polygon
namespace MyFem::Geo2D {

using namespace MyUtils;
using namespace LinAlg;

class ConvexPolygon { // this is more like a struct but whatever
 public:
  // the order of the points should be the order of their connections. It is not necessary to have a CCW or CW order, both will work (for now, maybe later CCW will be necessary)
  ConvexPolygon(std::vector<Vectord> pointsXY_) // maybe some checks here for convexity? Or maybe do that in the mesh algo idk
  : pointsXY(pointsXY_) {}
  
 public:
  std::vector<Vectord> pointsXY;
};

} // namespace MyFem::Geo2D
