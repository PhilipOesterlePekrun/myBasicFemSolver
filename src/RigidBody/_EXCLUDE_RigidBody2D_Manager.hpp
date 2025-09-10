#pragma once
#include <Global.hpp>

#include <myUtils.hpp>
#include <LinAlg.hpp>

#include "RadsymRigidBody.hpp"

namespace MyFem {

// Rigid Body
namespace RB {
  
using namespace LinAlg;

class Manager2D {
 public:
  Manager2D() {};
  
  // members
 private:
  bool wallsEnabled_ = false;
  double x0_, x1_, y0_, y1_;
  
  Array<RigidBody*> rigidBodies_;
  
 public:
  // i.e. rectangular boundary
  void enableWalls(double x0, double x1, double y0, double y1) {
    wallsEnabled_ = true;
    x0_ = x0;
    x1_ = x1;
    y0_ = y0;
    y1_ = y1;
  }
  // e.g. gravity (though doesnt have to be const wrt pos or time) //TODO: make this take a function object or something so you can set it outside idk--just do it however I did it with numerical integration
  Vectord globalAcceleration(double x, double y, double t) {
    return Vectord({0, -9.81});
  }
  void example1() {
    rigidBodies_.push_back(new)
  }
};

}

}
