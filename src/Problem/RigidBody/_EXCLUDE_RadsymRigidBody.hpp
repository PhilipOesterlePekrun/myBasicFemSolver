#pragma once

#include "mu.hpp"
#include "mu_core_LinAlg.hpp"

#include "RadsymRigidBody.hpp"

namespace MyFem {

// Rigid Body
namespace RB {

using namespace MyUtils;
using namespace LinAlg;

class RigidBody {
 public:
  RigidBody() {};
  
 public:
  std::vector<Vectord> positionXY_t;
 public:
  updateWithForce(double Fx, double Fy) {
    positionXY_t;
  }

};
  
}

}
