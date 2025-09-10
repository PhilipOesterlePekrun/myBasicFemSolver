#pragma once
#include <Global.hpp>

#include <myUtils.hpp>
#include <LinAlg.hpp>

#include "RadsymRigidBody.hpp"

namespace MyFem {

// Rigid Body
namespace RB {
  
using namespace LinAlg;

class RigidBody {
 public:
  RigidBody() {};
  
 public:
  Array<Vectord> positionXY_t;
 public:
  updateWithForce(double Fx, double Fy) {
    positionXY_t;
  }

};
  
}

}
