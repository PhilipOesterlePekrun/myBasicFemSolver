#pragma once
#include <Global.hpp>

#include <myUtils.hpp>

namespace Linear1D {

class Linear1DProblem {
public:
  Linear1DProblem() {};
  
  void runNoInputExample();
  
  
private:
  void readMeshTxt(std::string inputFilePath);
}; // class Linear1DProblem

} // namespace Linear1D
