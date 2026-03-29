#pragma once
#include "mu_core_GlobalCore.hpp"

// // PROJECT GLOBAL MACROS
#define STATUS(msg) std::cout<<"STATUS: "<<msg<<"\n";
#define THROW(msg) throw std::runtime_error(msg)

// // MANUAL SETTING OF MACROS RELATED TO MYUTILS
#undef myUtils_DbPr_ON

namespace MyFem {
  
// // PROJECT GLOBAL ALIASES
///auto stdtostr = [](auto x) { return std::to_string(x); };
  
// // PROJECT GLOBAL USINGS

} // namespace MyFem
