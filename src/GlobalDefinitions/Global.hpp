#pragma once

#include <string>
#include <iostream>

namespace db {

#define DEBUG_MODE true

// debug print
inline void pr(std::string s = "") {
  if(DEBUG_MODE)
    std::cout<<"-- DEBUG --\n"<<s<<"\n";
}

// throw an error
inline void throwAndExit(std::string msg = "") {
  std::cout<<"\n"<<msg<<"\n";
  std::exit(1);
}

}