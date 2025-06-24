#pragma once

#include <string>
#include <iostream>

// // MACRO DEFINITIONS
#define FOR(i, n) for(int i = 0; i < (n); ++i)

// // DEBUG FUNCTIONALITY AND UTILITIES
namespace db {

#define DEBUG_MODE true

// debug print
inline void pr(std::string s = "", int level = 0) {
  if(DEBUG_MODE) {
    std::string dbStr = "-- DEBUG --";
    if(level == 0)
      std::cout<<dbStr+"\n"<<s<<"\n";
    else {
      std::string lvlStr;
      FOR(i, level)
        lvlStr += "\t";
      std::cout<<lvlStr+dbStr+"\n";
      
      std::string tmpS = lvlStr;
      FOR(i, s.length()) {
        tmpS += s[i];
        if(s[i] == '\n')
          tmpS += lvlStr;
      }
      std::cout<<tmpS<<"\n";
    }
  }
}

// throw an error
inline void throwAndExit(std::string msg = "") {
  std::cout<<"\n"<<msg<<"\n";
  std::exit(1);
}

}
