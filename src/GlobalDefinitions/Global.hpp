#pragma once

#include <string>
#include <iostream>

// // MACRO DEFINITIONS
#define FOR(i, n) for(int i = 0; i < (n); ++i)

// // GLOBAL ALIASES
///auto stdtostr = [](auto x) { return std::to_string(x); };

namespace MyFem {

// // STRING UTILITIES INTENDED FOR EVENTUAL OUTPUT
inline std::string levelizeString(const std::string& s, int level) {
    if(level == 0)
      return s;
    // else
    std::string lvlStr;
    FOR(i, level)
      lvlStr += "\t";
    
    std::string tmpS = lvlStr;
    FOR(i, s.length()) {
      tmpS += s[i];
      if(s[i] == '\n' && i!=s.length()-1)// -1 because we don't want to turn the last "\n" into "\n\t", or else the next thing is affected
        tmpS += lvlStr;
    }
    return tmpS;
}

// // STANDARD WARNING OUTPUT
inline void warn(const std::string& s, int level = 0) {
  std::string wStr = "-- WARNING --\n"+s+"\n";
  std::cout<<levelizeString(wStr, level);
}
  
// // DEBUG FUNCTIONALITY
namespace db {

#define DEBUG_MODE true

// debug print
inline void pr(const std::string& s = "", int level = 0) { // TODO: make this shorter by using levelizeString()
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
inline void throwAndExit(const std::string& msg = "") {
  std::cout<<"\n"<<msg<<"\n";
  std::exit(1);
}

} // namespace db

} // namespace MyFem
