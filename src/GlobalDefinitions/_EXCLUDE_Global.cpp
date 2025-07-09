#include "Global.hpp"

#include <MyFem_Array_decl.hpp>

namespace MyFem {
  
Array<std::string> strToStrArray(const std::string& s) {
  auto out = Array<std::string>();
  if(s == "") return out;
  out.push_back("");
  size_t sLength = s.length();
  FOR(i, sLength) {
    if(i!=sLength-1 && s[i]=='\n')
      out.push_back("");
    else
      out(out.size()) += s[i];
  }
  return out;
}
std::string strArrayToStr(const Array<std::string>& a) {
  std::string out = "";
  FOR(i, a.size())
    out += a(i) + "\n";
  return out;
}

std::string alignStringAt(const std::string& s, std::string& alignerKey) {
  auto keyPositionsInLine = Array<Array<size_t>>();
  FOR(i, s.length()) {
    if(s[i] == '\n')
    ;
  }
  return "";
}

}
