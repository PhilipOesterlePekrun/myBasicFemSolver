#include <Global.hpp>

#include <myUtils.hpp>

#include <LinAlg.hpp>
#include <MyFem_Array_def.hpp>

using namespace MyFem;

int main(int argCount, char** args) {
  using namespace LinAlg;
  
  std::string testDelimiter = "----------------------------------------------------\n";
  
  std::cout<<"Start test main"<<"\n";
  
  {
    std::cout<<testDelimiter;
    std::cout<<"Print Matrix2d"<<"\n";
    auto m = Matrix2d(2, 3,
      {
        11, 12, 13,
        21, 22, 33
      }
    );
    m.print();
    std::cout<<"m.rowAt(1).print();\n";
    m.rowAt(1).print();
    std::cout<<"m.colAt(1).print();\n";
    m.colAt(1).print();
  }
  
  {
    std::cout<<testDelimiter;
    double ξ = 2.0;
    int א = 20;
    std::cout<<"Testing non-ASCII variable names:\n\tξ="<<ξ<<"\tא="<<א<<"\n";
  }
  
  {
    std::cout<<testDelimiter;
    std::cout<<"Testing LinAlg::Array\n";
    Array<int> arr({0, 2, 4, 6, 8, 10, 12});
    arr.print();
    arr.deleteIndices(Array<size_t>({2, 0}));
    arr.print();
    std::cout<<arr.toString();
  }
  
  {
    std::cout<<testDelimiter;
    std::cout<<"Testing Global.hpp string utils\n";
    Array<std::string> strArr({"0\n", "1"});
    std::cout<<"strArrayToStr(strArr):\n"<<strArrayToStr(strArr);
    std::cout<<"strToStrArray";
    std::cout<<strToStrArray(strArrayToStr(strArr)).toString();
    
    std::string testStr = "01234 |q|\n 56 | 78 | ww";
    std::cout<<alignStringAt(testStr, "|");
  }
  
  return 0;
}