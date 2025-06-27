#include <myUtils.hpp>

#include <LinAlg.hpp>

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
    Array<int> arr({0, 2, 4, 6, 8, 10});
    arr.print();
    arr.deleteIndices(Array<size_t>({2, 0}));
    arr.print();
  }
  
  return 0;
}