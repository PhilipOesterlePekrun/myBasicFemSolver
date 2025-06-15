#include <myUtils.hpp>

#include <LinAlg.hpp>

int main(int argCount, char** args) {
  using namespace LinAlg;
  
  std::string testDelimiter = "----------------------------------------------------\n";
  
  std::cout<<"Start test main"<<"\n";
  
  std::cout<<testDelimiter;
  std::cout<<"Print Matrix2d"<<"\n";
  auto m = Matrix2d(2, 3,
    {
      11, 12, 13,
      21, 22, 33
    }
  );
  m.print();
  m.rowAt(1).print();
  m.colAt(1).print();
  
  std::cout<<testDelimiter;
  double ξ = 2.0;
  int א = 20;
  std::cout<<"Testing non-ASCII variable names:\n\tξ="<<ξ<<"\tא="<<א<<"\n";
  
  return 0;
}