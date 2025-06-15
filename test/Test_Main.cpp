#include <myUtils.hpp>

#include <LinAlg.hpp>

int main(int argCount, char** args) {
  using namespace LinAlg;
  
  auto m = Matrix2d(2, 3,
    {
      11, 12, 13,
      21, 22, 33
    }
  );
  m.print();
  
  m.rowAt(1).print();
  m.colAt(1).print();
  
  double d;
  d++;
  std::cout<<"d="<<d<<"\n";
  std::cout<<"Main start\n";
  return 0;
}