#include <Global.hpp>

#include <iostream>

#include <myUtils.hpp>
#include <SFML/Graphics.hpp>

#include "Vis/Visualization.hpp"

int main2(int argCount, char** args) {
  db::pr("Test Global.hpp");
  sf::Font* timesNewRoman = new sf::Font("/home/oesterle/misc/myBasicFemSolver_Base/myBasicFemSolver/data//fonts/times.ttf");
  std::cout<<timesNewRoman->getInfo().family<<"\n";
  Visualization v = Visualization(2000, 1000, 60, sf::Color(120,120,120), sf::Color(200,200,0), sf::Color(0,0,0), timesNewRoman, 12);
  v.activate();
  while(v.active_) v.drawFrame();
  //Visualization v = Visualization();
  //v.activate();
  //while(v.active) v.drawFrame();
  std::cout<<"works2\n";
  std::string* lines = (std::string*)malloc(sizeof(std::string)*2);
  lines[0]="0000000";
  lines[1]="1111111";
  std::cout<<Utils::Strings::checkForIn("200", "1.200000",1)[0]<<"\n";
  //Utils::IO::writeFileLinesBinary("../testing.txt", lines, 2);
  
  free(lines);
  
  db::pr("Test Global.hpp");
  
  return 0;
}

#include "Problem/Linear1D_Manager.hpp"

int main(int argCount, char** args) {
  std::cout<<"Main start\n";
  //Linear1D::Linear1DProblem linProb = Linear1D::Linear1DProblem();
  //std::string dataDirPathAbs = "/home/oesterle/misc/myBasicFemSolver_Base/myBasicFemSolver/data/";
  //linProb.readMeshTxt(dataDirPathAbs+"Mesh/Linear1D/SimpleBar_linA_linC_3nodes.txt");
  
  Problem::Linear1D p;
  p.runNoInputExample();
  
  std::cout<<"Main end\n";
  return 0;
}