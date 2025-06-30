#include <Global.hpp>

#include <iostream>

#include <myUtils.hpp>
#include <SFML/Graphics.hpp>

#include "Vis/VisualizationBase.hpp"

using namespace MyFem;

int main1(int argCount, char** args) {
  db::pr("Test Global.hpp");
  sf::Font* timesNewRoman = new sf::Font("/home/oesterle/misc/myBasicFemSolver_Base/myBasicFemSolver/data//fonts/times.ttf");
  std::cout<<timesNewRoman->getInfo().family<<"\n";
  Vis::VisualizationBase v = Vis::VisualizationBase(2000, 1000, 60, sf::Color(120,120,120), sf::Color(200,200,0), sf::Color(0,0,0), timesNewRoman, 12);
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
#include "Problem/Linear2D_Manager.hpp"
#include "Vis/Visualization2D.hpp"

int main(int argCount, char** args) {
  std::cout<<"Main start\n";
  //Linear1D::Linear1DProblem linProb = Linear1D::Linear1DProblem();
  //std::string dataDirPathAbs = "/home/oesterle/misc/myBasicFemSolver_Base/myBasicFemSolver/data/";
  //linProb.readMeshTxt(dataDirPathAbs+"Mesh/Linear1D/SimpleBar_linA_linC_3nodes.txt");
  
  ///Problem::Linear1D p;
  ///p.runNoInputExample();
  
  Problem::Linear2D p2;
  //p2.runNoInputExample_SingleEle();
  p2.example_beam(2.0, 2.0, 20, 20);
  
  sf::Font* timesNewRoman = new sf::Font("/home/oesterle/misc/myBasicFemSolver_Base/myBasicFemSolver/data//fonts/times.ttf");
  std::cout<<timesNewRoman->getInfo().family<<"\n";
  
  std::shared_ptr<Problem::Linear2D> p2Shared(&p2);
  
  
  p2.getX_t().print();
  
  Vis::Visualization2D vis(2000, 1500, 2, sf::Color(200,200,200), sf::Color(200,200,0), sf::Color(0,0,0), timesNewRoman, 12,
  400, p2Shared);
  
  vis.activate();
  while(vis.active_) vis.drawFrame();
  
  std::cout<<"Main end\n";
  return 0;
}