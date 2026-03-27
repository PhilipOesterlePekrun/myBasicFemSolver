#include <iostream>

#include <SFML/Graphics.hpp>

#include "mu.hpp"
#include "mu_vis.hpp"

#include <Problem/SolidFem/Linear1D_Manager.hpp>
#include <Problem/SolidFem/Linear2D_Manager.hpp>
#include <Vis/Visualization_SolidFem2D.hpp>

namespace {
  using namespace MyUtils;
  using namespace MyFem;
}

int main(int argCount, char** args) {
  Timers::TimerRegistry::globalInstance().start();
  STATUS("Main start");
  //Linear1D::Linear1DProblem linProb = Linear1D::Linear1DProblem();
  //std::string dataDirPathAbs = "/home/oesterle/misc/myBasicFemSolver_Base/myBasicFemSolver/data/";
  //linProb.readMeshTxt(dataDirPathAbs+"Mesh/Linear1D/SimpleBar_linA_linC_3nodes.txt");
  
  ///Problem::Linear1D p;
  ///p.runNoInputExample();
  Problem::Linear2D p2;
  //p2.runNoInputExample_SingleEle();
  p2.example_beam(1.0, 0.1, 50, 5, 2000, 1e-8);
  //p2.example_beam(4.0, 0.5, 8, 2);
  //p2.example_torus(1, 1, 0.9, 1.1, 20, 2);
  
  sf::Font* timesNewRoman = new sf::Font("/home/oesterle/misc/myBasicFemSolver_Base/myBasicFemSolver/data//fonts/times.ttf");
  ///std::cout<<timesNewRoman->getInfo().family<<"\n";
  
  //std::shared_ptr<Problem::Linear2D> p2Shared(&p2);
  auto p2Shared = std::make_shared<Problem::Linear2D>(p2);
  
  std::cout<<Timers::TimerRegistry::globalInstance().timingReportStr()<<"\n";
  
  
  //p2.getX_t().print();
  
  MyFem::Vis::Visualization2D vis("2D vis", 1600, 1000, 20, sf::Color(200,200,200), sf::Color(200,200,0), sf::Color(0,0,0), timesNewRoman, 12,
  1200, p2Shared, 2);
  
  vis.activate();
  while(vis.isActive()) vis.drawFrame();
  
  STATUS("Main end");
  return 0;
}