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
  
  
  Problem::Linear2D p;
  {
    MyUtils::Timers::ScopedTimer timerLinSolve("p()");
    p.example_beam_dyn(1.0, 0.4, 20, 8, 200, 1e-8);
  }
  auto pShared = std::make_shared<Problem::Linear2D>(p);
  
  Problem::Linear2D p2;
  {
    MyUtils::Timers::ScopedTimer timerLinSolve("p2()");
    //p2.example_beam(1.0, 0.2, 20, 4, 10000, 1e-8, 0.05);
  }
  auto p2Shared = std::make_shared<Problem::Linear2D>(p2);
  
  sf::Font* timesNewRoman = new sf::Font("/home/oesterle/misc/myBasicFemSolver_Base/myBasicFemSolver/data//fonts/times.ttf");
  
  std::cout<<Timers::TimerRegistry::globalInstance().timingReportStr()<<"\n";
  
  
  //p2.getX_t().print();
  
  MyFem::Vis::Visualization2D vis("p", 1600, 1000, 1.0/p.deltaT_, sf::Color(200,200,200), sf::Color(200,200,0), sf::Color(0,0,0), timesNewRoman, 12, 1200, pShared, 2, false, 4000.0);
  //MyFem::Vis::Visualization2D vis2("p2", 1600, 1000, 20, sf::Color(200,200,200), sf::Color(200,200,0), sf::Color(0,0,0), timesNewRoman, 12, 1200, p2Shared, 2);
  
  vis.activate();
  //vis2.activate();
  while(vis.isActive()) {
    vis.drawFrame();
    //vis2.drawFrame();
  }
  
  STATUS("Main end");
  return 0;
}