#pragma once

#include <Global.hpp>

#include <myUtils.hpp>
#include <SFML/Graphics.hpp>

#include "VisualizationBase.hpp"

namespace MyFem {

namespace Vis::Objects {

// A graph object is more flexible than just a function. But the "Visualization" class is rather more like for a window with the whole UI, so graph will not derive from Visualization, instead it will reference the Visualization object; i.e., it is basically attached to a Visualization object but it controls itself so its like an autonomous object attached to the Visualization object.
class Graph {
public:
  Graph() = delete;
  Graph(VisualizationBase& vis)
    : vis_(vis), visWindow_(vis.renderWindow_) {};
  Graph(VisualizationBase& vis, uint width, uint height, uint posHorz, uint posVert)
    : vis_(vis), visWindow_(vis.renderWindow_), width_(width), height_(height), posHorz_(posHorz), posVert_(posVert) {};  

  VisualizationBase& vis_;
  sf::RenderWindow& visWindow_; // Like alias
  uint width_;
  uint height_;
  uint posHorz_;
  uint posVert_;
  
  std::vector<std::pair<double, double>> datXY;
  
  ///void operator()() const {}
  
  void draw();
  
};

} // namespace Vis::Objects

} // namespace MyFem
