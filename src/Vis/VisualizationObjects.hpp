#pragma once
#include <Global.hpp>

#include <myUtils.hpp>
#include <SFML/Graphics.hpp>

#include <MyFem_Array_decl.hpp>

namespace MyFem::Vis {
  
class VisualizationBase; // forward decl

class Object {
  friend class VisualizationBase;
 protected:
  Object() {};
  Object(uint posX, uint posY)
    : posX_(posX), posY_(posY) {};
    
  VisualizationBase* vis_;
  sf::RenderWindow* visWindow_; // Like alias
  uint posX_;
  uint posY_;
  
 private:
  bool attached_ = false; // == "activated"
 protected:
  void attach(VisualizationBase* vis);
  void deattach();
 public:
  inline bool isAttached() {
    return attached_;
  }
  void drawFrame() {
    if(attached_)
      draw();
  }
  
 protected:
  virtual void draw() const = 0;
};

// A graph object is more flexible than just a function. But the "Visualization" class is rather more like for a window with the whole UI, so graph will not derive from Visualization, instead it will reference the Visualization object; i.e., it is basically attached to a Visualization object but it controls itself so its like an autonomous object attached to the Visualization object.
class Graph : public Object {
 public:
  Graph();
  Graph(uint posX, uint posY, uint width, uint height, Array<Array<double>> graphDataXY)
    : Object(posX, posY), width_(width), height_(height), graphDataXY_(graphDataXY) {};

  uint width_;
  uint height_;
  
  Array<Array<double>> graphDataXY_; // inner array has size 2, so X, Y
  
 protected:
  virtual void draw() const override;
  
};

} // namespace MyFem::Vis
