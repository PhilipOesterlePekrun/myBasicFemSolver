#include "VisualizationObjects.hpp"

#include "VisualizationBase.hpp"
#include <MyFem_Array_def.hpp>

#include <SFML/Graphics.hpp>

namespace MyFem::Vis {
  
void Object::attach(VisualizationBase* vis) {
  vis_ = vis;
  visWindow_ = &vis->renderWindow_;
  attached_ = true;
}
void Object::deattach() {
  // I don't think it actually matters if the members continue to have dangling pointers
  attached_ = false;
}
  
void Graph::draw() const {
  visWindow_->draw(vis_->textConstructorXY("X-axis", posX_+(float)width_/2, posY_-10));
  sf::RectangleShape xAxis(sf::Vector2f(width_, 5));
  xAxis.setPosition(vis_->Vector2fInXY(posX_, posY_));
  xAxis.setFillColor(sf::Color(0, 0, 0));
  visWindow_->draw(xAxis);
  
  visWindow_->draw(vis_->textConstructorXY("Y-axis", posX_-50, posY_+(float)height_/2));
  sf::RectangleShape yAxis(sf::Vector2f(5, height_));
  yAxis.setPosition(vis_->Vector2fInXY(posX_, posY_+height_));
  yAxis.setFillColor(sf::Color(0, 0, 0));
  visWindow_->draw(yAxis);
  
  int n = graphDataXY_.size()-1;
  
  visWindow_->draw(vis_->textConstructorXY(std::to_string(graphDataXY_(0)(0)), posX_, posY_-10)); // x bottom
  visWindow_->draw(vis_->textConstructorXY(std::to_string(graphDataXY_(n)(0)), posX_+width_, posY_-10)); // x top
  
  visWindow_->draw(vis_->textConstructorXY(std::to_string(graphDataXY_(0)(1)), posX_-50, posY_)); // y bottom
  visWindow_->draw(vis_->textConstructorXY(std::to_string(graphDataXY_(n)(1)), posX_-50, posY_+height_)); // y top
  
  sf::VertexArray lines(sf::PrimitiveType::LineStrip, graphDataXY_.size());
  FOR(i, graphDataXY_.size()) {
    lines[i].position = vis_->Vector2fInXY(
      posX_+2+(double)width_*((graphDataXY_(i)(0)-graphDataXY_(0)(0))/(graphDataXY_(n)(0)-graphDataXY_(0)(0))),
      posY_+2+(double)height_*((graphDataXY_(i)(1)-graphDataXY_(0)(1))/(graphDataXY_(n)(1)-graphDataXY_(0)(1)))
    );
    lines[i].color = vis_->secondaryColor_;
  }
  visWindow_->draw(lines);
}

}
