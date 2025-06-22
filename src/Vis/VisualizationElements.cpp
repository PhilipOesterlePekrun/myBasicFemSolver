#include "VisualizationElements.hpp"

namespace VisElements {
  
void Graph::draw() {
  visWindow_.draw(vis_.textConstructorXY("X-axis", (float)width_/2, 1));
  sf::RectangleShape xAxis(sf::Vector2f(width_-2, 2));
  xAxis.setPosition(vis_.Vector2fInXY(2, 2));
  xAxis.setFillColor(sf::Color(255, 255, 255));
  visWindow_.draw(xAxis);
  
  visWindow_.draw(vis_.textConstructorXY("Y-axis", 1, (float)height_/2));
  sf::RectangleShape yAxis(sf::Vector2f(2, height_-2));
  yAxis.setPosition(vis_.Vector2fInXY(2, 2));
  yAxis.setFillColor(sf::Color(255, 255, 255));
  visWindow_.draw(yAxis);
}

}
