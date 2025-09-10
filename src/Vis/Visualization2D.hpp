#pragma once
#include <Global.hpp>

#include <myUtils.hpp>

#include "VisualizationBase.hpp"
#include "VisualizationObjects.hpp"
#include <Problem/Linear2D_Manager.hpp>

namespace MyFem {
  
namespace Vis {
  
class Visualization2D : public VisualizationBase {

 private:
  std::shared_ptr<Problem::Linear2D> problem_;
 public:
  double scaleFactor_;
  Array<float> XYoffset_ = Array<float>({200, 200});
  float outlineThickness_;

 public:
  Visualization2D(std::string title, int windowWidth, int windowHeight, int framerate, sf::Color baseColor, sf::Color secondaryColor, sf::Color defaultTextColor, sf::Font* defaultFont, int defaultFontSize,
  double scaleFactor, std::shared_ptr<Problem::Linear2D> problem, float outlineThickness)
  : VisualizationBase(title, windowWidth, windowHeight, framerate, baseColor, secondaryColor, defaultTextColor, defaultFont, defaultFontSize),
  scaleFactor_(scaleFactor), problem_(problem), outlineThickness_(outlineThickness) {}
  
  void setProblem(std::shared_ptr<Problem::Linear2D> problem) {
    problem_ = problem;
  }
  
 private:
  inline sf::Vector2f offsetScaledVectXY(float posX, float posY){
		return Vector2fInXY(XYoffset_(0)+scaleFactor_*posX, XYoffset_(1)+scaleFactor_*posY);
	}
  void drawFrameImplementation() override;
};

} // namespace Vis

} // namespace MyFem
