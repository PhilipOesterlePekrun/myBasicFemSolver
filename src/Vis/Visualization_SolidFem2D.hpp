#pragma once
#include "Global.hpp"

#include "mu.hpp"

#include "VisualizationBase.hpp"
#include "VisualizationObjects.hpp"
#include <Problem/SolidFem/Linear2D_Manager.hpp>

namespace MyFem {
  
namespace Vis {
  
class Visualization2D : public VisualizationBase {

 private:
  std::shared_ptr<Problem::Linear2D> problem_;
  
  double scaleFactor_;
  std::vector<float> XYoffset_ = std::vector<float>{{200, 200}};
  float outlineThickness_;
  bool showNodeIdLabels_;

 public:
  Visualization2D(std::string title, int windowWidth, int windowHeight, int framerate, sf::Color baseColor, sf::Color secondaryColor, sf::Color defaultTextColor, sf::Font* defaultFont, int defaultFontSize,
  double scaleFactor, std::shared_ptr<Problem::Linear2D> problem, float outlineThickness, bool showNodeIdLabels = false)
  : VisualizationBase(title, windowWidth, windowHeight, framerate, baseColor, secondaryColor, defaultTextColor, defaultFont, defaultFontSize),
  scaleFactor_(scaleFactor), problem_(problem), outlineThickness_(outlineThickness), showNodeIdLabels_(showNodeIdLabels) {}
  
  void setProblem(std::shared_ptr<Problem::Linear2D> problem) {
    problem_ = problem;
  }
  
 private:
  inline sf::Vector2f offsetScaledVectXY(float posX, float posY){
		return sf::Vector2f(XYoffset_[0]+scaleFactor_*posX, XYoffset_[1]+scaleFactor_*posY);
	}
  void drawFrameImplementation() override;
};

} // namespace Vis

} // namespace MyFem
