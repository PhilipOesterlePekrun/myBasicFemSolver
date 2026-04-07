#pragma once
#include "Global.hpp"

#include "mu_core.hpp"
#include "mu_db.hpp"
#include "mu_vis.hpp"

#include <Problem/SolidFem/Linear2D_Manager.hpp>

namespace MyFem {
  
namespace Vis {
  
class Visualization2D : public MyUtils::Vis::VisualizationBase {

 private:
  std::shared_ptr<Problem::Linear2D> problem_;
  
  double scaleFactor_;
  std::vector<float> XYoffset_ = std::vector<float>{{200, 200}};
  float outlineThickness_;
  bool showNodeIdLabels_;
  double exaggerationMultiplier_;

 public:
  Visualization2D(std::string title, int windowWidth, int windowHeight, int framerate, sf::Color baseColor, sf::Color secondaryColor, sf::Color defaultTextColor, sf::Font* defaultFont, int defaultFontSize,
  double scaleFactor, std::shared_ptr<Problem::Linear2D> problem, float outlineThickness, bool showNodeIdLabels = false, double exaggerationMultiplier = 1.0)
  : VisualizationBase(title, windowWidth, windowHeight, framerate, baseColor, secondaryColor, defaultTextColor, defaultFont, defaultFontSize),
  scaleFactor_(scaleFactor), problem_(problem), outlineThickness_(outlineThickness), showNodeIdLabels_(showNodeIdLabels), exaggerationMultiplier_(exaggerationMultiplier) {
    maxFrame_ = problem_->get_timeSteps();
  }
  
 private:
  inline sf::Vector2f offsetScaledVectXY(float posX, float posY){
		return sf::Vector2f(XYoffset_[0]+scaleFactor_*posX, XYoffset_[1]+scaleFactor_*posY);
	}
  void drawFrameImplementation() override;
};

} // namespace Vis

} // namespace MyFem
