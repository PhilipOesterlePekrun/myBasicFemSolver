#include "Visualization_SolidFem2D.hpp"

namespace MyFem::Vis {
  
void Visualization2D::drawFrameImplementation() {
  // // INPUTS
  
  //# Doesnt work correctly currently
  while (const auto event = renderWindow_.pollEvent()) {
    if(const auto* key = event->getIf<sf::Event::KeyPressed>()) {
      switch (key->code) {
        case sf::Keyboard::Key::Equal:
          exaggerationMultiplier_*=1.2;
          STATUS("exaggerationMultiplier_="+std::to_string(exaggerationMultiplier_));
          break;
        case sf::Keyboard::Key::Subtract:
          exaggerationMultiplier_/=1.2;
          STATUS("exaggerationMultiplier_="+std::to_string(exaggerationMultiplier_));
          break;
        default:
          break;
      }
    }
  }
  
  
  // // UI
  
  //if(!paused_) { // Something like this could be an option for efficiency somehow eventually (dont redraw when paused), but its literally bottom of priority list and probably not even worth it at all
  auto elements = problem_->getElements();
  int ndofn = problem_->get_ndofn();
  
  auto nodeIdLabels = std::vector<sf::Text>(0, sf::Text(*defaultFont_, "", defaultFontSize_));
  auto X_0 = problem_->get_X_0();
  
  if(currentFrame_==0) {
    
    FOR(e, elements.size()) {
      auto ele = elements[e];
      sf::ConvexShape tri_t;
      tri_t.setFillColor(sf::Color(100,100,100));
      tri_t.setPointCount(3);
      
      FOR(triNode, 3) {
        auto tmpPt = offsetScaledVectXY(X_0(ele->getGlobalDofIds()[2*triNode]), X_0(ele->getGlobalDofIds()[2*triNode+1]));
        tri_t.setPoint(triNode, Vector2fInXY(tmpPt.x, tmpPt.y));
        nodeIdLabels.push_back(textConstructorXY(std::to_string(ele->getGlobalDofIds()[triNode*2]/2), tmpPt.x, tmpPt.y));
      }
      
      tri_t.setOutlineThickness(outlineThickness_);
      tri_t.setOutlineColor(sf::Color(20,20,20));
      
      renderWindow_.draw(tri_t);
    }
  }
  else {
    ///auto X_t = problem_->get_X_t()[currentFrame()];
    auto U_t = problem_->get_U_t()[currentFrame()];
    
    // Deformed configuration (element triangles)
    FOR(e, elements.size()) {
      auto ele = elements[e];
      sf::ConvexShape tri_t;
      tri_t.setFillColor(sf::Color(100,100,100));
      tri_t.setPointCount(3);
      FOR(triNode, 3) {
        auto tmpPt = offsetScaledVectXY(
          X_0(ele->getGlobalDofIds()[2*triNode]) + exaggerationMultiplier_*U_t(ele->getGlobalDofIds()[2*triNode]),
          X_0(ele->getGlobalDofIds()[2*triNode+1]) + exaggerationMultiplier_*U_t(ele->getGlobalDofIds()[2*triNode+1])
        );
        
        tri_t.setPoint(triNode, Vector2fInXY(tmpPt.x, tmpPt.y));
        nodeIdLabels.push_back(textConstructorXY(std::to_string(ele->getGlobalDofIds()[triNode*2]/2), tmpPt.x, tmpPt.y));
      }
      
      tri_t.setOutlineThickness(outlineThickness_);
      tri_t.setOutlineColor(sf::Color(20,20,20));
      
      renderWindow_.draw(tri_t);
    }
    
    // Initial positions (dots)
    /*FOR(i, problem_->get_nnode()) {
      sf::CircleShape pt(5);
      pt.setFillColor(sf::Color(0,100,200));
      pt.setPosition(offsetScaledVectXY(X_0(ndofn*i), X_0(ndofn*i+1)));
      renderWindow_.draw(pt);
    }*/
  }
  if(showNodeIdLabels_)
    FOR(i, nodeIdLabels.size())
      renderWindow_.draw(nodeIdLabels[i]);
}
  
} // namespace MyFem::Vis
