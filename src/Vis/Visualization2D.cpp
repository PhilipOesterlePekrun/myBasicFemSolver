#include "Visualization2D.hpp"

#include <MyFem_Array_def.hpp>

namespace MyFem::Vis {
  
void Visualization2D::drawFrameImplementation() {
  //if(!paused_) { // Something like this could be an option for efficiency somehow eventually (dont redraw when paused), but its literally bottom of priority list and probably not even worth it at all
    auto elements = problem_->getElements();
    int ndofn = problem_->get_ndofn();
    
    if(currentFrame_==0) {
      auto X_0 = problem_->getX_0();
      
      FOR(e, elements.size()) {
        auto ele = elements(e);
        sf::ConvexShape tri_t;
        tri_t.setFillColor(sf::Color(100,100,100));
        tri_t.setPointCount(3);
        tri_t.setPoint(0, offsetScaledVectXY(X_0(ele->getGlobalDofIds()(0)), X_0(ele->getGlobalDofIds()(1))));
        tri_t.setPoint(1, offsetScaledVectXY(X_0(ele->getGlobalDofIds()(2)), X_0(ele->getGlobalDofIds()(3))));
        tri_t.setPoint(2, offsetScaledVectXY(X_0(ele->getGlobalDofIds()(4)), X_0(ele->getGlobalDofIds()(5))));
        
        tri_t.setOutlineThickness(4.f);
        tri_t.setOutlineColor(sf::Color(20,20,20));
        
        renderWindow_.draw(tri_t);
      }
    }
    else {
    
      auto X_t = problem_->getX_t();
      
      
      // Deformed configuration (element triangles)
      FOR(e, elements.size()) {
        auto ele = elements(e);
        sf::ConvexShape tri_t;
        tri_t.setFillColor(sf::Color(100,100,100));
        tri_t.setPointCount(3);
        tri_t.setPoint(0, offsetScaledVectXY(X_t(ele->getGlobalDofIds()(0)), X_t(ele->getGlobalDofIds()(1))));
        tri_t.setPoint(1, offsetScaledVectXY(X_t(ele->getGlobalDofIds()(2)), X_t(ele->getGlobalDofIds()(3))));
        tri_t.setPoint(2, offsetScaledVectXY(X_t(ele->getGlobalDofIds()(4)), X_t(ele->getGlobalDofIds()(5))));
        
        tri_t.setOutlineThickness(4.f);
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
  
}
  
}