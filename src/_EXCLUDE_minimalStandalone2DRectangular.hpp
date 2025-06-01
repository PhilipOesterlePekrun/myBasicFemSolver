#pragma once
#include <Global.hpp>

#include <iostream>

#include <myUtils/Strings.hpp>
#include <SFML/Graphics.hpp>

namespace minimal {
  
class minimalStandalone2DRectangular {
  public:
  
  minimalStandalone2DRectangular() {std::cout<<"works\n";}
  
  void draw() {
    if(renderWindow.isOpen())
    {
      sf::Event event;
      while (renderWindow.pollEvent(event))
      {
        if (event.type == sf::Event::Closed)
        {
          deactivate();
        }
      }

      renderWindow.clear(baseColor);
      drawFrameImplementation();
      renderWindow.display();
    }
  }
  
  
  
  sf::RenderWindow renderWindow;
  sf::Vector2i renderWindowPos=sf::Vector2i(50, 50); // default, but you can change it from externally so you can open multiple simulation windows perfectly next to eachother
  std::string windowName="Visualization Base";
  int windowWidth;
  int windowHeight;
  int antiAliasingLevel=0;
  sf::Color baseColor; // == background color
  sf::Color secondaryColor;
  sf::Color defaultTextColor;
  sf::Font defaultFont;
  int defaultFontSize; // might delete if it does not make sense

  unsigned int currentFrame=0; // or use time, idk
  bool active=false; // active = window open
  bool paused=true;
  bool playOnLoop=false;
};
  
}