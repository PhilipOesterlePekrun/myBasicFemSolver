#pragma once
#include <Global.hpp>

#include <myUtils.hpp>
#include <SFML/Graphics.hpp>

namespace MyFem {

namespace Vis {
  
class VisualizationBase {
  // UTILITIES
 protected:
	inline sf::Vector2f Vector2fInXY(float posX, float posY){
		return sf::Vector2f(posX,windowHeight_-posY);
	}
  
	sf::Text textConstructorXY(const std::string& textString, float posX, float posY, int fontSize, sf::Color fillColor, sf::Font* font);
  // Overloads
  inline sf::Text textConstructorXY(const std::string& textString, float posX, float posY, int fontSize, sf::Color fillColor) {
    return textConstructorXY(textString, posX, posY, fontSize, fillColor, defaultFont_);
  };
  inline sf::Text textConstructorXY(const std::string& textString, float posX, float posY, int fontSize) {
    return textConstructorXY(textString, posX, posY, fontSize, defaultTextColor_, defaultFont_);
  };
  inline sf::Text textConstructorXY(const std::string& textString, float posX, float posY) {
    return textConstructorXY(textString, posX, posY, defaultFontSize_, defaultTextColor_, defaultFont_);
  };
  
  // CONSTRUCTORS
 public:
	inline VisualizationBase(){} // default constructor
	VisualizationBase(int windowWidth, int windowHeight, int framerate, sf::Color baseColor, sf::Color secondaryColor, sf::Color defaultTextColor, sf::Font* defaultFont, int defaultFontSize);
	VisualizationBase(int windowWidth, int windowHeight, int framerate, sf::Color baseColor, sf::Color secondaryColor, sf::Color defaultTextColor, sf::Font* defaultFont, int defaultFontSize, int antiAliasingLevel); // (antiAliasingLevel = 0) == off
  
  // LOGICAL FUNCTIONS
	bool activate();
	bool deactivate();
	bool play();
	bool pause();
  bool goToFrame(int toFrame);
	bool goToTime(double toTime); // go to nearest time

  // DRAWING
 protected:
	void drawBaseUI(); //# should this be virtual or not idk?
 protected:
	virtual void drawFrameImplementation() {} // this has the logic to draw the frame, with the functions inserted into it
 private:
	bool pDown=false;
 public:
	void drawFrame(); // this just has all of the functions sequentially for drawing frame

// GETTERS
 protected:
  inline unsigned int getMaxFrame() {
    return framerate_*maxTime_.asSeconds();
  }
	inline const sf::Time currentTime() {
    return sf::Time(sf::seconds((float)currentFrame_/framerate_));
  }
  
/*static inline void loadFont(sf::Font *fontObj, std::string filePath) {
	if(!fontObj->loadFromFile(filePath)){
		std::cout<<"ERROR: Cannot load font"<<"\n";
	}
}*/
  
// MAIN MEMBER VARIABLES
 public:
  sf::RenderWindow renderWindow_;
	sf::Vector2i renderWindowPos_=sf::Vector2i(50, 50); // default, but you can change it from externally so you can open multiple simulation windows perfectly next to eachother
	std::string windowName_="Visualization";
	uint windowWidth_;
	uint windowHeight_;
	uint framerate_;
	uint antiAliasingLevel_=0;
	sf::Color baseColor_; // == background color
	sf::Color secondaryColor_;
	sf::Color defaultTextColor_;
	sf::Font* defaultFont_;
	int defaultFontSize_ = 12; // might delete if it does not make sense

	unsigned int currentFrame_=0; // or use time, idk
	bool active_=false; // active = window open
	bool paused_=true;
	bool playOnLoop_=false;
 protected:
	sf::Time maxTime_; // or maxFrame idk; in any case, it is taken from simulation in the constructor for everything except the Visualization base class
};

} // namespace Vis

} // namespace MyFem
