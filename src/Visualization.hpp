#pragma once
#include <Global.hpp>

#include <myUtils.hpp>
#include <SFML/Graphics.hpp>

class Visualization {
  // UTILITIES
protected:
	inline sf::Vector2f Vector2fInXY(float posX, float posY){
		return sf::Vector2f(posX,windowHeight-posY);
	}
	sf::Text textConstructorXY(std::string textString, sf::Font *font, int fontSize, sf::Color fillColor, float posX, float posY);
  
  // CONSTRUCTORS
public:
	inline Visualization(){} // default constructor
	Visualization(int windowWidth, int windowHeight, int framerate, sf::Color baseColor, sf::Color secondaryColor, sf::Color defaultTextColor, sf::Font* defaultFont, int defaultFontSize);
	Visualization(int windowWidth, int windowHeight, int framerate, sf::Color baseColor, sf::Color secondaryColor, sf::Color defaultTextColor, sf::Font* defaultFont, int defaultFontSize, int antiAliasingLevel); // (antiAliasingLevel = 0) == off
  
  // LOGICAL FUNCTIONS
	bool activate();
	bool deactivate();
	bool play();
	bool pause();
	bool goToTime(double toTime); // go to nearest time

  // DRAWING
protected:
	void drawBaseUI(); //# should this be virtual or not idk?
protected:
	virtual void drawFrameImplementation(); // this has the logic to draw the frame, with the functions inserted into it
private:
	bool pDown=false;
public:
	void drawFrame(); // this just has all of the functions sequentially for drawing frame

// GETTERS
protected:
  inline unsigned int getMaxFrame() {
    return framerate*maxTime.asSeconds();
  }
	inline const sf::Time currentTime() {
    return sf::Time(sf::seconds((float)currentFrame/framerate));
  }
  
/*static inline void loadFont(sf::Font *fontObj, std::string filePath) {
	if(!fontObj->loadFromFile(filePath)){
		std::cout<<"ERROR: Cannot load font"<<"\n";
	}
}*/
  
// MAIN MEMBER VARIABLES
public:
	sf::Vector2i renderWindowPos=sf::Vector2i(50, 50); // default, but you can change it from externally so you can open multiple simulation windows perfectly next to eachother
	std::string windowName="Visualization";
	int windowWidth;
	int windowHeight;
	int framerate;
	int antiAliasingLevel=0;
	sf::Color baseColor; // == background color
	sf::Color secondaryColor;
	sf::Color defaultTextColor;
	sf::Font* defaultFont;
	int defaultFontSize; // might delete if it does not make sense

	unsigned int currentFrame=0; // or use time, idk
	bool active=false; // active = window open
	bool paused=true;
	bool playOnLoop=false;
protected:
	sf::RenderWindow renderWindow;
	sf::Time maxTime; // or maxFrame idk; in any case, it is taken from simulation in the constructor for everything except the Visualization base class
};
