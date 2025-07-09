#include "VisualizationBase.hpp"

#include <MyFem_Array_def.hpp>

namespace MyFem::Vis {

void VisualizationBase::drawBaseUI() {
	// pause bar and symbol
	{
		sf::RectangleShape playPauseLoopBar(sf::Vector2f(200, 50));
		playPauseLoopBar.setPosition(Vector2fInXY((float)windowWidth_/2-100, 50));
		playPauseLoopBar.setFillColor(sf::Color(100, 100, 100));
		renderWindow_.draw(playPauseLoopBar);
		if(paused_)
		{
			sf::RectangleShape pauseL(sf::Vector2f(4, 40));
			pauseL.setPosition(Vector2fInXY(40+(float)windowWidth_/2-10, 40));
			pauseL.setFillColor(sf::Color(0,0,0));
			renderWindow_.draw(pauseL);
			sf::RectangleShape pauseR(sf::Vector2f(4, 40));
			pauseR.setPosition(Vector2fInXY(40+(float)windowWidth_/2+10, 40));
			pauseR.setFillColor(sf::Color(0,0,0));
			renderWindow_.draw(pauseR); // TODO: make better
		}
		else
		{
			sf::ConvexShape playTriangle;
			playTriangle.setFillColor(sf::Color(0,0,0));
			playTriangle.setPointCount(3);
			playTriangle.setPoint(0, Vector2fInXY(40+(float)windowWidth_/2+10, 25));
			playTriangle.setPoint(1, Vector2fInXY(40+(float)windowWidth_/2-10, 5));
			playTriangle.setPoint(2, Vector2fInXY(40+(float)windowWidth_/2-10, 50-5));
			renderWindow_.draw(playTriangle);
		}
	}

	// window time and frame count
	{
		std::string descString="Window:";
		std::string frameCountString="current frame = "+std::to_string(currentFrame_);
		std::string timeString="t = "+std::to_string(currentTime().asSeconds())+"s";
		sf::Text descText=textConstructorXY(descString, (float)windowWidth_/2-100, 50, 14);
		sf::Text frameCountText=textConstructorXY(frameCountString, (float)windowWidth_/2-90, 50-12, 12);
		sf::Text timeText=textConstructorXY(timeString, (float)windowWidth_/2-90, 50-25, 12);
		renderWindow_.draw(descText);
		renderWindow_.draw(frameCountText);
		renderWindow_.draw(timeText);
	}
}

// NOTE: Already in my XY coordinate system
sf::Text VisualizationBase::textConstructorXY(const std::string& textString, float posX, float posY, int fontSize, sf::Color fillColor, sf::Font* font) {
	sf::Text textObj(*font);
	textObj.setString(textString);
	textObj.setCharacterSize(fontSize);
	textObj.setFillColor(fillColor);
	textObj.setPosition(Vector2fInXY(posX, posY));
	return textObj;
}


// CONSTRUCTORS
VisualizationBase::VisualizationBase(std::string title, int windowWidth, int windowHeight, int framerate, sf::Color baseColor, sf::Color secondaryColor, sf::Color defaultTextColor, sf::Font* defaultFont, int defaultFontSize)
  : windowName_(title) {
	this->windowWidth_=windowWidth;
	this->windowHeight_=windowHeight;
	this->framerate_=framerate;
	this->baseColor_=baseColor;
	this->secondaryColor_=secondaryColor;
	this->defaultTextColor_=defaultTextColor; // TODO: use initializer list, idk why i did it with this-> first
  if(!defaultFont) {
    std::cerr << "ERROR: Null font pointer passed to Visualization\n";
    std::terminate();
  }
	this->defaultFont_=defaultFont;
  //this->defaultFont = new sf::Font("/home/oesterle/misc/myBasicFemSolver_Base/myBasicFemSolver/data//fonts/times.ttf");
	this->defaultFontSize_=defaultFontSize;
}
VisualizationBase::VisualizationBase(std::string title, int windowWidth, int windowHeight, int framerate, sf::Color baseColor, sf::Color secondaryColor, sf::Color defaultTextColor, sf::Font *defaultFont, int defaultFontSize, int antiAliasingLevel)
  : VisualizationBase(title, windowWidth, windowHeight, framerate, baseColor, secondaryColor, defaultTextColor, defaultFont, defaultFontSize)
{
	this->antiAliasingLevel_=antiAliasingLevel;
}

bool VisualizationBase::activate() {
	if(renderWindow_.isOpen())
    return 0;
    
	renderWindow_.create(sf::VideoMode({windowWidth_, windowHeight_}), windowName_);
	renderWindow_.setPosition(renderWindowPos_);
	renderWindow_.setFramerateLimit(framerate_);
  active_=true;
	return 1;
}

bool VisualizationBase::deactivate() {
	active_=false;
	renderWindow_.close();
	return 1;
}

bool VisualizationBase::play() {
	if(!renderWindow_.isOpen())
    return 0;

	paused_=false;
	return 1;
}

bool VisualizationBase::pause() { // TODO*: uhh do this toggle or something or add unpause() as well?
	if(!renderWindow_.isOpen())
    return 0;

	paused_=true;
	return 1;
}

bool VisualizationBase::goToFrame(int toFrame) {
	if(!renderWindow_.isOpen())
    return 0;

	currentFrame_=toFrame;
	return 1;
}

// go to the nearest time
bool VisualizationBase::goToTime(double toTime) {
	if(!renderWindow_.isOpen()||toTime>maxTime_.asSeconds())
    return 0;

	currentFrame_=toTime*framerate_;
	return 1;
}

void VisualizationBase::drawFrame() {
	if(renderWindow_.isOpen()) {
    while (const std::optional event = renderWindow_.pollEvent())
      if(event->is<sf::Event::Closed>())
      renderWindow_.close();

		if(sf::Keyboard::isKeyPressed(sf::Keyboard::Key::P)) {
			if(!pDown)
        paused_=!paused_;
			pDown=true;
		}
		else // TODO*: wtf is this? Does toggle not work without this else? Or with a new pause Toggle function that I can make....? Or how do I do this?s
      pDown=false;
      
    if(!paused_)
      currentFrame_++;
    
    if(sf::Keyboard::isKeyPressed(sf::Keyboard::Key::R)) {
			goToFrame(0);
		}

		renderWindow_.clear(baseColor_);
    drawBaseUI();
    FOR(i, objects_.size()) {
      objects_(i)->draw();
    }
		drawFrameImplementation();
		renderWindow_.display();
    
    if(sf::Keyboard::isKeyPressed(sf::Keyboard::Key::Q)) {
			deactivate();
		}
	}
}

} // namespace MyFem
