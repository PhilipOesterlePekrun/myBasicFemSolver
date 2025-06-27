#include "VisualizationBase.hpp"

#include "VisualizationElements.hpp"
using namespace MyFem::VisElements;

namespace MyFem {

void VisualizationBase::drawBaseUI() {
	// pause bar and symbol
	{
		sf::RectangleShape playPauseLoopBar(sf::Vector2f(400, 100));
		playPauseLoopBar.setPosition(Vector2fInXY((float)windowWidth_/2-200, 200));
		playPauseLoopBar.setFillColor(sf::Color(100, 100, 100));
		renderWindow_.draw(playPauseLoopBar);
		if(paused_)
		{
			sf::RectangleShape pauseL(sf::Vector2f(10, 80));
			pauseL.setPosition(Vector2fInXY((float)windowWidth_/2-5-20, 200-10));
			pauseL.setFillColor(sf::Color(0,0,0));
			renderWindow_.draw(pauseL);
			sf::RectangleShape pauseR(sf::Vector2f(10, 80));
			pauseR.setPosition(Vector2fInXY((float)windowWidth_/2-5+20, 200-10));
			pauseR.setFillColor(sf::Color(0,0,0));
			renderWindow_.draw(pauseR);
		}
		else
		{
			sf::ConvexShape playTriangle;
			playTriangle.setFillColor(sf::Color(0,0,0));
			playTriangle.setPointCount(3);
			playTriangle.setPoint(0, Vector2fInXY((float)windowWidth_/2+20, 200-50));
			playTriangle.setPoint(1, Vector2fInXY((float)windowWidth_/2-20, 200-10));
			playTriangle.setPoint(2, Vector2fInXY((float)windowWidth_/2-20, 200-90));
			renderWindow_.draw(playTriangle);
		}
	}

	// window time and frame count
	{
		std::string descString="Window:";
		std::string frameCountString="current frame = "+std::to_string(currentFrame_);
		std::string timeString="t = "+std::to_string(currentTime().asSeconds())+"s";
		sf::Text descText=textConstructorXY(descString, (float)windowWidth_/2-200, 200, 14);
		sf::Text frameCountText=textConstructorXY(frameCountString, (float)windowWidth_/2-190, 200-12, 12);
		sf::Text timeText=textConstructorXY(timeString, (float)windowWidth_/2-190, 200-24, 12);
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
VisualizationBase::VisualizationBase(int windowWidth, int windowHeight, int framerate, sf::Color baseColor, sf::Color secondaryColor, sf::Color defaultTextColor, sf::Font* defaultFont, int defaultFontSize) {
  db::pr();
	this->windowWidth_=windowWidth;
	this->windowHeight_=windowHeight;
	this->framerate_=framerate;
	this->baseColor_=baseColor;
	this->secondaryColor_=secondaryColor;
	this->defaultTextColor_=defaultTextColor;
  if(!defaultFont) {
    std::cerr << "ERROR: Null font pointer passed to Visualization\n";
    std::terminate();
  }
	this->defaultFont_=defaultFont;
  //this->defaultFont = new sf::Font("/home/oesterle/misc/myBasicFemSolver_Base/myBasicFemSolver/data//fonts/times.ttf");
	this->defaultFontSize_=defaultFontSize;
  db::pr();

	//maxTime=sf::Time(sf::seconds(5));
  db::pr();
}
VisualizationBase::VisualizationBase(int windowWidth, int windowHeight, int framerate, sf::Color baseColor, sf::Color secondaryColor, sf::Color defaultTextColor, sf::Font *defaultFont, int defaultFontSize, int antiAliasingLevel)
  : VisualizationBase(windowWidth, windowHeight, framerate, baseColor, secondaryColor, defaultTextColor, defaultFont, defaultFontSize)
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
	std::cout<<"activate() traversed"<<"\n"; //#
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

bool VisualizationBase::pause() {
	if(!renderWindow_.isOpen())
    return 0;

	paused_=true;
	return 1;
}

// go to the nearest time
bool VisualizationBase::goToTime(double toTime) {
	if(!renderWindow_.isOpen()||toTime>maxTime_.asSeconds())
    return 0;

	currentFrame_=toTime*framerate_;
	return 1;
}

void VisualizationBase::drawFrameImplementation() {
	drawBaseUI();
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
		else
      pDown=false;

		if(!paused_)
      currentFrame_++;

		renderWindow_.clear(baseColor_);
		drawFrameImplementation();
		renderWindow_.display();
	}
}

} // namespace MyFem
