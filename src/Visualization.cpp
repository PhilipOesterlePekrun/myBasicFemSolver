#include "Visualization.hpp"

void Visualization::drawBaseUI() {
	// pause bar and symbol
	{
		sf::RectangleShape playPauseLoopBar(sf::Vector2f(400, 100));
		playPauseLoopBar.setPosition(Vector2fInXY(windowWidth/2-200, 200));
		playPauseLoopBar.setFillColor(sf::Color(100, 100, 100));
		renderWindow.draw(playPauseLoopBar);
		if(paused)
		{
			sf::RectangleShape pauseL(sf::Vector2f(10, 80));
			pauseL.setPosition(Vector2fInXY(windowWidth/2-5-20, 200-10));
			pauseL.setFillColor(sf::Color(0,0,0));
			renderWindow.draw(pauseL);
			sf::RectangleShape pauseR(sf::Vector2f(10, 80));
			pauseR.setPosition(Vector2fInXY(windowWidth/2-5+20, 200-10));
			pauseR.setFillColor(sf::Color(0,0,0));
			renderWindow.draw(pauseR);
		}
		else
		{
			sf::ConvexShape playTriangle;
			playTriangle.setFillColor(sf::Color(0,0,0));
			playTriangle.setPointCount(3);
			playTriangle.setPoint(0, Vector2fInXY(windowWidth/2+20, 200-50));
			playTriangle.setPoint(1, Vector2fInXY(windowWidth/2-20, 200-10));
			playTriangle.setPoint(2, Vector2fInXY(windowWidth/2-20, 200-90));
			renderWindow.draw(playTriangle);
		}
	}

	// window time and frame count
	{
		std::string descString="Window:";
		std::string frameCountString="current frame = "+std::to_string(currentFrame);
		std::string timeString="t = "+std::to_string(currentTime().asSeconds())+"s";
		sf::Text descText=textConstructorXY(descString, defaultFont, 14, defaultTextColor, windowWidth/2-200, 200);
		sf::Text frameCountText=textConstructorXY(frameCountString, defaultFont, 12, defaultTextColor, windowWidth/2-190, 200-12);
		sf::Text timeText=textConstructorXY(timeString, defaultFont, 12, defaultTextColor, windowWidth/2-190, 200-24);
		renderWindow.draw(descText);
		renderWindow.draw(frameCountText);
		renderWindow.draw(timeText);
	}
}

// NOTE: Already in my XY coordinate system
sf::Text Visualization::textConstructorXY(std::string textString, sf::Font* font, int fontSize, sf::Color fillColor, float posX, float posY) {
	sf::Text textObj(*font);
	textObj.setString(textString);
	textObj.setCharacterSize(fontSize);
	textObj.setFillColor(fillColor);
	textObj.setPosition(Vector2fInXY(posX, posY));
	return textObj;
}


// CONSTRUCTORS
Visualization::Visualization(int windowWidth, int windowHeight, int framerate, sf::Color baseColor, sf::Color secondaryColor, sf::Color defaultTextColor, sf::Font* defaultFont, int defaultFontSize) {
  db::pr();
	this->windowWidth=windowWidth;
	this->windowHeight=windowHeight;
	this->framerate=framerate;
	this->baseColor=baseColor;
	this->secondaryColor=secondaryColor;
	this->defaultTextColor=defaultTextColor;
  if(!defaultFont) {
    std::cerr << "ERROR: Null font pointer passed to Visualization\n";
    std::terminate();
  }
	this->defaultFont=defaultFont;
  //this->defaultFont = new sf::Font("/home/oesterle/misc/myBasicFemSolver_Base/myBasicFemSolver/data//fonts/times.ttf");
	this->defaultFontSize=defaultFontSize;
  db::pr();

	//maxTime=sf::Time(sf::seconds(5));
  db::pr();
}
Visualization::Visualization(int windowWidth, int windowHeight, int framerate, sf::Color baseColor, sf::Color secondaryColor, sf::Color defaultTextColor, sf::Font *defaultFont, int defaultFontSize, int antiAliasingLevel)
  : Visualization(windowWidth, windowHeight, framerate, baseColor, secondaryColor, defaultTextColor, defaultFont, defaultFontSize)
{
	this->antiAliasingLevel=antiAliasingLevel;
}

bool Visualization::activate() {
	if(renderWindow.isOpen())
    return 0;
    
	renderWindow.create(sf::VideoMode({windowWidth, windowHeight}), windowName);
	renderWindow.setPosition(renderWindowPos);
	renderWindow.setFramerateLimit(framerate);
  active=true;
	std::cout<<"activate() traversed"<<"\n"; //#
	return 1;
}

bool Visualization::deactivate() {
	active=false;
	renderWindow.close();
	return 1;
}

bool Visualization::play() {
	if(!renderWindow.isOpen())
    return 0;

	paused=false;
	return 1;
}

bool Visualization::pause() {
	if(!renderWindow.isOpen())
    return 0;

	paused=true;
	return 1;
}

// go to the nearest time
bool Visualization::goToTime(double toTime) {
	if(!renderWindow.isOpen()||toTime>maxTime.asSeconds())
    return 0;

	currentFrame=toTime*framerate;
	return 1;
}

void Visualization::drawFrameImplementation() {
	drawBaseUI();
}

void Visualization::drawFrame() {
	if(renderWindow.isOpen()) {
    while (const std::optional event = renderWindow.pollEvent())
      if(event->is<sf::Event::Closed>())
      renderWindow.close();

		if(sf::Keyboard::isKeyPressed(sf::Keyboard::Key::P)) {
			if(!pDown)
        paused=!paused;
			pDown=true;
		}
		else
      pDown=false;

		if(!paused)
      currentFrame++;

		renderWindow.clear(baseColor);
		drawFrameImplementation();
		renderWindow.display();
	}
}
