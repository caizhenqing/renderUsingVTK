#pragma once
#include "algFoundation.h"

struct RenderBuffer
{
	ncs::GeometryModel geoModel;

	std::vector<ncs::Color> vertexColor;  // 直接给颜色
	std::vector<ncs::Color> lineColor;
	std::vector<ncs::Color> triColor;
	std::vector<ncs::Color> quadColor;
	std::vector<ncs::Color> solidColor;

	std::vector<double> vertexValue;     // 根据值计算颜色
	std::vector<double> lineValue;
	std::vector<double> triValue;
	std::vector<double> quadValue;
	std::vector<double> solidValue;
	double minValue = 0;
	double maxValue = 100;
	int colorNum = 9;
	std::string colorBarName = "color bar";

	double opacity = 0.3;
	bool ifShowColorBar = false;
	bool ifUsingColor = true;
};
typedef std::vector<RenderBuffer> Buffers;


enum class sceneView
{
	Main,
	Up,
	Down,
	Left,
	Right,
	Suitable,
	rotX,
	rotY,
	rotZ
};

struct globalInfoSet
{
	sceneView view;
	bool autoSceneSize = true;
	double SceneSize = 10;
};

struct ModelInfoSet
{
	int partID;
	bool ifPlotOutline;
	double lineWidth;
	double pointSize;
	bool ifParllerPerspective;
	bool ifPartVisible;
};


class RenderBase
{
public:
	virtual void init(bool independent = false, double sceneSize = 100) = 0;
	virtual void* getRenderWindow() = 0;
	virtual void renderStaticFlame(RenderBuffer& renderData, bool replace = true) = 0;
	virtual void renderDynamicFlames(Buffers& datas, double intervalTime=500) = 0;
	virtual void setFlag(bool ifPlotStModel = false, bool ifPlotDynModel = false,
		int flameID = 0) = 0;
	virtual void drawGround(double gridSize=10, double sceneSize=100, 
		bool modify = false, bool ifPlot = false) = 0;
	virtual void destroy() = 0;
	virtual void update() = 0;
	virtual void start() = 0;
	virtual void setPartVisibility(int partID, bool ifVisible) = 0;
	virtual void setPartOpcityAndColor(int partID, double *opcity = nullptr, 
		ncs::Color *color = nullptr) = 0;

	virtual void setModelPara(int partID, bool ifPlotOutline, double *lineWidth, double *pSize) = 0;

	virtual void setView(sceneView view, bool autoSceneSize = true, 
		double SceneSize = 10, bool paraView = false) = 0;
};