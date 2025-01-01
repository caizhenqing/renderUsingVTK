#pragma once
#include <math.h>
#include "renderBase.h"
#include "vtkNew.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include <vtkDataSetMapper.h>
#include <vtkOpenGLActor.h>
#include <vtkOpenGLRenderer.h>
#include <vtkNamedColors.h>
#include <vtkOutputWindow.h>
#include <vtkRenderWindow.h>
#include <vtkInteractorStyleSwitch.h>
#include <vtkProgrammableFilter.h>
#include <vtkCallbackCommand.h>
#include <vtkProperty.h>
#include <vtkCellData.h>
#include <vtkPointData.h>
#include <vtkAxesActor.h>
#include <vtkLine.h>
#include <vtkVertex.h>
#include <vtkQuad.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkPNGReader.h>
#include <vtkTexture.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkCamera.h>
#include <vtkCoordinate.h>
#include <vtkButtonWidget.h>
#include <vtkGenericOpenGLRenderWindow.h>
#include <vtkTexturedButtonRepresentation2D.h>
#include <vtkAxesActor.h>
#include <vtkInteractorStyleSwitch.h>
#include <vtkUnstructuredGrid.h>
#include <vtkDoubleArray.h>
#include <vtkTriangle.h>
#include <vtkUnsignedCharArray.h>
#include <vtkCallbackCommand.h>
#include <vtkDataSetMapper.h>
#include <vtkOpenGLActor.h>
#include <vtkCaptionActor2D.h>
#include <vtkTextProperty.h>
#include <vtkSphere.h>
#include <vtkTransform.h>
#include <vtkScalarBarActor.h>
#include <vtkColorTransferFunction.h>
#include <vtkContourFilter.h>
#include <vtkCylinderSource.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkAppendPolyData.h>
#include <vtkLineSource.h>
#include <vtkTubeFilter.h>
#include<vtkPolyLine.h>
#include <vtkParametricSpline.h>
#include <vtkMath.h>
#include <vtkSmartPointer.h>
#include <vtkPolyLine.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkArrowSource.h>
#include <vtkTransform.h>
#include <vtkTransformPolyDataFilter.h>
#include <vtkParametricFunctionSource.h>
#include <unordered_map>
#include <vtkLight.h>
#include <vtkTextActor.h>


class RenderUsingVTK: public RenderBase
{
public:
	virtual void init(bool independent = false, double sceneSize = 100) override;
	virtual void* getRenderWindow() override;
	virtual void drawGround(double gridSize = 10, double sceneSize = 100,
		bool modify = false, bool ifPlot = false) override;
	virtual void renderStaticFlame(RenderBuffer& renderData, bool replace = true) override;
	virtual void renderDynamicFlames(Buffers& datas, double intervalTime) override;
	virtual void setFlag(bool ifPlotStModel = false, bool ifPlotDynModel = false,
		int flameNum = 0) override;
	virtual void destroy() override;
	virtual void update() override;
	virtual void start() override;
	virtual void setPartVisibility(int partID, bool ifVisible) override;
	virtual void setPartOpcityAndColor(int partID, double* opcity = nullptr,
		ncs::Color* color = nullptr) override;
	virtual void setModelPara(int partID, bool ifPlotOutline,
		double* lineWidth, double* pSize) override;
	virtual void setView(sceneView view, bool autoSceneSize = true, 
		double modelSize = 10, bool paraView = false) override;

	void rendGL();
	
private:
	void rotateViewByXAxis(double angleRadians, double focalPoint[3]);
	void rotateViewByYAxis(double angleRadians, double focalPoint[3]);
	void rotateViewByZAxis(double angleRadians, double focalPoint[3]);
	

private:
	void renderBufferToUg(RenderBuffer& renderData, std::unordered_map<int, vtkUnstructuredGrid*>& ugData);
	vtkUnstructuredGrid* RenderBufferToUg(RenderBuffer& renderData);
	void createColorBar(vtkUnstructuredGrid* data, RenderBuffer& renderInfo);
	void setColorTransferFunction(RenderBuffer& renderInfo);
	void UpdateTitlePosition(vtkObject* caller, unsigned long eventId, void* callData);
	void showAxies();
	int flameID = 0;
	// double gridSize = 5;
	double sceneSize;
	bool isIndependent = false;
	void drawSpring(double point1[3], double point2[3], double color[3]);
	void DrawCurvedArrow(
		double radius,         // 弯曲的半径
		double startAngle,     // 起始角度（弧度）
		double endAngle,       // 结束角度（弧度）
		int numPoints,         // 曲线点数量
		double tubeRadius,     // 管径
		double arrowTipLength, // 箭头长度
		double arrowColor[3],
		int axis,
		double pos[3]);  // 箭头颜色 (RGB)
	void crossProduct(const double a[3], const double b[3], double result[3]);
	void drawLinearArrow(double start[3], double end[3], double color[3]);

private:
	void* renderWindow = nullptr;
	vtkOpenGLRenderer* renderer = nullptr;
	vtkGenericOpenGLRenderWindow* renWin = nullptr;
	vtkRenderWindow *renWinInde = nullptr;       // 独立窗口
	vtkNew<vtkRenderWindowInteractor> interactor;
	vtkSmartPointer<vtkDataSetMapper> staticMapper;
	vtkSmartPointer<vtkDataSetMapper> mapper;
	vtkOpenGLActor* dynamicActor = nullptr;
	vtkOpenGLActor* staticActor = nullptr;
	vtkPoints* points;		           // 存储顶点数据集
	vtkCellArray* cells;	           // 存储cell面单元
	vtkUnsignedCharArray* cellTypes;   // 标记每个面单元的类型
	// vtkUnstructuredGrid* ug;	       // 最终数据集，依据类型设置为非结构化网格模型
	std::vector<vtkUnstructuredGrid*> ugs;

	// 每个 part 对应一个数据和 actor
	std::unordered_map<int, vtkUnstructuredGrid*> partData;
	std::unordered_map<int, vtkOpenGLActor*> partActor;
	vtkUnstructuredGrid* groundData = nullptr;
	vtkOpenGLActor* groundActor = nullptr;
	
	vtkAxesActor* axes_actor = nullptr;
	bool ifPlotStatic = false;
	bool ifPlotDynamic = false;
	bool timerHasBeenGenerated = false;
	double dynamicOpacity = 1;

	// 颜色棒
	vtkSmartPointer<vtkColorTransferFunction> colorTransferFunction = nullptr;
	vtkSmartPointer<vtkContourFilter> contourFilter = nullptr;
	vtkSmartPointer<vtkDataSetMapper> contourMapper = nullptr;
	vtkSmartPointer<vtkActor> contourActor = nullptr;
	vtkSmartPointer<vtkScalarBarActor> scalarBar = nullptr;
	vtkSmartPointer<vtkTextActor> titleActor = nullptr;

private:
	bool hasBeenGenerate = false;
	bool visible = true;
};