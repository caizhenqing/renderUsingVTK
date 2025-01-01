#include "renderByVTK.h"

namespace {
	class TimerEventCallback : public vtkCommand
	{
	public:
		static TimerEventCallback* New()
		{
			return new TimerEventCallback;
		}
		TimerEventCallback() : widget(0)
		{
		}

		virtual void Execute(vtkObject* caller, unsigned long, void*)
		{
			widget->rendGL();
		}

		RenderUsingVTK* widget;
	};
} // namespace

void RenderUsingVTK::init(bool independent, double scSize)
{
	vtkOutputWindow::SetGlobalWarningDisplay(0);
	sceneSize = scSize;
	this->isIndependent = independent;
	if (renWin != nullptr) renWin->Delete();
	if (renderer != nullptr) renderer->Delete();

	renderer = vtkOpenGLRenderer::New();
#if 0
	renderer->SetGradientBackground(true);  // 启用渐变背景色
	// 设置渐变背景色的起始颜色和结束颜色
	renderer->SetBackground(0.31, 0.31, 0.51);  // 起始颜色（浅色）
	renderer->SetBackground2(0.949, 0.949, 0.949); // 结束颜色（深灰色）
#else
	// renderer->SetBackground(0.9, 0.9, 0.9);
	renderer->SetBackground(1, 1, 1);
	// renderer->UseShadowsOn();
	// renderer->SetAmbient(0.3);
#endif
	double position[3]{ sceneSize, sceneSize, sceneSize };
	// double position[3]{ -sceneSize, sceneSize*0.5, sceneSize };
	renderer->GetActiveCamera()->SetPosition(position);
	double focusPoint[3]{ 0, 0, 0 };

	renderer->GetActiveCamera()->SetFocalPoint(focusPoint);
	double upDirection[3]{ 0, 1, 0 };
	renderer->GetActiveCamera()->SetViewUp(upDirection);
	// renderer->GetActiveCamera()->ParallelProjectionOn();
	renderer->GetActiveCamera()->SetViewAngle(30);

	if (!independent)
	{
		renWin = vtkGenericOpenGLRenderWindow::New();
		renWin->AddRenderer(renderer);
		return;
	}

	renWinInde = vtkRenderWindow::New();
	renWinInde->SetWindowName("subWindow");
	renWinInde->AddRenderer(renderer);
	int* screenSizeX, * screenSizeY;
	screenSizeY = renWinInde->GetScreenSize();
	screenSizeX = screenSizeY++;
	renWinInde->SetSize(*screenSizeX * 0.6, *screenSizeY * 0.6);
	renWinInde->SetPosition(*screenSizeX * 0.2, *screenSizeY * 0.2);
	
	interactor->SetRenderWindow(renWinInde);
	//设置交互器操控模式
	vtkSmartPointer<vtkInteractorStyleSwitch> style = vtkSmartPointer<vtkInteractorStyleSwitch>::New();
	interactor->SetInteractorStyle(style);
	style->SetCurrentStyleToTrackballCamera();
	interactor.Get()->Initialize();
	// interactor.Get()->CreateRepeatingTimer(150);
	// renWinInde->Start();
	// interactor.Get()->Start();
	/*interactor.Get()->Initialize();
	interactor.Get()->CreateRepeatingTimer(150);*/
	// renWin->GetInteractor()->Start();
}

void* RenderUsingVTK::getRenderWindow()
{
	if (this->isIndependent) return this->renWinInde;
	return this->renWin;
}

void RenderUsingVTK::drawGround(double step, double sceneSize, bool modify, bool ifPlot)
{
	//double xOrigin[3] = { 1.0, 0.0, 0.0 };
	//double yOrigin[3] = { 0.0, 1.0, 0.0 };
	//double zOrigin[3] = { 0.0, 0.0, 1.0 };
	//double xEnd[3] = { 2.0, 0.0, 0.0 };
	//double yEnd[3] = { 0.0, 2.0, 0.0 };
	//double zEnd[3] = { 0.0, 0.0, 2.0 };
	//double xColor[3] = { 1.0, 0.0, 0.0 }; // 红色
	//double yColor[3] = { 0.0, 1.0, 0.0 }; // 绿色
	//double zColor[3] = { 0.0, 0.0, 1.0 }; // 蓝色
	//drawLinearArrow(xOrigin, xEnd, xColor);
	//drawLinearArrow(yOrigin, yEnd, yColor);
	//drawLinearArrow(zOrigin, zEnd, zColor);

	//double xPos[3] = { 1.0, 0.0, 0.0 }; // 绿色箭头
	//double yPos[3] = { 0.0, 1.0, 0.0 }; // 绿色箭头
	//double zPos[3] = { 0.0, 0.0, 1.0 }; // 绿色箭头
	//DrawCurvedArrow(0.38, 0, vtkMath::Pi(), 20, 0.02, 0.3, xColor, 0, xPos);
	//DrawCurvedArrow(0.38, 0, vtkMath::Pi(), 20, 0.02, 0.3, yColor, 1, yPos);
	//DrawCurvedArrow(0.38, 0, vtkMath::Pi(), 20, 0.02, 0.3, zColor, 2, zPos);

	//double springBegin[3] = { 0.0, 0.0, 0.0 }; // 绿色箭头
	//double springXEnd[3] = { 1.0, 0.0, 0.0 };  // 绿色箭头
	//double springYEnd[3] = { 0.0, 1.0, 0.0 };  // 绿色箭头
	//double springZEnd[3] = { 0.0, 0.0, 1.0 };  // 绿色箭头
	//drawSpring(springBegin, springXEnd, xColor);
	//drawSpring(springBegin, springYEnd, yColor);
	//drawSpring(springBegin, springZEnd, zColor);
	//return;

	if (this->axes_actor == nullptr) showAxies();
	int num = int(sceneSize / 2.0 / step);

	if (this->groundActor == nullptr || modify)
	{
		groundActor = vtkOpenGLActor::New();
		vtkSmartPointer<vtkDataSetMapper> groundMapper =
			vtkSmartPointer<vtkDataSetMapper>::New();

		if (this->groundData == nullptr)
		{
			this->groundData = vtkUnstructuredGrid::New();
		}
		else
		{
			this->groundData->Delete();
			this->groundData = vtkUnstructuredGrid::New();
		}

		vtkSmartPointer<vtkUnsignedCharArray> cellColors = vtkSmartPointer<vtkUnsignedCharArray>::New();
		cellColors->SetNumberOfComponents(3);  // RGB
		cellColors->SetName("Color");
		this->groundData->GetCellData()->SetScalars(cellColors);
		/*vtkDoubleArray* cellNormals = static_cast<vtkDoubleArray*>
				(ugModel[md->lines[i].partID]->GetCellData()->GetNormals());*/

		double boundrySize = step * num;
		vtkPoints* groundPoints = vtkPoints::New();

		vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
		groundPoints->InsertNextPoint(-boundrySize, 0, 0);
		groundPoints->InsertNextPoint(boundrySize, 0, 0);
		line->GetPointIds()->SetId(0, 0);
		line->GetPointIds()->SetId(1, 1);
		cellColors->InsertNextTuple3(255, 0, 0);
		this->groundData->InsertNextCell(VTK_LINE, line->GetPointIds());

		groundPoints->InsertNextPoint(0, 0, -boundrySize);
		groundPoints->InsertNextPoint(0, 0, boundrySize);
		line->GetPointIds()->SetId(0, 2);
		line->GetPointIds()->SetId(1, 3);
		cellColors->InsertNextTuple3(0, 0, 255);
		this->groundData->InsertNextCell(VTK_LINE, line->GetPointIds());

		int ndID = 4;
		for (int i = 1; i < num + 1; i++)
		{
			// x 方向
			groundPoints->InsertNextPoint(-boundrySize, 0, i * step);
			groundPoints->InsertNextPoint(boundrySize, 0, i * step);
			line->GetPointIds()->SetId(0, ndID);
			line->GetPointIds()->SetId(1, ndID + 1);
			cellColors->InsertNextTuple3(0, 0, 0);
			this->groundData->InsertNextCell(VTK_LINE, line->GetPointIds());
			ndID += 2;

			groundPoints->InsertNextPoint(-boundrySize, 0, -i * step);
			groundPoints->InsertNextPoint(boundrySize, 0, -i * step);
			line->GetPointIds()->SetId(0, ndID);
			line->GetPointIds()->SetId(1, ndID + 1);
			cellColors->InsertNextTuple3(0, 0, 0);
			this->groundData->InsertNextCell(VTK_LINE, line->GetPointIds());
			ndID += 2;

			// z 方向
			groundPoints->InsertNextPoint(i * step, 0, -boundrySize);
			groundPoints->InsertNextPoint(i * step, 0, boundrySize);
			line->GetPointIds()->SetId(0, ndID);
			line->GetPointIds()->SetId(1, ndID + 1);
			cellColors->InsertNextTuple3(0, 0, 0);
			this->groundData->InsertNextCell(VTK_LINE, line->GetPointIds());
			ndID += 2;

			groundPoints->InsertNextPoint(-i * step, 0, -boundrySize);
			groundPoints->InsertNextPoint(-i * step, 0, boundrySize);
			line->GetPointIds()->SetId(0, ndID);
			line->GetPointIds()->SetId(1, ndID + 1);
			cellColors->InsertNextTuple3(0, 0, 0);
			this->groundData->InsertNextCell(VTK_LINE, line->GetPointIds());
			ndID += 2;
		}
		
		groundMapper->SetScalarVisibility(true);
		groundMapper->SetScalarModeToUseCellData(); // 确保使用单元格数据
		//groundActor->GetProperty()->SetLineWidth(8);

		groundData->SetPoints(groundPoints);
		groundMapper->SetInputData(this->groundData);
		groundMapper->Update();
		groundActor->SetMapper(groundMapper);
		renderer->AddActor(groundActor);
	}

	if (this->groundActor != nullptr)
	{
		this->groundActor->SetVisibility(ifPlot);
	}

	if (this->isIndependent)
	{
		this->renWinInde->GetInteractor()->Render();
	}
	else
	{
		this->renWin->GetInteractor()->Render();
	}
}

void RenderUsingVTK::renderStaticFlame(RenderBuffer& data, bool replace)
{
	this->renderBufferToUg(data, this->partData);
	if (replace)
	{
		for (auto& pair : this->partActor)
		{
			vtkOpenGLActor* actor = pair.second;
			if (actor)
			{
				this->renderer->RemoveActor(actor);
				//actor->Delete();  // 调用VTK对象的Delete方法来释放内存
			}
		}
		this->partActor.clear();
	}

	for (int i=0; i< this->partData.size(); i++)
	{
		vtkUnstructuredGrid* dataUg = this->partData[i];
		if (this->partActor[i] == nullptr)
		{
			this->partActor[i] = vtkOpenGLActor::New();
			vtkSmartPointer<vtkDataSetMapper> dataMapper =
				vtkSmartPointer<vtkDataSetMapper>::New();
			dataMapper->SetScalarVisibility(true);
			dataMapper->SetInputData(dataUg);
			dataMapper->Update();
			this->partActor[i]->SetMapper(dataMapper);
			
			vtkProperty* property = this->partActor[i]->GetProperty();
			property->SetEdgeVisibility(true);
			// 设置透明度
			property->SetOpacity(data.opacity);

			// 设置点大小
			property->SetPointSize(8);

			// 设置材质属性以保持亮度不变
			// 环境光系数，使亮度在不同视角下更一致
			property->SetAmbient(0.3);

			// 漫反射系数，控制光照的散射效果
			property->SetDiffuse(0.6);

			// 镜面反射系数，控制高光效果
			property->SetSpecular(0.3);

			// 镜面反射指数，控制高光的集中程度
			property->SetSpecularPower(60);


			// 开启双面渲染，使模型看起来更自然
			property->SetBackfaceCulling(0);

			renderer->AddActor(this->partActor[i]);
		}
	}

	// 设置光照
	vtkSmartPointer<vtkLight> light = vtkSmartPointer<vtkLight>::New();
	light->SetPositional(1); // 设置为位置光
	light->SetPosition(renderer->GetActiveCamera()->GetPosition()); // 光源位置与相机位置相同
	light->SetFocalPoint(renderer->GetActiveCamera()->GetFocalPoint()); // 光源焦点与相机焦点相同
	light->SetColor(1.0, 1.0, 1.0); // 白色光
	light->SetIntensity(0.5); // 光照强度
	renderer->AddLight(light);

	if (data.ifShowColorBar) createColorBar(this->partData[0], data);
	if (this->isIndependent)
	{
		this->renWinInde->GetInteractor()->Render();
	}
	else
	{
		this->renWin->GetInteractor()->Render();
	}
}

void RenderUsingVTK::renderDynamicFlames(Buffers& datas, double intervalTime)
{
	ifPlotDynamic = true;
	this->ugs.clear();
	if (datas.size() == 0) return;
	this->dynamicOpacity = datas[0].opacity;
	// this->ifPlotDynamic = true;
	if (timerHasBeenGenerated == false)
	{
		TimerEventCallback* eventCallback = nullptr;
		if (this->isIndependent)
		{
			renWinInde->GetInteractor()->CreateRepeatingTimer(intervalTime);
			eventCallback = TimerEventCallback::New();
			eventCallback->widget = this;
			renWinInde->GetInteractor()->AddObserver(vtkCommand::TimerEvent, eventCallback);
		}
		else
		{
			renWin->GetInteractor()->CreateRepeatingTimer(intervalTime);
			eventCallback = TimerEventCallback::New();
			eventCallback->widget = this;
			renWin->GetInteractor()->AddObserver(vtkCommand::TimerEvent, eventCallback);
		}
		timerHasBeenGenerated = true;
	}

	
	bool createColorBar = false;
	for (auto& buffer : datas)
	{
		vtkUnstructuredGrid* ugData = RenderBufferToUg(buffer);
		if (ugData != nullptr) this->ugs.emplace_back(ugData);
		if (!createColorBar && buffer.ifShowColorBar)
		{
			this->createColorBar(ugData, datas[0]);
			createColorBar = true;
		}
	}

}

void RenderUsingVTK::setFlag(bool ifPlotStModel, bool ifPlotDynModel, int flameNum)
{
	this->ifPlotStatic = ifPlotStModel;
	this->ifPlotDynamic = ifPlotDynModel;
	this->flameID = flameNum;

	if (this->partActor.size() > 0)
	{
		if (!ifPlotStatic)
		{
			for (auto& pActor : this->partActor)
			{
				pActor.second->SetVisibility(false);
			}
		}
		else
		{
			for (auto& pActor : this->partActor)
			{
				pActor.second->SetVisibility(true);
			}
		}
	}


	if (this->staticActor != nullptr)
	{
		if (!ifPlotStatic)
		{
			this->staticActor->SetVisibility(false);
		}
		else
		{
			this->staticActor->SetVisibility(true);
		}
	}

	if (this->dynamicActor != nullptr)
	{
		if (!ifPlotDynModel)
		{
			this->dynamicActor->SetVisibility(false);
		}
		else
		{
			this->dynamicActor->SetVisibility(true);
		}
	}
	
}

void RenderUsingVTK::destroy()
{
}

void RenderUsingVTK::update()
{
	if (this->axes_actor == nullptr) showAxies();
}

void RenderUsingVTK::start()
{
	if (this->isIndependent)
	{
		renWinInde->GetInteractor()->Start();
	}
}

void RenderUsingVTK::setPartVisibility(int partID, bool ifVisible)
{
	this->partActor[partID]->SetVisibility(ifVisible);
	//this->renderer->Render();
	if (this->isIndependent)
	{
		this->renWinInde->Render();
	}
	else
	{
		this->renWin->Render();
	}
}

void RenderUsingVTK::setPartOpcityAndColor(int partID, double* opcity, ncs::Color* color)
{
	if (opcity != nullptr)
	{
		this->partActor[partID]->GetProperty()->SetOpacity(*opcity);
	}

	if (color != nullptr)
	{
		vtkMapper* maperTem= this->partActor[partID]->GetMapper();
		if (maperTem)
		{
			maperTem->SetScalarVisibility(0); // 禁用标量映射

			double partColor[3]{ color->r / 255.0, color->g / 255.0, color->b / 255.0 };
			this->partActor[partID]->GetProperty()->SetColor(partColor);
		}


	}

	if (this->isIndependent)
	{
		this->renWinInde->Render();
	}
	else
	{
		this->renWin->Render();
	}
}

void RenderUsingVTK::setModelPara(int partID, bool ifPlotOutline, 
	double* lineWidth, double* pSize)
{
	this->partActor[partID]->GetProperty()->SetEdgeVisibility(ifPlotOutline);

	if (lineWidth != nullptr)
	{
		this->partActor[partID]->GetProperty()->SetLineWidth(*lineWidth);
	}
	if (pSize != nullptr)
	{
		this->partActor[partID]->GetProperty()->SetPointSize(*lineWidth);
	}

	if (this->isIndependent)
	{
		this->renWinInde->Render();
	}
	else
	{
		this->renWin->Render();
	}

}

void RenderUsingVTK::setView(sceneView view, bool autoSceneSize, 
	double modelSize, bool paraView)
{
	double bounds[6]{0};
	if (autoSceneSize)
	{
		renderer->ComputeVisiblePropBounds(bounds);
		// 包围盒的尺寸 (对角线长度)
		modelSize = sqrt(pow(bounds[1] - bounds[0], 2) +
			pow(bounds[3] - bounds[2], 2) +
			pow(bounds[5] - bounds[4], 2));
	}
	
	auto camera = this->renderer->GetActiveCamera();
	// 获取相机视角、焦点和原始位置
	double viewAngle = camera->GetViewAngle(); // 当前视角
	double distance = modelSize / (2 * tan(vtkMath::RadiansFromDegrees(viewAngle / 2))); // 推荐距离
	double focalPoint[3];
	camera->GetFocalPoint(focalPoint);

	// 调整视角
	if (view == sceneView::Main) {
		// 主视图：从正前方向 (0, 0, Z) 观察
		// 设置相机位置，确保模型刚好填充视野
		camera->SetPosition(focalPoint[0], focalPoint[1], focalPoint[2] + distance);
		camera->SetViewUp(0, 1, 0); // Y轴正方向为上

		// 重置相机以填充视野
		renderer->ResetCameraClippingRange();
	}
	else if (view == sceneView::Up) {
		// 俯视图：从正上方向 (0, Y, 0) 观察
		camera->SetPosition(focalPoint[0], focalPoint[1] + distance, focalPoint[2]);
		camera->SetViewUp(0, 0, -1); // Z轴负方向为上
		// camera->SetParallelProjection(true);
	}
	else if (view == sceneView::Down) {
		// 仰视图：从正下方向 (0, -Y, 0) 观察
		camera->SetPosition(focalPoint[0], focalPoint[1] - distance, focalPoint[2]);
		camera->SetViewUp(0, 0, 1); // Z轴正方向为上
	}
	else if (view == sceneView::Left) {
		// 左视图：从左侧方向 (-X, 0, 0) 观察
		camera->SetPosition(focalPoint[0] - distance, focalPoint[1], focalPoint[2]);
		camera->SetViewUp(0, 1, 0); // Y轴正方向为上
	}
	else if (view == sceneView::Right) {
		// 右视图：从右侧方向 (X, 0, 0) 观察
		camera->SetPosition(focalPoint[0] + distance, focalPoint[1], focalPoint[2]);
		camera->SetViewUp(0, 1, 0); // Y轴正方向为上
	}
	else if (view == sceneView::Suitable) 
	{
		double center[3];
		// 计算模型中心点
		if (!autoSceneSize)
		{
			renderer->ComputeVisiblePropBounds(bounds);
		}
		center[0] = (bounds[0] + bounds[1]) / 2.0;
		center[1] = (bounds[2] + bounds[3]) / 2.0;
		center[2] = (bounds[4] + bounds[5]) / 2.0;
		// 自适应视图：调整相机位置和裁剪范围，使对象完全填满视口
		// 确定相机与模型中心的距离
		double distance = 0.5 * modelSize; // 可以根据需要调整距离系数
		// 设置相机位置，使其位于模型中心的一定距离处，并且稍微偏上
		camera->SetPosition(center[0], center[1], center[2] + distance);
		// 设置相机焦点为模型中心
		camera->SetFocalPoint(center[0], center[1], center[2]);
		// 设置相机的向上方向，通常是Y轴正方向
		camera->SetViewUp(0, 1, 0);
		// 调整相机的视角范围，确保模型填充视图
		renderer->ResetCameraClippingRange();
	}
	else if (view == sceneView::rotX)
	{
		double theta = -45.0; // 定义旋转角度
		rotateViewByXAxis(vtkMath::RadiansFromDegrees(theta), focalPoint);
	}
	else if (view == sceneView::rotY)
	{
		double theta = -45.0; // 定义旋转角度
		rotateViewByYAxis(vtkMath::RadiansFromDegrees(theta), focalPoint);
	}
	else if (view == sceneView::rotZ)
	{
		double theta = -45.0; // 定义旋转角度
		rotateViewByZAxis(vtkMath::RadiansFromDegrees(theta), focalPoint);
	}

	camera->SetParallelProjection(paraView);
	// 渲染更新
	this->renderer->ResetCameraClippingRange(); // 更新裁剪范围
	this->renderer->GetRenderWindow()->Render();
}

void RenderUsingVTK::rotateViewByXAxis(double angleRadians, double focalPoint[3])
{
	vtkCamera* camera = this->renderer->GetActiveCamera();

	// 获取当前相机参数
	double cameraPosition[3];
	double viewUp[3];
	camera->GetPosition(cameraPosition);
	camera->GetViewUp(viewUp);

	// 计算相机位置相对焦点的向量
	double dx = cameraPosition[0] - focalPoint[0];
	double dy = cameraPosition[1] - focalPoint[1];
	double dz = cameraPosition[2] - focalPoint[2];

	// 计算旋转矩阵
	double cosTheta = cos(angleRadians);
	double sinTheta = sin(angleRadians);

	// 旋转后的相机位置
	double newCameraPosition[3];
	newCameraPosition[0] = focalPoint[0] + dx;
	newCameraPosition[1] = focalPoint[1] + (dy * cosTheta - dz * sinTheta);
	newCameraPosition[2] = focalPoint[2] + (dy * sinTheta + dz * cosTheta);

	// 旋转后的“上方向”
	double newViewUp[3];
	newViewUp[0] = viewUp[0];
	newViewUp[1] = viewUp[1] * cosTheta - viewUp[2] * sinTheta;
	newViewUp[2] = viewUp[1] * sinTheta + viewUp[2] * cosTheta;

	// 设置新的相机参数
	camera->SetPosition(newCameraPosition);
	camera->SetViewUp(newViewUp);
}

void RenderUsingVTK::rotateViewByYAxis(double angleRadians, double focalPoint[3])
{
	vtkCamera* camera = this->renderer->GetActiveCamera();

	// 获取当前相机参数
	double cameraPosition[3];
	double viewUp[3];
	camera->GetPosition(cameraPosition);
	camera->GetViewUp(viewUp);

	// 计算相机位置相对焦点的向量
	double dx = cameraPosition[0] - focalPoint[0];
	double dy = cameraPosition[1] - focalPoint[1];
	double dz = cameraPosition[2] - focalPoint[2];

	// 计算旋转矩阵
	double cosTheta = cos(angleRadians);
	double sinTheta = sin(angleRadians);

	// 旋转后的相机位置
	double newCameraPosition[3];
	newCameraPosition[0] = focalPoint[0] + (dx * cosTheta + dz * sinTheta);
	newCameraPosition[1] = focalPoint[1] + dy;
	newCameraPosition[2] = focalPoint[2] + (-dx * sinTheta + dz * cosTheta);

	// 旋转后的“上方向”
	double newViewUp[3];
	newViewUp[0] = viewUp[0] * cosTheta + viewUp[2] * sinTheta;
	newViewUp[1] = viewUp[1];
	newViewUp[2] = -viewUp[0] * sinTheta + viewUp[2] * cosTheta;

	// 设置新的相机参数
	camera->SetPosition(newCameraPosition);
	camera->SetViewUp(newViewUp);
}

void RenderUsingVTK::rotateViewByZAxis(double angleRadians, double focalPoint[3])
{
	vtkCamera* camera = this->renderer->GetActiveCamera();
	// 获取当前相机参数
	double cameraPosition[3];
	double viewUp[3];
	camera->GetPosition(cameraPosition);
	camera->GetViewUp(viewUp);

	// 计算相机位置相对焦点的向量
	double dx = cameraPosition[0] - focalPoint[0];
	double dy = cameraPosition[1] - focalPoint[1];
	double dz = cameraPosition[2] - focalPoint[2];

	// 计算旋转矩阵
	double cosTheta = cos(angleRadians);
	double sinTheta = sin(angleRadians);

	// 旋转后的相机位置
	double newCameraPosition[3];
	newCameraPosition[0] = focalPoint[0] + (dx * cosTheta - dy * sinTheta);
	newCameraPosition[1] = focalPoint[1] + (dx * sinTheta + dy * cosTheta);
	newCameraPosition[2] = focalPoint[2] + dz;

	// 旋转后的“上方向”
	double newViewUp[3];
	newViewUp[0] = viewUp[0] * cosTheta - viewUp[1] * sinTheta;
	newViewUp[1] = viewUp[0] * sinTheta + viewUp[1] * cosTheta;
	newViewUp[2] = viewUp[2];

	// 设置新的相机参数
	camera->SetPosition(newCameraPosition);
	camera->SetViewUp(newViewUp);
}


void RenderUsingVTK::rendGL()
{
	if (this->flameID >= this->ugs.size())
	{
		return;
		this->flameID = 0;
	}
	if (!ifPlotDynamic) this->flameID = 0;
	if (this->ugs.size() == 0) return;
	vtkUnstructuredGrid* model = this->ugs[flameID];
	this->flameID++;

	if (dynamicActor == nullptr)
	{
		dynamicActor = vtkOpenGLActor::New();
		mapper = vtkSmartPointer<vtkDataSetMapper>::New();
		mapper->SetScalarVisibility(true);
		mapper->SetInputData(model);
		mapper->Update();
		dynamicActor->SetMapper(mapper);
		dynamicActor->GetProperty()->SetEdgeVisibility(true);
		dynamicActor->GetProperty()->SetOpacity(this->dynamicOpacity);
		dynamicActor->GetProperty()->SetPointSize(8);
		// 设置虚线模式，16位模式，这里使用0xF0F0作为示例
		//dynamicActor->GetProperty()->SetLineStipplePattern(0xFFFF);
		//dynamicActor->GetProperty()->SetLineStippleRepeatFactor(1); // 设置重复因子
		renderer->AddActor(dynamicActor);
	}
	else
	{
		mapper->SetInputData(model);
		mapper->Update();
	}

	if (this->isIndependent)
	{
		this->renWinInde->Render();
	}
	else
	{
		this->renWin->Render();
	}
	if (!ifPlotDynamic) this->flameID = 0;
}

void RenderUsingVTK::renderBufferToUg(RenderBuffer& data, std::unordered_map<int, vtkUnstructuredGrid*>& ugModel)
{
	for (auto& pair : ugModel) 
	{
		vtkUnstructuredGrid* grid = pair.second;
		if (grid) 
		{
			grid->Delete();  // 调用VTK对象的Delete方法来释放内存
		}
	}
	ugModel.clear();

	for (auto& part : data.geoModel.parts)
	{
		ugModel[part.first] = vtkUnstructuredGrid::New();

		// 颜色
		if (!ugModel[part.first]->GetCellData()->HasArray("Color"))
		{
			vtkSmartPointer<vtkUnsignedCharArray> cellColors = vtkSmartPointer<vtkUnsignedCharArray>::New();
			cellColors->SetNumberOfComponents(3);  // RGB
			cellColors->SetName("Color");
			ugModel[part.first]->GetCellData()->SetScalars(cellColors);
		}

		// 法向
		if (!ugModel[part.first]->GetCellData()->GetNormals())
		{
			vtkSmartPointer<vtkDoubleArray> cellNormals = vtkSmartPointer<vtkDoubleArray>::New();
			cellNormals->SetNumberOfComponents(3);  // 法向量为 3D 向量
			cellNormals->SetName("Normal");
			ugModel[part.first]->GetCellData()->SetNormals(cellNormals);
		}
	}

	// 所有 ug 共用一套节点
	vtkPoints* points = vtkPoints::New();

	// 颜色配置器
	/*if (colorTransferFunction == nullptr)
	{
		colorTransferFunction = vtkSmartPointer<vtkColorTransferFunction>::New();
	}
	colorTransferFunction->RemoveAllPoints();
	colorTransferFunction->AddRGBPoint(data.minValue, 0, 0, 1);
	double dValue = data.maxValue - data.minValue;
	colorTransferFunction->AddRGBPoint(data.minValue + 0.25 * dValue, 0, 1, 1);
	colorTransferFunction->AddRGBPoint(data.minValue + 0.5 * dValue, 0, 1, 0);
	colorTransferFunction->AddRGBPoint(data.minValue + 0.75 * dValue, 1, 1, 0);
	colorTransferFunction->AddRGBPoint(data.maxValue, 1, 0, 0);*/

	setColorTransferFunction(data);

	double rgbColor[3]; // 用于存储RGB颜色的数组

	ncs::GeometryModel* md = &data.geoModel;
	int ndID = 0;
	// 依次添加 vertex、line、tri、quad、color、normal
	if (data.vertexColor.size() == md->vertex.size())
	{
		for (int i = 0; i < md->vertex.size(); i++)
		{
			vtkSmartPointer<vtkVertex> vertex = vtkSmartPointer<vtkVertex>::New();
			points->InsertNextPoint(md->vertex[i].x, md->vertex[i].y, md->vertex[i].z);

			vtkUnsignedCharArray* cellColors = static_cast<vtkUnsignedCharArray*>
				(ugModel[md->vertex[i].partID]->GetCellData()->GetScalars());
			vtkDoubleArray* cellNormals = static_cast<vtkDoubleArray*>
				(ugModel[md->vertex[i].partID]->GetCellData()->GetNormals());

			// 将颜色绑定到 cellColors
			cellColors->InsertNextTuple3(data.vertexColor[i].r,
				data.vertexColor[i].g,
				data.vertexColor[i].b);

			cellNormals->InsertNextTuple3(0, 1, 0);

			vertex->GetPointIds()->SetId(0, ndID);
			ndID += 1;
			ugModel[md->vertex[i].partID]->InsertNextCell(VTK_VERTEX, vertex->GetPointIds());
		}
	}

	if (data.lineColor.size() == md->lines.size())
	{
		for (int i = 0; i < md->lines.size(); i++)
		{
			// 插入线段的几何信息
			vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
			points->InsertNextPoint(md->nodes[md->lines[i].nd1].x,
				md->nodes[md->lines[i].nd1].y,
				md->nodes[md->lines[i].nd1].z);
			points->InsertNextPoint(md->nodes[md->lines[i].nd2].x,
				md->nodes[md->lines[i].nd2].y,
				md->nodes[md->lines[i].nd2].z);

			vtkUnsignedCharArray* cellColors = static_cast<vtkUnsignedCharArray*>
				(ugModel[md->lines[i].partID]->GetCellData()->GetScalars());
			vtkDoubleArray* cellNormals = static_cast<vtkDoubleArray*>
				(ugModel[md->lines[i].partID]->GetCellData()->GetNormals());

			// 添加单元颜色（基于 cell）
			cellColors->InsertNextTuple3(data.lineColor[i].r,
				data.lineColor[i].g,
				data.lineColor[i].b);

			// 添加法向量（可选）
			cellNormals->InsertNextTuple3(0, 1, 0); // 对点1的法向量
			cellNormals->InsertNextTuple3(0, 1, 0); // 对点2的法向量

			// 关联线段的点
			line->GetPointIds()->SetId(0, ndID);
			line->GetPointIds()->SetId(1, ndID + 1);
			ndID += 2;

			ugModel[md->lines[i].partID]->InsertNextCell(VTK_LINE, line->GetPointIds());
		}
	}

	// std::vector<ncs::Triangle> trisForVtk;
	std::vector<ncs::Quad> quadsForVtk;
	if (md->quads.size() == data.quadColor.size())
	{
		quadsForVtk.assign(md->quads.begin(), md->quads.end());
	}


	int face[6][4] =
	{
		{0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4},
		{1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}
	};

	ncs::Triangle tri;
	ncs::Quad quad;

	if (data.solidColor.size() == md->solids.size() ||
		data.solidValue.size() == md->solids.size())
	{
		for (int i = 0; i < md->solids.size(); i++)
		{
			for (int faceID = 0; faceID < 6; faceID++)
			{
				quad.partID = md->solids[i].partID;
				quad.nds[0] = md->solids[i].nds[face[faceID][0]];
				quad.nds[1] = md->solids[i].nds[face[faceID][1]];
				quad.nds[2] = md->solids[i].nds[face[faceID][2]];
				quad.nds[3] = md->solids[i].nds[face[faceID][3]];
				quadsForVtk.emplace_back(quad);
				if (data.ifUsingColor)
				{
					data.quadColor.emplace_back(data.solidColor[i]);
				}
				else
				{
					colorTransferFunction->GetColor(data.solidValue[i], rgbColor);
					data.quadColor.emplace_back(ncs::Color(rgbColor[0] * 255,
						rgbColor[1] * 255, rgbColor[2] * 255));
				}
			}
		}
	}


	if (md->triangles.size() == data.triColor.size())
	{
		ncs::Node* nd;
		ncs::Node* ndO;
		int triID = 0;
		for (auto& tri : md->triangles) {

			vtkUnsignedCharArray* cellColors = static_cast<vtkUnsignedCharArray*>
				(ugModel[tri.partID]->GetCellData()->GetScalars());
			vtkDoubleArray* cellNormals = static_cast<vtkDoubleArray*>
				(ugModel[tri.partID]->GetCellData()->GetNormals());

			//插入数据到vtk数据结构中
			vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();

			ncs::Vec3 v1, v2, nml;
			ncs::Node* n1, * n2, * n3;

			n1 = &(md->nodes)[tri.nds[0]];
			n2 = &(md->nodes)[tri.nds[1]];
			n3 = &(md->nodes)[tri.nds[2]];
			v1.SetValue(n2->x - n1->x, n2->y - n1->y, n2->z - n1->z);
			v2.SetValue(n3->x - n1->x, n3->y - n1->y, n3->z - n1->z);
			nml = v1.cross(v2);
			nml.normalize();

			for (int i = 0; i < 3; i++)
			{
				points->InsertNextPoint(md->nodes[tri.nds[i]].x, md->nodes[tri.nds[i]].y,
					md->nodes[tri.nds[i]].z);
				triangle->GetPointIds()->SetId(i, ndID + i);
			}

			cellColors->InsertNextTuple3(data.triColor[triID].r, data.triColor[triID].g,
				data.triColor[triID].b);
			cellNormals->InsertNextTuple3(nml.x, nml.y, nml.z);

			ndID += 3;
			triID += 1;
			//插入数据
			ugModel[tri.partID]->InsertNextCell(VTK_TRIANGLE, triangle->GetPointIds());
		}
	}

	if (quadsForVtk.size() == data.quadColor.size())
	{
		int quadID = 0;
		for (const auto& qd : quadsForVtk) 
		{
			vtkUnsignedCharArray* cellColors = static_cast<vtkUnsignedCharArray*>
				(ugModel[qd.partID]->GetCellData()->GetScalars());
			vtkDoubleArray* cellNormals = static_cast<vtkDoubleArray*>
				(ugModel[qd.partID]->GetCellData()->GetNormals());

			vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();

			// 计算法向量
			ncs::Vec3 v1, v2, nml;
			ncs::Node* n1, * n2, * n3;
			n1 = &(md->nodes)[qd.nds[0]];
			n2 = &(md->nodes)[qd.nds[1]];
			n3 = &(md->nodes)[qd.nds[2]];
			v1.SetValue(n2->x - n1->x, n2->y - n1->y, n2->z - n1->z);
			v2.SetValue(n3->x - n1->x, n3->y - n1->y, n3->z - n1->z);
			nml = v1.cross(v2);
			nml.normalize();

			for (int i = 0; i < 4; i++)
			{
				points->InsertNextPoint(md->nodes[qd.nds[i]].x, md->nodes[qd.nds[i]].y,
					md->nodes[qd.nds[i]].z);
				quad->GetPointIds()->SetId(i, ndID + i);
			}
			ndID += 4;

			// 设置单元颜色和法向量
			cellColors->InsertNextTuple3(data.quadColor[quadID].r, data.quadColor[quadID].g,
				data.quadColor[quadID].b);
			cellNormals->InsertNextTuple3(nml.x, nml.y, nml.z);

			ugModel[qd.partID]->InsertNextCell(VTK_QUAD, quad->GetPointIds());
			quadID++;
		}
	}

	for (auto& ug : ugModel)
	{
		ug.second->SetPoints(points);
	}
}

vtkUnstructuredGrid* RenderUsingVTK::RenderBufferToUg(RenderBuffer& data)
{
#if 1
	ncs::GeometryModel* md = &data.geoModel;

	vtkPoints* points;
	vtkCellArray* cells;
	vtkUnsignedCharArray* cellTypes;
	vtkUnstructuredGrid* ug = vtkUnstructuredGrid::New();
	points = vtkPoints::New();
	cells = vtkCellArray::New();
	cellTypes = vtkUnsignedCharArray::New();
	cellTypes->SetNumberOfComponents(1);

	// 插入points的颜色数据到vtk数据结构中
	//if (!ug->GetPointData()->HasArray("Color"))
	//{
	//	vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
	//	colors->SetNumberOfComponents(3);
	//	colors->SetName("Color");
	//	ug->GetPointData()->AddArray(colors);
	//}
	if (!ug->GetCellData()->HasArray("Color"))
	{
		vtkSmartPointer<vtkUnsignedCharArray> cellColors = vtkSmartPointer<vtkUnsignedCharArray>::New();
		cellColors->SetNumberOfComponents(3);  // RGB
		cellColors->SetName("Color");
		ug->GetCellData()->SetScalars(cellColors);
	}
	vtkUnsignedCharArray* cellColors = static_cast<vtkUnsignedCharArray*>(ug->GetCellData()->GetScalars());

	////插入points的颜色数据到vtk数据结构中
	//if (ug->GetPointData()->GetNormals() == nullptr)
	//{
	//	vtkSmartPointer<vtkDoubleArray> normals = vtkSmartPointer<vtkDoubleArray>::New();
	//	normals->SetNumberOfComponents(3);
	//	normals->SetName("Normal");
	//	ug->GetPointData()->SetNormals(normals);
	//}
	//vtkDoubleArray* normals = static_cast<vtkDoubleArray*>(ug->GetPointData()->GetNormals());

	if (!ug->GetCellData()->GetNormals())
	{
		vtkSmartPointer<vtkDoubleArray> cellNormals = vtkSmartPointer<vtkDoubleArray>::New();
		cellNormals->SetNumberOfComponents(3);  // 法向量为 3D 向量
		cellNormals->SetName("Normal");
		ug->GetCellData()->SetNormals(cellNormals);
	}
	vtkDoubleArray* cellNormals = static_cast<vtkDoubleArray*>(ug->GetCellData()->GetNormals());

	// 颜色配置器
	if (colorTransferFunction == nullptr)
	{
		colorTransferFunction = vtkSmartPointer<vtkColorTransferFunction>::New();
	}
	colorTransferFunction->RemoveAllPoints();
	colorTransferFunction->AddRGBPoint(data.minValue, 0, 0, 1);
	double dValue = data.maxValue - data.minValue;
	colorTransferFunction->AddRGBPoint(data.minValue + 0.25 * dValue, 0, 1, 1);
	colorTransferFunction->AddRGBPoint(data.minValue + 0.5 * dValue, 0, 1, 0);
	colorTransferFunction->AddRGBPoint(data.minValue + 0.75 * dValue, 1, 1, 0);
	colorTransferFunction->AddRGBPoint(data.maxValue, 1, 0, 0);
	double rgbColor[3]; // 用于存储RGB颜色的数组

	int ndID = points->GetNumberOfPoints();
	// 依次添加 vertex、line、tri、quad、color、normal
	if (data.vertexColor.size() == md->vertex.size())
	{
		for (int i = 0; i < md->vertex.size(); i++)
		{
			vtkSmartPointer<vtkVertex> vertex = vtkSmartPointer<vtkVertex>::New();
			points->InsertNextPoint(md->vertex[i].x, md->vertex[i].y, md->vertex[i].z);

			// 将颜色绑定到 cellColors 而不是 pointColors
			cellColors->InsertNextTuple3(data.vertexColor[i].r,
				data.vertexColor[i].g,
				data.vertexColor[i].b);

			cellNormals->InsertNextTuple3(0, 0, 1);
			vertex.Get()->GetPointIds()->SetId(0, ndID);
			ndID += 1;
			cells->InsertNextCell(vertex);
			cellTypes->InsertNextTuple1(VTK_VERTEX);
		}
	}


	if (data.lineColor.size() == md->lines.size())
	{
		for (int i = 0; i < md->lines.size(); i++)
		{
			// 插入线段的几何信息
			vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
			points->InsertNextPoint(md->nodes[md->lines[i].nd1].x,
				md->nodes[md->lines[i].nd1].y,
				md->nodes[md->lines[i].nd1].z);
			points->InsertNextPoint(md->nodes[md->lines[i].nd2].x,
				md->nodes[md->lines[i].nd2].y,
				md->nodes[md->lines[i].nd2].z);

			// 添加单元颜色（基于 cell）
			cellColors->InsertNextTuple3(data.lineColor[i].r,
				data.lineColor[i].g,
				data.lineColor[i].b);

			// 添加法向量（可选）
			cellNormals->InsertNextTuple3(0, 1, 0); // 对点1的法向量
			cellNormals->InsertNextTuple3(0, 1, 0); // 对点2的法向量

			// 关联线段的点
			line->GetPointIds()->SetId(0, ndID);
			line->GetPointIds()->SetId(1, ndID + 1);
			ndID += 2;

			// 将线段加入到单元集合
			cells->InsertNextCell(line);
			cellTypes->InsertNextTuple1(VTK_LINE);
		}
	}



	// std::vector<ncs::Triangle> trisForVtk;
	std::vector<ncs::Quad> quadsForVtk;
	if (md->quads.size() == data.quadColor.size())
	{
		quadsForVtk.assign(md->quads.begin(), md->quads.end());
	}


	int face[6][4] =
	{
		{0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4},
		{1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}
	};

	ncs::Triangle tri;
	ncs::Quad quad;

	if (data.solidColor.size() == md->solids.size() ||
		data.solidValue.size() == md->solids.size())
	{
		for (int i = 0; i < md->solids.size(); i++)
		{
			for (int faceID = 0; faceID < 6; faceID++)
			{
				quad.nds[0] = md->solids[i].nds[face[faceID][0]];
				quad.nds[1] = md->solids[i].nds[face[faceID][1]];
				quad.nds[2] = md->solids[i].nds[face[faceID][2]];
				quad.nds[3] = md->solids[i].nds[face[faceID][3]];
				quadsForVtk.emplace_back(quad);
				if (data.ifUsingColor)
				{
					data.quadColor.emplace_back(data.solidColor[i]);
				}
				else
				{
					colorTransferFunction->GetColor(data.solidValue[i], rgbColor);
					data.quadColor.emplace_back(ncs::Color(rgbColor[0] * 255,
						rgbColor[1] * 255, rgbColor[2] * 255));
				}
			}
		}
	}


	if (md->triangles.size() == data.triColor.size())
	{
		ncs::Node* nd;
		ncs::Node* ndO;
		int triID = 0;
		for (auto& tri : md->triangles) {
			//插入数据到vtk数据结构中
			vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();

			ncs::Vec3 v1, v2, nml;
			ncs::Node* n1, * n2, * n3;

			n1 = &(md->nodes)[tri.nds[0]];
			n2 = &(md->nodes)[tri.nds[1]];
			n3 = &(md->nodes)[tri.nds[2]];
			v1.SetValue(n2->x - n1->x, n2->y - n1->y, n2->z - n1->z);
			v2.SetValue(n3->x - n1->x, n3->y - n1->y, n3->z - n1->z);
			nml = v1.cross(v2);
			nml.normalize();

			for (int i = 0; i < 3; i++)
			{
				points->InsertNextPoint(md->nodes[tri.nds[i]].x, md->nodes[tri.nds[i]].y,
					md->nodes[tri.nds[i]].z);

				cellColors->InsertNextTuple3(data.triColor[triID].r, data.triColor[triID].g,
					data.triColor[triID].b);
				cellNormals->InsertNextTuple3(nml.x, nml.y, nml.z);
				triangle->GetPointIds()->SetId(i, ndID + i);
			}
			ndID += 3;
			triID += 1;
			//插入数据
			cells->InsertNextCell(triangle);
			cellTypes->InsertNextTuple1(VTK_TRIANGLE);
		}
	}

	if (quadsForVtk.size() == data.quadColor.size())
	{
		int quadID = 0;
		for (const auto& qd : quadsForVtk) {
			vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();

			// 计算法向量
			ncs::Vec3 v1, v2, nml;
			ncs::Node* n1, * n2, * n3;
			n1 = &(md->nodes)[qd.nds[0]];
			n2 = &(md->nodes)[qd.nds[1]];
			n3 = &(md->nodes)[qd.nds[2]];
			v1.SetValue(n2->x - n1->x, n2->y - n1->y, n2->z - n1->z);
			v2.SetValue(n3->x - n1->x, n3->y - n1->y, n3->z - n1->z);
			nml = v1.cross(v2);
			nml.normalize();

			for (int i = 0; i < 4; i++)
			{
				points->InsertNextPoint(md->nodes[qd.nds[i]].x, md->nodes[qd.nds[i]].y,
					md->nodes[qd.nds[i]].z);
				quad->GetPointIds()->SetId(i, ndID + i);
			}
			ndID += 4;

			// 设置单元颜色和法向量
			cellColors->InsertNextTuple3(data.quadColor[quadID].r, data.quadColor[quadID].g,
				data.quadColor[quadID].b);
			cellNormals->InsertNextTuple3(nml.x, nml.y, nml.z);

			// 插入单元
			cells->InsertNextCell(quad);
			cellTypes->InsertNextTuple1(VTK_QUAD);
			quadID++;
		}

	}

	//组装模型
	ug->SetPoints(points);
	ug->SetCells(cellTypes, cells);
	ug->GetCellData()->SetScalars(cellColors);
	return ug;
#else
	ncs::GeometryModel* md = &data.geoModel;
	//if (md->nodes.size() == 0) return nullptr;

	vtkPoints* points;
	vtkCellArray* cells;
	vtkUnsignedCharArray* cellTypes;
	vtkUnstructuredGrid* ug = vtkUnstructuredGrid::New();
	points = vtkPoints::New();
	cells = vtkCellArray::New();
	cellTypes = vtkUnsignedCharArray::New();
	cellTypes->SetNumberOfComponents(1);

	// 插入points的颜色数据到vtk数据结构中
	if (!ug->GetPointData()->HasArray("Color"))
	{
		vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
		colors->SetNumberOfComponents(3);
		colors->SetName("Color");
		ug->GetPointData()->AddArray(colors);
	}
	vtkUnsignedCharArray* colors = static_cast<vtkUnsignedCharArray*>(ug->GetPointData()->GetArray("Color"));

	//插入points的颜色数据到vtk数据结构中
	if (ug->GetPointData()->GetNormals() == nullptr)
	{
		vtkSmartPointer<vtkDoubleArray> normals = vtkSmartPointer<vtkDoubleArray>::New();
		normals->SetNumberOfComponents(3);
		normals->SetName("Normal");
		ug->GetPointData()->SetNormals(normals);
	}
	vtkDoubleArray* normals = static_cast<vtkDoubleArray*>(ug->GetPointData()->GetNormals());

	// 颜色配置器
	if (colorTransferFunction == nullptr)
	{
		colorTransferFunction =  vtkSmartPointer<vtkColorTransferFunction>::New();
	}
	colorTransferFunction->RemoveAllPoints();
	colorTransferFunction->AddRGBPoint(data.minValue, 0, 0, 1);
	double dValue = data.maxValue - data.minValue;
	colorTransferFunction->AddRGBPoint(data.minValue + 0.25 * dValue, 0, 1, 1);
	colorTransferFunction->AddRGBPoint(data.minValue + 0.5 * dValue, 0, 1, 0);
	colorTransferFunction->AddRGBPoint(data.minValue + 0.75 * dValue, 1, 1, 0);
	colorTransferFunction->AddRGBPoint(data.maxValue, 1, 0, 0);
	double rgbColor[3]; // 用于存储RGB颜色的数组

	int ndID = points->GetNumberOfPoints();
	// 依次添加 vertex、line、tri、quad、color、normal
	if (data.vertexColor.size() == md->vertex.size())
	{
		for (int i = 0; i < md->vertex.size(); i++)
		{
			vtkSmartPointer<vtkVertex> vertex = vtkSmartPointer<vtkVertex>::New();
			points->InsertNextPoint(md->vertex[i].x, md->vertex[i].y, md->vertex[i].z);
			colors->InsertNextTuple3(data.vertexColor[i].r, data.vertexColor[i].g,
				data.vertexColor[i].b);
			normals->InsertNextTuple3(0, 0, 1);
			vertex.Get()->GetPointIds()->SetId(0, ndID);
			ndID += 1;
			cells->InsertNextCell(vertex);
			cellTypes->InsertNextTuple1(VTK_VERTEX);
		}
	}
	
	if (data.lineColor.size() == md->lines.size())
	{
		for (int i = 0; i < md->lines.size(); i++) 
		{
			//插入数据到vtk数据结构中
			vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
			points->InsertNextPoint(md->nodes[md->lines[i].nd1].x, md->nodes[md->lines[i].nd1].y,
				md->nodes[md->lines[i].nd1].z);
			points->InsertNextPoint(md->nodes[md->lines[i].nd2].x, md->nodes[md->lines[i].nd2].y,
				md->nodes[md->lines[i].nd2].z);
			colors->InsertNextTuple3(data.lineColor[i].r, data.lineColor[i].g,
				data.lineColor[i].b);
			colors->InsertNextTuple3(data.lineColor[i].r, data.lineColor[i].g,
				data.lineColor[i].b);
			normals->InsertNextTuple3(0, 1, 0);
			normals->InsertNextTuple3(0, 1, 0);

			line.Get()->GetPointIds()->SetId(0, ndID);
			line.Get()->GetPointIds()->SetId(1, ndID + 1);
			ndID += 2;
			//插入数据
			cells->InsertNextCell(line);
			cellTypes->InsertNextTuple1(VTK_LINE);
		}
	}
	

	// std::vector<ncs::Triangle> trisForVtk;
	std::vector<ncs::Quad> quadsForVtk;
	if (md->quads.size() == data.quadColor.size())
	{
		quadsForVtk.assign(md->quads.begin(), md->quads.end());
	}
	

	int face[6][4] =
	{
		{0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4},
		{1, 2, 6, 5}, {2, 3, 7, 6}, {3, 0, 4, 7}
	};

	ncs::Triangle tri;
	ncs::Quad quad;
	
	if (data.solidColor.size() == md->solids.size() || 
		data.solidValue.size() == md->solids.size())
	{
		for (int i = 0; i < md->solids.size(); i++)
		{
			for (int faceID = 0; faceID < 6; faceID++)
			{
				quad.nds[0] = md->solids[i].nds[face[faceID][0]];
				quad.nds[1] = md->solids[i].nds[face[faceID][1]];
				quad.nds[2] = md->solids[i].nds[face[faceID][2]];
				quad.nds[3] = md->solids[i].nds[face[faceID][3]];
				quadsForVtk.emplace_back(quad);
				if (data.ifUsingColor)
				{
					data.quadColor.emplace_back(data.solidColor[i]);
				}
				else
				{
					colorTransferFunction->GetColor(data.solidValue[i], rgbColor);
					data.quadColor.emplace_back(ncs::Color(rgbColor[0] * 255,
						rgbColor[1] * 255, rgbColor[2] * 255));
				}
			}
		}
	}
	

	if (md->triangles.size() == data.triColor.size())
	{
		ncs::Node* nd;
		ncs::Node* ndO;
		int triID = 0;
		for (auto& tri : md->triangles) {
			//插入数据到vtk数据结构中
			vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();

			ncs::Vec3 v1, v2, nml;
			ncs::Node* n1, * n2, * n3;

			n1 = &(md->nodes)[tri.nds[0]];
			n2 = &(md->nodes)[tri.nds[1]];
			n3 = &(md->nodes)[tri.nds[2]];
			v1.SetValue(n2->x - n1->x, n2->y - n1->y, n2->z - n1->z);
			v2.SetValue(n3->x - n1->x, n3->y - n1->y, n3->z - n1->z);
			nml = v1.cross(v2);
			nml.normalize();

			for (int i = 0; i < 3; i++)
			{
				points->InsertNextPoint(md->nodes[tri.nds[i]].x, md->nodes[tri.nds[i]].y,
					md->nodes[tri.nds[i]].z);

				colors->InsertNextTuple3(data.triColor[triID].r, data.triColor[triID].g,
					data.triColor[triID].b);
				normals->InsertNextTuple3(nml.x, nml.y, nml.z);
				triangle->GetPointIds()->SetId(i, ndID + i);
			}
			ndID += 3;
			triID += 1;
			//插入数据
			cells->InsertNextCell(triangle);
			cellTypes->InsertNextTuple1(VTK_TRIANGLE);
		}
	}
	
	if (quadsForVtk.size() == data.quadColor.size())
	{
		int quadID = 0;
		for (const auto& qd : quadsForVtk) {
			//插入数据到vtk数据结构中
			vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();

			ncs::Vec3 v1, v2, nml;
			ncs::Node* n1, * n2, * n3;

			n1 = &(md->nodes)[qd.nds[0]];
			n2 = &(md->nodes)[qd.nds[1]];
			n3 = &(md->nodes)[qd.nds[2]];
			v1.SetValue(n2->x - n1->x, n2->y - n1->y, n2->z - n1->z);
			v2.SetValue(n3->x - n1->x, n3->y - n1->y, n3->z - n1->z);
			nml = v1.cross(v2);
			nml.normalize();

			for (int i = 0; i < 4; i++)
			{
				points->InsertNextPoint(md->nodes[qd.nds[i]].x, md->nodes[qd.nds[i]].y,
					md->nodes[qd.nds[i]].z);
				colors->InsertNextTuple3(data.quadColor[quadID].r, data.quadColor[quadID].g,
					data.quadColor[quadID].b);
				normals->InsertNextTuple3(nml.x, nml.y, nml.z);
				quad->GetPointIds()->SetId(i, ndID + i);
			}
			ndID += 4;
			quadID++;
			//插入数据
			cells->InsertNextCell(quad);
			cellTypes->InsertNextTuple1(VTK_QUAD);
		}
	}
	
	//组装模型
	ug->SetPoints(points);
	ug->SetCells(cellTypes, cells);
	ug->GetPointData()->SetScalars(ug->GetPointData()->GetArray("Color"));
	return ug;
#endif
}

void RenderUsingVTK::createColorBar(vtkUnstructuredGrid* data, RenderBuffer& renderInfo)
{
	if (this->scalarBar == nullptr)
	{
		// contourFilter = vtkSmartPointer<vtkContourFilter>::New();
		contourMapper = vtkSmartPointer<vtkDataSetMapper>::New();
		contourActor = vtkSmartPointer<vtkActor>::New();
		scalarBar = vtkSmartPointer<vtkScalarBarActor>::New();
		renderer->AddActor2D(scalarBar);
		renderer->AddActor(contourActor);
		titleActor = vtkSmartPointer<vtkTextActor>::New();
		renderer->AddActor2D(titleActor);
	}
	// 颜色配置器
	setColorTransferFunction(renderInfo);

	contourMapper->SetLookupTable(colorTransferFunction);
	// 创建演员
	contourActor->SetMapper(contourMapper);
	// 添加演员到渲染器
	// renderer->AddActor(contourActor);
	// 创建色棒
	scalarBar->SetLookupTable(colorTransferFunction);

	scalarBar->SetNumberOfLabels(renderInfo.colorNum);

# if 1
	// 水平色带
	scalarBar.Get()->SetOrientationToHorizontal();
	scalarBar.Get()->SetPosition(0.25, 0.1);
	scalarBar->SetWidth(0.5);
	scalarBar->SetHeight(0.08);

	// 获取 scalarBar 的位置和大小
	double* scalarBarPos = scalarBar->GetPosition();       // 左下角位置
	double* scalarBarSize = scalarBar->GetPosition2();     // 宽度和高度

	// 计算自定义标签的位置
	// 例如：将自定义标签放置在颜色条的正下方
	double labelX = scalarBarPos[0] + scalarBarSize[0] / 2.0;   // 水平居中
	double labelY = scalarBarPos[1] - scalarBarSize[1] * 0.8;       // 在颜色条正下方，稍微向下偏移

	// 创建独立的文本标签
	// if(titleActor = nullptr) titleActor = vtkSmartPointer<vtkTextActor>::New();
	
	titleActor->SetInput((renderInfo.colorBarName + '\n').c_str());  // 设置标题内容
	titleActor->GetTextProperty()->SetFontSize(15);
	titleActor->GetTextProperty()->SetColor(0.0, 0.0, 0.0);  // 设置字体颜色
	titleActor->GetTextProperty()->SetFontFamilyToTimes();
	titleActor->GetTextProperty()->SetItalic(1);  // 设置字体为斜体
	titleActor->GetTextProperty()->SetBold(1);
	// 设置标题的位置到色带下方
	renderer->NormalizedDisplayToDisplay(labelX, labelY);
	titleActor->SetPosition(labelX, labelY);  // 根据具体渲染窗口的坐标调整
	//std::cout << "init: labrlX: " << labelX << "    labrlY: " << labelY << '\n';

	if (this->isIndependent)
	{
		this->renWinInde->GetInteractor()->AddObserver(vtkCommand::WindowResizeEvent, this,
			&RenderUsingVTK::UpdateTitlePosition);
	}
	else
	{
		this->renWin->GetInteractor()->AddObserver(vtkCommand::WindowResizeEvent, this,
			&RenderUsingVTK::UpdateTitlePosition);
	}

#else
	// 竖直色带
	scalarBar->SetTitle((renderInfo.colorBarName + '\n').c_str());
	scalarBar->SetWidth(0.06);
	scalarBar->SetHeight(1.0);
	vtkTextProperty* textProperty = scalarBar->GetTitleTextProperty();
	// textProperty->SetLineOffset(12);
	textProperty->SetFontSize(15); // 设置字体大小
	textProperty->SetColor(0.0, 0.0, 0.0); // 设置字体颜色
	textProperty->SetFontFamilyToTimes();
	textProperty->SetVerticalJustificationToBottom();
	textProperty->SetLineOffset(200);
#endif
	vtkTextProperty* contentProperty = scalarBar->GetLabelTextProperty();
	contentProperty->SetFontSize(12);  // 设置字体大小
	contentProperty->SetColor(0.0, 0.0, 0.0); // 设置字体颜色
	contentProperty->SetFontFamilyToTimes();
	contentProperty->SetJustificationToRight();  // 设置标签在色带下方

	scalarBar->SetLabelFormat("%.3f"); // 保留有效位数
	scalarBar->Modified();
	scalarBar->SetUnconstrainedFontSize(true);
	
}

void RenderUsingVTK::setColorTransferFunction(RenderBuffer& renderInfo)
{
	// 颜色配置器
	if (colorTransferFunction == nullptr)
	{
		colorTransferFunction = vtkSmartPointer<vtkColorTransferFunction>::New();

	}
	colorTransferFunction->RemoveAllPoints();
	double dValue = renderInfo.maxValue - renderInfo.minValue;

	colorTransferFunction->RemoveAllPoints();
	colorTransferFunction->AddRGBPoint(renderInfo.minValue, 0.0, 0.0, 1.0);  // 蓝色
	colorTransferFunction->AddRGBPoint(renderInfo.minValue + 0.2 * dValue, 0.0, 1.0, 1.0);  // 青色
	colorTransferFunction->AddRGBPoint(renderInfo.minValue + 0.4 * dValue, 0.0, 1.0, 0.0);  // 绿色
	colorTransferFunction->AddRGBPoint(renderInfo.minValue + 0.6 * dValue, 1.0, 1.0, 0.0);  // 黄色
	colorTransferFunction->AddRGBPoint(renderInfo.minValue + 0.8 * dValue, 1.0, 0.5, 0.0);  // 橙色
	colorTransferFunction->AddRGBPoint(renderInfo.maxValue, 1.0, 0.0, 0.0);  // 红色
}

void RenderUsingVTK::UpdateTitlePosition(vtkObject* caller, unsigned long eventId, void* callData)
{
	// 获取 scalarBar 的位置和大小
	double* scalarBarPos = scalarBar->GetPosition();       // 左下角位置
	double* scalarBarSize = scalarBar->GetPosition2();     // 宽度和高度

	// 计算自定义标签的位置
	// 例如：将自定义标签放置在颜色条的正下方
	double labelX = scalarBarPos[0] + scalarBarSize[0] / 2.0;   // 水平居中
	double labelY = scalarBarPos[1] - scalarBarSize[1] * 0.8;       // 在颜色条正下方，稍微向下偏移

	// 设置标题的位置到色带下方
	renderer->NormalizedDisplayToDisplay(labelX, labelY);
	titleActor->SetPosition(labelX, labelY);  // 根据具体渲染窗口的坐标调整
	//std::cout << "adjust: labrlX: " << labelX << "    labrlY: " << labelY << '\n';
}


void RenderUsingVTK::showAxies()
{
#if 0
	vtkAxesActor* axes_actor = vtkAxesActor::New();
	axes_actor->SetAxisLabels(0);
	axes_actor->SetTotalLength(0.2, 0.2, 0.2);
	//// 获取坐标轴名称的文本属性
	//vtkTextProperty* labelXProperty = axes_actor->GetXAxisCaptionActor2D()->GetCaptionTextProperty();
	//vtkTextProperty* labelYProperty = axes_actor->GetYAxisCaptionActor2D()->GetCaptionTextProperty();
	//vtkTextProperty* labelZProperty = axes_actor->GetZAxisCaptionActor2D()->GetCaptionTextProperty();

	//// labelXProperty->SetFontFamilyToArial(); // 设置字体家族（可选）
	//labelXProperty->SetFontSize(1);    // 设置字体大小
	//labelYProperty->SetFontSize(1);    // 设置字体大小
	//labelZProperty->SetFontSize(1);    // 设置字体大小

	//// 设置字体颜色
	//double textColor[3] = { 0.0, 0.0, 0.0 }; // 红色
	//labelXProperty->SetColor(textColor);
	//labelYProperty->SetColor(textColor);
	//labelZProperty->SetColor(textColor);

	//// 设置字体为 Times New Roman
	//labelXProperty->SetFontFamily(VTK_TIMES);
	//labelYProperty->SetFontFamily(VTK_TIMES);
	//labelZProperty->SetFontFamily(VTK_TIMES);

	//// 设置为斜体（示例）
	//labelXProperty->ItalicOn();
	//labelYProperty->ItalicOn();
	//labelZProperty->ItalicOn();

	//vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
	//transform->Translate(10, 10, 10); // 设置新的原点位置

	double position[3]{ 1, 1, 1 };
	axes_actor->SetPosition(position);

	this->renderer->AddActor(axes_actor);
	//this->renderer->Render();

#else
	axes_actor = vtkAxesActor::New();
	axes_actor->SetAxisLabels(1);
	// axes_actor->SetTotalLength(1, 1, 1);
	vtkOrientationMarkerWidget* markerWidget = vtkOrientationMarkerWidget::New();
	markerWidget->SetOrientationMarker(axes_actor);
	markerWidget->SetOutlineColor(0.93, 0.57, 0.13);

	axes_actor->GetXAxisShaftProperty()->SetLineWidth(3.0);
	axes_actor->GetYAxisShaftProperty()->SetLineWidth(3.0);
	axes_actor->GetZAxisShaftProperty()->SetLineWidth(3.0);

	if (this->isIndependent)
	{
		markerWidget->SetInteractor(this->renWinInde->GetInteractor());
	}
	else
	{
		markerWidget->SetInteractor(this->renWin->GetInteractor());
	}
	
	markerWidget->SetViewport(0, 0, 0.25, 0.25);
	// markerWidget->SetViewport(0.25, 0.25, 0.5, 0.5);
	markerWidget->SetZoom(1);
	markerWidget->SetEnabled(true);
	markerWidget->InteractiveOff();
#endif
}

void RenderUsingVTK::drawSpring(double point1[3], double point2[3], double color[3])
{
	//// 给定两个端点
	//double point1[3] = { 0.0, 0.0, 0.0 }; // 起点
	//double point2[3] = { 1.0, 1.0, 1.0 }; // 终点

	// 弹簧参数
	double springRadius = 0.1;    // 弹簧半径
	int numberOfTurns = 12;        // 螺旋圈数
	int numPointsPerTurn = 100;   // 每圈的点数
	double lineThickness = 0.01;
	double springHeight = vtkMath::Distance2BetweenPoints(point1, point2); // 弹簧的总高度
	double pitch = springHeight / numberOfTurns; // 每圈的螺距

	// 计算方向向量
	double direction[3];
	vtkMath::Subtract(point2, point1, direction);
	vtkMath::Normalize(direction); // 归一化方向向量

	// 构建一个正交基，保证局部坐标的螺旋平面与方向向量正交
	double arbitraryVector[3];
	if (fabs(direction[0]) > fabs(direction[1]) && fabs(direction[0]) > fabs(direction[2])) {
		arbitraryVector[0] = 0.0;
		arbitraryVector[1] = 1.0;
		arbitraryVector[2] = 0.0;
	}
	else {
		arbitraryVector[0] = 1.0;
		arbitraryVector[1] = 0.0;
		arbitraryVector[2] = 0.0;
	}

	// 计算正交基
	double tangent[3];
	vtkMath::Cross(arbitraryVector, direction, tangent);
	vtkMath::Normalize(tangent);

	double binormal[3];
	vtkMath::Cross(direction, tangent, binormal);
	vtkMath::Normalize(binormal);

	// 创建螺旋线的点集
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	for (int i = 0; i <= numberOfTurns * numPointsPerTurn; ++i) {
		double theta = 2 * vtkMath::Pi() * i / numPointsPerTurn; // 当前点的角度
		double localX = springRadius * cos(theta);              // 螺旋截面局部 x
		double localY = springRadius * sin(theta);              // 螺旋截面局部 y
		double localZ = pitch * i / numPointsPerTurn;           // 局部 z（沿方向向量的偏移）

		// 将局部坐标转换到全局坐标
		double globalPoint[3] = {
			point1[0] + tangent[0] * localX + binormal[0] * localY + direction[0] * localZ,
			point1[1] + tangent[1] * localX + binormal[1] * localY + direction[1] * localZ,
			point1[2] + tangent[2] * localX + binormal[2] * localY + direction[2] * localZ
		};
		points->InsertNextPoint(globalPoint);
	}

	// 创建螺旋线的多段线数据
	vtkSmartPointer<vtkPolyLine> polyLine = vtkSmartPointer<vtkPolyLine>::New();
	int totalPoints = numberOfTurns * numPointsPerTurn + 1;
	polyLine->GetPointIds()->SetNumberOfIds(totalPoints);
	for (int i = 0; i < totalPoints; ++i) {
		polyLine->GetPointIds()->SetId(i, i);
	}

	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	lines->InsertNextCell(polyLine);

	// 创建螺旋线数据
	vtkSmartPointer<vtkPolyData> helixData = vtkSmartPointer<vtkPolyData>::New();
	helixData->SetPoints(points);
	helixData->SetLines(lines);

	// 将螺旋线转换为管状结构
	vtkSmartPointer<vtkTubeFilter> tubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
	tubeFilter->SetInputData(helixData);
	tubeFilter->SetRadius(lineThickness); // 管状结构的半径
	tubeFilter->SetNumberOfSides(10); // 管状结构的侧面数

	// 映射和渲染
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(tubeFilter->GetOutputPort());

	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetColor(color[0], color[1], color[2]); // 设置颜色

	renderer->AddActor(actor);
}

void RenderUsingVTK::DrawCurvedArrow(
	double radius,         // 弯曲的半径
	double startAngle,     // 起始角度（弧度）
	double endAngle,       // 结束角度（弧度）
	int numPoints,         // 曲线点数量
	double tubeRadius,     // 管径
	double arrowTipLength, // 箭头长度
	double arrowColor[3],
	int axis,
	double pos[3])  // 箭头颜色 (RGB)
{
	// 生成弯曲箭头的点
	vtkSmartPointer<vtkPoints> arrowPoints = vtkSmartPointer<vtkPoints>::New();
	double angleStep = (endAngle - startAngle) / (numPoints - 1);
	for (int i = 0; i < numPoints; ++i) {
		double angle = startAngle + i * angleStep;

		double x = 0;
		double y = 0;
		double z = 0;
		if (axis == 0)
		{
			x = pos[0];
			y = -radius * cos(angle) + pos[1];
			z = -radius * sin(angle) + pos[2];
		}
		else if (axis == 1)
		{
			x = -radius * sin(angle) + pos[0];
			y = pos[1];
			z = -radius * cos(angle) + pos[2];
		}
		else
		{
			x = -radius * sin(angle) + pos[0];
			y = radius * cos(angle) + pos[1];
			z = pos[2];
		}
		
		arrowPoints->InsertNextPoint(x, y, z);
	}

	// 创建箭头的曲线
	vtkSmartPointer<vtkPolyLine> polyLine = vtkSmartPointer<vtkPolyLine>::New();
	polyLine->GetPointIds()->SetNumberOfIds(numPoints);
	for (int i = 0; i < numPoints; ++i) {
		polyLine->GetPointIds()->SetId(i, i);
	}

	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	lines->InsertNextCell(polyLine);

	vtkSmartPointer<vtkPolyData> curveData = vtkSmartPointer<vtkPolyData>::New();
	curveData->SetPoints(arrowPoints);
	curveData->SetLines(lines);

	// 转换为管状结构
	vtkSmartPointer<vtkTubeFilter> tubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
	tubeFilter->SetInputData(curveData);
	tubeFilter->SetRadius(tubeRadius);
	tubeFilter->SetNumberOfSides(20);

	// 映射并创建曲线箭头部分的Actor
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(tubeFilter->GetOutputPort());

	vtkSmartPointer<vtkActor> arrowBodyActor = vtkSmartPointer<vtkActor>::New();
	arrowBodyActor->SetMapper(mapper);
	arrowBodyActor->GetProperty()->SetColor(arrowColor[0], arrowColor[1], arrowColor[2]);

	// 添加箭头的尾部部分到渲染器
	renderer->AddActor(arrowBodyActor);

	// 创建箭头的头部（尖端）
	double arrowTipPosition[3];
	double arrowDirection[3];
	if (axis == 0)
	{
		arrowTipPosition[0] = pos[0];
		arrowTipPosition[1] = -radius * cos(endAngle) + pos[1];
		arrowTipPosition[2] = -radius * sin(endAngle) + pos[2];

		arrowDirection[0] = 0.0;
		arrowDirection[1] = -sin(endAngle);
		arrowDirection[2] = cos(endAngle);
	}
	else if (axis == 1)
	{
		arrowTipPosition[0] = -radius * sin(endAngle) + pos[0];
		arrowTipPosition[1] = pos[1];
		arrowTipPosition[2] = -radius * cos(endAngle) + pos[2];

		arrowDirection[0] = cos(endAngle);
		arrowDirection[1] = 0;
		arrowDirection[2] = -sin(endAngle);
	}
	else
	{
		arrowTipPosition[0] = -radius * sin(endAngle) + pos[0];
		arrowTipPosition[1] = radius * cos(endAngle) + pos[1];
		arrowTipPosition[2] = pos[2];

		arrowDirection[0] = cos(endAngle);
		arrowDirection[1] = sin(endAngle);
		arrowDirection[2] = 0.0;
	}
	vtkMath::Normalize(arrowDirection);

	vtkSmartPointer<vtkArrowSource> arrowTip = vtkSmartPointer<vtkArrowSource>::New();
	arrowTip->SetTipLength(arrowTipLength);
	arrowTip->SetTipRadius(tubeRadius * 2);
	arrowTip->SetShaftRadius(tubeRadius);

	// 创建箭头尖端的变换
	vtkSmartPointer<vtkTransform> arrowTransform = vtkSmartPointer<vtkTransform>::New();
	arrowTransform->Translate(arrowTipPosition);
	if (axis == 0)
	{
		arrowTransform->RotateWXYZ(vtkMath::DegreesFromRadians(90.0/180*vtkMath::Pi()), 0, -1, 0);
		// arrowTransform->RotateWXYZ(vtkMath::DegreesFromRadians(90.0 / 180 * vtkMath::Pi()), -1, 0, 0);
	}
	else if (axis == 1)
	{
		// arrowTransform->RotateWXYZ(vtkMath::DegreesFromRadians(atan2(-arrowDirection[2], arrowDirection[0])), 0, 1, 0);
	}
	else
	{
		// arrowTransform->RotateWXYZ(vtkMath::DegreesFromRadians(atan2(arrowDirection[1], arrowDirection[0])), 0, 0, 1);
	}
	arrowTransform->Scale(0.3, 1.0, 1.0);

	vtkSmartPointer<vtkTransformPolyDataFilter> arrowTransformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	arrowTransformFilter->SetTransform(arrowTransform);
	arrowTransformFilter->SetInputConnection(arrowTip->GetOutputPort());

	vtkSmartPointer<vtkPolyDataMapper> arrowTipMapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	arrowTipMapper->SetInputConnection(arrowTransformFilter->GetOutputPort());

	vtkSmartPointer<vtkActor> arrowTipActor = vtkSmartPointer<vtkActor>::New();
	arrowTipActor->SetMapper(arrowTipMapper);
	arrowTipActor->GetProperty()->SetColor(arrowColor[0], arrowColor[1], arrowColor[2]);

	// 添加箭头尖端到渲染器
	renderer->AddActor(arrowTipActor);
}



// 辅助函数，计算两个向量的叉积
void RenderUsingVTK::crossProduct(const double a[3], const double b[3], double result[3]) {
	result[0] = a[1] * b[2] - a[2] * b[1];
	result[1] = a[2] * b[0] - a[0] * b[2];
	result[2] = a[0] * b[1] - a[1] * b[0];
}

void RenderUsingVTK::drawLinearArrow(double start[3], double end[3], double color[3])
{
	// 创建箭头源
	vtkSmartPointer<vtkArrowSource> arrowSource = vtkSmartPointer<vtkArrowSource>::New();

	// 计算方向向量
	double direction[3];
	for (int i = 0; i < 3; ++i) {
		direction[i] = end[i] - start[i];
	}

	// 计算方向向量的长度
	double length = sqrt(direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2]);

	if (length < 1e-6) {
		// 如果起点和终点太接近，不绘制箭头
		return;
	}

	// 创建变换
	vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
	transform->Translate(start); // 平移到起点

	// 计算旋转到目标方向的角度和轴
	double zAxis[3] = { 1.0, 0.0, 0.0 }; // 默认箭头方向为 Z 轴
	double rotationAxis[3];
	vtkMath::Cross(zAxis, direction, rotationAxis); // 计算旋转轴
	double rotationAngle = vtkMath::DegreesFromRadians(acos(vtkMath::Dot(zAxis, direction) / length)); // 计算旋转角度

	if (vtkMath::Norm(rotationAxis) > 1e-6) {
		// 如果旋转轴非零，应用旋转
		vtkMath::Normalize(rotationAxis);
		transform->RotateWXYZ(rotationAngle, rotationAxis);
	}

	// 缩放到目标长度
	transform->Scale(length * 0.6, length, length);

	// 应用变换
	vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	transformFilter->SetTransform(transform);
	transformFilter->SetInputConnection(arrowSource->GetOutputPort());
	transformFilter->Update();

	// 映射器
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(transformFilter->GetOutputPort());

	// Actor
	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetColor(color);

	// 添加到渲染器
	renderer->AddActor(actor);
}
