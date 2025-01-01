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
	renderer->SetGradientBackground(true);  // ���ý��䱳��ɫ
	// ���ý��䱳��ɫ����ʼ��ɫ�ͽ�����ɫ
	renderer->SetBackground(0.31, 0.31, 0.51);  // ��ʼ��ɫ��ǳɫ��
	renderer->SetBackground2(0.949, 0.949, 0.949); // ������ɫ�����ɫ��
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
	//���ý������ٿ�ģʽ
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
	//double xColor[3] = { 1.0, 0.0, 0.0 }; // ��ɫ
	//double yColor[3] = { 0.0, 1.0, 0.0 }; // ��ɫ
	//double zColor[3] = { 0.0, 0.0, 1.0 }; // ��ɫ
	//drawLinearArrow(xOrigin, xEnd, xColor);
	//drawLinearArrow(yOrigin, yEnd, yColor);
	//drawLinearArrow(zOrigin, zEnd, zColor);

	//double xPos[3] = { 1.0, 0.0, 0.0 }; // ��ɫ��ͷ
	//double yPos[3] = { 0.0, 1.0, 0.0 }; // ��ɫ��ͷ
	//double zPos[3] = { 0.0, 0.0, 1.0 }; // ��ɫ��ͷ
	//DrawCurvedArrow(0.38, 0, vtkMath::Pi(), 20, 0.02, 0.3, xColor, 0, xPos);
	//DrawCurvedArrow(0.38, 0, vtkMath::Pi(), 20, 0.02, 0.3, yColor, 1, yPos);
	//DrawCurvedArrow(0.38, 0, vtkMath::Pi(), 20, 0.02, 0.3, zColor, 2, zPos);

	//double springBegin[3] = { 0.0, 0.0, 0.0 }; // ��ɫ��ͷ
	//double springXEnd[3] = { 1.0, 0.0, 0.0 };  // ��ɫ��ͷ
	//double springYEnd[3] = { 0.0, 1.0, 0.0 };  // ��ɫ��ͷ
	//double springZEnd[3] = { 0.0, 0.0, 1.0 };  // ��ɫ��ͷ
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
			// x ����
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

			// z ����
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
		groundMapper->SetScalarModeToUseCellData(); // ȷ��ʹ�õ�Ԫ������
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
				//actor->Delete();  // ����VTK�����Delete�������ͷ��ڴ�
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
			// ����͸����
			property->SetOpacity(data.opacity);

			// ���õ��С
			property->SetPointSize(8);

			// ���ò��������Ա������Ȳ���
			// ������ϵ����ʹ�����ڲ�ͬ�ӽ��¸�һ��
			property->SetAmbient(0.3);

			// ������ϵ�������ƹ��յ�ɢ��Ч��
			property->SetDiffuse(0.6);

			// ���淴��ϵ�������Ƹ߹�Ч��
			property->SetSpecular(0.3);

			// ���淴��ָ�������Ƹ߹�ļ��г̶�
			property->SetSpecularPower(60);


			// ����˫����Ⱦ��ʹģ�Ϳ���������Ȼ
			property->SetBackfaceCulling(0);

			renderer->AddActor(this->partActor[i]);
		}
	}

	// ���ù���
	vtkSmartPointer<vtkLight> light = vtkSmartPointer<vtkLight>::New();
	light->SetPositional(1); // ����Ϊλ�ù�
	light->SetPosition(renderer->GetActiveCamera()->GetPosition()); // ��Դλ�������λ����ͬ
	light->SetFocalPoint(renderer->GetActiveCamera()->GetFocalPoint()); // ��Դ���������������ͬ
	light->SetColor(1.0, 1.0, 1.0); // ��ɫ��
	light->SetIntensity(0.5); // ����ǿ��
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
			maperTem->SetScalarVisibility(0); // ���ñ���ӳ��

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
		// ��Χ�еĳߴ� (�Խ��߳���)
		modelSize = sqrt(pow(bounds[1] - bounds[0], 2) +
			pow(bounds[3] - bounds[2], 2) +
			pow(bounds[5] - bounds[4], 2));
	}
	
	auto camera = this->renderer->GetActiveCamera();
	// ��ȡ����ӽǡ������ԭʼλ��
	double viewAngle = camera->GetViewAngle(); // ��ǰ�ӽ�
	double distance = modelSize / (2 * tan(vtkMath::RadiansFromDegrees(viewAngle / 2))); // �Ƽ�����
	double focalPoint[3];
	camera->GetFocalPoint(focalPoint);

	// �����ӽ�
	if (view == sceneView::Main) {
		// ����ͼ������ǰ���� (0, 0, Z) �۲�
		// �������λ�ã�ȷ��ģ�͸պ������Ұ
		camera->SetPosition(focalPoint[0], focalPoint[1], focalPoint[2] + distance);
		camera->SetViewUp(0, 1, 0); // Y��������Ϊ��

		// ��������������Ұ
		renderer->ResetCameraClippingRange();
	}
	else if (view == sceneView::Up) {
		// ����ͼ�������Ϸ��� (0, Y, 0) �۲�
		camera->SetPosition(focalPoint[0], focalPoint[1] + distance, focalPoint[2]);
		camera->SetViewUp(0, 0, -1); // Z�Ḻ����Ϊ��
		// camera->SetParallelProjection(true);
	}
	else if (view == sceneView::Down) {
		// ����ͼ�������·��� (0, -Y, 0) �۲�
		camera->SetPosition(focalPoint[0], focalPoint[1] - distance, focalPoint[2]);
		camera->SetViewUp(0, 0, 1); // Z��������Ϊ��
	}
	else if (view == sceneView::Left) {
		// ����ͼ������෽�� (-X, 0, 0) �۲�
		camera->SetPosition(focalPoint[0] - distance, focalPoint[1], focalPoint[2]);
		camera->SetViewUp(0, 1, 0); // Y��������Ϊ��
	}
	else if (view == sceneView::Right) {
		// ����ͼ�����Ҳ෽�� (X, 0, 0) �۲�
		camera->SetPosition(focalPoint[0] + distance, focalPoint[1], focalPoint[2]);
		camera->SetViewUp(0, 1, 0); // Y��������Ϊ��
	}
	else if (view == sceneView::Suitable) 
	{
		double center[3];
		// ����ģ�����ĵ�
		if (!autoSceneSize)
		{
			renderer->ComputeVisiblePropBounds(bounds);
		}
		center[0] = (bounds[0] + bounds[1]) / 2.0;
		center[1] = (bounds[2] + bounds[3]) / 2.0;
		center[2] = (bounds[4] + bounds[5]) / 2.0;
		// ����Ӧ��ͼ���������λ�úͲü���Χ��ʹ������ȫ�����ӿ�
		// ȷ�������ģ�����ĵľ���
		double distance = 0.5 * modelSize; // ���Ը�����Ҫ��������ϵ��
		// �������λ�ã�ʹ��λ��ģ�����ĵ�һ�����봦��������΢ƫ��
		camera->SetPosition(center[0], center[1], center[2] + distance);
		// �����������Ϊģ������
		camera->SetFocalPoint(center[0], center[1], center[2]);
		// ������������Ϸ���ͨ����Y��������
		camera->SetViewUp(0, 1, 0);
		// ����������ӽǷ�Χ��ȷ��ģ�������ͼ
		renderer->ResetCameraClippingRange();
	}
	else if (view == sceneView::rotX)
	{
		double theta = -45.0; // ������ת�Ƕ�
		rotateViewByXAxis(vtkMath::RadiansFromDegrees(theta), focalPoint);
	}
	else if (view == sceneView::rotY)
	{
		double theta = -45.0; // ������ת�Ƕ�
		rotateViewByYAxis(vtkMath::RadiansFromDegrees(theta), focalPoint);
	}
	else if (view == sceneView::rotZ)
	{
		double theta = -45.0; // ������ת�Ƕ�
		rotateViewByZAxis(vtkMath::RadiansFromDegrees(theta), focalPoint);
	}

	camera->SetParallelProjection(paraView);
	// ��Ⱦ����
	this->renderer->ResetCameraClippingRange(); // ���²ü���Χ
	this->renderer->GetRenderWindow()->Render();
}

void RenderUsingVTK::rotateViewByXAxis(double angleRadians, double focalPoint[3])
{
	vtkCamera* camera = this->renderer->GetActiveCamera();

	// ��ȡ��ǰ�������
	double cameraPosition[3];
	double viewUp[3];
	camera->GetPosition(cameraPosition);
	camera->GetViewUp(viewUp);

	// �������λ����Խ��������
	double dx = cameraPosition[0] - focalPoint[0];
	double dy = cameraPosition[1] - focalPoint[1];
	double dz = cameraPosition[2] - focalPoint[2];

	// ������ת����
	double cosTheta = cos(angleRadians);
	double sinTheta = sin(angleRadians);

	// ��ת������λ��
	double newCameraPosition[3];
	newCameraPosition[0] = focalPoint[0] + dx;
	newCameraPosition[1] = focalPoint[1] + (dy * cosTheta - dz * sinTheta);
	newCameraPosition[2] = focalPoint[2] + (dy * sinTheta + dz * cosTheta);

	// ��ת��ġ��Ϸ���
	double newViewUp[3];
	newViewUp[0] = viewUp[0];
	newViewUp[1] = viewUp[1] * cosTheta - viewUp[2] * sinTheta;
	newViewUp[2] = viewUp[1] * sinTheta + viewUp[2] * cosTheta;

	// �����µ��������
	camera->SetPosition(newCameraPosition);
	camera->SetViewUp(newViewUp);
}

void RenderUsingVTK::rotateViewByYAxis(double angleRadians, double focalPoint[3])
{
	vtkCamera* camera = this->renderer->GetActiveCamera();

	// ��ȡ��ǰ�������
	double cameraPosition[3];
	double viewUp[3];
	camera->GetPosition(cameraPosition);
	camera->GetViewUp(viewUp);

	// �������λ����Խ��������
	double dx = cameraPosition[0] - focalPoint[0];
	double dy = cameraPosition[1] - focalPoint[1];
	double dz = cameraPosition[2] - focalPoint[2];

	// ������ת����
	double cosTheta = cos(angleRadians);
	double sinTheta = sin(angleRadians);

	// ��ת������λ��
	double newCameraPosition[3];
	newCameraPosition[0] = focalPoint[0] + (dx * cosTheta + dz * sinTheta);
	newCameraPosition[1] = focalPoint[1] + dy;
	newCameraPosition[2] = focalPoint[2] + (-dx * sinTheta + dz * cosTheta);

	// ��ת��ġ��Ϸ���
	double newViewUp[3];
	newViewUp[0] = viewUp[0] * cosTheta + viewUp[2] * sinTheta;
	newViewUp[1] = viewUp[1];
	newViewUp[2] = -viewUp[0] * sinTheta + viewUp[2] * cosTheta;

	// �����µ��������
	camera->SetPosition(newCameraPosition);
	camera->SetViewUp(newViewUp);
}

void RenderUsingVTK::rotateViewByZAxis(double angleRadians, double focalPoint[3])
{
	vtkCamera* camera = this->renderer->GetActiveCamera();
	// ��ȡ��ǰ�������
	double cameraPosition[3];
	double viewUp[3];
	camera->GetPosition(cameraPosition);
	camera->GetViewUp(viewUp);

	// �������λ����Խ��������
	double dx = cameraPosition[0] - focalPoint[0];
	double dy = cameraPosition[1] - focalPoint[1];
	double dz = cameraPosition[2] - focalPoint[2];

	// ������ת����
	double cosTheta = cos(angleRadians);
	double sinTheta = sin(angleRadians);

	// ��ת������λ��
	double newCameraPosition[3];
	newCameraPosition[0] = focalPoint[0] + (dx * cosTheta - dy * sinTheta);
	newCameraPosition[1] = focalPoint[1] + (dx * sinTheta + dy * cosTheta);
	newCameraPosition[2] = focalPoint[2] + dz;

	// ��ת��ġ��Ϸ���
	double newViewUp[3];
	newViewUp[0] = viewUp[0] * cosTheta - viewUp[1] * sinTheta;
	newViewUp[1] = viewUp[0] * sinTheta + viewUp[1] * cosTheta;
	newViewUp[2] = viewUp[2];

	// �����µ��������
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
		// ��������ģʽ��16λģʽ������ʹ��0xF0F0��Ϊʾ��
		//dynamicActor->GetProperty()->SetLineStipplePattern(0xFFFF);
		//dynamicActor->GetProperty()->SetLineStippleRepeatFactor(1); // �����ظ�����
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
			grid->Delete();  // ����VTK�����Delete�������ͷ��ڴ�
		}
	}
	ugModel.clear();

	for (auto& part : data.geoModel.parts)
	{
		ugModel[part.first] = vtkUnstructuredGrid::New();

		// ��ɫ
		if (!ugModel[part.first]->GetCellData()->HasArray("Color"))
		{
			vtkSmartPointer<vtkUnsignedCharArray> cellColors = vtkSmartPointer<vtkUnsignedCharArray>::New();
			cellColors->SetNumberOfComponents(3);  // RGB
			cellColors->SetName("Color");
			ugModel[part.first]->GetCellData()->SetScalars(cellColors);
		}

		// ����
		if (!ugModel[part.first]->GetCellData()->GetNormals())
		{
			vtkSmartPointer<vtkDoubleArray> cellNormals = vtkSmartPointer<vtkDoubleArray>::New();
			cellNormals->SetNumberOfComponents(3);  // ������Ϊ 3D ����
			cellNormals->SetName("Normal");
			ugModel[part.first]->GetCellData()->SetNormals(cellNormals);
		}
	}

	// ���� ug ����һ�׽ڵ�
	vtkPoints* points = vtkPoints::New();

	// ��ɫ������
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

	double rgbColor[3]; // ���ڴ洢RGB��ɫ������

	ncs::GeometryModel* md = &data.geoModel;
	int ndID = 0;
	// ������� vertex��line��tri��quad��color��normal
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

			// ����ɫ�󶨵� cellColors
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
			// �����߶εļ�����Ϣ
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

			// ��ӵ�Ԫ��ɫ������ cell��
			cellColors->InsertNextTuple3(data.lineColor[i].r,
				data.lineColor[i].g,
				data.lineColor[i].b);

			// ��ӷ���������ѡ��
			cellNormals->InsertNextTuple3(0, 1, 0); // �Ե�1�ķ�����
			cellNormals->InsertNextTuple3(0, 1, 0); // �Ե�2�ķ�����

			// �����߶εĵ�
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

			//�������ݵ�vtk���ݽṹ��
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
			//��������
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

			// ���㷨����
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

			// ���õ�Ԫ��ɫ�ͷ�����
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

	// ����points����ɫ���ݵ�vtk���ݽṹ��
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

	////����points����ɫ���ݵ�vtk���ݽṹ��
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
		cellNormals->SetNumberOfComponents(3);  // ������Ϊ 3D ����
		cellNormals->SetName("Normal");
		ug->GetCellData()->SetNormals(cellNormals);
	}
	vtkDoubleArray* cellNormals = static_cast<vtkDoubleArray*>(ug->GetCellData()->GetNormals());

	// ��ɫ������
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
	double rgbColor[3]; // ���ڴ洢RGB��ɫ������

	int ndID = points->GetNumberOfPoints();
	// ������� vertex��line��tri��quad��color��normal
	if (data.vertexColor.size() == md->vertex.size())
	{
		for (int i = 0; i < md->vertex.size(); i++)
		{
			vtkSmartPointer<vtkVertex> vertex = vtkSmartPointer<vtkVertex>::New();
			points->InsertNextPoint(md->vertex[i].x, md->vertex[i].y, md->vertex[i].z);

			// ����ɫ�󶨵� cellColors ������ pointColors
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
			// �����߶εļ�����Ϣ
			vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
			points->InsertNextPoint(md->nodes[md->lines[i].nd1].x,
				md->nodes[md->lines[i].nd1].y,
				md->nodes[md->lines[i].nd1].z);
			points->InsertNextPoint(md->nodes[md->lines[i].nd2].x,
				md->nodes[md->lines[i].nd2].y,
				md->nodes[md->lines[i].nd2].z);

			// ��ӵ�Ԫ��ɫ������ cell��
			cellColors->InsertNextTuple3(data.lineColor[i].r,
				data.lineColor[i].g,
				data.lineColor[i].b);

			// ��ӷ���������ѡ��
			cellNormals->InsertNextTuple3(0, 1, 0); // �Ե�1�ķ�����
			cellNormals->InsertNextTuple3(0, 1, 0); // �Ե�2�ķ�����

			// �����߶εĵ�
			line->GetPointIds()->SetId(0, ndID);
			line->GetPointIds()->SetId(1, ndID + 1);
			ndID += 2;

			// ���߶μ��뵽��Ԫ����
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
			//�������ݵ�vtk���ݽṹ��
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
			//��������
			cells->InsertNextCell(triangle);
			cellTypes->InsertNextTuple1(VTK_TRIANGLE);
		}
	}

	if (quadsForVtk.size() == data.quadColor.size())
	{
		int quadID = 0;
		for (const auto& qd : quadsForVtk) {
			vtkSmartPointer<vtkQuad> quad = vtkSmartPointer<vtkQuad>::New();

			// ���㷨����
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

			// ���õ�Ԫ��ɫ�ͷ�����
			cellColors->InsertNextTuple3(data.quadColor[quadID].r, data.quadColor[quadID].g,
				data.quadColor[quadID].b);
			cellNormals->InsertNextTuple3(nml.x, nml.y, nml.z);

			// ���뵥Ԫ
			cells->InsertNextCell(quad);
			cellTypes->InsertNextTuple1(VTK_QUAD);
			quadID++;
		}

	}

	//��װģ��
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

	// ����points����ɫ���ݵ�vtk���ݽṹ��
	if (!ug->GetPointData()->HasArray("Color"))
	{
		vtkSmartPointer<vtkUnsignedCharArray> colors = vtkSmartPointer<vtkUnsignedCharArray>::New();
		colors->SetNumberOfComponents(3);
		colors->SetName("Color");
		ug->GetPointData()->AddArray(colors);
	}
	vtkUnsignedCharArray* colors = static_cast<vtkUnsignedCharArray*>(ug->GetPointData()->GetArray("Color"));

	//����points����ɫ���ݵ�vtk���ݽṹ��
	if (ug->GetPointData()->GetNormals() == nullptr)
	{
		vtkSmartPointer<vtkDoubleArray> normals = vtkSmartPointer<vtkDoubleArray>::New();
		normals->SetNumberOfComponents(3);
		normals->SetName("Normal");
		ug->GetPointData()->SetNormals(normals);
	}
	vtkDoubleArray* normals = static_cast<vtkDoubleArray*>(ug->GetPointData()->GetNormals());

	// ��ɫ������
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
	double rgbColor[3]; // ���ڴ洢RGB��ɫ������

	int ndID = points->GetNumberOfPoints();
	// ������� vertex��line��tri��quad��color��normal
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
			//�������ݵ�vtk���ݽṹ��
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
			//��������
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
			//�������ݵ�vtk���ݽṹ��
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
			//��������
			cells->InsertNextCell(triangle);
			cellTypes->InsertNextTuple1(VTK_TRIANGLE);
		}
	}
	
	if (quadsForVtk.size() == data.quadColor.size())
	{
		int quadID = 0;
		for (const auto& qd : quadsForVtk) {
			//�������ݵ�vtk���ݽṹ��
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
			//��������
			cells->InsertNextCell(quad);
			cellTypes->InsertNextTuple1(VTK_QUAD);
		}
	}
	
	//��װģ��
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
	// ��ɫ������
	setColorTransferFunction(renderInfo);

	contourMapper->SetLookupTable(colorTransferFunction);
	// ������Ա
	contourActor->SetMapper(contourMapper);
	// �����Ա����Ⱦ��
	// renderer->AddActor(contourActor);
	// ����ɫ��
	scalarBar->SetLookupTable(colorTransferFunction);

	scalarBar->SetNumberOfLabels(renderInfo.colorNum);

# if 1
	// ˮƽɫ��
	scalarBar.Get()->SetOrientationToHorizontal();
	scalarBar.Get()->SetPosition(0.25, 0.1);
	scalarBar->SetWidth(0.5);
	scalarBar->SetHeight(0.08);

	// ��ȡ scalarBar ��λ�úʹ�С
	double* scalarBarPos = scalarBar->GetPosition();       // ���½�λ��
	double* scalarBarSize = scalarBar->GetPosition2();     // ��Ⱥ͸߶�

	// �����Զ����ǩ��λ��
	// ���磺���Զ����ǩ��������ɫ�������·�
	double labelX = scalarBarPos[0] + scalarBarSize[0] / 2.0;   // ˮƽ����
	double labelY = scalarBarPos[1] - scalarBarSize[1] * 0.8;       // ����ɫ�����·�����΢����ƫ��

	// �����������ı���ǩ
	// if(titleActor = nullptr) titleActor = vtkSmartPointer<vtkTextActor>::New();
	
	titleActor->SetInput((renderInfo.colorBarName + '\n').c_str());  // ���ñ�������
	titleActor->GetTextProperty()->SetFontSize(15);
	titleActor->GetTextProperty()->SetColor(0.0, 0.0, 0.0);  // ����������ɫ
	titleActor->GetTextProperty()->SetFontFamilyToTimes();
	titleActor->GetTextProperty()->SetItalic(1);  // ��������Ϊб��
	titleActor->GetTextProperty()->SetBold(1);
	// ���ñ����λ�õ�ɫ���·�
	renderer->NormalizedDisplayToDisplay(labelX, labelY);
	titleActor->SetPosition(labelX, labelY);  // ���ݾ�����Ⱦ���ڵ��������
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
	// ��ֱɫ��
	scalarBar->SetTitle((renderInfo.colorBarName + '\n').c_str());
	scalarBar->SetWidth(0.06);
	scalarBar->SetHeight(1.0);
	vtkTextProperty* textProperty = scalarBar->GetTitleTextProperty();
	// textProperty->SetLineOffset(12);
	textProperty->SetFontSize(15); // ���������С
	textProperty->SetColor(0.0, 0.0, 0.0); // ����������ɫ
	textProperty->SetFontFamilyToTimes();
	textProperty->SetVerticalJustificationToBottom();
	textProperty->SetLineOffset(200);
#endif
	vtkTextProperty* contentProperty = scalarBar->GetLabelTextProperty();
	contentProperty->SetFontSize(12);  // ���������С
	contentProperty->SetColor(0.0, 0.0, 0.0); // ����������ɫ
	contentProperty->SetFontFamilyToTimes();
	contentProperty->SetJustificationToRight();  // ���ñ�ǩ��ɫ���·�

	scalarBar->SetLabelFormat("%.3f"); // ������Чλ��
	scalarBar->Modified();
	scalarBar->SetUnconstrainedFontSize(true);
	
}

void RenderUsingVTK::setColorTransferFunction(RenderBuffer& renderInfo)
{
	// ��ɫ������
	if (colorTransferFunction == nullptr)
	{
		colorTransferFunction = vtkSmartPointer<vtkColorTransferFunction>::New();

	}
	colorTransferFunction->RemoveAllPoints();
	double dValue = renderInfo.maxValue - renderInfo.minValue;

	colorTransferFunction->RemoveAllPoints();
	colorTransferFunction->AddRGBPoint(renderInfo.minValue, 0.0, 0.0, 1.0);  // ��ɫ
	colorTransferFunction->AddRGBPoint(renderInfo.minValue + 0.2 * dValue, 0.0, 1.0, 1.0);  // ��ɫ
	colorTransferFunction->AddRGBPoint(renderInfo.minValue + 0.4 * dValue, 0.0, 1.0, 0.0);  // ��ɫ
	colorTransferFunction->AddRGBPoint(renderInfo.minValue + 0.6 * dValue, 1.0, 1.0, 0.0);  // ��ɫ
	colorTransferFunction->AddRGBPoint(renderInfo.minValue + 0.8 * dValue, 1.0, 0.5, 0.0);  // ��ɫ
	colorTransferFunction->AddRGBPoint(renderInfo.maxValue, 1.0, 0.0, 0.0);  // ��ɫ
}

void RenderUsingVTK::UpdateTitlePosition(vtkObject* caller, unsigned long eventId, void* callData)
{
	// ��ȡ scalarBar ��λ�úʹ�С
	double* scalarBarPos = scalarBar->GetPosition();       // ���½�λ��
	double* scalarBarSize = scalarBar->GetPosition2();     // ��Ⱥ͸߶�

	// �����Զ����ǩ��λ��
	// ���磺���Զ����ǩ��������ɫ�������·�
	double labelX = scalarBarPos[0] + scalarBarSize[0] / 2.0;   // ˮƽ����
	double labelY = scalarBarPos[1] - scalarBarSize[1] * 0.8;       // ����ɫ�����·�����΢����ƫ��

	// ���ñ����λ�õ�ɫ���·�
	renderer->NormalizedDisplayToDisplay(labelX, labelY);
	titleActor->SetPosition(labelX, labelY);  // ���ݾ�����Ⱦ���ڵ��������
	//std::cout << "adjust: labrlX: " << labelX << "    labrlY: " << labelY << '\n';
}


void RenderUsingVTK::showAxies()
{
#if 0
	vtkAxesActor* axes_actor = vtkAxesActor::New();
	axes_actor->SetAxisLabels(0);
	axes_actor->SetTotalLength(0.2, 0.2, 0.2);
	//// ��ȡ���������Ƶ��ı�����
	//vtkTextProperty* labelXProperty = axes_actor->GetXAxisCaptionActor2D()->GetCaptionTextProperty();
	//vtkTextProperty* labelYProperty = axes_actor->GetYAxisCaptionActor2D()->GetCaptionTextProperty();
	//vtkTextProperty* labelZProperty = axes_actor->GetZAxisCaptionActor2D()->GetCaptionTextProperty();

	//// labelXProperty->SetFontFamilyToArial(); // ����������壨��ѡ��
	//labelXProperty->SetFontSize(1);    // ���������С
	//labelYProperty->SetFontSize(1);    // ���������С
	//labelZProperty->SetFontSize(1);    // ���������С

	//// ����������ɫ
	//double textColor[3] = { 0.0, 0.0, 0.0 }; // ��ɫ
	//labelXProperty->SetColor(textColor);
	//labelYProperty->SetColor(textColor);
	//labelZProperty->SetColor(textColor);

	//// ��������Ϊ Times New Roman
	//labelXProperty->SetFontFamily(VTK_TIMES);
	//labelYProperty->SetFontFamily(VTK_TIMES);
	//labelZProperty->SetFontFamily(VTK_TIMES);

	//// ����Ϊб�壨ʾ����
	//labelXProperty->ItalicOn();
	//labelYProperty->ItalicOn();
	//labelZProperty->ItalicOn();

	//vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
	//transform->Translate(10, 10, 10); // �����µ�ԭ��λ��

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
	//// ���������˵�
	//double point1[3] = { 0.0, 0.0, 0.0 }; // ���
	//double point2[3] = { 1.0, 1.0, 1.0 }; // �յ�

	// ���ɲ���
	double springRadius = 0.1;    // ���ɰ뾶
	int numberOfTurns = 12;        // ����Ȧ��
	int numPointsPerTurn = 100;   // ÿȦ�ĵ���
	double lineThickness = 0.01;
	double springHeight = vtkMath::Distance2BetweenPoints(point1, point2); // ���ɵ��ܸ߶�
	double pitch = springHeight / numberOfTurns; // ÿȦ���ݾ�

	// ���㷽������
	double direction[3];
	vtkMath::Subtract(point2, point1, direction);
	vtkMath::Normalize(direction); // ��һ����������

	// ����һ������������֤�ֲ����������ƽ���뷽����������
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

	// ����������
	double tangent[3];
	vtkMath::Cross(arbitraryVector, direction, tangent);
	vtkMath::Normalize(tangent);

	double binormal[3];
	vtkMath::Cross(direction, tangent, binormal);
	vtkMath::Normalize(binormal);

	// ���������ߵĵ㼯
	vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();

	for (int i = 0; i <= numberOfTurns * numPointsPerTurn; ++i) {
		double theta = 2 * vtkMath::Pi() * i / numPointsPerTurn; // ��ǰ��ĽǶ�
		double localX = springRadius * cos(theta);              // ��������ֲ� x
		double localY = springRadius * sin(theta);              // ��������ֲ� y
		double localZ = pitch * i / numPointsPerTurn;           // �ֲ� z���ط���������ƫ�ƣ�

		// ���ֲ�����ת����ȫ������
		double globalPoint[3] = {
			point1[0] + tangent[0] * localX + binormal[0] * localY + direction[0] * localZ,
			point1[1] + tangent[1] * localX + binormal[1] * localY + direction[1] * localZ,
			point1[2] + tangent[2] * localX + binormal[2] * localY + direction[2] * localZ
		};
		points->InsertNextPoint(globalPoint);
	}

	// ���������ߵĶ��������
	vtkSmartPointer<vtkPolyLine> polyLine = vtkSmartPointer<vtkPolyLine>::New();
	int totalPoints = numberOfTurns * numPointsPerTurn + 1;
	polyLine->GetPointIds()->SetNumberOfIds(totalPoints);
	for (int i = 0; i < totalPoints; ++i) {
		polyLine->GetPointIds()->SetId(i, i);
	}

	vtkSmartPointer<vtkCellArray> lines = vtkSmartPointer<vtkCellArray>::New();
	lines->InsertNextCell(polyLine);

	// ��������������
	vtkSmartPointer<vtkPolyData> helixData = vtkSmartPointer<vtkPolyData>::New();
	helixData->SetPoints(points);
	helixData->SetLines(lines);

	// ��������ת��Ϊ��״�ṹ
	vtkSmartPointer<vtkTubeFilter> tubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
	tubeFilter->SetInputData(helixData);
	tubeFilter->SetRadius(lineThickness); // ��״�ṹ�İ뾶
	tubeFilter->SetNumberOfSides(10); // ��״�ṹ�Ĳ�����

	// ӳ�����Ⱦ
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(tubeFilter->GetOutputPort());

	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetColor(color[0], color[1], color[2]); // ������ɫ

	renderer->AddActor(actor);
}

void RenderUsingVTK::DrawCurvedArrow(
	double radius,         // �����İ뾶
	double startAngle,     // ��ʼ�Ƕȣ����ȣ�
	double endAngle,       // �����Ƕȣ����ȣ�
	int numPoints,         // ���ߵ�����
	double tubeRadius,     // �ܾ�
	double arrowTipLength, // ��ͷ����
	double arrowColor[3],
	int axis,
	double pos[3])  // ��ͷ��ɫ (RGB)
{
	// ����������ͷ�ĵ�
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

	// ������ͷ������
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

	// ת��Ϊ��״�ṹ
	vtkSmartPointer<vtkTubeFilter> tubeFilter = vtkSmartPointer<vtkTubeFilter>::New();
	tubeFilter->SetInputData(curveData);
	tubeFilter->SetRadius(tubeRadius);
	tubeFilter->SetNumberOfSides(20);

	// ӳ�䲢�������߼�ͷ���ֵ�Actor
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(tubeFilter->GetOutputPort());

	vtkSmartPointer<vtkActor> arrowBodyActor = vtkSmartPointer<vtkActor>::New();
	arrowBodyActor->SetMapper(mapper);
	arrowBodyActor->GetProperty()->SetColor(arrowColor[0], arrowColor[1], arrowColor[2]);

	// ��Ӽ�ͷ��β�����ֵ���Ⱦ��
	renderer->AddActor(arrowBodyActor);

	// ������ͷ��ͷ������ˣ�
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

	// ������ͷ��˵ı任
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

	// ��Ӽ�ͷ��˵���Ⱦ��
	renderer->AddActor(arrowTipActor);
}



// �����������������������Ĳ��
void RenderUsingVTK::crossProduct(const double a[3], const double b[3], double result[3]) {
	result[0] = a[1] * b[2] - a[2] * b[1];
	result[1] = a[2] * b[0] - a[0] * b[2];
	result[2] = a[0] * b[1] - a[1] * b[0];
}

void RenderUsingVTK::drawLinearArrow(double start[3], double end[3], double color[3])
{
	// ������ͷԴ
	vtkSmartPointer<vtkArrowSource> arrowSource = vtkSmartPointer<vtkArrowSource>::New();

	// ���㷽������
	double direction[3];
	for (int i = 0; i < 3; ++i) {
		direction[i] = end[i] - start[i];
	}

	// ���㷽�������ĳ���
	double length = sqrt(direction[0] * direction[0] + direction[1] * direction[1] + direction[2] * direction[2]);

	if (length < 1e-6) {
		// ��������յ�̫�ӽ��������Ƽ�ͷ
		return;
	}

	// �����任
	vtkSmartPointer<vtkTransform> transform = vtkSmartPointer<vtkTransform>::New();
	transform->Translate(start); // ƽ�Ƶ����

	// ������ת��Ŀ�귽��ĽǶȺ���
	double zAxis[3] = { 1.0, 0.0, 0.0 }; // Ĭ�ϼ�ͷ����Ϊ Z ��
	double rotationAxis[3];
	vtkMath::Cross(zAxis, direction, rotationAxis); // ������ת��
	double rotationAngle = vtkMath::DegreesFromRadians(acos(vtkMath::Dot(zAxis, direction) / length)); // ������ת�Ƕ�

	if (vtkMath::Norm(rotationAxis) > 1e-6) {
		// �����ת����㣬Ӧ����ת
		vtkMath::Normalize(rotationAxis);
		transform->RotateWXYZ(rotationAngle, rotationAxis);
	}

	// ���ŵ�Ŀ�곤��
	transform->Scale(length * 0.6, length, length);

	// Ӧ�ñ任
	vtkSmartPointer<vtkTransformPolyDataFilter> transformFilter = vtkSmartPointer<vtkTransformPolyDataFilter>::New();
	transformFilter->SetTransform(transform);
	transformFilter->SetInputConnection(arrowSource->GetOutputPort());
	transformFilter->Update();

	// ӳ����
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(transformFilter->GetOutputPort());

	// Actor
	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetColor(color);

	// ��ӵ���Ⱦ��
	renderer->AddActor(actor);
}
