#pragma once
#include "algFoundation.h"

namespace ncs
{
	namespace alg
	{
		ALGBASE_API void calSolidEtVolume(SolidElement& solid, std::vector<Node>& nds);
		ALGBASE_API bool geometryToNcbuffers(GeometryModel& originModel, GeometryModel& model,
			Buffer& buffer);
		ALGBASE_API bool getSolidExtent(std::vector<Node>& nodes, SolidElement& et, 
			Extent& extent);
		ALGBASE_API bool getTriangleExtent(const std::vector<Node>& nodes, const Triangle& tri,
			Extent& extent);

		ALGBASE_API bool pointIsUpPlane(double x, double y, double z,
			ncs::Node n1, ncs::Node n2, ncs::Node n3);

		// HSV转RGB算法
		void hsvToRgb(float h, float s, float v, ncs::Color& color);
		ALGBASE_API void getPartColor(int partID, int totalPart, 
			ncs::Color& color);

		ALGBASE_API void getResultColor(double minV, double maxV, int myValue,
			ncs::Color& color);

		/*
		* 
		* 计算单元集合的质量并将质心移动到原点
		* model 中须包含solid、nodes、mat、part
		* 计算修改model的nodes坐标
		* 
		*/
		ALGBASE_API void calcBodyMassAndInertia(ncs::GeometryModel& model,
			std::unordered_map<int, ncs::Material>& mat,
			double& volume, double& mass, ncs::Vec3& Inertia, double& length,
			double& ratio, ncs::Vec3* massCenPos = nullptr);

		/*
		* 有单元质量，只计算转动惯量
		*/
		ALGBASE_API void calcBodyInertia(const ncs::VNode& nodes,
			ncs::VSolidElement& solids, ncs::Vec3& Inertia);
		/*
		* 有单元质量后，直接计算转动惯量
		*/
		ALGBASE_API void calcInertia(ncs::VNode& nodes, ncs::VSolidElement& solids, ncs::Vec3& Inertia);
		
		/*
		* vector to matrix
		*/
		ALGBASE_API void vectorToMatrix(double v[3], double m[3][3]);


		ALGBASE_API bool lineInterPlane(ncs::Vec3& origin, ncs::Vec3& dir, ncs::Vec3 A, ncs::Vec3 B,
			ncs::Vec3 C, double& dis);
		

	}
	
}