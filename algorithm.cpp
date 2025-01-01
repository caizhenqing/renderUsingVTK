#include "algorithm.h"

void ncs::alg::calSolidEtVolume(SolidElement& solid, std::vector<Node>& nds)
{
	double a, b, c, d, e, f;
	double tmp0, tmp1, tmp2;
	solid.volume = 0;

	// 4,5,7,0
	ncs::Vec3 n0, n1, n2, n3, n4, n5, n6, n7;
	ncs::Node* nd;
	nd = &nds[solid.nds[0]];
	n0 = ncs::Vec3(nd->x, nd->y, nd->z);
	nd = &nds[solid.nds[1]];
	n1 = ncs::Vec3(nd->x, nd->y, nd->z);
	nd = &nds[solid.nds[2]];
	n2 = ncs::Vec3(nd->x, nd->y, nd->z);
	nd = &nds[solid.nds[3]];
	n3 = ncs::Vec3(nd->x, nd->y, nd->z);
	nd = &nds[solid.nds[4]];
	n4 = ncs::Vec3(nd->x, nd->y, nd->z);
	nd = &nds[solid.nds[5]];
	n5 = ncs::Vec3(nd->x, nd->y, nd->z);
	nd = &nds[solid.nds[6]];
	n6 = ncs::Vec3(nd->x, nd->y, nd->z);
	nd = &nds[solid.nds[7]];
	n7 = ncs::Vec3(nd->x, nd->y, nd->z);

	a = (n4 - n5).magnitude();
	b = (n4 - n7).magnitude();
	c = (n4 - n0).magnitude();
	d = (n0 - n7).magnitude();
	e = (n0 - n5).magnitude();
	f = (n5 - n7).magnitude();

	tmp0 = b * b + c * c - d * d;
	tmp1 = a * a + c * c - e * e;
	tmp2 = a * a + b * b - f * f;
	double value;
	value = 4 * a * a * b * b * c * c - a * a * tmp0 * tmp0 - b * b * tmp1 * tmp1
		- c * c * tmp2 * tmp2 + tmp0 * tmp1 * tmp2;
	if (value > 0.0)
	{
		solid.volume += sqrt(value) / 12;
	}
	

	// 1,2,5,0
	a = (n1 - n2).magnitude();
	b = (n1 - n5).magnitude();
	c = (n1 - n0).magnitude();
	d = (n0 - n5).magnitude();
	e = (n0 - n2).magnitude();
	f = (n2 - n5).magnitude();

	tmp0 = b * b + c * c - d * d;
	tmp1 = a * a + c * c - e * e;
	tmp2 = a * a + b * b - f * f;
	// solid.volume += sqrt(4 * a * a * b * b * c * c - a * a * tmp0 * tmp0 - b * b * tmp1 * tmp1 - c * c * tmp2 * tmp2 + tmp0 * tmp1 * tmp2) / 12;
	value = 4 * a * a * b * b * c * c - a * a * tmp0 * tmp0 - b * b * tmp1 * tmp1
		- c * c * tmp2 * tmp2 + tmp0 * tmp1 * tmp2;
	if (value > 0.0)
	{
		solid.volume += sqrt(value) / 12;
	}


	// 6,5,2,7
	a = (n6 - n5).magnitude();
	b = (n6 - n2).magnitude();
	c = (n6 - n7).magnitude();
	d = (n7 - n2).magnitude();
	e = (n7 - n5).magnitude();
	f = (n5 - n2).magnitude();

	tmp0 = b * b + c * c - d * d;
	tmp1 = a * a + c * c - e * e;
	tmp2 = a * a + b * b - f * f;
	// solid.volume += sqrt(4 * a * a * b * b * c * c - a * a * tmp0 * tmp0 - b * b * tmp1 * tmp1 - c * c * tmp2 * tmp2 + tmp0 * tmp1 * tmp2) / 12;
	value = 4 * a * a * b * b * c * c - a * a * tmp0 * tmp0 - b * b * tmp1 * tmp1
		- c * c * tmp2 * tmp2 + tmp0 * tmp1 * tmp2;
	if (value > 0.0)
	{
		solid.volume += sqrt(value) / 12;
	}

	// 3,7,2,0
	a = (n3 - n7).magnitude();
	b = (n3 - n2).magnitude();
	c = (n3 - n0).magnitude();
	d = (n0 - n2).magnitude();
	e = (n0 - n7).magnitude();
	f = (n7 - n2).magnitude();

	tmp0 = b * b + c * c - d * d;
	tmp1 = a * a + c * c - e * e;
	tmp2 = a * a + b * b - f * f;
	// solid.volume += sqrt(4 * a * a * b * b * c * c - a * a * tmp0 * tmp0 - b * b * tmp1 * tmp1 - c * c * tmp2 * tmp2 + tmp0 * tmp1 * tmp2) / 12;
	value = 4 * a * a * b * b * c * c - a * a * tmp0 * tmp0 - b * b * tmp1 * tmp1
		- c * c * tmp2 * tmp2 + tmp0 * tmp1 * tmp2;
	if (value > 0.0)
	{
		solid.volume += sqrt(value) / 12;
	}

	// 5,7,0,2
	a = (n5 - n7).magnitude();
	b = (n5 - n0).magnitude();
	c = (n5 - n2).magnitude();
	d = (n2 - n0).magnitude();
	e = (n2 - n7).magnitude();
	f = (n7 - n0).magnitude();

	tmp0 = b * b + c * c - d * d;
	tmp1 = a * a + c * c - e * e;
	tmp2 = a * a + b * b - f * f;
	// solid.volume += sqrt(4 * a * a * b * b * c * c - a * a * tmp0 * tmp0 - b * b * tmp1 * tmp1 - c * c * tmp2 * tmp2 + tmp0 * tmp1 * tmp2) / 12;
	value = 4 * a * a * b * b * c * c - a * a * tmp0 * tmp0 - b * b * tmp1 * tmp1
		- c * c * tmp2 * tmp2 + tmp0 * tmp1 * tmp2;
	if (value > 0.0)
	{
		solid.volume += sqrt(value) / 12;
	}

}

bool ncs::alg::geometryToNcbuffers(GeometryModel& originModel, GeometryModel& model,
	Buffer& buffer)
{
	if (model.nodes.size() <= 0) return false;

	ncs::Node* nd;
	ncs::Node* ndO;
	ncs::Triangle* tri;
	double vmin, vmax;
	for (int i = 0; i < model.nodes.size(); i++)
	{
		nd = &model.nodes[i];
		ndO = &originModel.nodes[i];
		double v = abs(ndO->y - nd->y);
		if (i == 0)
		{
			vmin = v;
			vmax = v;
			continue;
		}

		if (vmin > v)vmin = v;
		if (vmax < v)vmax = v;
	}

	ncs::Color posyColor;

	for (int i = 0; i < model.triangles.size(); i++)
	{
		/*tri = &(model.triangles)[i];*/
		tri = &model.triangles[i];

		ncs::Vec3 v1, v2, normal;
		ncs::Node* n1, * n2, * n3;
		n1 = &(model.nodes)[tri->nds[0]];
		n2 = &(model.nodes)[tri->nds[1]];
		n3 = &(model.nodes)[tri->nds[2]];
		v1.SetValue(n2->x - n1->x, n2->y - n1->y, n2->z - n1->z);
		v2.SetValue(n3->x - n1->x, n3->y - n1->y, n3->z - n1->z);

		normal = v1.cross(v2);
		normal.normalize();

		for (int j = 0; j < 3; j++)
		{
			nd = &(model.nodes)[tri->nds[j]];
			buffer.nds.emplace_back(std::vector<double>{ nd->x, nd->y, nd->z });
			posyColor.setColor(vmin, vmax, abs(nd->y - originModel.nodes[tri->nds[j]].y));
			buffer.colors.emplace_back(std::vector<int>{posyColor.r, posyColor.g, posyColor.b});
			// buffer.colors.emplace_back(std::vector<int>{255, 0, 0});
			buffer.normals.emplace_back(std::vector<double>{ normal.x, normal.y, normal.z});
		}
		buffer.tris.emplace_back(std::vector<int>{i * 3, i * 3 + 1, i * 3 + 2});
	}
	return false;
}

bool ncs::alg::getSolidExtent(std::vector<Node>& nodes, SolidElement& et, Extent& extent)
{
	for (int i = 0; i < 8; i++)
	{
		if (i == 0)
		{
			extent.xmin = nodes[et.nds[i]].x;
			extent.xmax = nodes[et.nds[i]].x;

			extent.ymin = nodes[et.nds[i]].y;
			extent.ymax = nodes[et.nds[i]].y;

			extent.zmin = nodes[et.nds[i]].z;
			extent.zmax = nodes[et.nds[i]].z;
			continue;
		}

		if (extent.xmin > nodes[et.nds[i]].x)extent.xmin = nodes[et.nds[i]].x;
		if (extent.xmax < nodes[et.nds[i]].x)extent.xmax = nodes[et.nds[i]].x;

		if (extent.ymin > nodes[et.nds[i]].y)extent.ymin = nodes[et.nds[i]].y;
		if (extent.ymax < nodes[et.nds[i]].y)extent.ymax = nodes[et.nds[i]].y;

		if (extent.zmin > nodes[et.nds[i]].z)extent.zmin = nodes[et.nds[i]].z;
		if (extent.zmax < nodes[et.nds[i]].z)extent.zmax = nodes[et.nds[i]].z;
	}
	return true;
}

bool ncs::alg::getTriangleExtent(const std::vector<Node>& nodes, const Triangle& tri, Extent& extent)
{
	if (tri.thickness <= 1e-10) return false;
	unsigned int ndID;
	ndID = tri.nds[0];

	extent.xmin = nodes[ndID].x;
	extent.xmax = nodes[ndID].x;
	extent.ymin = nodes[ndID].y;
	extent.ymax = nodes[ndID].y;
	extent.zmin = nodes[ndID].z;
	extent.zmax = nodes[ndID].z;
	double tem1, tem2;
	for (int i = 0; i < 3; i++)
	{
		ndID = tri.nds[i];
		if (extent.xmin > nodes[ndID].x) extent.xmin = nodes[ndID].x;
		if (extent.xmax < nodes[ndID].x) extent.xmax = nodes[ndID].x;
		if (extent.ymin > nodes[ndID].y) extent.ymin = nodes[ndID].y;
		if (extent.ymax < nodes[ndID].y) extent.ymax = nodes[ndID].y;
		if (extent.zmin > nodes[ndID].z) extent.zmin = nodes[ndID].z;
		if (extent.zmax < nodes[ndID].z) extent.zmax = nodes[ndID].z;

		ncs::Node ndNormal;
		ndNormal.x = nodes[ndID].x + tri.normal.x * tri.thickness;
		ndNormal.y = nodes[ndID].y + tri.normal.y * tri.thickness;
		ndNormal.z = nodes[ndID].z + tri.normal.z * tri.thickness;
		if (extent.xmin > ndNormal.x) extent.xmin = ndNormal.x;
		if (extent.xmax < ndNormal.x) extent.xmax = ndNormal.x;
		if (extent.ymin > ndNormal.y) extent.ymin = ndNormal.y;
		if (extent.ymax < ndNormal.y) extent.ymax = ndNormal.y;
		if (extent.zmin > ndNormal.z) extent.zmin = ndNormal.z;
		if (extent.zmax < ndNormal.z) extent.zmax = ndNormal.z;
	}

	return false;
}

bool ncs::alg::pointIsUpPlane(double x, double y, double z, ncs::Node n1, ncs::Node n2, ncs::Node n3)
{
	ncs::Vec3 v1, v2, v3, normal;
	v1.x = n2.x - n1.x;
	v1.y = n2.y - n1.y;
	v1.z = n2.z - n1.z;

	v2.x = n3.x - n1.x;
	v2.y = n3.y - n1.y;
	v2.z = n3.z - n1.z;

	v3.x = x - n1.x;
	v3.y = y - n1.y;
	v3.z = z - n1.z;

	normal = v1.cross(v2);
	return normal.dot(v3) > 0;
}

// HSV转RGB算法
void ncs::alg::hsvToRgb(float h, float s, float v, ncs::Color& color) {
    int hi = static_cast<int>(std::floor(h / 60)) % 6;
    float f = h / 60 - hi;
    float p = v * (1 - s);
    float q = v * (1 - f * s);
    float t = v * (1 - (1 - f) * s);

    switch (hi) {
        case 0: 
			// 红到黄
			color.r = int(v * 255);
			color.g = int(t * 255);
			color.b = int(p * 255);
			break;
        case 1:
			// 黄到绿
			color.r = int(q * 255);
			color.g = int(v * 255);
			color.b = int(p * 255);
			break;
        case 2: 
			// 绿到青
			color.r = int(p * 255);
			color.g = int(v * 255);
			color.b = int(t * 255);
			break;
        case 3: 
			// 青到蓝
			color.r = int(p * 255);
			color.g = int(q * 255);
			color.b = int(v * 255);
			break;
        case 4: 
			// 蓝到紫
			color.r = int(t * 255);
			color.g = int(p * 255);
			color.b = int(v * 255);
			break;
        case 5: 
			// 紫到红
			color.r = int(v * 255);
			color.g = int(p * 255);
			color.b = int(q * 255);
			break;
    }
}

void ncs::alg::getPartColor(int partID, int totalPart,
	ncs::Color& color)
{
	// 计算颜色的色相
	float hue = static_cast<float>(partID) / totalPart;
	// 设置饱和度和亮度，这里可以根据需要调整
	float saturation = 1.0f;
	float value = 1.0f;

	// 将HSV转换为RGB
	hsvToRgb(hue * 360, saturation, value, color);
}

void ncs::alg::getResultColor(double vMin, double vMax, int value, ncs::Color& color)
{
	vMax += 1e-7;
	double rate = (value - vMin) / (vMax - vMin);
	if (rate > 1.0)
		rate = 1.0;
	double rc, gc, bc;
	if (rate <= 0.25)
	{
		rc = 0.0;
		gc = rate / 0.25;
		bc = 1.0;
	}
	else if (rate > 0.25 && rate <= 0.5)
	{
		rc = 0.0;
		gc = 1.0;
		bc = 1.0 - (rate - 0.25) / 0.25;
	}
	else if (rate > 0.5 && rate <= 0.75)
	{
		rc = (rate - 0.5) / 0.25;
		gc = 1.0;
		bc = 0.0;
	}
	else
	{
		rc = 1.0;
		gc = 1.0 - (rate - 0.75) / 0.25;
		bc = 0.0;
	}
	color.r = int(rc * 255);
	color.g = int(gc * 255);
	color.b = int(bc * 255);
}

void ncs::alg::calcBodyMassAndInertia(ncs::GeometryModel& model,
	std::unordered_map<int, ncs::Material>& mat,
	double& volume, double& mass, ncs::Vec3& Inertia, double& length,
	double& ratio, ncs::Vec3 *massCenPos)
{
	volume = 0;
	mass = 0;
	Inertia.SetValue(0, 0, 0);

	ncs::Vec3 massCenter(0, 0, 0);
	for (auto& et : model.solids)
	{
		et.density = mat[model.parts[et.partID].matID].density;
		calSolidEtVolume(et, model.nodes);
		et.mass = et.volume * et.density;
		volume += et.volume;
		mass += et.mass;

		et.centerPos.SetValue(0, 0, 0);
		for (int i = 0; i < 8; i++)
		{
			et.centerPos.x += model.nodes[et.nds[i]].x;
			et.centerPos.y += model.nodes[et.nds[i]].y;
			et.centerPos.z += model.nodes[et.nds[i]].z;
		}
		et.centerPos /= 8.0;

		massCenter += et.centerPos * et.mass;
	}
	massCenter /= mass;

	// translate the model to (0, 0, 0)
	double ymax = -1e30;
	double ymin = 1e30;
	for (auto& nd : model.nodes)
	{
		if (nd.y > ymax) ymax = nd.y;
		if (nd.y < ymin) ymin = nd.y;
		nd.x -= massCenter.x;
		nd.y -= massCenter.y;
		nd.z -= massCenter.z;
	}
	length = ymax - ymin;
	ratio = (massCenter.y - ymin) / length;

	if (massCenPos != nullptr)
	{
		massCenPos->x = massCenter.x;
		massCenPos->y = massCenter.y;
		massCenPos->z = massCenter.z;
	}

	// 计算转动惯量
	for (auto& et : model.solids)
	{
		et.centerPos -= massCenter;
		double a = pow(et.volume, 1.0 / 3.0);
		double selfInertia = et.mass * a * a / 6.0;

		ncs::Vec3 dis(et.centerPos.x, et.centerPos.y, et.centerPos.z);

		Inertia.x += et.mass * (dis.y * dis.y + dis.z * dis.z) + selfInertia;
		Inertia.y += et.mass * (dis.x * dis.x + dis.z * dis.z) + selfInertia;
		Inertia.z += et.mass * (dis.x * dis.x + dis.y * dis.y) + selfInertia;
	}
}

void ncs::alg::calcBodyInertia(const ncs::VNode& nodes, ncs::VSolidElement& solids, ncs::Vec3& Inertia)
{
	Inertia.SetValue(0, 0, 0);

	// 先计算质量中心
	ncs::Vec3 wdCenter;
	double mass = 0;
	for (auto& et : solids)
	{
		et.centerPos.SetValue(0, 0, 0);
		for (int i = 0; i < 8; i++)
		{
			et.centerPos.x += nodes[et.nds[i]].x;
			et.centerPos.y += nodes[et.nds[i]].y;
			et.centerPos.z += nodes[et.nds[i]].z;
		}
		et.centerPos /= 8.0;
		mass += et.mass;
		wdCenter += et.centerPos * et.mass;
	}
	wdCenter /= mass;

	// 计算转动惯量
	for (auto& et : solids)
	{
		double a = pow(et.volume, 1.0 / 3.0);
		double selfInertia = et.mass * a * a / 6.0;

		ncs::Vec3 dis = et.centerPos - wdCenter;

		Inertia.x += et.mass * (dis.y * dis.y + dis.z * dis.z) + selfInertia;
		Inertia.y += et.mass * (dis.x * dis.x + dis.z * dis.z) + selfInertia;
		Inertia.z += et.mass * (dis.x * dis.x + dis.y * dis.y) + selfInertia;
	}
}

void ncs::alg::calcInertia(ncs::VNode& nodes, ncs::VSolidElement& solids, ncs::Vec3& Inertia)
{
	for (auto& et : solids)
	{
		et.centerPos.SetValue(0, 0, 0);
		for (int i = 0; i < 8; i++)
		{
			et.centerPos.x += nodes[et.nds[i]].x;
			et.centerPos.y += nodes[et.nds[i]].y;
			et.centerPos.z += nodes[et.nds[i]].z;
		}
		et.centerPos /= 8.0;

		double a = pow(et.volume, 1.0 / 3.0);
		double selfInertia = et.mass * a * a / 6.0;

		ncs::Vec3 dis(et.centerPos.x, et.centerPos.y, et.centerPos.z);

		Inertia.x += et.mass * (dis.y * dis.y + dis.z * dis.z) + selfInertia;
		Inertia.y += et.mass * (dis.x * dis.x + dis.z * dis.z) + selfInertia;
		Inertia.z += et.mass * (dis.x * dis.x + dis.y * dis.y) + selfInertia;
	}
}

void ncs::alg::vectorToMatrix(double v[3], double m[3][3])
{
	m[0][0] = 0.0;
	m[0][1] = -v[2];
	m[0][2] = v[1];
	m[1][0] = v[2];
	m[1][1] = 0.0;
	m[1][2] = -v[0];
	m[2][0] = -v[1];
	m[2][1] = v[0];
	m[2][2] = 0.0;
}

bool ncs::alg::lineInterPlane(ncs::Vec3& origin, ncs::Vec3& dir, ncs::Vec3 A, ncs::Vec3 B, 
	ncs::Vec3 C, double& dis)
{
	ncs::Vec3 D, DE;
	D = origin;
	DE = dir;
	ncs::Vec3 DA = A - D;
	ncs::Vec3 AB = B - A;
	ncs::Vec3 AC = C - A;

	ncs::Vec3 normalABC = AB.cross(AC);
	if (DA.dot(DE) <= 0) return false;

	double d;
	d = DA.dot(normalABC) / (DE.dot(normalABC));


	ncs::Vec3 p0;
	p0.x = D.x + d * DE.x;
	p0.y = D.y + d * DE.y;
	p0.z = D.z + d * DE.z;

	double a, a1, a2, a3, eps;
	eps = 1e-6;

	a = AB.AreaWithOtherVector(AC);

	ncs::Vec3 pa, pb;
	pa = A - p0;
	pb = B - p0;
	a1 = pa.AreaWithOtherVector(pb);

	pa = B - p0;
	pb = C - p0;
	a2 = pa.AreaWithOtherVector(pb);

	pa = C - p0;
	pb = A - p0;
	a3 = pa.AreaWithOtherVector(pb);

	if (abs(a1 + a2 + a3 - a) > eps) return false;
	dis = (p0 - D).magnitude();
	return true;
}
