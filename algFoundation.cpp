#include <vector>
#include "algFoundation.h"
//#include <iostream>

ncs::Node::~Node()
{
}

ncs::Node::Node()
{
	this->id = 0;
	this->x = 0.0;
	this->y = 0.0;
	this->z = 0.0;
}

ncs::Node::Node(double x_, double y_, double z_)
{
	this->id = 0;
	this->x = x_;
	this->y = y_;
	this->z = z_;
}

inline ncs::Vec3::Vec3()
{
	this->x = 0.0;
	this->y = 0.0;
	this->z = 0.0;
}

inline ncs::Vec3::Vec3(double x_, double y_, double z_)
{
	this->x = x_;
	this->y = y_;
	this->z = z_;
}

inline ncs::Vec3::Vec3(const Vec3& other)
{
	this->x = other.x;
	this->y = other.y;
	this->z = other.z;
}

inline ncs::Vec3::Vec3(double a)
{
	this->x = a;
	this->y = a;
	this->z = a;
}

inline ncs::Vec3::Vec3(ZERO)
{
	x = 0.0;
	y = 0.0;
	z = 0.0;
}

ncs::Vec3 ncs::Vec3::operator-(const Vec3& other) const
{
	return Vec3(this->x - other.x, this->y - other.y, this->z - other.z);
}

ncs::Vec3 ncs::Vec3::operator+(const Vec3& other) const
{
	return Vec3(this->x + other.x, this->y + other.y, this->z + other.z);
}

ncs::Vec3& ncs::Vec3::operator+=(const Vec3& other)
{
	this->x += other.x;
	this->y += other.y;
	this->z += other.z;
	return *this;

}

ncs::Vec3& ncs::Vec3::operator-=(const Vec3& other)
{
	this->x -= other.x;
	this->y -= other.y;
	this->z -= other.z;
	return *this;
}

ncs::Vec3& ncs::Vec3::operator*=(double d)
{
	x *= d;
	y *= d;
	z *= d;
	return *this;
	// TODO: 在此处插入 return 语句
}

ncs::Vec3& ncs::Vec3::operator/=(double d)
{
	d = 1.0 / d;
	x *= d;
	y *= d;
	z *= d;
	return *this;
	// TODO: 在此处插入 return 语句
}

ncs::Vec3 ncs::Vec3::operator*(double d) const
{
	return Vec3(this->x * d, this->y * d, this->z * d);
}

ncs::Vec3 ncs::Vec3::operator/(double d) const
{
	d = 1.0 / d;
	return Vec3(this->x * d, this->y * d, this->z * d);
}

inline ncs::Vec3 ncs::Vec3::operator=(const Vec3& other)
{
	this->x = other.x;
	this->y = other.y;
	this->z = other.z;
	return *this;
}

double& ncs::Vec3::operator[](unsigned int index)
{
	return reinterpret_cast<double*>(this)[index];
	// TODO: 在此处插入 return 语句
}

const double& ncs::Vec3::operator[](unsigned int index) const
{
	return reinterpret_cast<const double*>(this)[index];
	// TODO: 在此处插入 return 语句
}

bool ncs::Vec3::operator==(const Vec3& other) const
{
	return x == other.x && y == other.y && z == other.z;
}

bool ncs::Vec3::operator!=(const Vec3& other) const
{
	return x != other.x || y != other.y || z != other.z;
}

bool ncs::Vec3::isZero() const
{
	return x == 0.0 && y == 0.0 && z == 0.0;
}

double ncs::Vec3::magnitudeSquared() const
{
	return x * x + y * y + z * z;
}

double ncs::Vec3::magnitude() const
{
	double moduls = magnitudeSquared();
	if (moduls > 0.0)
	{
		return sqrt(magnitudeSquared());
	}
	return 0;
}

ncs::Vec3 ncs::Vec3::operator-() const
{
	return Vec3(-x, -y, -z);
}

double ncs::Vec3::dot(const Vec3& other) const
{
	return x * other.x + y * other.y + z * other.z;
}

ncs::Vec3 ncs::Vec3::multiply(const Vec3& other) const
{
	return ncs::Vec3(x * other.x, y * other.y, z * other.z);
}

ncs::Vec3 ncs::Vec3::cross(const Vec3& other) const
{
	return Vec3(y * other.z - z * other.y, z * other.x - x * other.z, x * other.y - y * other.x);
}


double ncs::Vec3::AreaWithOtherVector(const Vec3& other) const
{
	return this->cross(other).magnitude() * 0.5;
}

void ncs::Vec3::normalize()
{
	const double m = magnitude();
	if (m > double(0.0))
		*this /= m;
}

ncs::Vec3 ncs::Vec3::getNormalized() const
{
	const double m = magnitude();
	return m > double(0.0) ? *this / m : ncs::Vec3(double(0));
}

void ncs::Vec3::SetValue(double fx, double fy, double fz)
{
	x = fx;
	y = fy;
	z = fz;
}

double ncs::Vec3::minElement() const
{
	return min(x, min(y, z));
}

double ncs::Vec3::maxElement() const
{
	return max(x, max(y, z));
}

ncs::Vec3 ncs::Vec3::absolute() const
{
	return ncs::Vec3(abs(x), abs(y), abs(z));
}


ncs::RotMatrix::RotMatrix()
{
	posx = 0;
	posy = 0;
	posz = 0;
	angx = 0;
	angy = 0;
	angz = 0;
}

void ncs::RotMatrix::setAttitude(const Vec3& pos, const Vec3& ang)
{
	posx = pos.x;
	posy = pos.y;
	posz = pos.z;
	angx = ang.x;
	angy = ang.y;
	angz = ang.z;

	this->matrixGenerate();
}

void ncs::RotMatrix::setAttitude(double px, double py, double pz, double rotx, double roty, double rotz)
{
	posx = px;
	posy = py;
	posz = pz;
	angx = rotx;
	angy = roty;
	angz = rotz;

	this->matrixGenerate();

}

ncs::Vec3 ncs::RotMatrix::getByRotate(const Vec3& in)
{
	Vec3 tem;
	tem.x = in.x * rotMatrix[0][0] + in.y * rotMatrix[0][1] + in.z * rotMatrix[0][2];
	tem.y = in.x * rotMatrix[1][0] + in.y * rotMatrix[1][1] + in.z * rotMatrix[1][2];
	tem.z = in.x * rotMatrix[2][0] + in.y * rotMatrix[2][1] + in.z * rotMatrix[2][2];
	return tem;
}

ncs::Vec3 ncs::RotMatrix::getByRotateAndOffset(const Vec3& in)
{
	Vec3 tem;
	tem.x = in.x * rotMatrix[0][0] + in.y * rotMatrix[0][1] + in.z * rotMatrix[0][2] + posx;
	tem.y = in.x * rotMatrix[1][0] + in.y * rotMatrix[1][1] + in.z * rotMatrix[1][2] + posy;
	tem.z = in.x * rotMatrix[2][0] + in.y * rotMatrix[2][1] + in.z * rotMatrix[2][2] + posz;
	return tem;
}

ncs::Node ncs::RotMatrix::getByRotate(const Node& in)
{
	Node tem;
	tem.x = in.x * rotMatrix[0][0] + in.y * rotMatrix[0][1] + in.z * rotMatrix[0][2];
	tem.y = in.x * rotMatrix[1][0] + in.y * rotMatrix[1][1] + in.z * rotMatrix[1][2];
	tem.z = in.x * rotMatrix[2][0] + in.y * rotMatrix[2][1] + in.z * rotMatrix[2][2];
	return tem;
}

ncs::Node ncs::RotMatrix::getByRotateAndOffset(const Node& in)
{
	Node tem;
	tem.x = in.x * rotMatrix[0][0] + in.y * rotMatrix[0][1] + in.z * rotMatrix[0][2] + posx;
	tem.y = in.x * rotMatrix[1][0] + in.y * rotMatrix[1][1] + in.z * rotMatrix[1][2] + posy;
	tem.z = in.x * rotMatrix[2][0] + in.y * rotMatrix[2][1] + in.z * rotMatrix[2][2] + posz;
	return tem;
}

double ncs::RotMatrix::A2R(double angle)
{
	return (angle / 180.0) * Pi;
}

void ncs::RotMatrix::matrixGenerate()
{
	double rotX[3][3];//roll
	double rotY[3][3];//pitch
	double rotZ[3][3];//yaw

	double alpha = A2R(angx);
	double beta = A2R(angy);
	double gama = A2R(angz);

	rotX[0][0] = 1.0;
	rotX[1][0] = 0.0;
	rotX[2][0] = 0.0; //
	rotX[0][1] = 0.0;
	rotX[1][1] = cos(alpha);
	rotX[2][1] = sin(alpha); //
	rotX[0][2] = 0.0;
	rotX[1][2] = -sin(alpha);
	rotX[2][2] = cos(alpha);

	rotY[0][0] = cos(beta);
	rotY[1][0] = 0.0;
	rotY[2][0] = -sin(beta);
	rotY[0][1] = 0.0;
	rotY[1][1] = 1.0;
	rotY[2][1] = 0.0;
	rotY[0][2] = sin(beta);
	rotY[1][2] = 0.0;
	rotY[2][2] = cos(beta);

	rotZ[0][0] = cos(gama);
	rotZ[1][0] = sin(gama);
	rotZ[2][0] = 0.0;
	rotZ[0][1] = -sin(gama);
	rotZ[1][1] = cos(gama);
	rotZ[2][1] = 0.0;
	rotZ[0][2] = 0.0;
	rotZ[1][2] = 0.0;
	rotZ[2][2] = 1.0;

	double tmp[3][3];

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			tmp[i][j] = rotY[i][0] * rotX[0][j] + rotY[i][1] * rotX[1][j] + rotY[i][2] * rotX[2][j];
		}
	}

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			rotMatrix[i][j] = rotZ[i][0] * tmp[0][j] + rotZ[i][1] * tmp[1][j] + rotZ[i][2] * tmp[2][j];
		}
	}
}

ncs::Color::Color()
{
	this->r = 0;
	this->g = 0;
	this->b = 0;
}

ncs::Color::Color(double r_, double g_, double b_)
{
	r = r_;
	g = g_;
	b = b_;
};

void ncs::Color::setColor(double vMin, double vMax, double value)
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
	this->r = int(rc * 255);
	this->g = int(gc * 255);
	this->b = int(bc * 255);
}

ncs::Line::Line()
{
	this->nd1 = 0;
	this->nd2 = 0;
}

ncs::Line::Line(int n1, int n2)
{
	this->nd1 = n1;
	this->nd2 = n2;
}

ncs::Part::Part()
{
}

ncs::Part::Part(int tp, int mID)
{
	this->type = tp;
	this->matID = mID;
}

//void ncs::Material::calSectionStrength(Threshold& th, double len, double h, double b)
//{
//	double rRatio = rcRatioc + rcRatiot;
//	double fai = this->getFactor(len / b);
//
//	double area = h * b;
//	// 最大压力
//	th.fCompressive = 0.9 * fai * (rRatio * fy + (1 - rRatio) * fc) * area;
//	// 最大拉力
//	th.fTensile = rRatio * fs * area;
//	// 最大剪力, 抗拉强度取抗压强度 1/10
//	double h0 = h - this->preDis;
//	th.fShear = 0.7 * this->fc * 0.1 * h0 * b;
//	th.fShear += fy * this->n * Pi * stRadius * stRadius / s * h0;
//	// 拉伸刚度
//	th.linearStiffness = this->Ec * h * b / len;
//	// 剪切刚度
//	th.shearStiffness = th.linearStiffness;
//	// 最大弯矩，沿截面 b 方向
//	double fb, x;
//	fb = 0.8 / (1 + this->fy / (0.003 * this->Es));
//	x = fb * h0 / 2;
//	th.moment1 = b * x * (h0 - x / 2.0) * this->fc +
//		rcRatioc * (h0 - this->preDis) * b * h * fy;
//	// 最大弯矩，沿截面 h 方向
//	h0 = b - this->preDis;
//	x = fb * h0 / 2;
//	th.moment2 = h * x * (h0 - x / 2.0) * this->fc +
//		rcRatioc * (h0 - this->preDis) * h * b * fy;
//	// 弯曲刚度
//
//	// 扭刚度
//
//	// th.strain = 0.3;
//}

//void ncs::ncMaterial::calSectionStrength(ncsThread& th)
//{
//	double rRatio = rcRatioc + rcRatiot;
//	double hoPerh = 0.18;
//
//	th.tensile = rRatio * fs;
//	th.compressive = 0.9 * (rRatio * fy + (1 - rRatio) * fc);
//	// 抗拉强度取抗压强度1/10
//	th.shear = 0.7 * this->fc * 0.1 * hoPerh;
//	th.shear += fy * this->n * ncPI * stRadius * stRadius / s;
//
//	double fb, x;
//	fb = 0.8 / (1 + this->fy / (0.003 * this->Es));
//	x = fb * h0 / 2;
//	th.bend = x * (h0 - x / 2.0) * this->fc / h +
//		rcRatioc * (h0 - this->preDis) * fy;
//
//	th.strain = 0.3;
//}

//double ncs::Material::getFactor(double lb)
//{
//	std::vector<double> x{ 8,10,12,14,16,18,20,22,24,26,28,
//		30,32,34,36,38,40,42,44,46,48,50 };
//
//	std::vector<double> y{ 1,0.98,0.95,0.92,0.87,0.81,0.75,0.70,0.65,
//		0.60,0.56,0.52,0.48,0.44,0.40,0.36,0.32,0.29,0.26,0.23,0.21,0.19 };
//
//	for (int i = 0; i < x.size(); i++)
//	{
//		if (lb <= x[i]) return y[i];
//	}
//	return 0.19;
//}

double ncs::min(double a, double b)
{
	return a < b ? a : b;
}

double ncs::max(double a, double b)
{
	return a < b ? b : a;
}

ncs::Extent::Extent()
{
	xmin = 0;
	xmax = 0;
	ymin = 0;
	ymax = 0;
	zmin = 0;
	zmax = 0;
}
