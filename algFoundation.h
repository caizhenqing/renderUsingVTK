#pragma once
#include <vector>
#include <unordered_map>
#include <fstream>
#include "dllExport.h"
#include <math.h>


namespace ncs
{
	double min(double a, double b);
	double max(double a, double b);
	static const float Pi = float(3.141592653589793);
	static const float angToRad = float(0.017453292519943295);
	static const float HalfPi = float(1.57079632679489661923);
	static const float TwoPi = float(6.28318530717958647692);
	static const float InvPi = float(0.31830988618379067154);
	static const float InvTwoPi = float(0.15915494309189533577);
	static const float PiDivFour = float(0.78539816339744830962);
	static const float Sqrt2 = float(1.4142135623730951);
	static const float InvSqrt2 = float(0.7071067811865476);

	enum EMPTY
	{
		Empty
	};

	enum ZERO
	{
		Zero
	};

	enum ONE
	{
		One
	};

	enum IDENTITY
	{
		Identity
	};

	struct ALGBASE_API Node
	{
		 ~Node();
		int id;
		double x;
		double y;
		double z;
		Node();
		Node(double x_, double y_, double z_);
	};
	typedef std::vector<Node> VNode;

	struct ALGBASE_API Line
	{
		int nd1, nd2;
		int partID;
		Line();
		Line(int n1, int n2);
	};

	struct ALGBASE_API Vec3
	{
		inline Vec3();
		inline Vec3(double a);
		inline Vec3(ZERO);
		inline Vec3(double x_, double y_, double z_);
		inline Vec3(const Vec3& other);

		inline Vec3 operator = (const Vec3& other);
		double& operator[](unsigned int index);
		const double& operator[](unsigned int index) const;

		bool operator == (const Vec3& other) const;
		bool operator != (const Vec3& other) const;

		bool isZero() const;
		double magnitudeSquared() const;
		double magnitude() const;

		Vec3 operator -() const;
		Vec3 operator + (const Vec3& other) const;
		Vec3 operator - (const Vec3& other) const;

		Vec3 operator * (double d) const;
		Vec3 operator / (double d) const;
		Vec3& operator += (const Vec3& other);
		Vec3& operator -= (const Vec3& other);
		Vec3& operator *= (double d);
		Vec3& operator /= (double d);

		double dot(const Vec3& other) const;
		Vec3 multiply(const Vec3& other) const;
		Vec3 cross(const Vec3& other) const;
		double AreaWithOtherVector(const Vec3& other) const;
		void normalize();
		Vec3 getNormalized() const;
		void SetValue(double fx, double fy, double fz);
		double minElement() const;
		double maxElement() const;
		Vec3 absolute() const;
		double x;
		double y;
		double z;
	};
	typedef std::vector<Vec3> VVec3;

	struct ALGBASE_API Curve
	{
		std::string curveTitle;
		std::string paraXTitle;
		std::vector<double> paraX;
		std::vector<std::string> paraYTitle;
		std::vector<std::vector<double>> paraY;
	};

	struct ALGBASE_API Triangle
	{
		int partID;
		int triID;
		int nds[3];
		double area;
		double thickness;
		Vec3 normal;
		Vec3 centerPos;
		Triangle()
		{
			area = 0;
			thickness = 0;
		}
		Triangle(int n1, int n2, int n3)
		{
			area = 0;
			thickness = 0;
			nds[0] = n1;
			nds[1] = n2;
			nds[2] = n3;
		}
	};
	typedef std::vector<Triangle> VTriangle;

	struct ALGBASE_API Quad
	{
		Quad(){};
		Quad(int n1, int n2, int n3, int n4)
		{
			nds[0] = n1;
			nds[1] = n2;
			nds[2] = n3;
			nds[3] = n4;
		};
		int partID;
		int quadID;
		int nds[4]{ 0 };
		double area;
		Vec3 normal;
	};
	typedef std::vector<Quad> VQuad;

	struct Face
	{
		// 是否是自由面
		bool isFreeface = false;
		// 如果不是自由面，邻接面ID
		int adjacencyEtID = 0;
		// 此面的四个节点
		int nds[4];
	};
	typedef std::vector<Face> VFace;

	struct SolidElement
	{
		int partID;
		int etID = 0;
		int nds[8];
		Face faces[6];

		double mass = 0;
		double volume = 0;
		double density = 0;
		Vec3 centerPos;
	};
	typedef std::vector<SolidElement> VSolidElement;

	struct DemElement
	{
		int partID;
		int etID;
		double dimX;
		double dimY;
		double dimZ;
		double posX;
		double posY;
		double posZ;
		double rotX;
		double rotY;
		double rotZ;
	};

	struct Threshold
	{
		double fCompressive;    // 最大压力
		double fTensile;        // 最大拉力
		double fShear;          // 最大剪力
		// linear spring
		double linearStiffness;  // 拉伸刚度
		double shearStiffness;   // 剪切刚度

		double moment1;          // 最大弯矩，沿截面 b 方向
		double moment2;          // 最大弯矩，沿截面 h 方向
		double torque;           // 最大扭矩

		// angular spring
		double bendingStiffness;  // 弯曲刚度
		double torqueStiffness;   // 扭刚度
	};

	struct RcMaterial
	{
		double prxy;
		double destroyStrain;    // 破坏应变

		// 钢筋
		double rEx;             // 钢筋弹性模量
		double rDensity;        // 钢筋密度
		double rYeildstress;    // 钢筋屈服强度
		double logRadious;      // 纵筋截面半径
		int logNumber;          // 纵筋肢数，一般为2或者4
		double stripRadious;    // 箍筋截面半径
		double stripDis;        // 箍筋间距

		// 混凝土
		double cDensity;        // 混凝土密度
		double cEx;             // 混凝土弹性模量
		double cFc;             // 混凝土抗压强度
	};



	struct ALGBASE_API Material
	{
		int matID;
		double density = 0;                // 密度
		double yeidstress = 0;          // 屈服强度
		double shearMoudle = 0;         // 剪切模量
		double Ex = 0;                  // 弹性模量
		double fc = 0;                  // 抗压强度
		double prxy;

		// rebar
		//double Es;                  // 钢筋弹性模量
		//double fy;                  // 钢筋屈服强度
		//double fs;                  // 钢筋极限强度
		//double dstrain;             // 钢筋破坏应变
		//double rcRatiot;            // 纵筋受拉区配筋率
		//double rcRatioc;            // 纵筋受压区配筋率
		//double stRadius;            // 箍筋截面半径
		//double s;                   // 长度方向箍筋距离
		
		// concrete
		//double density;             // 混凝土材料密度
		//double Ec;                  // 混凝土弹性模量
		//double fc;                  // 混凝土抗压强度
		//double preDis;              // 混凝土保护尺寸

		// int n = 2;
		/*void calSectionStrength(Threshold& th, double len, double h, double b);
		double getFactor(double lb);*/
	};

	struct ALGBASE_API Part
	{
		int type;   // 构件类型
		int matID;
		Part();
		// type and matID
		Part(int tp, int mID);
	};

	struct ALGBASE_API Particle
	{
		int ndID;
		int partID;
		double mass;
		Vec3 v0;
		double dens;
		double pressure;
		double smoothLen;
		double diameter;
	};
	typedef std::vector<Particle> VParticle;

	struct ALGBASE_API Vertex
	{
		int id;
		int partID;
		double x;
		double y;
		double z;
		Vertex();
		Vertex(double x_, double y_, double z_);
	};
	typedef std::vector<Vertex> VVertex;

	struct ALGBASE_API GeometryModel
	{
		std::vector<Node> nodes;           // 每个节点的坐标信息
		std::vector<Vertex> vertex;
		std::vector<Particle> particles;
		std::vector<Line> lines;
		std::vector<Curve> curves;
		std::vector<Triangle> triangles;
		std::vector<Quad> quads;
		std::vector<SolidElement> solids;  // 离散元单元信息
		std::vector<DemElement> demSolids;
		std::unordered_map<int, ncs::Part> parts;
	};


	class ALGBASE_API RotMatrix
	{
	public:
		RotMatrix();
		void setAttitude(const Vec3& pos, const Vec3& ang);
		void setAttitude(double px, double py, double pz, double rotx, double roty, double rotz);

		Vec3 getByRotate(const Vec3& in);
		Vec3 getByRotateAndOffset(const Vec3& in);

		Node getByRotate(const Node& in);
		Node getByRotateAndOffset(const Node& in);
		double rotMatrix[3][3]{ {0} };

	private:
		double posx, posy, posz;
		double angx, angy, angz;

		double A2R(double angle);
		void matrixGenerate();
	};


	struct TimeSetting
	{
		double timeStep;
		double calcTime;
		int IntervalStepNum;
		TimeSetting()
		{
			timeStep = 0.0;
			calcTime = 0.0;
			IntervalStepNum = 0;
		}
	};


	struct Buffer
	{
		std::vector<std::vector<double>> nds;
		std::vector<std::vector<int>> tris;
		std::vector<std::vector<int>> colors;
		std::vector<std::vector<double>> textures;
		std::vector<std::vector<double>> normals;
	};
	typedef std::vector<Buffer> buffers;



	struct ALGBASE_API Color
	{
		Color();
		Color(double r_, double g_, double b_);
		void setColor(double vMin, double vMax, double value);

		int r;
		int g;
		int b;
	};

	struct ALGBASE_API  Extent
	{
		double xmin, xmax, ymin, ymax, zmin, zmax;
		Extent();
	};
}



