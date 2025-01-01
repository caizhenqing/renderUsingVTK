#pragma once
#include"algFoundation.h"
namespace ncs 
{
	struct ALGBASE_API ncMatrix3
	{
		ncMatrix3();
		ncMatrix3(Vec3& v);
		ncMatrix3(double a[3][3]);

		ncMatrix3 vectorToMatrix(Vec3& v) const;
		ncMatrix3 inverseMatrix() const;
		ncMatrix3 transposition() const;
		Vec3 mutltVec(const Vec3& v) const;
		void setValue(double matrix[3][3]);
		double det() const;
		ncMatrix3 operator * (ncMatrix3& other) const;
		Vec3 operator * (Vec3& v) const;
		ncMatrix3 operator + (ncMatrix3& other) const;
		ncMatrix3 operator - (ncMatrix3& other) const;
		double getQuaternion(double& w, double& x, double& y, double& z);
		void operator = (const ncMatrix3& other);

		double m[3][3];
	};
}
