#include "matrix.h"


ncs::ncMatrix3::ncMatrix3()
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			this->m[i][j] = 0;
		}
	}
	m[0][0] = 1;
	m[1][1] = 1;
	m[2][2] = 1;
}

ncs::ncMatrix3::ncMatrix3(Vec3& v)
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			this->m[i][j] = 0;
		}
	}
	this->m[0][0] = v.x;
	this->m[1][1] = v.y;
	this->m[2][2] = v.z;
}

ncs::ncMatrix3::ncMatrix3(double a[3][3])
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			this->m[i][j] = a[i][j];
		}
	}
}

ncs::ncMatrix3 ncs::ncMatrix3::vectorToMatrix(Vec3& v) const
{
	ncs::ncMatrix3 matrix;
	matrix.m[0][1] = -v.z;
	matrix.m[0][2] = v.y;
	matrix.m[1][0] = v.z;
	matrix.m[1][2] = -v.x;
	matrix.m[2][0] = -v.y;
	matrix.m[2][1] = v.x;
	return matrix;
}

ncs::ncMatrix3 ncs::ncMatrix3::inverseMatrix() const
{
	ncMatrix3 a;
	a = *this;

	int n = 3;
	int* is, * js, i, j, k;
	double d, p;
	is = new int[n];
	js = new int[n];
	for (k = 0; k <= n - 1; k++)
	{
		d = 0.0;
		for (i = k; i <= n - 1; i++)
		{
			for (j = k; j <= n - 1; j++)
			{
				p = fabs(a.m[i][j]);
				if (p > d)
				{
					d = p;
					is[k] = i;
					js[k] = j;
				}
			}
		}
		if (d + 1.0 == 1.0)
		{
			delete[] is, js;
		}

		if (is[k] != k)
		{
			for (j = 0; j <= n - 1; j++)
			{
				p = a.m[k][j];
				a.m[k][j] = a.m[is[k]][j];
				a.m[is[k]][j] = p;
			}
		}
		if (js[k] != k)
		{
			for (i = 0; i <= n - 1; i++)
			{
				p = a.m[i][k];
				a.m[i][k] = a.m[i][js[k]];
				a.m[i][js[k]] = p;
			}
		}
		a.m[k][k] = 1.0 / a.m[k][k];
		for (j = 0; j <= n - 1; j++)
		{
			if (j != k)
				a.m[k][j] = a.m[k][j] * a.m[k][k];
		}
		for (i = 0; i <= n - 1; i++)
		{
			if (i != k)
			{
				for (j = 0; j <= n - 1; j++)
				{
					if (j != k)
						a.m[i][j] = a.m[i][j] - a.m[i][k] * a.m[k][j];
				}
			}
		}
		for (i = 0; i <= n - 1; i++)
		{
			if (i != k)
				a.m[i][k] = -a.m[i][k] * a.m[k][k];
		}
	}
	for (k = n - 1; k >= 0; k--)
	{
		if (js[k] != k)
		{
			for (j = 0; j <= n - 1; j++)
			{
				p = a.m[k][j];
				a.m[k][j] = a.m[js[k]][j];
				a.m[js[k]][j] = p;
			}
		}
		if (is[k] != k)
		{
			for (i = 0; i <= n - 1; i++)
			{
				p = a.m[i][k];
				a.m[i][k] = a.m[i][is[k]];
				a.m[i][is[k]] = p;
			}
		}
	}
	delete[] is, js;
	return a;
}

ncs::ncMatrix3 ncs::ncMatrix3::transposition() const
{
	ncs::ncMatrix3 mOut;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			mOut.m[j][i] = this->m[i][j];
		}
	}

	return mOut;
}

ncs::Vec3 ncs::ncMatrix3::mutltVec(const Vec3& v) const
{
	Vec3 out;
	for (int i = 0; i < 3; i++)
	{
		out[i] = this->m[i][0] * v[0] + this->m[i][1] * v[1] + this->m[i][2] * v[2];
	}
	return out;
}

void ncs::ncMatrix3::setValue(double matrix[3][3])
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			this->m[i][j] = matrix[i][j];
		}
	}

		
}

double ncs::ncMatrix3::det() const
{
	double res = 0;

	res += m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1]);
	res -= m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0]);
	res += m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0]);
	return res;
}

ncs::ncMatrix3 ncs::ncMatrix3::operator*(ncMatrix3& mB) const
{
	ncs::ncMatrix3 mOut;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			mOut.m[i][j] = this->m[i][0] * mB.m[0][j] + this->m[i][1] * mB.m[1][j] + this->m[i][2] * mB.m[2][j];
		}
	}
	return mOut;
}

ncs::Vec3 ncs::ncMatrix3::operator*(Vec3& v) const
{
	ncs::Vec3 vOut;
	vOut.x = this->m[0][0] * v.x + this->m[0][1] * v.y + this->m[0][2] * v.z;
	vOut.y = this->m[1][0] * v.x + this->m[1][1] * v.y + this->m[1][2] * v.z;
	vOut.z = this->m[2][0] * v.x + this->m[2][1] * v.y + this->m[2][2] * v.z;
	return vOut;
}

ncs::ncMatrix3 ncs::ncMatrix3::operator+(ncMatrix3& other) const
{
	ncMatrix3 matrix;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			matrix.m[i][j] = this->m[i][j] + other.m[i][j];
		}
	}
	return matrix;
}

ncs::ncMatrix3 ncs::ncMatrix3::operator-(ncMatrix3& other) const
{
	ncMatrix3 matrix;
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			matrix.m[i][j] = this->m[i][j] - other.m[i][j];
		}
	}
	return matrix;
}

double ncs::ncMatrix3::getQuaternion(double& w, double& x, double& y, double& z)
{


	return 0.0;
}

void ncs::ncMatrix3::operator=(const ncMatrix3& other)
{
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			this->m[i][j] = other.m[i][j];
		}
	}
}