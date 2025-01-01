#pragma once
#include "algFoundation.h"
// #include <iostream>
#include <iomanip>

ALGBASE_API void geoModelToKfile(const ncs::GeometryModel& model, std::string fileName);

ALGBASE_API void curveToTxt(const ncs::Curve& curve, std::string fileName);

// ending = true,Çå¿ÕÄÚ´æ
ALGBASE_API void matrixWriteOut(double x, double matrix[3][3], const char* fileName, 
	bool ending = false);