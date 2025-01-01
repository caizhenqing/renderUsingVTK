#include "writeOut.h"
#include <fstream>
#include <stdio.h>

static FILE* fileOut = nullptr;

void geoModelToKfile(const ncs::GeometryModel& model, std::string fileName)
{
	std::ofstream file;
	if (fileName.find('.') == std::string::npos)
	{
		fileName += ".k";
	}
	file.open(fileName, std::ios::out);

	file << "*KEYWORD" << '\n';

	// node
	file << "*NODE" << '\n';
	int ndID = 1;
	for (auto& nd: model.nodes)
	{
		file << ndID << "," << nd.x << "," << nd.y << "," << nd.z << '\n';
		ndID++;
	}

	if (model.vertex.size() > 0)
	{
		file << "*VERTEX" << '\n';
	}
	int vID = 1;
	for (auto& v : model.vertex)
	{
		file << vID << "," << v.id +100 << "," << v.x << "," << v.y << "," << v.z << '\n';
		ndID++;
	}

	// tri
	int etID = 1;
	if (model.triangles.size() > 0) file << "*ELEMENT_SHELL" << '\n';
	for (auto &tri: model.triangles)
	{
		file << etID << "," << tri.partID << "," << tri.nds[0] + 1 << ",";
		file << tri.nds[1] + 1 << "," << tri.nds[2] + 1 << ",";
		file << tri.nds[2] + 1 << '\n';
		etID++;
	}

	// quad
	for (auto& quad : model.quads)
	{
		file << etID << "," << quad.partID << "," << quad.nds[0] + 1 << ",";
		file << quad.nds[1] + 1 << "," << quad.nds[2] + 1 << ",";
		file << quad.nds[3] + 1 << std::endl;
		etID++;
	}

	// solid
	if (model.solids.size() > 0) file << "*ELEMENT_SOLID" << '\n';
	for (auto &solid: model.solids)
	{
		file << etID << "," << solid.partID;
		for (int i = 0; i < 8; i++)
		{
			int ndID = 4 + i;
			if (i > 3) ndID = i - 4;

			file << "," << solid.nds[ndID] + 1;
		}
		file << std::endl;
		etID++;
	}
	file << "*END";
}

void curveToTxt(const ncs::Curve& curve, std::string fName)
{
	if (curve.paraX.size() == 0) return;
	if (curve.paraX.size() != curve.paraY.size()) return;

	std::ofstream file;
	std::string fileName = fName;
	if (fileName.find('.') == std::string::npos)
	{
		fileName += ".txt";
	}
	file.open(fileName, std::ios::out);
	if (curve.curveTitle.size()>0)
	{
		file << curve.curveTitle<<'\n';
	}
	if (curve.paraXTitle.size() > 0)
	{
		file << curve.paraXTitle;
	}
	for (unsigned int i = 0; i < curve.paraYTitle.size(); i++)
	{
		file << ", " << curve.paraYTitle[i];
	}
	file << '\n';

	for (int i = 0; i < curve.paraX.size(); i++)
	{
		int num = 12;
		file << std::setw(num) << std::right << std::fixed << std::setprecision(5);
		file << curve.paraX[i];
		
		for (int j = 0; j < curve.paraY[0].size(); j++)
		{
			file << ',';
			file << std::setw(num) << std::right << std::fixed << std::setprecision(5);
			file << curve.paraY[i][j];
		}
		file << '\n';
	}
	file.close();
}

void matrixWriteOut(double x, double matrix[3][3], const char* fileName, bool ending)
{
	// 如果文件指针是空的，打开文件以写入模式
	if (fileOut == nullptr)
	{
		fileOut = fopen(fileName, "w");
		if (fileOut == nullptr) {
			return;
		}
		// 写入表头
		fprintf(fileOut, "%12s%12s%12s%12s%12s%12s%12s%12s%12s%12s\n",
			"t", "m[0][0]", "m[0][1]", "m[0][2]",
			"m[1][0]", "m[1][1]", "m[1][2]",
			"m[2][0]", "m[2][1]", "m[2][2]");
	}

	fprintf(fileOut, "%12f", x);
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			fprintf(fileOut, "%12f", matrix[i][j]);
		}
	}
	fprintf(fileOut, "\n");

	if (ending)
	{
		fclose(fileOut);
		fileOut = nullptr;
	}
	
	
}
