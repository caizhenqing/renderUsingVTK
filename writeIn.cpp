#include "writeIn.h"
#include <sstream>

bool kFileToGeoModel(ncs::GeometryModel& model, std::string fileName, 
	std::unordered_map<int, ncs::RcMaterial> &mats)
{
	std::ifstream files;
	files.open(fileName, std::ios::in);
	std::string tem;
	if (!files.is_open())
	{
		printf("No k file!");
		return false;
	}

	int sign = 0;
	ncs::Curve* curve = nullptr;

	while (getline(files, tem))
	{
		// tem = std::string::strim(tem);
		if (tem == "*NODE" || tem == "*node" || tem == "*Node")
		{
			sign = 1;
			continue;
		}
		else if (tem == "*ELEMENT_SOLID")
		{
			sign = 2;
			continue;
		}
		else if (tem == "*BIM_ELEMENT")
		{
			sign = 3;
			continue;
		}
		else if (tem == "*DAMAGE")
		{
			sign = 4;
			model.curves.emplace_back(ncs::Curve());
			curve = &model.curves.back();
			continue;
		}
		else if (tem == "*BIM_PART")
		{
			sign = 5;
			continue;
		}
		else if (tem == "*BIM_MATERIAL")
		{
			sign = 6;
			continue;
		}
		else if (tem == "*SET_SEGMENT")
		{
			sign = 7;
			continue;
		}
		else if (tem == "*PART" || tem == "*part")
		{
			sign = 8;
			continue;
		}
		else if (tem == "*ELEMENT_SHELL")
		{
			sign = 9;
			continue;
		}

		if (tem[0] == '$') continue;
		if (tem == "" || tem == " ") continue;
		if (tem[0] == '*')
		{
			sign = 0;
			continue;
		}

		int i = 0;
		int etID = 0;

		ncs::Node nd;
		ncs::SolidElement solid;
		ncs::DemElement demEt;
		ncs::Part part;
		ncs::RcMaterial mat;
		ncs::Quad quad;

		std::vector<double> content;
		if (tem.find(',') != std::string::npos)
		{
			std::stringstream ss(tem);
			std::string s;
			while (getline(ss, s, ','))
			{
				std::stringstream converter(s);
				double val;
				if (converter >> val) {
					// 如果转换成功，将值添加到向量中
					content.emplace_back(val);
				}
				else
				{
					content.clear();
					break;
				}
				// content.emplace_back(stod(s));
			}
		}
		else
		{
			std::istringstream ss(tem);
			std::string s;
			while (ss >> s)
			{
				std::stringstream converter(s);
				double val;
				if (converter >> val) {
					// 如果转换成功，将值添加到向量中
					content.emplace_back(val);
				}
				else
				{
					content.clear();
					break;
				}
				// content.emplace_back(stod(s));
			}
		}
		if (content.size() == 0) continue;

		switch (sign)
		{
		case 1:
			nd.id = int(content[0]);
			nd.x = content[1];
			nd.y = content[2];
			nd.z = content[3];
			model.nodes.emplace_back(nd);
			break;
		case 2:
			solid.etID = int(content[0]);
			solid.partID = int(content[1]);
			for (int ndID = 0; ndID < 8; ndID++)
			{
				solid.nds[ndID] = int(content[ndID + 2]) - 1;
			}
			model.solids.emplace_back(solid);
			break;
		case 3:
			demEt.etID = int(content[0]);
			demEt.partID = int(content[1]);
			demEt.dimX = content[2];
			demEt.dimY = content[3];
			demEt.dimZ = content[4];
			demEt.posX = content[5];
			demEt.posY = content[6];
			demEt.posZ = content[7];
			demEt.rotX = content[8];
			demEt.rotY = content[9];
			demEt.rotZ = content[10];
			model.demSolids.emplace_back(demEt);
			break;
		case 4:
			curve->paraX.emplace_back(content[0]);
			curve->paraY.emplace_back(std::vector<double>());
			for (int i = 0; i < content.size() - 1; i++)
			{
				curve->paraY.back().emplace_back(content[i + 1]);
			}
			break;
		case 5:
			part.type = int(content[1]);
			part.matID = int(content[2]);
			model.parts[int(content[0])] = part;
			break;
		case 6:
			if (content.size() != 14) break;
			// mat.fc = content[1];
			/*mat.preDis = content[2];
			mat.fy = content[3];
			mat.fs = content[4];
			mat.Es = content[5];
			mat.rcRatiot = content[6];
			mat.rcRatioc = content[7];
			mat.stRadius = content[8];
			mat.s = content[9];
			mat.dstrain = content[10];
			mat.density = content[11];
			mat.Ec = content[12];
			mat.n = content[13];*/
			// mats[int(content[0])] = mat;
			break;
		case 7:
			if (content.size() < 4) break;
			quad.partID = 1;
			quad.nds[0] = int(content[0]) - 1;
			quad.nds[1] = int(content[1]) - 1;
			quad.nds[2] = int(content[2]) - 1;
			quad.nds[3] = int(content[3]) - 1;
			model.quads.emplace_back(quad);
			break;
		case 8:
		{
			part.type = int(content[1]);
			part.matID = int(content[2]);
			//model.parts.emplace(int(content[0]), part);
			int x1 = content[0];
			model.parts[x1] = part;
			break;
		}
		case 9:
		{
			quad.partID = int(content[1]);
			for (int i = 0; i < 4; i++)
			{
				quad.nds[i] = int(content[i + 2] - 1);
			}
			model.quads.emplace_back(quad);
		}

		default:
			break;
		}
	}

	files.close();
	return true;
}
