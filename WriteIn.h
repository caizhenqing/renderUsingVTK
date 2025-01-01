#pragma once
#include "algFoundation.h"
#include "dllExport.h"

ALGBASE_API bool kFileToGeoModel(ncs::GeometryModel& model, std::string fileName,
	std::unordered_map<int, ncs::RcMaterial> &mats);
