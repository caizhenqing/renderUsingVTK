#pragma once
#include "renderBase.h"

#ifdef RENDERBASE_EXPORTS
#define RENDERBASE_API _declspec(dllexport)
#else
#define RENDERBASE_API _declspec(dllimport)
#endif

RENDERBASE_API RenderBase* CreateRenderUsingVTK();


