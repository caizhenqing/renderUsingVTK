#pragma once

#ifdef _WIN32
#ifdef ALGBASE_EXPORTS
#define ALGBASE_API _declspec(dllexport)
#else
#define ALGBASE_API _declspec(dllimport)
#endif
#else
#define ALGBASE_API
#endif