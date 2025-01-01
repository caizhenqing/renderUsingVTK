#pragma once

#include "algFoundation.h"
using namespace ncs;

namespace ncs
{
    struct particlePair
    {
        unsigned int p1;
        unsigned int p2;

        double w;
        Vec3 dw;

        double w1;
        Vec3 dw1;
        Vec3 secW;
    };
    typedef std::vector <particlePair> VparticlePair;

    struct FluidBlocks
    {
        int id;
        Extent box;
        double dens;
        double mass;
        double pressure;
        double particleRadius;
        double smoothLen;

        Vec3 vel;
    };

    struct sphSetting
    {
        Extent box;

        Vec3 gravity; // 重力
        int kernelID = 3;
        int dimension;
        double particleRadius;
        double supportR;
        double dens0;
        double stiffness;
        double exponent;

        double time;
        double stepSize;
        int aglID;
        int densCalType;
        std::vector<FluidBlocks> fBlocks;

        int maxNeighbor = 80;

        double Hmin = 0.2;
        double Hmax = 2.0;
        double CSLH = 1.2;  // 应用于质点光滑长度的常量

        // 计算时间控制
        double dt = 5e-4;
        double tFactor = 0.2;
        int nstep = 500;
        int outputInterval = 3;

        // 粘性项
        bool viscous = true;
        double viscousFactor = 0.01;   // 粘性系数

        int dnesCalType = 0;          // 0: 密度求和法计算密度； 1: 积分法计算密度
        double restiCoeff = 0.9;      // 反弹系数
        double fricCoeff = 0.2;       // 摩擦系数
        double surfaceTension = 0.01;
        double threshold = 7.065;
        double artSoundVel;
        double h0;
    };

    struct sphInput
    {
        sphSetting config;
        GeometryModel model;
        std::unordered_map<int, Material> mat;

    };
}


