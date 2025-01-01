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

        Vec3 gravity; // ����
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
        double CSLH = 1.2;  // Ӧ�����ʵ�⻬���ȵĳ���

        // ����ʱ�����
        double dt = 5e-4;
        double tFactor = 0.2;
        int nstep = 500;
        int outputInterval = 3;

        // ճ����
        bool viscous = true;
        double viscousFactor = 0.01;   // ճ��ϵ��

        int dnesCalType = 0;          // 0: �ܶ���ͷ������ܶȣ� 1: ���ַ������ܶ�
        double restiCoeff = 0.9;      // ����ϵ��
        double fricCoeff = 0.2;       // Ħ��ϵ��
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


