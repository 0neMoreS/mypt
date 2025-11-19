#pragma once
#include "material.h"
#include "light.h"


struct IntersectInfo
{
    float t; // distance to the hit point
    SurfaceInfo surfaceInfo;
    std::shared_ptr<BxDF> bxdf;
    std::shared_ptr<Light> light;
    // const Primitive* hitPrimitive;

    IntersectInfo& operator=(const IntersectInfo& other)
    {
        t = other.t;
        surfaceInfo = other.surfaceInfo;
        bxdf = other.bxdf;
        light = other.light;
        return *this;
    }

    bool hasAreaLight() const { return !light->isNull(); }

    bool hasAreaLight(const std::shared_ptr<Light>& light) const
    {
        return !light->isNull() && light == light;
    }

    // return emission
    Vec3f Le(const Vec3f& dir) const
    {
        return light->Le(surfaceInfo, dir);
    }

    BxDFType getBxDFType() const { return bxdf->getType(); }

    Vec3f evaluateBxDF(const Vec3f& wo, const Vec3f& wi, const TransportDirection& mode) const
    {
        // world to local transform
        const Vec3f wo_l = worldToLocal(wo, surfaceInfo.dpdu, surfaceInfo.shadingNormal, surfaceInfo.dpdv);
        const Vec3f wi_l = worldToLocal(wi, surfaceInfo.dpdu, surfaceInfo.shadingNormal, surfaceInfo.dpdv);

        return bxdf->evaluate(wo_l, wi_l, mode);
    }

    // sample direction by BxDF
    // its pdf is propotional to the shape od BxDF
    Vec3f sampleBxDF(const Vec3f& wo, const TransportDirection& mode, Sampler& sampler, Vec3f& wi, float& pdf) const
    {
        // world to local transform
        const Vec3f wo_l = worldToLocal(wo, surfaceInfo.dpdu, surfaceInfo.shadingNormal, surfaceInfo.dpdv);

        // sample direction in tangent space
        Vec3f wi_l;
        const Vec3f f = bxdf->sampleDirection(wo_l, mode, sampler, wi_l, pdf);

        // local to world transform
        wi = localToWorld(wi_l, surfaceInfo.dpdu, surfaceInfo.shadingNormal, surfaceInfo.dpdv);

        return f;
    }

    // get all samplable direction
    std::vector<DirectionPair> sampleAllBxDF(const Vec3f& wo, const TransportDirection& mode) const
    {
        // world to local transform
        const Vec3f wo_l = worldToLocal(wo, surfaceInfo.dpdu, surfaceInfo.shadingNormal, surfaceInfo.dpdv);

        // sample all direction in tangent space
        std::vector<DirectionPair> dir_pairs = bxdf->sampleAllDirection(wo_l, mode);

        // local to world transform
        for (auto& dp : dir_pairs)
        {
            dp.first = localToWorld(dp.first, surfaceInfo.dpdu, surfaceInfo.shadingNormal, surfaceInfo.dpdv);
        }

        return dir_pairs;
    }

    SurfaceInfo sampleLightPoint(Sampler& sampler, float& pdf) const
    {
        if (light != nullptr)
        {
            return light->samplePoint(sampler, pdf);
        }
        return SurfaceInfo();
    }
};