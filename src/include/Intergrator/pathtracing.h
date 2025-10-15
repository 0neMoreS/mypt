#pragma once
#include "integrator.h"
// implementation of path tracing
class PathTracing : public Integrator
{
private:
    const int maxDepth;

public:
    PathTracing(int maxDepth = 100) : maxDepth(maxDepth) {}

    void build(const Scene& scene, Sampler& sampler) override {}

    Vec3f integrate(const Ray& ray_in, const Scene& scene, Sampler& sampler) const override
    {
        Vec3f radiance(0);
        Ray ray = ray_in;
        Vec3f throughput(1, 1, 1);
        for (int k = 0; k < maxDepth; ++k)
        {
            IntersectInfo info;
            if (scene.intersect(ray, info))
            {
                if (k >= maxDepth / 2) {
                    // russian roulette
                    const float russian_roulette_prob = std::min(std::max(throughput[0], std::max(throughput[1], throughput[2])), 1.0f);
                    if (sampler.getNext1D() >= russian_roulette_prob)
                    {
                        break;
                    }
                    throughput /= russian_roulette_prob;
                }

                // Le
                if (info.hitPrimitive->hasAreaLight())
                {
                    radiance += throughput * info.hitPrimitive->Le(info.surfaceInfo, -ray.direction);
                    break;
                }

                // sample direction by BxDF
                Vec3f dir;
                float pdf_dir;
                Vec3f f = info.hitPrimitive->sampleBxDF(-ray.direction, info.surfaceInfo, TransportDirection::FROM_CAMERA, sampler, dir, pdf_dir);

                // update throughput and ray
                throughput *= f * cosTerm(-ray.direction, dir, info.surfaceInfo, TransportDirection::FROM_CAMERA) / pdf_dir;
                ray = Ray(info.surfaceInfo.position, dir);
            }
            else
            {
                break;
            }
        }

        return radiance;
    }
};