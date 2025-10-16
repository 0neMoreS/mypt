#pragma once
#include "integrator.h"
// implementation of path tracing with next event estimation
class PathTracingNEE : public Integrator
{
private:
    const int maxDepth;

public:
    PathTracingNEE(int maxDepth = 100) : maxDepth(maxDepth) {}

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
                // russian roulette
                if (k >= maxDepth / 2) {
                    const float russian_roulette_prob = std::min(std::max(throughput[0], std::max(throughput[1], throughput[2])), 1.0f);
                    if (sampler.getNext1D() >= russian_roulette_prob)
                    {
                        break;
                    }
                    throughput /= russian_roulette_prob;
                }

                // ============ Direct Lighting ============

                // NEE
                float scene_lights_pdf;
                std::shared_ptr<Light> light = scene.sampleLight(sampler, scene_lights_pdf);

                float light_pdf;
                SurfaceInfo light_surf = light->samplePoint(sampler, light_pdf);
                const Vec3f to_light = normalize(light_surf.position - info.surfaceInfo.position);

                IntersectInfo nee_info;
                if (scene.intersect(Ray{ info.surfaceInfo.position, to_light }, nee_info) && nee_info.hitPrimitive->hasAreaLight()) {
                    const float cos_theta_light = std::max(0.0f, dot(-to_light, light_surf.shadingNormal));
                    const float cos_theta_surface = std::max(0.0f, dot(to_light, info.surfaceInfo.shadingNormal));
                    const float dist2 = length(light_surf.position - info.surfaceInfo.position);

                    // evaluate BSDF
                    const Vec3f f = info.hitPrimitive->evaluateBxDF(-ray.direction, to_light, info.surfaceInfo, TransportDirection::FROM_CAMERA);
                    if (f[0] > 0 || f[1] > 0 || f[2] > 0) {
                        const Vec3f Le = light->Le(light_surf, -to_light);
                        radiance += throughput * f * Le * cos_theta_light * cos_theta_surface / (dist2 * scene_lights_pdf * light_pdf);
                    }
                }

                // ============ Direct Lighting ============

                // ============ Indirect Lighting ============

                // Hit Light
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

                // ============ Indirect Lighting ============
            }
            else
            {
                break;
            }
        }

        return radiance * 0.5f;
    }
};