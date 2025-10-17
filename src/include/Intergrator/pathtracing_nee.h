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
        Ray wo = ray_in;
        Vec3f throughput(1, 1, 1);
        for (int k = 0; k < maxDepth; ++k)
        {
            IntersectInfo info;
            if (scene.intersect(wo, info))
            {
                wo = -wo;

                if (k == 0 && info.hitPrimitive->hasAreaLight()) {
                    // directly visible light source
                    radiance += throughput * info.hitPrimitive->Le(info.surfaceInfo, wo.direction);
                    // pdf here is camera's pdf, which is divided in the main loop
                }

                if (k >= maxDepth / 2) {
                    // russian roulette
                    const float russian_roulette_prob = std::min(std::max(throughput[0], std::max(throughput[1], throughput[2])), 1.0f);
                    if (sampler.getNext1D() >= russian_roulette_prob)
                    {
                        break;
                    }
                    throughput /= russian_roulette_prob;
                }

                // NEE from this shading point
                float scene_lights_pdf;
                std::shared_ptr<Light> light = scene.sampleLight(sampler, scene_lights_pdf);
                radiance += (throughput * NEE_MIS(wo, info, light, scene, sampler) / scene_lights_pdf);


                // sample direction by BxDF, change the NEE shading point
                Vec3f dir;
                float pdf_dir;
                Vec3f f = info.hitPrimitive->sampleBxDF(wo.direction, info.surfaceInfo, TransportDirection::FROM_CAMERA, sampler, dir, pdf_dir);

                // update throughput and ray
                throughput *= f * cosTerm(wo.direction, dir, info.surfaceInfo, TransportDirection::FROM_CAMERA) / pdf_dir;
                wo = Ray(info.surfaceInfo.position, dir);
            }
            else
            {
                break;
            }
        }

        return radiance;
    }

    Vec3f NEE_MIS(const Ray& wo, const IntersectInfo& info, const std::shared_ptr<Light>& light, const Scene& scene, Sampler& sampler) const
    {
        Vec3f radiance(0);

        // Get pdfs
        Vec3f dir;
        float dir_pdf;
        Vec3f f_brdf = info.hitPrimitive->sampleBxDF(wo.direction, info.surfaceInfo, TransportDirection::FROM_CAMERA, sampler, dir, dir_pdf);

        if (info.hitPrimitive->getBxDFType() == BxDFType::DIFFUSE) {
            float light_pdf;
            SurfaceInfo light_surf = light->samplePoint(sampler, light_pdf);
            const Vec3f to_light = normalize(light_surf.position - info.surfaceInfo.position);
            IntersectInfo light_info;

            // ============ Sample Light ===========

            if (scene.intersect(Ray{ info.surfaceInfo.position, to_light }, light_info) && light_info.hitPrimitive->hasAreaLight(light)) {
                const float cos_theta_light = std::abs(dot(-to_light, light_info.surfaceInfo.shadingNormal));
                const float dist2 = std::pow(length(light_info.surfaceInfo.position - info.surfaceInfo.position), 2.f);
                const float w_light_pdf = light_pdf * dist2 / cos_theta_light;
                const float weight_light = w_light_pdf / (w_light_pdf + dir_pdf); // balance heuristic

                const Vec3f f_light = info.hitPrimitive->evaluateBxDF(wo.direction, to_light, info.surfaceInfo, TransportDirection::FROM_CAMERA);

                radiance += weight_light * (f_light * cosTerm(wo.direction, to_light, info.surfaceInfo, TransportDirection::FROM_CAMERA)) * light_info.hitPrimitive->Le(light_info.surfaceInfo, -to_light) / w_light_pdf;
            }

            // ============ Sample Light ============

            // ============ Sample BRDF ============

            if (scene.intersect(Ray{ info.surfaceInfo.position, dir }, light_info) && light_info.hitPrimitive->hasAreaLight(light))
            {
                const float cos_theta_light = std::abs(dot(-dir, light_info.surfaceInfo.shadingNormal));
                const float dist2 = std::pow(length(info.surfaceInfo.position - light_info.surfaceInfo.position), 2.f);
                const float w_light_pdf = light_pdf * dist2 / cos_theta_light;
                const float weight_light = dir_pdf / (w_light_pdf + dir_pdf); // balance heuristic

                radiance += weight_light * (f_brdf * cosTerm(wo.direction, dir, info.surfaceInfo, TransportDirection::FROM_CAMERA)) * light_info.hitPrimitive->Le(light_info.surfaceInfo, -dir) / dir_pdf;
            }

            // ============ Sample BRDF ============
        }
        else {
            IntersectInfo light_info;
            if (scene.intersect(Ray{ info.surfaceInfo.position, dir }, light_info) && light_info.hitPrimitive->hasAreaLight(light))
            {
                radiance += (f_brdf * cosTerm(wo.direction, dir, info.surfaceInfo, TransportDirection::FROM_CAMERA)) * light->Le(light_info.surfaceInfo, -dir) / dir_pdf;
            }
        }
        return radiance;
    }
};