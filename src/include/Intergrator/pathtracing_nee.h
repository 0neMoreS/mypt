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
        Vec3f throughput(1, 1, 1); // the throughput of previous shading point, because the different cos in Sample Light and Sample BRDF
        IntersectInfo info; // current shading point
        IntersectInfo info_next; // next shading point
        bool hit_next = scene.intersect(wo, info_next);
        for (int k = 0; k < maxDepth; ++k)
        {
            if (hit_next)
            {
                // update info for current shading point
                wo = -wo;
                info = info_next;
                hit_next = false;

                // directly visible light source
                if (k == 0 && info.hasAreaLight()) {
                    radiance += throughput * info.Le(wo.direction);
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

                // Get pdfs
                Vec3f dir;
                float pdf_dir;
                Vec3f f_brdf = info.sampleBxDF(wo.direction, TransportDirection::FROM_CAMERA, sampler, dir, pdf_dir);

                // sample light from the scene
                float scene_lights_pdf;
                std::shared_ptr<Light> light = scene.sampleLight(sampler, scene_lights_pdf);

                // NEE from this shading point
                if (info.getBxDFType() == BxDFType::DIFFUSE) {
                    float light_pdf;
                    SurfaceInfo light_surf = light->samplePoint(sampler, light_pdf);
                    const Vec3f to_light = normalize(light_surf.position - info.surfaceInfo.position);

                    // ============ Sample Light ===========

                    if (scene.intersect(Ray{ info.surfaceInfo.position, to_light }, info_next) && info_next.hasAreaLight(light)) {
                        const float cos_theta_light = std::abs(dot(-to_light, info_next.surfaceInfo.shadingNormal));
                        const float dist2 = std::pow(length(info_next.surfaceInfo.position - info.surfaceInfo.position), 2.f);
                        const float w_light_pdf = light_pdf * dist2 / cos_theta_light;
                        const float weight_light = w_light_pdf / (w_light_pdf + pdf_dir); // balance heuristic

                        const Vec3f f_light = info.evaluateBxDF(wo.direction, to_light, TransportDirection::FROM_CAMERA);

                        radiance += throughput * weight_light * (f_light * cosTerm(wo.direction, to_light, info.surfaceInfo, TransportDirection::FROM_CAMERA)) * info_next.Le(-to_light) / (w_light_pdf * scene_lights_pdf);
                    }

                    // ============ Sample Light ============

                    // ============ Sample BRDF ============

                    // if hit, info_next will update to cover the previous info_next, if not hit, hit_next remains false and loop will break
                    if (scene.intersect(Ray{ info.surfaceInfo.position, dir }, info_next))
                    {
                        hit_next = true;
                        if (info_next.hasAreaLight(light)) {

                            const float cos_theta_light = std::abs(dot(-dir, info_next.surfaceInfo.shadingNormal));
                            const float dist2 = std::pow(length(info.surfaceInfo.position - info_next.surfaceInfo.position), 2.f);
                            const float w_light_pdf = light_pdf * dist2 / cos_theta_light;
                            const float weight_light = pdf_dir / (w_light_pdf + pdf_dir); // balance heuristic

                            radiance += throughput * weight_light * (f_brdf * cosTerm(wo.direction, dir, info.surfaceInfo, TransportDirection::FROM_CAMERA)) * info_next.Le(-dir) / (pdf_dir * scene_lights_pdf);
                        }
                    }

                    // ============ Sample BRDF ============
                } else if (info.getBxDFType() == BxDFType::SPECULAR) {
                    if (scene.intersect(Ray{ info.surfaceInfo.position, dir }, info_next))
                    {
                        hit_next = true;
                        if (info_next.hasAreaLight(light)) {
                            radiance += throughput * (f_brdf * cosTerm(wo.direction, dir, info.surfaceInfo, TransportDirection::FROM_CAMERA)) * light->Le(info_next.surfaceInfo, -dir) / (pdf_dir * scene_lights_pdf);
                        }
                    }
                }

                // update ray and throughput
                throughput *= f_brdf * cosTerm(wo.direction, dir, info.surfaceInfo, TransportDirection::FROM_CAMERA) / pdf_dir;
                wo = Ray(info.surfaceInfo.position, dir);
            } else
            {
                break;
            }
        }

        return radiance;
    }
};