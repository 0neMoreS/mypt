#pragma once
#include "Intergrator/integrator.h"
#include "Intergrator/photo_mapping/photo_mapping_core.h"

// implementation of photon mapping
class PhotonMapping : public Integrator
{
private:
    // number of photons used for making global photon map
    const int nPhotonsGlobal;

    // number of photons used for radiance estimation by global photon map
    const int nEstimationGlobal;

    // number of photons for making caustics photon map
    const int nPhotonsCaustics;

    // number of photons used for radiance estimation by caustics photon map
    const int nEstimationCaustics;

    // maximum depth to estimate radiance by final gathering
    const int finalGatheringDepth;

    // maximum depth of photon tracing, eye tracing
    const int maxDepth;

    PhotonMap globalPhotonMap;
    PhotonMap causticsPhotonMap;

    // compute reflected radiance with global photon map
    Vec3f computeRadianceWithPhotonMap(const Vec3f& wo,
        const IntersectInfo& info) const {
        // get nearby photons
        float max_dist2;
        const std::vector<int> photon_indices = globalPhotonMap.queryKNearestPhotons(info.surfaceInfo.position, nEstimationGlobal, max_dist2);

        Vec3f Lo;
        for (const int photon_idx : photon_indices) {
            const Photon& photon = globalPhotonMap.getIthPhoton(photon_idx);
            const Vec3f f = info.hitPrimitive->evaluateBxDF(wo, photon.wi, info.surfaceInfo, TransportDirection::FROM_CAMERA);
            Lo += f * photon.throughput;
        }
        if (photon_indices.size() > 0) {
            Lo /= (nPhotonsGlobal * PI * max_dist2);
        }
        return Lo;
    }

    // compute reflected radiance with caustics photon map
    Vec3f computeCausticsWithPhotonMap(const Vec3f& wo, const IntersectInfo& info) const {
        // get nearby photons
        float max_dist2;
        const std::vector<int> photon_indices = causticsPhotonMap.queryKNearestPhotons(info.surfaceInfo.position, nEstimationGlobal, max_dist2);

        Vec3f Lo;
        for (const int photon_idx : photon_indices) {
            const Photon& photon = causticsPhotonMap.getIthPhoton(photon_idx);
            const Vec3f f = info.hitPrimitive->evaluateBxDF(wo, photon.wi, info.surfaceInfo, TransportDirection::FROM_CAMERA);
            Lo += f * photon.throughput;
        }
        if (photon_indices.size() > 0) {
            Lo /= (nPhotonsCaustics * PI * max_dist2);
        }

        return Lo;
    }

    // sample initial ray from light and compute initial throughput
    Ray sampleRayFromLight(const Scene& scene, Sampler& sampler, Vec3f& throughput) {
        // sample light
        float light_choose_pdf;
        const std::shared_ptr<Light> light = scene.sampleLight(sampler, light_choose_pdf);

        // sample point on light
        float light_pos_pdf;
        const SurfaceInfo light_surf = light->samplePoint(sampler, light_pos_pdf);

        // sample direction on light
        float light_dir_pdf;
        const Vec3f dir = light->sampleDirection(light_surf, sampler, light_dir_pdf);

        // spawn ray
        Ray ray(light_surf.position, dir);
        throughput = light->Le(light_surf, dir) / (light_choose_pdf * light_pos_pdf * light_dir_pdf) * std::abs(dot(dir, light_surf.shadingNormal));

        return ray;
    }

public:
    PhotonMapping(int nPhotonsGlobal, int nEstimationGlobal, float nPhotonsCausticsMultiplier, int nEstimationCaustics, int strictCalcDepth, int maxDepth) :
        nPhotonsGlobal(nPhotonsGlobal),
        nEstimationGlobal(nEstimationGlobal),
        nPhotonsCaustics(nPhotonsGlobal* nPhotonsCausticsMultiplier),
        nEstimationCaustics(nEstimationCaustics),
        finalGatheringDepth(strictCalcDepth),
        maxDepth(maxDepth) {
    }

    const PhotonMap* getPhotonMapPtr() const { return &globalPhotonMap; }

    // photon tracing and build photon map
    void build(const Scene& scene, Sampler& sampler) override {
        std::vector<Photon> thread_photons[omp_get_max_threads()];
        std::vector<Photon> photons;

        // init sampler for each thread
        std::vector<std::unique_ptr<Sampler>> samplers(omp_get_max_threads());
        for (int i = 0; i < samplers.size(); ++i) {
            samplers[i] = sampler.clone();
            samplers[i]->setSeed(samplers[i]->getSeed() * (i + 1));
        }

        // build global photon map
        // photon tracing
#pragma omp parallel for
        for (int i = 0; i < nPhotonsGlobal; ++i) {
            auto& sampler_per_thread = *samplers[omp_get_thread_num()];

            // sample initial ray from light and set initial throughput
            Vec3f throughput;
            Ray ray = sampleRayFromLight(scene, sampler_per_thread, throughput);

            // trace photons
            // whener hitting diffuse surface, add photon to the photon array
            // recursively tracing photon with russian roulette
            for (int k = 0; k < maxDepth; ++k) {
                if (std::isnan(throughput[0]) || std::isnan(throughput[1]) || std::isnan(throughput[2])) {
                    break;
                } else if (throughput[0] < 0 || throughput[1] < 0 || throughput[2] < 0) {
                    break;
                }

                IntersectInfo info;
                if (scene.intersect(ray, info)) {
                    const BxDFType bxdf_type = info.hitPrimitive->getBxDFType();
                    if (bxdf_type == BxDFType::DIFFUSE) {
                        {
                            thread_photons[omp_get_thread_num()].emplace_back(throughput, info.surfaceInfo.position, -ray.direction);
                        }
                    }

                    // russian roulette
                    if (k > 0) {
                        const float russian_roulette_prob = std::min(std::max(throughput[0], std::max(throughput[1], throughput[2])), 1.0f);
                        if (sampler_per_thread.getNext1D() >= russian_roulette_prob) {
                            break;
                        }
                        throughput /= russian_roulette_prob;
                    }

                    // sample direction by BxDF
                    Vec3f dir;
                    float pdf_dir;
                    const Vec3f f = info.hitPrimitive->sampleBxDF(-ray.direction, info.surfaceInfo, TransportDirection::FROM_LIGHT, sampler_per_thread, dir, pdf_dir);

                    // update throughput and ray
                    throughput *= f * cosTerm(-ray.direction, dir, info.surfaceInfo, TransportDirection::FROM_LIGHT) / pdf_dir;
                    ray = Ray(info.surfaceInfo.position, dir);
                } else {
                    // photon goes to the sky
                    break;
                }
            }
        }

        for (int i = 0; i < omp_get_max_threads(); ++i) {
            for (const auto& photon : thread_photons[i]) {
                photons.push_back(photon);
            }
        }

        // build photon map
        globalPhotonMap.setPhotons(photons);
        globalPhotonMap.build();

        // build caustics photon map
        if (finalGatheringDepth > 0) {
            photons.clear();
            for (int i = 0; i < omp_get_max_threads(); ++i) {
                thread_photons[i].clear();
            }
            // photon tracing
#pragma omp parallel for
            for (int i = 0; i < nPhotonsCaustics; ++i) {
                auto& sampler_per_thread = *samplers[omp_get_thread_num()];

                // sample initial ray from light and set initial throughput
                Vec3f throughput;
                Ray ray = sampleRayFromLight(scene, sampler_per_thread, throughput);

                // when hitting diffuse surface after specular, add photon to the photon
                // array
                bool prev_specular = false;
                for (int k = 0; k < maxDepth; ++k) {
                    if (std::isnan(throughput[0]) || std::isnan(throughput[1]) ||
                        std::isnan(throughput[2])) {
                        break;
                    } else if (throughput[0] < 0 || throughput[1] < 0 ||
                        throughput[2] < 0) {
                        break;
                    }

                    IntersectInfo info;
                    if (scene.intersect(ray, info)) {
                        const BxDFType bxdf_type = info.hitPrimitive->getBxDFType();

                        // break when hitting diffuse surface without previous specular
                        if (!prev_specular && bxdf_type == BxDFType::DIFFUSE) {
                            break;
                        }

                        // add photon when hitting diffuse surface after specular
                        if (prev_specular && bxdf_type == BxDFType::DIFFUSE) {
                            {
                                thread_photons[omp_get_thread_num()].emplace_back(throughput, info.surfaceInfo.position, -ray.direction);
                            }
                            break;
                        }

                        prev_specular = (bxdf_type == BxDFType::SPECULAR);

                        // russian roulette
                        if (k > 0) {
                            const float russian_roulette_prob =
                                std::min(std::max(throughput[0],
                                    std::max(throughput[1], throughput[2])),
                                    1.0f);
                            if (sampler_per_thread.getNext1D() >= russian_roulette_prob) {
                                break;
                            }
                            throughput /= russian_roulette_prob;
                        }

                        // sample direction by BxDF
                        Vec3f dir;
                        float pdf_dir;
                        const Vec3f f = info.hitPrimitive->sampleBxDF(-ray.direction, info.surfaceInfo, TransportDirection::FROM_LIGHT, sampler_per_thread, dir, pdf_dir);

                        // update throughput and ray
                        throughput *= f * cosTerm(-ray.direction, dir, info.surfaceInfo, TransportDirection::FROM_LIGHT) / pdf_dir;
                        ray = Ray(info.surfaceInfo.position, dir);
                    } else {
                        // photon goes to the sky
                        break;
                    }
                }
            }

            for (int i = 0; i < omp_get_max_threads(); ++i) {
                for (const auto& photon : thread_photons[i]) {
                    photons.push_back(photon);
                }
            }

            causticsPhotonMap.setPhotons(photons);
            causticsPhotonMap.build();
        }
    }

    Vec3f integrate(const Ray& ray_in, const Scene& scene, Sampler& sampler) const override {
        Vec3f radiance(0);
        Ray wo = ray_in;
        Vec3f throughput(1, 1, 1); // the throughput of previous shading point, because the different cos in Sample Light and Sample BRDF
        IntersectInfo info; // current shading point
        for (int k = 0; k < maxDepth; ++k)
        {
            if (scene.intersect(wo, info))
            {
                // update info for current shading point
                wo = -wo;

                // directly visible light source
                if (info.hitPrimitive->hasAreaLight()) {
                    radiance += throughput * info.hitPrimitive->Le(info.surfaceInfo, wo.direction);
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

                // NEE for DIFFUSE
                if (info.hitPrimitive->getBxDFType() == BxDFType::DIFFUSE) {
                    if (k > finalGatheringDepth) {
                        radiance += throughput * computeRadianceWithPhotonMap(wo.direction, info);
                        break;
                    }

                    radiance += throughput * computeDirectIllumination(scene, wo.direction, info, sampler);
                    radiance += throughput * computeIndirectWithPhotonMap(scene, wo.direction, info, sampler);
                    radiance += throughput * computeCausticsWithPhotonMap(wo.direction, info);

                    break;
                } else if (info.hitPrimitive->getBxDFType() == BxDFType::SPECULAR) {
                    // Get pdfs
                    Vec3f dir;
                    float pdf_dir;
                    Vec3f f_brdf = info.hitPrimitive->sampleBxDF(wo.direction, info.surfaceInfo, TransportDirection::FROM_CAMERA, sampler, dir, pdf_dir);

                    // update ray and throughput
                    throughput *= f_brdf * cosTerm(wo.direction, dir, info.surfaceInfo, TransportDirection::FROM_CAMERA) / pdf_dir;
                    wo = Ray(info.surfaceInfo.position, dir);
                }
            } else
            {
                break;
            }
        }

        return radiance;
    }

private:

    // compute direct illumination with explicit light sampling(NEE)
    Vec3f computeDirectIllumination(const Scene& scene, const Vec3f& wo, const IntersectInfo& info, Sampler& sampler) const {
        Vec3f radiance(0);

        // sample light from the scene
        float scene_lights_pdf;
        std::shared_ptr<Light> light = scene.sampleLight(sampler, scene_lights_pdf);

        // sample light position
        float light_pdf;
        SurfaceInfo light_surf = light->samplePoint(sampler, light_pdf);
        const Vec3f to_light = normalize(light_surf.position - info.surfaceInfo.position);

        // sample brdf
        Vec3f dir;
        float pdf_dir;
        const Vec3f f_brdf = info.hitPrimitive->sampleBxDF(wo, info.surfaceInfo, TransportDirection::FROM_CAMERA, sampler, dir, pdf_dir);

        IntersectInfo info_next;

        // ============ Sample Light ===========

        if (scene.intersect(Ray{ info.surfaceInfo.position, to_light }, info_next) && info_next.hitPrimitive->hasAreaLight(light)) {
            const float cos_theta_light = std::abs(dot(-to_light, info_next.surfaceInfo.shadingNormal));
            const float dist2 = std::pow(length(info_next.surfaceInfo.position - info.surfaceInfo.position), 2.f);
            const float w_light_pdf = light_pdf * dist2 / cos_theta_light;
            const float weight_light = w_light_pdf / (w_light_pdf + pdf_dir); // balance heuristic

            const Vec3f f_light = info.hitPrimitive->evaluateBxDF(wo, to_light, info.surfaceInfo, TransportDirection::FROM_CAMERA);

            radiance += weight_light * (f_light * cosTerm(wo, to_light, info.surfaceInfo, TransportDirection::FROM_CAMERA)) * info_next.hitPrimitive->Le(info_next.surfaceInfo, -to_light) / (w_light_pdf * scene_lights_pdf);
        }

        // ============ Sample Light ============

        // ============ Sample BRDF ============

        // if hit, info_next will update to cover the previous info_next, if not hit, hit_next remains false and loop will break
        if (scene.intersect(Ray{ info.surfaceInfo.position, dir }, info_next) && info_next.hitPrimitive->hasAreaLight(light))
        {
            const float cos_theta_light = std::abs(dot(-dir, info_next.surfaceInfo.shadingNormal));
            const float dist2 = std::pow(length(info.surfaceInfo.position - info_next.surfaceInfo.position), 2.f);
            const float w_light_pdf = light_pdf * dist2 / cos_theta_light;
            const float weight_light = pdf_dir / (w_light_pdf + pdf_dir); // balance heuristic

            radiance += weight_light * (f_brdf * cosTerm(wo, dir, info.surfaceInfo, TransportDirection::FROM_CAMERA)) * info_next.hitPrimitive->Le(info_next.surfaceInfo, -dir) / (pdf_dir * scene_lights_pdf);
        }

        return radiance;
    }

    Vec3f computeIndirectWithPhotonMap(const Scene& scene, const Vec3f& ray_in, const IntersectInfo& info, Sampler& sampler) const {
        Vec3f radiance(0);
        Vec3f throughput(1, 1, 1);
        IntersectInfo info_current = info;
        Ray wo = Ray(info_current.surfaceInfo.position, ray_in);

        for (int k = 0; k < maxDepth; ++k) {
            Vec3f dir;
            float pdf_dir;
            const Vec3f f = info_current.hitPrimitive->sampleBxDF(wo.direction, info_current.surfaceInfo, TransportDirection::FROM_CAMERA, sampler, dir, pdf_dir);
            IntersectInfo info_next;

            if (scene.intersect(Ray{ info_current.surfaceInfo.position, dir }, info_next)) {

                // RR
                if (k >= maxDepth / 2) {
                    // russian roulette
                    const float russian_roulette_prob = std::min(std::max(throughput[0], std::max(throughput[1], throughput[2])), 1.0f);
                    if (sampler.getNext1D() >= russian_roulette_prob)
                    {
                        break;
                    }
                    throughput /= russian_roulette_prob;
                }

                if (info_next.hitPrimitive->getBxDFType() == BxDFType::DIFFUSE) {
                    // when hitting diffuse, compute radiance with photon map
                    radiance += throughput * f * cosTerm(wo.direction, dir, info_current.surfaceInfo, TransportDirection::FROM_CAMERA) * computeRadianceWithPhotonMap(-dir, info_next) / pdf_dir;
                    break;
                } else if (info_next.hitPrimitive->getBxDFType() == BxDFType::SPECULAR) {
                    // when hitting specular, next hit
                    wo = Ray(info_next.surfaceInfo.position, -dir);
                    throughput *= f * cosTerm(wo.direction, dir, info_current.surfaceInfo, TransportDirection::FROM_CAMERA) / pdf_dir;
                    info_current = info_next;
                }
            } else {
                break;
            }
        }

        return radiance;
    }
};