#pragma once
#include <omp.h>

#include <optional>

#include "core.h"
#include "scene.h"

class Integrator
{
public:
  // do preliminary jobs before calling integrate
  virtual void build(const Scene &scene, Sampler &sampler) = 0;

  // compute radiance coming from the given ray
  virtual Vec3f integrate(const Ray &ray, const Scene &scene,
                          Sampler &sampler) const = 0;

  // compute cosine term
  // NOTE: need to account for the asymmetry of BSDF when photon tracing
  // https://pbr-book.org/3ed-2018/Light_Transport_III_Bidirectional_Methods/The_Path-Space_Measurement_Equation#x3-Non-symmetryDuetoShadingNormals
  // Veach, Eric. Robust Monte Carlo methods for light transport simulation.
  // Stanford University, 1998. Section 5.3
  static float cosTerm(const Vec3f &wo, const Vec3f &wi,
                       const SurfaceInfo &surfaceInfo,
                       const TransportDirection &transport_dir)
  {
    const float wi_ns = dot(wi, surfaceInfo.shadingNormal);
    const float wi_ng = dot(wi, surfaceInfo.geometricNormal);
    const float wo_ns = dot(wo, surfaceInfo.shadingNormal);
    const float wo_ng = dot(wo, surfaceInfo.geometricNormal);

    // prevent light leaks
    if (wi_ng * wi_ns <= 0 || wo_ng * wo_ns <= 0)
    {
      return 0;
    }

    if (transport_dir == TransportDirection::FROM_CAMERA)
    {
      return std::abs(wi_ns);
    }
    else if (transport_dir == TransportDirection::FROM_LIGHT)
    {
      return std::abs(wo_ns) * std::abs(wi_ng) / std::abs(wo_ng);
    }
    else
    {
      std::exit(EXIT_FAILURE);
    }
  }
};

// implementation of path tracing
// NOTE: for reference purpose
class PathTracing : public Integrator
{
private:
  const int maxDepth;

public:
  PathTracing(int maxDepth = 100) : maxDepth(maxDepth) {}

  void build(const Scene &scene, Sampler &sampler) override {}

  Vec3f integrate(const Ray &ray_in, const Scene &scene,
                  Sampler &sampler) const override
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
        if (k > 0)
        {
          const float russian_roulette_prob = std::min(
              std::max(throughput[0], std::max(throughput[1], throughput[2])),
              1.0f);
          if (sampler.getNext1D() >= russian_roulette_prob)
          {
            break;
          }
          throughput /= russian_roulette_prob;
        }

        // Le
        if (info.hitPrimitive->hasAreaLight())
        {
          radiance += throughput *
                      info.hitPrimitive->Le(info.surfaceInfo, -ray.direction);
        }

        // sample direction by BxDF
        Vec3f dir;
        float pdf_dir;
        Vec3f f = info.hitPrimitive->sampleBxDF(
            -ray.direction, info.surfaceInfo, TransportDirection::FROM_CAMERA,
            sampler, dir, pdf_dir);

        // update throughput and ray
        throughput *= f *
                      cosTerm(-ray.direction, dir, info.surfaceInfo,
                              TransportDirection::FROM_CAMERA) /
                      pdf_dir;
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