#pragma once
#include <omp.h>

#include <optional>

#include "core.h"
#include "scene.h"

class Integrator
{
public:
  // do preliminary jobs before calling integrate
  virtual void build(const Scene& scene, Sampler& sampler) = 0;

  // compute radiance coming from the given ray
  virtual Vec3f integrate(const Ray& ray, const Scene& scene, Sampler& sampler) const = 0;

  // compute cosine term
  // NOTE: need to account for the asymmetry of BSDF when photon tracing
  // https://pbr-book.org/3ed-2018/Light_Transport_III_Bidirectional_Methods/The_Path-Space_Measurement_Equation#x3-Non-symmetryDuetoShadingNormals
  // Veach, Eric. Robust Monte Carlo methods for light transport simulation.
  // Stanford University, 1998. Section 5.3
  static float cosTerm(const Vec3f& wo, const Vec3f& wi, const SurfaceInfo& surfaceInfo, const TransportDirection& transport_dir)
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