
#include <chrono>
#include <omp.h>

#include "camera.h"
#include "image.h"
#include "Intergrator/pathtracing.h"
#include "Intergrator/pathtracing_nee.h"
#include "scene.h"

int main()
{

  const int width = 512;
  const int height = 512;
  const int n_samples = 1024;
  const int max_depth = 8;


  auto start = std::chrono::high_resolution_clock::now();

  Image image(width, height);
  Camera camera(Vec3f(0, 1.f, 3), Vec3f(0, 1.f, -1), Vec3f(0, 1.f, 0), 50.f, float(width) / float(height));

  Scene scene;
  scene.loadModel("../models/CornellBox-Original.obj");
  // scene.loadModel("../models/cornellbox-water2.obj");
  scene.build();

  // photon tracing and build photon map
  PathTracing integrator(max_depth);

#pragma omp parallel for collapse(2) schedule(dynamic, 1)
  for (int i = 0; i < height; ++i)
  {
    for (int j = 0; j < width; ++j)
    {
      fprintf(stderr, "\rRendering (%d spp) %5.2f%%", n_samples, 100. * i / (height - 1));
      UniformSampler sampler(j + width * i);

      for (int k = 0; k < n_samples; ++k)
      {
        Ray ray;
        float pdf;
        float dx = sampler.getNext1D();
        float dy = sampler.getNext1D();
        camera.sampleRay(Vec2f(float(j + dx) / width, float(i + dy) / height), ray, pdf);
        const Vec3f radiance = integrator.integrate(ray, scene, sampler) / pdf;

        if (std::isnan(radiance[0]) || std::isnan(radiance[1]) ||
          std::isnan(radiance[2]))
        {
          continue;
        }

        image.addPixel((height - i - 1), (width - j - 1), radiance);
      }
    }
  }

  // take average
  image.divide(n_samples);

  image.gammaCorrection(2.2f);
  image.writePPM("output_test_path_tracing.ppm");

  auto end = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = end - start;
  printf("\nTotal rendering time: %.2f hours (%.2f seconds)\n", elapsed.count() / 3600.0, elapsed.count());
}