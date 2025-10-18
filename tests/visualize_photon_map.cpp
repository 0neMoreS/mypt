#include <omp.h>

#include "camera.h"
#include "image.h"
#include "Intergrator/photo_mapping/photon_mapping_integrator.h"
#include "scene.h"

int main() {
  const int width = 512;
  const int height = 512;
  const int n_photons = 100000;
  const int max_depth = 8;

  Image image(width, height);
  Camera camera(Vec3f(0, 1.f, 3), Vec3f(0, 1.f, -1), Vec3f(0, 1.f, 0), 50.f, float(width) / float(height));

  Scene scene;
  scene.loadModel("../models/cornellbox-water2.obj");
  scene.build();

  // photon tracing and build photon map
  PhotonMapping integrator(n_photons, 1, 0, 0, 0, max_depth);
  UniformSampler sampler;
  integrator.build(scene, sampler);

  const PhotonMap* photon_map = integrator.getPhotonMapPtr();

#pragma omp parallel for collapse(2) schedule(dynamic, 1)
  for (int i = 0; i < height; ++i) {
    for (int j = 0; j < width; ++j) {
      Ray ray;
      float pdf;
      camera.sampleRay(Vec2f((float(j) / width), (float(i) / height)), ray, pdf);
      IntersectInfo info;
      if (scene.intersect(ray, info)) {
        // query photon map
        float r2;
        const std::vector<int> photon_indices = photon_map->queryKNearestPhotons(info.surfaceInfo.position, 1, r2);
        const int photon_idx = photon_indices[0];

        // if distance to the photon is small enough, write photon's
        // throughput to the image
        if (r2 < 0.001f) {
          const Photon& photon = photon_map->getIthPhoton(photon_idx);
          image.setPixel((height - i - 1), (width - j - 1), photon.throughput);
        }
      } else {
        image.setPixel((height - i - 1), (width - j - 1), Vec3f(0));
      }

    }
  }

  image.writePPM("output_visualize_photon_map.ppm");

  return 0;
}