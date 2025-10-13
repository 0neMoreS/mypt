#include "camera.h"
#include "image.h"
#include "scene.h"

int main()
{
  const int w = 400, h = 400;
  const int n_samples = 32;
  Image image(w, h);
  Camera camera(Vec3f(0, 0, 5), Vec3f(0, 0, 0), Vec3f(0, 1, 0), 45.0, float(w) / float(h), 1.0, 100.0);

  Scene scene;
  scene.loadModel("../models/CornellBox-Original.obj");
  scene.build();

#pragma omp parallel for schedule(dynamic, 1)
  for (unsigned short y = 0; y < h; y++)
  {
    fprintf(stderr, "\rRendering (%d spp) %5.2f%%", n_samples, 100. * y / (h - 1));
    for (unsigned short x = 0; x < w; x++)
    {
      UniformSampler sampler((y * w + x) * 9781 + 1);
      for (unsigned short s = 0; s < n_samples; s++)
      {
        double dx = sampler.getNext1D();
        double dy = sampler.getNext1D();
        Ray r;
        camera.getRay(Vec2f((x + dx) / w, (y + dy) / h), r);
        IntersectInfo info;
        if (scene.intersect(r, info))
        {
          image.setPixel(y, x, 0.5f * (info.surfaceInfo.shadingNormal + 1.0f));
        }
        else
        {
          image.setPixel(y, x, Vec3f(0.0f));
        }
      }
    }
  }

  image.writePPM("output.ppm");

  return 0;
}