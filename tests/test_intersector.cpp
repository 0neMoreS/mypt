#include "camera.h"
#include "image.h"
#include "scene.h"

int main()
{
  const int width = 400, height = 400;
  const int n_samples = 32;
  Image image(width, height);
  // Camera camera(Vec3f(0, 1.f, 2.f), Vec3f(0, 1.f, -1.f), Vec3f(0, 1, 0), 114.f, float(width) / float(height), 0.1f, 100.f);
  Camera camera(Vec3f(0, 1, 7), Vec3f(0, 0, -1), 0.25f * PI);

  Scene scene;
  scene.loadModel("../models/CornellBox-Original.obj");
  scene.build();

  // #pragma omp parallel for schedule(dynamic, 1)
  //   for (unsigned short y = 0; y < h; y++)
  //   {
  //     fprintf(stderr, "\rRendering (%d spp) %5.2f%%", n_samples, 100. * y / (h - 1));
  //     for (unsigned short x = 0; x < w; x++)
  //     {
  //       UniformSampler sampler((y * w + x) * 9781 + 1);
  //       for (unsigned short s = 0; s < n_samples; s++)
  //       {
  //         double dx = sampler.getNext1D();
  //         double dy = sampler.getNext1D();
  //         Ray r;
  //         camera.sampleRay(Vec2f((x + dx) / w, (y + dy) / h), r);
  //         IntersectInfo info;
  //         if (scene.intersect(r, info))
  //         {
  //           image.setPixel(y, x, 0.5f * (info.surfaceInfo.shadingNormal + 1.0f));
  //         }
  //         else
  //         {
  //           image.setPixel(y, x, Vec3f(0.0f));
  //         }
  //       }
  //     }
  //   }

  for (int i = 0; i < height; ++i)
  {
    for (int j = 0; j < width; ++j)
    {
      const float u = (2.0f * j - width) / height;
      const float v = (2.0f * i - height) / height;

      Ray ray;
      float pdf;
      if (camera.sampleRay(Vec2f(u, v), ray, pdf))
      {
        IntersectInfo info;
        if (scene.intersect(ray, info))
        {
          image.setPixel(i, j, 0.5f * (info.surfaceInfo.shadingNormal + 1.0f));
        }
      }
      else
      {
        image.setPixel(i, j, Vec3f(0));
      }
    }
  }

  image.writePPM("output.ppm");

  return 0;
}