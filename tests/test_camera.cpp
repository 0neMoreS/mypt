#include "camera.h"
#include "image.h"
#include "core.h"
#include "sampler.h"

#ifndef OLD
int main()
{
  const int height = 512, width = 512;
  const int n_samples = 32;
  Image image(width, height);
  Camera camera(Vec3f(0, 0, 0), Vec3f(0, 0, -1.f), Vec3f(0, 1.f, 0), 90.f, float(width) / float(height));

#pragma omp parallel for schedule(dynamic, 1)
  for (unsigned short i = 0; i < height; i++)
  {
    fprintf(stderr, "\rRendering (%d spp) %5.2f%%", n_samples, 100. * i / (height - 1));
    for (unsigned short j = 0; j < width; j++)
    {
      Ray ray;
      camera.sampleRay(Vec2f(float(j) / width, float(i) / height), ray);

      image.setPixel(i, j, 0.5f * (ray.direction + 1.0f));
    }
  }

  image.writePPM("output_test_camera.ppm");

  return 0;
}
#endif

#ifdef OLD
#include "camera.h"
#include "image.h"

int main()
{
  const int width = 512;
  const int height = 512;

  Camera camera(Vec3f(0), Vec3f(0, 0, -1));

  Image image(width, height);
  for (int i = 0; i < height; ++i)
  {
    for (int j = 0; j < width; ++j)
    {
      Ray ray;
      float pdf;

      const float u = (2.0f * j - width) / width;
      const float v = (2.0f * i - height) / height;

      if (camera.sampleRay(Vec2f(u, v), ray, pdf))
      {
        image.setPixel(i, j, 0.5f * (ray.direction + 1.0f));
      }
      else
      {
        image.setPixel(i, j, Vec3f(0));
      }
    }
  }

  image.writePPM("output_test_camera_old.ppm");

  return 0;
}
#endif