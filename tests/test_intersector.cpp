#include "camera.h"
#include "image.h"
#include "scene.h"

#ifndef OLD
int main()
{
  const int width = 512, height = 512;
  const int n_samples = 32;
  Image image(width, height);

  Scene scene;
  scene.loadModel("../models/CornellBox-Original.obj");
  scene.build();

  Camera camera(Vec3f(0, 1.f, 7.f), Vec3f(0, 1.f, -1.f), Vec3f(0, 1.f, 0), 25.f, float(width) / float(height));

#pragma omp parallel for schedule(dynamic, 1)
  for (unsigned short i = 0; i < height; i++)
  {
    fprintf(stderr, "\rRendering (%d spp) %5.2f%%", n_samples, 100. * i / (height - 1));
    for (unsigned short j = 0; j < width; j++)
    {
      Ray ray;
      float pdf;
      camera.sampleRay(Vec2f(float(j) / width, float(i) / height), ray, pdf);
      IntersectInfo info;
      if (scene.intersect(ray, info))
        image.setPixel((height - i - 1), (width - j - 1), 0.5f * (info.surfaceInfo.shadingNormal + 1.0f));
      else
        image.setPixel((height - i - 1), (width - j - 1), Vec3f(0));
    }
  }

  image.writePPM("output_test_intersector.ppm");
  return 0;
}
#endif

#ifdef OLD
int main()
{
  const int width = 512;
  const int height = 512;

  Scene scene;
  scene.loadModel("../models/CornellBox-Original.obj");
  scene.build();

  Camera camera(Vec3f(0, 1, 7), Vec3f(0, 0, -1), 0.25f * PI);

  Image image(width, height);
  for (int i = 0; i < height; ++i)
  {
    for (int j = 0; j < width; ++j)
    {
      const float u = (2.0f * j - width) / width;
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

  image.writePPM("output_test_intersector_old.ppm");

  return 0;
}
#endif