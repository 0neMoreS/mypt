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

  Camera camera(Vec3f(0, 1.f, 7.f), Vec3f(0, 1.f, -1.f), Vec3f(0, 1, 0), 45.f, float(width) / float(height), 0.1f, 100.f);

  for (int i = 0; i < height; ++i)
  {
    for (int j = 0; j < width; ++j)
    {
      const float u = (2.0f * j - width) / height;
      const float v = (2.0f * i - height) / height;

      Ray ray;
      camera.sampleRay(Vec2f(u, v), ray);
      IntersectInfo info;
      if (scene.intersect(ray, info))
      {
        image.setPixel(i, j, 0.5f * (info.surfaceInfo.shadingNormal + 1.0f));
      }

      image.setPixel(i, j, Vec3f(0));
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

  image.writePPM("output_test_intersector_old.ppm");

  return 0;
}
#endif