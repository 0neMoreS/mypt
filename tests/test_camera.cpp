// #include "camera.h"
// #include "image.h"
// #include "core.h"
// #include "sampler.h"

// int main()
// {
//   const int w = 400, h = 400;
//   const int n_samples = 32;
//   Image image(w, h);
//   Camera camera(Vec3f(0, 0, 5), Vec3f(0, 0, 0), Vec3f(0, 1, 0), 45.0, float(w) / float(h), 1.0, 100.0);

// #pragma omp parallel for schedule(dynamic, 1)
//   for (unsigned short y = 0; y < h; y++)
//   {
//     fprintf(stderr, "\rRendering (%d spp) %5.2f%%", n_samples, 100. * y / (h - 1));
//     for (unsigned short x = 0; x < w; x++)
//     {
//       UniformSampler sampler((y * w + x) * 9781 + 1);
//       for (unsigned short s = 0; s < n_samples; s++)
//       {
//         float dx = sampler.getNext1D();
//         float dy = sampler.getNext1D();
//         Ray r;
//         camera.getRay(Vec2f(float(x) / w, float(y) / h), r);

//         image.setPixel(y, x, 0.5f * (r.direction + 1.0f));
//       }
//     }
//   }

//   image.writePPM("output.ppm");

//   return 0;
// }

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
      const float u = (2.0f * j - width) / height;
      const float v = (2.0f * i - height) / height;

      Ray ray;
      float pdf;
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

  image.writePPM("output.ppm");

  return 0;
}