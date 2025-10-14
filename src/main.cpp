#include <core.h>
#include <camera.h>
#include <image.h>
#include <sampler.h>

int main(int argc, char *argv[])
{
    const int w = 400, h = 400;
    const int n_samples = 32;
    Image image(w, h);
    // Camera camera(Vec3f(0, 0, 5), Vec3f(0, 0, 0), Vec3f(0, 1, 0), 45.0, float(w) / float(h), 1.0, 100.0);
    // Camera camera(Vec3f(0), Vec3f(0, 0, -1));

    return 0;
}