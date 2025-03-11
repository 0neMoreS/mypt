#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2

// g++ -O3 -fopenmp mypt.cpp -o mypt

struct Vec
{                   // Usage: time ./smallpt 5000 && xv image.ppm
    double x, y, z; // position, also color (r,g,b)
    Vec(double x_ = 0, double y_ = 0, double z_ = 0)
    {
        x = x_;
        y = y_;
        z = z_;
    }
    Vec operator+(const Vec &b) const { return Vec(x + b.x, y + b.y, z + b.z); }
    Vec operator-(const Vec &b) const { return Vec(x - b.x, y - b.y, z - b.z); }
    Vec operator*(double b) const { return Vec(x * b, y * b, z * b); }
    Vec mult(const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }
    Vec &norm() { return *this = *this * (1 / sqrt(x * x + y * y + z * z)); }
    double dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; }
    Vec operator%(Vec &b) { return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); } // cross:
};
struct Ray
{
    Vec o, d;
    Ray(Vec o_, Vec d_) : o(o_), d(d_) {}
};
enum Refl_t
{
    DIFF,
    SPEC,
    REFR
}; // material types, used in radiance()
struct Sphere
{
    double rad;  // radius
    Vec p, e, c; // position, emission, color
    Refl_t refl; // reflection type (DIFFuse, SPECular, REFRactive)
    Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_) : rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}
    double intersect(const Ray &r) const
    { // returns distance, 0 if nohit TODO
        double a = r.d.dot(r.d);
        double b = 2 * (r.d.dot(r.o) - r.d.dot(p));
        double c = r.o.dot(r.o) + p.dot(p) - rad * rad;
        double delta = b * b - 4 * a * c;
        if (delta < 0)
        {
            return 0;
        }
        return (-b - sqrt(delta)) / 2 * a;
        // Vec op = p - r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
        // double t, eps = 1e-4, b = op.dot(r.d), det = b * b - op.dot(op) + rad * rad;
        // if (det < 0)
        //     return 0;
        // else
        //     det = sqrt(det);
        // return (t = b - det) > eps ? t : ((t = b + det) > eps ? t : 0);
    }
};
Sphere spheres[] = {
    // Scene: radius, position, emission, color, material
    Sphere(1e5, Vec(1e5 + 1, 40.8, 81.6), Vec(), Vec(.75, .25, .25), DIFF),   // Left
    Sphere(1e5, Vec(-1e5 + 99, 40.8, 81.6), Vec(), Vec(.25, .25, .75), DIFF), // Rght
    Sphere(1e5, Vec(50, 40.8, 1e5), Vec(), Vec(.75, .75, .75), DIFF),         // Back
    Sphere(1e5, Vec(50, 40.8, -1e5 + 170), Vec(), Vec(), DIFF),               // Frnt
    Sphere(1e5, Vec(50, 1e5, 81.6), Vec(), Vec(.75, .75, .75), DIFF),         // Botm
    Sphere(1e5, Vec(50, -1e5 + 81.6, 81.6), Vec(), Vec(.75, .75, .75), DIFF), // Top
    Sphere(16.5, Vec(27, 16.5, 47), Vec(), Vec(1, 1, 1) * .999, SPEC),        // Mirr
    Sphere(16.5, Vec(73, 16.5, 78), Vec(), Vec(1, 1, 1) * .999, REFR),        // Glas
    Sphere(600, Vec(50, 681.6 - .27, 81.6), Vec(12, 12, 12), Vec(), DIFF)     // Lite
};

inline double clamp(double x) { return x < 0 ? 0 : x > 1 ? 1
                                                         : x; }

inline int toInt(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); }

Vec radiance(const Ray &r, int depth, unsigned short *Xi)
{
    int min_dis = __INT_MAX__;
    Vec color{0, 0, 0};
    for (int i = 0; i < 9; i++)
    {
        if (spheres[i].intersect(r) != 0 && spheres[i].intersect(r) < min_dis)
        {
            min_dis = spheres[i].intersect(r);
            color = spheres[i].c;
        }
    }
    return color;
}

// int main(int argc, char *argv[])
// {
//     int w = 1024 / 8, h = 768 / 8, samps = argc == 2 ? atoi(argv[1]) / 4 : 30;
//     Ray cam(Vec(50, 52, 295.6), Vec()); // cam pos, dir
//     FILE *f = fopen("image.ppm", "w");  // Write image to PPM file.
//     fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
//     Vec *c = new Vec[w * h];

//     for (unsigned short x = 0; x < w; x++)
//     {
//         for (unsigned short y = 0; y < h; y++)
//         {
//             Vec r{0, 0, 0};
//             for (unsigned short s = 0; s < samps; s++)
//             {
//                 unsigned short Xi[3] = {0, 0, (unsigned short)(x * y)};
//                 double rx = erand48(Xi);
//                 double ry = erand48(Xi);
//                 Vec d = {x + rx, y + ry, -1.0};
//                 r = r + radiance(Ray{cam.o, d.norm()}, 2, Xi) * (1.0 / samps);
//             }
//         }
//     }

//     for (int i = 0; i < w * h; i++)
//         fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
// }

int main(int argc, char *argv[])
{
    int w = 1024 / 8, h = 768 / 8, samps = argc == 2 ? atoi(argv[1]) / 4 : 30;
    Ray cam(Vec(50, 52, 295.6), Vec(0, -0.042612, -1).norm()); // cam pos, dir
    Vec cx = Vec(w * .5135 / h), cy = (cx % cam.d).norm() * .5135, r, *c = new Vec[w * h];
    for (int y = 0; y < h; y++)
    { // Loop over image rows
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps * 4, 100. * y / (h - 1));
        for (unsigned short x = 0, Xi[3] = {0, 0, (unsigned short)(y * y * y)}; x < w; x++)
            for (int sy = 0, i = (h - y - 1) * w + x; sy < 2; sy++) // 2x2 subpixel rows
                for (int sx = 0; sx < 2; sx++, r = Vec())
                { // 2x2 subpixel cols
                    for (int s = 0; s < samps; s++)
                    {
                        double r1 = 2 * erand48(Xi), dx = r1 < 1 ? sqrt(r1) - 1 : 1 - sqrt(2 - r1);
                        double r2 = 2 * erand48(Xi), dy = r2 < 1 ? sqrt(r2) - 1 : 1 - sqrt(2 - r2);
                        Vec d = cx * (((sx + .5 + dx) / 2 + x) / w - .5) +
                                cy * (((sy + .5 + dy) / 2 + y) / h - .5) + cam.d;
                        r = r + radiance(Ray(cam.o + d * 140, d.norm()), 0, Xi) * (1. / samps);
                    } // Camera rays are pushed ^^^^^ forward to start in interior
                    c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z)) * .25;
                }
    }
    FILE *f = fopen("image.ppm", "w"); // Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    for (int i = 0; i < w * h; i++)
        fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
}