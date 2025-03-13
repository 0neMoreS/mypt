#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2
#define SPHERES 9
#define MAXDEPTH 5
#define LAMBERTALBEDO 0.7
// #define ORIGINCODE

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
    bool operator==(const Vec &b) const { return x == b.x && y == b.y && z == b.z; }
    Vec mult(const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }
    Vec &norm() { return *this = *this * (1 / sqrt(x * x + y * y + z * z)); }
    double dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; }
    Vec operator%(Vec &b) { return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); } // cross:
};
struct Ray
{
    Vec o, d;
    Ray(Vec o_, Vec d_) : o(o_), d(d_.norm()) {}
};
enum Refl_t
{
    DIFFUSE,
    SPECULAR,
    REFRACTIVE,
    LIGHT
}; // material types, used in radiance()
struct Sphere
{
    double rad;  // radius
    Vec p, e, c; // position, emission, color
    Refl_t refl; // reflection type (DIFFUSEuse, SPECULARular, REFRACTIVEactive)
    Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_) : rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}
    double intersect(const Ray &r) const
    { // returns distance(in root from, not real distance), 0 if not hit
        Vec op = r.o - p;
        double a = r.d.dot(r.d);
        double b = 2 * r.d.dot(op);
        double c = op.dot(op) - rad * rad;
        double delta = b * b - 4 * a * c;
        if (delta < 0)
        {
            return 0;
        }
        double root1 = (-b - sqrt(delta)) / 2 * a;
        double root2 = (-b + sqrt(delta)) / 2 * a;
        return root1 > 0 ? root1 : (root2 > 0 ? root2 : 0); // consider the light source located in the sphere
    }
};
Sphere spheres[SPHERES] = {
    // Scene: radius, position, emission, color, material
    // [-50, 50] ^ 3
    Sphere(1e5, Vec(1e5 - 50, 0, 0), Vec(), Vec(.75, .25, .25), DIFFUSE),   // Left
    Sphere(1e5, Vec(-1e5 + 50, 0, 0), Vec(), Vec(.25, .25, .75), DIFFUSE),  // Rght
    Sphere(1e5, Vec(0, 0, 1e5 - 50), Vec(), Vec(.25, .75, .75), DIFFUSE),   // Back
    Sphere(1e5, Vec(0, 0, -1e5 + 50), Vec(), Vec(.9, .2, .5), DIFFUSE),     // Frnt
    Sphere(1e5, Vec(0, 1e5 - 50, 0), Vec(), Vec(.75, .25, .75), DIFFUSE),   // Botm
    Sphere(1e5, Vec(0, -1e5 + 50, 0), Vec(), Vec(.75, .75, .75), DIFFUSE),  // Top
    Sphere(16.5, Vec(-20, -18, -14), Vec(), Vec(1, 1, 1) * .999, SPECULAR), // Mirr
    Sphere(16.5, Vec(13, 8, 8), Vec(), Vec(1, 1, 1) * .999, SPECULAR),      // Glas
    Sphere(1e3, Vec(0, 1e3 + 50 - 0.5, 15), Vec(12, 12, 12), Vec(), LIGHT)  // Lite
};

inline double clamp(double x) { return x < 0 ? 0 : x > 1 ? 1
                                                         : x; } // limit x to [0,1]

inline int toInt(double x) { return int(pow(clamp(x), 1 / 2.2) * 255 + .5); } // gamma

#ifdef ORIGINCODE
inline bool intersect(const Ray &r, double &t, int &id)
{
    double n = sizeof(spheres) / sizeof(Sphere), d, inf = t = 1e20;
    for (int i = int(n); i--;)
        if ((d = spheres[i].intersect(r)) && d < t)
        {
            t = d;
            id = i;
        }
    return t < inf;
}

Vec radiance(const Ray &r, int depth, unsigned short *Xi)
{
    double t;   // distance to intersection
    int id = 0; // id of intersected object
    if (!intersect(r, t, id))
        return Vec();                // if miss, return black
    const Sphere &obj = spheres[id]; // the hit object
    Vec x = r.o + r.d * t, n = (x - obj.p).norm(), nl = n.dot(r.d) < 0 ? n : n * -1, f = obj.c;
    double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y
                                                        : f.z; // max refl
    if (++depth > 5)
        if (erand48(Xi) < p)
            f = f * (1 / p);
        else
            return obj.e; // R.R.
    if (obj.refl == DIFFUSE)
    { // Ideal DIFFUSEUSE reflection
        double r1 = 2 * M_PI * erand48(Xi), r2 = erand48(Xi), r2s = sqrt(r2);
        Vec w = nl, u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm(), v = w % u;
        Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
        return obj.e + f.mult(radiance(Ray(x, d), depth, Xi));
    }
    else if (obj.refl == SPECULAR) // Ideal SPECULARULAR reflection
        return obj.e + f.mult(radiance(Ray(x, r.d - n * 2 * n.dot(r.d)), depth, Xi));
    Ray reflRay(x, r.d - n * 2 * n.dot(r.d)); // Ideal dielectric REFRACTIVEACTION
    bool into = n.dot(nl) > 0;                // Ray from outside going in?
    double nc = 1, nt = 1.5, nnt = into ? nc / nt : nt / nc, ddn = r.d.dot(nl), cos2t;
    if ((cos2t = 1 - nnt * nnt * (1 - ddn * ddn)) < 0) // Total internal reflection
        return obj.e + f.mult(radiance(reflRay, depth, Xi));
    Vec tdir = (r.d * nnt - n * ((into ? 1 : -1) * (ddn * nnt + sqrt(cos2t)))).norm();
    double a = nt - nc, b = nt + nc, R0 = a * a / (b * b), c = 1 - (into ? -ddn : tdir.dot(n));
    double Re = R0 + (1 - R0) * c * c * c * c * c, Tr = 1 - Re, P = .25 + .5 * Re, RP = Re / P, TP = Tr / (1 - P);
    return obj.e + f.mult(depth > 2 ? (erand48(Xi) < P ? // Russian roulette
                                           radiance(reflRay, depth, Xi) * RP
                                                       : radiance(Ray(x, tdir), depth, Xi) * TP)
                                    : radiance(reflRay, depth, Xi) * Re + radiance(Ray(x, tdir), depth, Xi) * Tr);
}
#endif

#ifndef ORIGINCODE
Vec radiance(const Ray &r, int depth, unsigned short *Xi)
{
    // ========================== intersect ==========================
    int id = -1;
    int min_t = __INT_MAX__;
    for (int i = 0; i < SPHERES; i++)
    {
        if (spheres[i].intersect(r) != 0 && spheres[i].intersect(r) < min_t)
        {
            min_t = spheres[i].intersect(r);
            id = i;
        }
    }
    if (id == -1)
    {
        return Vec();
    }
    Sphere &obj = spheres[id];
    // ========================== init ==========================
    // try no R.R. first
    if (++depth >= MAXDEPTH)
    {
        return obj.e;
    }
    Vec f = obj.c; // rgb can standard radiance
    Vec x = r.o + r.d * min_t;
    Vec n = (x - obj.p).norm();
    Vec nl = n.dot(r.d) < 0 ? n : n * -1;
    // ========================== specular ==========================
    if (obj.refl == SPECULAR)
    {
        return obj.e + f.mult(radiance(Ray{x, r.d - nl * 2 * n.dot(r.d)}, depth, Xi)) * abs(r.d.dot(nl));
    }
    // ========================== diffuse ==========================
    else if (obj.refl == DIFFUSE)
    { // Ideal DIFFUSEUSE reflection
        Vec random_dir{2 * erand48(Xi) - 1, 2 * erand48(Xi) - 1, 2 * erand48(Xi) - 1};
        random_dir.norm();
        Vec output_dir = x + nl + random_dir;
        return obj.e + f.mult(radiance(Ray{x, output_dir - x}, depth, Xi)) * LAMBERTALBEDO;
    }
    // ========================== light ==========================
    else if (obj.refl == LIGHT)
    {
        return obj.e;
    }
    return obj.e;
}
#endif

int main(int argc, char *argv[])
{
    int w = 800, h = 400, samps = argc == 2 ? atoi(argv[1]) : 30;
    FILE *f = fopen("image.ppm", "w"); // Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    Vec *c = new Vec[w * h];
    Vec cam{0.0, 0.0, 50.0};
    Vec lower_left_corner{2.0, 1.0, 49.0};
    Vec horizontal{-4.0, 0.0, 0.0};
    Vec vertical{0.0, -2.0, 0.0};
    Vec r;
#pragma omp parallel for schedule(dynamic, 1) private(r)
    for (unsigned short y = 0; y < h; y++)
    {
        fprintf(stderr, "\rRendering (%d spp) %5.2f%%", samps, 100. * y / (h - 1));
        for (unsigned short x = 0; x < w; x++)
        {
            unsigned short Xi[3] = {x, y, (unsigned short)(x * y)};
            r = Vec();
            for (unsigned short s = 0; s < samps; s++)
            {
                double dx = erand48(Xi);
                double dy = erand48(Xi);
                Vec pixel = horizontal * ((double)(x + dx) / (double)w) + vertical * ((double)(y + dy) / (double)h) + lower_left_corner;
                Vec ro = cam;
                Vec rd = pixel - cam;
                // printf("ro: %f, %f, %f \nrd: %f, %f, %f \n", ro.x, ro.y, ro.z, rd.x, rd.y, rd.z);
                r = r + radiance(Ray{ro, rd}, 0, Xi) * (1.0 / samps);
            }
            int i = y * w + (w - x);
            c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z));
        }
    }

    for (int i = 0; i < w * h; i++)
        fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
}