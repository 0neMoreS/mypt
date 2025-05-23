#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2
#define SPHERES 8
#define MAXDEPTH 5
#define LAMBERTALBEDO 0.5
#define SCALE 1
#define NAIR 1.0
#define NGLASS 1.8

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
    Vec operator/(double b) const { return Vec(x / b, y / b, z / b); }
    bool operator==(const Vec &b) const { return x == b.x && y == b.y && z == b.z; }
    Vec mult(const Vec &b) const { return Vec(x * b.x, y * b.y, z * b.z); }
    double moudle() { return sqrt(x * x + y * y + z * z); }
    Vec &norm() { return *this = *this * (1 / moudle()); }
    double dot(const Vec &b) const { return x * b.x + y * b.y + z * b.z; }
    Vec operator%(const Vec &b) const { return Vec(y * b.z - z * b.y, z * b.x - x * b.z, x * b.y - y * b.x); } // cross:
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
    Refl_t refl; // reflection type (DIFFUSEuse, SPECULARular, REFRACTIONactive)
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
    Sphere(1e5, Vec(1e5 - 50, 0, 0), Vec(), Vec(.75, .25, .25), DIFFUSE),  // Left
    Sphere(1e5, Vec(-1e5 + 50, 0, 0), Vec(), Vec(.25, .25, .75), DIFFUSE), // Rght
    Sphere(1e5, Vec(0, 0, 1e5 - 50), Vec(), Vec(.25, .75, .25), DIFFUSE),  // Back
    // Sphere(1e5, Vec(0, 0, -1e5 + 50), Vec(), Vec(.75, .25, .75), DIFFUSE), // Frnt
    // Sphere(1e5, Vec(0, 0, -1e5 + 50), Vec(), Vec(.0, .0, .0), DIFFUSE),    // Frnt
    Sphere(1e5, Vec(0, 1e5 - 50, 0), Vec(), Vec(.25, .75, .75), DIFFUSE),  // Botm
    Sphere(1e5, Vec(0, -1e5 + 50, 0), Vec(), Vec(.75, .75, .25), DIFFUSE), // Top
    Sphere(16.5, Vec(-20, -18, -14), Vec(), Vec(1, 1, 1), SPECULAR),       // Mirr
    Sphere(16.5, Vec(13, 8, 8), Vec(), Vec(1, 1, 1), REFRACTIVE),          // Glas
    Sphere(1e3, Vec(0, 1e3 + 50 - 0.5, 15), Vec(16, 16, 16), Vec(), LIGHT) // Lite
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
    {
        if (depth > 16)
        {
            return obj.e;
        }
        if (erand48(Xi) < p)
            f = f * (1 / p);
        else
            return obj.e; // R.R.
    }
    if (obj.refl == DIFFUSE)
    { // Ideal DIFFUSEUSE reflection
        double r1 = 2 * M_PI * erand48(Xi), r2 = erand48(Xi), r2s = sqrt(r2);
        Vec w = nl, u = ((fabs(w.x) > .1 ? Vec(0, 1) : Vec(1)) % w).norm(), v = w % u;
        Vec d = (u * cos(r1) * r2s + v * sin(r1) * r2s + w * sqrt(1 - r2)).norm();
        return obj.e + f.mult(radiance(Ray(x, d), depth, Xi));
    }
    else if (obj.refl == SPECULAR) // Ideal SPECULARULAR reflection
        return obj.e + f.mult(radiance(Ray(x, r.d - n * 2 * n.dot(r.d)), depth, Xi));
    Ray reflRay(x, r.d - n * 2 * n.dot(r.d)); // Ideal dielectric REFRACTIONACTION
    bool into = n.dot(r.d) < 0;               // Ray from outside going in?
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
    double min_t = 1e20;
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
    Vec f = obj.c; // rgb can standard radiance
    Vec x = r.o + r.d * min_t;
    Vec n = (x - obj.p).norm();
    Vec nl = n.dot(r.d) < 0 ? n : n * -1;
    // f = f * -(nl.dot(r.d)); // cos in rendering equation
    double p = f.x > f.y && f.x > f.z ? f.x : f.y > f.z ? f.y
                                                        : f.z;
    if (++depth >= MAXDEPTH)
    {
        if (depth > 8)
        {
            return obj.e;
        }
        if (erand48(Xi) < p)
        {

            f = f / p;
        }
        else
        {
            return obj.e; // R.R.
        }
    }
    // ========================== specular ==========================
    if (obj.refl == SPECULAR)
    {
        return obj.e + f.mult(radiance(Ray{x, r.d - nl * 2 * nl.dot(r.d)}, depth, Xi));
    }
    // ========================== diffuse ==========================
    else if (obj.refl == DIFFUSE)
    { // Ideal DIFFUSEUSE reflection
        double random_theta = 2 * M_PI * erand48(Xi);
        double random_phi = acos(2 * erand48(Xi) - 1);
        double random_r = pow(erand48(Xi), 1.0 / 3);
        Vec random_dir{random_r * sin(random_theta) * cos(random_phi), random_r * sin(random_theta) * sin(random_phi), random_r * cos(random_theta)};
        Vec output_dir = nl + random_dir;
        // return obj.e + f.mult(radiance(Ray{x, output_dir}, depth, Xi)) * LAMBERTALBEDO / (3 * random_r * random_r * sin(random_phi) / 4 * M_PI);
        return obj.e + f.mult(radiance(Ray{x, output_dir}, depth, Xi)) * LAMBERTALBEDO * -(nl.dot(r.d));
    }
    // ========================== light ==========================
    else if (obj.refl == LIGHT)
    {
        return obj.e;
    }
    // ========================== refraction ==========================
    else if (obj.refl == REFRACTIVE)
    {
        bool air_to_glass = n.dot(r.d) < 0; // N = nl
        double n1 = air_to_glass ? NAIR : NGLASS;
        double n2 = air_to_glass ? NGLASS : NAIR;
        Vec reflect_output_dir = r.d - n * 2 * n.dot(r.d);
        double cos_theta1 = nl.dot(r.d);
        double discriminant = 1.0 - ((n1 * n1) / (n2 * n2)) * (1 - cos_theta1 * cos_theta1);
        if (discriminant < 0)
        {
            return obj.e + f.mult(radiance(Ray{x, reflect_output_dir}, depth, Xi));
        }
        // Vec refract_output_dir = (L - nl * cos_theta1) * n1 / n2 - nl * sqrt(discriminant);
        Vec refract_output_dir = r.d * (n1 / n2) - nl * cos_theta1 * (n1 / n2) - nl * sqrt(discriminant);
        double R0 = ((n1 - n2) * (n1 - n2)) / ((n1 + n2) * (n1 + n2));
        // double c = 1 - (air_to_glass ? -cos_theta1 : refract_output_dir.dot(n));
        double Re = R0 + (1 - R0) * pow(1 + cos_theta1, 5.0);
        if (depth < 2)
        {
            return obj.e + f.mult(radiance(Ray{x, refract_output_dir}, depth, Xi) * (1 - Re) + radiance(Ray{x, reflect_output_dir}, depth, Xi) * Re);
        }
        else
        {
            double PR = .25 + .5 * Re;
            if (erand48(Xi) < PR)
            {

                return obj.e + f.mult(radiance(Ray{x, reflect_output_dir}, depth, Xi) * Re) / PR;
            }
            else
            {
                return obj.e + f.mult(radiance(Ray{x, refract_output_dir}, depth, Xi) * (1 - Re)) / (1 - PR);
            }
        }
    }
    return Vec{};
}
#endif

int main(int argc, char *argv[])
{
    const int w = 400 * SCALE, h = 400 * SCALE, samps = argc == 2 ? atoi(argv[1]) : 30;
    FILE *f = fopen("image.ppm", "w"); // Write image to PPM file.
    fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
    Vec *c = new Vec[w * h];
    const double near = 50.0;
    const double far = 10.0;
    // const double far = 49.0;
    const double vfov = 90 * M_PI / 180;
    const double w_ratio_h = 4 / 4;
    const double half_height = tan(vfov / 2) * (near - far);
    const double half_width = half_height * w_ratio_h;
    const Vec vup{0.0, 1.0, 0.0};
    const Vec look_at{10.0, -10.0, 0.0};
    const Vec cam{0.0, 0.0, near + 10.0};
    const Vec cam_z = (cam - look_at).norm();
    const Vec cam_x = (vup % cam_z).norm();
    const Vec cam_y = cam_z % cam_x;
    const Vec lower_left_corner = cam - cam_x * half_width - cam_y * half_height - cam_z;
    const Vec horizontal = cam_x * 2 * half_width;
    const Vec vertical = cam_y * 2 * half_height;
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
                // translate pixel pos on film plane to world pos
                Vec pixel = horizontal * ((double)(x + dx) / (double)w) + vertical * ((double)(y + dy) / (double)h) + lower_left_corner;
                Vec ro = cam;
                Vec rd = pixel - cam;
                // printf("ro: %f, %f, %f \nrd: %f, %f, %f \n", ro.x, ro.y, ro.z, rd.x, rd.y, rd.z);
                r = r + radiance(Ray{pixel, cam_z * -1.0}, 0, Xi) / samps;
                // r = r + radiance(Ray{ro, rd}, 0, Xi) / samps;
            }
            int i = (h - y - 1) * w + x;
            c[i] = c[i] + Vec(clamp(r.x), clamp(r.y), clamp(r.z));
        }
    }
    for (int i = 0; i < w * h; i++)
        fprintf(f, "%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
}