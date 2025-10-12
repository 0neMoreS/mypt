#pragma once
#include "core.h"
#include "ray.h"

struct Camera
{
    Vec origin;
    Vec lower_left_corner;
    Vec horizontal;
    Vec vertical;
    Vec u, v, w;
    double lens_radius;

    Camera(
        Vec lookfrom, Vec lookat, Vec vup,
        double vfov,
        double aspect,
        double near, double far,
        double aperture = 0.0)
    {
        lens_radius = aperture / 2;
        double theta = vfov * M_PI / 180;
        double half_height = tan(theta / 2) * (near - far);
        double half_width = aspect * half_height;

        origin = lookfrom;
        w = (lookfrom - lookat).norm();
        u = (vup % w).norm();
        v = w % u;

        lower_left_corner = origin - u * half_width - v * half_height - w;
        horizontal = u * 2 * half_width;
        vertical = v * 2 * half_height;
    }

    Ray getRay(double s, double t, unsigned short *Xi)
    {
        Vec rd = Vec();
        Vec offset = Vec();
        return Ray(origin + offset,
                   (lower_left_corner + horizontal * s + vertical * t - origin - offset).norm());
    }
};