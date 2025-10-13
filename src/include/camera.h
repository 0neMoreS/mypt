#pragma once
#include "core.h"

struct Camera
{
    Vec3f origin;
    Vec3f lower_left_corner;
    Vec3f horizontal;
    Vec3f vertical;
    Vec3f u, v, w;

    Camera(
        Vec3f lookfrom, Vec3f lookat, Vec3f vup,
        double vfov,
        double aspect,
        double near, double far)
    {
        double theta = vfov * M_PI / 180;
        double half_height = tan(theta / 2) * (near - far);
        double half_width = aspect * half_height;

        origin = lookfrom;
        w = normalize(lookfrom - lookat);
        u = normalize(cross(vup, w));
        v = cross(w, u);

        lower_left_corner = origin - u * half_width - v * half_height - w;
        horizontal = u * 2 * half_width;
        vertical = v * 2 * half_height;
    }

    void getRay(const Vec2f &uv, Ray &ray)
    {
        ray = Ray(origin, normalize(lower_left_corner + uv[0] * horizontal + uv[1] * vertical - origin));
    }
};