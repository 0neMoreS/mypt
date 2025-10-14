#pragma once
#include "core.h"

// class Camera
// {
// public:
//     Vec3f origin;
//     Vec3f lower_left_corner;
//     Vec3f horizontal;
//     Vec3f vertical;
//     Vec3f u, v, w;

//     Camera(
//         Vec3f lookfrom, Vec3f lookat, Vec3f vup,
//         double vfov,
//         double aspect,
//         double near, double far)
//     {
//         double theta = vfov * M_PI / 180;
//         double half_height = tan(theta / 2) * (near - far);
//         double half_width = aspect * half_height;

//         origin = lookfrom;
//         w = normalize(lookat - lookfrom);
//         u = normalize(cross(w, vup));
//         v = cross(u, w);

//         lower_left_corner = origin - u * half_width - v * half_height - w;
//         horizontal = u * 2 * half_width;
//         vertical = v * 2 * half_height;

//         std::cout << std::fixed << std::setprecision(3);

//         std::cout << "[Camera] position: ("
//                   << origin[0] << ", " << origin[1] << ", " << origin[2] << ")\n";

//         std::cout << "[Camera] forward: ("
//                   << w[0] << ", " << w[1] << ", " << w[2] << ")\n";

//         std::cout << "[Camera] right: ("
//                   << u[0] << ", " << u[1] << ", " << u[2] << ")\n";

//         std::cout << "[Camera] up: ("
//                   << v[0] << ", " << v[1] << ", " << v[2] << ")\n";
//     }

//     void getRay(const Vec2f &uv, Ray &ray)
//     {
//         ray = Ray(origin, normalize(lower_left_corner + uv[0] * horizontal + uv[1] * vertical - origin));
//     }
// };

class Camera
{
private:
    Vec3f position;
    Vec3f forward;
    Vec3f right;
    Vec3f up;

    float FOV;
    float focal_length;

public:
    Camera(const Vec3f &position, const Vec3f &forward, float FOV = 0.5f * PI)
        : position(position), forward(forward)
    {
        right = normalize(cross(forward, Vec3f(0, 1, 0)));
        up = normalize(cross(right, forward));

        // compute focal length from FOV
        focal_length = 1.0f / std::tan(0.5f * FOV);
    }

    // sample ray emitting from the given sensor coordinate
    // NOTE: uv: [-1, -1] x [1, 1], sensor coordinate
    bool sampleRay(const Vec2f &uv, Ray &ray, float &pdf) const
    {
        const Vec3f pinholePos = position + focal_length * forward;
        const Vec3f sensorPos = position + uv[0] * right + uv[1] * up;
        ray = Ray(sensorPos, normalize(pinholePos - sensorPos));
        pdf = 1.0f;
        return true;
    }
};