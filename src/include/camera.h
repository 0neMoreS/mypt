#pragma once
#include "core.h"

#ifndef OLD
class Camera
{
public:
    Vec3f position;
    Vec3f lower_left_corner;
    Vec3f horizontal;
    Vec3f vertical;
    Vec3f right, up, forward;

    Camera(Vec3f lookfrom, Vec3f lookat, Vec3f vup, double vfov, double aspect)
    {
        double theta = vfov * M_PI / 180;
        double half_height = tan(theta / 2);
        double half_width = aspect * half_height;

        position = lookfrom;
        forward = normalize(lookat - lookfrom);
        right = normalize(cross(forward, vup));
        up = normalize(cross(right, forward));

        lower_left_corner = position - right * half_width - up * half_height + forward;
        horizontal = right * 2 * half_width;
        vertical = up * 2 * half_height;

        std::cout << std::fixed << std::setprecision(3);

        std::cout << "[Camera] position: ("
            << position[0] << ", " << position[1] << ", " << position[2] << ")\n";

        std::cout << "[Camera] forward: ("
            << forward[0] << ", " << forward[1] << ", " << forward[2] << ")\n";

        std::cout << "[Camera] right: ("
            << right[0] << ", " << right[1] << ", " << right[2] << ")\n";

        std::cout << "[Camera] up: ("
            << up[0] << ", " << up[1] << ", " << up[2] << ")\n";

        std::cout << "[Camera] vfov: " << vfov << " degrees\n";

        std::cout << "[Camera] half_height: " << half_height << "\n";

        std::cout << "[Camera] half_width: " << half_width << "\n";

        std::cout << "[Camera] lower_left_corner: ("
            << lower_left_corner[0] << ", " << lower_left_corner[1] << ", " << lower_left_corner[2] << ")\n";
    }

    void sampleRay(const Vec2f& uv, Ray& ray, float& pdf) const
    {
        ray = Ray(position, normalize(lower_left_corner + uv[0] * horizontal + uv[1] * vertical - position));
        pdf = 1.0f;
    }
};
#endif

#ifdef OLD
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
    Camera(const Vec3f& position, const Vec3f& forward, float FOV = 0.5f * PI)
        : position(position), forward(forward)
    {
        right = normalize(cross(forward, Vec3f(0, 1, 0)));
        up = normalize(cross(right, forward));

        // compute focal length from FOV
        focal_length = 1.0f / std::tan(0.5f * FOV);

        std::cout << std::fixed << std::setprecision(3);
        std::cout << "[Camera] position: ("
            << position[0] << ", " << position[1] << ", " << position[2] << ")\n";
        std::cout << "[Camera] forward: ("
            << forward[0] << ", " << forward[1] << ", " << forward[2] << ")\n";
        std::cout << "[Camera] right: ("
            << right[0] << ", " << right[1] << ", " << right[2] << ")\n";
        std::cout << "[Camera] up: ("
            << up[0] << ", " << up[1] << ", " << up[2] << ")\n";
        std::cout << "[Camera] FOV: " << FOV * 180.0f / PI << " degrees\n";
    }

    // sample ray emitting from the given sensor coordinate
    // NOTE: uv: [-1, -1] x [1, 1], sensor coordinate
    bool sampleRay(const Vec2f& uv, Ray& ray, float& pdf) const
    {
        const Vec3f pinholePos = position + focal_length * forward;
        const Vec3f sensorPos = position + uv[0] * right + uv[1] * up;
        ray = Ray(sensorPos, normalize(pinholePos - sensorPos));
        pdf = 1.0f;
        return true;
    }
};
#endif