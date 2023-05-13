#pragma once

#include "ray.hpp"
#include "utils.hpp"

class Camera
{
private:
    Point3 origin;
    Point3 top_left_corner;
    Vec3f horizontal;
    Vec3f vertical;
    Vec3f u, v, w; // u, v, w are the basis vectors of the camera
    float lens_radius; // the radius of the lens
    float time0, time1; // shutter open/close times

public:
    Camera(
        Point3 lookfrom,
        Point3 lookat,
        Vec3f vup, // the up vector of the camera
        float vfov, // vertical field-of-view in degrees
        float aspect_ratio, // the ratio of the width to the height of the viewport
        float aperture, // the size of the lens
        float focus_dist, // the distance between the projection point and the plane where everything is in perfect focus
        float _time0 = 0.0f, float _time1 = 0.0f
    )
    {
        auto theta = degrees_to_radians(vfov);
        auto h = tan(theta / 2);
        // std::cout << "theta: " << theta << " h: " << h << std::endl;
        auto viewport_height = 2.0 * h;
        auto viewport_width = aspect_ratio * viewport_height;

        w = (lookfrom - lookat).normalize();
        u = cross(vup, w).normalize();
        v = cross(w, u);

        origin = lookfrom;
        horizontal = focus_dist * viewport_width * u; // the horizontal vector of the plane
        vertical = focus_dist * viewport_height * v; // the vertical vector of the plane
        top_left_corner = origin - horizontal / 2 + vertical / 2 - focus_dist * w; // the lower left corner of the plane

        lens_radius = aperture / 2;
        time0 = _time0;
        time1 = _time1;
    }

    Ray get_ray(float s, float t) const
    {
        auto rd = lens_radius * random_in_unit_disk(); // simulate the lens
        auto offset = u * rd[0] + v * rd[1];
        // return the ray from the origin to the point on the viewport
        // std::cout<<"cam";
        return Ray(
            origin + offset,
            top_left_corner + s * horizontal - t * vertical - origin - offset,
            random_float(time0, time1)
        );
    }
};

