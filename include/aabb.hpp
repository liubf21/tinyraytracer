#pragma once

#include "utils.hpp"

class AABB
{
    Point3 minimum;
    Point3 maximum;

public:
    AABB() {}
    AABB(const Point3 &a, const Point3 &b)
    {
        minimum = a;
        maximum = b;
    }

    Point3 min() const { return minimum; }
    Point3 max() const { return maximum; }

    bool ray_intersect(const Ray &ray, float t_min, float t_max) const
    {
        for (int a = 0; a < 3; a++)
        {
            auto invD = 1.0f / ray.direction()[a];
            auto t0 = (min()[a] - ray.origin()[a]) * invD;
            auto t1 = (max()[a] - ray.origin()[a]) * invD;
            if (invD < 0.0f)
                std::swap(t0, t1);
            t_min = t0 > t_min ? t0 : t_min;
            t_max = t1 < t_max ? t1 : t_max;
            if (t_max <= t_min)
                return false;
        }
        return true;

        // Vec3f hit1 = minimum - ray.origin();
        // Vec3f hit2 = maximum - ray.origin();
        // float tenter = std::numeric_limits<float>::min(), texit = std::numeric_limits<float>::max();
        // for (size_t i = 0; i < 3; i++)
        // {
        //     if (ray.direction()[i] != 0)
        //     {
        //         float t1 = hit1[i] / ray.direction()[i];
        //         float t2 = hit2[i] / ray.direction()[i];
        //         if (t1 > t2)
        //             std::swap(t1, t2); // make sure t1 < t2, direction is important
        //         tenter = std::max(tenter, t1);
        //         texit = std::min(texit, t2);
        //     }
        // }
        // return tenter < texit && texit > 0;
    }
};

inline AABB surrounding_box(AABB box0, AABB box1) // return a box that contains both box0 and box1
{
    Point3 small(fmin(box0.min().x, box1.min().x),
                 fmin(box0.min().y, box1.min().y),
                 fmin(box0.min().z, box1.min().z));

    Point3 big(fmax(box0.max().x, box1.max().x),
               fmax(box0.max().y, box1.max().y),
               fmax(box0.max().z, box1.max().z));

    return AABB(small, big);
}
