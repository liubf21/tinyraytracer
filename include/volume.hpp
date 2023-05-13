#pragma once

#include "primitive.hpp"
#include "material.hpp"
#include "texture.hpp"

class ConstantMedium : public Primitive
{
    shared_ptr<Primitive> boundary; // the boundary of the medium
    float neg_inv_density; // negative inverse density
    shared_ptr<Material> phase_function; // the material of the medium, which is isotropic

public:
    ConstantMedium(shared_ptr<Primitive> b, float d, shared_ptr<Texture> a) : boundary(b), neg_inv_density(-1 / d), phase_function(make_shared<Isotropic>(a)) {} // d is in the range of [0, 1]
    ConstantMedium(shared_ptr<Primitive> b, float d, Color c) : boundary(b), neg_inv_density(-1 / d), phase_function(make_shared<Isotropic>(c)) {}

    virtual bool ray_intersect(const Ray &ray, float t_min, float t_max, Intersection &intersection) const override
    {
        bool db = (random_float() < 0.0001); // debug
        db = false;
        Intersection intersection1, intersection2; // the two intersections of the ray with the boundary

        // make sure the ray intersects the boundary
        if (!boundary->ray_intersect(ray, -infinity, infinity, intersection1)) // find the first intersection
            return false;

        if (!boundary->ray_intersect(ray, intersection1.t + 0.0001, infinity, intersection2)) // find the second intersection
            return false;

        if (db)
        {
            std::cerr << "\nt0 t1 " << intersection1.t << " " << intersection2.t << "\n";
        }

        if (intersection1.t < t_min)
            intersection1.t = t_min;
        if (intersection2.t > t_max)
            intersection2.t = t_max;

        if (intersection1.t >= intersection2.t) // if the first intersection is after the second intersection
            return false;

        if (intersection1.t < 0)
            intersection1.t = 0;

        const auto distance_inside_boundary = (intersection2.t - intersection1.t); // because ray.direction() is a unit vector
        const auto hit_distance = neg_inv_density * log(random_float()); // the distance the ray travels inside the medium, using the exponential distribution

        if (hit_distance > distance_inside_boundary)
            return false;

        intersection.t = intersection1.t + hit_distance;
        intersection.position = ray.at(intersection.t);

        if (db)
        {
            std::cerr << "hit_distance = " << hit_distance << "\n"
                      << "intersection.t = " << intersection.t << "\n"
                      << "intersection.position = " << intersection.position << "\n";
        }

        intersection.normal = Vec3f(1, 0, 0); // arbitrary, not used
        intersection.material = phase_function;

        return true;
    }

    virtual bool bounding_box(float t0, float t1, AABB &output_box) const override
    {
        return boundary->bounding_box(t0, t1, output_box);
    }
};
