#pragma once

#include "utils.hpp"
#include "primitive.hpp"

class BVH_Node : public Primitive
{
    // The children pointers are to generic hittables. They can be other bvh_nodes, or spheres, or any other hittable.
    std::shared_ptr<Primitive> left;
    std::shared_ptr<Primitive> right;
    AABB box;

public:
    BVH_Node() {}
    BVH_Node(std::vector<std::shared_ptr<Primitive>> &primitives, float time0, float time1) : BVH_Node(primitives, 0, primitives.size(), time0, time1) {}
    BVH_Node(std::vector<std::shared_ptr<Primitive>> &primitives, size_t start, size_t end, float time0, float time1);

    virtual bool ray_intersect(const Ray &ray, float t_min, float t_max, Intersection &intersection) const override;
    virtual bool bounding_box(float time0, float time1, AABB &output_box) const override;

};

inline bool box_compare(const std::shared_ptr<Primitive> &a, const std::shared_ptr<Primitive> &b, int axis)
{
    AABB box_a, box_b;
    if (!a->bounding_box(0, 0, box_a) || !b->bounding_box(0, 0, box_b))
        std::cerr << "No bounding box in BVH_Node constructor.\n";
    return box_a.min()[axis] < box_b.min()[axis];
}

inline bool box_x_compare(const std::shared_ptr<Primitive> &a, const std::shared_ptr<Primitive> &b)
{
    return box_compare(a, b, 0);
}

inline bool box_y_compare(const std::shared_ptr<Primitive> &a, const std::shared_ptr<Primitive> &b)
{
    return box_compare(a, b, 1);
}

inline bool box_z_compare(const std::shared_ptr<Primitive> &a, const std::shared_ptr<Primitive> &b)
{
    return box_compare(a, b, 2);
}


class Box : public Primitive
{
    Point3 box_min, box_max;
    // shared_ptr<Primitive> bvh_sides;
    std::vector<shared_ptr<Primitive>> sides;

public:
    Box() {}
    Box(const Point3 &p0, const Point3 &p1, shared_ptr<Material> m) : box_min(p0), box_max(p1)
    {
        // std::vector<shared_ptr<Primitive>> sides;
        sides.emplace_back(make_shared<Plane>(Point3(p0.x, p0.y, p0.z), Point3(p1.x, p0.y, p0.z), Point3(p1.x, p1.y, p0.z), Point3(p0.x, p1.y, p0.z), m)); // plane at z = p0.z
        sides.emplace_back(make_shared<Plane>(Point3(p0.x, p0.y, p1.z), Point3(p1.x, p0.y, p1.z), Point3(p1.x, p1.y, p1.z), Point3(p0.x, p1.y, p1.z), m)); // plane at z = p1.z
        sides.emplace_back(make_shared<Plane>(Point3(p0.x, p0.y, p0.z), Point3(p0.x, p0.y, p1.z), Point3(p0.x, p1.y, p1.z), Point3(p0.x, p1.y, p0.z), m)); // plane at x = p0.x
        sides.emplace_back(make_shared<Plane>(Point3(p1.x, p0.y, p0.z), Point3(p1.x, p0.y, p1.z), Point3(p1.x, p1.y, p1.z), Point3(p1.x, p1.y, p0.z), m)); // plane at x = p1.x
        sides.emplace_back(make_shared<Plane>(Point3(p0.x, p0.y, p0.z), Point3(p1.x, p0.y, p0.z), Point3(p1.x, p0.y, p1.z), Point3(p0.x, p0.y, p1.z), m)); // plane at y = p0.y
        sides.emplace_back(make_shared<Plane>(Point3(p0.x, p1.y, p0.z), Point3(p1.x, p1.y, p0.z), Point3(p1.x, p1.y, p1.z), Point3(p0.x, p1.y, p1.z), m)); // plane at y = p1.y
        // bvh_sides = make_shared<BVH_Node>(sides, 0, 1);
    }

    virtual bool ray_intersect(const Ray &ray, float t_min, float t_max, Intersection &intersection) const override
    {
        // return bvh_sides->ray_intersect(ray, t_min, t_max, intersection);
        // not elegant, repeated code
        bool hit_anything = false;
        float closest_so_far = t_max;
        for (const auto &side : sides)
        {
            if (side->ray_intersect(ray, t_min, closest_so_far, intersection))
            {
                hit_anything = true;
                closest_so_far = intersection.t;
            }
        }
        return hit_anything;
    }

    virtual bool bounding_box(float time0, float time1, AABB &output_box) const override
    {
        output_box = AABB(box_min, box_max);
        return true;
    }

};
