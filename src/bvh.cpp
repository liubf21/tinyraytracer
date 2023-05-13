#include "bvh.hpp"

BVH_Node::BVH_Node(std::vector<std::shared_ptr<Primitive>> &primitives, size_t start, size_t end, float time0, float time1)
{
    int axis = random_int(0, 2);
    // used to sort the primitives
    auto comparator = (axis == 0) ? box_x_compare
                    : (axis == 1) ? box_y_compare
                    : box_z_compare;

    size_t object_span = end - start;

    if (object_span == 1)
    {
        left = right = primitives[start];
    }
    else if (object_span == 2)
    {
        if (comparator(primitives[start], primitives[start + 1]))
        {
            left = primitives[start];
            right = primitives[start + 1];
        }
        else
        {
            left = primitives[start + 1];
            right = primitives[start];
        }
    }
    else
    {
        std::sort(primitives.begin() + start, primitives.begin() + end, comparator);

        auto mid = start + object_span / 2;
        left = std::make_shared<BVH_Node>(primitives, start, mid, time0, time1); // recursive call
        right = std::make_shared<BVH_Node>(primitives, mid, end, time0, time1);
    }

    AABB box_left, box_right;

    if (!left->bounding_box(time0, time1, box_left) || !right->bounding_box(time0, time1, box_right))
        std::cerr << "No bounding box in BVH_Node constructor.\n";

    box = surrounding_box(box_left, box_right);
}

bool BVH_Node::ray_intersect(const Ray &ray, float t_min, float t_max, Intersection &intersection) const
{
    if (!box.ray_intersect(ray, t_min, t_max))
        return false;

    bool hit_left = left->ray_intersect(ray, t_min, t_max, intersection); // recursive call
    bool hit_right = right->ray_intersect(ray, t_min, hit_left ? intersection.t : t_max, intersection);

    return hit_left || hit_right;
}

bool BVH_Node::bounding_box(float time0, float time1, AABB &output_box) const
{
    output_box = box;
    return true;
}
