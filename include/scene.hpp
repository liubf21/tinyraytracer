#pragma once

#include "utils.hpp"
#include "bvh.hpp"
#include "model.hpp"
#include "volume.hpp"

enum RenderAlgo
{
    RAY_TRACING,
    PATH_TRACING,
    PHONTON_MAPPING
};

class Scene
{
    // creating the scene (adding objects and lights)
    std::vector<shared_ptr<Primitive>> primitives;
    std::vector<shared_ptr<Primitive>> lights;
    
    Vec3f backgroundColor = Vec3f(0.0f, 0.0f, 0.0f);
    int envmap_width, envmap_height;
    std::vector<Vec3f> envmap;

public:
    // setting up options
    int width;
    int height;
    float aspect_ratio;
    float fov;
    int samples; // number of samples per pixel
    int maxDepth;
    int iterations = 60;
    float radius = 2.0f;

    RenderAlgo algo = PATH_TRACING;

    Scene(int w = 600, int h = 600, float f = 90, int s = 20, int d = 100) : width(w), height(h), fov(f), samples(s), maxDepth(d)
    {
        aspect_ratio = (float)width / (float)height;
    }
    Scene(std::vector<shared_ptr<Primitive>> &primitives, std::vector<shared_ptr<Primitive>> &lights) : primitives(primitives), lights(lights) {}

    void clear_scene()
    {
        primitives.clear();
        lights.clear();
        envmap.clear();
    }

    void Add(shared_ptr<Primitive> object) { primitives.push_back(std::move(object)); }
    void Add_light(shared_ptr<Primitive> light) { lights.push_back(std::move(light)); }

    const std::vector<shared_ptr<Primitive>> &get_primitives() const { return primitives; }
    const std::vector<shared_ptr<Primitive>> &get_lights() const { return lights; }

    void naive_scene();
    void path_tracing_scene();
    void random_scene();
    void final_scene();
    void cornell_box();
    void cornell_smoke();
    bool scene_intersect(const Ray &ray, float t_min, float t_max, Intersection &intersection) const;
    bool bounding_box(float time0, float time1, AABB &output_box) const;
    void add_model(Model &&model, shared_ptr<Material> material, Vec3f translate = Vec3f(), float rotate_y = .0);
    void load_envmap(const char *filename);
    Vec3f background(const Vec3f &dir) const;
};
