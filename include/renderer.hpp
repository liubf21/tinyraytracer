#pragma once

#include <fstream>

#include "camera.hpp"
#include "scene.hpp"
#include "kdtree.hpp"


class Renderer
{
    std::vector<Vec3f> framebuffer;

public:
    void render(const Camera &camera, const Scene &scene, RenderAlgo algo = PATH_TRACING);
    Vec3f ray_tracing(const Ray &ray, const Scene &scene, int depth = 0);
    Vec3f path_tracing(const Ray &ray, const Scene &scene, int depth = 0);
    std::vector<PPMnode> pm_backtrace(const Ray &ray, const Scene &scene, int depth, int index, Vec3f pref = Vec3f(1, 1, 1), float prob = 1.0);
    void pm_forward(const Ray &ray, const Scene &scene, int depth, Vec3f col, IMGbuf *c, KDTree *kdt, float prob = 1.0);

    void write_to_ppm(int width, int height, float gamma = 2.2, std::string filename = std::string());
    void write_to_jpg(int width, int height, float gamma = 2.2, std::string filename = std::string());
};
