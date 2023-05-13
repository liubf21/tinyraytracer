#include "renderer.hpp"

void Renderer::render(const Camera &camera, const Scene &scene, RenderAlgo algo)
{
    framebuffer.resize(scene.width * scene.height);

    switch (scene.algo)
    {
    case RAY_TRACING:    // ray traing
#pragma omp parallel for // multi-threading
        for (int j = 0; j < scene.height; ++j)
        {
            std::cerr << "\rRendering " << 100. * j / scene.height << "%" << std::flush;
            for (int i = 0; i < scene.width; ++i)
            {
                // normalize the pixel coordinates to the range [-1, 1]
                // j = 0, i = 0 in the image is the bottom left corner
                auto u = (i + 0.5) / scene.width;
                auto v = (j + 0.5) / scene.height;
                Ray r = camera.get_ray(u, v);
                framebuffer[i + j * scene.width] = ray_tracing(r, scene);
            }
        }
        write_to_jpg(scene.width, scene.height);
        break;
    case PATH_TRACING:   // path tracing
#pragma omp parallel for // multi-threading
        for (int j = 0; j < scene.height; ++j)
        {
            std::cerr << "\rRendering (" << scene.samples << " spp)" << 100. * j / scene.height << "%" << std::flush;
            for (int i = 0; i < scene.width; ++i)
            {
                for (int s = 0; s < scene.samples; ++s)
                {
                    for (int _ = 0; _ < 4; ++_) // 4 samples per pixel (2x2 subpixel sampling)
                    {
                        auto u = (i + 0.5 * (s % 2) + random_float() * 0.5) / scene.width;
                        auto v = (j + 0.5 * (s / 2) + random_float() * 0.5) / scene.height;
                        Ray r = camera.get_ray(u, v);
                        framebuffer[i + j * scene.width] = framebuffer[i + j * scene.width] + path_tracing(r, scene);
                    }
                }
                framebuffer[i + j * scene.width] = framebuffer[i + j * scene.width] * (1. / scene.samples);
            }
        }
        write_to_jpg(scene.width, scene.height);
        break;
    case PHONTON_MAPPING: // photon mapping
        int nth = omp_get_num_procs();
        IMGbuf **c = new IMGbuf *[nth];
        for (int i = 0; i < nth; i++)
        {
            c[i] = new IMGbuf[scene.width * scene.height];
        }
        IMGbuf *final = new IMGbuf[scene.width * scene.height];
        IMGbuf *now = new IMGbuf[scene.width * scene.height];
        std::vector<PPMnode> ball[nth]; // used for multi-threading
        KDTree tree;
        int iter = scene.iterations;
        float samp = scene.samples;
        float rad = scene.radius, alpha = 0.8;
        for (int _ = 1; _ <= iter; fprintf(stderr, "\riter %d done!\n", _), ++_)
        {
            if (_ < 3 || _ % 10 == 0) // if not enough photons, or every 10 iterations
            {
                if (_ > 1)
                {
                    // samp /= sqrt(alpha);
                    samp /= alpha * alpha;
                    rad *= alpha;
                }
                framebuffer.resize(scene.width * scene.height);
#pragma omp parallel for num_threads(nth) schedule(dynamic, 1)
                for (int y = 0; y < scene.height; ++y)
                {
                    int num = omp_get_thread_num();
                    fprintf(stderr, "\rbuiding kd-tree %5.2f%%", 100. * y / scene.height);
                    for (int x = 0; x < scene.width; ++x)
                    {
                        for (int s = 0; s < 4; ++s) // 4 samples per pixel
                        {
                            auto u = (x + 0.5 * (s % 2) + random_float() * 0.5) / scene.width;
                            auto v = (y + 0.5 * (s / 2) + random_float() * 0.5) / scene.height;
                            Ray r = camera.get_ray(u, v);
                            std::vector<PPMnode> tmp = pm_backtrace(r, scene, 0, y * scene.width + x); // backtrace and store the result in tmp
                            for (auto &i : tmp)
                                if (i.index >= 0)
                                {
                                    i.r = rad;
                                    ball[num].push_back(i);
                                    // std::cout << i.col << std::endl;
                                    framebuffer[y * scene.width + x] = ball[num].back().col;
                                }
                        }
                    }
                }
                write_to_jpg(scene.width, scene.height, 1.0, "../output/pm");
                // return;

                std::vector<PPMnode> totball; // total ball
                fprintf(stderr, "\rbuild tree ...");
                for (auto &i : ball)
                    totball.insert(totball.end(), i.begin(), i.end());
                std::cout << "totball size = " << totball.size() << std::endl;
                tree.init(totball);
                fprintf(stderr, "done!\n");
            }
            fprintf(stderr, "rad = %f samp = %.0f\n", rad, samp);
            int per = samp / nth + 1;
#pragma omp parallel for num_threads(nth) schedule(dynamic, 1)
            for (int t = 0; t < nth; ++t)
            {
                int num = 0;
                for (int __ = 0; __ < per; __++)
                {
                    if (num == 0 && __ % 1000 == 0)
                        fprintf(stderr, "\rtracing %5.2f%%", 100. * __ / per);
                    auto vec = scene.get_lights()[0]->random();
                    Ray r = Ray(vec, random_in_unit_sphere()); // generate a ray from light source
                    Vec3f col = Vec3f(1, 1, 1) + .4;
                    tree.query(PPMnode(r.origin(), col, r.direction()), c[num]); // query the tree about the photons before sending the ray
                    pm_forward(r, scene, 0, col, c[num], &tree);
                    // std::cout << "c[num] " << c[num][0].getcol() << ":" << c[num][0].n << std::endl;
                }
            }
            memset(now, 0, sizeof(now));
            for (int i = 0; i < nth; memset(c[i++], 0, sizeof(c[i])))
            {
                for (int j = 0; j < scene.width * scene.height; ++j)
                {
                    now[j] = now[j] + c[i][j];
                    // if (c[i][j].n > 0)
                    // {
                    // std::cout << "c[i] " << c[i][j].getcol() << ":" << c[i][j].n << std::endl;
                    // std::cout << "now " << now[j].getcol() <<":"<< now[j].n << std::endl;
                    // }
                }
            }
            framebuffer.resize(scene.width * scene.height);
            for (int i = 0; i < scene.width * scene.height; ++i)
            {
                // std::cout << final[i].getcol() << std::endl;
                final[i] = final[i] + now[i] / now[i].n; // * alpha;
                framebuffer[i] = final[i].get();
                if (final[i].n > 0)
                {
                    // std::cout << "final " << final[i].getcol() << ":" << final[i].n << std::endl;
                    // std::cout << "get " << final[i].get() << std::endl;
                    // std::cout << "framebuffer " << framebuffer[i] << std::endl;
                }
            }
            if (_ == 1 || _ % 10 == 0)
            {
                write_to_jpg(scene.width, scene.height);
                std::cout << "write to jpg" << std::endl;
                // return;
            }
        }
        break;
    }
}

Vec3f Renderer::ray_tracing(const Ray &ray, const Scene &scene, int depth)
{
    Intersection intersection;
    if (depth > 4 || !scene.scene_intersect(ray, epsilon, max_float, intersection))
    {
        return scene.background(ray.direction());
    }
    Vec3f point = intersection.position, N = intersection.normal;
    // std::cout << "point: " << point << std::endl;
    shared_ptr<BlinnPhong> material = std::dynamic_pointer_cast<BlinnPhong>(intersection.material); // BlinnPhong material
    Vec3f reflect_color = Vec3f();
    if (material->getRelf() > epsilon)
    {
        Vec3f reflect_dir = reflect(ray.direction(), N).normalize();
        Vec3f reflect_orig = reflect_dir * N < 0 ? point - N * 1e-3 : point + N * 1e-3; // offset the original point to avoid occlusion by the object itself
        reflect_color = ray_tracing(Ray(reflect_orig, reflect_dir), scene, depth + 1);
    }
    Vec3f refract_color = Vec3f();
    if (material->getRefr() > epsilon)
    {
        Vec3f refract_dir = refract(ray.direction(), N, material->getRefractiveIndex()).normalize();
        Vec3f refract_orig = refract_dir * N < 0 ? point - N * 1e-3 : point + N * 1e-3; // offset the original point to avoid occlusion by the object itself
        refract_color = ray_tracing(Ray(refract_orig, refract_dir), scene, depth + 1);
    }
    float ambient_light_intensity = 0.2, diffuse_light_intensity = 0, specular_light_intensity = 0;

    // auto light = std::dynamic_pointer_cast<Light>(scene.get_lights()[0]); // only one light source
    for (const auto &l : scene.get_lights())
    {
        auto light = std::dynamic_pointer_cast<Light>(l);
        Vec3f light_dir = (light->getPosition() - point).normalize(); // from the hit point to the light, which is opposite to the direction of the light
        float light_distance = (light->getPosition() - point).norm();
        // shadow, there is only a small subtlety: I perturb the point by moving it in the direction of normal
        Vec3f shadow_orig = light_dir * N < 0 ? point - N * 1e-3 : point + N * 1e-3;
        // following code will prevent the shadow if the light is behind the hit point
        Intersection shadow_intersection;
        if (scene.scene_intersect(Ray(shadow_orig, light_dir), epsilon, max_float, shadow_intersection) && (shadow_intersection.position - shadow_orig).norm() < light_distance)
            continue;
        diffuse_light_intensity += light->getIntensity() * std::max(0.f, light_dir * N);
        // specular_light_intensity += powf(std::max(0.f, -reflect(-light_dir, N) * dir), material.getSpecularExponent()) * light->getIntensity();
        Vec3f h = (light_dir - ray.direction()).normalize();
        specular_light_intensity += powf(std::max(0.f, h * N), material->getSpecularExponent()) * light->getIntensity();
    }
    // return Vec3f(.8, .3, .3);
    return material->getDiffuseColor() * (ambient_light_intensity * material->getAlbedo()[0] + diffuse_light_intensity * material->getAlbedo()[1]) +
           Vec3f(1., 1., 1.) * specular_light_intensity * material->getAlbedo()[2] + reflect_color * material->getRelf() + refract_color * material->getRefr();
}

Color Renderer::path_tracing(const Ray &ray, const Scene &scene, int depth)
{
    Intersection intersection;
    if (depth > scene.maxDepth)
        return Color();
    if (!scene.scene_intersect(ray, epsilon, max_float, intersection))
    {
        return scene.background(ray.direction());
    }
    ScatterRecord srec;
    Color emission = intersection.material->emitted(intersection.u, intersection.v, intersection.position);
    if (!intersection.material->scatter(ray, intersection, srec))
        return emission;

    // emission = emission + wiseProduct(attenuation, path_tracing(scattered, scene, depth + 1)) * intersection.material->scattering_pdf(ray, intersection, scattered) / pdf;
    // return emission;
    if (srec.is_specular) // specular, no direct lighting
    {
        return wiseProduct(srec.attenuation, path_tracing(srec.scattered, scene, depth + 1));
    }
    // Russian roulette
    float p = std::max(srec.attenuation[0], std::max(srec.attenuation[1], srec.attenuation[2]));
    if (random_float() < p)
    {
        srec.attenuation = srec.attenuation / p;
    }
    else
    {
        return emission;
    }
    // direct lighting
    if (random_float() < 0.5 && !scene.get_lights().empty()) // randomly choose direct lighting or indirect lighting
    {
        auto light = scene.get_lights()[0];
        // auto light = scene.get_lights()[random_int(0, scene.get_lights().size() - 1)];
        srec.scattered = Ray(intersection.position, light->random(intersection.position), ray.time()); // sample a light
        srec.pdf = light->pdf_value(srec.scattered.origin(), srec.scattered.direction());
        if (srec.pdf < epsilon) // avoid division by zero
            return emission;
    }
    // hard-coded direct lighting
    /*
    auto on_light = Point3(random_float(-5, 5), 30, random_float(-40, -30));
    auto to_light = on_light - intersection.position;
    auto distance_squared = to_light.norm2();
    to_light.normalize();
    if (to_light * intersection.normal < 0) // light is behind the hit point
        return emission;
    auto light_area = 100.0;
    auto light_cosine = fabs(to_light.y);
    if (light_cosine < 0.0001) // avoid division by zero
        return emission;
    pdf = distance_squared / (light_cosine * light_area);
    scattered = Ray(intersection.position, to_light, ray.time());
    */

    return emission + wiseProduct(srec.attenuation, path_tracing(srec.scattered, scene, depth + 1)) *
                          intersection.material->scattering_pdf(ray, intersection, srec.scattered) / srec.pdf;
}

// trace a ray from camera, stop when hit diffuse surface, store information in PPMnode
std::vector<PPMnode> Renderer::pm_backtrace(const Ray &ray, const Scene &scene, int depth, int index, Vec3f pref, float prob)
{
    std::vector<PPMnode> result, tmp;
    if (pref[0] < epsilon && pref[1] < epsilon && pref[2] < epsilon || prob < epsilon)
        return result;
    Intersection intersection;
    if (!scene.scene_intersect(ray, epsilon, max_float, intersection))
    {
        return result;
    }
    Color emission = intersection.material->emitted(intersection.u, intersection.v, intersection.position);
    // std::cout << "emission: " << emission << std::endl;
    ScatterRecord srec;
    if (!intersection.material->scatter(ray, intersection, srec))
    {
        result.push_back(PPMnode(intersection.position, emission, intersection.normal, index, prob, 1.0));
        return result;
    }
    if (srec.is_specular) // specular, no direct lighting
    {
        tmp = pm_backtrace(srec.scattered, scene, depth + 1, index, wiseProduct(pref, srec.attenuation) + emission, prob);
        result.insert(result.end(), tmp.begin(), tmp.end());
    }
    else
    {
        result.push_back(PPMnode(intersection.position, wiseProduct(pref, srec.attenuation) + emission, intersection.normal, index, prob, 1.0));
        // tmp = pm_backtrace(srec.scattered, scene, depth + 1, index, wiseProduct(pref, srec.attenuation), prob);
        // result.insert(result.end(), tmp.begin(), tmp.end());
    }
    return result;
}

// cast photon from light source, when hit diffuse surface, query KDTree and store information in IMGbuf
void Renderer::pm_forward(const Ray &ray, const Scene &scene, int depth, Vec3f col, IMGbuf *c, KDTree *kdt, float prob)
{
    if (col[0] < epsilon && col[1] < epsilon && col[2] < epsilon || prob < epsilon)
        return;
    Intersection intersection;
    ScatterRecord srec;
    if (!scene.scene_intersect(ray, epsilon, max_float, intersection) || !intersection.material->scatter(ray, intersection, srec))
    {
        return;
    }
    Color emission = intersection.material->emitted(intersection.u, intersection.v, intersection.position); // photon
    if (srec.attenuation[0] < epsilon && srec.attenuation[1] < epsilon && srec.attenuation[2] < epsilon)
    {
        kdt->query(PPMnode(intersection.position, emission, intersection.normal), c);
        return;
    }
    if (depth > scene.maxDepth)
        return;
    if (srec.is_specular) // specular, no direct lighting
    {
        pm_forward(srec.scattered, scene, depth + 1, wiseProduct(col, srec.attenuation) + emission, c, kdt, prob);
    }
    else
    {
        kdt->query(PPMnode(intersection.position, wiseProduct(col, srec.attenuation) + emission, intersection.normal), c); // query the kd-tree
        pm_forward(srec.scattered, scene, depth + 1, wiseProduct(col, srec.attenuation) + emission, c, kdt, prob);
    }
}

void Renderer::write_to_ppm(int width, int height, float gamma, std::string filename)
{
    if (filename.empty())
    {
        filename = "../output/out";
        int i = 0;
        while (std::ifstream(filename + ".ppm").good()) // check if the file exists
        {
            filename = "../output/out" + std::to_string(++i);
        }
    }
    std::cout << "filename: " << filename << std::endl;
    std::ofstream ofs(filename + ".ppm"); // save the framebuffer to file, ppm is a simple image format
    ofs << "P6\n"
        << width << " " << height << "\n255\n";
    for (Vec3f &c : framebuffer)
    {
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max > 1)
            c = c * (1. / max);
        for (int j = 0; j < 3; ++j)
        {
            if (c[j] != c[j])
                c[j] = 0;                                                // remove NaN
            ofs << (char)(255.999 * clamp(powf(c[j], 1. / 2.), 0., 1.)); // type conversion from float to char (gamma correction)
        }
    }
}

void Renderer::write_to_jpg(int width, int height, float gamma, std::string filename)
{
    // save as jpg
    if (filename.empty())
    {
        filename = "../output/out";
        int i = 0;
        while (std::ifstream(filename + ".jpg").good()) // check if the file exists
        {
            filename = "../output/out" + std::to_string(++i);
        }
    }
    std::cout << "filename: " << filename << std::endl;
    std::vector<unsigned char> pixmap(height * width * 3);
    for (int i = 0; i < height * width; ++i)
    {
        Vec3f &c = framebuffer[i];
        float max = std::max(c[0], std::max(c[1], c[2]));
        if (max > 1)
            c = c * (1. / max);
        for (int j = 0; j < 3; j++)
        {
            if (c[j] != c[j])
                c[j] = 0; // remove NaN
            pixmap[i * 3 + j] = (unsigned char)(255.999 * clamp(powf(framebuffer[i][j], 1. / gamma), 0., 1.));
        }
    }
    stbi_write_jpg((filename + ".jpg").c_str(), width, height, 3, pixmap.data(), 100);
}
