#include "scene.hpp"

void Scene::naive_scene()
{
    shared_ptr<Material> ivory = make_shared<BlinnPhong>(Vec3f(0.1, 0.6, 0.4), Vec3f(0.5, 0.5, 0.4), 50.);
    shared_ptr<Material> red_rubber = make_shared<BlinnPhong>(Vec3f(0.1, 0.9, 0.1), Vec3f(0.8, 0.1, 0.1), 10.);
    shared_ptr<Material> mirror = make_shared<BlinnPhong>(Vec3f(0, 0.9, 0.9), Vec3f(0., 0., 0.), 1425., 0, 0.8, 0);
    shared_ptr<Material> glass = make_shared<BlinnPhong>(Vec3f(0, 0.7, 0.8), Vec3f(0., 0., 0.), 800., 1.5, 0.2, 0.8);

    load_envmap("../img/envmap.jpg");
    Add_light(make_shared<Light>(Vec3f(-20, 20, 20), 1.5));
    Add_light(make_shared<Light>(Vec3f(30, 50, -25), 1.8));

    std::vector<shared_ptr<Primitive>> primitives;
    primitives.push_back(make_shared<Sphere>(Vec3f(0, 0, -10), 1, ivory));
    primitives.push_back(make_shared<Sphere>(Vec3f(1.5, -0.5, -18), 3, red_rubber));
    primitives.push_back(make_shared<Sphere>(Vec3f(-1.0, -1.5, -12), 2, glass));
    primitives.push_back(make_shared<Sphere>(Vec3f(-2, -1, -4), 1, glass));
    primitives.push_back(make_shared<Sphere>(Vec3f(7, 5, -18), 4, mirror));
    primitives.push_back(make_shared<Sphere>(Vec3f(-7, 5, -18), 4, mirror));
    // primitives.push_back(make_shared<Plane>(Vec3f(0, -3, 0), Vec3f(0, -2, 0), ivory));

    add_model(Model("../models/duck.obj"), glass);
    // add_model(Model("../models/bunny/bunny_big.obj"), glass);
    Add(make_shared<BVH_Node>(primitives, 0, 1.0));
}

void Scene::path_tracing_scene()
{
    shared_ptr<Material> ivory = make_shared<Lambertian>(Vec3f(0.7, 0.7, 0.6));
    shared_ptr<Material> red_rubber = make_shared<Lambertian>(Vec3f(0.9, 0.1, 0.1));
    shared_ptr<Material> mirror = make_shared<Metal>(Vec3f(1, 1, 1), 0.0);
    shared_ptr<Material> glass = make_shared<Dielectric>(Vec3f(1, 1, 1), 1.5);
    shared_ptr<Material> checker = make_shared<Lambertian>(make_shared<CheckerTexture>(make_shared<SolidColor>(0.2, 0.3, 0.1), make_shared<SolidColor>(0.9, 0.9, 0.9)));
    shared_ptr<Material> noise = make_shared<Lambertian>(make_shared<NoiseTexture>(4));
    shared_ptr<Material> img = make_shared<Lambertian>(make_shared<ImageTexture>("../img/t2.jpg"));

    load_envmap("../img/envmap.jpg");

    std::vector<shared_ptr<Primitive>> primitives;
    primitives.push_back(make_shared<Sphere>(Vec3f(0, 0, -10), 1, ivory));
    primitives.push_back(make_shared<Sphere>(Vec3f(1.5, -0.5, -18), 3, red_rubber));
    primitives.push_back(make_shared<Sphere>(Vec3f(-1.0, -1.5, -12), 2, glass));
    // primitives.push_back(make_shared<Sphere>(Vec3f(-2, -1, -4), 1, glass));
    primitives.push_back(make_shared<Sphere>(Vec3f(7, 5, -18), 4, mirror));
    // primitives.push_back(make_shared<Sphere>(Vec3f(-7, 5, -18), 4, mirror));
    // primitives.push_back(make_shared<Plane>(Vec3f(0, -3, 0), Vec3f(0, -2, 0), ivory));

    primitives.push_back(make_shared<Sphere>(Vec3f(-7, 5, -18), 4, noise));
    primitives.push_back(make_shared<Sphere>(Vec3f(-2, -1, -4), 1, img));
    primitives.push_back(make_shared<MovingSphere>(Vec3f(-4, 0, -10), Vec3f(-4, -0.3, -10), 0, 1, 0.3, ivory));

    add_model(Model("../models/duck.obj"), glass);
    // add_model(Model("../models/bunny/bunny_big.obj"), glass, Vec3f(-400, -100, -400));
    Add(make_shared<BVH_Node>(primitives, 0, 1.0));
}

void Scene::random_scene()
{
    std::vector<shared_ptr<Primitive>> spheres;

    auto ground_material = std::make_shared<Lambertian>(Color(0.5, 0.5, 0.5));
    spheres.push_back(std::make_shared<Sphere>(Point3(0, -1002, 0), 1000, ground_material));
    auto light = std::make_shared<DiffuseLight>(Color(15, 15, 15));

    for (int a = -11; a < 11; ++a)
    {
        for (int b = -24; b < -0; ++b)
        {
            auto choose_mat = random_float();
            Point3 center(a + 0.9 * random_float(), -1.8, b + 0.9 * random_float());

            if ((center - Point3(4, 0.2, 0)).norm() > 0.9)
            {
                shared_ptr<Material> sphere_material;

                if (choose_mat < 0.8)
                {
                    // diffuse
                    auto albedo = wiseProduct(random_vector3(), random_vector3());
                    sphere_material = std::make_shared<Lambertian>(albedo);
                    spheres.push_back(std::make_shared<Sphere>(center, 0.2, sphere_material));
                }
                else if (choose_mat < 0.95)
                {
                    // metal
                    auto albedo = random_vector3(0.5, 1);
                    auto fuzz = random_float(0, 0.5);
                    sphere_material = std::make_shared<Metal>(albedo, fuzz);
                    spheres.push_back(std::make_shared<Sphere>(center, 0.2, sphere_material));
                }
                else
                {
                    // glass
                    sphere_material = std::make_shared<Dielectric>(1.5);
                    // spheres.push_back(std::make_shared<Sphere>(center, 0.2, sphere_material));
                    spheres.push_back(std::make_shared<Sphere>(center, 0.2, light));
                    Add_light(std::make_shared<Sphere>(center, 0.2, light));
                }
            }
        }
    }

    auto material1 = std::make_shared<Dielectric>(1.5);
    spheres.push_back(std::make_shared<Sphere>(Point3(0, -1, -5), 1.0, material1));

    auto material2 = std::make_shared<Lambertian>(Color(0.4, 0.2, 0.1));
    spheres.push_back(std::make_shared<Sphere>(Point3(-4, -1, -5), 1.0, material2));

    auto material3 = std::make_shared<Metal>(Color(0.7, 0.6, 0.5), 0.0);
    spheres.push_back(std::make_shared<Sphere>(Point3(4, -1, -5), 1.0, light));
    Add_light(std::make_shared<Sphere>(Point3(4, -1, -5), 1.0, light));

    Add(make_shared<BVH_Node>(spheres, 0, 1));
}

void Scene::final_scene()
{
    std::vector<shared_ptr<Primitive>> box1;
    auto ground = make_shared<Lambertian>(Color(0.48, 0.83, 0.53));

    const int boxes_per_side = 20;
    for (int j = 0; j < boxes_per_side; j++)
    {
        for (int i = 0; i < boxes_per_side; i++)
        {
            auto w = 3.0;
            auto x0 = -30.0 + i * w;
            auto z0 = -63.0 + j * w;
            auto y0 = -8.0;
            auto x1 = x0 + w;
            auto y1 = random_float(-7, -5);
            auto z1 = z0 + w;

            box1.push_back(make_shared<Box>(Point3(x0, y0, z0), Point3(x1, y1, z1), ground));
        }
    }

    std::vector<shared_ptr<Primitive>> primitives;

    primitives.push_back(make_shared<BVH_Node>(box1, 0, 1));

    auto light = make_shared<DiffuseLight>(Color(15, 15, 15));
    // primitives.push_back(make_shared<Plane>(Vec3f(-8, 28, -30), Vec3f(-8, 28, -45), Vec3f(8, 28, -45), Vec3f(8, 28, -30), light));
    Add(make_shared<Plane>(Vec3f(-10, 25, -10), Vec3f(-10, 25, -36), Vec3f(10, 25, -36), Vec3f(10, 25, -10), light));       // light
    Add_light(make_shared<Plane>(Vec3f(-10, 25, -10), Vec3f(-10, 25, -36), Vec3f(10, 25, -36), Vec3f(10, 25, -10), light)); // light

    auto center1 = Point3(-8, 0, -10);
    auto center2 = center1 + Vec3f(0.5, 0, 0);
    auto moving_sphere_material = make_shared<Lambertian>(Color(0.7, 0.3, 0.1));
    primitives.push_back(make_shared<MovingSphere>(center1, center2, 0.0, 1.0, 0.8, moving_sphere_material));

    primitives.push_back(make_shared<Sphere>(Point3(-1, 0, -5), 1, make_shared<Dielectric>(1.5)));
    primitives.push_back(make_shared<Sphere>(Point3(-3, 2, -5), 0.8, make_shared<Metal>(Color(0.8, 0.8, 0.9), 1.0)));

    auto boundary = make_shared<Sphere>(Point3(0, -2, -5), 1, make_shared<Dielectric>(1.5));
    primitives.push_back(boundary);
    primitives.push_back(make_shared<ConstantMedium>(boundary, 0.2, Color(0.2, 0.4, 0.9)));
    boundary = make_shared<Sphere>(Point3(0, 0, 0), 500, make_shared<Dielectric>(1.5));
    primitives.push_back(make_shared<ConstantMedium>(boundary, 0.001, Color(1, 1, 1)));
    boundary = make_shared<Sphere>(Point3(0, -1, -4), 2, make_shared<Dielectric>(1.5));
    primitives.push_back(make_shared<ConstantMedium>(boundary, 0.5, Color(0.3, 0.1, 0.3)));
    boundary = make_shared<Sphere>(Point3(0, -1, -4), 1, make_shared<Dielectric>(1.5));
    primitives.push_back(make_shared<ConstantMedium>(boundary, 0.8, Color(1, 0.1, 0.1)));

    auto emat = make_shared<Lambertian>(make_shared<ImageTexture>("../img/earthmap.jpg"));
    primitives.push_back(make_shared<Sphere>(Point3(-3, -2, -5), 1.2, emat));
    auto pertext = make_shared<NoiseTexture>(4);
    primitives.push_back(make_shared<Sphere>(Point3(4, -3, -6), 1.2, make_shared<Lambertian>(pertext)));

    std::vector<shared_ptr<Primitive>> boxes2;
    auto white = make_shared<Lambertian>(Color(0.73, 0.73, 0.73));
    int ns = 1000;
    for (int j = 0; j < ns; j++)
    {
        boxes2.push_back(make_shared<Sphere>(random_vector3(0, 8), 0.5, white));
    }
    primitives.push_back(make_shared<Translate>(make_shared<RotateY>(make_shared<BVH_Node>(boxes2, 0.0, 1.0), 15), Vec3f(4, 0, -28)));

    Add(make_shared<BVH_Node>(primitives, 0, 1));
}

void Scene::cornell_box()
{
    // backgroundColor = Color(1, 1, 1);
    std::vector<shared_ptr<Primitive>> primitives;

    auto red = make_shared<Lambertian>(Color(0.65, 0.05, 0.05));
    auto white = make_shared<Lambertian>(Color(0.73, 0.73, 0.73));
    auto green = make_shared<Lambertian>(Color(0.12, 0.45, 0.15));
    auto light = make_shared<DiffuseLight>(Color(15, 15, 15));

    Add(make_shared<Plane>(Vec3f(-30, -20, 0), Vec3f(-30, -20, -50), Vec3f(-30, 30, -50), Vec3f(-30, 30, 0), green));   // left
    Add(make_shared<Plane>(Vec3f(30, -20, 0), Vec3f(30, -20, -50), Vec3f(30, 30, -50), Vec3f(30, 30, 0), red));         // right
    Add(make_shared<Plane>(Vec3f(-30, -20, 0), Vec3f(30, -20, 0), Vec3f(30, -20, -50), Vec3f(-30, -20, -50), white));   // bottom
    Add(make_shared<Plane>(Vec3f(-30, 30, 0), Vec3f(30, 30, 0), Vec3f(30, 30, -50), Vec3f(-30, 30, -50), white));       // top
    Add(make_shared<Plane>(Vec3f(-30, -20, -50), Vec3f(30, -20, -50), Vec3f(30, 30, -50), Vec3f(-30, 30, -50), white)); // back
    Add(make_shared<Plane>(Vec3f(-5, 30, -30), Vec3f(-5, 30, -40), Vec3f(5, 30, -40), Vec3f(5, 30, -30), light));       // light
    Add_light(make_shared<Plane>(Vec3f(-5, 30, -30), Vec3f(-5, 30, -40), Vec3f(5, 30, -40), Vec3f(5, 30, -30), light)); // light, just for sampling

    auto mirror = make_shared<Metal>(Color(0.8, 0.8, 0.8), 0.0);
    shared_ptr<Primitive> box1 = make_shared<Box>(Vec3f(-15, -20, -45), Vec3f(0, 10, -35), mirror);
    box1 = make_shared<RotateY>(box1, 15);
    box1 = make_shared<Translate>(box1, Vec3f(8, 0, 1));

    // shared_ptr<Primitive> box2 = make_shared<Box>(Vec3f(-2, -20, -32), Vec3f(12, -10, -20), white);
    // box2 = make_shared<RotateY>(box2, -18);
    // box2 = make_shared<Translate>(box2, Vec3f(-8, 0, -5));

    auto glass = make_shared<Dielectric>(1.5);
    shared_ptr<Primitive> box2 = make_shared<Sphere>(Point3(5, -15, -26), 5, glass);

    // Add(std::make_shared<BVH_Node>(primitives, 0.0, 1.0));
    Add(box1);
    Add(box2);
}

void Scene::cornell_smoke()
{
    std::vector<shared_ptr<Primitive>> primitives;

    auto red = make_shared<Lambertian>(Color(0.65, 0.05, 0.05));
    auto white = make_shared<Lambertian>(Color(0.73, 0.73, 0.73));
    auto green = make_shared<Lambertian>(Color(0.12, 0.45, 0.15));
    auto light = make_shared<DiffuseLight>(Color(15, 15, 15));

    Add(make_shared<Plane>(Vec3f(-30, -20, 0), Vec3f(-30, -20, -50), Vec3f(-30, 30, -50), Vec3f(-30, 30, 0), green));   // left
    Add(make_shared<Plane>(Vec3f(30, -20, 0), Vec3f(30, -20, -50), Vec3f(30, 30, -50), Vec3f(30, 30, 0), red));         // right
    Add(make_shared<Plane>(Vec3f(-30, -20, 0), Vec3f(30, -20, 0), Vec3f(30, -20, -50), Vec3f(-30, -20, -50), white));   // bottom
    Add(make_shared<Plane>(Vec3f(-30, 30, 0), Vec3f(30, 30, 0), Vec3f(30, 30, -50), Vec3f(-30, 30, -50), white));       // top
    Add(make_shared<Plane>(Vec3f(-30, -20, -50), Vec3f(30, -20, -50), Vec3f(30, 30, -50), Vec3f(-30, 30, -50), white)); // back
    Add(make_shared<Plane>(Vec3f(-5, 30, -30), Vec3f(-5, 30, -40), Vec3f(5, 30, -40), Vec3f(5, 30, -30), light));       // light
    Add_light(make_shared<Plane>(Vec3f(-5, 28, -30), Vec3f(-5, 28, -40), Vec3f(5, 28, -40), Vec3f(5, 28, -30), light)); // light

    // Add(std::make_shared<ConstantMedium>(std::make_shared<Box>(Vec3f(-15, -20, -45), Vec3f(0, 10, -35), white), 0.01, Color(0.5, 0.5, 0.5)));
    shared_ptr<Primitive> box1 = make_shared<Box>(Vec3f(-15, -20, -45), Vec3f(0, 10, -35), white);
    box1 = make_shared<RotateY>(box1, 15);
    box1 = make_shared<Translate>(box1, Vec3f(8, 0, 1));
    box1 = make_shared<ConstantMedium>(box1, 0.3, Color(0., 0., 0.));

    shared_ptr<Primitive> box2 = make_shared<Box>(Vec3f(-2, -20, -32), Vec3f(13, -10, -22), white);
    box2 = make_shared<RotateY>(box2, -18);
    box2 = make_shared<Translate>(box2, Vec3f(-8, 0, -5));
    box2 = make_shared<ConstantMedium>(box2, 0.3, Color(0.93, 0.93, 0.93));

    Add(box1);
    Add(box2);
}

bool Scene::scene_intersect(const Ray &ray, float t_min, float t_max, Intersection &intersection) const
{
    // traverse all the primitives and find the closest one that intersects with the ray (should accelerate with BVH)
    bool hit_anything = false;
    auto closest_so_far = t_max;
    for (const auto &primitive : primitives) // primitive may be BVH_Node
    {
        if (primitive->ray_intersect(ray, t_min, closest_so_far, intersection))
        {
            hit_anything = true;
            closest_so_far = intersection.t;
        }
    }
    return hit_anything;
}

bool Scene::bounding_box(float time0, float time1, AABB &output_box) const
{
    if (primitives.empty())
        return false;

    AABB temp_box;
    bool first_box = true;

    for (const auto &object : primitives)
    {
        if (!object->bounding_box(time0, time1, temp_box))
            return false;
        output_box = first_box ? temp_box : surrounding_box(output_box, temp_box);
        first_box = false;
    }

    return true;
}

void Scene::add_model(Model &&model, shared_ptr<Material> material, Vec3f translate, float rotate_y)
{
    std::vector<shared_ptr<Primitive>> model_primitives;
    for (int i = 0; i < model.nfaces(); ++i)
    {
        model_primitives.push_back(make_shared<Triangle>(model.point(model.vert(i, 0)),
                                                         model.point(model.vert(i, 1)),
                                                         model.point(model.vert(i, 2)),
                                                         material));
    }
    // Add(make_shared<BVH_Node>(model_primitives, 0, 1));
    Add(make_shared<Translate>(make_shared<RotateY>(make_shared<BVH_Node>(model_primitives, 0, 1), rotate_y), translate));
    // for (const auto &primitive : model_primitives)
    // {
    //     Add(primitive);
    // }
}

void Scene::load_envmap(const char *filename)
{
    int n = -1;
    unsigned char *pixmap = stbi_load(filename, &envmap_width, &envmap_height, &n, 0);
    if (!pixmap || n != 3)
    {
        std::cerr << "ERROR: could not load environment map" << filename << "\n";
        exit(1);
    }
    envmap.resize(envmap_width * envmap_height);
    for (int j = envmap_height - 1; j >= 0; j--)
    {
        for (int i = 0; i < envmap_width; i++)
        {
            int pixel_index = (j * envmap_width + i) * n;
            envmap[i + j * envmap_width] = Vec3f(powf(pixmap[pixel_index] / 255.0, 2.2),
                                                 powf(pixmap[pixel_index + 1] / 255.0, 2.2),
                                                 powf(pixmap[pixel_index + 2] / 255.0, 2.2)); // gamma correction
        }
    }
    stbi_image_free(pixmap);
}

Vec3f Scene::background(const Vec3f &dir) const
{
    if (envmap.empty())
    {
        return backgroundColor;
        auto t = 0.5 * (dir.y + 1.0);
        return (1.0 - t) * Color(1.0, 1.0, 1.0) + t * Color(0.5, 0.7, 1.0);
    }
    float phi = atan2(dir.z, dir.x);                // [-pi, pi]
    float theta = asin(dir.y);                      // [-pi/2, pi/2]
    int x = (phi + pi) / (2 * pi) * envmap_width;   // [0, width]
    int y = (-theta + pi / 2) / pi * envmap_height; // [0, height]
    return envmap[y * envmap_width + x];
}