#include <fstream>
#include <vector>

#include "utils.hpp"
#include "camera.hpp"
#include "scene.hpp"
#include "renderer.hpp"

int main(int argc, char **argv)
{
    if (argc < 2)
    {
        std::cout << "use -h or --help to get help" << std::endl;
        return 0;
    }
    Scene scene;
    if (strcmp(argv[1], "-f") == 0)
    {
        scene.load_volume(argv[2]);
        // scene.naive_scene();
        scene.algo = RAY_TRACING;
    }
    else
        for (int i = 1; i < argc; ++i)
        {
            // algorithm, scene, width, height, samples
            if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)
            {
                std::cout << "usage: ./tinyraytracer -r <render algorithm> -sc <scene>(optional) -w <width>(optional) -h <height>(optional) -s <samples>(optional)" << std::endl;
                std::cout << "render algorithm: 1. naive 2. path tracing 3. progressive photon mapping" << std::endl;
                std::cout << "scene: 1. cornell box 2. cornell smoke 3. random scene 4. final scene" << std::endl;
                std::cout << "(naive algorithm can not specify scene)" << std::endl;
                std::cout << "if using progressive photon mapping, you can specify the number of iterations by adding -i <iterations> and the radius by adding -ra <radius>" << std::endl;
                return 0;
            }
            else if (strcmp(argv[i], "-r") == 0)
            {
                if (i + 1 < argc)
                {
                    if (strcmp(argv[i + 1], "1") == 0)
                    {
                        std::cout << "naive render" << std::endl;
                        scene.algo = RAY_TRACING;
                        scene.naive_scene();
                    }
                    else if (strcmp(argv[i + 1], "2") == 0)
                    {
                        std::cout << "path tracing" << std::endl;
                        scene.algo = PATH_TRACING;
                        scene.path_tracing_scene();
                    }
                    else if (strcmp(argv[i + 1], "3") == 0)
                    {
                        std::cout << "progressive photon mapping" << std::endl;
                        scene.algo = PHONTON_MAPPING;
                        scene.path_tracing_scene();
                        scene.samples = 8000;
                    }
                    else
                    {
                        std::cout << "invalid render algorithm" << std::endl;
                        return -1;
                    }
                }
                else
                {
                    std::cout << "invalid render algorithm" << std::endl;
                    return -1;
                }
            }
            else if (strcmp(argv[i], "-sc") == 0)
            {
                if (scene.algo == RAY_TRACING)
                {
                    std::cout << "have not implemented yet" << std::endl;
                    return -1;
                }
                if (i + 1 < argc)
                {
                    scene.clear_scene();
                    if (strcmp(argv[i + 1], "1") == 0)
                    {
                        std::cout << "cornell box" << std::endl;
                        scene.cornell_box();
                    }
                    else if (strcmp(argv[i + 1], "2") == 0)
                    {
                        std::cout << "cornell smoke" << std::endl;
                        scene.cornell_smoke();
                    }
                    else if (strcmp(argv[i + 1], "3") == 0)
                    {
                        std::cout << "random scene" << std::endl;
                        scene.random_scene();
                    }
                    else if (strcmp(argv[i + 1], "4") == 0)
                    {
                        std::cout << "final scene" << std::endl;
                        scene.final_scene();
                    }
                    else
                    {
                        std::cout << "invalid scene" << std::endl;
                        return -1;
                    }
                }
                else
                {
                    std::cout << "invalid scene" << std::endl;
                    return -1;
                }
            }
            else if (strcmp(argv[i], "-w") == 0)
            {
                if (i + 1 < argc)
                {
                    scene.width = atoi(argv[i + 1]);
                }
                else
                {
                    std::cout << "invalid width" << std::endl;
                    return -1;
                }
            }
            else if (strcmp(argv[i], "-h") == 0)
            {
                if (i + 1 < argc)
                {
                    scene.height = atoi(argv[i + 1]);
                }
                else
                {
                    std::cout << "invalid height" << std::endl;
                    return -1;
                }
            }
            else if (strcmp(argv[i], "-s") == 0)
            {
                if (i + 1 < argc)
                {
                    scene.samples = atoi(argv[i + 1]);
                }
                else
                {
                    std::cout << "invalid samples" << std::endl;
                    return -1;
                }
            }
            else if (strcmp(argv[i], "-i") == 0)
            {
                if (scene.algo != PHONTON_MAPPING)
                {
                    std::cout << "invalid argument" << std::endl;
                    return -1;
                }
                if (i + 1 < argc)
                {
                    scene.iterations = atoi(argv[i + 1]);
                }
                else
                {
                    std::cout << "invalid iterations" << std::endl;
                    return -1;
                }
            }
            else if (strcmp(argv[i], "-ra") == 0)
            {
                if (scene.algo != PHONTON_MAPPING)
                {
                    std::cout << "invalid argument" << std::endl;
                    return -1;
                }
                if (i + 1 < argc)
                {
                    scene.radius = atof(argv[i + 1]);
                }
                else
                {
                    std::cout << "invalid radius" << std::endl;
                    return -1;
                }
            }
        }

    Camera camera(Vec3f(0, 0, 40), Vec3f(20, 20, -1), Vec3f(0, 1, 0), scene.fov, scene.aspect_ratio, 0, 10, 0, 1);
    Renderer renderer;

    renderer.render(camera, scene);

    return 0;
}