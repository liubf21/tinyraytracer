#pragma once

#include <limits> // for std::numeric_limits
#include <memory> // for std::shared_ptr
#include <omp.h>  // for OpenMP
#include <random>    // for random number generator
#include <algorithm> // for std::sort

#include "stb_image_write.h"
#include "stb_image.h"

#include "vec.hpp"
#include "ray.hpp"

// Usings
using std::make_shared;
using std::shared_ptr;
using std::sqrt;

// Constants
const float epsilon = 0.0001f;
const float max_float = std::numeric_limits<float>::max();
const float min_float = std::numeric_limits<float>::min();
const float infinity = std::numeric_limits<float>::infinity();
const float pi = M_PI;

// Utility Functions
inline float degrees_to_radians(float degrees) // inline: to avoid multiple definitions
{
    return degrees * pi / 180.0;
}

inline float random_float()
{
    // Returns a random real in (0,1).
    static std::uniform_real_distribution<float> distribution(0.0, 1.0);
    // static std::mt19937 generator;
    static std::random_device rd; // Will be used to obtain a seed for the random number engine
    static std::mt19937 generator(rd()); // Standard mersenne_twister_engine seeded with rd()
    return distribution(generator);
}

inline float random_float(float min, float max)
{
    // Returns a random real in [min,max).
    return min + (max - min) * random_float();
}


inline float clamp(float x, float min, float max)
{
    if (x < min)
        return min;
    if (x > max)
        return max;
    return x;
}

inline int random_int(int min, int max)
{
    // Returns a random integer in [min,max].
    return static_cast<int>(clamp(random_float(min, max + 1), min, max));
}

inline Vec3f random_in_unit_sphere()
{
    while (true)
    {
        auto p = Vec3f(random_float(-1, 1), random_float(-1, 1), random_float(-1, 1));
        if (p.norm2() >= 1)
            continue;
        return p;
    }
}

inline Vec3f random_unit_vector()
{
    return random_in_unit_sphere().normalize();
}

inline Vec3f random_in_hemisphere(const Vec3f &normal)
{
    Vec3f in_unit_sphere = random_in_unit_sphere();
    if (dot(in_unit_sphere, normal) > 0.0) // In the same hemisphere as the normal
        return in_unit_sphere;
    else
        return -in_unit_sphere;
}

inline Vec3f random_in_unit_disk()
{
    while (true)
    {
        auto p = Vec3f(random_float(-1, 1), random_float(-1, 1), 0);
        if (p.norm2() >= 1)
            continue;
        return p;
    }
}

inline bool near_zero(const Vec3f &v)
{
    // Return true if the vector is close to zero in all dimensions.
    return (fabs(v.x) < epsilon) && (fabs(v.y) < epsilon) && (fabs(v.z) < epsilon);
}

inline Vec3f random_vector3(float min = 0.f, float max = 1.f)
{
    return Vec3f(random_float(min, max), random_float(min, max), random_float(min, max));
}

// I is the incident vector, N is the normal vector(unit vector), return the reflected vector in the direction of I
inline Vec3f reflect(const Vec3f &I, const Vec3f &N)
{
    return (I - N * 2.f * (I * N)).normalize();
}

// I is the incident vector, N is the normal vector(unit vector), return the refracted vector in the direction of I
// eta_t is the refractive index of the medium that the ray is entering, eta_i is the refractive index of the medium that the ray is leaving
// 折射进入玻璃可能的问题: 对于球，可以通过法向量方向判断光线射入还是射出，但是对于三角形，法向量方向不一定能判断光线射入还是射出
// 解决方法: 保证法向量指向物体外部
// 为了保证获得正确的法向量，模型导出的obj文件中，对每个三角面的描述里的三个顶点的顺序是固定的。并且按照“(顶点2-顶点0)×(顶点1-顶点0)”就能得到正确的法向量。而通常来说，这种定义好的法向量都会指向模型的外部。
inline Vec3f refract(const Vec3f &I, const Vec3f &N, const float &eta_t, const float &eta_i = 1.f)
{
    float cosi = -std::max(-1.f, std::min(1.f, I * N));
    if (cosi < 0)
        return refract(I, -N, eta_i, eta_t); // if the ray comes from the inside the object, swap the air and the media
    float eta = eta_i / eta_t;
    float cost_sq = 1 - eta * eta * (1 - cosi * cosi);
    return cost_sq < 0 ? reflect(I, N) : (I * eta + N * (eta * cosi - sqrtf(cost_sq))).normalize(); // : I*eta + N*(eta*cosi - sqrtf(cost_sq)); // k<0 = total reflection, no ray to refract
}
