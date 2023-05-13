#pragma once

#include "utils.hpp"
#include "primitive.hpp"
#include "texture.hpp"

struct ScatterRecord
{
    bool is_specular;
    Ray scattered;
    Color attenuation;
    float pdf;
};

class Material
{
public:
    virtual ~Material() {} // virtual destructor, so that derived classes can be deleted properly
    virtual bool scatter(const Ray &ray, const Intersection &intersection, ScatterRecord& srec) const { return false; }; // for path tracing
    // virtual Vec3f shade(const Ray &ray, const Intersection &intersection, const std::vector<shared_ptr<Light>> &lights, const std::vector<shared_ptr<Primitive>> &primitives, const int &depth) const = 0;
    virtual Color emitted(const float &u, const float &v, const Point3 &p) const { return Color(0.f, 0.f, 0.f); } // for light sources
    virtual float scattering_pdf(const Ray &ray, const Intersection &intersection, const Ray &scattered) const { return 0.f; } // return the pdf of the scattered ray
};

class BlinnPhong : public Material // used in ray tracing
{
    Vec3f albedo = Vec3f(0.f, 1.f, 1.f); // for Blinn-Phong shading
    Vec3f diffuse_color = Vec3f();       // color of the object
    float specular_exponent = 100.f;     // Phong specular exponent
    float refractive_index = 1.f;        // refractive index
    float Relf = 0.f;                    // reflection factor
    float Refr = 0.f;                    // refraction factor
public:
    BlinnPhong() {}
    BlinnPhong(const Vec3f &a, const Vec3f &dc, const float &se) : albedo(a), diffuse_color(dc), specular_exponent(se) {}                                       // for Blinn-Phong shading
    BlinnPhong(const float &se, const float &ri, const float &Relf, const float &Refr) : specular_exponent(se), refractive_index(ri), Relf(Relf), Refr(Refr) {} // for glass or mirror
    BlinnPhong(const Vec3f &a, const Vec3f &dc, const float &se, const float &ri, const float &Relf, const float &Refr) : albedo(a), diffuse_color(dc), specular_exponent(se), refractive_index(ri), Relf(Relf), Refr(Refr) {}

    const Vec3f &getAlbedo() const { return albedo; }
    const Vec3f &getDiffuseColor() const { return diffuse_color; }
    const float &getSpecularExponent() const { return specular_exponent; }
    const float &getRefractiveIndex() const { return refractive_index; }
    const float &getRelf() const { return Relf; }
    const float &getRefr() const { return Refr; }

};

class Lambertian : public Material
{
    // Color albedo; // reflectance in three channels
    shared_ptr<Texture> albedo; // Texture albedo instead of Color albedo
public:
    Lambertian(const Color &a) : albedo(make_shared<SolidColor>(a)) {}
    Lambertian(shared_ptr<Texture> a) : albedo(a) {}

    virtual bool scatter(const Ray &ray, const Intersection &intersection, ScatterRecord& srec) const override
    {
        // auto scatter_direction = intersection.normal + random_unit_vector();
        auto scatter_direction = random_in_hemisphere(intersection.normal);
        // auto scatter_direction = bump(intersection) + random_unit_vector();

        // Catch degenerate scatter direction
        // if (near_zero(scatter_direction))
        //     scatter_direction = intersection.normal;

        srec.scattered = Ray(intersection.position, scatter_direction, ray.time());
        srec.attenuation = albedo->value(intersection.u, intersection.v, intersection.position);
        // pdf = dot(intersection.normal, scattered.direction()) / pi; // sample from cosine distribution
        srec.pdf = 0.5f / pi;
        srec.is_specular = false;
        return true;
    }

    virtual float scattering_pdf(const Ray &ray, const Intersection &intersection, const Ray &scattered) const override
    {
        auto cosine = dot(intersection.normal, scattered.direction());
        return cosine < 0 ? 0 : cosine / pi;
    }

    Vec3f bump(const Intersection &intersection) const
    {
        // test bump mapping
        auto kh = 0.2, kn = 0.2;
        auto n = intersection.normal;
        auto x = n.x, y = n.y, z = n.z;
        // std::cout << "n: " << n << std::endl;
        // printf("n: %f %f %f\n", n.x, n.y, n.z);
        auto t = Vec3f(x * y / sqrt(x * x + z * z), sqrt(x * x + z * z), z * y / sqrt(x * x + z * z)); // tangent
        auto b = cross(intersection.normal, t);
        auto T = Vec3f(t.x, t.y, t.z);
        auto B = Vec3f(b.x, b.y, b.z);
        auto N = Vec3f(n.x, n.y, n.z);
        // std::cout << "T: " << T << std::endl;
        auto dU = kh * kn * (albedo->value(intersection.u + 0.01, intersection.v, intersection.position).x - albedo->value(intersection.u, intersection.v, intersection.position).x);
        auto dV = kh * kn * (albedo->value(intersection.u, intersection.v + 0.01, intersection.position).x - albedo->value(intersection.u, intersection.v, intersection.position).x);
        auto ln = Vec3f(-dU, -dV, 1.f);
        // std::cout << "ln: " << ln << std::endl;
        auto new_normal = Vec3f(dot(T, ln), dot(B, ln), dot(N, ln)).normalize();
        // auto new_normal = intersection.normal + Vec3f(random_float(0, 0.9), random_float(0, 0.9), random_float(0, 0.9));
        return new_normal.normalize();
    }
};

class Metal : public Material
{
    Color albedo; // reflectance in three channels
    float fuzz; // fuzziness of the metal, if 0, it is a perfect mirror
public:
    Metal(const Color &a, float f) : albedo(a), fuzz(f < 1 ? f : 1) {}

    virtual bool scatter(const Ray &ray, const Intersection &intersection, ScatterRecord& srec) const override
    {
        Vec3f reflected = reflect(ray.direction(), intersection.normal);
        srec.scattered = Ray(intersection.position + (reflected * intersection.normal > 0 ? 1e-3 * intersection.normal : -1e-3 * intersection.normal
            ), reflected + fuzz * random_in_unit_sphere(), ray.time());
        srec.attenuation = albedo;
        srec.is_specular = true;
        srec.pdf = 1;
        return (dot(srec.scattered.direction(), intersection.normal) > 0); // if the scattered ray is not in the same hemisphere as the normal, it is not reflected
    }
};

class Dielectric : public Material
{
    Color albedo; // reflectance in three channels
    float ir; // Index of Refraction
    static float reflectance(float cosine, float ref_idx)
    {
        // Use Schlick's approximation for reflectance.
        auto r0 = (1 - ref_idx) / (1 + ref_idx);
        r0 = r0 * r0;
        return r0 + (1 - r0) * pow((1 - cosine), 5);
    }
public:
    Dielectric(const Color &a, float index_of_refraction) : albedo(a), ir(index_of_refraction) {}
    Dielectric(float index_of_refraction) : albedo(Color(1, 1, 1)), ir(index_of_refraction) {}

    virtual bool scatter(const Ray &ray, const Intersection &intersection, ScatterRecord& srec) const override
    {
        srec.attenuation = albedo;
        // float refraction_ratio = intersection.front_face ? (1.0 / ir) : ir;

        Vec3f unit_direction = ray.direction();
        float cos_theta = fmin(dot(-unit_direction, intersection.normal), 1.0);
        // float sin_theta = sqrt(1.0 - cos_theta * cos_theta);
        // bool cannot_refract = refraction_ratio * sin_theta > 1.0;
        Vec3f direction;

        // if (reflectance(cos_theta, ir) > random_float())
            // direction = reflect(unit_direction, intersection.normal);
        // else
            direction = refract(unit_direction, intersection.normal, ir);

        srec.scattered = Ray(intersection.position + (direction * intersection.normal > 0 ? 1e-3 * intersection.normal : -1e-3 * intersection.normal
            ), direction, ray.time());
        srec.is_specular = true;
        srec.pdf = 1;
        return true;
    }
};

class DiffuseLight : public Material
{
    shared_ptr<Texture> emit; // Texture albedo instead of Color albedo
public:
    DiffuseLight(shared_ptr<Texture> a) : emit(a) {}
    DiffuseLight(Color c) : emit(make_shared<SolidColor>(c)) {}

    virtual bool scatter(const Ray &ray, const Intersection &intersection, ScatterRecord& srec) const override
    {
        return false;
    }

    virtual Color emitted(const float &u, const float &v, const Point3 &p) const override
    {
        return emit->value(u, v, p);
    }
};

class Isotropic : public Material // for fog
{
    shared_ptr<Texture> albedo;

public:
    Isotropic(shared_ptr<Texture> a) : albedo(a) {}
    Isotropic(Color c) : albedo(make_shared<SolidColor>(c)) {}

    virtual bool scatter(const Ray &ray, const Intersection &intersection, ScatterRecord& srec) const override
    {
        srec.scattered = Ray(intersection.position, random_in_unit_sphere(), ray.time());
        srec.attenuation = albedo->value(intersection.u, intersection.v, intersection.position);
        srec.is_specular = false;
        srec.pdf = 1. / (4. * M_PI);
        return true;
    }

    virtual float scattering_pdf(const Ray &ray, const Intersection &intersection, const Ray &scattered) const override
    {
        return 1. / (4. * M_PI);
    }
};

