#pragma once

#include "utils.hpp"
#include "aabb.hpp"

class Material; // forward declaration

struct Intersection
{
    Point3 position;
    Vec3f normal; // should be normalized, always pointing outwards
    shared_ptr<Material> material;
    float t;
    float u, v; // texture coordinates
};

class Primitive
{
protected:
    shared_ptr<Material> material;

public:
    Primitive() {}
    Primitive(shared_ptr<Material> m) : material(m) {}
    virtual ~Primitive() {}

    virtual bool ray_intersect(const Ray &ray, float t_min, float t_max, Intersection &intersection) const = 0; // t_min and t_max are used to limit the range of t
    // shared_ptr<Material> getMaterial() const { return material; }
    // virtual Vec3f getNormal(const Vec3f& hitpoint) const = 0; // the direction should be checked

    virtual bool bounding_box(float time0, float time1, AABB &output_box) const = 0;
    virtual float pdf_value(const Point3 &o, const Vec3f &v) const { return 0.0; } // v is the direction of the ray, must be normalized
    virtual Vec3f random() const { return Vec3f(1, 0, 0); }                        // return a random point on the primitive
    virtual Vec3f random(const Vec3f &o) const { return Vec3f(1, 0, 0); }          // return a random direction from the primitive to the origin
};

class Light : public Primitive // point light
{
    Vec3f position;
    float intensity;

public:
    Light(const Vec3f &p, const float &i) : position(p), intensity(i) {}
    const Vec3f &getPosition() const { return position; }
    const float &getIntensity() const { return intensity; }
    virtual bool ray_intersect(const Ray &ray, float t_min, float t_max, Intersection &intersection) const override
    {
        return false;
    }
    virtual bool bounding_box(float time0, float time1, AABB &output_box) const override
    {
        return false;
    }
    virtual float pdf_value(const Point3 &o, const Vec3f &v) const override
    {
        if ((position - o).normalize() * v > 0.999)
            return 1.0;
        else
            return 0.0;
    }
    virtual Vec3f random(const Vec3f &o) const override { return position - o; }
};

class Sphere : public Primitive
{
    Vec3f center;
    float radius;
    static void get_sphere_uv(const Point3 &p, float &u, float &v) // p is the normal vector
    {
        // p: a given point on the sphere of radius one, centered at the origin.
        // u: returned value [0,1] of angle around the Y axis from X=-1.
        // v: returned value [0,1] of angle from Y=-1 to Y=+1.
        //     <1 0 0> yields <0.50 0.50>       <-1  0  0> yields <0.00 0.50>
        //     <0 1 0> yields <0.50 1.00>       < 0 -1  0> yields <0.50 0.00>
        //     <0 0 1> yields <0.25 0.50>       < 0  0 -1> yields <0.75 0.50>

        float theta = acosf(-p.y);
        float phi = atan2f(-p.z, p.x) + pi;

        u = phi / (2 * pi);
        v = theta / pi;
    }

public:
    Sphere(const Vec3f &c, const float &r, shared_ptr<Material> m) : Primitive(m), center(c), radius(r) {}

    virtual bool ray_intersect(const Ray &ray, float t_min, float t_max, Intersection &intersection) const override
    {
        Vec3f L = center - ray.origin();
        float tca = L * ray.direction(); // ray direction should be normalized
        float d2 = L * L - tca * tca;
        if (d2 > radius * radius)
            return false;
        float thc = sqrtf(radius * radius - d2);
        auto t0 = tca - thc;
        if (t0 < t_min || t0 > t_max)
            t0 = tca + thc;
        if (t0 < t_min || t0 > t_max)
            return false; // return false if both t0 and t1 are out of range
        intersection.t = t0;
        intersection.position = ray.origin() + ray.direction() * t0;
        intersection.normal = (intersection.position - center).normalize(); // should be normalized, always pointing outwards
        intersection.material = material;
        get_sphere_uv(intersection.normal, intersection.u, intersection.v);
        // std::cout << intersection.normal << std::endl;

        return true;
    }

    virtual bool bounding_box(float time0, float time1, AABB &output_box) const override
    {
        output_box = AABB(center - Vec3f(radius, radius, radius), center + Vec3f(radius, radius, radius));
        return true;
    }
};

class MovingSphere : public Primitive
{
    Vec3f center0, center1;
    float time0, time1;
    float radius;

public:
    MovingSphere(const Vec3f &c0, const Vec3f &c1, const float &t0, const float &t1, const float &r, shared_ptr<Material> m) : Primitive(m), center0(c0), center1(c1), time0(t0), time1(t1), radius(r) {}

    virtual bool ray_intersect(const Ray &ray, float t_min, float t_max, Intersection &intersection) const override
    {
        Vec3f L = center(ray.time()) - ray.origin();
        float tca = L * ray.direction(); // ray direction should be normalized
        float d2 = L * L - tca * tca;
        if (d2 > radius * radius)
            return false;
        float thc = sqrtf(radius * radius - d2);
        auto t0 = tca - thc;
        if (t0 < t_min || t0 > t_max)
            t0 = tca + thc;
        if (t0 < t_min || t0 > t_max)
            return false; // return false if both t0 and t1 are out of range
        intersection.t = t0;
        intersection.position = ray.origin() + ray.direction() * t0;
        intersection.normal = (intersection.position - center(ray.time())).normalize(); // should be normalized, always pointing outwards
        intersection.material = material;
        return true;
    }

    virtual bool bounding_box(float time0, float time1, AABB &output_box) const override
    {
        AABB box0(center0 - Vec3f(radius, radius, radius), center0 + Vec3f(radius, radius, radius));
        AABB box1(center1 - Vec3f(radius, radius, radius), center1 + Vec3f(radius, radius, radius));
        output_box = surrounding_box(box0, box1);
        return true;
    }

    Vec3f center(float time) const
    {
        return center0 + ((time - time0) / (time1 - time0)) * (center1 - center0);
    }
};

class Plane : public Primitive
{
    Vec3f normal; // normal vector of the plane, pointing from the origin
    Point3 p[4];  // four points on the plane
    bool isInf = false;

public:
    Plane(Vec3f n, Point3 p0, shared_ptr<Material> m) : Primitive(m), normal(n.normalize()), p{p0, p0, p0, p0} { isInf = true; }
    Plane(Point3 p0, Point3 p1, Point3 p2, Point3 p3, shared_ptr<Material> m) : Primitive(m), p{p0, p1, p2, p3} { normal = cross(p1 - p0, p2 - p0).normalize(); } // hasn't insure the points are on the same plane

    virtual bool ray_intersect(const Ray &ray, float t_min, float t_max, Intersection &intersection) const override
    {
        float denom = normal * ray.direction(); // dot product of the normal and the direction of the ray
        if (denom > 1e-6 || denom < -1e-6)
        {
            auto t0 = (normal * p[0] - normal * ray.origin()) / denom; // dot product of the normal and the vector from the origin to the point on the plane
            if (t0 >= t_min && t0 <= t_max)
            {
                Point3 temp = ray.origin() + ray.direction() * t0;
                if (!isInf)
                {
                    // bool inside = true;
                    int tmp = 1;
                    for (int i = 0; i < 4; i++)
                    {
                        Vec3f v1 = p[(i + 1) % 4] - p[i]; // the order of the points matters
                        Vec3f v2 = temp - p[i];
                        Vec3f out = cross(v1, v2);
                        tmp *= dot(out, normal) > 0 ? 1 : -1;
                        if (tmp < 0)
                        {
                            return false;
                        }
                    }
                }
                intersection.t = t0;
                intersection.position = temp;
                intersection.normal = dot(ray.direction(), normal) < 0 ? normal : normal * -1; // note: direction of the normal should be determined by the direction of the ray
                intersection.material = material;
                return true;
            }
        }
        return false;
    }

    virtual bool bounding_box(float time0, float time1, AABB &output_box) const override
    {
        if (isInf)
            return false;
        float min_[3], max_[3];
        for (int i = 0; i < 3; i++)
        {
            min_[i] = fmin(p[0][i], fmin(p[1][i], fmin(p[2][i], p[3][i])));
            max_[i] = fmax(p[0][i], fmax(p[1][i], fmax(p[2][i], p[3][i])));
        }
        output_box = AABB(Point3(min_[0], min_[1], min_[2]), Point3(max_[0], max_[1], max_[2]));
        return true;
    }

    virtual float pdf_value(const Point3 &o, const Vec3f &v) const override
    {
        Intersection intersection;
        if (this->ray_intersect(Ray(o, v), epsilon, max_float, intersection))
        {
            float area = cross(p[1] - p[0], p[2] - p[0]).norm();
            float distance_squared = intersection.t * intersection.t * v.norm2();
            float cosine = fabs(dot(v, intersection.normal) / v.norm());
            return distance_squared / (cosine * area);
        }
        else
            // std::cout << "!";
            return 0;
    }

    virtual Vec3f random() const override // random point on the plane
    {
        // interpolate the point on the plane
        auto p1 = p[0] + random_float() * (p[1] - p[0]);
        auto p2 = p[3] + random_float() * (p[2] - p[3]);
        return p1 + random_float() * (p2 - p1);
    }

    virtual Vec3f random(const Vec3f &o) const override
    {
        return random() - o;
    }
};

class Triangle : public Primitive
{
    Vec3f v0, v1, v2;
    Vec3f normal;
    Vec3f edge1, edge2;

public:
    Triangle(const Vec3f &v0, const Vec3f &v1, const Vec3f &v2, shared_ptr<Material> m) : Primitive(m), v0(v0), v1(v1), v2(v2)
    {
        edge1 = v1 - v0;
        edge2 = v2 - v0;
        normal = cross(edge1, edge2).normalize(); // normal vector of the triangle, pointing outwards
    }

    virtual bool ray_intersect(const Ray &ray, float t_min, float t_max, Intersection &intersection) const override
    {
        Vec3f pvec = cross(ray.direction(), edge2);
        float det = edge1 * pvec;
        if (det < 1e-5)
            return false;

        Vec3f tvec = ray.origin() - v0;
        float u = tvec * pvec;
        if (u < 0 || u > det)
            return false;

        Vec3f qvec = cross(tvec, edge1);
        float v = ray.direction() * qvec;
        if (v < 0 || u + v > det)
            return false;

        auto tnear = edge2 * qvec * (1. / det);

        if (tnear < t_min || tnear > t_max)
            return false;

        intersection.t = tnear;
        intersection.position = ray.origin() + ray.direction() * tnear;
        intersection.normal = normal;
        intersection.material = material;
        return true;
    }

    virtual bool bounding_box(float time0, float time1, AABB &output_box) const override
    {
        Vec3f min = Vec3f(std::min(v0.x, std::min(v1.x, v2.x)), std::min(v0.y, std::min(v1.y, v2.y)), std::min(v0.z, std::min(v1.z, v2.z)));
        Vec3f max = Vec3f(std::max(v0.x, std::max(v1.x, v2.x)), std::max(v0.y, std::max(v1.y, v2.y)), std::max(v0.z, std::max(v1.z, v2.z)));
        output_box = AABB(min, max);
        return true;
    }
};

class Translate : public Primitive // move the object
{
    shared_ptr<Primitive> primitive;
    Vec3f offset;

public:
    Translate(shared_ptr<Primitive> p, const Vec3f &displacement) : primitive(p), offset(displacement) {}

    virtual bool ray_intersect(const Ray &ray, float t_min, float t_max, Intersection &intersection) const override
    {
        Ray moved_ray(ray.origin() - offset, ray.direction(), ray.time()); // move the ray in the opposite direction of the object
        if (!primitive->ray_intersect(moved_ray, t_min, t_max, intersection))
            return false;

        intersection.position = intersection.position + offset; // move the intersection point in the direction of the object

        return true;
    }

    virtual bool bounding_box(float time0, float time1, AABB &output_box) const override
    {
        if (!primitive->bounding_box(time0, time1, output_box))
            return false;

        output_box = AABB(output_box.min() + offset, output_box.max() + offset);

        return true;
    }
};

class RotateY : public Primitive
{
    shared_ptr<Primitive> primitive;
    float sin_theta;
    float cos_theta;
    bool has_box;
    AABB bbox;

public:
    RotateY(shared_ptr<Primitive> p, float angle) : primitive(p)
    {
        float radians = degrees_to_radians(angle);
        sin_theta = sin(radians);
        cos_theta = cos(radians);
        has_box = primitive->bounding_box(0, 1, bbox);
        Point3 min(infinity, infinity, infinity);
        Point3 max(-infinity, -infinity, -infinity);
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 3; ++j)
            {
                float x = i * bbox.max().x + (1 - i) * bbox.min().x;
                float y = j * bbox.max().y + (1 - j) * bbox.min().y;
                float z = i * bbox.max().z + (1 - i) * bbox.min().z;
                float newx = cos_theta * x + sin_theta * z;
                float newz = -sin_theta * x + cos_theta * z;
                Vec3f tester(newx, y, newz); // rotate the bounding box
                for (int c = 0; c < 3; c++)
                {
                    min[c] = fmin(min[c], tester[c]);
                    max[c] = fmax(max[c], tester[c]);
                }
            }
        }
        bbox = AABB(min, max);
    }

    virtual bool ray_intersect(const Ray &ray, float t_min, float t_max, Intersection &intersection) const override
    {
        auto origin = ray.origin();
        auto direction = ray.direction();

        origin[0] = cos_theta * ray.origin()[0] - sin_theta * ray.origin()[2]; // rotate the ray
        origin[2] = sin_theta * ray.origin()[0] + cos_theta * ray.origin()[2];

        direction[0] = cos_theta * ray.direction()[0] - sin_theta * ray.direction()[2];
        direction[2] = sin_theta * ray.direction()[0] + cos_theta * ray.direction()[2];

        Ray rotated_ray(origin, direction, ray.time());

        if (!primitive->ray_intersect(rotated_ray, t_min, t_max, intersection))
            return false;

        auto position = intersection.position;
        auto normal = intersection.normal;

        position[0] = cos_theta * intersection.position[0] + sin_theta * intersection.position[2]; // rotate the intersection point
        position[2] = -sin_theta * intersection.position[0] + cos_theta * intersection.position[2];

        normal[0] = cos_theta * intersection.normal[0] + sin_theta * intersection.normal[2]; // rotate the normal
        normal[2] = -sin_theta * intersection.normal[0] + cos_theta * intersection.normal[2];

        intersection.position = position;

        return true;
    }

    virtual bool bounding_box(float time0, float time1, AABB &output_box) const override
    {
        output_box = bbox;
        return has_box;
    }
};
