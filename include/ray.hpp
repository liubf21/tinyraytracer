#pragma once

#include "vec.hpp"

class Ray {
private:
    Point3 orig;
    Vec3f dir; // should be normalized
    float tm; // time, for motion blur
public:
    Ray() {}
    Ray(const Point3 &origin, Vec3f direction, float time = 0.0f) : orig(origin), dir(direction.normalize()), tm(time) {//std::cout<<"ray";
    }
    
    const Point3& origin() const { return orig; }
    const Vec3f& direction() const { return dir; }
    float time() const { return tm; }
    Point3 at(float t) const { return orig + dir*t; }
    
};
