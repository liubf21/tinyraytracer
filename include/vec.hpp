#pragma once
#include <cmath>
#include <cassert>
#include <iostream>

template <size_t DIM, typename T> // DIM: dimension, T: type
class vec
{
    T data_[DIM];

public:
    vec()
    {
        for (size_t i = DIM; i--; data_[i] = T())
            ;
    }
    T &operator[](const size_t i)
    {
        assert(i < DIM);
        return data_[i];
    }
    const T &operator[](const size_t i) const
    {
        assert(i < DIM);
        return data_[i];
    }
};

// Type aliases
using Vec2f = vec<2, float>;
using Vec3f = vec<3, float>;
using Vec3i = vec<3, int>;
using Vec4f = vec<4, float>;

using Point3 = Vec3f; // 3D point
using Color = Vec3f;  // RGB color

template <typename T>
class vec<2, T>
{
public:
    vec() : x(T()), y(T()) {}
    vec(T X, T Y) : x(X), y(Y) {}
    template <class U>
    vec<2, T>(const vec<2, U> &v) : x(T(v.x)), y(T(v.y)) {} // casting constructor
    T &operator[](const size_t i)
    {
        assert(i < 2);
        return i <= 0 ? x : y;
    }
    const T &operator[](const size_t i) const
    {
        assert(i < 2);
        return i <= 0 ? x : y;
    }
    T x, y;
};

template <typename T>
class vec<3, T>
{
public:
    vec() : x(T()), y(T()), z(T()) {}
    vec(T X, T Y, T Z) : x(X), y(Y), z(Z) {}
    T &operator[](const size_t i)
    {
        assert(i < 3);
        return i <= 0 ? x : (1 == i ? y : z);
    }
    const T &operator[](const size_t i) const
    {
        assert(i < 3);
        return i <= 0 ? x : (1 == i ? y : z);
    }
    T norm2() const { return x * x + y * y + z * z; }
    T norm() const { return sqrt(norm2()); } // length of the vector
    vec<3, T> &normalize(T l = 1.0)
    { // std::cout<<*this<<":"<<norm()<<std::endl;
        *this = (*this) * (l / norm());
        return *this;
    }
    T x, y, z;
};

template <typename T>
class vec<4, T>
{
public:
    vec() : x(T()), y(T()), z(T()), w(T()) {}
    vec(T X, T Y, T Z, T W) : x(X), y(Y), z(Z), w(W) {}
    T &operator[](const size_t i)
    {
        assert(i < 4);
        return i <= 0 ? x : (1 == i ? y : (2 == i ? z : w));
    }
    const T &operator[](const size_t i) const
    {
        assert(i < 4);
        return i <= 0 ? x : (1 == i ? y : (2 == i ? z : w));
    }
    T x, y, z, w;
};

// functions
template <size_t DIM, typename T>
T operator*(const vec<DIM, T> &lhs, const vec<DIM, T> &rhs)
{ // dot product
    T ret = T();
    for (size_t i = DIM; i--; ret += lhs[i] * rhs[i])
        ;
    return ret;
}

template <size_t DIM, typename T>
vec<DIM, T> operator+(vec<DIM, T> lhs, const vec<DIM, T> &rhs)
{
    for (size_t i = DIM; i--; lhs[i] += rhs[i])
        ;
    return lhs;
}

template <size_t DIM, typename T>
vec<DIM, T> operator-(vec<DIM, T> lhs, const vec<DIM, T> &rhs)
{
    for (size_t i = DIM; i--; lhs[i] -= rhs[i])
        ;
    return lhs;
}

template <size_t DIM, typename T, typename U>
vec<DIM, T> operator+(vec<DIM, T> lhs, const U &rhs)
{
    for (size_t i = DIM; i--; lhs[i] += rhs)
        ;
    return lhs;
}

template <size_t DIM, typename T, typename U>
vec<DIM, T> operator-(vec<DIM, T> lhs, const U &rhs)
{
    for (size_t i = DIM; i--; lhs[i] -= rhs)
        ;
    return lhs;
}

template <size_t DIM, typename T, typename U>
vec<DIM, T> operator*(const vec<DIM, T> &lhs, const U &rhs)
{ // scalar product
    vec<DIM, T> ret;
    for (size_t i = DIM; i--; ret[i] = lhs[i] * rhs)
        ;
    return ret;
}

template <size_t DIM, typename T, typename U>
vec<DIM, T> operator*(const U &lhs, const vec<DIM, T> &rhs)
{
    return rhs * lhs;
}

template <size_t DIM, typename T>
vec<DIM, T> operator-(const vec<DIM, T> &lhs)
{
    return lhs * T(-1);
}

template <size_t DIM, typename T, typename U>
vec<DIM, T> operator/(const vec<DIM, T> &lhs, const U &rhs)
{
    vec<DIM, T> ret;
    for (size_t i = DIM; i--; ret[i] = lhs[i] / rhs)
        ;
    return ret;
}

template <typename T>
T dot(const vec<3, T> &v1, const vec<3, T> &v2)
{ // dot product
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;
}

template <typename T>
vec<3, T> cross(vec<3, T> v1, vec<3, T> v2)
{ // cross product
    return vec<3, T>(v1.y * v2.z - v1.z * v2.y, v1.z * v2.x - v1.x * v2.z, v1.x * v2.y - v1.y * v2.x);
}

template <size_t DIM, typename T>
vec<DIM, T> wiseProduct(const vec<DIM, T> &v1, const vec<DIM, T> &v2)
{ // element-wise product
    vec<DIM, T> ret;
    for (size_t i = DIM; i--; ret[i] = v1[i] * v2[i])
        ;
    return ret;
}

template <size_t DIM, typename T>
vec<DIM, T> min(const vec<DIM, T> &v1, const vec<DIM, T> &v2)
{ // element-wise min
    vec<DIM, T> ret;
    for (size_t i = DIM; i--; ret[i] = std::min(v1[i], v2[i]))
        ;
    return ret;
}

template <size_t DIM, typename T>
vec<DIM, T> max(const vec<DIM, T> &v1, const vec<DIM, T> &v2)
{ // element-wise max
    vec<DIM, T> ret;
    for (size_t i = DIM; i--; ret[i] = std::max(v1[i], v2[i]))
        ;
    return ret;
}

template <size_t DIM, typename T>
std::ostream &operator<<(std::ostream &out, const vec<DIM, T> &v)
{
    for (unsigned int i = 0; i < DIM; i++)
    {
        out << v[i] << " ";
    }
    return out;
}
