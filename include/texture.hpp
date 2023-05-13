#pragma once

#include "utils.hpp"

class Texture
{
public:
    virtual ~Texture() {}
    virtual Color value(float u, float v, const Point3 &p) const = 0;
};

class SolidColor : public Texture
{
private:
    Color color_value;

public:
    SolidColor() {}
    SolidColor(Color c) : color_value(c) {}
    SolidColor(float red, float green, float blue)
        : SolidColor(Color(red, green, blue)) {}

    virtual Color value(float u, float v, const Point3 &p) const override
    {
        return color_value;
    }
};

class CheckerTexture : public Texture
{
private:
    shared_ptr<Texture> odd;
    shared_ptr<Texture> even;

public:
    CheckerTexture() {}
    CheckerTexture(shared_ptr<Texture> t0, shared_ptr<Texture> t1) : odd(t0), even(t1) {}

    virtual Color value(float u, float v, const Point3 &p) const override
    {
        auto sines = sin(10 * p.x) * sin(10 * p.y) * sin(10 * p.z); // a 3D checker pattern
        if (sines < 0)
            return odd->value(u, v, p);
        else
            return even->value(u, v, p);
    }
};

class ImageTexture : public Texture
{
private:
    unsigned char *data;
    int width, height;
    int bytes_per_scanline;

public:
    const static int bytes_per_pixel = 3; // RGB
    ImageTexture() : data(nullptr), width(0), height(0), bytes_per_scanline(0) {}
    ImageTexture(const char *filename)
    {
        auto components_per_pixel = bytes_per_pixel;

        // data stores the image data in a flat array, with the pixel at (i, j) being at data[bytes_per_scanline * i + bytes_per_pixel * j + k]
        data = stbi_load(filename, &width, &height, &components_per_pixel, components_per_pixel);

        if (!data)
        {
            std::cerr << "ERROR: Could not load texture image file '" << filename << "'.\n";
            width = height = 0;
        }

        bytes_per_scanline = bytes_per_pixel * width;
    }
    ~ImageTexture()
    {
        // delete data;
        stbi_image_free(data);
    }

    virtual Color value(float u, float v, const Point3 &p) const override
    {
        // If we have no texture data, then return solid cyan as a debugging aid.
        if (data == nullptr)
            return Color(0, 1, 1);

        // Clamp input texture coordinates to [0,1] x [1,0]
        u = clamp(u, 0.0f, 1.0f);
        v = 1.0f - clamp(v, 0.0f, 1.0f); // Flip V to image coordinates

        int i = static_cast<int>(u * width);
        int j = static_cast<int>(v * height);

        // Clamp integer mapping, since actual coordinates should be less than 1.0
        if (i >= width)
            i = width - 1;
        if (j >= height)
            j = height - 1;

        const float color_scale = 1.0f / 255.0f;
        auto pixel = data + j * bytes_per_scanline + i * bytes_per_pixel; // pointer to the first byte of the pixel

        return Color(color_scale * pixel[0], color_scale * pixel[1], color_scale * pixel[2]); // to [0,1]
    }
};

class Perlin
{
    static const int point_count = 256;
    Vec3f *ranvec; // random float array
    int *perm_x;   // permutation array
    int *perm_y;
    int *perm_z;

    static int *perlin_generate_perm()
    {
        auto p = new int[point_count];
        for (int i = 0; i < Perlin::point_count; ++i) // fill the array with 0, 1, ..., point_count - 1
            p[i] = i;
        permute(p, point_count); // shuffle the array
        return p;
    }

    static void permute(int *p, int n) // generate a permutation
    {
        for (int i = n - 1; i > 0; --i)
        {
            int target = random_int(0, i); // random index in [0, i]
            int tmp = p[i];                // swap p[i] and p[target]
            p[i] = p[target];
            p[target] = tmp;
        }
    }

    static float perlin_interp(Vec3f c[2][2][2], float u, float v, float w) // trilinear interpolation, c is the 2*2*2 cube of the grid points
    {
        auto uu = u * u * (3 - 2 * u); // smooth the interpolation
        auto vv = v * v * (3 - 2 * v);
        auto ww = w * w * (3 - 2 * w);
        auto accum = 0.0f;
        for (int i = 0; i < 2; ++i)
            for (int j = 0; j < 2; ++j)
                for (int k = 0; k < 2; ++k)
                {
                    Vec3f weight_v(u - i, v - j, w - k); // weight vector for the interpolation
                    accum += (i * uu + (1 - i) * (1 - uu)) *
                             (j * vv + (1 - j) * (1 - vv)) *
                             (k * ww + (1 - k) * (1 - ww)) * dot(c[i][j][k], weight_v);
                }

        return accum;
    }

public:
    Perlin()
    {
        ranvec = new Vec3f[point_count];
        for (int i = 0; i < point_count; ++i)
            ranvec[i] = random_vector3(-1, 1).normalize(); // use random unit vectors
        perm_x = perlin_generate_perm();
        perm_y = perlin_generate_perm();
        perm_z = perlin_generate_perm();
    }
    ~Perlin()
    {
        delete[] ranvec;
        delete[] perm_x;
        delete[] perm_y;
        delete[] perm_z;
    }

    float noise(const Point3 &p) const
    {
        // auto i = static_cast<int>(4 * p.x) & 255; // 1111 1111
        // auto j = static_cast<int>(4 * p.y) & 255;
        // auto k = static_cast<int>(4 * p.z) & 255;

        // return ranfloat[perm_x[i] ^ perm_y[j] ^ perm_z[k]];

        auto u = p.x - floor(p.x);
        auto v = p.y - floor(p.y);
        auto w = p.z - floor(p.z);
        // u = u * u * (3 - 2 * u); // use hermite cubic to smooth the value
        // v = v * v * (3 - 2 * v);
        // w = w * w * (3 - 2 * w);

        auto i = static_cast<int>(floor(p.x));
        auto j = static_cast<int>(floor(p.y));
        auto k = static_cast<int>(floor(p.z));
        Vec3f c[2][2][2];

        for (int di = 0; di < 2; ++di)
            for (int dj = 0; dj < 2; ++dj)
                for (int dk = 0; dk < 2; ++dk)
                    c[di][dj][dk] = ranvec[perm_x[(i + di) & 255] ^
                                           perm_y[(j + dj) & 255] ^
                                           perm_z[(k + dk) & 255]]; // get the 2*2*2 cube of the grid points

        return perlin_interp(c, u, v, w); // use trilinear interpolation to get the final value
    }

    float turb(const Point3 &p, int depth = 7) const // sum of noise
    {
        auto accum = 0.0f;
        auto temp_p = p;
        auto weight = 1.0f;

        for (int i = 0; i < depth; ++i)
        {
            accum += weight * noise(temp_p);
            weight *= 0.5;
            temp_p = temp_p * 2;
        }

        return fabs(accum);
    }
};

class NoiseTexture : public Texture
{
private:
    Perlin noise;
    float scale;

public:
    NoiseTexture() {}
    NoiseTexture(float sc) : scale(sc) {} // scale the input point

    virtual Color value(float u, float v, const Point3 &p) const override
    {
        // return Color(1, 1, 1) * 0.5 * (1 + noise.turb(scale * p)); // cast perlin output to [0,1]
        // return Color(1, 1, 1) * noise.turb(scale * p);             // turbulence
        return Color(1, 1, 1) * 0.5 * (1 + sin(scale * p.z + 10 * noise.turb(p))); // marble
    }
};
