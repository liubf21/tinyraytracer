#pragma once

#include <vector>
#include <string>
#include "ray.hpp"

class Model 
{
private:
    std::vector<Vec3f> verts; // list of vertices
    std::vector<Vec3i> faces; // indices of vertices for each triangle
    Vec3f m_min, m_max; // bounding box, computed in the constructor
    // void get_bbox(Vec3f &min, Vec3f &max); // bounding box for all the vertices, including isolated ones
public:
    Model(const char *filename);

    int nverts() const;                          // number of vertices
    int nfaces() const;                          // number of triangles

    const Vec3f &point(int i) const;                   // coordinates of the vertex i
    Vec3f &point(int i);                   // coordinates of the vertex i
    int vert(int fi, int li) const;              // index of the vertex for the triangle fi and local index li

};

std::ostream& operator<<(std::ostream& out, Model &m);

