# tinyraytracer

Welcome to TinyRayTracer, a ray tracing project written in C++. This ray tracer is capable of rendering various objects such as spheres, planes, cubes, and triangles with diffuse and specular lighting. It also supports advanced features such as shadows, reflections, and refractions.

## Features

- Algorithm: Whitted ray tracing algorithm and Monte Carlo path tracing algorithm
- Primitives: sphere, plane, triangle, cube
- Materials: Lambertian, metal, dielectric
- Textures: checkerboard, image, perlin noise
- Lighting: point light, directional light, ambient light
- Acceleration structure: bounding box, BVH
- Other features: multi-threading using OpenMP, gobal illumination, soft shadows, anti-aliasing, depth of field, motion blur, smoke etc.

## Installation Instructions

### Prerequisites

- C++ compiler (gcc, g++, clang, etc.)
- CMake 3.12 or above

### Build and Run the Program

1. Clone the project repository to your local machine using Git.
```
git clone https://github.com/liubf21/tinyraytracer.git
```

2. Navigate to the project directory.
```
cd tinyraytracer
```

3. Build the project using CMake.
```
mkdir build
cd build
cmake ..
```

4. Compile the source code using Make.
```
make
```

5. Run the program.
```
./tinyraytracer
```

## Usage

use -h or --help to get help information

## Examples

![image1](output/pt1.jpg)
(more examples in output folder)

## References
[GAMES101](https://sites.cs.ucsb.edu/~lingqi/teaching/games101.html)

[smallpt](http://www.kevinbeason.com/smallpt/)

[tinyraytracer](https://github.com/ssloy/tinyraytracer/wiki)

[raytracing in one weekend](https://raytracing.github.io)

[minilight](https://www.hxa.name/minilight/)
