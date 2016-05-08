# Yocto/GL: C99 Single File Libraries for Physically-Based Graphics

Yocto/GL is a collection of single-file libraries for building
physically-based graphics applications.
Yocto/GL is written in C99 without dependencies and compiles on
OSX (clang/gcc), Linux (clang/gcc) and Windows (cl).

## Main Libraries

- **yocto_obj.h** - Wavefront OBJ/MTL loader that supports arbitrary polygons (with/without triangulation), lines, and points. Includes optionals extensions for per-vertex color and radius, camera and environment map. Optionally depends on `stb_image.h` for texture loading.
- **yocto_bvh.h** - Ray casting and closet point queries of points, lines and triangles accelerated by a two-level bounding volume hierarchy.
- **yocto_shape.h** - Utilities for manipulating shapes composed of points, lines, triangles or quads. Includes parametric shape generation, uniform tesselation, normal computation, uniform shape sampling.
- **yocto_trace.h** - Path tracer with support for point, line or triangle geometry, mesh area lights and environment maps, materials with either GGX or Phong (only opaque for now). Support both incremental and offline computation on single- or multi-core machines.
- **yocto_rigid.h** - Rigid body solver supporting convex and concave triangle meshes based on Sequential Impulses (aka Projected Gauss-Sidel).

## Support Libraries

- **yocto_cmdline.h** - Utilities for writing command line applications. Includes in  particular a command line parsing library that support options and arguments of ints, floats, strings, enums.
- **yocto_glu.h** - Quick and dirty rendering of images and shapes in OpenGL, useful to create interactive viewers.
- **yocto_math.h** - A few vector math routines helpful to write test apps. This is not complete and not suitable for others use.

## Examples

This repository contains Yocto/GL applications written to test the libraries.

- **ytestgen.c**: Creates various test cases for the path tracer and GL viewer.
- **yview.c**: HDR/PNG/JPG image viewer with exposure/gamma tone mapping.
- **yshade.c**: OpenGL viewer.
- **ytrace.c**: Interactive path-tracer, that can also run in offline mode.
- **ytrace.c**: Simple rigid body demo code. Does not work on offline mode.

A few screenshots from **ytrace** are included here for demonstration.

![](images/sh03.path.png)
![](images/ls02.direct.png)
![](images/cb01.path.png)
![](images/rs02.path.png)

A screenshotted movie from **ysym** is included here for demonstration.

![](images/rb02.ysym.gif)

## Possible Future Development

- Likely moving to C++ for internal implementation, but still maintaining C99 interfaces.
    - While we like the pure C99 code, the lack of operator overloaing makes the code very hard to read. This is a drawback for people that want to understand the code.
    - Also a judicious use of templates would remove a lot of redundand code.
- Likely moving a common math library that for now is not shared between YOCTO/GL. 
    - This is for readability, but each single YOCTO/GL library will remain independent from the others.
- Likely removing support for quads and moving to simplex-based geometry representation.
- Particle-based simluation coming soon.

## License

Yocto/GL libraries are released under the permissive MIT license, while the
example apps are released under the 2-clause BSD (to include warranty for binary distribution).
