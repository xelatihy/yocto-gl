# Yocto/GL: Tiny C++ Libraries for Data-Driven Physically-based Graphics

[![Build Status](https://travis-ci.org/xelatihy/yocto-gl.svg?branch=master)](https://travis-ci.org/xelatihy/yocto-gl) [![Build status](https://ci.appveyor.com/api/projects/status/rkqw7a8cenl877m6/branch/master?svg=true)](https://ci.appveyor.com/project/xelatihy/yocto-gl/branch/master)

Yocto/GL is a collection of utility C++17 libraries for building 
physically-based graphics algorithms released under the MIT license.
Yocto/GL is written in a deliberatly data-driven style for ease of
development and use. 

Features include:
- convenience math functions for graphics
- static length vectors for 2, 3, 4 length of arbitrary type
- static length matrices for 2x2, 3x3, 4x4 of arbitrary type
- static length rigid transforms (frames), specialized for 2d and 3d space
- linear algebra operations and transforms
- axis aligned bounding boxes
- rays and ray-primitive intersection
- point-primitive distance and overlap tests
- normal and tangent computation for meshes and lines
- generation of tesselated meshes
- mesh refinement with linear tesselation and Catmull-Cark subdivision
- keyframed animation, skinning and morphing
- random number generation via PCG32
- simple image data structure and a few image operations
- simple scene format
- generation of image examples
- generation of scene examples
- procedural sun and sky HDR
- procedural Perlin noise
- BVH for intersection and closest point query
- Python-like path operations
- immediate mode command line parser
- path tracer supporting surfaces and hairs, GGX and MIS
- support for loading and saving Wavefront OBJ and Khronos glTF
- fast, hackable, proprietary JSON format

Yocto/GL is written in C++17 and compiles on OSX (clang from Xcode 10+),
Linux (gcc 6+, clang 6+) and Windows (MSVC 2015, MSVC 2017). For compilation
options, check the individual libraries.

Here are two images rendered with the builtin path tracer, where the
scenes are crated with the test generator.

![Yocto/GL](images/shapes.png)

![Yocto/GL](images/lines.png)


## Credits

This library includes code from the PCG random number generator,
boost hash_combine, base64 encode/decode by René Nyffenegger and 
public domain code from github.com/sgorsten/linalg, 
gist.github.com/badboy/6267743 and github.com/nothings/stb_perlin.h.
Other external libraries are included with their own license.


## Libraries

Yocto/GL is split into two small libraries to make code navigation easier.
See each header file for documentation.

- `yocto/yocto_math.{h}`: fixed size vectors, matrices, frames, rays, 
   bounding boxes, transforms
- `yocto/yocto_random.{h}`: random number generation, Perlin noise, Monte Carlo
   utilities
- `yocto/yocto_shape.{h}`: geometry utilities, shape manipulation, 
   procedural shapes
- `yocto/yocto_bvh.{h}`: ray intersection and closest point queries 
   using a two-level bounding volume hierarchy
- `yocto/yocto_image.{h,cpp}`: color utilities, image manipulation, 
   procedural images, procedural sun-sky, image input/output
- `yocto/yocto_scene.{h,cpp}`: simple scene storage, evaluation of scene 
   properties
- `yocto/yocto_trace.{h,cpp}`: path tracing
- `yocto/yocto_utils.{h}`: printing and parsing values, path utlities, file io,
   command line parsing
- `yocto/yocto_obj.{h}`: OBJ parser based on callbacks (SAX-like)
- `yocto/yocto_pbrt.{h}`: pbrt parser based on callbacks (SAX-like)
- `yocto/yocto_imageio.{h,cpp}`: image loading and saving
- `yocto/yocto_sceneio.{h,cpp}`: scene loading and saving


## Example Applications

You can see Yocto/GL in action in the following applications written to
test the library:

- `apps/yview.cpp`: simple OpenGL viewer
- `apps/ytrace.cpp`: offline path-tracer
- `apps/yitrace.cpp`: interactive path-tracer
- `apps/yscnproc.cpp`: scene manipulation and conversion to/from OBJ and glTF
- `apps/yimview.cpp`: HDR/PNG/JPG image viewer with exposure/gamma tone mapping
- `apps/yimproc.cpp`: offline image manipulation.

You can build the example applications using CMake with
    `mkdir build; cd build; cmake ..; cmake --build`


## Compilation

This library requires a C++17 compiler and is know to compiled on 
OsX (Xcode >= 10), Windows (MSVC 2017) and Linux (gcc >= 7, clang >= 4).

For image loading and saving, Yocto/GL depends on `stb_image.h`,
`stb_image_write.h`, `stb_image_resize.h` and `tinyexr.h`.
To support Khronos glTF, Yocto/GL depends on `json.hpp`. 
All dependencies are included in the distribution.

OpenGL utilities include the OpenGL libraries, use GLEW on Windows/Linux,
GLFW for windows handling and Dear ImGui for UI support.
Since OpenGL is quite onerous and hard to link, its support can be disabled
by defining YGL_OPENGL to 1 before including this file. If you use any of
the OpenGL calls, make sure to properly link to the OpenGL libraries on
your system. OpenGL extensions use `glad.{h, cpp}` For ImGUI, build with the 
libraries `imgui.cpp`, `imgui_draw.cpp`, `imgui_impl_glfw_gl3.cpp`.
For raytracing, we optionally link to Intel's Embree if `YGL_EMBREE` is 
defined at build time.


## Design Considerations

Yocto/GL follows a "data-driven programming model" that makes data explicit.
Data is stored in simple structs and access with free functions or directly.
All data is public, so we make no attempt at encapsulation.
All objects is Yocto/GL have value semantic and we do not use pointers
in data structure but indices. This means that everything can be trivially
serialized and there is no need for memory management.

In terms of code style we prefer a functional approach rather than an
object oriented one, favoring free functions to class methods. All functions
and data are defined in the `yocto` namespace so any library can call all
others. We do this to make it as easy as possible to extend the library simply
by extending the `yocto` namespace.

We do this since this makes Yocto/GL easier to extend and quicker to learn,
with a more explicit data flow that is easier whrn writing parallel code.
Since Yocto/GL is mainly used for research and teaching,
explicit data is both more hackable and easier to understand.

The use of templates in Yocto was the reason for many refactorings, going
from no template to heavy template use. After many changes, we settled
on using templates just like the C++ STL. We do this to avoid code duplication
and since in some cases templates are the natural way of modeling some types.

We use exception for error reporting to reduce code size and make it seimpler to 
write mroe robust io code. This follows the stardard practice in the C++ STL.

Finally, we import math symbols from the standard library rather than
using the `std::name` pattern into the `yocto` namespace. This makes math code 
easier to read, and allows us to override functions or types if desired.
