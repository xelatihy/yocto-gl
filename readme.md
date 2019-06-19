[Image rendered with Yocto/GL path tracer](images/island.png)

# Yocto/GL: Tiny C++ Libraries for Data-Driven Physically-based Graphics

[![Build Status](https://travis-ci.org/xelatihy/yocto-gl.svg?branch=master)](https://travis-ci.org/xelatihy/yocto-gl) [![Build status](https://ci.appveyor.com/api/projects/status/rkqw7a8cenl877m6/branch/master?svg=true)](https://ci.appveyor.com/project/xelatihy/yocto-gl/branch/master)

Yocto/GL is a collection of small C++17 libraries for building 
physically-based graphics algorithms released under the MIT license.
Yocto/GL is written in a deliberatly data-driven style for ease of
development and use.
Yocto/GL is split into two small libraries to make code navigation easier.
See each header file for documentation.

- `yocto/yocto_math.{h}`: fixed-size vectors, matrices, rigid frames, rays, 
   bounding boxes, transforms, timer
- `yocto/yocto_random.{h}`: random number generation, Perlin noise, Monte Carlo
   integration utilities
- `yocto/yocto_shape.{h,cpp}`:  various utilities for manipulating 
   triangle meshes, quads meshes and line sets, computation of normals and 
   tangents, linear and Catmull-Clark subdivision, mesh loading and saving, 
   procedural shapes generation, geometry utilities 
- `yocto/yocto_bvh.{h,cpp}`: ray intersection and closest point queries of 
   triangle meshes, quads meshes, line sets and instances scenes using a 
   two-level bounding volume hierarchy
- `yocto/yocto_image.{h,cpp}`: simple image data type, image resizing, 
   tonemapping, color correction, image loading and saving, 
   procedural images, procedural sun-sky, color conversion utilities
- `yocto/yocto_trace.{h,cpp}`: path tracing of surfaces and hairs supporting
   area and environment illumination, microfacet GGX and subsurface scattering,
   multiple importance sampling
- `yocto/yocto_scene.{h,cpp}`: simple scene storage, evaluation of scene 
   properties
- `yocto/yocto_obj.{h,cpp}`: OBJ parser based on callbacks (SAX-like)
- `yocto/yocto_pbrt.{h,cpp}`: pbrt parser based on callbacks (SAX-like)
- `yocto/yocto_sceneio.{h,cpp}`: scene loading and saving of OBJ, pbrt, glTF,
   and a custom, hackable, YAML format

You can see Yocto/GL in action in the following applications written to
test the library:

- `apps/yscntrace.cpp`: command-line path-tracer
- `apps/yscnitrace.cpp`: interactive path-tracer
- `apps/yscnproc.cpp`: command-line scene manipulation and conversion
- `apps/yimgview.cpp`: HDR/PNG/JPG image viewer with tonemapping and color grading
- `apps/yimgproc.cpp`: command-line image manipulation
- `apps/yscnview.cpp`: simple OpenGL viewer

Here are some test images rendered with the path tracer. More images are 
included in the [gallery](gallery.md).

![Example materials: matte, plastic, metal, glass, subsurface, normal mapping](images/features1.png)

![Example shapes: procedural shapes, Catmull-Clark subdivision, hairs, displacement mapping](images/features2.png)

## Design Considerations

Yocto/GL follows a "data-driven programming model" that makes data explicit.
Data is stored in simple structs and accessed with free functions or directly.
All data is public, so we make no attempt at encapsulation.
All objects is Yocto/GL have value semantic and we do not use pointers
in data structure but indices. This means that everything can be trivially
serialized and there is no need for memory management.

We do this since this makes Yocto/GL easier to extend and quicker to learn,
with a more explicit data flow that is easier when writing parallel code.
Since Yocto/GL is mainly used for research and teaching,
explicit data is both more hackable and easier to understand.

In terms of code style we prefer a functional approach rather than an
object oriented one, favoring free functions to class methods. All functions
and data are defined in the `yocto` namespace so any library can call all
others. We do this to make it as easy as possible to extend the library simply
by extending the `yocto` namespace.

The use of templates in Yocto was the reason for many refactorings, going
from no template to heavy template use. After many changes, we settled
on using few templates for readability.

We use exception for error reporting to reduce code size and make it seimpler to 
write and more robust io code. This follows the stardard practice in the C++ STL.

## Credits

Main contributors:
  - Fabio Pellacini (lead developer): [web](http://pellacini.di.uniroma1.it), [github](https://github.com/xelatihy) 
  - Edoardo Carra: [linkedin](https://www.linkedin.com/in/edoardocarra/)
  - Giacomo Nazzaro: [github](https://github.com/giacomonazzaro)

This library includes code from the [PCG random number generator](http://www.pcg-random.org),
boost `hash_combine`, and public domain code from `github.com/sgorsten/linalg`, 
`gist.github.com/badboy/6267743` and `github.com/nothings/stb_perlin.h`.
Other external libraries are included with their own license.

## Compilation

This library requires a C++17 compiler and is know to compiled on 
OsX (Xcode >= 10), Windows (MSVC 2017) and Linux (gcc >= 7, clang >= 4).

You can build the example applications using CMake with
    `mkdir build; cd build; cmake ..; cmake --build`

Yocto/GL depends on `stb_image.h`, `stb_image_write.h`, `stb_image_resize.h` and
`tinyexr.h` for image loading, saving and resizing,  `happly.hpp` and `cgltf.h` 
for PLY and glTF support, and `filesystem.hpp` to support C++17 filesystem API 
when missing. All dependencies are included in the distribution.

OpenGL utilities include the OpenGL libraries, use GLEW on Windows/Linux,
GLFW for windows handling and Dear ImGui for UI support.
Since OpenGL is quite onerous and hard to link, its support can be disabled
by defining YGL_OPENGL to 1 before including this file. If you use any of
the OpenGL calls, make sure to properly link to the OpenGL libraries on
your system. OpenGL extensions use `glad.{h, cpp}` For ImGUI, build with the 
libraries `imgui.cpp`, `imgui_draw.cpp`, `imgui_impl_glfw_gl3.cpp`.
For raytracing, we optionally link to Intel's Embree if `YGL_EMBREE` is 
defined at build time.
