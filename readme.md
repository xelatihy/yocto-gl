# Yocto/GL: Tiny C++ Libraries for Physically-based Graphics

Yocto/GL is a collection of utility libraries for building physically-based 
graphics algorithms implemented released under the MIT license.
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
- simple logger
- path tracer supporting surfaces and hairs, GGX and MIS
- support for loading and saving Wavefront OBJ and Khronos glTF
- OpenGL utilities to manage textures, buffers and prograrms
- OpenGL shader for image viewing and GGX microfacet and hair rendering

The current version is 0.5.0.

Yocto/GL is written in C++14 and compiles on OSX (clang from Xcode 9+),
Linux (gcc 6+, clang 4+) and Windows (MSVC 2015, MSVC 2017). For compilation
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

Yocto/GL is split into many small libraries to make code navigation easier.

- `yocto_math.h`: fixed-size vectors, matrices, frames, transforms, bounding
  boxes, rays
- `yocto_image.{h,cpp}`: image container, image loading and saving, 
  tone mapping, procedural sun-sky, image examples, image utlities
- `yocto_shape.{h,cpp}`: 3D shapes utilities, shape tesselation, Catmull-Clark 
   subdivision surfaces, primitive shape generation, procedural hair
- `yocto_bvh.{h,cpp}`: ray-scene intersection and closest point queries
  accelerated by a two-level bounding volume hierarchy
- `yocto_scene.{h,cpp}`: simple scene data structure for demos and path tracing
- `yocto_trace.{h,cpp}`: path tracing with GGX material, hair and point shading,
  mesh lights, environment illumination, all with MIS sampling
- `yocto_obj.{h,cpp}`: Wavefront OBJ loading and saving
- `yocto_gltf.{h,cpp}`: Khronos glTF loading and saving
- `yocto_glutils.{h,cpp}`: OpenGL utilities, GLFW and dear ImGui UIs wrapper
- `yocto_utils.h`: utilities for command line applications


## Example Applications

You can see Yocto/GL in action in the following applications written to
test the library:

- `yview.cpp`: simple OpenGL viewer for OBJ and glTF scenes
- `ytrace.cpp`: offline path-tracer
- `yitrace.cpp`: interactive path-tracer
- `yscnproc.cpp`: scene manipulation and conversion to/from OBJ and glTF
- `ytestgen.cpp`: creates test cases for the path tracer and GL viewer
- `yimview.cpp`: HDR/PNG/JPG image viewer with exposure/gamma tone mapping
- `yimproc.cpp`: offline image manipulation.

You can build the example applications using CMake with
    `mkdir build; cd build; cmake ..; cmake --build`


## Compilation

For image loading and saving, Yocto/GL depends on `stb_image.h`,
`stb_image_write.h`, `stb_image_resize.h` and `tinyexr.h`. These features
can be disabled by defining YGL_IMAGEIO to 0 before including this file.
If these features are useful, then the implementation files need to
included in the manner described by the respective libraries. To simplify
builds, we provide a file that builds these libraries, `stb_image.cpp`.

To support Khronos glTF, Yocto/GL depends on `json.hpp`. This feature can
be disabled by defining YGL_GLTF to 0 before including this file.

To support SVG, Yocto/GL depends on `nanosvg.h`. This feature can
be disabled by defining YGL_SVG to 0 before including this file.

OpenGL utilities include the OpenGL libraries, use GLEW on Windows/Linux,
GLFW for windows handling and Dear ImGui for UI support.
Since OpenGL is quite onerous and hard to link, its support can be disabled
by defining YGL_OPENGL to 1 before including this file. If you use any of
the OpenGL calls, make sure to properly link to the OpenGL libraries on
your system. For ImGUI, build with the libraries `imgui.cpp`,
`imgui_draw.cpp`, `imgui_impl_glfw_gl3.cpp`.


## Design Considerations

Yocto/GL tries to follow a simple programming model inspired by C but with
heavy use of operator overloading for math readability. We attempt to make
the code easy to use use rather than as performant as possible.
We adopt a functional style and only rarely use classes and methods.
Using a function style makes the code easier to extend, more explicit in
the function requirements, and easier to write parallel-friendly APIs.
I guess you could call this "data-driven programming".

The use of templates in Yocto was the reason for many refactorings, going
from no template to heavy templates use. At this time, templates are used
in the basic types to make the codebase shorter and reduce bugs,
at the price of accessibility for beginners. The truth is that modern C++,
a tenant of Yocto, is heavily templated anyway, so being able to read
template code is necessary no matter how Yocto does things.

We make use of exception for error reporting. This makes the code
much cleaner and more in line with the expectation of most other programming
languages.

Finally, we often import symbols from the standard library rather than
using the `std::name` pattern. We found that this improves consistency
when using math functions, and is more readable with templates. We realize
this is not standard, but the imports are hidden within the ygl namespace,
so library users do not have to be concern about it.
