# Yocto/GL: Tiny C++ Libraries for Physically-based Graphics

Yocto/GL is a collection of utility C++14 libraries for building 
physically-based graphics algorithms released under the MIT license.
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

Yocto/GL is written in C++14 and compiles on OSX (clang from Xcode 9+),
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

- `yocto/yocto_gl.{h,cpp}`: math library, image and shape operations, ray-scene 
  intesection with multi-level BVH, simple scene format, path tracing
- `yocto/yocto_glio.{h,cpp}`: file and path operations, image loading and saving, 
  scene loading and saving, immediate-mode command line parsing


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

This library requires a C++14 compiler and is know to compiled on 
OsX (Xcode >= 8), Windows (MSVC 2017) and Linux (gcc >= 6, clang >= 4).

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

Yocto/GL tries to follow a "data-driven programming model" that makes data
explicit. Data is stored in simple structs and access with free functions
or directly. All data is public, so we make no attempt at encapsulation.
We do this since this makes Yocto/GL easier to extend and quicker to learn,
which a more explicit data flow that is easier whrn writing parallel code.
Since Yocto/GL is mainly used for research and teaching,
explicit data is both more hackable and easier to understand.

The use of templates in Yocto was the reason for many refactorings, going
from no template to heavy template use. After many changes, we settled
on using templates following the established convention in the C++ standard 
library.

We make use of exception for error reporting. This makes the code
cleaner and more in line with the expectation of most other programming
languages. At the same time, exception are not as easy to use in C++
and are disabled in many libraries. For this reasons, this will likely 
change in the future.

Finally, we import math symbols from the standard library rather than
using the `std::name` pattern into the `ygl` namespace. This makes math code 
easier to read, and allows us to override come function implementation when
desired.
