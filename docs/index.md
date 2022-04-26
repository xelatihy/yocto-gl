# Yocto/GL: Tiny C++ Libraries for Data-Oriented Physically-based Graphics

Yocto/GL is a collection of small C++17 libraries for building
physically-based graphics algorithms released under the MIT license.
Yocto/GL is written in a deliberately data-oriented style for ease of
development and use.

## Libraries

Yocto/GL is split into small libraries to make code navigation easier.
Here is a list of the main libraries.

- [Yocto/Math](yocto/yocto_math.md): fixed-size vectors, matrices, rigid frames,
  transforms
- [Yocto/Color](yocto/yocto_color.md): color conversion, color adjustment,
  tone mapping functions, color grading, color maps, color spaces
- [Yocto/Geometry](yocto/yocto_geometry.md): rays, bounding boxes,
  geometry functions, ray-primitive intersection, point-primitive overlap
- [Yocto/Noise](yocto/yocto_noise.md): Perlin noise
- [Yocto/Sampling](yocto/yocto_sampling.md): random number generation,
  generation of points and directions, Monte Carlo utilities
- [Yocto/Shading](yocto/yocto_shading.md): evaluation and sampling of fresnel
  functions, bsdf lobes, transmittance lobes, phase functions
- [Yocto/Image](yocto/yocto_image.md): simple image data type, image resizing,
  tonemapping, color correction, procedural images, procedural sun-sky
- [Yocto/Shape](yocto/yocto_shape.md): simple shape data structure, utilities 
  for manipulating triangle meshes, quads meshes and line sets, computation of 
  normals and tangents, linear and Catmull-Clark subdivision, 
  procedural shapes generation, ray intersection and closest point queries
- [Yocto/Scene](yocto/yocto_scene.md): scene representation and properties
  evaluation
- [Yocto/Bvh](yocto/yocto_bvh.md): ray intersection and closest point queries
  of triangle meshes, quads meshes, line sets and shape instances using a
  two-level bounding volume hierarchy
- [Yocto/Trace](yocto/yocto_trace.md): path tracing of surfaces and hairs
  supporting area and environment illumination, microfacet GGX and subsurface
  scattering, multiple importance sampling
- [Yocto/SceneIO](yocto/yocto_sceneio.md): image, shape and scene serialization
- [Yocto/ModelIO](yocto/yocto_modelio.md): low-level parsing and writing for
  Ply, Obj, Stl formats
- [Yocto/PbrtIO](yocto/yocto_pbrtio.md): low-level parsing and writing for
  Pbrt format
- [Yocto/Cli](yocto/yocto_cli.md): printing utilities and command line parsing
- [Yocto/Parallel](yocto/yocto_parallel.md): concurrency utilities (deprecated)

## Example Applications

You can see Yocto/GL in action in the following applications written to
test the library:

- `apps/ytonemap.cpp`: image conversion and viewing
- `apps/ycolorgrade.cpp`: image color grading
- `apps/yconvert.cpp`: scene conversion
- `apps/yconverts.cpp`: shape conversion
- `apps/ytrace.cpp`: offline and interactive scene rendering
- `apps/ycutrace.cpp`: offline and interactive scene rendering with CUDA
- `apps/yview.cpp`: interactive scene viewing

Here are some test images rendered with the path tracer. More images are
included in the [project site](https://xelatihy.github.io/yocto-gl/).

![Example materials: matte, plastic, metal, glass, subsurface, normal mapping](images/features1.jpg)

![Example shapes: procedural shapes, Catmull-Clark subdivision, hairs, displacement mapping](images/features2.jpg)

![Image rendered with Yocto/GL path tracer. Model by Disney Animation Studios.](images/island.jpg)

## Design Considerations

Yocto/GL follows a "data-oriented programming model" that makes data explicit.
Data is stored in simple structs and accessed with free functions or directly.
All data is public, so we make no attempt at encapsulation.
We do this since this makes Yocto/GL easier to extend and quicker to learn,
with a more explicit data flow that is easier when writing parallel code.
Since Yocto/GL is mainly used for research and teaching,
explicit data is both more hackable and easier to understand.

Nearly all objects in Yocto/GL have value semantic. This means that everything
can be trivially copied and serialized and there is no need for memory management. 
While this has the drawback of potentially introducing spurious copies, it does
have the benefit of ensuring that no memory corruption can occur, which
turned out was a major issue for novice C++ users, even in a very small
library like this one. 

In terms of code style we prefer a functional approach rather than an
object oriented one, favoring free functions to class methods. All functions
and data are defined in the `yocto` namespace so libraries can call each others
trivially.

The use of templates in Yocto was the reason for many refactoring, going
from no template to heavy template use. At this point, Yocto uses some templates
for readability. In the future, we will increase the use of templates in math
code, while keeping many APIs explicitly typed.

For error handling in IO we either return status object or an interface that
uses boolean flags and error strings. Internally exceptions are used when used
by external libraries, but otherwise no exception are used. At the moment,
exceptions are only used to report "programmer errors", namely when 
preconditions or post conditions are violated in functions, just lime the
standard library does.

## License

The library is released under the MIT license. We include various external 
dependencies in the distribution that each have thir own license, compatible
with the chosen one.

## Compilation

This library requires a C++17 compiler and is know to compiled on
OsX (Xcode >= 11), Windows (MSVC 2019) and Linux (gcc >= 9, clang >= 9).

You can build the example applications using CMake with
`mkdir build; cd build; cmake ..; cmake --build .`

Yocto/GL required dependencies are included in the distribution and do not
need to be installed separately.

Yocto/GL optionally supports building OpenGL demos. OpenGL support is enabled 
by defining the cmake option `YOCTO_OPENGL` and contained in the `yocto_gui` 
library. OpenGL dependencies are included in this repo.

Yocto/GL optionally supports the use of Intel's Embree for ray casting.
See the main CMake file for how to link to it. Embree support is enabled by
defining the cmake option `YOCTO_EMBREE`. Embree needs to be installed separately.

Yocto/GL optionally supports the use of Intel's Open Image Denoise for denoising.
See the main CMake file for how to link to it. Open Image Denoise support
is enabled by defining the cmake option `YOCTO_DENOISE`. 
OIDN needs to be installed separately.
