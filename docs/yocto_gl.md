# Yocto/GL: Tiny C++ Library for Physically-based Graphics

Yocto/GL is a collection utiliies for building physically-based graphics
algorithms implemented as a two-file library (`yocto_gl.h`, `yocto_gl.cpp`),
and released under the MIT license. Features include:

- convenience math functions for graphics
- static length vectors for 2, 3, 4 length and int and float type
- static length matrices for 2x2, 3x3, 4x4 and float type
- static length rigid transforms (frames), specialized for 2d and 3d space
- linear algebra operations and transforms
- axis aligned bounding boxes
- rays and ray-primitive intersection
- point-primitive distance and overlap tests
- normal and tangent computation for meshes and lines
- generation of tesselated meshes
- mesh refinement with linear tesselation and Catmull-Cark subdivision
- random number generation via PCG32
- simple image data structure and a few image operations
- simple scene format
- generation of image examples
- generation of scene examples
- procedural sun and sky HDR
- procedural Perlin noise
- BVH for intersection and closest point query
- Python-like iterators, string, path and container operations
- utilities to load and save entire text and binary files
- immediate mode command line parser
- simple logger and thread pool
- path tracer supporting surfaces and hairs, GGX and MIS
- support for loading and saving Wavefront OBJ and Khronos glTF
- OpenGL utilities to manage textures, buffers and prograrms
- OpenGL shader for image viewing and GGX microfacet and hair rendering

The current version is 0.1.0. You can access the previous multi-file version
with tag "v0.0.1" in this repository.

## Credits

This library includes code from the PCG random number generator,
the LLVM thread pool, boost hash_combine, Pixar multijittered sampling,
code from "Real-Time Collision Detection" by Christer Ericson, base64
encode/decode by René Nyffenegger and public domain code from
github.com/sgorsten/linalg, gist.github.com/badboy/6267743 and
github.com/nothings/stb_perlin.h.

This library imports many symbols from std for three reasons: avoid
verbosity , esnuring better conventions when calling math functions and
allowing easy overriding of std containers if desired. Just do not
flatten this namespace into yours if this is a concern.

For most components of the library, the use should be relatively easy to
understand if you are familiar with 3d computer graphics. For more complex
components, we follow the usage below.


## Design Considerations

Yocto/GL tries to follow a simple programming model inspired by C but with
heavy use of operator overloading for math readability. We attempt tp make
the code weasy to use rather than as performant as possible. The APIs
attempt to make using the code as little error prone as possible, sometimes
at the price of some slowdown. We adopt a functional style and only rarely
use classes and methods. Using a function style makes the code easier to
extend, more explicit in the requirements, and easier to write
parallel-friendly APIs. I guess you could call this "data-driven
programming". We use templates very little now, after a major refactoring,
to improve error reporting, reduce compilation times and make the codebase
more accessible to beginners. This lead to a small increase in copied code
that we deem ok at this time. Finally, we often import symbols from the
standard library rather than using the `std::name` pattern. We found that
this improves consistency, especially when using math functions, is
significantly more readable when using templates and allows to to more
easily switch STL implementation if desired.


## Compilation

Yocto/GL is written in C++14, with compilation supported on C++11, and
compiles on OSX (clang from Xcode 9+), Linux (gcc 6+, clang 4+)
and Windows (MSVC 2017).

For image loading and saving, Yocto/GL depends on `stb_image.h`,
`stb_image_write.h`, `stb_image_resize.h` and `tinyexr.h`. These features
can be disabled by defining YGL_IMAGEIO to 0 before including this file.
If these features are useful, then the implementation files need to
included in the manner described by the respective libraries. To simplify
builds, we provice a file that builds these libraries, `stb_image.cpp`.

To support Khronos glTF, Yocto/GL depends on `json.hpp`. These feature can
be disabled by defining YGL_GLTF to 0 before including this file.

OpenGL utilities include the OpenGL libaries, use GLEW on Windows/Linux,
GLFW for windows handling and Dear ImGui for UI support.
Since OpenGL is quite onerous and hard to link, its support is disabled by
default. You can enable it by defining YGL_OPENGL to 1 before including
this file. If you use any of the OpenGL calls, make sure to properly link to
the OpenGL libraries on your system. For ImGUI, build with the libraries
`imgui.cpp`, `imgui_draw.cpp`, `imgui_impl_glfw_gl3.cpp`.


## Example Applications

You can see Yocto/GL in action in the following applications written to
test the library:

- `yview.cpp`: simple OpenGL viewer for OBJ and glTF scenes
- `ytrace.cpp`: offline path-tracer
- `yitrace.cpp.cpp`: interactive path-tracer
- `yscnproc.cpp`: scene manipulation and conversion to/from OBJ and glTF
- `ytestgen.cpp`: creates test cases for the path tracer and GL viewer
- `yimview.cpp`: HDR/PNG/JPG image viewer with exposure/gamma tone mapping
- `yimproc.cpp`: offline image manipulation.

You can build the example applications using CMake with
    `mkdir build; cd build; cmake ..; cmake --build`

Here are two images rendered with the buildin path tracer, where the
scenes are crated with the test generator.

![Yocto/GL](images/shapes.png)

![Yocto/GL](images/lines.png)


## Usage

To use the library simply include this file and setup the compilation
option as described above.
All library features are documented at the definition and should be
relatively easy to use if you are familiar with writing graphics code.
You can find the extracted documentation at `yocto_gl.md`.
Here we give an overview of some of the main features.


### Small Vectors and Matrices, Frames, Bounding Boxes and Transforms

We provide common operations for small vectors and matrices typically used
in graphics. In particular, we support 2-4 dimensional float vectors
`vec2f`, `vec3f`, `vec4f`, 2-4 dimensional int vectors `vec2i`, `vec3i`,
`vec4i` and a 4 dimensional byte vector `vec4b`. The float vectors
support most arithmetic and vector operations.

We support 2-4 dimensional float matrices `mat2f`, `mat3f`, `mat4f`, with
matrix-matrix and matrix-vector products, trasposes and inverses. Matrices
are stored in column-major ordered and are accessed and constructed by
column.

To represent transformations, most of the library facilities prefer the use
cooordinate frames, aka rigid transforms, represented as `frame3f`.
The structure store three coodinate axis and the frame origin. This is
equivalenent to a rigid transform written as a column-major affine
matrix. Transform operations are better behaved with this representation.

We represent coordinate bounds with axis-aligned bounding boxes in 1-4
dimensions: `bbox1f`, `bbox2f`, `bbox3f`, `bbox4f`. These types support
expansion operation, union and containment. We provide operations to
compute bounds for points, lines, triangles and quads.

For all basic types we support iteration with `begin()`/`end()` pairs
and stream inout and output.

For both matrices and frames we support transform operations for points,
vectors and directions (`trasform_point()`, `trasform_vector()`,
`trasform_direction()`). For frames we also the support inverse operations
(`transform_xxx_inverse()`). Transform matrices and frames can be
constructed from basic translation, rotation and scaling, e.g. with
`translation_mat4f()` or `translation_frame3f()` repsectively, etc. For
rotation we support axis-angle and quaternions, with slerp.


### Random Number Generation, Noise, Hashing and Monte Carlo support

This library supportds many facitlities helpful in writing sampling
functions targeting path tracing and shape generations.

1. Random number generation with PCG32:
    1. initialize the random number generator with `init_rng()`
    2. advance the random number state with `advance_rng()`
    3. if necessary, you can reseed the rng with `seed_rng()`
    4. generate random integers in an interval with `next_rand1i()`
    5. generate random floats and double in the [0,1) range with
       `next_rand1f()`, `next_rand2f()`, `next_rand3f()`, `next_rand1d()`
    6. you can skip random numbers with `advance_rng()` and get the skipped
       length with `rng_distance()`
    7. generate random shaffled sequences with `rng_shuffle()`
2. Perlin noise: `perlin_noise()` to generate Perlin noise with optional
   wrapping, with fractal variations `perlin_ridge_noise()`,
   `perlin_fbm_noise()`, `perlin_turbulence_noise()`
3. Integer hashing: public domain hash functions for integer values as
   `hash_permute()`, `hash_uint32()`, `hash_uint64()`, `hash_uint64_32()`
   and `hash_combine()`.
4. Monte Carlo support: warp functions from [0,1)^k domains to domains
   commonly used in path tracing. In particular, use `sample_hemisphere()`,
   `sample_sphere()`, `sample_hemisphere_cosine()`,
   `sample_hemisphere_cospower()`. `sample_disk()`. `sample_cylinder()`.
   `sample_triangle()`. For each warp, you can compute the PDF with
   `sample_xxx_pdf()`.


### Python-like container operations and iterators

To make the code more readable, we adopt Python-like iterations and
container operations extensively throughout Yocto/GL. These operations
are mostly for internal use but could also be used externally.

1. Python iterators with `range()` and `enumerate()`
2. Python operators for containers: support for + and += for `std::vector`
3. Check for containment with `contains`  similarly to `in` in Python


### Shape Utilities

The library contains a few function to help with typically geometry
manipulation useful to support scene viewing and path tracing.

1. compute line tangents, and triangle and quad areas and normals
2. compute barycentric interpolation with `eval_barycentric_line()`,
   `eval_barycentric_triangle()` and `eval_barycentric_quad()`
3. evaluate Bezier curve and derivatives with `eval_bezier_cubic()` and
   `eval_bezier_cubic_derivative()`
4. compute smooth normals and tangents with `compute_normals()`
5. compute tangent frames from texture coordinates with
   `compute_tangent_space()`
6. compute skinning with `compute_skinning()` and
   `compute_matrix_skinning()`
6. shape creation with `make_points()`, `make_lines()`, `make_uvgrid()`
7. element merging with `marge_elems()`
8. facet elements with `facet_elems()`
9. shape sampling with `sample_points()`, `sample_lines()`,
   `sample_triangles()`; initialize the sampling CDFs with
   `sample_points_cdf()`, `sample_lines_cdf()`, `sample_triangles_cdf()`
10. samnple a could of point over a surface with `sample_triangles_points()`
11. get edges and boundaries with `get_edges()` and `get_boundary_edges()`
12. convert quads to triangles with `convert_quads_to_triangles()`
13. convert face varying to vertex shared representations with
    `convert_face_varying()`
14. subdivide elements by edge splits with `subdivide_elems()` and
    `subdivide_vert()`
15. Catmull-Clark subdivision surface with `subdivide_catmullclark()` with
    support for edge and vertex creasing
16. example shapes: `make_cube()`, `make_uvsphere()`, `make_uvhemisphere()`,
    `make_uvquad()`, `make_uvcube()`, `make_fvcube()`, `make_hair()`,
    `make_suzanne()`


### Image and color

We support simple containers for either 4-byte per pixel sRGB images
`image4b`, or 4-float per pixel HDR images `image4f`.

1. convert between byte and float images with `srgb_to_linear()` and
   `linear_to_srgb()`
2. color conversion with `hsv_to_rgb()`, `xyz_to_rgb()` and `rgb_to_xyz()`
3. exposure-gamma tonemapping, with optional filmic curve, with
   `tonemap_image()`
4. compositing support with `image_over()`
5. example image generation with `m,ake_grid_image()`,
   `make_checker_image()`, `make_bumpdimple_image()`, `make_ramp_image()`,
   `make_gammaramp_image()`, `make_gammaramp_imagef()`, `make_uv_image()`,
   `make_uvgrid_image()`, `make_recuvgrid_image()`
6. bump to normal mapping with `bump_to_normal_map()`
7. HDR sun-sky with `m ake_sunsky_image()`
8. various noise images with `make_noise_image()`, `make_fbm_image()`,
   `make_ridge_image()`, `make_turbulence_image()`
9. image loading and saving with `load_image4b()`, `load_image4f()`,
   `save_image4b()`, `save_image4f()`
10. image resizing with `resize_image()`


### Ray Intersection and Point Overlap Queries

We support ray-scene intersection for points, lines and triangles
accelerated by a simple BVH data structure.  Our BVH is written for minimal
code and not maximum speed, but still gives reasonable results. We suggest
the use of Intel's Embree as a fast alternative.

1. use `ray3f` to represent rays
2. build the BVH with `build_points_bvh()`, `build_points_bvh()` or
  `build_points_bvh()`
3. perform ray-element intersection with `intersect_points_bvh()`,
  `intersect_lines_bvh()` and `intersect_triangles_bvh()`
4. perform point overlap queries with `overlap_points_bvh()`,
  `overlap_lines_bvh()` and `overlap_triangles_bvh()`
5. to support custom elements, use `buid_bvh()`, `intersect_bvh()` and
  `overlap_bvh()` and provide them with proper callbacks
6. we also experimentally support quads with the `xxx_quads_xxx()` functions


### Simple scene

We support a simple scene model used to quickly write demos that lets you
load/save Wavefront OBJ and Khronos glTF and perform several simple scene
manipulation including ray-scene intersection and closest point queries.

The geometry model is comprised of a set of shapes, which are indexed
collections of points, lines, triangles and quads. Each shape may contain
only one element type. Shapes are organized into a scene by creating shape
instances, each its own transform. Materials are specified like in glTF and
include emission, base-metallic and diffuse-specular parametrization,
normal, occlusion and displacement mapping. Finally, the scene containes
caemras and environement maps. Quad support in shapes is experimental and
mostly supported for loading and saving.

For low-level access to OBJ/glTF formats, you are best accssing the formats
directly with Yocto/Obj and Yocto/glTF. This components provides a
simplified high-level access to each format which is sufficient for most
applications and tuned for quick creating viewers, renderers and simulators.

1. load a scene with `load_scene()` and save it with `save_scene()`.
2. add missing data with `add_elements()`
3. use `compute_bounds()` to compute element bounds
4. can merge scene together with `merge_into()`
5. make example scenes with `make_test_scene()`

Ray-intersection and closet-point routines supporting points,
lines and triangles accelerated by a two-level bounding volume
hierarchy (BVH). Quad support is experimental.

1. build the bvh with `build_bvh()`
2. perform ray-interseciton tests with `intersect_ray()`
    - use early_exit=false if you want to know the closest hit point
    - use early_exit=false if you only need to know whether there is a hit
    - for points and lines, a radius is required
    - for triangles, the radius is ignored
2. perform point overlap tests with `overlap_point()` to check whether
   a point overlaps with an element within a maximum distance
    - use early_exit as above
    - for all primitives, a radius is used if defined, but should
      be very small compared to the size of the primitive since the radius
      overlap is approximate
3. perform instance overlap queries with `overlap_instance_bounds()`
4. use `refit_bvh()` to recompute the bvh bounds if transforms or vertices
   are changed (you should rebuild the bvh for large changes)

Notes: Quads are internally handled as a pair of two triangles v0,v1,v3 and
v2,v3,v1, with the u/v coordinates of the second triangle corrected as 1-u
and 1-v to produce a quad parametrization where u and v go from 0 to 1. This
is equivalent to Intel's Embree.


### Pathtracing

We supply a path tracer implementation with support for textured mesh
lights, GGX/Phong materials, environment mapping. The interface supports
progressive parallel execution. The path tracer takes as input a scene
and update pixels in image with traced samples. We use a straightfoward
path tracer with MIS and also a few simpler shaders for debugging or
quick image generation.

Materials are represented as sums of an emission term, a diffuse term and
a specular microfacet term (GGX or Phong). Only opaque for now. We pick
a proper material type for each shape element type (points, lines,
triangles).

Lights are defined as any shape with a material emission term. Additionally
one can also add environment maps. But even if you can, you might want to
add a large triangle mesh with inward normals instead. The latter is more
general (you can even more an arbitrary shape sun). For now only the first
env is used.

1. build the ray-tracing acceleration structure with `build_bvh()`
2. prepare lights for rendering `update_lights()`
3. define rendering params with the `trace_params` structure
4. render blocks of samples with `trace_block()`

The code can also run in fully asynchronous mode to preview images in a
window.

1. build the ray-tracing acceleration structure with `build_bvh()`
2. prepare lights for rendering `update_lights()`
3. define rendering params with the `trace_params` structure
4. initialize the prograssive rendering buffers
5. start the progressive renderer with `trace_async_start()`
7. stop the progressive renderer with `trace_async_stop()`


### Wavefront OBJ

Wavefront OBJ/MTL loader and writer with support for points,
lines, triangles and general polygons and all materials properties.
Contains also a few extensions to easily create demos such as per-vertex
color and radius, cameras, environment maps and instances.
Can use either a low-level OBJ representation, from this files,
or a high level flattened representation included in Yocto/Scn.

Both in reading and writing, OBJ has no clear convention on the orientation
of textures Y axis. So in many cases textures appears flipped. To handle
that, use the option to flip textures coordinates on either saving or
loading. By default texture coordinates are flipped since this seems
the convention found on test cases collected on the web. The value Tr
has similar problems, since its relation to opacity is software specific.
Again we let the user chose the convension and set the default to the
one found on the web.

In the high level interface, shapes are indexed meshes and are described
by arrays of vertex indices for points/lines/triangles and arrays for vertex
positions, normals, texcoords, color and radius. The latter two as
extensions. Since OBJ is a complex formats that does not match well with
current GPU rendering / path tracing algorithms, we adopt a simplification
similar to other single file libraries:
1. vertex indices are unique, as in OpenGL and al standard indexed triangle
  meshes data structures, and not OBJ triplets; YOCTO_OBJ ensures that no
  vertex dusplication happens thought for same triplets
2. we split shapes on changes to groups and materials, instead of keeping
  per-face group/material data; this makes the data usable right away in
  a GPU viewer; this is not a major limitation if we accept the previous
  point that already changes shapes topology.

1. load a obj data with `load_obj()`; can load also textues
2. look at the `obj_XXX` data structures for access to individual elements
3. use obj back to disk with `save_obj()`; can also save textures
4. use get_shape() to get a flattened shape version that contains only
   triangles, lines or points


### Khronos glTF

Khronos GLTF loader and writer for Khronos glTF format. Supports
all the glTF spec and the Khronos extensions. All parsing and writing code
is autogenerated form the schema. Supports glTF version 2.0 and the
following extensions: `KHR_binary_glTF` and `KHR_specular_glossiness`.

This component depends on `json.hpp` and, for image loading and saving,
it depends on `stb_image.h`, `stb_image_write.h`, `stb_image_resize.h` and
`tinyexr.h`. This feature can be disabled as before.

The library provides a low  level interface that is a direct
C++ translation of the glTF schemas and should be used if one wants
complete control over the fromat or an application wants to have their
own scene code added. A higher-level interface is provided by the scene
or by `yocto_gltf.h`.

glTF is a very complex file format and was designed mainly with untyped
languages in mind. We attempt to match the glTF low-level interface
to C++ as best as it can. Since the code is generated from the schema, we
follow glTF naming conventions and typing quite well. To simplify adoption
and keep the API relatively simple we use vector as arrays and use
pointers to reference to all glTF objects. While this makes it less effcient
than it might have been, glTF heavy use of optional values makes this
necessary. At the same time, we do not keep track of set/unset values
for basic types (int, float, bool) as a compromise for efficieny.

glTF uses integer indices to access objects.
While writing code ourselves we found that we add signiicant problems
since we would use an index to access the wriong type of scene objects.
For this reasons, we use an explit index `glTFid<T>` that can only access
an object of type T. Internally this is just the same old glTF index. But
this can used to access the scene data with `glTF::get<T>(index)`.

1. load a glTF model with `load_gltf()`
2. look at the `glTFXXX` data structures for access to individual elements
3. save glTF back to disk with `save_gltf()`


### OpenGL support

We include a set of utilities to draw on screen with OpenGL 3.3, manage
windows with GLFW and draw immediate-mode widgets with ImGui.

1. texture and buffer objects with `gl_texture` and `gl_buffer`
    - create textures/buffers with appropriate constructors
    - check validity wiht `is_valid()`
    - update textures/buffers with `update()` functions
    - delete textures/buffers with `clear()`
    - bind/unbind textures/buffers with `bind()`/`unbind()`
    - draw elements with `gl_buffer::draw_elems()`
2. program objects with `gl_program`
    - program creation with constructor
    - check validity wiht `is_valid()`
    - delete with `clear()`
    - uniforms with `set_program_uniform()`
    - vertex attrib with `set_program_vertattr()`
    - draw elements with `gl_buffer::draw_elems()`
3. image viewing with `gl_stdimage_program`, with support for tone mapping.
4. draw surfaces and hair with GGX/Kayjia-Kay with `gl_stdsurface_program`
    - initialize the program with constructor
    - check validity wiht `is_valid()`
    - start/end each frame with `begin_frame()`, `end_frame()`
    - define lights with `set_lights()`
    - start/end each shape with `begin_shape()`, `end_shape()`
    - define material Parameters with `set_material()`
    - define vertices with `set_vert()`
    - draw elements with `draw_elems()`
5. draw yocto scenes using the above shader
    - initialize the rendering state with `init_stdprogram_state()`
    - load/update meshes and textures with `update_stdprogram_state()`
    - setup draw params using a `gl_stdsurface_params` struct
    - draw scene with `draw_stdprogram_scene()`
6. also includes other utlities for quick OpenGL hacking
7. GLFW window with `gl_window`
    - create with constructor
    - delete with `clear()`
    - set callbacks with `set_callbacks()`
    - includes carious utiliies to query window, mouse and keyboard
8. immediate mode widgets using ImGui
    - init with `init_widget()`
    - use the various widget calls to draw the widget and handle events


### Other Utilities

We include additional utilities for writing command line applications and
manipulating files.

1. Python-like string opeations: `startswith()`, `endswith()`, `contains()`,
   `splitlines()`, `partition()`, `split()`, `splitlines()`, `strip()`,
   `rstrip()`, `lstrip()`, `join()`, `lower()`, `upper()`, `isspace()`,
   `replace()`
2. Path-like path operations: `path_dirname()`, `path_extension()`,
   `path_basename()`, `path_filename()`, `replace_path_extension()`,
   `prepend_path_extension()`, `split_path()`
3. Python-like format strings (only support for position arguments and no
   formatting commands): `format()`, `print()`
5. load/save entire files: `load_binfile()`, `load_txtfile()`,
   `save_binfile()` and `save_binfile()`
4. simple logger with support for console and file streams:
    1. create a `logger`
    2. add more streams with `add_console_stream()` or `add_file_stream()`
    3. write log messages with `log_msg()` and its variants
    4. you can also use a global default logger with the free functions
       `log_XXX()`
5. thead pool for concurrent execution (waiting the standard to catch up):
    1. either create a `thread_pool` or use the global one
    2. run tasks in parallel `parallel_for()`
    3. run tasks asynchronously `async()`
6. timer for simple access to `std::chrono`:
    1. create a `timer`
    2. start and stop the clock with `start()` and `stop()`
    3. get time with `elapsed_time()`


### Command Line Parsing

The library includes a simple command line parser that parses commands in
immediate mode, i.e. when an option is declared. The parser supports options
and unnamed arguments with generic types parsed using C++ stream. The
parser autogenerates its own documentation. This allows to write complex
command lines with a tiny amount of implementation code on both the library
and user end.

1. create a `cmdline` parser object by passing `argc, argv, name, help`
    - an option for printing help is automatically added
2. for each option, parse it calling the functions `parse_opt()`
    - options are parsed on the fly and a comprehensive help is
      automatically generated
    - supports bool (flags), int, float, double, string, enums
    - options names are "--longname" for longname and "-s" for short
    - command line format is "--longname value", "-s v" for all but flags
    - values are parsed with `iostream <<` operators
    - for general use `opt = parse_opt<type>()`
    - for boolean flags is `parse_flag()`
    - for enums use `parse_opte()`
3. for each unnamed argument, parse it calling the functions parse_arg()
    - names are only used for help
    - supports types as above
    - for general use `arg = parse_arg<type>()`
    - to parse all remaining values use `args = parse_arga<type>(...)`
4. end cmdline parsing with `check_parsing()` to check for unsued values,
   missing arguments
5. to check for error use `should_exit()` and to print the message use
   `get_message()`
6. since arguments are parsed immediately, one can easily implement
   subcommands by just branching the command line code based on a read
   argument without any need for complex syntax


## History

Here we mark only major features added to the library. Small refactorings
and bug fixes are reported here.

- v 0.1.0: initial release after refactoring

## API Documentation

#### Function Alias sqrt()

~~~ .cpp
using std::sqrt;
~~~

sqrt

#### Function Alias pow()

~~~ .cpp
using std::pow;
~~~

pow

#### Function Alias exp()

~~~ .cpp
using std::exp;
~~~

pow

#### Function Alias log()

~~~ .cpp
using std::log;
~~~

log

#### Function Alias log10()

~~~ .cpp
using std::log10;
~~~

log10

#### Function Alias sin()

~~~ .cpp
using std::sin;
~~~

sin

#### Function Alias cos()

~~~ .cpp
using std::cos;
~~~

cos

#### Function Alias tan()

~~~ .cpp
using std::tan;
~~~

tan

#### Function Alias asin()

~~~ .cpp
using std::asin;
~~~

asin

#### Function Alias acos()

~~~ .cpp
using std::acos;
~~~

acos

#### Function Alias atan()

~~~ .cpp
using std::atan;
~~~

atan

#### Function Alias atan2()

~~~ .cpp
using std::atan2;
~~~

atan2

#### Function Alias abs()

~~~ .cpp
using std::abs;
~~~

absolute value

#### Function Alias fabs()

~~~ .cpp
using std::fabs;
~~~

floating point absolute value

#### Function Alias floor()

~~~ .cpp
using std::floor;
~~~

floor

#### Function Alias ceil()

~~~ .cpp
using std::ceil;
~~~

ceil

#### Function Alias round()

~~~ .cpp
using std::round;
~~~

round

#### Function Alias isfinite()

~~~ .cpp
using std::isfinite;
~~~

isfinate

#### Function Alias string()

~~~ .cpp
using std::string;
~~~

string

#### Function Alias vector()

~~~ .cpp
using std::vector;
~~~

vector

#### Function Alias array()

~~~ .cpp
using std::array;
~~~

array

#### Function Alias map()

~~~ .cpp
using std::map;
~~~

map

#### Function Alias unordered_map()

~~~ .cpp
using std::unordered_map;
~~~

unordered map

#### Function Alias unordered_set()

~~~ .cpp
using std::unordered_set;
~~~

unordered set

#### Function Alias pair()

~~~ .cpp
using std::pair;
~~~

pair

#### Function Alias tuple()

~~~ .cpp
using std::tuple;
~~~

tuple

#### Function Alias unique_ptr()

~~~ .cpp
using std::unique_ptr;
~~~

unique pointer

#### Function Alias function()

~~~ .cpp
using std::function;
~~~

function

#### Namespace using std::string_literals;

string literals

#### Function Alias numeric_limits()

~~~ .cpp
using std::numeric_limits;
~~~

numeric limits

#### Function Alias initializer_list()

~~~ .cpp
using std::initializer_list;
~~~

initializer list

#### Function Alias ostream()

~~~ .cpp
using std::ostream;
~~~

output stream

#### Function Alias istream()

~~~ .cpp
using std::istream;
~~~

input stream

#### Function Alias stringstream()

~~~ .cpp
using std::stringstream;
~~~

string stream

#### Function Alias fstream()

~~~ .cpp
using std::fstream;
~~~

file stream

#### Function Alias runtime_error()

~~~ .cpp
using std::runtime_error;
~~~

runtime error

#### Function Alias exception()

~~~ .cpp
using std::exception;
~~~

exception

#### Function Alias ios_base()

~~~ .cpp
using std::ios_base;
~~~

ios base

#### Function Alias find()

~~~ .cpp
using std::find;
~~~

find algorithms

#### Function Alias swap()

~~~ .cpp
using std::swap;
~~~

swap algorithms

#### Function Alias getline()

~~~ .cpp
using std::getline;
~~~

get line from streams

#### Function Alias to_string()

~~~ .cpp
using std::to_string;
~~~

convert to string

#### Function Alias cout()

~~~ .cpp
using std::cout;
~~~

cout object for printing

#### Typedef byte

~~~ .cpp
using byte = unsigned char;
~~~

convenient typedef for bytes

#### Typedef uint

~~~ .cpp
using uint = unsigned int;
~~~

convenient typedef for bytes

#### Constant pif

~~~ .cpp
const auto pif = 3.14159265f;
~~~

pi (float)

#### Constant pi

~~~ .cpp
const auto pi = 3.1415926535897932384626433832795;
~~~

pi (double)

#### Constant flt_max

~~~ .cpp
const auto flt_max = numeric_limits<float>::max();
~~~

shortcat for float max value

#### Constant flt_min

~~~ .cpp
const auto flt_min = numeric_limits<float>::lowest();
~~~

shortcat for float min value

#### Constant flt_eps

~~~ .cpp
const auto flt_eps = numeric_limits<float>::epsilon();
~~~

shortcat for float epsilon

#### Constant int_max

~~~ .cpp
const auto int_max = numeric_limits<int>::max();
~~~

shortcat for int max value

#### Constant int_min

~~~ .cpp
const auto int_min = numeric_limits<int>::min();
~~~

shortcat for int min value

#### Function min()

~~~ .cpp
inline int min(int x, int y);
~~~

Safe minimum value.

#### Function min()

~~~ .cpp
inline float min(float x, float y);
~~~

Safe minimum value.

#### Function min()

~~~ .cpp
inline int min(initializer_list<int> vs);
~~~

Safe minimum value.

#### Function min()

~~~ .cpp
inline float min(initializer_list<float> vs);
~~~

Safe minimum value.

#### Function max()

~~~ .cpp
inline int max(int x, int y);
~~~

Safe maximum value.

#### Function max()

~~~ .cpp
inline float max(float x, float y);
~~~

Safe maximum value.

#### Function max()

~~~ .cpp
inline int max(initializer_list<int> vs);
~~~

Safe maximum value.

#### Function max()

~~~ .cpp
inline float max(initializer_list<float> vs);
~~~

Safe maximum value.

#### Function clamp()

~~~ .cpp
inline int clamp(int x, int min_, int max_);
~~~

Clamp a value between a minimum and a maximum.

#### Function clamp()

~~~ .cpp
inline float clamp(float x, float min_, float max_);
~~~

Clamp a value between a minimum and a maximum.

#### Function lerp()

~~~ .cpp
inline float lerp(float a, float b, float t);
~~~

Linear interpolation.

#### Function bilerp()

~~~ .cpp
inline float bilerp(float aa, float ba, float ab, float bb, float s, float t);
~~~

bilinear interpolation

#### Function pow2()

~~~ .cpp
inline int pow2(int x);
~~~

Integer power of two

#### Function fastfloor()

~~~ .cpp
inline int fastfloor(float x);
~~~

Fast floor

#### Function float_to_byte()

~~~ .cpp
inline byte float_to_byte(float x);
~~~

Safe float to byte conversion

#### Function byte_to_float()

~~~ .cpp
inline float byte_to_float(byte x);
~~~

Safe byte to float conversion

#### Struct vec2f

~~~ .cpp
struct vec2f {
    vec2f(); 
    explicit vec2f(float vv); 
    vec2f(float x, float y); 
    float& operator[](int i); 
    const float& operator[](int i) const; 
    float* data(); 
    const float* data() const; 
    float x;
    float y;
}
~~~

Vector of 2 float elements.

- Members:
    - vec2f():      default constructor
    - vec2f():      element constructor
    - vec2f():      element constructor
    - operator[]():      element access
    - operator[]():      element access
    - data():      data access
    - data():      data access
    - x:      element data
    - y:      element data


#### Struct vec3f

~~~ .cpp
struct vec3f {
    vec3f(); 
    explicit vec3f(float vv); 
    vec3f(float x, float y, float z); 
    float& operator[](int i); 
    const float& operator[](int i) const; 
    float* data(); 
    const float* data() const; 
    float x;
    float y;
    float z;
}
~~~

Vector of 3 float elements.

- Members:
    - vec3f():      default constructor
    - vec3f():      element constructor
    - vec3f():      element constructor
    - operator[]():      element access
    - operator[]():      element access
    - data():      data access
    - data():      data access
    - x:      element data
    - y:      element data
    - z:      element data


#### Struct vec4f

~~~ .cpp
struct vec4f {
    vec4f(); 
    explicit vec4f(float vv); 
    vec4f(float x, float y, float z, float w); 
    vec4f(const vec3f& xyz, float w); 
    float& operator[](int i); 
    const float& operator[](int i) const; 
    float* data(); 
    const float* data() const; 
    vec3f& xyz(); 
    const vec3f& xyz() const; 
    float x;
    float y;
    float z;
    float w;
}
~~~

Vector of 4 float elements.

- Members:
    - vec4f():      default constructor
    - vec4f():      element constructor
    - vec4f():      element constructor
    - vec4f():      constructor from smaller vector
    - operator[]():      element access
    - operator[]():      element access
    - data():      data access
    - data():      data access
    - xyz():      access xyz components
    - xyz():      access xyz components
    - x:      element data
    - y:      element data
    - z:      element data
    - w:      element data


#### Struct vec2i

~~~ .cpp
struct vec2i {
    vec2i(); 
    explicit vec2i(int vv); 
    vec2i(int x, int y); 
    int& operator[](int i); 
    const int& operator[](int i) const; 
    int* data(); 
    const int* data() const; 
    int x;
    int y;
}
~~~

Vector of 2 int elements.

- Members:
    - vec2i():      default constructor
    - vec2i():      element constructor
    - vec2i():      element constructor
    - operator[]():      element access
    - operator[]():      element access
    - data():      data access
    - data():      data access
    - x:      element data
    - y:      element data


#### Struct vec3i

~~~ .cpp
struct vec3i {
    vec3i(); 
    explicit vec3i(int vv); 
    vec3i(int x, int y, int z); 
    int& operator[](int i); 
    const int& operator[](int i) const; 
    int* data(); 
    const int* data() const; 
    int x;
    int y;
    int z;
}
~~~

Vector of 3 int elements.

- Members:
    - vec3i():      default constructor
    - vec3i():      element constructor
    - vec3i():      element constructor
    - operator[]():      element access
    - operator[]():      element access
    - data():      data access
    - data():      data access
    - x:      element data
    - y:      element data
    - z:      element data


#### Struct vec4i

~~~ .cpp
struct vec4i {
    vec4i(); 
    explicit vec4i(int vv); 
    vec4i(int x, int y, int z, int w); 
    vec4i(const vec3i& xyz, int w); 
    int& operator[](int i); 
    const int& operator[](int i) const; 
    int* data(); 
    const int* data() const; 
    vec3i& xyz(); 
    const vec3i& xyz() const; 
    int x;
    int y;
    int z;
    int w;
}
~~~

Vector of 4 int elements.

- Members:
    - vec4i():      default constructor
    - vec4i():      element constructor
    - vec4i():      element constructor
    - vec4i():      constructor from smaller vector
    - operator[]():      element access
    - operator[]():      element access
    - data():      data access
    - data():      data access
    - xyz():      access xyz components
    - xyz():      access xyz components
    - x:      element data
    - y:      element data
    - z:      element data
    - w:      element data


#### Struct vec3b

~~~ .cpp
struct vec3b {
    vec3b(); 
    explicit vec3b(int vv); 
    vec3b(byte x, byte y, byte z); 
    byte& operator[](int i); 
    const byte& operator[](int i) const; 
    byte* data(); 
    const byte* data() const; 
    byte x;
    byte y;
    byte z;
}
~~~

Vector of 3 byte elements.

- Members:
    - vec3b():      default constructor
    - vec3b():      element constructor
    - vec3b():      element constructor
    - operator[]():      element access
    - operator[]():      element access
    - data():      data access
    - data():      data access
    - x:      element data
    - y:      element data
    - z:      element data


#### Struct vec4b

~~~ .cpp
struct vec4b {
    vec4b(); 
    explicit vec4b(byte vv); 
    vec4b(byte x, byte y, byte z, byte w); 
    vec4b(const vec3b& xyz, byte w); 
    byte& operator[](int i); 
    const byte& operator[](int i) const; 
    byte* data(); 
    const byte* data() const; 
    vec3b& xyz(); 
    const vec3b& xyz() const; 
    byte x;
    byte y;
    byte z;
    byte w;
}
~~~

Vector of 4 byte elements.

- Members:
    - vec4b():      default constructor
    - vec4b():      element constructor
    - vec4b():      element constructor
    - vec4b():      constructor from smaller vector
    - operator[]():      element access
    - operator[]():      element access
    - data():      data access
    - data():      data access
    - xyz():      access xyz components
    - xyz():      access xyz components
    - x:      element data
    - y:      element data
    - z:      element data
    - w:      element data


#### Constant zero2f

~~~ .cpp
const auto zero2f = vec2f();
~~~

2-dimensional float zero vector

#### Constant zero3f

~~~ .cpp
const auto zero3f = vec3f();
~~~

3-dimensional float zero vector

#### Constant zero4f

~~~ .cpp
const auto zero4f = vec4f();
~~~

4-dimensional float zero vector

#### Constant zero2i

~~~ .cpp
const auto zero2i = vec2i();
~~~

2-dimensional int zero vector

#### Constant zero3i

~~~ .cpp
const auto zero3i = vec3i();
~~~

3-dimensional int zero vector

#### Constant zero4i

~~~ .cpp
const auto zero4i = vec4i();
~~~

4-dimensional int zero vector

#### Constant zero4b

~~~ .cpp
const auto zero4b = vec4b();
~~~

4-dimensional byte zero vector

#### Function begin()

~~~ .cpp
inline int* begin(vec2i& a);
~~~

iteration support

#### Function begin()

~~~ .cpp
inline const int* begin(const vec2i& a);
~~~

iteration support

#### Function end()

~~~ .cpp
inline int* end(vec2i& a);
~~~

iteration support

#### Function end()

~~~ .cpp
inline const int* end(const vec2i& a);
~~~

iteration support

#### Function begin()

~~~ .cpp
inline int* begin(vec3i& a);
~~~

iteration support

#### Function begin()

~~~ .cpp
inline const int* begin(const vec3i& a);
~~~

iteration support

#### Function end()

~~~ .cpp
inline int* end(vec3i& a);
~~~

iteration support

#### Function end()

~~~ .cpp
inline const int* end(const vec3i& a);
~~~

iteration support

#### Function begin()

~~~ .cpp
inline int* begin(vec4i& a);
~~~

iteration support

#### Function begin()

~~~ .cpp
inline const int* begin(const vec4i& a);
~~~

iteration support

#### Function end()

~~~ .cpp
inline int* end(vec4i& a);
~~~

iteration support

#### Function end()

~~~ .cpp
inline const int* end(const vec4i& a);
~~~

iteration support

#### Function operator==()

~~~ .cpp
inline bool operator==(const vec2f& a, const vec2f& b);
~~~

vector operator ==

#### Function operator!=()

~~~ .cpp
inline bool operator!=(const vec2f& a, const vec2f& b);
~~~

vector operator !=

#### Function operator==()

~~~ .cpp
inline bool operator==(const vec3f& a, const vec3f& b);
~~~

vector operator ==

#### Function operator!=()

~~~ .cpp
inline bool operator!=(const vec3f& a, const vec3f& b);
~~~

vector operator !=

#### Function operator==()

~~~ .cpp
inline bool operator==(const vec4f& a, const vec4f& b);
~~~

vector operator ==

#### Function operator!=()

~~~ .cpp
inline bool operator!=(const vec4f& a, const vec4f& b);
~~~

vector operator !=

#### Function operator==()

~~~ .cpp
inline bool operator==(const vec2i& a, const vec2i& b);
~~~

vector operator ==

#### Function operator!=()

~~~ .cpp
inline bool operator!=(const vec2i& a, const vec2i& b);
~~~

vector operator !=

#### Function operator==()

~~~ .cpp
inline bool operator==(const vec3i& a, const vec3i& b);
~~~

vector operator ==

#### Function operator!=()

~~~ .cpp
inline bool operator!=(const vec3i& a, const vec3i& b);
~~~

vector operator !=

#### Function operator==()

~~~ .cpp
inline bool operator==(const vec4i& a, const vec4i& b);
~~~

vector operator ==

#### Function operator!=()

~~~ .cpp
inline bool operator!=(const vec4i& a, const vec4i& b);
~~~

vector operator !=

#### Function operator <()

~~~ .cpp
inline bool operator<(const vec2i& a, const vec2i& b);
~~~

vector operator < (lexicographic order - useful for map)

#### Function operator <()

~~~ .cpp
inline bool operator<(const vec3i& a, const vec3i& b);
~~~

vector operator < (lexicographic order - useful for map)

#### Function operator <()

~~~ .cpp
inline bool operator<(const vec4i& a, const vec4i& b);
~~~

vector operator < (lexicographic order - useful for map)

#### Function operator+()

~~~ .cpp
inline vec2f operator+(const vec2f& a);
~~~

vector operator +

#### Function operator-()

~~~ .cpp
inline vec2f operator-(const vec2f& a);
~~~

vector operator -

#### Function operator+()

~~~ .cpp
inline vec2f operator+(const vec2f& a, const vec2f& b);
~~~

vector operator +

#### Function operator-()

~~~ .cpp
inline vec2f operator-(const vec2f& a, const vec2f& b);
~~~

vector operator -

#### Function operator*()

~~~ .cpp
inline vec2f operator*(const vec2f& a, const vec2f& b);
~~~

vector operator *

#### Function operator*()

~~~ .cpp
inline vec2f operator*(const vec2f& a, float b);
~~~

vector operator *

#### Function operator*()

~~~ .cpp
inline vec2f operator*(float a, const vec2f& b);
~~~

vector operator *

#### Function operator/()

~~~ .cpp
inline vec2f operator/(const vec2f& a, const vec2f& b);
~~~

vector operator /

#### Function operator/()

~~~ .cpp
inline vec2f operator/(const vec2f& a, float b);
~~~

vector operator /

#### Function operator/()

~~~ .cpp
inline vec2f operator/(float a, const vec2f& b);
~~~

vector operator /

#### Function operator+()

~~~ .cpp
inline vec3f operator+(const vec3f& a);
~~~

vector operator +

#### Function operator-()

~~~ .cpp
inline vec3f operator-(const vec3f& a);
~~~

vector operator -

#### Function operator+()

~~~ .cpp
inline vec3f operator+(const vec3f& a, const vec3f& b);
~~~

vector operator +

#### Function operator-()

~~~ .cpp
inline vec3f operator-(const vec3f& a, const vec3f& b);
~~~

vector operator -

#### Function operator*()

~~~ .cpp
inline vec3f operator*(const vec3f& a, const vec3f& b);
~~~

vector operator *

#### Function operator*()

~~~ .cpp
inline vec3f operator*(const vec3f& a, float b);
~~~

vector operator *

#### Function operator*()

~~~ .cpp
inline vec3f operator*(float a, const vec3f& b);
~~~

vector operator *

#### Function operator/()

~~~ .cpp
inline vec3f operator/(const vec3f& a, const vec3f& b);
~~~

vector operator /

#### Function operator/()

~~~ .cpp
inline vec3f operator/(const vec3f& a, float b);
~~~

vector operator /

#### Function operator/()

~~~ .cpp
inline vec3f operator/(float a, const vec3f& b);
~~~

vector operator /

#### Function operator+()

~~~ .cpp
inline vec4f operator+(const vec4f& a);
~~~

vector operator +

#### Function operator-()

~~~ .cpp
inline vec4f operator-(const vec4f& a);
~~~

vector operator -

#### Function operator+()

~~~ .cpp
inline vec4f operator+(const vec4f& a, const vec4f& b);
~~~

vector operator +

#### Function operator-()

~~~ .cpp
inline vec4f operator-(const vec4f& a, const vec4f& b);
~~~

vector operator -

#### Function operator*()

~~~ .cpp
inline vec4f operator*(const vec4f& a, const vec4f& b);
~~~

vector operator *

#### Function operator*()

~~~ .cpp
inline vec4f operator*(const vec4f& a, float b);
~~~

vector operator *

#### Function operator*()

~~~ .cpp
inline vec4f operator*(float a, const vec4f& b);
~~~

vector operator *

#### Function operator/()

~~~ .cpp
inline vec4f operator/(const vec4f& a, const vec4f& b);
~~~

vector operator /

#### Function operator/()

~~~ .cpp
inline vec4f operator/(const vec4f& a, float b);
~~~

vector operator /

#### Function operator/()

~~~ .cpp
inline vec4f operator/(float a, const vec4f& b);
~~~

vector operator /

#### Function operator+=()

~~~ .cpp
inline vec2f& operator+=(vec2f& a, const vec2f& b);
~~~

vector operator +=

#### Function operator-=()

~~~ .cpp
inline vec2f& operator-=(vec2f& a, const vec2f& b);
~~~

vector operator -=

#### Function operator*=()

~~~ .cpp
inline vec2f& operator*=(vec2f& a, const vec2f& b);
~~~

vector operator *=

#### Function operator*=()

~~~ .cpp
inline vec2f& operator*=(vec2f& a, float b);
~~~

vector operator *=

#### Function operator/=()

~~~ .cpp
inline vec2f& operator/=(vec2f& a, const vec2f& b);
~~~

vector operator /=

#### Function operator/=()

~~~ .cpp
inline vec2f& operator/=(vec2f& a, float b);
~~~

vector operator /=

#### Function operator+=()

~~~ .cpp
inline vec3f& operator+=(vec3f& a, const vec3f& b);
~~~

vector operator +=

#### Function operator-=()

~~~ .cpp
inline vec3f& operator-=(vec3f& a, const vec3f& b);
~~~

vector operator -=

#### Function operator*=()

~~~ .cpp
inline vec3f& operator*=(vec3f& a, const vec3f& b);
~~~

vector operator *=

#### Function operator*=()

~~~ .cpp
inline vec3f& operator*=(vec3f& a, float b);
~~~

vector operator *=

#### Function operator/=()

~~~ .cpp
inline vec3f& operator/=(vec3f& a, const vec3f& b);
~~~

vector operator /=

#### Function operator/=()

~~~ .cpp
inline vec3f& operator/=(vec3f& a, float b);
~~~

vector operator /=

#### Function operator+=()

~~~ .cpp
inline vec4f& operator+=(vec4f& a, const vec4f& b);
~~~

vector operator +=

#### Function operator-=()

~~~ .cpp
inline vec4f& operator-=(vec4f& a, const vec4f& b);
~~~

vector operator -=

#### Function operator*=()

~~~ .cpp
inline vec4f& operator*=(vec4f& a, const vec4f& b);
~~~

vector operator *=

#### Function operator*=()

~~~ .cpp
inline vec4f& operator*=(vec4f& a, float b);
~~~

vector operator *=

#### Function operator/=()

~~~ .cpp
inline vec4f& operator/=(vec4f& a, const vec4f& b);
~~~

vector operator /=

#### Function operator/=()

~~~ .cpp
inline vec4f& operator/=(vec4f& a, float b);
~~~

vector operator /=

#### Function operator+()

~~~ .cpp
inline vec2i operator+(const vec2i& a);
~~~

vector operator +

#### Function operator-()

~~~ .cpp
inline vec2i operator-(const vec2i& a);
~~~

vector operator -

#### Function operator+()

~~~ .cpp
inline vec2i operator+(const vec2i& a, const vec2i& b);
~~~

vector operator +

#### Function operator-()

~~~ .cpp
inline vec2i operator-(const vec2i& a, const vec2i& b);
~~~

vector operator -

#### Function operator+()

~~~ .cpp
inline vec3i operator+(const vec3i& a);
~~~

vector operator +

#### Function operator-()

~~~ .cpp
inline vec3i operator-(const vec3i& a);
~~~

vector operator -

#### Function operator+()

~~~ .cpp
inline vec3i operator+(const vec3i& a, const vec3i& b);
~~~

vector operator +

#### Function operator-()

~~~ .cpp
inline vec3i operator-(const vec3i& a, const vec3i& b);
~~~

vector operator -

#### Function operator+()

~~~ .cpp
inline vec4i operator+(const vec4i& a);
~~~

vector operator +

#### Function operator-()

~~~ .cpp
inline vec4i operator-(const vec4i& a);
~~~

vector operator -

#### Function operator+()

~~~ .cpp
inline vec4i operator+(const vec4i& a, const vec4i& b);
~~~

vector operator +

#### Function operator-()

~~~ .cpp
inline vec4i operator-(const vec4i& a, const vec4i& b);
~~~

vector operator -

#### Function operator+=()

~~~ .cpp
inline vec2i& operator+=(vec2i& a, const vec2i& b);
~~~

vector operator +=

#### Function operator-=()

~~~ .cpp
inline vec2i& operator-=(vec2i& a, const vec2i& b);
~~~

vector operator -=

#### Function operator+=()

~~~ .cpp
inline vec3i& operator+=(vec3i& a, const vec3i& b);
~~~

vector operator +=

#### Function operator-=()

~~~ .cpp
inline vec3i& operator-=(vec3i& a, const vec3i& b);
~~~

vector operator -=

#### Function operator+=()

~~~ .cpp
inline vec4i& operator+=(vec4i& a, const vec4i& b);
~~~

vector operator +=

#### Function operator-=()

~~~ .cpp
inline vec4i& operator-=(vec4i& a, const vec4i& b);
~~~

vector operator -=

#### Function dot()

~~~ .cpp
inline float dot(const vec2f& a, const vec2f& b);
~~~

vector dot product

#### Function dot()

~~~ .cpp
inline float dot(const vec3f& a, const vec3f& b);
~~~

vector dot product

#### Function dot()

~~~ .cpp
inline float dot(const vec4f& a, const vec4f& b);
~~~

vector dot product

#### Function cross()

~~~ .cpp
inline float cross(const vec2f& a, const vec2f& b);
~~~

vector cross product

#### Function cross()

~~~ .cpp
inline vec3f cross(const vec3f& a, const vec3f& b);
~~~

vector cross product

#### Function length()

~~~ .cpp
inline float length(const vec2f& a);
~~~

vector length

#### Function length()

~~~ .cpp
inline float length(const vec3f& a);
~~~

vector length

#### Function length()

~~~ .cpp
inline float length(const vec4f& a);
~~~

vector length

#### Function normalize()

~~~ .cpp
inline vec2f normalize(const vec2f& a);
~~~

vector normalization

#### Function normalize()

~~~ .cpp
inline vec3f normalize(const vec3f& a);
~~~

vector normalization

#### Function normalize()

~~~ .cpp
inline vec4f normalize(const vec4f& a);
~~~

vector normalization

#### Function uangle()

~~~ .cpp
inline float uangle(const vec3f& a, const vec3f& b);
~~~

angle between normalized vectors

#### Function uangle()

~~~ .cpp
inline float uangle(const vec4f& a, const vec4f& b);
~~~

angle between normalized vectors

#### Function angle()

~~~ .cpp
inline float angle(const vec3f& a, const vec3f& b);
~~~

angle between vectors

#### Function lerp()

~~~ .cpp
inline vec2f lerp(const vec2f& a, const vec2f& b, float t);
~~~

vector linear interpolation

#### Function lerp()

~~~ .cpp
inline vec3f lerp(const vec3f& a, const vec3f& b, float t);
~~~

vector linear interpolation

#### Function lerp()

~~~ .cpp
inline vec4f lerp(const vec4f& a, const vec4f& b, float t);
~~~

vector linear interpolation

#### Function bilerp()

~~~ .cpp
inline vec3f bilerp(const vec3f& aa, const vec3f& ba, const vec3f& ab,
    const vec3f& bb, float s, float t);
~~~

vector bilinear interpolation

#### Function nlerp()

~~~ .cpp
inline vec3f nlerp(const vec3f& a, const vec3f& b, float t);
~~~

vector normalized linear interpolation

#### Function slerp()

~~~ .cpp
inline vec3f slerp(const vec3f& a, const vec3f& b, float t);
~~~

vector spherical linear interpolation (vectors have to be normalized)

#### Function nlerp()

~~~ .cpp
inline vec4f nlerp(const vec4f& a, const vec4f& b, float t);
~~~

vector normalized linear interpolation

#### Function slerp()

~~~ .cpp
inline vec4f slerp(const vec4f& a, const vec4f& b, float t);
~~~

vector spherical linear interpolation (vectors have to be normalized)

#### Function orthogonal()

~~~ .cpp
inline vec3f orthogonal(const vec3f& v);
~~~

orthogonal vector

#### Function orthonormalize()

~~~ .cpp
inline vec3f orthonormalize(const vec3f& a, const vec3f& b);
~~~

orthonormalize two vectors

#### Function clamp()

~~~ .cpp
inline vec2f clamp(const vec2f& x, float min, float max);
~~~

vector component-wise clamp

#### Function clamp()

~~~ .cpp
inline vec3f clamp(const vec3f& x, float min, float max);
~~~

vector component-wise clamp

#### Function clamp()

~~~ .cpp
inline vec4f clamp(const vec4f& x, float min, float max);
~~~

vector component-wise clamp

#### Function clamplen()

~~~ .cpp
inline vec2f clamplen(const vec2f& x, float max);
~~~

clamp the length of a vector

#### Function clamplen()

~~~ .cpp
inline vec3f clamplen(const vec3f& x, float max);
~~~

clamp the length of a vector

#### Function clamplen()

~~~ .cpp
inline vec4f clamplen(const vec4f& x, float max);
~~~

clamp the length of a vector

#### Function min_element()

~~~ .cpp
inline pair<int, float> min_element(const vec2f& a);
~~~

min vector element

#### Function min_element()

~~~ .cpp
inline pair<int, float> min_element(const vec3f& a);
~~~

min vector element

#### Function min_element()

~~~ .cpp
inline pair<int, float> min_element(const vec4f& a);
~~~

min vector element

#### Function _max_element()

~~~ .cpp
inline pair<int, float> _max_element(int N, const float* a);
~~~

index of the max vector element

#### Function max_element()

~~~ .cpp
inline pair<int, float> max_element(const vec2f& a);
~~~

index of the min vector element

#### Function max_element()

~~~ .cpp
inline pair<int, float> max_element(const vec3f& a);
~~~

index of the min vector element

#### Function max_element()

~~~ .cpp
inline pair<int, float> max_element(const vec4f& a);
~~~

index of the min vector element

#### Function float_to_byte()

~~~ .cpp
inline vec3b float_to_byte(const vec3f& a);
~~~

Element-wise conversion

#### Function byte_to_float()

~~~ .cpp
inline vec3f byte_to_float(const vec3b& a);
~~~

Element-wise conversion

#### Function float_to_byte()

~~~ .cpp
inline vec4b float_to_byte(const vec4f& a);
~~~

Element-wise conversion

#### Function byte_to_float()

~~~ .cpp
inline vec4f byte_to_float(const vec4b& a);
~~~

Element-wise conversion

#### Function operator < <()

~~~ .cpp
inline ostream& operator<<(ostream& os, const vec2f& a);
~~~

stream write

#### Function operator < <()

~~~ .cpp
inline ostream& operator<<(ostream& os, const vec3f& a);
~~~

stream write

#### Function operator < <()

~~~ .cpp
inline ostream& operator<<(ostream& os, const vec4f& a);
~~~

stream write

#### Function operator < <()

~~~ .cpp
inline ostream& operator<<(ostream& os, const vec2i& a);
~~~

stream write

#### Function operator < <()

~~~ .cpp
inline ostream& operator<<(ostream& os, const vec3i& a);
~~~

stream write

#### Function operator < <()

~~~ .cpp
inline ostream& operator<<(ostream& os, const vec4i& a);
~~~

stream write

#### Function operator \> \>()

~~~ .cpp
inline istream& operator>>(istream& is, vec2f& a);
~~~

stream read

#### Function operator \> \>()

~~~ .cpp
inline istream& operator>>(istream& is, vec3f& a);
~~~

stream read

#### Function operator \> \>()

~~~ .cpp
inline istream& operator>>(istream& is, vec4f& a);
~~~

stream read

#### Function operator \> \>()

~~~ .cpp
inline istream& operator>>(istream& is, vec2i& a);
~~~

stream read

#### Function operator \> \>()

~~~ .cpp
inline istream& operator>>(istream& is, vec3i& a);
~~~

stream read

#### Function operator \> \>()

~~~ .cpp
inline istream& operator>>(istream& is, vec4i& a);
~~~

stream read

#### Struct hash <ygl::vec2i \>

~~~ .cpp
template <>
struct hash<ygl::vec2i> {
~~~

Hash functor for vector for use with unordered_map

#### Struct hash <ygl::vec3i \>

~~~ .cpp
template <>
struct hash<ygl::vec3i> {
~~~

Hash functor for vector for use with unordered_map

#### Struct hash <ygl::vec4i \>

~~~ .cpp
template <>
struct hash<ygl::vec4i> {
~~~

Hash functor for vector for use with unordered_map

#### Struct mat2f

~~~ .cpp
struct mat2f {
    mat2f(); 
    explicit mat2f(float vv); 
    mat2f(const vec2f& x, const vec2f& y); 
    vec2f& operator[](int i); 
    const vec2f& operator[](int i) const; 
    vec2f* data(); 
    const vec2f* data() const; 
    vec2f x;
    vec2f y;
}
~~~

Matrix of 2x2 elements stored in column major format.
Colums access via operator[].

- Members:
    - mat2f():      default constructor
    - mat2f():      diagonal constructor
    - mat2f():      list constructor
    - operator[]():      element access
    - operator[]():      element access
    - data():      data access
    - data():      data access
    - x:      element data
    - y:      element data


#### Struct mat3f

~~~ .cpp
struct mat3f {
    mat3f(); 
    explicit mat3f(float vv); 
    mat3f(const vec3f& x, const vec3f& y, const vec3f& z); 
    vec3f& operator[](int i); 
    const vec3f& operator[](int i) const; 
    vec3f* data(); 
    const vec3f* data() const; 
    vec3f x;
    vec3f y;
    vec3f z;
}
~~~

Matrix of 3x3 elements stored in column major format.
Colums access via operator[].

- Members:
    - mat3f():      default constructor
    - mat3f():      diagonal constructor
    - mat3f():      list constructor
    - operator[]():      element access
    - operator[]():      element access
    - data():      data access
    - data():      data access
    - x:      element data
    - y:      element data
    - z:      element data


#### Struct mat4f

~~~ .cpp
struct mat4f {
    mat4f(); 
    explicit mat4f(float vv); 
    mat4f(const vec4f& x, const vec4f& y, const vec4f& z, const vec4f& w); 
    vec4f& operator[](int i); 
    const vec4f& operator[](int i) const; 
    vec4f* data(); 
    const vec4f* data() const; 
    vec4f x;
    vec4f y;
    vec4f z;
    vec4f w;
}
~~~

Matrix of 4x4 elements stored in column major format.
Colums access via operator[].

- Members:
    - mat4f():      default constructor
    - mat4f():      diagonal constructor
    - mat4f():      list constructor
    - operator[]():      element access
    - operator[]():      element access
    - data():      data access
    - data():      data access
    - x:      element data
    - y:      element data
    - z:      element data
    - w:      element data


#### Constant identity_mat2f

~~~ .cpp
const auto identity_mat2f = mat2f();
~~~

2-dimensional float identity matrix

#### Constant identity_mat3f

~~~ .cpp
const auto identity_mat3f = mat3f();
~~~

3-dimensional float identity matrix

#### Constant identity_mat4f

~~~ .cpp
const auto identity_mat4f = mat4f();
~~~

4-dimensional float identity matrix

#### Function operator==()

~~~ .cpp
inline bool operator==(const mat2f& a, const mat2f& b);
~~~

matrix operator ==

#### Function operator!=()

~~~ .cpp
inline bool operator!=(const mat2f& a, const mat2f& b);
~~~

matrix operator !=

#### Function operator==()

~~~ .cpp
inline bool operator==(const mat3f& a, const mat3f& b);
~~~

matrix operator ==

#### Function operator!=()

~~~ .cpp
inline bool operator!=(const mat3f& a, const mat3f& b);
~~~

matrix operator !=

#### Function operator==()

~~~ .cpp
inline bool operator==(const mat4f& a, const mat4f& b);
~~~

matrix operator ==

#### Function operator!=()

~~~ .cpp
inline bool operator!=(const mat4f& a, const mat4f& b);
~~~

matrix operator !=

#### Function operator+()

~~~ .cpp
inline mat2f operator+(const mat2f& a, const mat2f& b);
~~~

matrix operator +

#### Function operator+()

~~~ .cpp
inline mat3f operator+(const mat3f& a, const mat3f& b);
~~~

matrix operator +

#### Function operator+()

~~~ .cpp
inline mat4f operator+(const mat4f& a, const mat4f& b);
~~~

matrix operator +

#### Function operator*()

~~~ .cpp
inline mat2f operator*(const mat2f& a, float b);
~~~

matrix scalar multiply

#### Function operator/()

~~~ .cpp
inline mat2f operator/(const mat2f& a, float b);
~~~

matrix scalar division

#### Function operator*()

~~~ .cpp
inline vec2f operator*(const mat2f& a, const vec2f& b);
~~~

matrix-vector right multiply

#### Function operator*()

~~~ .cpp
inline vec2f operator*(const vec2f& a, const mat2f& b);
~~~

matrix-vector left multiply

#### Function operator*()

~~~ .cpp
inline mat2f operator*(const mat2f& a, const mat2f& b);
~~~

matrix-matrix multiply

#### Function operator*()

~~~ .cpp
inline mat3f operator*(const mat3f& a, float b);
~~~

matrix scalar multiply

#### Function operator/()

~~~ .cpp
inline mat3f operator/(const mat3f& a, float b);
~~~

matrix scalar division

#### Function operator*()

~~~ .cpp
inline mat4f operator*(const mat4f& a, float b);
~~~

matrix scalar multiply

#### Function operator/()

~~~ .cpp
inline mat4f operator/(const mat4f& a, float b);
~~~

matrix scalar division

#### Function operator*()

~~~ .cpp
inline vec3f operator*(const mat3f& a, const vec3f& b);
~~~

matrix-vector right multiply

#### Function operator*()

~~~ .cpp
inline vec3f operator*(const vec3f& a, const mat3f& b);
~~~

matrix-vector left multiply

#### Function operator*()

~~~ .cpp
inline mat3f operator*(const mat3f& a, const mat3f& b);
~~~

matrix-matrix multiply

#### Function operator*()

~~~ .cpp
inline vec4f operator*(const mat4f& a, const vec4f& b);
~~~

matrix-vector right multiply

#### Function operator*()

~~~ .cpp
inline vec4f operator*(const vec4f& a, const mat4f& b);
~~~

matrix-vector left multiply

#### Function operator*()

~~~ .cpp
inline mat4f operator*(const mat4f& a, const mat4f& b);
~~~

matrix-matrix multiply

#### Function operator+=()

~~~ .cpp
inline mat2f& operator+=(mat2f& a, const mat2f& b);
~~~

matrix sum assignment

#### Function operator*=()

~~~ .cpp
inline mat2f& operator*=(mat2f& a, const mat2f& b);
~~~

matrix-matrix multiply assignment

#### Function operator*=()

~~~ .cpp
inline mat2f& operator*=(mat2f& a, float b);
~~~

matrix scaling assignment

#### Function operator/=()

~~~ .cpp
inline mat2f& operator/=(mat2f& a, float b);
~~~

matrix scaling assignment

#### Function operator+=()

~~~ .cpp
inline mat3f& operator+=(mat3f& a, const mat3f& b);
~~~

matrix sum assignment

#### Function operator*=()

~~~ .cpp
inline mat3f& operator*=(mat3f& a, const mat3f& b);
~~~

matrix-matrix multiply assignment

#### Function operator*=()

~~~ .cpp
inline mat3f& operator*=(mat3f& a, float b);
~~~

matrix scaling assignment

#### Function operator/=()

~~~ .cpp
inline mat3f& operator/=(mat3f& a, float b);
~~~

matrix scaling assignment

#### Function operator+=()

~~~ .cpp
inline mat4f& operator+=(mat4f& a, const mat4f& b);
~~~

matrix sum assignment

#### Function operator*=()

~~~ .cpp
inline mat4f& operator*=(mat4f& a, const mat4f& b);
~~~

matrix-matrix multiply assignment

#### Function operator*=()

~~~ .cpp
inline mat4f& operator*=(mat4f& a, float b);
~~~

matrix scaling assignment

#### Function operator/=()

~~~ .cpp
inline mat4f& operator/=(mat4f& a, float b);
~~~

matrix scaling assignment

#### Function mat_diagonal()

~~~ .cpp
inline vec2f mat_diagonal(const mat2f& a);
~~~

matrix diagonal

#### Function mat_diagonal()

~~~ .cpp
inline vec3f mat_diagonal(const mat3f& a);
~~~

matrix diagonal

#### Function mat_diagonal()

~~~ .cpp
inline vec4f mat_diagonal(const mat4f& a);
~~~

matrix diagonal

#### Function transpose()

~~~ .cpp
inline mat2f transpose(const mat2f& a);
~~~

matrix transpose

#### Function transpose()

~~~ .cpp
inline mat3f transpose(const mat3f& a);
~~~

matrix transpose

#### Function transpose()

~~~ .cpp
inline mat4f transpose(const mat4f& a);
~~~

matrix transpose

#### Function adjugate()

~~~ .cpp
inline mat2f adjugate(const mat2f& a);
~~~

matrix adjugate

#### Function adjugate()

~~~ .cpp
inline mat3f adjugate(const mat3f& a);
~~~

matrix adjugate

#### Function adjugate()

~~~ .cpp
inline mat4f adjugate(const mat4f& a);
~~~

matrix adjugate

#### Function determinant()

~~~ .cpp
inline float determinant(const mat2f& a);
~~~

matrix determinant

#### Function determinant()

~~~ .cpp
inline float determinant(const mat3f& a);
~~~

matrix determinant

#### Function determinant()

~~~ .cpp
inline float determinant(const mat4f& a);
~~~

matrix determinant

#### Function inverse()

~~~ .cpp
inline mat2f inverse(const mat2f& a);
~~~

matrix inverse (uses adjugate and determinant)

#### Function inverse()

~~~ .cpp
inline mat3f inverse(const mat3f& a);
~~~

matrix inverse (uses adjugate and determinant)

#### Function inverse()

~~~ .cpp
inline mat4f inverse(const mat4f& a);
~~~

matrix inverse (uses adjugate and determinant)

#### Function operator < <()

~~~ .cpp
inline ostream& operator<<(ostream& os, const mat2f& a);
~~~

stream write

#### Function operator \> \>()

~~~ .cpp
inline istream& operator>>(istream& is, mat2f& a);
~~~

stream read

#### Function operator < <()

~~~ .cpp
inline ostream& operator<<(ostream& os, const mat3f& a);
~~~

stream write

#### Function operator \> \>()

~~~ .cpp
inline istream& operator>>(istream& is, mat3f& a);
~~~

stream read

#### Function operator < <()

~~~ .cpp
inline ostream& operator<<(ostream& os, const mat4f& a);
~~~

stream write

#### Function operator \> \>()

~~~ .cpp
inline istream& operator>>(istream& is, mat4f& a);
~~~

stream read

#### Struct frame3f

~~~ .cpp
struct frame3f {
    static const int N = 3;
    frame3f(); 
    frame3f(const vec3f& x, const vec3f& y, const vec3f& z, const vec3f& o); 
    frame3f(const mat3f& m, const vec3f& t); 
    frame3f(const mat4f& m); 
    vec3f& operator[](int i); 
    const vec3f& operator[](int i) const; 
    vec3f* data(); 
    const vec3f* data() const; 
    vec3f& pos(); 
    const vec3f& pos() const; 
    mat3f& rot(); 
    const mat3f& rot() const; 
    vec3f x;
    vec3f y;
    vec3f z;
    vec3f o;
}
~~~

Rigid transforms stored as a column-major affine matrix.
In memory, this representation is equivalent to storing an 3x3 rotation
followed by a 3x1 translation. Viewed this way, the representation allows
also to retrive the axis of the coordinate frame as the first 3 columns and
the translation as the 4th column. Colums access via operator[].
Access rotation and position with pos() and rot().

- Members:
    - N:      size
    - frame3f():      default constructor
    - frame3f():      element constructor
    - frame3f():      element constructor
    - frame3f():      conversion from matrix (assumes the matrix is a frame, so dangerous!)
    - operator[]():      element access
    - operator[]():      element access
    - data():      data access
    - data():      data access
    - pos():      access position
    - pos():      access position
    - rot():      access rotation
    - rot():      access rotation
    - x:      element data
    - y:      element data
    - z:      element data
    - o:      element data


#### Constant identity_frame3f

~~~ .cpp
const auto identity_frame3f =
    frame3f{ {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0} };
~~~

indentity frame

#### Function begin()

~~~ .cpp
inline vec3f* begin(frame3f& a);
~~~

iteration support

#### Function begin()

~~~ .cpp
inline const vec3f* begin(const frame3f& a);
~~~

iteration support

#### Function end()

~~~ .cpp
inline vec3f* end(frame3f& a);
~~~

iteration support

#### Function end()

~~~ .cpp
inline const vec3f* end(const frame3f& a);
~~~

iteration support

#### Function to_mat4f()

~~~ .cpp
inline mat4f to_mat4f(const frame3f& a);
~~~

frame to matrix conversion

#### Function to_frame3f()

~~~ .cpp
inline frame3f to_frame3f(const mat4f& a);
~~~

matrix to frame conversion

#### Function operator==()

~~~ .cpp
inline bool operator==(const frame3f& a, const frame3f& b);
~~~

vector operator ==

#### Function operator!=()

~~~ .cpp
inline bool operator!=(const frame3f& a, const frame3f& b);
~~~

vector operator !=

#### Function operator*()

~~~ .cpp
inline frame3f operator*(const frame3f& a, const frame3f& b);
~~~

frame composition (equivalent to affine matrix multiply)

#### Function inverse()

~~~ .cpp
inline frame3f inverse(const frame3f& a);
~~~

frame inverse (equivalent to rigid affine inverse)

#### Function operator < <()

~~~ .cpp
inline ostream& operator<<(ostream& os, const frame3f& a);
~~~

stream write

#### Function operator \> \>()

~~~ .cpp
inline istream& operator>>(istream& is, frame3f& a);
~~~

stream read

#### Struct quat4f

~~~ .cpp
struct quat4f {
    quat4f(); 
    explicit quat4f(const vec4f& vv); 
    explicit operator vec4f() const; 
    float& operator[](int i); 
    const float& operator[](int i) const; 
    float* data(); 
    const float* data() const; 
    float x;
    float y;
    float z;
    float w;
}
~~~

Quaternions implemented as a vec<T,4>. Data access via operator[].
Quaterions are xi + yj + zk + w.

- Members:
    - quat4f():      default constructor
    - quat4f():      conversion from vec
    - operator vec4f():      conversion to vec
    - operator[]():      element access
    - operator[]():      element access
    - data():      data access
    - data():      data access
    - x:      data
    - y:      data
    - z:      data
    - w:      data


#### Constant identity_quat4f

~~~ .cpp
const auto identity_quat4f = quat4f{0, 0, 0, 1};
~~~

float identity quaterion

#### Function operator==()

~~~ .cpp
inline bool operator==(const quat4f& a, const quat4f& b);
~~~

vector operator ==

#### Function operator!=()

~~~ .cpp
inline bool operator!=(const quat4f& a, const quat4f& b);
~~~

vector operator !=

#### Function operator*()

~~~ .cpp
inline quat4f operator*(const quat4f& a, const quat4f& b);
~~~

quaterion multiply

#### Function operator*()

~~~ .cpp
inline quat4f operator*(const quat4f& a, float b);
~~~

quaterion multiply

#### Function operator/()

~~~ .cpp
inline quat4f operator/(const quat4f& a, float b);
~~~

quaterion division

#### Function conjugate()

~~~ .cpp
inline quat4f conjugate(const quat4f& v);
~~~

quaterion conjugate

#### Function inverse()

~~~ .cpp
inline quat4f inverse(const quat4f& v);
~~~

quaterion inverse

#### Function normalize()

~~~ .cpp
inline quat4f normalize(const quat4f& v);
~~~

quaterion inverse

#### Function nlerp()

~~~ .cpp
inline quat4f nlerp(const quat4f& a, const quat4f& b, float t);
~~~

quaterion normalized linear interpolation

#### Function slerp()

~~~ .cpp
inline quat4f slerp(const quat4f& a, const quat4f& b, float t);
~~~

quaterion spherical linear interpolation

#### Function operator < <()

~~~ .cpp
inline ostream& operator<<(ostream& os, const quat4f& a);
~~~

stream write

#### Function operator \> \>()

~~~ .cpp
inline istream& operator>>(istream& is, quat4f& a);
~~~

stream read

#### Struct bbox1f

~~~ .cpp
struct bbox1f {
    bbox1f(); 
    bbox1f(float m, float M); 
    float& operator[](int i); 
    const float& operator[](int i) const; 
    float min;
    float max;
}
~~~

Axis aligned bounding box represented as a min/max vector pair.

- Members:
    - bbox1f():      initializes an invalid bbox
    - bbox1f():      list constructor
    - operator[]():      element access
    - operator[]():      element access
    - min:      element data
    - max:      element data


#### Struct bbox2f

~~~ .cpp
struct bbox2f {
    bbox2f(); 
    bbox2f(const vec2f& m, const vec2f& M); 
    vec2f& operator[](int i); 
    const vec2f& operator[](int i) const; 
    vec2f min;
    vec2f max;
}
~~~

Axis aligned bounding box represented as a min/max vector pair.

- Members:
    - bbox2f():      initializes an invalid bbox
    - bbox2f():      list constructor
    - operator[]():      element access
    - operator[]():      element access
    - min:      element data
    - max:      element data


#### Struct bbox3f

~~~ .cpp
struct bbox3f {
    bbox3f(); 
    bbox3f(const vec3f& m, const vec3f& M); 
    vec3f& operator[](int i); 
    const vec3f& operator[](int i) const; 
    vec3f min;
    vec3f max;
}
~~~

Axis aligned bounding box represented as a min/max vector pair.

- Members:
    - bbox3f():      initializes an invalid bbox
    - bbox3f():      list constructor
    - operator[]():      element access
    - operator[]():      element access
    - min:      element data
    - max:      element data


#### Struct bbox4f

~~~ .cpp
struct bbox4f {
    bbox4f(); 
    bbox4f(const vec4f& m, const vec4f& M); 
    vec4f& operator[](int i); 
    const vec4f& operator[](int i) const; 
    vec4f min;
    vec4f max;
}
~~~

Axis aligned bounding box represented as a min/max vector pair.

- Members:
    - bbox4f():      initializes an invalid bbox
    - bbox4f():      list constructor
    - operator[]():      element access
    - operator[]():      element access
    - min:      element data
    - max:      element data


#### Constant invalid_bbox1f

~~~ .cpp
const auto invalid_bbox1f = bbox1f();
~~~

1-dimensional float empty bbox

#### Constant invalid_bbox2f

~~~ .cpp
const auto invalid_bbox2f = bbox2f();
~~~

2-dimensional float empty bbox

#### Constant invalid_bbox3f

~~~ .cpp
const auto invalid_bbox3f = bbox3f();
~~~

3-dimensional float empty bbox

#### Constant invalid_bbox4f

~~~ .cpp
const auto invalid_bbox4f = bbox4f();
~~~

4-dimensional float empty bbox

#### Function operator==()

~~~ .cpp
inline bool operator==(const bbox1f& a, const bbox1f& b);
~~~

bbox operator ==

#### Function operator!=()

~~~ .cpp
inline bool operator!=(const bbox1f& a, const bbox1f& b);
~~~

bbox operator !=

#### Function operator==()

~~~ .cpp
inline bool operator==(const bbox2f& a, const bbox2f& b);
~~~

bbox operator ==

#### Function operator!=()

~~~ .cpp
inline bool operator!=(const bbox2f& a, const bbox2f& b);
~~~

bbox operator !=

#### Function operator==()

~~~ .cpp
inline bool operator==(const bbox3f& a, const bbox3f& b);
~~~

bbox operator ==

#### Function operator!=()

~~~ .cpp
inline bool operator!=(const bbox3f& a, const bbox3f& b);
~~~

bbox operator !=

#### Function operator==()

~~~ .cpp
inline bool operator==(const bbox4f& a, const bbox4f& b);
~~~

bbox operator ==

#### Function operator!=()

~~~ .cpp
inline bool operator!=(const bbox4f& a, const bbox4f& b);
~~~

bbox operator !=

#### Function bbox_center()

~~~ .cpp
inline float bbox_center(const bbox1f& a);
~~~

computes the center of a bbox

#### Function bbox_diagonal()

~~~ .cpp
inline float bbox_diagonal(const bbox1f& a);
~~~

computes the diagonal of a bbox

#### Function bbox_center()

~~~ .cpp
inline vec2f bbox_center(const bbox2f& a);
~~~

computes the center of a bbox

#### Function bbox_diagonal()

~~~ .cpp
inline vec2f bbox_diagonal(const bbox2f& a);
~~~

computes the diagonal of a bbox

#### Function bbox_center()

~~~ .cpp
inline vec3f bbox_center(const bbox3f& a);
~~~

computes the center of a bbox

#### Function bbox_diagonal()

~~~ .cpp
inline vec3f bbox_diagonal(const bbox3f& a);
~~~

computes the diagonal of a bbox

#### Function bbox_center()

~~~ .cpp
inline vec4f bbox_center(const bbox4f& a);
~~~

computes the center of a bbox

#### Function bbox_diagonal()

~~~ .cpp
inline vec4f bbox_diagonal(const bbox4f& a);
~~~

computes the diagonal of a bbox

#### Function expand()

~~~ .cpp
inline bbox1f expand(const bbox1f& a, float b);
~~~

expands a bounding box with a point

#### Function expand()

~~~ .cpp
inline bbox2f expand(const bbox2f& a, const vec2f& b);
~~~

expands a bounding box with a point

#### Function expand()

~~~ .cpp
inline bbox3f expand(const bbox3f& a, const vec3f& b);
~~~

expands a bounding box with a point

#### Function expand()

~~~ .cpp
inline bbox4f expand(const bbox4f& a, const vec4f& b);
~~~

expands a bounding box with a point

#### Function expand()

~~~ .cpp
inline bbox1f expand(const bbox1f& a, const bbox1f& b);
~~~

expands a bounding box with a bounding box

#### Function expand()

~~~ .cpp
inline bbox2f expand(const bbox2f& a, const bbox2f& b);
~~~

expands a bounding box with a bounding box

#### Function expand()

~~~ .cpp
inline bbox3f expand(const bbox3f& a, const bbox3f& b);
~~~

expands a bounding box with a bounding box

#### Function expand()

~~~ .cpp
inline bbox4f expand(const bbox4f& a, const bbox4f& b);
~~~

expands a bounding box with a bounding box

#### Function contains()

~~~ .cpp
inline bool contains(const bbox3f& a, const vec3f& b);
~~~

check if a bounding box contains a point

#### Function contains()

~~~ .cpp
inline bool contains(const bbox3f& a, const bbox3f& b);
~~~

check if a bounding box contains a bounding box

#### Function operator+=()

~~~ .cpp
inline bbox1f& operator+=(bbox1f& a, float b);
~~~

assign to expand()

#### Function operator+=()

~~~ .cpp
inline bbox1f& operator+=(bbox1f& a, const bbox1f& b);
~~~

assign to expand()

#### Function operator+=()

~~~ .cpp
inline bbox2f& operator+=(bbox2f& a, const vec2f& b);
~~~

assign to expand()

#### Function operator+=()

~~~ .cpp
inline bbox2f& operator+=(bbox2f& a, const bbox2f& b);
~~~

assign to expand()

#### Function operator+=()

~~~ .cpp
inline bbox3f& operator+=(bbox3f& a, const vec3f& b);
~~~

assign to expand()

#### Function operator+=()

~~~ .cpp
inline bbox3f& operator+=(bbox3f& a, const bbox3f& b);
~~~

assign to expand()

#### Function operator+=()

~~~ .cpp
inline bbox4f& operator+=(bbox4f& a, const vec4f& b);
~~~

assign to expand()

#### Function operator+=()

~~~ .cpp
inline bbox4f& operator+=(bbox4f& a, const bbox4f& b);
~~~

assign to expand()

#### Function make_bbox()

~~~ .cpp
inline bbox1f make_bbox(int count, const float* v);
~~~

initialize a bonding box from a list of points

#### Function make_bbox()

~~~ .cpp
inline bbox2f make_bbox(int count, const vec2f* v);
~~~

initialize a bonding box from a list of points

#### Function make_bbox()

~~~ .cpp
inline bbox3f make_bbox(int count, const vec3f* v);
~~~

initialize a bonding box from a list of points

#### Function make_bbox()

~~~ .cpp
inline bbox4f make_bbox(int count, const vec4f* v);
~~~

initialize a bonding box from a list of points

#### Function make_bbox()

~~~ .cpp
inline bbox3f make_bbox(const initializer_list<vec3f>& v);
~~~

initialize a bonding box from a list of points

#### Function operator < <()

~~~ .cpp
inline ostream& operator<<(ostream& os, const bbox1f& a);
~~~

stream write

#### Function operator \> \>()

~~~ .cpp
inline istream& operator>>(istream& is, bbox1f& a);
~~~

stream read

#### Function operator < <()

~~~ .cpp
inline ostream& operator<<(ostream& os, const bbox2f& a);
~~~

stream write

#### Function operator \> \>()

~~~ .cpp
inline istream& operator>>(istream& is, bbox2f& a);
~~~

stream read

#### Function operator < <()

~~~ .cpp
inline ostream& operator<<(ostream& os, const bbox3f& a);
~~~

stream write

#### Function operator \> \>()

~~~ .cpp
inline istream& operator>>(istream& is, bbox3f& a);
~~~

stream read

#### Function operator < <()

~~~ .cpp
inline ostream& operator<<(ostream& os, const bbox4f& a);
~~~

stream write

#### Function operator \> \>()

~~~ .cpp
inline istream& operator>>(istream& is, bbox4f& a);
~~~

stream read

#### Function point_bbox()

~~~ .cpp
inline bbox3f point_bbox(const vec3f& p, float r = 0);
~~~

Point bounds

#### Function line_bbox()

~~~ .cpp
inline bbox3f line_bbox(
    const vec3f& v0, const vec3f& v1, float r0 = 0, float r1 = 0);
~~~

Line bounds

#### Function triangle_bbox()

~~~ .cpp
inline bbox3f triangle_bbox(const vec3f& v0, const vec3f& v1, const vec3f& v2);
~~~

Triangle bounds

#### Function quad_bbox()

~~~ .cpp
inline bbox3f quad_bbox(
    const vec3f& v0, const vec3f& v1, const vec3f& v2, const vec3f& v3);
~~~

Quad bounds

#### Function tetrahedron_bbox()

~~~ .cpp
inline bbox3f tetrahedron_bbox(
    const vec3f& v0, const vec3f& v1, const vec3f& v2, const vec3f& v3);
~~~

Tetrahedron bounds

#### Struct ray3f

~~~ .cpp
struct ray3f {
    static const int N = 3;
    using T = float;
    vec3f o;
    vec3f d;
    float tmin;
    float tmax;
    ray3f(); 
    ray3f(const vec3f& o, const vec3f& d, float tmin = 0, float tmax = flt_max); 
}
~~~

Rays with origin, direction and min/max t value.

- Members:
    - N:      size
    - T:      type
    - o:      origin
    - d:      direction
    - tmin:      minimum distance
    - tmax:      maximum distance
    - ray3f():      default constructor
    - ray3f():      initializes a ray from its elements


#### Function operator < <()

~~~ .cpp
inline ostream& operator<<(ostream& os, const ray3f& a);
~~~

stream write

#### Function operator \> \>()

~~~ .cpp
inline istream& operator>>(istream& is, ray3f& a);
~~~

stream read

#### Function transform_point()

~~~ .cpp
inline vec3f transform_point(const mat4f& a, const vec3f& b);
~~~

transforms a point by a matrix

#### Function transform_vector()

~~~ .cpp
inline vec3f transform_vector(const mat4f& a, const vec3f& b);
~~~

transforms a vector by a matrix

#### Function transform_direction()

~~~ .cpp
inline vec3f transform_direction(const mat4f& a, const vec3f& b);
~~~

transforms a direction by a matrix

#### Function transform_point()

~~~ .cpp
inline vec3f transform_point(const frame3f& a, const vec3f& b);
~~~

transforms a point by a frame (rigid affine transform)

#### Function transform_vector()

~~~ .cpp
inline vec3f transform_vector(const frame3f& a, const vec3f& b);
~~~

transforms a vector by a frame (rigid affine transform)

#### Function transform_direction()

~~~ .cpp
inline vec3f transform_direction(const frame3f& a, const vec3f& b);
~~~

transforms a direction by a frame (rigid affine transform)

#### Function transform_frame()

~~~ .cpp
inline frame3f transform_frame(const frame3f& a, const frame3f& b);
~~~

transforms a frame by a frame (rigid affine transform)

#### Function transform_point_inverse()

~~~ .cpp
inline vec3f transform_point_inverse(const frame3f& a, const vec3f& b);
~~~

inverse transforms a point by a frame (rigid affine transform)

#### Function transform_vector_inverse()

~~~ .cpp
inline vec3f transform_vector_inverse(const frame3f& a, const vec3f& b);
~~~

inverse transforms a vector by a frame (rigid affine transform)

#### Function transform_direction_inverse()

~~~ .cpp
inline vec3f transform_direction_inverse(const frame3f& a, const vec3f& b);
~~~

inverse transforms a direction by a frame (rigid affine transform)

#### Function transform_ray()

~~~ .cpp
inline ray3f transform_ray(const mat4f& a, const ray3f& b);
~~~

transforms a ray by a matrix (direction is not normalized after)

#### Function transform_bbox()

~~~ .cpp
inline bbox3f transform_bbox(const mat4f& a, const bbox3f& b);
~~~

transforms a bbox by a matrix

#### Function transform_ray()

~~~ .cpp
inline ray3f transform_ray(const frame3f& a, const ray3f& b);
~~~

transforms a ray by a frame (rigid affine transform)

#### Function transform_bbox()

~~~ .cpp
inline bbox3f transform_bbox(const frame3f& a, const bbox3f& b);
~~~

transforms a bbox by a frame (rigid affine transform)

#### Function transform_ray_inverse()

~~~ .cpp
inline ray3f transform_ray_inverse(const frame3f& a, const ray3f& b);
~~~

inverse transforms a ray by a frame (rigid affine transform)

#### Function transform_bbox_inverse()

~~~ .cpp
inline bbox3f transform_bbox_inverse(const frame3f& a, const bbox3f& b);
~~~

inverse transforms a bbox by a frame (rigid affine transform)

#### Function rotation_mat3f()

~~~ .cpp
inline mat3f rotation_mat3f(const vec3f& axis, float angle);
~~~

rotation matrix from axis-angle

#### Function translation_frame3f()

~~~ .cpp
inline frame3f translation_frame3f(const vec3f& a);
~~~

translation frame

#### Function translation_mat4f()

~~~ .cpp
inline mat4f translation_mat4f(const vec3f& a);
~~~

translation matrix

#### Function scaling_frame3f()

~~~ .cpp
inline frame3f scaling_frame3f(const vec3f& a);
~~~

scaling frame (this is not rigid and here for symmatry of API)

#### Function scaling_mat4f()

~~~ .cpp
inline mat4f scaling_mat4f(const vec3f& a);
~~~

scaling matrix

#### Function rotation_frame3f()

~~~ .cpp
inline frame3f rotation_frame3f(const vec3f& axis, float angle);
~~~

rotation frame

#### Function rotation_mat4f()

~~~ .cpp
inline mat4f rotation_mat4f(const mat3f& rot);
~~~

rotation matrix

#### Function rotation_mat4f()

~~~ .cpp
inline mat4f rotation_mat4f(const vec3f& axis, float angle);
~~~

rotation matrix

#### Function rotation_axisangle4()

~~~ .cpp
inline vec4f rotation_axisangle4(const quat4f& a);
~~~

quaternion axis-angle conversion

#### Function rotation_quat4f()

~~~ .cpp
inline quat4f rotation_quat4f(const vec4f& axis_angle);
~~~

axis-angle to quaternion

#### Function rotation_mat3f()

~~~ .cpp
inline mat3f rotation_mat3f(const quat4f& v);
~~~

quaterion to matrix conversion

#### Function rotation_mat4f()

~~~ .cpp
inline mat4f rotation_mat4f(const quat4f& v);
~~~

rotation matrix

#### Function rotation_quat4f()

~~~ .cpp
inline quat4f rotation_quat4f(const mat3f& m_);
~~~

matrix to quaternion

#### Function lookat_frame3f()

~~~ .cpp
inline frame3f lookat_frame3f(const vec3f& eye, const vec3f& center,
    const vec3f& up, bool inv_xz = false);
~~~

OpenGL lookat frame

#### Function lookat_mat4f()

~~~ .cpp
inline mat4f lookat_mat4f(
    const vec3f& eye, const vec3f& center, const vec3f& up);
~~~

OpenGL lookat matrix

#### Function frustum_mat4f()

~~~ .cpp
inline mat4f frustum_mat4f(
    float l, float r, float b, float t, float n, float f);
~~~

OpenGL frustum matrix

#### Function ortho_mat4f()

~~~ .cpp
inline mat4f ortho_mat4f(float l, float r, float b, float t, float n, float f);
~~~

OpenGL orthographic matrix

#### Function ortho2d_mat4f()

~~~ .cpp
inline mat4f ortho2d_mat4f(float left, float right, float bottom, float top);
~~~

OpenGL orthographic 2D matrix

#### Function ortho_mat4f()

~~~ .cpp
inline mat4f ortho_mat4f(float xmag, float ymag, float near, float far);
~~~

OpenGL/GLTF orthographic matrix

#### Function perspective_mat4f()

~~~ .cpp
inline mat4f perspective_mat4f(
    float fovy, float aspect, float near, float far);
~~~

OpenGL/GLTF perspective matrix

#### Function perspective_mat4f()

~~~ .cpp
inline mat4f perspective_mat4f(float fovy, float aspect, float near);
~~~

OpenGL/GLTF infinite perspective matrix

#### Function decompose_mat4f()

~~~ .cpp
inline void decompose_mat4f(
    const mat4f& m, vec3f& translation, mat3f& rotation, vec3f& scale);
~~~

Decompose an affine matrix into translation, rotation, scale.
Assumes there is no shear and the matrix is affine.

#### Function to_quat4f()

~~~ .cpp
inline quat4f to_quat4f(const mat3f& a);
~~~

Convert a rotation matrix to a quaternion

#### Function decompose_mat4f()

~~~ .cpp
inline void decompose_mat4f(
    const mat4f& m, vec3f& translation, quat4f& rotation, vec3f& scale);
~~~

Decompose an affine matrix into translation, rotation, scale.
Assumes there is no shear and the matrix is affine.

#### Function compose_mat4f()

~~~ .cpp
inline mat4f compose_mat4f(
    const vec3f& translation, const mat3f& rotation, const vec3f& scale);
~~~

Decompose an affine matrix into translation, rotation, scale.
Assumes there is no shear and the matrix is affine.

#### Function compose_mat4f()

~~~ .cpp
inline mat4f compose_mat4f(
    const vec3f& translation, const quat4f& rotation, const vec3f& scale);
~~~

Decompose an affine matrix into translation, rotation, scale.
Assumes there is no shear and the matrix is affine.

#### Function camera_turntable()

~~~ .cpp
inline void camera_turntable(vec3f& from, vec3f& to, vec3f& up,
    const vec3f& rotate, float dolly, const vec3f& pan);
~~~

Turntable for UI navigation from a from/to/up parametrization of the
camera.

#### Function camera_turntable()

~~~ .cpp
inline void camera_turntable(frame3f& frame, float& focus, const vec2f& rotate,
    float dolly, const vec2f& pan);
~~~

Turntable for UI navigation for a frame/distance parametrization of the
camera.

#### Function camera_fps()

~~~ .cpp
inline void camera_fps(
    frame3f& frame, const vec3f& transl, const vec2f& rotate);
~~~

FPS camera for UI navigation for a frame parametrization.
https://gamedev.stackexchange.com/questions/30644/how-to-keep-my-quaternion-using-fps-camera-from-tilting-and-messing-up

#### Struct rng_pcg32

~~~ .cpp
struct rng_pcg32 {
    uint64_t state = 0x853c49e6748fea9bULL;
    uint64_t inc = 0xda3e39cb94b95bdbULL;
}
~~~

PCG random numbers. A family of random number generators that supports
multiple sequences. In our code, we allocate one sequence for each sample.
PCG32 from http://www.pcg-random.org/

- Members:
    - state:      RNG state.
    - inc:      RNG sequence. Must be odd.


#### Function advance_rng()

~~~ .cpp
inline uint32_t advance_rng(rng_pcg32& rng);
~~~

Next random number

#### Function advance_rng()

~~~ .cpp
inline void advance_rng(rng_pcg32& rng, uint64_t delta);
~~~

Multi-step advance function (jump-ahead, jump-back).

#### Function advance_rng()

~~~ .cpp
inline void advance_rng(rng_pcg32& rng, int64_t delta);
~~~

Multi-step advance function (jump-ahead, jump-back).

#### Function seed_rng()

~~~ .cpp
inline void seed_rng(rng_pcg32& rng, uint64_t state, uint64_t seq = 1);
~~~

Seeds a random number generator with a state state from the sequence seq.

#### Function init_rng()

~~~ .cpp
inline rng_pcg32 init_rng(uint64_t state, uint64_t seq = 1);
~~~

Init a random number generator with a state state from the sequence seq.

#### Function next_rand1i()

~~~ .cpp
inline uint32_t next_rand1i(rng_pcg32& rng, uint32_t n);
~~~

Next random uint in [0,n) range with proper weighting

#### Function next_rand1f()

~~~ .cpp
inline float next_rand1f(rng_pcg32& rng);
~~~

Next random float in [0,1).

#### Function next_rand1f()

~~~ .cpp
inline float next_rand1f(rng_pcg32& rng, float a, float b);
~~~

Next random float in [a,b).

#### Function next_rand2f()

~~~ .cpp
inline vec2f next_rand2f(rng_pcg32& rng);
~~~

Next random float2 in [0,1)x[0,1).

#### Function next_rand2f()

~~~ .cpp
inline vec2f next_rand2f(rng_pcg32& rng, const vec2f& a, const vec2f& b);
~~~

Next random float in [a.x,b.x)x[a.y,b.y).

#### Function next_rand3f()

~~~ .cpp
inline vec3f next_rand3f(rng_pcg32& rng);
~~~

Next random float3 in [0,1)x[0,1)x[0,1).

#### Function next_rand2f()

~~~ .cpp
inline vec3f next_rand2f(rng_pcg32& rng, const vec3f& a, const vec3f& b);
~~~

Next random float in [a.x,b.x)x[a.y,b.y)x[a.z,b.z).

#### Function next_rand1d()

~~~ .cpp
inline double next_rand1d(rng_pcg32& rng);
~~~

Next random double in [0, 1). Only 32 mantissa bits are filled, but still
better than float that uses 23.

#### Function rng_distance()

~~~ .cpp
inline int64_t rng_distance(const rng_pcg32& a, const rng_pcg32& b);
~~~

Distance between random number generators

#### Function rng_shuffle()

~~~ .cpp
template <typename Iterator>
inline void rng_shuffle(rng_pcg32& rng, Iterator begin, Iterator end);
~~~

Random shuffle of a sequence.

#### Function rng_shuffle()

~~~ .cpp
template <typename T>
inline void rng_shuffle(rng_pcg32& rng, T* vals, int num);
~~~

Random shuffle of a sequence.

#### Function rng_shuffle()

~~~ .cpp
template <typename T>
inline void rng_shuffle(rng_pcg32& rng, vector<T>& vals);
~~~

Random shuffle of a sequence.

#### Function operator==()

~~~ .cpp
inline bool operator==(const rng_pcg32& a, const rng_pcg32& b);
~~~

Equality operator

#### Function operator!=()

~~~ .cpp
inline bool operator!=(const rng_pcg32& a, const rng_pcg32& b);
~~~

Inequality operator

#### Function sample_hemisphere()

~~~ .cpp
inline vec3f sample_hemisphere(const vec2f& ruv);
~~~

sample hemispherical direction with uniform distribution

#### Function sample_hemisphere_pdf()

~~~ .cpp
inline float sample_hemisphere_pdf(const vec3f& w);
~~~

pdf for hemispherical direction with uniform distribution

#### Function sample_sphere()

~~~ .cpp
inline vec3f sample_sphere(const vec2f ruv);
~~~

spherical direction with uniform distribution

#### Function sample_sphere_pdf()

~~~ .cpp
inline float sample_sphere_pdf(const vec3f& w);
~~~

pdf for spherical direction with uniform distribution

#### Function sample_hemisphere_cosine()

~~~ .cpp
inline vec3f sample_hemisphere_cosine(const vec2f& ruv);
~~~

hemispherical direction with cosine distribution

#### Function sample_hemisphere_cosine_pdf()

~~~ .cpp
inline float sample_hemisphere_cosine_pdf(const vec3f& w);
~~~

pdf for hemispherical direction with cosine distribution

#### Function sample_hemisphere_cospower()

~~~ .cpp
inline vec3f sample_hemisphere_cospower(const vec2f& ruv, float n);
~~~

hemispherical direction with cosine power distribution

#### Function sample_hemisphere_cospower_pdf()

~~~ .cpp
inline float sample_hemisphere_cospower_pdf(const vec3f& w, float n);
~~~

pdf for hemispherical direction with cosine power distribution

#### Function sample_disk()

~~~ .cpp
inline vec3f sample_disk(const vec2f& ruv);
~~~

uniform disk

#### Function sample_disk_pdf()

~~~ .cpp
inline float sample_disk_pdf();
~~~

pdf for uniform disk

#### Function sample_cylinder()

~~~ .cpp
inline vec3f sample_cylinder(const vec2f& ruv);
~~~

uniform cylinder

#### Function sample_cylinder_pdf()

~~~ .cpp
inline float sample_cylinder_pdf();
~~~

pdf for uniform cylinder

#### Function sample_triangle()

~~~ .cpp
inline vec2f sample_triangle(const vec2f& ruv);
~~~

uniform triangle

#### Function sample_triangle()

~~~ .cpp
inline vec3f sample_triangle(
    const vec2f& ruv, const vec3f& v0, const vec3f& v1, const vec3f& v2);
~~~

uniform triangle

#### Function sample_triangle_pdf()

~~~ .cpp
inline float sample_triangle_pdf(
    const vec3f& v0, const vec3f& v1, const vec3f& v2);
~~~

pdf for uniform triangle (triangle area)

#### Function sample_index()

~~~ .cpp
inline int sample_index(float r, int size);
~~~

index with uniform distribution

#### Function sample_index_pdf()

~~~ .cpp
inline float sample_index_pdf(int size);
~~~

pdf for index with uniform distribution

#### Function hash_permute()

~~~ .cpp
inline uint32_t hash_permute(uint32_t i, uint32_t n, uint32_t key);
~~~

Computes the i-th term of a permutation of l values keyed by p.
From Correlated Multi-Jittered Sampling by Kensler @ Pixar

#### Function hash_randfloat()

~~~ .cpp
inline float hash_randfloat(uint32_t i, uint32_t key);
~~~

Computes a float value by hashing i with a key p.
From Correlated Multi-Jittered Sampling by Kensler @ Pixar

#### Function hash_uint32()

~~~ .cpp
inline uint32_t hash_uint32(uint64_t a);
~~~

32 bit integer hash. Public domain code.

#### Function hash_uint64()

~~~ .cpp
inline uint64_t hash_uint64(uint64_t a);
~~~

64 bit integer hash. Public domain code.

#### Function hash_uint64_32()

~~~ .cpp
inline uint32_t hash_uint64_32(uint64_t a);
~~~

64-to-32 bit integer hash. Public domain code.

#### Function hash_combine()

~~~ .cpp
inline size_t hash_combine(size_t a, size_t b);
~~~

Combines two 64 bit hashes as in boost::hash_combine

#### Function perlin_noise()

~~~ .cpp
float perlin_noise(const vec3f& p, const vec3i& wrap = zero3i);
~~~

Compute the revised Pelin noise function. Wrap provides a wrapping noise
but must be power of two (wraps at 256 anyway). For octave based noise,
good values are obtained with octaves=6 (numerber of noise calls),
lacunarity=~2.0 (spacing between successive octaves: 2.0 for warpping
output), gain=0.5 (relative weighting applied to each successive octave),
offset=1.0 (used to invert the ridges).

#### Function perlin_ridge_noise()

~~~ .cpp
float perlin_ridge_noise(const vec3f& p, float lacunarity = 2.0f,
    float gain = 0.5f, float offset = 1.0f, int octaves = 6,
    const vec3i& wrap = zero3i);
~~~

Ridge noise function

#### Function perlin_fbm_noise()

~~~ .cpp
float perlin_fbm_noise(const vec3f& p, float lacunarity = 2.0f,
    float gain = 0.5f, int octaves = 6, const vec3i& wrap = zero3i);
~~~

Fractal brownian motion noise - see perlin_noise() for params.

#### Function perlin_turbulence_noise()

~~~ .cpp
float perlin_turbulence_noise(const vec3f& p, float lacunarity = 2.0f,
    float gain = 0.5f, int octaves = 6, const vec3i& wrap = zero3i);
~~~

Fractal turbulence noise - see perlin_noise() for params.

#### Function range()

~~~ .cpp
inline range_generator range(int max);
~~~

Python-like range

#### Function range()

~~~ .cpp
inline range_generator range(int min, int max, int step = 1);
~~~

Python-like range

#### Function enumerate()

~~~ .cpp
template <typename T>
inline enumerate_generator<const T> enumerate(const vector<T>& vv);
~~~

Python-like range

#### Function enumerate()

~~~ .cpp
template <typename T>
inline enumerate_generator<T> enumerate(vector<T>& vv);
~~~

Python-like range

#### Function operator+()

~~~ .cpp
template <typename T>
inline vector<T> operator+(const vector<T>& v, const T& vv);
~~~

Append an element to a vector

#### Function operator+=()

~~~ .cpp
template <typename T>
inline vector<T>& operator+=(vector<T>& v, const T& vv);
~~~

Append an element to a vector

#### Function operator+()

~~~ .cpp
template <typename T, typename ET>
inline vector<T> operator+(const vector<T>& v, const ET& vv);
~~~

Append an element to a vector

#### Function operator+=()

~~~ .cpp
template <typename T, typename ET>
inline vector<T>& operator+=(vector<T>& v, const ET& vv);
~~~

Append an element to a vector

#### Function operator+()

~~~ .cpp
template <typename T>
inline vector<T> operator+(const vector<T>& v, const vector<T>& vv);
~~~

Append a vector to a vector

#### Function operator+=()

~~~ .cpp
template <typename T>
inline vector<T>& operator+=(vector<T>& v, const vector<T>& vv);
~~~

Append a vector to a vector

#### Function get_key()

~~~ .cpp
template <typename K, typename V>
inline K get_key(const std::vector<pair<K, V>>& kvs, const V& v);
~~~

Get a key

#### Function get_value()

~~~ .cpp
template <typename K, typename V>
inline V get_value(const std::vector<pair<K, V>>& kvs, const K& k);
~~~

Get a value

#### Function find_idx()

~~~ .cpp
template <typename T>
inline int find_idx(const vector<T>& v, const T& vv);
~~~

Find a value in an array

#### Function contains()

~~~ .cpp
template <typename T>
inline bool contains(const vector<T>& v, const T& vv);
~~~

Checks if a containers contains a value

#### Function contains()

~~~ .cpp
template <typename K, typename V>
inline bool contains(const map<K, V>& v, const K& vv);
~~~

Checks if a containers contains a value

#### Function contains()

~~~ .cpp
template <typename K, typename V>
inline bool contains(const unordered_map<K, V>& v, const K& vv);
~~~

Checks if a containers contains a value

#### Function line_tangent()

~~~ .cpp
inline vec3f line_tangent(const vec3f& v0, const vec3f& v1);
~~~

line tangent

#### Function line_length()

~~~ .cpp
inline float line_length(const vec3f& v0, const vec3f& v1);
~~~

line length

#### Function triangle_normal()

~~~ .cpp
inline vec3f triangle_normal(
    const vec3f& v0, const vec3f& v1, const vec3f& v2);
~~~

triangle normal

#### Function triangle_area()

~~~ .cpp
inline float triangle_area(const vec3f& v0, const vec3f& v1, const vec3f& v2);
~~~

triangle area

#### Function quad_area()

~~~ .cpp
inline float quad_area(
    const vec3f& v0, const vec3f& v1, const vec3f& v2, const vec3f& v3);
~~~

quad area

#### Function tetrahedron_volume()

~~~ .cpp
inline float tetrahedron_volume(
    const vec3f& v0, const vec3f& v1, const vec3f& v2, const vec3f& v3);
~~~

tetrahedron volume

#### Function eval_barycentric_point()

~~~ .cpp
template <typename T>
inline T eval_barycentric_point(const vector<T>& vals, const int& p, float w);
~~~

line barycentric interpolation

#### Function eval_barycentric_line()

~~~ .cpp
template <typename T>
inline T eval_barycentric_line(
    const vector<T>& vals, const vec2i& l, const vec2f& w);
~~~

line barycentric interpolation

#### Function eval_barycentric_triangle()

~~~ .cpp
template <typename T>
inline T eval_barycentric_triangle(
    const vector<T>& vals, const vec3i& t, const vec3f& w);
~~~

triangle barycentric interpolation

#### Function eval_barycentric_tetra()

~~~ .cpp
template <typename T>
inline T eval_barycentric_tetra(
    const vector<T>& vals, const vec4i& t, const vec4f& w);
~~~

tetrahedron barycentric interpolation

#### Function eval_barycentric_quad()

~~~ .cpp
template <typename T>
inline T eval_barycentric_quad(
    const vector<T>& vals, const vec4i& t, const vec4f& w);
~~~

quad interpolation based on the two-triangle representation

#### Function eval_bernstein()

~~~ .cpp
template <typename T>
inline T eval_bernstein(T u, int i, int degree);
~~~

bernstein polynomials (for Bezier)

#### Function eval_bernstein_derivative()

~~~ .cpp
template <typename T>
inline T eval_bernstein_derivative(T u, int i, int degree);
~~~

bernstein polynomials (for Bezier)

#### Function eval_bezier_cubic()

~~~ .cpp
template <typename T, typename T1>
inline T eval_bezier_cubic(
    const T& v0, const T& v1, const T& v2, const T& v3, T1 t);
~~~

eval bezier

#### Function eval_bezier_cubic()

~~~ .cpp
template <typename T, typename T1>
inline T eval_bezier_cubic(const vector<T>& vals, const vec4i& b, T1 t);
~~~

eval bezier

#### Function eval_bezier_cubic_derivative()

~~~ .cpp
template <typename T, typename T1>
inline T eval_bezier_cubic_derivative(
    const T& v0, const T& v1, const T& v2, const T& v3, T1 t);
~~~

eval bezier derivative

#### Function eval_bezier_cubic_derivative()

~~~ .cpp
template <typename T, typename T1>
inline T eval_bezier_cubic_derivative(
    const vector<T>& vals, const vec4i& b, T1 t);
~~~

eval bezier derivative

#### Function compute_normals()

~~~ .cpp
inline vector<vec3f> compute_normals(const vector<vec2i>& lines,
    const vector<vec3i>& triangles, const vector<vec4i>& quads,
    const vector<vec3f>& pos, bool weighted = true);
~~~

Compute per-vertex normals/tangents for lines, triangles and quads with
positions pos. Weighted indicated whether the normals/tangents are
weighted by line length.

#### Function compute_tangent_frames()

~~~ .cpp
inline vector<vec4f> compute_tangent_frames(const vector<vec3i>& triangles,
    const vector<vec3f>& pos, const vector<vec3f>& norm,
    const vector<vec2f>& texcoord, bool weighted = true);
~~~

Compute per-vertex tangent frame for triangle meshes.
Tangent space is defined by a four component vector.
The first three components are the tangent with respect to the U texcoord.
The fourth component is the sign of the tangent wrt the V texcoord.
Tangent frame is useful in normal mapping.

#### Function compute_skinning()

~~~ .cpp
inline void compute_skinning(const vector<vec3f>& pos,
    const vector<vec3f>& norm, const vector<vec4f>& weights,
    const vector<vec4i>& joints, const vector<mat4f>& xforms,
    vector<vec3f>& skinned_pos, vector<vec3f>& skinned_norm);
~~~

Apply skinning

#### Function compute_skinning()

~~~ .cpp
inline void compute_skinning(const vector<vec3f>& pos,
    const vector<vec3f>& norm, const vector<vec4f>& weights,
    const vector<vec4i>& joints, const vector<frame3f>& xforms,
    vector<vec3f>& skinned_pos, vector<vec3f>& skinned_norm);
~~~

Apply skinning

#### Function compute_matrix_skinning()

~~~ .cpp
inline void compute_matrix_skinning(const vector<vec3f>& pos,
    const vector<vec3f>& norm, const vector<vec4f>& weights,
    const vector<vec4i>& joints, const vector<mat4f>& xforms,
    vector<vec3f>& skinned_pos, vector<vec3f>& skinned_norm);
~~~

Apply skinning as specified in Khronos glTF

#### Function get_edges()

~~~ .cpp
inline vector<vec2i> get_edges(const vector<vec2i>& lines,
    const vector<vec3i>& triangles, const vector<vec4i>& quads);
~~~

Create an array of edges.

#### Function get_boundary_edges()

~~~ .cpp
inline vector<vec2i> get_boundary_edges(const vector<vec2i>& lines,
    const vector<vec3i>& triangles, const vector<vec4i>& quads);
~~~

Create an array of boundary edges. Lines are always considered boundaries.

#### Function get_verts()

~~~ .cpp
inline vector<int> get_verts(const vector<vec2i>& lines,
    const vector<vec3i>& triangles, const vector<vec4i>& quads);
~~~

Get a list of all unique vertices.

#### Function get_boundary_verts()

~~~ .cpp
inline vector<int> get_boundary_verts(const vector<vec2i>& lines,
    const vector<vec3i>& triangles, const vector<vec4i>& quads);
~~~

Create an array of boundary vertices. Lines are always considered
boundaries.

#### Function convert_quads_to_triangles()

~~~ .cpp
inline vector<vec3i> convert_quads_to_triangles(const vector<vec4i>& quads);
~~~

Convert quads to triangles

#### Function convert_quads_to_triangles()

~~~ .cpp
inline vector<vec3i> convert_quads_to_triangles(
    const vector<vec4i>& quads, int row_length);
~~~

Convert quads to triangles with a diamond-like topology.
Quads have to be consecutive one row after another.

#### Function convert_face_varying()

~~~ .cpp
inline tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>
convert_face_varying(const vector<vec4i>& quads_pos,
    const vector<vec4i>& quads_norm, const vector<vec4i>& quads_texcoord,
    const vector<vec3f>& pos, const vector<vec3f>& norm,
    const vector<vec2f>& texcoord);
~~~

Convert face varying data to single primitives. Returns the quads indices
and filled vectors for pos, norm and texcoord.

#### Function subdivide_elems()

~~~ .cpp
inline tuple<vector<vec2i>, vector<vec3i>, vector<vec4i>, vector<vec2i>,
    vector<vec4i>>
subdivide_elems(const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads, int nverts);
~~~

Tesselate lines, triangles and quads by spolitting edges.
Returns the tesselated elements and dictionaries for vertex calculations.

#### Function subdivide_vert()

~~~ .cpp
template <typename T>
inline vector<T> subdivide_vert(const vector<T>& vert,
    const vector<vec2i>& edges, const vector<vec4i>& faces,
    bool normalized = false);
~~~

Subdivide vertex properties given the maps

#### Function subdivide_catmullclark()

~~~ .cpp
template <typename T>
inline vector<T> subdivide_catmullclark(const vector<vec4i>& quads,
    const vector<T>& vert, const vector<vec2i>& crease_tlines,
    const vector<int>& crease_tpoints, bool normalized = false);
~~~

Performs the smoothing step of Catmull-Clark. Start with a tesselate quad
mesh obtained with subdivide_elems() and subdivide_vert(). To handle open
meshes with boundary, get the boundary from make_boundary_edge() and pass it
as crease_lines. To fix the boundary entirely, just get the boundary
vertices and pass it as creases.

#### Function make_uvquads()

~~~ .cpp
inline tuple<vector<vec4i>, vector<vec2f>> make_uvquads(int usteps, int vsteps,
    bool uwrap = false, bool vwrap = false, bool vpole0 = false,
    bool vpole1 = false);
~~~

Generate a rectangular grid of usteps x vsteps uv values for parametric
surface generation.

#### Function make_uvlines()

~~~ .cpp
inline tuple<vector<vec2i>, vector<vec2f>> make_uvlines(int num, int usteps);
~~~

Generate parametric num lines of usteps segments.

#### Function make_uvpoints()

~~~ .cpp
inline tuple<vector<int>, vector<vec2f>> make_uvpoints(int num);
~~~

Generate a parametric point set. Mostly here for completeness.

#### Function merge_elems()

~~~ .cpp
inline tuple<vector<vec2i>, vector<vec3i>, vector<vec4i>> merge_elems(
    int nverts, const vector<vec2i>& lines1, const vector<vec3i>& triangles1,
    const vector<vec4i>& quads1, const vector<vec2i>& lines2,
    const vector<vec3i>& triangles2, const vector<vec4i>& quads2);
~~~

Merge elements between shapes. The elements are merged by increasing the
array size of the second array by the number of vertices of the first.
Vertex data can then be concatenated successfully.

#### Function facet_elems()

~~~ .cpp
inline tuple<vector<vec2i>, vector<vec3i>, vector<vec4i>, vector<int>>
facet_elems(const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads);
~~~

Unshare shape data by duplicating all vertex data for each element,
giving a faceted look. Note that faceted tangents are not computed.

#### Function facet_vert()

~~~ .cpp
template <typename T>
inline vector<T> facet_vert(const vector<T>& vert, const vector<int>& vmap);
~~~

Unshare vertices for faceting

#### Function sample_points()

~~~ .cpp
inline int sample_points(int npoints, float re);
~~~

Pick a point

#### Function sample_points_cdf()

~~~ .cpp
inline vector<float> sample_points_cdf(int npoints);
~~~

Compute a distribution for sampling points uniformly

#### Function sample_points()

~~~ .cpp
inline int sample_points(const vector<float>& cdf, float re);
~~~

Pick a point

#### Function sample_lines_cdf()

~~~ .cpp
inline vector<float> sample_lines_cdf(
    const vector<vec2i>& lines, const vector<vec3f>& pos);
~~~

Compute a distribution for sampling lines uniformly

#### Function sample_lines()

~~~ .cpp
inline pair<int, vec2f> sample_lines(
    const vector<float>& cdf, float re, float ruv);
~~~

Pick a point on lines

#### Function sample_triangles_cdf()

~~~ .cpp
inline vector<float> sample_triangles_cdf(
    const vector<vec3i>& triangles, const vector<vec3f>& pos);
~~~

Compute a distribution for sampling triangle meshes uniformly

#### Function sample_triangles()

~~~ .cpp
inline pair<int, vec3f> sample_triangles(
    const vector<float>& cdf, float re, const vec2f& ruv);
~~~

Pick a point on a triangle mesh

#### Function sample_quads_cdf()

~~~ .cpp
inline vector<float> sample_quads_cdf(
    const vector<vec4i>& quads, const vector<vec3f>& pos);
~~~

Compute a distribution for sampling quad meshes uniformly

#### Function sample_quads()

~~~ .cpp
inline pair<int, vec4f> sample_quads(
    const vector<float>& cdf, float re, const vec2f& ruv);
~~~

Pick a point on a quad mesh

#### Function sample_triangles_points()

~~~ .cpp
inline tuple<vector<vec3f>, vector<vec3f>, vector<vec2f>>
sample_triangles_points(const vector<vec3i>& triangles,
    const vector<vec3f>& pos, const vector<vec3f>& norm,
    const vector<vec2f>& texcoord, int npoints, uint64_t seed = 0);
~~~

Samples a set of points over a triangle mesh uniformly. The rng function
takes the point index and returns vec3f numbers uniform directibuted in
[0,1]^3. unorm and texcoord are optional.

#### Function make_sphere()

~~~ .cpp
tuple<vector<vec3i>, vector<vec3f>> make_sphere(int level);
~~~

Make a sphere. Returns quads, pos.

#### Function make_geodesicsphere()

~~~ .cpp
tuple<vector<vec3i>, vector<vec3f>> make_geodesicsphere(int level);
~~~

Make a geodesic sphere. Returns quads, pos.

#### Function make_cube()

~~~ .cpp
tuple<vector<vec4i>, vector<vec3f>> make_cube();
~~~

Make a cube with unique vertices. This is watertight but has no
texture coordinates or normals. Returns quads, pos.

#### Function make_uvsphere()

~~~ .cpp
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>> make_uvsphere(
    int level, bool flipped = false);
~~~

Make a sphere. This is not watertight. Returns quads, pos, norm, texcoord.

#### Function make_uvhemisphere()

~~~ .cpp
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>
make_uvhemisphere(int level, bool flipped = false);
~~~

Make a sphere. This is not watertight. Returns quads, pos, norm, texcoord.

#### Function make_uvquad()

~~~ .cpp
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>> make_uvquad(
    int level);
~~~

Make a quad. Returns quads, pos, norm, texcoord.

#### Function make_fvsphere()

~~~ .cpp
tuple<vector<vec4i>, vector<vec3f>, vector<vec4i>, vector<vec3f>, vector<vec4i>,
    vector<vec2f>>
make_fvsphere();
~~~

Make a facevarying sphere with unique vertices but different texture
coordinates. Returns (quads, pos), (quads, norm), (quads, texcoord).

#### Function make_fvcube()

~~~ .cpp
tuple<vector<vec4i>, vector<vec3f>, vector<vec4i>, vector<vec3f>, vector<vec4i>,
    vector<vec2f>>
make_fvcube();
~~~

Make a facevarying cube with unique vertices but different texture
coordinates. Returns (quads, pos), (quads, norm), (quads, texcoord).

#### Function make_suzanne()

~~~ .cpp
tuple<vector<vec4i>, vector<vec3f>> make_suzanne();
~~~

Make a suzanne monkey model for testing. Note that some quads are
degenerate. Returns quads, pos.

#### Function make_uvcube()

~~~ .cpp
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>> make_uvcube(
    int level);
~~~

Make a cube with uv. This is not watertight. Returns quads, pos, norm,
texcoord.

#### Function make_uvspherecube()

~~~ .cpp
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>
make_uvspherecube(int level);
~~~

Make a sphere from a cube. This is not watertight. Returns quads, pos, norm,
texcoord.

#### Function make_uvspherizedcube()

~~~ .cpp
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>
make_uvspherizedcube(int level, float radius);
~~~

Make a cube than stretch it towards a sphere. This is not watertight.
Returns quads, pos, norm, texcoord.

#### Function make_uvflipcapsphere()

~~~ .cpp
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>
make_uvflipcapsphere(int level, float z, bool flipped = false);
~~~

Make a flipped sphere. This is not watertight. Returns quads, pos, norm,
texcoord.

#### Function make_uvcutsphere()

~~~ .cpp
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>
make_uvcutsphere(int level, float z, bool flipped = false);
~~~

Make a cutout sphere. This is not watertight. Returns quads, pos, norm,
texcoord.

#### Function make_hair()

~~~ .cpp
tuple<vector<vec2i>, vector<vec3f>, vector<vec3f>, vector<vec2f>, vector<float>>
make_hair(int num, int level, const vec2f& len, const vec2f& rad,
    const vector<vec3i>& striangles, const vector<vec4i>& squads,
    const vector<vec3f>& spos, const vector<vec3f>& snorm,
    const vector<vec2f>& stexcoord, const vec2f& noise = zero2f,
    const vec2f& clump = zero2f, const vec2f& rotation = zero2f,
    uint32_t seed = 0);
~~~

Make a hair ball around a shape. Returns lines, pos, norm, texcoord, radius.

#### Struct image4f

~~~ .cpp
struct image4f {
    image4f(); 
    image4f(int w, int h, const vec4f& v = zero4f); 
    image4f(int w, int h, const vec4f* v); 
    int width() const; 
    int height() const; 
    vec2i size() const; 
    bool empty() const; 
    explicit operator bool() const; 
    void resize(int w, int h, const vec4f& v = zero4f); 
    void assign(int w, int h, const vec4f& v); 
    void set(const vec4f& v); 
    vec4f& operator[](const vec2i& ij); 
    const vec4f& operator[](const vec2i& ij) const; 
    vec4f& at(const vec2i& ij); 
    const vec4f& at(const vec2i& ij) const; 
    vec4f& at(int i, int j); 
    const vec4f& at(int i, int j) const; 
    vec4f* data(); 
    const vec4f* data() const; 
}
~~~

HDR image

- Members:
    - image4f():      empty image constructor
    - image4f():      image constructor
    - image4f():      image constructor
    - width():      width
    - height():      height
    - size():      size
    - empty():      check for empty
    - operator bool():      check for empty
    - resize():      reallocate memory
    - assign():      reallocate memory
    - set():      set values
    - operator[]():      element access
    - operator[]():      element access
    - at():      element access
    - at():      element access
    - at():      element access
    - at():      element access
    - data():      data access
    - data():      data access


#### Struct image4b

~~~ .cpp
struct image4b {
    image4b(); 
    image4b(int w, int h, const vec4b& v = zero4b); 
    image4b(int w, int h, const vec4b* v); 
    int width() const; 
    int height() const; 
    vec2i size() const; 
    bool empty() const; 
    explicit operator bool() const; 
    void resize(int w, int h, const vec4b& v = zero4b); 
    void assign(int w, int h, const vec4b& v); 
    void set(const vec4b& v); 
    vec4b& operator[](const vec2i& ij); 
    const vec4b& operator[](const vec2i& ij) const; 
    vec4b& at(const vec2i& ij); 
    const vec4b& at(const vec2i& ij) const; 
    vec4b& at(int i, int j); 
    const vec4b& at(int i, int j) const; 
    vec4b* data(); 
    const vec4b* data() const; 
}
~~~

LDR image

- Members:
    - image4b():      empty image constructor
    - image4b():      image constructor
    - image4b():      image constructor
    - width():      width
    - height():      height
    - size():      size
    - empty():      check for empty
    - operator bool():      check for empty
    - resize():      reallocate memory
    - assign():      reallocate memory
    - set():      set values
    - operator[]():      element access
    - operator[]():      element access
    - at():      element access
    - at():      element access
    - at():      element access
    - at():      element access
    - data():      data access
    - data():      data access


#### Function srgb_to_linear()

~~~ .cpp
inline vec3f srgb_to_linear(const vec3b& srgb);
~~~

Approximate conversion from srgb.

#### Function srgb_to_linear()

~~~ .cpp
inline vec4f srgb_to_linear(const vec4b& srgb);
~~~

Approximate conversion from srgb.

#### Function linear_to_srgb()

~~~ .cpp
inline vec3b linear_to_srgb(const vec3f& lin);
~~~

Approximate conversion to srgb.

#### Function linear_to_srgb()

~~~ .cpp
inline vec4b linear_to_srgb(const vec4f& lin);
~~~

Approximate conversion to srgb.

#### Function xyz_to_xyY()

~~~ .cpp
inline vec3f xyz_to_xyY(const vec3f& xyz);
~~~

Convert between CIE XYZ and xyY

#### Function xyY_to_xyz()

~~~ .cpp
inline vec3f xyY_to_xyz(const vec3f& xyY);
~~~

Convert between CIE XYZ and xyY

#### Function xyz_to_rgb()

~~~ .cpp
inline vec3f xyz_to_rgb(const vec3f& xyz);
~~~

Convert between CIE XYZ and RGB

#### Function rgb_to_xyz()

~~~ .cpp
inline vec3f rgb_to_xyz(const vec3f& rgb);
~~~

Convert between CIE XYZ and RGB

#### Function tonemap_filmic()

~~~ .cpp
inline float tonemap_filmic(float hdr);
~~~

Tone map with a fitted filmic curve.

Implementation from
https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/

#### Function tonemap_image()

~~~ .cpp
inline image4b tonemap_image(
    const image4f& hdr, float exposure, float gamma, bool filmic = false);
~~~

Tone mapping HDR to LDR images.

#### Function image_over()

~~~ .cpp
inline void image_over(
    vec4f* img, int width, int height, int nlayers, vec4f** layers);
~~~

Image over operator

#### Function image_over()

~~~ .cpp
inline void image_over(
    vec4b* img, int width, int height, int nlayers, vec4b** layers);
~~~

Image over operator

#### Function hsv_to_rgb()

~~~ .cpp
inline vec4b hsv_to_rgb(const vec4b& hsv);
~~~

Convert HSV to RGB

Implementatkion from
http://stackoverflow.com/questions/3018313/algorithm-to-convert-rgb-to-hsv-and-hsv-to-rgb-in-range-0-255-for-both

#### Function make_grid_image()

~~~ .cpp
inline image4b make_grid_image(int width, int height, int tile = 64,
    const vec4b& c0 =;
~~~

Make a grid image

#### Function make_checker_image()

~~~ .cpp
inline image4b make_checker_image(int width, int height, int tile = 64,
    const vec4b& c0 =;
~~~

Make a checkerboard image

#### Function make_bumpdimple_image()

~~~ .cpp
inline image4b make_bumpdimple_image(int width, int height, int tile = 64);
~~~

Make an image with bumps and dimples.

#### Function make_ramp_image()

~~~ .cpp
inline image4b make_ramp_image(int width, int height, const vec4b& c0,
    const vec4b& c1, bool srgb = false);
~~~

Make a uv colored grid

#### Function make_gammaramp_image()

~~~ .cpp
inline image4b make_gammaramp_image(int width, int height);
~~~

Make a gamma ramp image

#### Function make_gammaramp_imagef()

~~~ .cpp
inline image4f make_gammaramp_imagef(int width, int height);
~~~

Make a gamma ramp image

#### Function make_uv_image()

~~~ .cpp
inline image4b make_uv_image(int width, int height);
~~~

Make an image color with red/green in the [0,1] range. Helpful to visualize
uv texture coordinate application.

#### Function make_uvgrid_image()

~~~ .cpp
inline image4b make_uvgrid_image(
    int width, int height, int tile = 64, bool colored = true);
~~~

Make a uv colored grid

#### Function make_recuvgrid_image()

~~~ .cpp
inline image4b make_recuvgrid_image(
    int width, int height, int tile = 64, bool colored = true);
~~~

Make a uv recusive colored grid

#### Function bump_to_normal_map()

~~~ .cpp
inline image4b bump_to_normal_map(const image4b& img, float scale = 1);
~~~

Comvert a bump map to a normal map.

#### Function make_sunsky_image()

~~~ .cpp
image4f make_sunsky_image(int res, float thetaSun, float turbidity = 3,
    bool has_sun = false, bool has_ground = true);
~~~

Make a sunsky HDR model with sun at theta elevation in [0,pi/2], turbidity
in [1.7,10] with or without sun.

#### Function make_noise_image()

~~~ .cpp

image4b make_noise_image(int resx, int resy, float scale = 1, bool wrap = true);
~~~

Compute the revised Pelin noise function. Wrap provides a wrapping noise
but must be power of two (wraps at 256 anyway). For octave based noise,
good values are obtained with octaves=6 (numerber of noise calls),
lacunarity=~2.0 (spacing between successive octaves: 2.0 for warpping
output), gain=0.5 (relative weighting applied to each successive octave),
offset=1.0 (used to invert the ridges).
Make a noise image. Wrap works only if both resx and resy are powers of two.

#### Function make_fbm_image()

~~~ .cpp
image4b make_fbm_image(int resx, int resy, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, int octaves = 6, bool wrap = true);
~~~

Make a noise image. Wrap works only if both resx and resy are powers of two.

#### Function make_ridge_image()

~~~ .cpp
image4b make_ridge_image(int resx, int resy, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, float offset = 1.0f,
    int octaves = 6, bool wrap = true);
~~~

Make a noise image. Wrap works only if both resx and resy are powers of two.

#### Function make_turbulence_image()

~~~ .cpp
image4b make_turbulence_image(int resx, int resy, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, int octaves = 6, bool wrap = true);
~~~

Make a noise image. Wrap works only if both resx and resy are powers of two.

#### Function is_hdr_filename()

~~~ .cpp
bool is_hdr_filename(const string& filename);
~~~

Check if an image is HDR based on filename

#### Function load_image4b()

~~~ .cpp
image4b load_image4b(const string& filename);
~~~

Loads an ldr image.

#### Function load_image4f()

~~~ .cpp
image4f load_image4f(const string& filename);
~~~

Loads an hdr image.

#### Function save_image4b()

~~~ .cpp
bool save_image4b(const string& filename, const image4b& img);
~~~

Saves an ldr image.

#### Function save_image4f()

~~~ .cpp
bool save_image4f(const string& filename, const image4f& img);
~~~

Saves an hdr image.

#### Function load_imagef()

~~~ .cpp
vector<float> load_imagef(
    const string& filename, int& width, int& height, int& ncomp);
~~~

Loads an image

#### Function load_image()

~~~ .cpp
vector<byte> load_image(
    const string& filename, int& width, int& height, int& ncomp);
~~~

Loads an image

#### Function load_imagef_from_memory()

~~~ .cpp
vector<float> load_imagef_from_memory(const string& filename, const byte* data,
    int length, int& width, int& height, int& ncomp);
~~~

Loads an image from memory.

#### Function load_image_from_memory()

~~~ .cpp
vector<byte> load_image_from_memory(const string& filename, const byte* data,
    int length, int& width, int& height, int& ncomp);
~~~

Loads an image from memory.

#### Function save_imagef()

~~~ .cpp
bool save_imagef(
    const string& filename, int width, int height, int ncomp, const float* hdr);
~~~

Saves an image

#### Function save_image()

~~~ .cpp
bool save_image(
    const string& filename, int width, int height, int ncomp, const byte* ldr);
~~~

Saves an image

#### Function save_image()

~~~ .cpp
inline bool save_image(const string& filename, const image4f& hdr,
    float exposure, float gamma, bool filmic = false);
~~~

Save an HDR or LDR image with tonemapping based on filename

#### Enum resize_filter

~~~ .cpp
enum struct resize_filter {
    def = 0,
    box = 1,
    triangle = 2,
    cubic_spline = 3,
    catmull_rom = 4,
    mitchell = 5

}
~~~

Filter for resizing

- Values:
    - def:      default
    - box:      box filter
    - triangle:      triangle filter
    - cubic_spline:      cubic spline
    - catmull_rom:      Catmull-Rom interpolating sline
    - mitchell:      Mitchel-Netrevalli filter with B=1/3, C=1/3


#### Enum resize_edge

~~~ .cpp
enum struct resize_edge {
    def = 0,
    clamp = 1,
    reflect = 2,
    wrap = 3,
    zero = 4

}
~~~

Edge mode for resizing

- Values:
    - def:      default
    - clamp:      clamp
    - reflect:      reflect
    - wrap:      wrap
    - zero:      zero


#### Function resize_image()

~~~ .cpp
void resize_image(const image4f& img, image4f& res_img,
    resize_filter filter = resize_filter::def,
    resize_edge edge = resize_edge::def, bool premultiplied_alpha = false);
~~~

Resize image.

#### Function resize_image()

~~~ .cpp
void resize_image(const image4b& img, image4b& res_img,
    resize_filter filter = resize_filter::def,
    resize_edge edge = resize_edge::def, bool premultiplied_alpha = false);
~~~

Resize image.

#### Function intersect_point()

~~~ .cpp
inline bool intersect_point(
    const ray3f& ray, const vec3f& p, float r, float& ray_t);
~~~

Intersect a ray with a point (approximate)

Parameters:
- ray: ray origin and direction, parameter min, max range
- p: point position
- r: point radius

Out Parameters:
- ray_t: ray parameter at the intersection point
- euv: primitive uv ( {0,0} for points )

Returns:
- whether the intersection occurred

Iplementation Notes:
- out Parameters and only writtent o if an intersection occurs
- algorithm finds the closest point on the ray segment to the point and
   test their distance with the point radius
- based on http://geomalgorithms.com/a02-lines.html.

#### Function intersect_line()

~~~ .cpp
inline bool intersect_line(const ray3f& ray, const vec3f& v0, const vec3f& v1,
    float r0, float r1, float& ray_t, vec2f& euv);
~~~

Intersect a ray with a line

Parameters:
- ray: ray origin and direction, parameter min, max range
- v0, v1: line segment points
- r0, r1: line segment radia

Out Parameters:
- ray_t: ray parameter at the intersection point
- euv: euv.x is the line parameter at the intersection ( euv.y is zero )

Returns:
- whether the intersection occurred

Notes:
- out Parameters and only writtent o if an intersection occurs
- algorithm find the closest points on line and ray segment and test
  their distance with the line radius at that location
- based on http://geomalgorithms.com/a05-intersect-1.html
- based on http://geomalgorithms.com/a07-distance.html#
    dist3D_Segment_to_Segment

#### Function intersect_triangle()

~~~ .cpp
inline bool intersect_triangle(const ray3f& ray, const vec3f& v0,
    const vec3f& v1, const vec3f& v2, float& ray_t, vec3f& euv);
~~~

Intersect a ray with a triangle

Parameters:
- ray: ray origin and direction, parameter min, max range
- v0, v1, v2: triangle vertices

Out Parameters:
- ray_t: ray parameter at the intersection point
- euv: baricentric coordinates of the intersection

Returns:
- whether the intersection occurred

Notes:
- out Parameters and only writtent o if an intersection occurs
- algorithm based on Muller-Trombone intersection test

#### Function intersect_quad()

~~~ .cpp
inline bool intersect_quad(const ray3f& ray, const vec3f& v0, const vec3f& v1,
    const vec3f& v2, const vec3f& v3, float& ray_t, vec4f& euv);
~~~

Intersect a ray with a quad represented as two triangles (0,1,3) and
(2,3,1), with the uv coordinates of the second triangle corrected by u =
1-u' and v = 1-v' to produce a quad parametrization where u and v go from 0
to 1. This is equivalent to Intel's Embree. The external user does not have
to be concerned about the parametrization and can just use the euv as
specified.

Parameters:
- ray: ray origin and direction, parameter min, max range
- v0, v1, v2, v3: quad vertices

Out Parameters:
- ray_t: ray parameter at the intersection point
- euv: baricentric coordinates of the intersection

Returns:
- whether the intersection occurred

#### Function intersect_tetrahedron()

~~~ .cpp
inline bool intersect_tetrahedron(const ray3f& ray_, const vec3f& v0,
    const vec3f& v1, const vec3f& v2, const vec3f& v3, float& ray_t,
    vec4f& euv);
~~~

Intersect a ray with a tetrahedron. Note that we consider only
intersection wiht the tetrahedra surface and discount intersction with
the interior.

Parameters:
- ray: ray to intersect with
- v0, v1, v2: triangle vertices

Out Parameters:
- ray_t: ray parameter at the intersection point
- euv: baricentric coordinates of the intersection

Returns:
- whether the intersection occurred

TODO: check order
TODO: uv

#### Function intersect_check_bbox()

~~~ .cpp
inline bool intersect_check_bbox(const ray3f& ray, const bbox3f& bbox);
~~~

Intersect a ray with a axis-aligned bounding box

Parameters:
- ray: ray to intersect with
- bbox: bounding box min/max bounds

Returns:
- whether the intersection occurred

#### Function _safemin()

~~~ .cpp
static inline const float& _safemin(const float& a, const float& b);
~~~

Min/max used in BVH traversal. Copied here since the traversal code
relies on the specific behaviour wrt NaNs.

#### Function _safemax()

~~~ .cpp
static inline const float& _safemax(const float& a, const float& b);
~~~

Min/max used in BVH traversal. Copied here since the traversal code
relies on the specific behaviour wrt NaNs.

#### Function intersect_check_bbox()

~~~ .cpp
inline bool intersect_check_bbox(const ray3f& ray, const vec3f& ray_dinv,
    const vec3i& ray_dsign, const bbox3f& bbox_);
~~~

Intersect a ray with a axis-aligned bounding box

Parameters:
- ray_o, ray_d: ray origin and direction
- ray_tmin, ray_tmax: ray parameter min, max range
- ray_dinv: ray inverse direction
- ray_dsign: ray direction sign
- bbox_min, bbox_max: bounding box min/max bounds

Returns:
- whether the intersection occurred

Implementation Notes:
- based on "Robust BVH Ray Traversal" by T. Ize published at
http://jcgt.org/published/0002/02/02/paper.pdf

#### Struct bvh_node

~~~ .cpp
struct bvh_node {
    bbox3f bbox;
    uint32_t start;
    uint16_t count;
    uint8_t isleaf;
    uint8_t axis;
}
~~~

BVH tree node containing its bounds, indices to the BVH arrays of either
sorted primitives or internal nodes, whether its a leaf or an internal node,
and the split axis. Leaf and internal nodes are identical, except that
indices refer to primitives for leaf nodes or other nodes for internal
nodes. See bvh_tree for more details.

This is an internal data structure.

- Members:
    - bbox:      bounding box
    - start:      index to the first sorted primitive/node
    - count:      number of primitives/nodes
    - isleaf:      whether it is a leaf
    - axis:      slit axis


#### Struct bvh_tree

~~~ .cpp
struct bvh_tree {
    vector<bvh_node> nodes;
    vector<int> sorted_prim;
}
~~~

BVH tree, stored as a node array. The tree structure is encoded using array
indices instead of pointers, both for speed but also to simplify code.
BVH nodes indices refer to either the node array, for internal nodes,
or a primitive array, for leaf nodes. BVH trees may contain only one type
of geometric primitive, like points, lines, triangle or shape other BVHs.
To handle multiple primitive types and transformed primitices, build
a two-level hierarchy with the outer BVH, the scene BVH, containing inner
BVHs, shape BVHs, each of which of a uniform primitive type.

This is an internal data structure.

- Members:
    - nodes:      sorted array of internal nodes
    - sorted_prim:      sorted elements


#### Function build_bvh()

~~~ .cpp
inline bvh_tree* build_bvh(
    int nprims, bool equalsize, const function<bbox3f(int)>& elem_bbox);
~~~

Build a BVH from a set of primitives.

#### Function build_triangles_bvh()

~~~ .cpp
inline bvh_tree* build_triangles_bvh(const vector<vec3i>& triangles,
    const vector<vec3f>& pos, bool equal_size = true);
~~~

Build a triangles BVH.

#### Function build_quads_bvh()

~~~ .cpp
inline bvh_tree* build_quads_bvh(const vector<vec4i>& quads,
    const vector<vec3f>& pos, bool equal_size = true);
~~~

Build a quads BVH.

#### Function build_lines_bvh()

~~~ .cpp
inline bvh_tree* build_lines_bvh(const vector<vec2i>& lines,
    const vector<vec3f>& pos, const vector<float>& radius,
    bool equal_size = true);
~~~

Build a lines BVH.

#### Function build_points_bvh()

~~~ .cpp
inline bvh_tree* build_points_bvh(const vector<int>& points,
    const vector<vec3f>& pos, const vector<float>& radius,
    bool equal_size = true);
~~~

Build a points BVH.

#### Function build_points_bvh()

~~~ .cpp
inline bvh_tree* build_points_bvh(const vector<vec3f>& pos,
    const vector<float>& radius, bool equal_size = true);
~~~

Build a points BVH.

#### Function refit_bvh()

~~~ .cpp
inline void refit_bvh(
    bvh_tree* bvh, int nodeid, const function<bbox3f(int)>& elem_bbox);
~~~

Recursively recomputes the node bounds for a shape bvh

#### Function refit_triangles_bvh()

~~~ .cpp
inline void refit_triangles_bvh(
    bvh_tree* bvh, const vec3i* triangles, const vec3f* pos);
~~~

Refit triangles bvh

#### Function refit_triangles_bvh()

~~~ .cpp
inline void refit_triangles_bvh(
    bvh_tree* bvh, const vector<vec3i>& triangles, const vector<vec3f>& pos);
~~~

Refit triangles bvh

#### Function refit_quads_bvh()

~~~ .cpp
inline void refit_quads_bvh(
    bvh_tree* bvh, const vec4i* quads, const vec3f* pos);
~~~

Refit quads bvh

#### Function refit_quads_bvh()

~~~ .cpp
inline void refit_quads_bvh(
    bvh_tree* bvh, const vector<vec4i>& quads, const vector<vec3f>& pos);
~~~

Refit quads bvh

#### Function refit_lines_bvh()

~~~ .cpp
inline void refit_lines_bvh(
    bvh_tree* bvh, const vec2i* lines, const vec3f* pos, const float* radius);
~~~

Refit lines bvh

#### Function refit_lines_bvh()

~~~ .cpp
inline void refit_lines_bvh(bvh_tree* bvh, const vector<vec2i>& lines,
    const vector<vec3f>& pos, const vector<float>& radius);
~~~

Refit lines bvh

#### Function refit_points_bvh()

~~~ .cpp
inline void refit_points_bvh(
    bvh_tree* bvh, const int* points, const vec3f* pos, const float* radius);
~~~

Refit points bvh

#### Function refit_points_bvh()

~~~ .cpp
inline void refit_points_bvh(bvh_tree* bvh, const vector<int>& points,
    const vector<vec3f>& pos, const vector<float>& radius);
~~~

Refit points bvh

#### Function refit_points_bvh()

~~~ .cpp
inline void refit_points_bvh(
    bvh_tree* bvh, const vec3f* pos, const float* radius);
~~~

Refit points bvh

#### Function refit_points_bvh()

~~~ .cpp
inline void refit_points_bvh(
    bvh_tree* bvh, const vector<vec3f>& pos, const vector<float>& radius);
~~~

Refit lines bvh

#### Function intersect_bvh()

~~~ .cpp
inline bool intersect_bvh(const bvh_tree* bvh, const ray3f& ray_,
    bool early_exit, float& ray_t, int& eid,
    const function<bool(int, const ray3f&, float&)>& intersect_elem);
~~~

Intersect ray with a bvh.

#### Function overlap_bvh()

~~~ .cpp
inline bool overlap_bvh(const bvh_tree* bvh, const vec3f& pos, float max_dist,
    bool early_exit, float& dist, int& eid,
    const function<bool(int, const vec3f&, float, float&)>& overlap_elem);
~~~

Finds the closest element with a bvh.

#### Function intersect_triangles_bvh()

~~~ .cpp
inline bool intersect_triangles_bvh(const bvh_tree* bvh, const vec3i* triangles,
    const vec3f* pos, const ray3f& ray, bool early_exit, float& ray_t, int& eid,
    vec3f& euv);
~~~

Intersect a triangle BVH

#### Function intersect_triangles_bvh()

~~~ .cpp
inline bool intersect_triangles_bvh(const bvh_tree* bvh,
    const vector<vec3i>& triangles, const vector<vec3f>& pos, const ray3f& ray,
    bool early_exit, float& ray_t, int& eid, vec3f& euv);
~~~

Intersect a triangle BVH

#### Function intersect_quads_bvh()

~~~ .cpp
inline bool intersect_quads_bvh(const bvh_tree* bvh, const vec4i* quads,
    const vec3f* pos, const ray3f& ray, bool early_exit, float& ray_t, int& eid,
    vec4f& euv);
~~~

Intersect a quad BVH

#### Function intersect_quads_bvh()

~~~ .cpp
inline bool intersect_quads_bvh(const bvh_tree* bvh, const vector<vec4i>& quads,
    const vector<vec3f>& pos, const ray3f& ray, bool early_exit, float& ray_t,
    int& eid, vec4f& euv);
~~~

Intersect a quad BVH

#### Function intersect_lines_bvh()

~~~ .cpp
inline bool intersect_lines_bvh(const bvh_tree* bvh, const vec2i* lines,
    const vec3f* pos, const float* radius, const ray3f& ray, bool early_exit,
    float& ray_t, int& eid, vec2f& euv);
~~~

Intersect a line BVH

#### Function intersect_lines_bvh()

~~~ .cpp
inline bool intersect_lines_bvh(const bvh_tree* bvh, const vector<vec2i>& lines,
    const vector<vec3f>& pos, const vector<float>& radius, const ray3f& ray,
    bool early_exit, float& ray_t, int& eid, vec2f& euv);
~~~

Intersect a line BVH

#### Function intersect_points_bvh()

~~~ .cpp
inline bool intersect_points_bvh(const bvh_tree* bvh, const int* points,
    const vec3f* pos, const float* radius, const ray3f& ray, bool early_exit,
    float& ray_t, int& eid);
~~~

Intersect a point BVH

#### Function intersect_points_bvh()

~~~ .cpp
inline bool intersect_points_bvh(const bvh_tree* bvh, const vector<int>& points,
    const vector<vec3f>& pos, const vector<float>& radius, const ray3f& ray,
    bool early_exit, float& ray_t, int& eid);
~~~

Intersect a point BVH

#### Function intersect_points_bvh()

~~~ .cpp
inline bool intersect_points_bvh(const bvh_tree* bvh, const vec3f* pos,
    const float* radius, const ray3f& ray, bool early_exit, float& ray_t,
    int& eid);
~~~

Intersect a point BVH

#### Function intersect_points_bvh()

~~~ .cpp
inline bool intersect_points_bvh(const bvh_tree* bvh, const vector<vec3f>& pos,
    const vector<float>& radius, const ray3f& ray, bool early_exit,
    float& ray_t, int& eid);
~~~

Intersect a point BVH

#### Function overlap_triangles_bvh()

~~~ .cpp
inline bool overlap_triangles_bvh(const bvh_tree* bvh, const vec3i* triangles,
    const vec3f* pos, const float* radius, const vec3f& pt, float max_dist,
    bool early_exit, float& dist, int& eid, vec3f& euv);
~~~

Intersect a triangle BVH

#### Function overlap_triangles_bvh()

~~~ .cpp
inline bool overlap_triangles_bvh(const bvh_tree* bvh,
    const vector<vec3i>& triangles, const vector<vec3f>& pos,
    const vector<float>& radius, const vec3f& pt, float max_dist,
    bool early_exit, float& dist, int& eid, vec3f& euv);
~~~

Intersect a triangle BVH

#### Function overlap_quads_bvh()

~~~ .cpp
inline bool overlap_quads_bvh(const bvh_tree* bvh, const vec4i* quads,
    const vec3f* pos, const float* radius, const vec3f& pt, float max_dist,
    bool early_exit, float& dist, int& eid, vec4f& euv);
~~~

Intersect a quad BVH

#### Function overlap_quads_bvh()

~~~ .cpp
inline bool overlap_quads_bvh(const bvh_tree* bvh, const vector<vec4i>& quads,
    const vector<vec3f>& pos, const vector<float>& radius, const vec3f& pt,
    float max_dist, bool early_exit, float& dist, int& eid, vec4f& euv);
~~~

Intersect a quad BVH

#### Function overlap_lines_bvh()

~~~ .cpp
inline bool overlap_lines_bvh(const bvh_tree* bvh, const vec2i* lines,
    const vec3f* pos, const float* radius, const vec3f& pt, float max_dist,
    bool early_exit, float& dist, int& eid, vec2f& euv);
~~~

Intersect a line BVH

#### Function overlap_lines_bvh()

~~~ .cpp
inline bool overlap_lines_bvh(const bvh_tree* bvh, const vector<vec2i>& lines,
    const vector<vec3f>& pos, const vector<float>& radius, const vec3f& pt,
    float max_dist, bool early_exit, float& dist, int& eid, vec2f& euv);
~~~

Intersect a line BVH

#### Function overlap_points_bvh()

~~~ .cpp
inline bool overlap_points_bvh(const bvh_tree* bvh, const int* points,
    const vec3f* pos, const float* radius, const vec3f& pt, float max_dist,
    bool early_exit, float& dist, int& eid);
~~~

Intersect a point BVH

#### Function overlap_points_bvh()

~~~ .cpp
inline bool overlap_points_bvh(const bvh_tree* bvh, const vector<int>& points,
    const vector<vec3f>& pos, const vector<float>& radius, const vec3f& pt,
    float max_dist, bool early_exit, float& dist, int& eid);
~~~

Intersect a point BVH

#### Function overlap_points_bvh()

~~~ .cpp
inline bool overlap_points_bvh(const bvh_tree* bvh, const vec3f* pos,
    const float* radius, const vec3f& pt, float max_dist, bool early_exit,
    float& dist, int& eid);
~~~

Intersect a point BVH

#### Function overlap_points_bvh()

~~~ .cpp
inline bool overlap_points_bvh(const bvh_tree* bvh, const vector<vec3f>& pos,
    const vector<float>& radius, const vec3f& pt, float max_dist,
    bool early_exit, float& dist, int& eid);
~~~

Intersect a point BVH

#### Function overlap_bvh_elems()

~~~ .cpp
template <typename OverlapElem>
void overlap_bvh_elems(const bvh_tree* bvh1, const bvh_tree* bvh2,
    bool skip_duplicates, bool skip_self, vector<vec2i>& overlaps,
    const OverlapElem& overlap_elems);
~~~

Finds the overlap between BVH leaf nodes.

#### Struct texture

~~~ .cpp
struct texture {
    string name;
    string path;
    image4b ldr;
    image4f hdr;
    int width() const; 
    int height() const; 
}
~~~

Scene Texture

- Members:
    - name:      name
    - path:      path
    - ldr:      if loaded, ldr image
    - hdr:      if loaded, hdr image
    - width():      get texture width
    - height():      get texture height


#### Struct texture_info

~~~ .cpp
struct texture_info {
    texture* txt = nullptr;
    bool wrap_s = true;
    bool wrap_t = true;
    bool linear = true;
    bool mipmap = true;
    float scale = 1;
    operator bool() const; 
}
~~~

Scene Texture Additional Information

- Members:
    - txt:      texture pointer
    - wrap_s:      wrap s coordinate
    - wrap_t:      wrap t coordinate
    - linear:      linear interpolation
    - mipmap:      mipmaping
    - scale:      texture strength (occlusion and normal)
    - operator bool():      check whether the texture if present


#### Enum material_type

~~~ .cpp
enum struct material_type {
    specular_roughness = 0,
    metallic_roughness = 1,
    specular_glossiness = 2,
}
~~~

Material type

- Values:
    - specular_roughness:      Microfacet material type (OBJ)
    - metallic_roughness:      Base and metallic material (metallic-roughness in glTF)
    - specular_glossiness:      Diffuse and specular material (specular-glossness in glTF)


#### Struct material

~~~ .cpp
struct material {
    string name;
    bool double_sided = false;
    material_type mtype = material_type::specular_roughness;
    vec3f ke = {0, 0, 0};
    vec3f kd = {0, 0, 0};
    vec3f ks = {0, 0, 0};
    vec3f kr = {0, 0, 0};
    vec3f kt = {0, 0, 0};
    float rs = 0.0001;
    float op = 1;
    texture_info ke_txt = {};
    texture_info kd_txt = {};
    texture_info ks_txt = {};
    texture_info kr_txt = {};
    texture_info kt_txt = {};
    texture_info rs_txt = {};
    texture_info bump_txt = {};
    texture_info disp_txt = {};
    texture_info norm_txt = {};
    texture_info occ_txt = {};
}
~~~

Scene Material

- Members:
    - name:      material name
    - double_sided:      double-sided rendering
    - mtype:      material type
    - ke:      emission color
    - kd:      diffuse color / base color
    - ks:      specular color / metallic factor
    - kr:      clear coat reflection
    - kt:      transmission color
    - rs:      roughness
    - op:      opacity
    - ke_txt:      emission texture
    - kd_txt:      diffuse texture
    - ks_txt:      specular texture
    - kr_txt:      reflection texture
    - kt_txt:      transmission texture
    - rs_txt:      roughness texture
    - bump_txt:      bump map texture (heighfield)
    - disp_txt:      displacement map texture (heighfield)
    - norm_txt:      normal texture
    - occ_txt:      occlusion texture


#### Struct shape

~~~ .cpp
struct shape {
    string name = "";
    string path = "";
    material* mat = nullptr;
    vector<int> points;
    vector<vec2i> lines;
    vector<vec3i> triangles;
    vector<vec4i> quads;
    vector<vec4i> quads_pos;
    vector<vec4i> quads_norm;
    vector<vec4i> quads_texcoord;
    vector<vec3f> pos;
    vector<vec3f> norm;
    vector<vec2f> texcoord;
    vector<vec2f> texcoord1;
    vector<vec4f> color;
    vector<float> radius;
    vector<vec4f> tangsp;
    vector<float> elem_cdf;
    bvh_tree* bvh = nullptr;
    bbox3f bbox = invalid_bbox3f;
}
~~~

Shape data represented as an indexed array.
May contain only one of the points/lines/triangles/quads.

- Members:
    - name:      shape name
    - path:      path (used for saving in glTF)
    - mat:      shape material
    - points:      points
    - lines:      lines
    - triangles:      triangles
    - quads:      quads
    - quads_pos:      face-varying indices for position
    - quads_norm:      face-varying indices for normal
    - quads_texcoord:      face-varying indices for texcoord
    - pos:      per-vertex position (3 float)
    - norm:      per-vertex normals (3 float)
    - texcoord:      per-vertex texcoord (2 float)
    - texcoord1:      per-vertex second texcoord (2 float)
    - color:      per-vertex color (4 float)
    - radius:      per-vertex radius (1 float)
    - tangsp:      per-vertex tangent space (4 float)
    - elem_cdf:      element CDF for sampling
    - bvh:      BVH
    - bbox:      bounding box (needs to be updated explicitly)


#### Struct instance

~~~ .cpp
struct instance {
    frame3f frame = identity_frame3f;
    shape* shp = nullptr;
    bbox3f bbox = invalid_bbox3f;
    mat4f xform() const; 
}
~~~

Shape instance.

- Members:
    - frame:      transform frame
    - shp:      shape instance
    - bbox:      bounding box (needs to be updated explicitly)
    - xform():      instance transform as matrix


#### Struct camera

~~~ .cpp
struct camera {
    string name;
    frame3f frame = identity_frame3f;
    bool ortho = false;
    float yfov = 2;
    float aspect = 16.0f / 9.0f;
    float focus = 1;
    float aperture = 0;
    float near = 0.01f;
    float far = 10000;
}
~~~

Scene Camera

- Members:
    - name:      name
    - frame:      transform frame
    - ortho:      ortho cam
    - yfov:      vertical field of view
    - aspect:      aspect ratio
    - focus:      focus distance
    - aperture:      lens aperture
    - near:      near plane distance
    - far:      far plane distance


#### Struct environment

~~~ .cpp
struct environment {
    string name;
    frame3f frame = identity_frame3f;
    vec3f ke = {0, 0, 0};
    texture_info ke_txt = {};
}
~~~

Envinonment map

- Members:
    - name:      name
    - frame:      transform frame
    - ke:      emission coefficient
    - ke_txt:      emission texture


#### Struct light

~~~ .cpp
struct light {
    instance* ist = nullptr;
    environment* env = nullptr;
}
~~~

Light, either an instance or an environment.
This is only used internally to avoid looping over all objects every time.

- Members:
    - ist:      instance
    - env:      environment


#### Struct scene

~~~ .cpp
struct scene {
    vector<shape*> shapes;
    vector<instance*> instances;
    vector<material*> materials;
    vector<texture*> textures;
    vector<camera*> cameras;
    vector<environment*> environments;
    vector<light*> lights;
    bvh_tree* bvh = nullptr;
    bbox3f bbox = invalid_bbox3f;
    ~scene(); 
}
~~~

Scene

- Members:
    - shapes:      shape array
    - instances:      instance array
    - materials:      material array
    - textures:      texture array
    - cameras:      camera array
    - environments:      environment array
    - lights:      light array
    - bvh:      BVH
    - bbox:      bounding box (needs to be updated explicitly)
    - ~scene():      cleanup


#### Function eval_barycentric()

~~~ .cpp
template <typename T>
inline T eval_barycentric(
    const shape* shp, const vector<T>& vals, int eid, const vec4f& euv);
~~~

Shape value interpolated using barycentric coordinates

#### Function eval_pos()

~~~ .cpp
inline vec3f eval_pos(const shape* shp, int eid, const vec4f& euv);
~~~

Shape position interpolated using barycentric coordinates

#### Function eval_norm()

~~~ .cpp
inline vec3f eval_norm(const shape* shp, int eid, const vec4f& euv);
~~~

Shape normal interpolated using barycentric coordinates

#### Function eval_texcoord()

~~~ .cpp
inline vec2f eval_texcoord(const shape* shp, int eid, const vec4f& euv);
~~~

Shape texcoord interpolated using barycentric coordinates

#### Function eval_color()

~~~ .cpp
inline vec4f eval_color(const shape* shp, int eid, const vec4f& euv);
~~~

Shape texcoord interpolated using barycentric coordinates

#### Function eval_tangsp()

~~~ .cpp
inline vec4f eval_tangsp(const shape* shp, int eid, const vec4f& euv);
~~~

Shape tangent space interpolated using barycentric coordinates

#### Function eval_pos()

~~~ .cpp
inline vec3f eval_pos(const instance* ist, int eid, const vec4f& euv);
~~~

Instance position interpolated using barycentric coordinates

#### Function eval_norm()

~~~ .cpp
inline vec3f eval_norm(const instance* ist, int eid, const vec4f& euv);
~~~

Instance normal interpolated using barycentric coordinates

#### Function eval_texture()

~~~ .cpp
inline vec4f eval_texture(const texture_info& info, const vec2f& texcoord,
    bool srgb = true, const vec4f& def =;
~~~

Evaluate a texture

#### Function subdivide_shape()

~~~ .cpp
inline void subdivide_shape(shape* shp, bool subdiv = false);
~~~

Subdivides shape elements. Apply subdivision surface rules if subdivide
is true.

#### Function facet_shape()

~~~ .cpp
inline void facet_shape(shape* shp);
~~~

Facet a shape. Supports only non-face0varying shapes

#### Function tesselate_shape()

~~~ .cpp
inline void tesselate_shape(shape* shp);
~~~

Tesselate a shape into basic primitives

#### Function tesselate_shapes()

~~~ .cpp
inline void tesselate_shapes(scene* scn);
~~~

Tesselate scene shapes and update pointers

#### Struct load_options

~~~ .cpp
struct load_options {
    bool load_textures = true;
    bool skip_missing = true;
    bool obj_flip_texcoord = true;
    bool obj_facet_non_smooth = false;
    bool obj_flip_tr = true;
    bool preserve_quads = false;
    bool preserve_facevarying = false;
}
~~~

Loading options

- Members:
    - load_textures:      Whether to load textures
    - skip_missing:      Skip missing files without giving and error
    - obj_flip_texcoord:      Whether to flip the v coordinate in OBJ
    - obj_facet_non_smooth:      Duplicate vertices if smoothing off in OBJ
    - obj_flip_tr:      Whether to flip tr in OBJ
    - preserve_quads:      whether to preserve quads
    - preserve_facevarying:      whether to preserve face-varying faces


#### Function load_scene()

~~~ .cpp
scene* load_scene(const string& filename, const load_options& opts =;
~~~

Loads a scene. For now OBJ or glTF are supported.
Throws an exception if an error occurs.

#### Struct save_options

~~~ .cpp
struct save_options {
    bool save_textures = true;
    bool skip_missing = true;
    bool obj_flip_texcoord = true;
    bool obj_flip_tr = true;
    bool gltf_separate_buffers = false;
}
~~~

Save options

- Members:
    - save_textures:      Whether to save textures
    - skip_missing:      Skip missing files without giving and error
    - obj_flip_texcoord:      Whether to flip the v coordinate in OBJ
    - obj_flip_tr:      Whether to flip tr in OBJ
    - gltf_separate_buffers:      Whether to use separate buffers in gltf


#### Function save_scene()

~~~ .cpp
void save_scene(
    const string& filename, const scene* scn, const save_options& opts);
~~~

Saves a scene. For now OBJ and glTF are supported.
Throws an exception if an error occurs.

#### Struct add_elements_options

~~~ .cpp
struct add_elements_options {
    bool smooth_normals = true;
    float pointline_radius = 0;
    bool tangent_space = true;
    bool texture_data = true;
    bool shape_instances = true;
    bool default_camera = true;
    bool default_environment = false;
    bool default_names = true;
    bool default_paths = true;
    static add_elements_options none(); 
}
~~~

Add elements options

- Members:
    - smooth_normals:      Add missing normal
    - pointline_radius:      Add missing radius for points and lines (<=0 for no adding)
    - tangent_space:      Add missing trangent space
    - texture_data:      texture data
    - shape_instances:      Add instances
    - default_camera:      Add default camera
    - default_environment:      Add an empty default environment
    - default_names:      Add default names
    - default_paths:      Add default paths
    - none():      initialize to no element


#### Function add_elements()

~~~ .cpp
void add_elements(scene* scn, const add_elements_options& opts =;
~~~

Add elements

#### Function merge_into()

~~~ .cpp
void merge_into(scene* merge_into, scene* merge_from);
~~~

Merge scene into one another. Note that the objects are _moved_ from
merge_from to merged_into, so merge_from will be empty after this function.

#### Function update_bounds()

~~~ .cpp
inline void update_bounds(shape* shp);
~~~

Computes a shape bounding box (quick computation that ignores radius)

#### Function update_bounds()

~~~ .cpp
inline void update_bounds(instance* ist, bool do_shape = true);
~~~

Updates the instance bounding box

#### Function update_bounds()

~~~ .cpp
inline void update_bounds(scene* scn, bool do_shapes = true);
~~~

Updates the scene and scene's instances bounding boxes

#### Function flatten_instances()

~~~ .cpp
inline void flatten_instances(scene* scn);
~~~

Flatten scene instances into separate meshes.

#### Function update_lights()

~~~ .cpp
void update_lights(
    scene* scn, bool include_env = false, bool sampling_cdf = false);
~~~

Initialize the lights

#### Function print_info()

~~~ .cpp
void print_info(const scene* scn);
~~~

Print scene information (call update bounds bes before)

#### Function build_bvh()

~~~ .cpp
inline void build_bvh(shape* shp, bool equalsize = true);
~~~

Build a shape BVH

#### Function build_bvh()

~~~ .cpp
inline void build_bvh(
    scene* scn, bool equalsize = true, bool do_shapes = true);
~~~

Build a scene BVH

#### Function refit_bvh()

~~~ .cpp
inline void refit_bvh(shape* shp);
~~~

Refits a scene BVH

#### Function refit_bvh()

~~~ .cpp
inline void refit_bvh(scene* scn, bool do_shapes = true);
~~~

Refits a scene BVH

#### Function intersect_ray()

~~~ .cpp
inline bool intersect_ray(const shape* shp, const ray3f& ray, bool early_exit,
    float& ray_t, int& eid, vec4f& euv);
~~~

Intersect the shape with a ray. Find any interstion if early_exit,
otherwise find first intersection.

- Parameters:
    - scn: scene to intersect
    - ray: ray to be intersected
    - early_exit: whether to stop at the first found hit
    - ray_t: ray distance at intersection
    - eid: shape element index
    - euv: element barycentric coordinates
- Returns:
    - whether it intersected

#### Function intersect_ray()

~~~ .cpp
inline bool intersect_ray(const instance* ist, const ray3f& ray,
    bool early_exit, float& ray_t, int& eid, vec4f& euv);
~~~

Intersect the instance with a ray. Find any interstion if early_exit,
otherwise find first intersection.

- Parameters:
    - scn: scene to intersect
    - ray: ray to be intersected
    - early_exit: whether to stop at the first found hit
    - ray_t: ray distance at intersection
    - eid: shape element index
    - euv: element barycentric coordinates
- Returns:
    - whether it intersected

#### Function intersect_ray()

~~~ .cpp
inline bool intersect_ray(const scene* scn, const ray3f& ray, bool early_exit,
    float& ray_t, int& iid, int& eid, vec4f& euv);
~~~

Intersect the scene with a ray. Find any interstion if early_exit,
otherwise find first intersection.

- Parameters:
    - scn: scene to intersect
    - ray: ray to be intersected
    - early_exit: whether to stop at the first found hit
    - ray_t: ray distance at intersection
    - iid: instance index
    - eid: shape element index
    - euv: element barycentric coordinates
- Returns:
    - whether it intersected

#### Struct intersection_point

~~~ .cpp
struct intersection_point {
    float dist = 0;
    int iid = -1;
    int eid = -1;
    vec4f euv = zero4f;
    operator bool() const; 
}
~~~

Surface point.

- Members:
    - dist:      distance of the hit along the ray or from the point
    - iid:      instance index
    - eid:      shape element index
    - euv:      shape barycentric coordinates
    - operator bool():      check if intersection is valid


#### Function intersect_ray()

~~~ .cpp
inline intersection_point intersect_ray(
    const scene* scn, const ray3f& ray, bool early_exit);
~~~

Intersect the scene with a ray. Find any interstion if early_exit,
otherwise find first intersection.

- Parameters:
    - scn: scene to intersect
    - ray: ray to be intersected
    - early_exit: whether to stop at the first found hit
- Returns:
    - intersection record

#### Function overlap_point()

~~~ .cpp
inline bool overlap_point(const shape* shp, const vec3f& pos, float max_dist,
    bool early_exit, float& dist, int& eid, vec4f& euv);
~~~

Finds the closest element that overlaps a point within a given distance.

- Parameters:
    - scn: scene to intersect
    - pos: point position
    - max_dist: maximu valid distance
    - early_exit: whether to stop at the first found hit
    - dist: distance at intersection
    - eid: shape element index
    - euv: element barycentric coordinates
- Returns:
    - whether it intersected

#### Function overlap_point()

~~~ .cpp
inline bool overlap_point(const instance* ist, const vec3f& pos, float max_dist,
    bool early_exit, float& dist, int& eid, vec4f& euv);
~~~

Finds the closest element that overlaps a point within a given distance.

- Parameters:
    - scn: scene to intersect
    - pos: point position
    - max_dist: maximu valid distance
    - early_exit: whether to stop at the first found hit
    - dist: distance at intersection
    - eid: shape element index
    - euv: element barycentric coordinates
- Returns:
    - whether it intersected

#### Function overlap_point()

~~~ .cpp
inline bool overlap_point(const scene* scn, const vec3f& pos, float max_dist,
    bool early_exit, float& dist, int& iid, int& eid, vec4f& euv);
~~~

Finds the closest element that overlaps a point within a given distance.

- Parameters:
    - scn: scene to intersect
    - pos: point position
    - max_dist: maximu valid distance
    - early_exit: whether to stop at the first found hit
    - dist: distance at intersection
    - iid: instance index
    - eid: shape element index
    - euv: element barycentric coordinates
- Returns:
    - whether it intersected

#### Function overlap_instance_bounds()

~~~ .cpp
inline void overlap_instance_bounds(const scene* scn1, const scene* scn2,
    bool skip_duplicates, bool skip_self, vector<vec2i>& overlaps);
~~~

Find the list of overlaps between instance bounds.

#### Function make_cornell_box_scene()

~~~ .cpp
scene* make_cornell_box_scene();
~~~

Makes the Cornell Box scene

#### Enum test_scene_type

~~~ .cpp
enum struct test_scene_type {
~~~

Test scene enumeration.

#### Function test_scene_names()

~~~ .cpp
inline const vector<pair<string, test_scene_type>>& test_scene_names();
~~~

Names for enumeration

#### Function make_test_scene()

~~~ .cpp
scene* make_test_scene(test_scene_type stype);
~~~

Makes a test scene

#### Enum trace_shader_type

~~~ .cpp
enum struct trace_shader_type {
    pathtrace = 0,
    eyelight,
    direct,
    pathtrace_nomis,
    debug_normal,
    debug_albedo,
    debug_texcoord,
}
~~~

Type of rendering algorithm (shader)

- Values:
    - pathtrace:      pathtrace
    - eyelight:      eye hight for quick previews
    - direct:      direct illumination
    - pathtrace_nomis:      pathtrace without MIS (usedful ony for debugging)
    - debug_normal:      debug normal
    - debug_albedo:      debug albedo
    - debug_texcoord:      debug texcoord


#### Function trace_shader_names()

~~~ .cpp
inline const vector<pair<string, trace_shader_type>>& trace_shader_names();
~~~

Names for enumeration

#### Enum trace_rng_type

~~~ .cpp
enum struct trace_rng_type {
    uniform = 0,
    stratified,
}
~~~

Random number generator type

- Values:
    - uniform:      uniform random numbers
    - stratified:      stratified random numbers


#### Function trace_rng_names()

~~~ .cpp
inline const vector<pair<string, trace_rng_type>>& trace_rng_names();
~~~

Names for enumeration

#### Enum trace_filter_type

~~~ .cpp
enum struct trace_filter_type {
    box = 1,
    triangle = 2,
    cubic = 3,
    catmull_rom = 4,
    mitchell = 5

}
~~~

Filter type

- Values:
    - box:      box filter
    - triangle:      hat filter
    - cubic:      cubic spline
    - catmull_rom:      Catmull-Rom spline
    - mitchell:      Mitchell-Netrevalli


#### Function trace_filter_names()

~~~ .cpp
inline const vector<pair<string, trace_filter_type>>& trace_filter_names();
~~~

Names for enumeration

#### Struct trace_params

~~~ .cpp
struct trace_params {
    int camera_id = 0;
    int width = 360;
    int height = 360;
    int nsamples = 256;
    trace_shader_type stype = trace_shader_type::pathtrace;
    bool shadow_notransmission = false;
    trace_rng_type rtype = trace_rng_type::stratified;
    trace_filter_type ftype = trace_filter_type::box;
    vec3f ambient = {0, 0, 0};
    bool envmap_invisible = false;
    int min_depth = 3;
    int max_depth = 8;
    float pixel_clamp = 10;
    float ray_eps = 1e-4f;
    bool parallel = true;
    uint32_t seed = 0;
    int block_size = 32;
}
~~~

Rendering params

- Members:
    - camera_id:      camera id
    - width:      image width
    - height:      image height
    - nsamples:      number of samples
    - stype:      sampler type
    - shadow_notransmission:      wheter to test transmission in shadows
    - rtype:      random number generation type
    - ftype:      filter type
    - ambient:      ambient lighting
    - envmap_invisible:      view environment map
    - min_depth:      minimum ray depth
    - max_depth:      maximum ray depth
    - pixel_clamp:      final pixel clamping
    - ray_eps:      ray intersection epsilon
    - parallel:      parallel execution
    - seed:      seed for the random number generators
    - block_size:      block size for parallel batches (probably leave it as is)


#### Function trace_blocks()

~~~ .cpp
inline vector<pair<vec2i, vec2i>> trace_blocks(const trace_params& params);
~~~

Make image blocks

#### Function trace_rngs()

~~~ .cpp
inline vector<rng_pcg32> trace_rngs(const trace_params& params);
~~~

Make a 2D array of random number generators for parallelization

#### Function trace_block()

~~~ .cpp
void trace_block(const scene* scn, image4f& img, const vec2i& block_min,
    const vec2i& block_max, int samples_min, int samples_max,
    vector<rng_pcg32>& rngs, const trace_params& params);
~~~

Renders a block of samples

Notes: It is safe to call the function in parallel on different blocks.
But two threads should not access the same pixels at the same time. If
the same block is rendered with different samples, samples have to be
sequential.

- Parameters:
    - scn: trace scene
    - img: pixel data in RGBA format (width/height in params)
    - block: range of pixels to render
    - samples_min, samples_max: range of samples to render
    - params: trace params

#### Function trace_samples()

~~~ .cpp
void trace_samples(const scene* scn, image4f& img, int samples_min,
    int samples_max, vector<rng_pcg32>& rngs, const trace_params& params);
~~~

Trace the next samples in [samples_min, samples_max) range.
Samples have to be traced consecutively.

#### Function trace_block_filtered()

~~~ .cpp
void trace_block_filtered(const scene* scn, image4f& img, image4f& acc,
    image4f& weight, const vec2i& block_min, const vec2i& block_max,
    int samples_min, int samples_max, vector<rng_pcg32>& rngs,
    std::mutex& image_mutex, const trace_params& params);
~~~

Renders a filtered block of samples

Notes: It is safe to call the function in parallel on different blocks.
But two threads should not access the same pixels at the same time. If
the same block is rendered with different samples, samples have to be
sequential.

- Parameters:
    - scn: trace scene
    - img: pixel data in RGBA format (width/height in params)
    - acc: accumulation buffer in RGBA format (width/height in params)
    - weight: weight buffer in float format (width/height in params)
    - block: range of pixels to render
    - samples_min, samples_max: range of samples to render
    - image_mutex: mutex for locking
    - params: trace params

#### Function trace_filtered_samples()

~~~ .cpp
void trace_filtered_samples(const scene* scn, image4f& img, image4f& acc,
    image4f& weight, int samples_min, int samples_max, vector<rng_pcg32>& rngs,
    const trace_params& params);
~~~

Trace the next samples in [samples_min, samples_max) range.
Samples have to be traced consecutively.

#### Function trace_image()

~~~ .cpp
inline image4f trace_image(const scene* scn, const trace_params& params);
~~~

Trace the whole image

#### Function trace_async_start()

~~~ .cpp
void trace_async_start(const scene* scn, image4f& img, vector<rng_pcg32>& rngs,
    const trace_params& params, thread_pool* pool,
    const function<void(int)>& callback);
~~~

Starts an anyncrhounous renderer with a maximum of 256 samples.

#### Function trace_async_stop()

~~~ .cpp
void trace_async_stop(thread_pool* pool);
~~~

Stop the asynchronous renderer.

#### Struct obj_vertex

~~~ .cpp
struct obj_vertex {
    int pos;
    int texcoord;
    int norm;
    int color;
    int radius;
    obj_vertex(int pos = -1, int texcoord = -1, int norm = -1, int color = -1, int radius = -1); 
}
~~~

Face vertex

- Members:
    - pos:      position
    - texcoord:      texcoord
    - norm:      normal
    - color:      color [extension]
    - radius:      radius [extension]
    - obj_vertex():      Constructor (copies members initializing missing ones to -1)


#### Enum obj_element_type : uint16_t

~~~ .cpp
enum struct obj_element_type : uint16_t {
    point = 1,
    line = 2,
    face = 3,
    tetra = 4,
}
~~~

element type

- Values:
    - point:      lists of points
    - line:      polylines
    - face:      polygon faces
    - tetra:      tetrahedrons


#### Struct obj_element

~~~ .cpp
struct obj_element {
    uint32_t start;
    obj_element_type type;
    uint16_t size;
}
~~~

Element vertex indices

- Members:
    - start:      starting vertex index
    - type:      element type
    - size:      number of vertices


#### Struct obj_group

~~~ .cpp
struct obj_group {
    string matname;
    string groupname;
    bool smoothing = true;
    vector<obj_vertex> verts;
    vector<obj_element> elems;
}
~~~

Element group

- Members:
    - matname:      material name
    - groupname:      group name
    - smoothing:      smoothing
    - verts:      element vertices
    - elems:      element faces


#### Struct obj_object

~~~ .cpp
struct obj_object {
    string name;
    vector<obj_group> groups;
}
~~~

Obj object

- Members:
    - name:      object name
    - groups:      element groups


#### Struct obj_texture_info

~~~ .cpp
struct obj_texture_info {
    string path = "";
    bool clamp = false;
    float scale = 1;
    unordered_map<string, vector<string>> unknown_props;
}
~~~

Texture information for OBJ

- Members:
    - path:      the texture path
    - clamp:      whether to clamp tp th edge
    - scale:      the scale for bump and displacement
    - unknown_props:      the rest of the unknown properties


#### Struct obj_texture

~~~ .cpp
struct obj_texture {
    string path;
    int width = 0;
    int height = 0;
    int ncomp = 0;
    vector<uint8_t> datab;
    vector<float> dataf;
}
~~~

OBJ texture. Texture data is loaded only if desired.

- Members:
    - path:      texture path
    - width:      Width
    - height:      Height
    - ncomp:      Number of Channels
    - datab:      Buffer data for 8-bit images
    - dataf:      Buffer data for float images


#### Struct obj_material

~~~ .cpp
struct obj_material {
    string name;
    int illum = 0;
    vec3f ke = {0, 0, 0};
    vec3f ka = {0, 0, 0};
    vec3f kd = {0, 0, 0};
    vec3f ks = {0, 0, 0};
    vec3f kr = {0, 0, 0};
    vec3f kt = {0, 0, 0};
    float ns = 1;
    float ior = 1;
    float op = 1;
    obj_texture_info ke_txt;
    obj_texture_info ka_txt;
    obj_texture_info kd_txt;
    obj_texture_info ks_txt;
    obj_texture_info kr_txt;
    obj_texture_info kt_txt;
    obj_texture_info ns_txt;
    obj_texture_info op_txt;
    obj_texture_info ior_txt;
    obj_texture_info bump_txt;
    obj_texture_info disp_txt;
    obj_texture_info norm_txt;
    unordered_map<string, vector<string>> unknown_props;
}
~~~

OBJ material

- Members:
    - name:      material name
    - illum:      MTL illum mode
    - ke:      emission color
    - ka:      ambient color
    - kd:      diffuse color
    - ks:      specular color
    - kr:      reflection color
    - kt:      transmision color
    - ns:      phong exponent for ks
    - ior:      index of refraction
    - op:      opacity
    - ke_txt:      emission texture
    - ka_txt:      ambient texture
    - kd_txt:      diffuse texture
    - ks_txt:      specular texture
    - kr_txt:      reflection texture
    - kt_txt:      transmission texture
    - ns_txt:      specular exponent texture
    - op_txt:      opacity texture
    - ior_txt:      index of refraction
    - bump_txt:      bump map texture (heighfield)
    - disp_txt:      displacement map texture (heighfield)
    - norm_txt:      normal map texture
    - unknown_props:      unknown string props


#### Struct obj_camera

~~~ .cpp
struct obj_camera {
    string name;
    frame3f frame = identity_frame3f;
    bool ortho = false;
    float yfov = 2 * atan(0.5f);
    float aspect = 16.0f / 9.0f;
    float aperture = 0;
    float focus = 1;
}
~~~

Camera [extension]

- Members:
    - name:      camera name
    - frame:      transform frame (affine matrix)
    - ortho:      orthografic camera
    - yfov:      vertical field of view
    - aspect:      aspect ratio
    - aperture:      lens aperture
    - focus:      focus distance


#### Struct obj_environment

~~~ .cpp
struct obj_environment {
    string name;
    frame3f frame = identity_frame3f;
    string matname;
}
~~~

Environment [extension]

- Members:
    - name:      environment name
    - frame:      transform frame (affine matrix)
    - matname:      material name


#### Struct obj_instance

~~~ .cpp
struct obj_instance {
    string name;
    frame3f frame = identity_frame3f;
    string objname;
}
~~~

Instance [extension]

- Members:
    - name:      instance name
    - frame:      transform frame (affine matrix)
    - objname:      object name


#### Struct obj_scene

~~~ .cpp
struct obj_scene {
    vector<vec3f> pos;
    vector<vec3f> norm;
    vector<vec2f> texcoord;
    vector<vec4f> color;
    vector<float> radius;
    vector<obj_object*> objects;
    vector<obj_material*> materials;
    vector<obj_texture*> textures;
    vector<obj_camera*> cameras;
    vector<obj_environment*> environments;
    vector<obj_instance*> instances;
    ~obj_scene(); 
}
~~~

OBJ asset

- Members:
    - pos:      vertex positions
    - norm:      vertex normals
    - texcoord:      vertex texcoord
    - color:      vertex color [extension]
    - radius:      vertex radius [extension]
    - objects:      objects
    - materials:      materials
    - textures:      textures
    - cameras:      cameras [extension]
    - environments:      env maps [extension]
    - instances:      instances [extension]
    - ~obj_scene():      cleanup


#### Function load_obj()

~~~ .cpp
obj_scene* load_obj(const string& filename, bool load_textures = false,
    bool skip_missing = false, bool flip_texcoord = true, bool flip_tr = true);
~~~

Load OBJ

- Parameters:
    - filename: filename
    - load_texture: whether to load textures
    - skip_missing: whether to skip missing files
    - flip_texcoord: whether to flip the v coordinate
    - flip_tr: whether to flip the Tr value
- Return:
    - obj (nullptr on error)

#### Function save_obj()

~~~ .cpp
void save_obj(const string& filename, const obj_scene* model,
    bool save_textures = false, bool skip_missing = false,
    bool flip_texcoord = true, bool flip_tr = true);
~~~

Save OBJ

- Parameters:
    - filename: filename
    - model: obj data to save
    - save_textures: whether to save textures
    - skip_missing: whether to skip missing files
    - flip_texcoord: whether to flip the v coordinate
    - flip_tr: whether to flip the Tr value
- Returns:
    - whether an error occurred

#### Struct obj_shape

~~~ .cpp
struct obj_shape {
    string name = "";
    string matname = "";
    vector<int> points;
    vector<vec2i> lines;
    vector<vec3i> triangles;
    vector<vec4i> tetras;
    vector<vec3f> pos;
    vector<vec3f> norm;
    vector<vec2f> texcoord;
    vector<vec4f> color;
    vector<float> radius;
}
~~~

Shape. May contain only one of the points/lines/triangles.

- Members:
    - name:      name of the group that enclosed it
    - matname:      name of the material
    - points:      points
    - lines:      lines
    - triangles:      triangles
    - tetras:      tetrahedrons
    - pos:      per-vertex position (3 float)
    - norm:      per-vertex normals (3 float)
    - texcoord:      per-vertex texcoord (2 float)
    - color:      [extension] per-vertex color (4 float)
    - radius:      [extension] per-vertex radius (1 float)


#### Struct obj_mesh

~~~ .cpp
struct obj_mesh {
    vector<obj_shape> shapes;
    ~obj_mesh(); 
}
~~~

Mesh

- Members:
    - shapes:      primitives
    - ~obj_mesh():      cleanup


#### Function get_mesh()

~~~ .cpp
obj_mesh* get_mesh(
    const obj_scene* model, const obj_object& oobj, bool facet_non_smooth);
~~~

Gets a mesh from an OBJ object.

#### Typedef buffer_data

~~~ .cpp
using buffer_data = vector<unsigned char>;
~~~

Generic buffer data.

#### Struct image_data

~~~ .cpp
struct image_data {
    int width = 0;
    int height = 0;
    int ncomp = 0;
    vector<uint8_t> datab;
    vector<float> dataf;
}
~~~

Generic image data.

- Members:
    - width:      Width
    - height:      Height
    - ncomp:      Number of Channels
    - datab:      Buffer data for 8-bit images
    - dataf:      Buffer data for float images


#### Struct glTFid

~~~ .cpp
template <typename T>
struct glTFid {
    glTFid(); 
    explicit glTFid(int id); 
    explicit operator int() const; 
    bool is_valid() const; 
    explicit operator bool() const; 
}
~~~

glTFid

- Members:
    - glTFid():      defaoult constructor to an invalid id
    - glTFid():      explicit conversion from integer
    - operator int():      explicit convcersion to integer
    - is_valid():      check if it is valid
    - operator bool():      check if it is valid


#### Struct glTFProperty

~~~ .cpp
struct glTFProperty {
    map<string, nlohmann::json> extensions = {};
    nlohmann::json extras = {};
}
~~~

Generic glTF object

- Members:
    - extensions:      Extensions.
    - extras:      Extra data.


#### Struct glTFChildOfRootProperty

~~~ .cpp
struct glTFChildOfRootProperty : glTFProperty {
    string name = "";
}
~~~

Generic glTF named object

- Members:
    - name:      The user-defined name of this object.


#### Enum class glTFAccessorSparseIndicesComponentType

~~~ .cpp
enum class glTFAccessorSparseIndicesComponentType {
    NotSet = -1,
}
~~~

Values for glTFAccessorSparseIndices::componentType

- Values:
    - NotSet:      Not set


#### Struct glTFAccessorSparseIndices

~~~ .cpp
struct glTFAccessorSparseIndices : glTFProperty {
    glTFid<glTFBufferView> bufferView = {};
    int byteOffset = 0;
    glTFAccessorSparseIndicesComponentType componentType =
        glTFAccessorSparseIndicesComponentType::NotSet;
}
~~~

Indices of those attributes that deviate from their initialization value.

- Members:
    - bufferView:      The index of the bufferView with sparse indices. Referenced bufferView
     can't have ARRAY_BUFFER or ELEMENT_ARRAY_BUFFER target. [required]
    - byteOffset:      The offset relative to the start of the bufferView in bytes. Must be
     aligned.
    - componentType:      The indices data type. [required]


#### Struct glTFAccessorSparseValues

~~~ .cpp
struct glTFAccessorSparseValues : glTFProperty {
    glTFid<glTFBufferView> bufferView = {};
    int byteOffset = 0;
}
~~~

Array of size `accessor.sparse.count` times number of components storing the
displaced accessor attributes pointed by `accessor.sparse.indices`.

- Members:
    - bufferView:      The index of the bufferView with sparse values. Referenced bufferView
     can't have ARRAY_BUFFER or ELEMENT_ARRAY_BUFFER target. [required]
    - byteOffset:      The offset relative to the start of the bufferView in bytes. Must be
     aligned.


#### Struct glTFAccessorSparse

~~~ .cpp
struct glTFAccessorSparse : glTFProperty {
    int count = 0;
    glTFAccessorSparseIndices* indices = nullptr;
    glTFAccessorSparseValues* values = nullptr;
    ~glTFAccessorSparse(); 
}
~~~

Sparse storage of attributes that deviate from their initialization value.

- Members:
    - count:      Number of entries stored in the sparse array. [required]
    - indices:      Index array of size `count` that points to those accessor attributes
     that deviate from their initialization value. Indices must strictly
     increase. [required]
    - values:      Array of size `count` times number of components, storing the displaced
     accessor attributes pointed by `indices`. Substituted values must have
     the same `componentType` and number of components as the base accessor.
     [required]
    - ~glTFAccessorSparse():      Cleanup


#### Enum class glTFAccessorComponentType

~~~ .cpp
enum class glTFAccessorComponentType {
    NotSet = -1,
}
~~~

Values for glTFAccessor::componentType

- Values:
    - NotSet:      Not set


#### Enum class glTFAccessorType

~~~ .cpp
enum class glTFAccessorType {
    NotSet = -1,
}
~~~

Values for glTFAccessor::type

- Values:
    - NotSet:      Not set


#### Struct glTFAccessor

~~~ .cpp
struct glTFAccessor : glTFChildOfRootProperty {
    glTFid<glTFBufferView> bufferView = {};
    int byteOffset = 0;
    glTFAccessorComponentType componentType = glTFAccessorComponentType::NotSet;
    bool normalized = false;
    int count = 0;
    glTFAccessorType type = glTFAccessorType::NotSet;
    vector<float> max = {};
    vector<float> min = {};
    glTFAccessorSparse* sparse = nullptr;
    ~glTFAccessor(); 
}
~~~

A typed view into a bufferView.  A bufferView contains raw binary data.  An
accessor provides a typed view into a bufferView or a subset of a bufferView
similar to how WebGL's `vertexAttribPointer()` defines an attribute in a
buffer.

- Members:
    - bufferView:      The index of the bufferView.
    - byteOffset:      The offset relative to the start of the bufferView in bytes.
    - componentType:      The datatype of components in the attribute. [required]
    - normalized:      Specifies whether integer data values should be normalized.
    - count:      The number of attributes referenced by this accessor. [required]
    - type:      Specifies if the attribute is a scalar, vector, or matrix. [required]
    - max:      Maximum value of each component in this attribute.
    - min:      Minimum value of each component in this attribute.
    - sparse:      Sparse storage of attributes that deviate from their initialization
     value.
    - ~glTFAccessor():      Cleanup


#### Enum class glTFAnimationChannelTargetPath

~~~ .cpp
enum class glTFAnimationChannelTargetPath {
    NotSet = -1,
}
~~~

Values for glTFAnimationChannelTarget::path

- Values:
    - NotSet:      Not set


#### Struct glTFAnimationChannelTarget

~~~ .cpp
struct glTFAnimationChannelTarget : glTFProperty {
    glTFid<glTFNode> node = {};
    glTFAnimationChannelTargetPath path =
        glTFAnimationChannelTargetPath::NotSet;
}
~~~

The index of the node and TRS property that an animation channel targets.

- Members:
    - node:      The index of the node to target. [required]
    - path:      The name of the node's TRS property to modify, or the "weights" of the
     Morph Targets it instantiates. [required]


#### Struct glTFAnimationChannel

~~~ .cpp
struct glTFAnimationChannel : glTFProperty {
    glTFid<glTFAnimationSampler> sampler = {};
    glTFAnimationChannelTarget* target = nullptr;
    ~glTFAnimationChannel(); 
}
~~~

Targets an animation's sampler at a node's property.

- Members:
    - sampler:      The index of a sampler in this animation used to compute the value for
     the target. [required]
    - target:      The index of the node and TRS property to target. [required]
    - ~glTFAnimationChannel():      Cleanup


#### Enum class glTFAnimationSamplerInterpolation

~~~ .cpp
enum class glTFAnimationSamplerInterpolation {
    NotSet = -1,
}
~~~

Values for glTFAnimationSampler::interpolation

- Values:
    - NotSet:      Not set


#### Struct glTFAnimationSampler

~~~ .cpp
struct glTFAnimationSampler : glTFProperty {
    glTFid<glTFAccessor> input = {};
    glTFAnimationSamplerInterpolation interpolation =
        glTFAnimationSamplerInterpolation::Linear;
    glTFid<glTFAccessor> output = {};
}
~~~

Combines input and output accessors with an interpolation algorithm to
define a keyframe graph (but not its target).

- Members:
    - input:      The index of an accessor containing keyframe input values, e.g., time.
     [required]
    - interpolation:      Interpolation algorithm.
    - output:      The index of an accessor, containing keyframe output values. [required]


#### Struct glTFAnimation

~~~ .cpp
struct glTFAnimation : glTFChildOfRootProperty {
    vector<glTFAnimationChannel*> channels = {};
    vector<glTFAnimationSampler*> samplers = {};
    glTFAnimationChannel* get(const glTFid<glTFAnimationChannel>& id) const; 
    glTFAnimationSampler* get(const glTFid<glTFAnimationSampler>& id) const; 
    ~glTFAnimation(); 
}
~~~

A keyframe animation.

- Members:
    - channels:      An array of channels, each of which targets an animation's sampler at a
     node's property. Different channels of the same animation can't have
     equal targets. [required]
    - samplers:      An array of samplers that combines input and output accessors with an
     interpolation algorithm to define a keyframe graph (but not its target).
     [required]
    - get():      typed access for nodes
    - get():      typed access for nodes
    - ~glTFAnimation():      Cleanup


#### Struct glTFAsset

~~~ .cpp
struct glTFAsset : glTFProperty {
    string copyright = "";
    string generator = "";
    string version = "";
    string minVersion = "";
}
~~~

Metadata about the glTF asset.

- Members:
    - copyright:      A copyright message suitable for display to credit the content creator.
    - generator:      Tool that generated this glTF model.  Useful for debugging.
    - version:      The glTF version that this asset targets. [required]
    - minVersion:      The minimum glTF version that this asset targets.


#### Struct glTFBuffer

~~~ .cpp
struct glTFBuffer : glTFChildOfRootProperty {
    string uri = "";
    int byteLength = 0;
    buffer_data data = {};
}
~~~

A buffer points to binary geometry, animation, or skins.

- Members:
    - uri:      The uri of the buffer.
    - byteLength:      The length of the buffer in bytes. [required]
    - data:      Stores buffer content after loading. [required]


#### Enum class glTFBufferViewTarget

~~~ .cpp
enum class glTFBufferViewTarget {
    NotSet = -1,
}
~~~

Values for glTFBufferView::target

- Values:
    - NotSet:      Not set


#### Struct glTFBufferView

~~~ .cpp
struct glTFBufferView : glTFChildOfRootProperty {
    glTFid<glTFBuffer> buffer = {};
    int byteOffset = 0;
    int byteLength = 0;
    int byteStride = 0;
    glTFBufferViewTarget target = glTFBufferViewTarget::NotSet;
}
~~~

A view into a buffer generally representing a subset of the buffer.

- Members:
    - buffer:      The index of the buffer. [required]
    - byteOffset:      The offset into the buffer in bytes.
    - byteLength:      The length of the bufferView in bytes. [required]
    - byteStride:      The stride, in bytes.
    - target:      The target that the GPU buffer should be bound to.


#### Struct glTFCameraOrthographic

~~~ .cpp
struct glTFCameraOrthographic : glTFProperty {
    float xmag = 0;
    float ymag = 0;
    float zfar = 0;
    float znear = 0;
}
~~~

An orthographic camera containing properties to create an orthographic
projection matrix.

- Members:
    - xmag:      The floating-point horizontal magnification of the view. [required]
    - ymag:      The floating-point vertical magnification of the view. [required]
    - zfar:      The floating-point distance to the far clipping plane. `zfar` must be
     greater than `znear`. [required]
    - znear:      The floating-point distance to the near clipping plane. [required]


#### Struct glTFCameraPerspective

~~~ .cpp
struct glTFCameraPerspective : glTFProperty {
    float aspectRatio = 0;
    float yfov = 0;
    float zfar = 0;
    float znear = 0;
}
~~~

A perspective camera containing properties to create a perspective
projection matrix.

- Members:
    - aspectRatio:      The floating-point aspect ratio of the field of view.
    - yfov:      The floating-point vertical field of view in radians. [required]
    - zfar:      The floating-point distance to the far clipping plane.
    - znear:      The floating-point distance to the near clipping plane. [required]


#### Enum class glTFCameraType

~~~ .cpp
enum class glTFCameraType {
    NotSet = -1,
}
~~~

Values for glTFCamera::type

- Values:
    - NotSet:      Not set


#### Struct glTFCamera

~~~ .cpp
struct glTFCamera : glTFChildOfRootProperty {
    glTFCameraOrthographic* orthographic = nullptr;
    glTFCameraPerspective* perspective = nullptr;
    glTFCameraType type = glTFCameraType::NotSet;
    ~glTFCamera(); 
}
~~~

A camera's projection.  A node can reference a camera to apply a transform
to place the camera in the scene.

- Members:
    - orthographic:      An orthographic camera containing properties to create an orthographic
     projection matrix.
    - perspective:      A perspective camera containing properties to create a perspective
     projection matrix.
    - type:      Specifies if the camera uses a perspective or orthographic projection.
     [required]
    - ~glTFCamera():      Cleanup


#### Enum class glTFImageMimeType

~~~ .cpp
enum class glTFImageMimeType {
    NotSet = -1,
}
~~~

Values for glTFImage::mimeType

- Values:
    - NotSet:      Not set


#### Struct glTFImage

~~~ .cpp
struct glTFImage : glTFChildOfRootProperty {
    string uri = "";
    glTFImageMimeType mimeType = glTFImageMimeType::NotSet;
    glTFid<glTFBufferView> bufferView = {};
    image_data data = {};
}
~~~

Image data used to create a texture. Image can be referenced by URI or
`bufferView` index. `mimeType` is required in the latter case.

- Members:
    - uri:      The uri of the image.
    - mimeType:      The image's MIME type.
    - bufferView:      The index of the bufferView that contains the image. Use this instead of
     the image's uri property.
    - data:      Stores image content after loading.


#### Struct glTFTextureInfo

~~~ .cpp
struct glTFTextureInfo : glTFProperty {
    glTFid<glTFTexture> index = {};
    int texCoord = 0;
}
~~~

Reference to a texture.

- Members:
    - index:      The index of the texture. [required]
    - texCoord:      The set index of texture's TEXCOORD attribute used for texture
     coordinate mapping.


#### Struct glTFTexture

~~~ .cpp
struct glTFTexture : glTFChildOfRootProperty {
    glTFid<glTFSampler> sampler = {};
    glTFid<glTFImage> source = {};
}
~~~

A texture and its sampler.

- Members:
    - sampler:      The index of the sampler used by this texture. When undefined, a sampler
     with repeat wrapping and auto filtering should be used.
    - source:      The index of the image used by this texture.


#### Struct glTFMaterialNormalTextureInfo

~~~ .cpp
struct glTFMaterialNormalTextureInfo : glTFTextureInfo {
    float scale = 1;
}
~~~

Normal texture information.

- Members:
    - scale:      The scalar multiplier applied to each normal vector of the normal
     texture.


#### Struct glTFMaterialOcclusionTextureInfo

~~~ .cpp
struct glTFMaterialOcclusionTextureInfo : glTFTextureInfo {
    float strength = 1;
}
~~~

Occlusion texture information.

- Members:
    - strength:      A scalar multiplier controlling the amount of occlusion applied.


#### Struct glTFMaterialPbrMetallicRoughness

~~~ .cpp
struct glTFMaterialPbrMetallicRoughness : glTFProperty {
    vec4f baseColorFactor = {1, 1, 1, 1};
    glTFTextureInfo* baseColorTexture = nullptr;
    float metallicFactor = 1;
    float roughnessFactor = 1;
    glTFTextureInfo* metallicRoughnessTexture = nullptr;
    ~glTFMaterialPbrMetallicRoughness(); 
}
~~~

A set of parameter values that are used to define the metallic-roughness
material model from Physically-Based Rendering (PBR) methodology.

- Members:
    - baseColorFactor:      The material's base color factor.
    - baseColorTexture:      The base color texture.
    - metallicFactor:      The metalness of the material.
    - roughnessFactor:      The roughness of the material.
    - metallicRoughnessTexture:      The metallic-roughness texture.
    - ~glTFMaterialPbrMetallicRoughness():      Cleanup


#### Struct glTFMaterialPbrSpecularGlossiness

~~~ .cpp
struct glTFMaterialPbrSpecularGlossiness : glTFProperty {
    vec4f diffuseFactor = {1, 1, 1, 1};
    glTFTextureInfo* diffuseTexture = nullptr;
    vec3f specularFactor = {1, 1, 1};
    float glossinessFactor = 1;
    glTFTextureInfo* specularGlossinessTexture = nullptr;
    ~glTFMaterialPbrSpecularGlossiness(); 
}
~~~

glTF extension that defines the specular-glossiness material model from
Physically-Based Rendering (PBR) methodology.

- Members:
    - diffuseFactor:      The reflected diffuse factor of the material.
    - diffuseTexture:      The diffuse texture.
    - specularFactor:      The specular RGB color of the material.
    - glossinessFactor:      The glossiness or smoothness of the material.
    - specularGlossinessTexture:      The specular-glossiness texture.
    - ~glTFMaterialPbrSpecularGlossiness():      Cleanup


#### Enum class glTFMaterialAlphaMode

~~~ .cpp
enum class glTFMaterialAlphaMode {
    NotSet = -1,
}
~~~

Values for glTFMaterial::alphaMode

- Values:
    - NotSet:      Not set


#### Struct glTFMaterial

~~~ .cpp
struct glTFMaterial : glTFChildOfRootProperty {
    glTFMaterialPbrMetallicRoughness* pbrMetallicRoughness = nullptr;
    glTFMaterialPbrSpecularGlossiness* pbrSpecularGlossiness = nullptr;
    glTFMaterialNormalTextureInfo* normalTexture = nullptr;
    glTFMaterialOcclusionTextureInfo* occlusionTexture = nullptr;
    glTFTextureInfo* emissiveTexture = nullptr;
    vec3f emissiveFactor = {0, 0, 0};
    glTFMaterialAlphaMode alphaMode = glTFMaterialAlphaMode::Opaque;
    float alphaCutoff = 0.5;
    bool doubleSided = false;
    ~glTFMaterial(); 
}
~~~

The material appearance of a primitive.

- Members:
    - pbrMetallicRoughness:      A set of parameter values that are used to define the metallic-roughness
     material model from Physically-Based Rendering (PBR) methodology. When
     not specified, all the default values of `pbrMetallicRoughness` apply.
    - pbrSpecularGlossiness:      A set of parameter values that are used to define the
     specular-glossiness material model from Physically-Based Rendering (PBR)
     methodology. When not specified, all the default values of
     `pbrMetallicRoughness` apply.
    - normalTexture:      The normal map texture.
    - occlusionTexture:      The occlusion map texture.
    - emissiveTexture:      The emissive map texture.
    - emissiveFactor:      The emissive color of the material.
    - alphaMode:      The alpha rendering mode of the material.
    - alphaCutoff:      The alpha cutoff value of the material.
    - doubleSided:      Specifies whether the material is double sided.
    - ~glTFMaterial():      Cleanup


#### Enum class glTFMeshPrimitiveMode

~~~ .cpp
enum class glTFMeshPrimitiveMode {
    NotSet = -1,
}
~~~

Values for glTFMeshPrimitive::mode

- Values:
    - NotSet:      Not set


#### Struct glTFMeshPrimitive

~~~ .cpp
struct glTFMeshPrimitive : glTFProperty {
    map<string, glTFid<glTFAccessor>> attributes = {};
    glTFid<glTFAccessor> indices = {};
    glTFid<glTFMaterial> material = {};
    glTFMeshPrimitiveMode mode = glTFMeshPrimitiveMode::Triangles;
    vector<map<string, glTFid<glTFAccessor>>> targets = {};
}
~~~

Geometry to be rendered with the given material.

- Members:
    - attributes:      A dictionary object, where each key corresponds to mesh attribute
     semantic and each value is the index of the accessor containing
     attribute's data. [required]
    - indices:      The index of the accessor that contains the indices.
    - material:      The index of the material to apply to this primitive when rendering.
    - mode:      The type of primitives to render.
    - targets:      An array of Morph Targets, each  Morph Target is a dictionary mapping
     attributes (only `POSITION`, `NORMAL`, and `TANGENT` supported) to their
     deviations in the Morph Target.


#### Struct glTFMesh

~~~ .cpp
struct glTFMesh : glTFChildOfRootProperty {
    vector<glTFMeshPrimitive*> primitives = {};
    vector<float> weights = {};
    ~glTFMesh(); 
}
~~~

A set of primitives to be rendered.  A node can contain one mesh.  A node's
transform places the mesh in the scene.

- Members:
    - primitives:      An array of primitives, each defining geometry to be rendered with a
     material. [required]
    - weights:      Array of weights to be applied to the Morph Targets.
    - ~glTFMesh():      Cleanup


#### Struct glTFNode

~~~ .cpp
struct glTFNode : glTFChildOfRootProperty {
    glTFid<glTFCamera> camera = {};
    vector<glTFid<glTFNode>> children = {};
    glTFid<glTFSkin> skin = {};
    mat4f matrix = { {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1} };
    glTFid<glTFMesh> mesh = {};
    quat4f rotation = {0, 0, 0, 1};
    vec3f scale = {1, 1, 1};
    vec3f translation = {0, 0, 0};
    vector<float> weights = {};
}
~~~

A node in the node hierarchy.  When the node contains `skin`, all
`mesh.primitives` must contain `JOINTS_0` and `WEIGHTS_0` attributes.  A
node can have either a `matrix` or any combination of
`translation`/`rotation`/`scale` (TRS) properties. TRS properties are
converted to matrices and postmultiplied in the `T * R * S` order to compose
the transformation matrix; first the scale is applied to the vertices, then
the rotation, and then the translation. If none are provided, the transform
is the identity. When a node is targeted for animation (referenced by an
animation.channel.target), only TRS properties may be present; `matrix` will
not be present.

- Members:
    - camera:      The index of the camera referenced by this node.
    - children:      The indices of this node's children.
    - skin:      The index of the skin referenced by this node.
    - matrix:      A floating-point 4x4 transformation matrix stored in column-major order.
    - mesh:      The index of the mesh in this node.
    - rotation:      The node's unit quaternion rotation in the order (x, y, z, w), where w
     is the scalar.
    - scale:      The node's non-uniform scale.
    - translation:      The node's translation.
    - weights:      The weights of the instantiated Morph Target. Number of elements must
     match number of Morph Targets of used mesh.


#### Enum class glTFSamplerMagFilter

~~~ .cpp
enum class glTFSamplerMagFilter {
    NotSet = -1,
}
~~~

Values for glTFSampler::magFilter

- Values:
    - NotSet:      Not set


#### Enum class glTFSamplerMinFilter

~~~ .cpp
enum class glTFSamplerMinFilter {
    NotSet = -1,
}
~~~

Values for glTFSampler::minFilter

- Values:
    - NotSet:      Not set


#### Enum class glTFSamplerWrapS

~~~ .cpp
enum class glTFSamplerWrapS {
    NotSet = -1,
}
~~~

glTFSampler::wrapS

- Values:
    - NotSet:      Not set


#### Enum class glTFSamplerWrapT

~~~ .cpp
enum class glTFSamplerWrapT {
    NotSet = -1,
}
~~~

glTFSampler::wrapT

- Values:
    - NotSet:      Not set


#### Struct glTFSampler

~~~ .cpp
struct glTFSampler : glTFChildOfRootProperty {
    glTFSamplerMagFilter magFilter = glTFSamplerMagFilter::NotSet;
    glTFSamplerMinFilter minFilter = glTFSamplerMinFilter::NotSet;
    glTFSamplerWrapS wrapS = glTFSamplerWrapS::Repeat;
    glTFSamplerWrapT wrapT = glTFSamplerWrapT::Repeat;
}
~~~

Texture sampler properties for filtering and wrapping modes.

- Members:
    - magFilter:      Magnification filter.
    - minFilter:      Minification filter.
    - wrapS:      s wrapping mode.
    - wrapT:      t wrapping mode.


#### Struct glTFScene

~~~ .cpp
struct glTFScene : glTFChildOfRootProperty {
    vector<glTFid<glTFNode>> nodes = {};
}
~~~

The root nodes of a scene.

- Members:
    - nodes:      The indices of each root node.


#### Struct glTFSkin

~~~ .cpp
struct glTFSkin : glTFChildOfRootProperty {
    glTFid<glTFAccessor> inverseBindMatrices = {};
    glTFid<glTFNode> skeleton = {};
    vector<glTFid<glTFNode>> joints = {};
}
~~~

Joints and matrices defining a skin.

- Members:
    - inverseBindMatrices:      The index of the accessor containing the floating-point 4x4 inverse-bind
     matrices.  The default is that each matrix is a 4x4 identity matrix,
     which implies that inverse-bind matrices were pre-applied.
    - skeleton:      The index of the node used as a skeleton root. When undefined, joints
     transforms resolve to scene root.
    - joints:      Indices of skeleton nodes, used as joints in this skin. [required]


#### Struct glTF

~~~ .cpp
struct glTF : glTFProperty {
    vector<string> extensionsUsed = {};
    vector<string> extensionsRequired = {};
    vector<glTFAccessor*> accessors = {};
    vector<glTFAnimation*> animations = {};
    glTFAsset* asset = nullptr;
    vector<glTFBuffer*> buffers = {};
    vector<glTFBufferView*> bufferViews = {};
    vector<glTFCamera*> cameras = {};
    vector<glTFImage*> images = {};
    vector<glTFMaterial*> materials = {};
    vector<glTFMesh*> meshes = {};
    vector<glTFNode*> nodes = {};
    vector<glTFSampler*> samplers = {};
    glTFid<glTFScene> scene = {};
    vector<glTFScene*> scenes = {};
    vector<glTFSkin*> skins = {};
    vector<glTFTexture*> textures = {};
    glTFAccessor* get(const glTFid<glTFAccessor>& id) const; 
    glTFAnimation* get(const glTFid<glTFAnimation>& id) const; 
    glTFBuffer* get(const glTFid<glTFBuffer>& id) const; 
    glTFBufferView* get(const glTFid<glTFBufferView>& id) const; 
    glTFCamera* get(const glTFid<glTFCamera>& id) const; 
    glTFImage* get(const glTFid<glTFImage>& id) const; 
    glTFMaterial* get(const glTFid<glTFMaterial>& id) const; 
    glTFMesh* get(const glTFid<glTFMesh>& id) const; 
    glTFNode* get(const glTFid<glTFNode>& id) const; 
    glTFSampler* get(const glTFid<glTFSampler>& id) const; 
    glTFScene* get(const glTFid<glTFScene>& id) const; 
    glTFSkin* get(const glTFid<glTFSkin>& id) const; 
    glTFTexture* get(const glTFid<glTFTexture>& id) const; 
    ~glTF(); 
}
~~~

The root object for a glTF asset.

- Members:
    - extensionsUsed:      Names of glTF extensions used somewhere in this asset.
    - extensionsRequired:      Names of glTF extensions required to properly load this asset.
    - accessors:      An array of accessors.
    - animations:      An array of keyframe animations.
    - asset:      Metadata about the glTF asset. [required]
    - buffers:      An array of buffers.
    - bufferViews:      An array of bufferViews.
    - cameras:      An array of cameras.
    - images:      An array of images.
    - materials:      An array of materials.
    - meshes:      An array of meshes.
    - nodes:      An array of nodes.
    - samplers:      An array of samplers.
    - scene:      The index of the default scene.
    - scenes:      An array of scenes.
    - skins:      An array of skins.
    - textures:      An array of textures.
    - get():      typed access for nodes
    - get():      typed access for nodes
    - get():      typed access for nodes
    - get():      typed access for nodes
    - get():      typed access for nodes
    - get():      typed access for nodes
    - get():      typed access for nodes
    - get():      typed access for nodes
    - get():      typed access for nodes
    - get():      typed access for nodes
    - get():      typed access for nodes
    - get():      typed access for nodes
    - get():      typed access for nodes
    - ~glTF():      Cleanup


#### Function load_gltf()

~~~ .cpp
glTF* load_gltf(const string& filename, bool load_bin = true,
    bool load_img = false, bool skip_missing = false);
~~~

Loads a gltf file from disk

- Parameters:
    - filename: scene filename
    - load_bin/load_img: load binary data
    - skip_missing: do not throw an exception if a file is missing
- Returns:
    - gltf data loaded (nullptr on error)

#### Function load_binary_gltf()

~~~ .cpp
glTF* load_binary_gltf(const string& filename, bool load_bin = true,
    bool load_img = false, bool skip_missing = false);
~~~

Loads a binary gltf file from disk

- Parameters:
    - filename: scene filename
    - other params as above
- Returns:
    - gltf data loaded (nullptr on error)

#### Function save_gltf()

~~~ .cpp
void save_gltf(const string& filename, const glTF* gltf, bool save_bin = true,
    bool save_images = false);
~~~

Saves a scene to disk

- Parameters:
    - filename: scene filename
    - gltf: data to save
    - save_bin/save_images: save binary data

#### Function save_binary_gltf()

~~~ .cpp
void save_binary_gltf(const string& filename, const glTF* gltf,
    bool save_bin = true, bool save_images = false);
~~~

Saves a scene to disk

- Parameters:
    - filename: scene filename
    - gltf: data to save
    - save_bin/save_images: save binary data

#### Function node_transform()

~~~ .cpp
inline mat4f node_transform(const glTFNode* node);
~~~

Computes the local node transform and its inverse.

#### Struct accessor_view

~~~ .cpp
struct accessor_view {
    accessor_view(const glTF* gltf, const glTFAccessor* accessor); 
    int size() const; 
    int count() const; 
    int ncomp() const; 
    bool valid() const; 
    vec2f getv2f(int idx, const vec2f& def =; 
    vec3f getv3f(int idx, const vec3f& def =; 
    vec4f getv4f(int idx, const vec4f& def =; 
    mat4f getm4f(int idx) const; 
    float get(int idx, int c = 0) const; 
    vec2i getv2i(int idx, const vec2i& def =; 
    vec3i getv3i(int idx, const vec3i& def =; 
    vec4i getv4i(int idx, const vec4i& def =; 
    int geti(int idx, int c = 0) const; 
}
~~~

A view for gltf array buffers that allows for typed access.

- Members:
    - accessor_view():      construct a view from an accessor
    - size():      number of elements in the view
    - count():      number of elements in the view
    - ncomp():      number of components per element
    - valid():      check whether the view is valid
    - getv2f():      get the idx-th element of fixed length width default values
    - getv3f():      get the idx-th element of fixed length width default values
    - getv4f():      get the idx-th element of fixed length width default values
    - getm4f():      get the idx-th element of fixed length as a matrix
    - get():      get the c-th component of the idx-th element
    - getv2i():      get the idx-th element as integer with fixed length
    - getv3i():      get the idx-th element as integer with fixed length
    - getv4i():      get the idx-th element as integer with fixed length
    - geti():      get the c-th component of the idx-th element as integer


#### Function startswith()

~~~ .cpp
inline bool startswith(const string& str, const string& substr);
~~~

Checks if a string starts with a prefix.

#### Function endswith()

~~~ .cpp
inline bool endswith(const string& str, const string& substr);
~~~

Checks if a string ends with a prefix.

#### Function contains()

~~~ .cpp
inline bool contains(const string& str, const string& substr);
~~~

Check is a string contains a substring.

#### Function splitlines()

~~~ .cpp
inline vector<string> splitlines(const string& str, bool keep_newline = false);
~~~

Splits a string into lines at the '\n' character. The line
terminator is kept if keep_newline. This function does not work on
Window if keep_newline is true.

#### Function partition()

~~~ .cpp
inline vector<string> partition(const string& str, const string& split);
~~~

Partition the string.

#### Function split()

~~~ .cpp
inline vector<string> split(const string& str);
~~~

Splits the string.

#### Function split()

~~~ .cpp
inline vector<string> split(const string& str, const string& substr);
~~~

Splits the string.

#### Function split()

~~~ .cpp
inline vector<string> split(const string& str, char substr);
~~~

Splits the string.

#### Function rstrip()

~~~ .cpp
inline string rstrip(const string& str);
~~~

Strip the string.

#### Function lstrip()

~~~ .cpp
inline string lstrip(const string& str);
~~~

Strip the string.

#### Function strip()

~~~ .cpp
inline string strip(const string& str);
~~~

Strip the string.

#### Function join()

~~~ .cpp
inline string join(const vector<string>& strs, const string& sep);
~~~

Joins a list of string with a string as separator.

#### Function lower()

~~~ .cpp
inline string lower(const string& str);
~~~

Converts an ASCII string to lowercase.

#### Function upper()

~~~ .cpp
inline string upper(const string& str);
~~~

Converts an ASCII string to uppercase.

#### Function isspace()

~~~ .cpp
inline bool isspace(const string& str);
~~~

Check if a string is space.

#### Function replace()

~~~ .cpp
inline string replace(const string& str, const string& s1, const string& s2);
~~~

Replace s1 with s2 in str.

#### Function path_dirname()

~~~ .cpp
inline string path_dirname(const string& filename);
~~~

Get directory name (including '/').

#### Function path_extension()

~~~ .cpp
inline string path_extension(const string& filename);
~~~

Get extension (including '.').

#### Function path_basename()

~~~ .cpp
inline string path_basename(const string& filename);
~~~

Get file basename.

#### Function path_filename()

~~~ .cpp
inline string path_filename(const string& filename);
~~~

Get filename without directory (equiv to get_basename() +
get_extension()).

#### Function replace_path_extension()

~~~ .cpp
inline string replace_path_extension(
    const string& filename, const string& ext);
~~~

Replace extension.

#### Function prepend_path_extension()

~~~ .cpp
inline string prepend_path_extension(
    const string& filename, const string& prep);
~~~

Prepend a string to the extension.

#### Function split_path()

~~~ .cpp
inline void split_path(
    const string& filename, string& dirname, string& basename, string& ext);
~~~

Splits a path calling the above functions.

#### Function format()

~~~ .cpp
inline string format(const string& fmt, const vector<string>& args);
~~~

Really-minimal Python like string format. The implementation is not fast
nor memory efficient. But it is good enough for some needs.

#### Function format()

~~~ .cpp
template <typename... Args>
inline string format(const string& fmt, const Args&... args);
~~~

Really-minimal Python like string format. Internally uses streams for
generality and supports for now only the '{}' operator. The implementation
is not fast nor memory efficient. But it is good enough for some needs.

#### Function print()

~~~ .cpp
template <typename... Args>
inline void print(const string& fmt, const Args&... args);
~~~

Wrapper for the above function that prints to stdout.

#### Function println()

~~~ .cpp
template <typename... Args>
inline void println(const string& fmt, const Args&... args);
~~~

Wrapper for the above function that prints to stdout with endline.

#### Function load_binfile()

~~~ .cpp
inline vector<unsigned char> load_binfile(const string& filename);
~~~

Loads the contents of a binary file in an in-memory array.

#### Function load_txtfile()

~~~ .cpp
inline string load_txtfile(const string& filename);
~~~

Loads the contents of a text file into a string.

#### Function save_binfile()

~~~ .cpp
inline void save_binfile(
    const string& filename, const vector<unsigned char>& data);
~~~

Saves binary data to a file.

#### Function save_txtfile()

~~~ .cpp
inline void save_txtfile(const string& filename, const string& str);
~~~

Saves a string to a text file.

#### Struct cmdline_parser

~~~ .cpp
struct cmdline_parser;
~~~

Immediate mode command line parser (opaque type)

#### Struct cmdline_parser

~~~ .cpp
struct cmdline_parser {
~~~

Immediate mode command line parser

#### Function should_exit()

~~~ .cpp
inline bool should_exit(cmdline_parser& parser);
~~~

check unused arguments

#### Function get_usage()

~~~ .cpp
inline string get_usage(const cmdline_parser& parser);
~~~

returns the usage string

#### Function parse_flag()

~~~ .cpp
inline bool parse_flag(cmdline_parser& parser, const string& name,
    const string& flag, const string& help, bool def = false,
    bool req = false);
~~~

parse a flag from the command line

#### Function parse_opt()

~~~ .cpp
template <typename T>
inline T parse_opt(cmdline_parser& parser, const string& name,
    const string& flag, const string& help, const T& def =;
~~~

parse an option from the command line

#### Function parse_opt()

~~~ .cpp
template <typename T>
inline T parse_opt(cmdline_parser& parser, const string& name,
    const string& flag, const string& help,
    const vector<pair<string, T>>& key_values, const T& def, bool req = false,
    const vector<T>& choices =;
~~~

parse an enum option from the command line

#### Function make_parser()

~~~ .cpp
inline cmdline_parser make_parser(
    int argc, char** argv, const string& prog, const string& help);
~~~

initialize the command line

#### Struct logger

~~~ .cpp
struct logger {
    bool _verbose = true;
    bool _console = true;
    FILE* _file = nullptr;
}
~~~

Logger object. A logger can output messages to multiple streams.
Use add streams commands for it.

- Members:
    - _verbose:      whether to output verbose
    - _console:      whether to output to console
    - _file:      file stream for stream output


#### Function make_logger()

~~~ .cpp
inline logger* make_logger(bool console = true, bool verbose = true);
~~~

Make a logger with an optional console stream and a verbosity level

#### Function add_file_stream()

~~~ .cpp
inline void add_file_stream(logger* lgr, const string& filename, bool append);
~~~

Add a file stream to a logger.

- Parameters:
    - lgr: logger
    - filename: filename
    - append: append or write open mode for file logger
    - short_message: whether to use a short message version
    - output_level: output level
    - flush_level: output level
- Returns:
    - true if ok

#### Function get_default_logger()

~~~ .cpp
inline logger* get_default_logger();
~~~

Get default logger.
By default a non-verbose stdout logger is creater.

#### Function log_info()

~~~ .cpp
template <typename... Args>
inline void log_info(logger* lgr, const string& msg, const Args&... args);
~~~

Log a info message

#### Function log_warning()

~~~ .cpp
template <typename... Args>
inline void log_warning(logger* lgr, const string& msg, const Args&... args);
~~~

Log a info message

#### Function log_error()

~~~ .cpp
template <typename... Args>
inline void log_error(logger* lgr, const string& msg, const Args&... args);
~~~

Log an error message

#### Function log_fatal()

~~~ .cpp
template <typename... Args>
inline void log_fatal(logger* lgr, const string& msg, const Args&... args);
~~~

Log a fatal message and exit

#### Function add_file_stream()

~~~ .cpp
inline void add_file_stream(const string& filename, bool append);
~~~

Adds a file stream to the default logger

#### Function log_info()

~~~ .cpp
template <typename... Args>
inline void log_info(const string& msg, const Args&... args);
~~~

Logs a message to the default loggers

#### Function log_error()

~~~ .cpp
template <typename... Args>
inline void log_error(const string& msg, const Args&... args);
~~~

Logs a message to the default loggers

#### Function log_fatal()

~~~ .cpp
template <typename... Args>
inline void log_fatal(const string& msg, const Args&... args);
~~~

Logs a message to the default loggers

#### Struct thread_pool

~~~ .cpp
struct thread_pool {
~~~

Thread pool for concurrency. This code is derived from LLVM ThreadPool

#### Function make_pool()

~~~ .cpp
inline thread_pool* make_pool(
    int nthreads = std::thread::hardware_concurrency());
~~~

Makes a thread pool

#### Function run_async()

~~~ .cpp
inline std::shared_future<void> run_async(
    thread_pool* pool, const function<void()>& task);
~~~

Runs a task asynchronously onto the global thread pool

#### Function wait_pool()

~~~ .cpp
inline void wait_pool(thread_pool* pool);
~~~

Wait for all jobs to finish on the global thread pool

#### Function clear_pool()

~~~ .cpp
inline void clear_pool(thread_pool* pool);
~~~

Clear all jobs on the global thread pool

#### Function parallel_for()

~~~ .cpp
inline void parallel_for(
    thread_pool* pool, int count, const function<void(int idx)>& task);
~~~

Parallel for implementation on the global thread pool

#### Function get_global_pool()

~~~ .cpp
inline thread_pool* get_global_pool();
~~~

Global pool

#### Function run_async()

~~~ .cpp
inline std::shared_future<void> run_async(const function<void()>& task);
~~~

Runs a task asynchronously onto the global thread pool

#### Function wait_pool()

~~~ .cpp
inline void wait_pool();
~~~

Wait for all jobs to finish on the global thread pool

#### Function clear_pool()

~~~ .cpp
inline void clear_pool();
~~~

Clear all jobs on the global thread pool

#### Function parallel_for()

~~~ .cpp
inline void parallel_for(int count, const function<void(int idx)>& task);
~~~

Parallel for implementation on the global thread pool

#### Struct timer

~~~ .cpp
struct timer {
    timer(bool autostart = true); 
    void start(); 
    void stop(); 
    double elapsed(); 
}
~~~

A simple wrapper for std::chrono.

- Members:
    - timer():      initialize a timer and start it if necessary
    - start():      start a timer
    - stop():      stops a timer
    - elapsed():      elapsed time


#### Enum gl_etype : int

~~~ .cpp
enum struct gl_etype : int {
    point = 1,
    line = 2,
    triangle = 3,
    quad = 4,
}
~~~

Shape types

- Values:
    - point:      points
    - line:      lines
    - triangle:      triangles
    - quad:      quads


#### Enum gl_ltype : int

~~~ .cpp
enum struct gl_ltype : int {
    point = 0,
    directional = 1,
}
~~~

Light types

- Values:
    - point:      point lights
    - directional:      directional lights


#### Function gl_check_error()

~~~ .cpp
bool gl_check_error(bool print = true);
~~~

Checks for GL error and then prints

#### Function gl_clear_buffers()

~~~ .cpp
void gl_clear_buffers(const vec4f& background =;
~~~

Clear window

#### Function gl_enable_depth_test()

~~~ .cpp
void gl_enable_depth_test(bool enabled);
~~~

Enable/disable depth test

#### Function gl_enable_culling()

~~~ .cpp
void gl_enable_culling(bool enabled);
~~~

Enable/disable culling

#### Function gl_enable_wireframe()

~~~ .cpp
void gl_enable_wireframe(bool enabled);
~~~

Enable/disable wireframe

#### Function gl_enable_edges()

~~~ .cpp
void gl_enable_edges(bool enabled, float tolerance = 0.9999f);
~~~

Enable/disable edges. Attempts to avoid z-fighting but the method is not
robust.

#### Function gl_enable_blending()

~~~ .cpp
void gl_enable_blending(bool enabled);
~~~

Enable/disable blending

#### Function gl_set_blend_over()

~~~ .cpp
void gl_set_blend_over();
~~~

Set blending to over operator

#### Function gl_line_width()

~~~ .cpp
void gl_line_width(float w);
~~~

Line width

#### Function gl_set_viewport()

~~~ .cpp
void gl_set_viewport(const vec4i& v);
~~~

Set viewport

#### Function gl_set_viewport()

~~~ .cpp
void gl_set_viewport(const vec2i& v);
~~~

Set viewport

#### Struct gl_texture

~~~ .cpp
struct gl_texture {
~~~

Opengl texture object

#### Function make_texture()

~~~ .cpp
inline gl_texture make_texture(int w, int h, int nc, const float* pixels,
    bool linear, bool mipmap, bool as_float);
~~~

Creates a texture with pixels values of size w, h with nc number of
components (1-4).
Internally use float if as_float and filtering if filter.
Returns the texture id.

#### Function make_texture()

~~~ .cpp
inline gl_texture make_texture(int w, int h, int nc,
    const unsigned char* pixels, bool linear, bool mipmap, bool as_srgb);
~~~

Creates a texture with pixels values of size w, h with nc number of
components (1-4).
Internally use srgb lookup if as_srgb and filtering if filter.
Returns the texture id.

#### Function make_texture()

~~~ .cpp
inline gl_texture make_texture(
    const image4f& img, bool linear, bool mipmap, bool as_float);
~~~

Creates a texture from an image.
Internally use float if as_float and filtering if filter.
Returns the texture id.

#### Function make_texture()

~~~ .cpp
inline gl_texture make_texture(
    const image4b& img, bool linear, bool mipmap, bool as_srgb);
~~~

Creates a texture from an image.
Internally use srgb lookup if as_srgb and filtering if filter.
Returns the texture id.

#### Function update_texture()

~~~ .cpp
inline void update_texture(
    gl_texture& txt, int w, int h, int nc, const float* pixels);
~~~

Updates the texture tid with new image data.

#### Function update_texture()

~~~ .cpp
inline void update_texture(
    gl_texture& txt, int w, int h, int nc, const unsigned char* pixels);
~~~

Updates the texture tid with new image data.

#### Function update_texture()

~~~ .cpp
inline void update_texture(gl_texture& txt, const image4f& img);
~~~

Updates the texture tid with new image data.

#### Function update_texture()

~~~ .cpp
inline void update_texture(gl_texture& txt, const image4b& img);
~~~

Updates the texture tid with new image data.

#### Function bind_texture()

~~~ .cpp
void bind_texture(const gl_texture& txt, uint unit);
~~~

Binds a texture to a texture unit

#### Function unbind_texture()

~~~ .cpp
void unbind_texture(const gl_texture& txt, uint unit);
~~~

Unbinds

#### Function get_texture_id()

~~~ .cpp
inline uint get_texture_id(const gl_texture& txt);
~~~

Get id

#### Function is_texture_valid()

~~~ .cpp
inline bool is_texture_valid(const gl_texture& txt);
~~~

Check if defined

#### Function clear_texture()

~~~ .cpp
void clear_texture(gl_texture& txt);
~~~

Destroys the texture tid.

#### Enum gl_texture_wrap

~~~ .cpp
enum struct gl_texture_wrap {
    not_set = 0,
    repeat = 1,
    clamp = 2,
    mirror = 3,
}
~~~

Wrap values for texture

- Values:
    - not_set:      not set
    - repeat:      repeat
    - clamp:      clamp to edge
    - mirror:      mirror


#### Enum gl_texture_filter

~~~ .cpp
enum struct gl_texture_filter {
    not_set = 0,
    linear = 1,
    nearest = 2,
    linear_mipmap_linear = 3,
    nearest_mipmap_nearest = 4,
    linear_mipmap_nearest = 5,
    nearest_mipmap_linear = 6,
}
~~~

Filter values for texture

- Values:
    - not_set:      not set
    - linear:      linear
    - nearest:      nearest
    - linear_mipmap_linear:      mip-mapping
    - nearest_mipmap_nearest:      mip-mapping
    - linear_mipmap_nearest:      mip-mapping
    - nearest_mipmap_linear:      mip-mapping


#### Struct gl_texture_info

~~~ .cpp
struct gl_texture_info {
    gl_texture txt = {};
    int texcoord = 0;
    float scale = 1;
    gl_texture_wrap wrap_s = gl_texture_wrap::not_set;
    gl_texture_wrap wrap_t = gl_texture_wrap::not_set;
    gl_texture_filter filter_mag = gl_texture_filter::not_set;
    gl_texture_filter filter_min = gl_texture_filter::not_set;
    gl_texture_info(); 
    gl_texture_info(const gl_texture& tid); 
}
~~~

Texture information for parameter setting.

- Members:
    - txt:      texture
    - texcoord:      texture coordinate set
    - scale:      texture strength/scale (used by some models)
    - wrap_s:      wrap mode
    - wrap_t:      wrap mode
    - filter_mag:      filter mode
    - filter_min:      filter mode
    - gl_texture_info():      default constructor
    - gl_texture_info():      constructor from texture id only


#### Struct gl_vertex_buffer

~~~ .cpp
struct gl_vertex_buffer {
~~~

OpenGL vertex/element buffer

#### Function make_vertex_buffer()

~~~ .cpp
inline gl_vertex_buffer make_vertex_buffer(
    int num, int ncomp, const float* values, bool dynamic = false);
~~~

Creates a buffer.

#### Function make_vertex_buffer()

~~~ .cpp
inline gl_vertex_buffer make_vertex_buffer(
    int num, int ncomp, const int* values, bool dynamic = false);
~~~

Creates a buffer.

#### Function make_vertex_buffer()

~~~ .cpp
inline gl_vertex_buffer make_vertex_buffer(
    const vector<float>& values, bool dynamic = false);
~~~

Creates a buffer.

#### Function make_vertex_buffer()

~~~ .cpp
inline gl_vertex_buffer make_vertex_buffer(
    const vector<vec2f>& values, bool dynamic = false);
~~~

Creates a buffer.

#### Function make_vertex_buffer()

~~~ .cpp
inline gl_vertex_buffer make_vertex_buffer(
    const vector<vec3f>& values, bool dynamic = false);
~~~

Creates a buffer.

#### Function make_vertex_buffer()

~~~ .cpp
inline gl_vertex_buffer make_vertex_buffer(
    const vector<vec4f>& values, bool dynamic = false);
~~~

Creates a buffer.

#### Function make_vertex_buffer()

~~~ .cpp
inline gl_vertex_buffer make_vertex_buffer(
    const vector<int>& values, bool dynamic = false);
~~~

Creates a buffer.

#### Function make_vertex_buffer()

~~~ .cpp
inline gl_vertex_buffer make_vertex_buffer(
    const vector<vec2i>& values, bool dynamic = false);
~~~

Creates a buffer.

#### Function make_vertex_buffer()

~~~ .cpp
inline gl_vertex_buffer make_vertex_buffer(
    const vector<vec3i>& values, bool dynamic = false);
~~~

Creates a buffer.

#### Function make_vertex_buffer()

~~~ .cpp
inline gl_vertex_buffer make_vertex_buffer(
    const vector<vec4i>& values, bool dynamic = false);
~~~

Creates a buffer.

#### Function update_vertex_buffer()

~~~ .cpp
inline void update_vertex_buffer(
    gl_vertex_buffer& buf, int num, int ncomp, const float* values);
~~~

Updates the buffer with new data.

#### Function update_vertex_buffer()

~~~ .cpp
inline void update_vertex_buffer(
    gl_vertex_buffer& buf, int num, int ncomp, const int* values);
~~~

Updates the buffer with new data.

#### Function update_vertex_buffer()

~~~ .cpp
inline void update_vertex_buffer(
    gl_vertex_buffer& buf, const vector<float>& values);
~~~

Updates the buffer bid with new data.

#### Function update_vertex_buffer()

~~~ .cpp
inline void update_vertex_buffer(
    gl_vertex_buffer& buf, const vector<vec2f>& values);
~~~

Updates the buffer bid with new data.

#### Function update_vertex_buffer()

~~~ .cpp
inline void update_vertex_buffer(
    gl_vertex_buffer& buf, const vector<vec3f>& values);
~~~

Updates the buffer bid with new data.

#### Function update_vertex_buffer()

~~~ .cpp
inline void update_vertex_buffer(
    gl_vertex_buffer& buf, const vector<vec4f>& values);
~~~

Updates the buffer bid with new data.

#### Function update_vertex_buffer()

~~~ .cpp
inline void update_vertex_buffer(
    gl_vertex_buffer& buf, const vector<int>& values);
~~~

Updates the buffer bid with new data.

#### Function update_vertex_buffer()

~~~ .cpp
inline void update_vertex_buffer(
    gl_vertex_buffer& buf, const vector<vec2i>& values);
~~~

Updates the buffer bid with new data.

#### Function update_vertex_buffer()

~~~ .cpp
inline void update_vertex_buffer(
    gl_vertex_buffer& buf, const vector<vec3i>& values);
~~~

Updates the buffer bid with new data.

#### Function update_vertex_buffer()

~~~ .cpp
inline void update_vertex_buffer(
    gl_vertex_buffer& buf, const vector<vec4i>& values);
~~~

Updates the buffer bid with new data.

#### Function bind_vertex_buffer()

~~~ .cpp
void bind_vertex_buffer(const gl_vertex_buffer& buf, uint vattr);
~~~

Bind the buffer at a particular attribute location

#### Function unbind_vertex_buffer()

~~~ .cpp
void unbind_vertex_buffer(const gl_vertex_buffer& buf, uint vattr);
~~~

Unbind the buffer

#### Function unbind_vertex_buffer()

~~~ .cpp
void unbind_vertex_buffer(uint vattr);
~~~

Unbind the buffer

#### Function get_vertex_buffer_id()

~~~ .cpp
inline uint get_vertex_buffer_id(const gl_vertex_buffer& buf);
~~~

Get id

#### Function is_vertex_buffer_valid()

~~~ .cpp
inline bool is_vertex_buffer_valid(const gl_vertex_buffer& buf);
~~~

Check if defined

#### Function clear_vertex_buffer()

~~~ .cpp
void clear_vertex_buffer(gl_vertex_buffer& buf);
~~~

Destroys the buffer

#### Struct gl_element_buffer

~~~ .cpp
struct gl_element_buffer {
~~~

OpenGL vertex/element buffer

#### Function make_element_buffer()

~~~ .cpp
inline gl_element_buffer make_element_buffer(
    int num, int ncomp, const int* values, bool dynamic = false);
~~~

Creates a buffer.

#### Function make_element_buffer()

~~~ .cpp
inline gl_element_buffer make_element_buffer(
    const vector<int>& values, bool dynamic = false);
~~~

Creates a buffer.

#### Function make_element_buffer()

~~~ .cpp
inline gl_element_buffer make_element_buffer(
    const vector<vec2i>& values, bool dynamic = false);
~~~

Creates a buffer.

#### Function make_element_buffer()

~~~ .cpp
inline gl_element_buffer make_element_buffer(
    const vector<vec3i>& values, bool dynamic = false);
~~~

Creates a buffer.

#### Function make_element_buffer()

~~~ .cpp
inline gl_element_buffer make_element_buffer(
    const vector<vec4i>& values, bool dynamic = false);
~~~

Creates a buffer.

#### Function update_element_buffer()

~~~ .cpp
inline void update_element_buffer(
    gl_element_buffer& buf, int num, int ncomp, const int* values);
~~~

Updates the buffer with new data.

#### Function update_element_buffer()

~~~ .cpp
inline void update_element_buffer(
    gl_element_buffer& buf, const vector<int>& values);
~~~

Updates the buffer bid with new data.

#### Function update_element_buffer()

~~~ .cpp
inline void update_element_buffer(
    gl_element_buffer& buf, const vector<vec2i>& values);
~~~

Updates the buffer bid with new data.

#### Function update_element_buffer()

~~~ .cpp
inline void update_element_buffer(
    gl_element_buffer& buf, const vector<vec3i>& values);
~~~

Updates the buffer bid with new data.

#### Function update_element_buffer()

~~~ .cpp
inline void update_element_buffer(
    gl_element_buffer& buf, const vector<vec4i>& values);
~~~

Updates the buffer bid with new data.

#### Function draw_elems()

~~~ .cpp
void draw_elems(const gl_element_buffer& buf);
~~~

Draws elements.

#### Function get_element_buffer_id()

~~~ .cpp
inline uint get_element_buffer_id(const gl_element_buffer& buf);
~~~

Get id

#### Function is_element_buffer_valid()

~~~ .cpp
inline bool is_element_buffer_valid(const gl_element_buffer& buf);
~~~

Check if defined

#### Function clear_element_buffer()

~~~ .cpp
void clear_element_buffer(gl_element_buffer& buf);
~~~

Destroys the buffer

#### Struct gl_program

~~~ .cpp
struct gl_program {
~~~

OpenGL program

#### Function make_program()

~~~ .cpp
gl_program make_program(const string& vertex, const string& fragment);
~~~

Creates and OpenGL program from vertex and fragment code. Returns the
program id. Optionally return vertex and fragment shader ids. A VAO is
created.

#### Function clear_program()

~~~ .cpp
void clear_program(gl_program& prog);
~~~

Destroys the program pid and optionally the sahders vid and fid.

#### Function get_program_uniform_location()

~~~ .cpp
int get_program_uniform_location(const gl_program& prog, const string& name);
~~~

Get uniform location (simple GL wrapper that avoids GL includes)

#### Function get_program_attrib_location()

~~~ .cpp
int get_program_attrib_location(const gl_program& prog, const string& name);
~~~

Get uniform location (simple GL wrapper that avoids GL includes)

#### Function get_program_uniforms_names()

~~~ .cpp
vector<pair<string, int>> get_program_uniforms_names(const gl_program& prog);
~~~

Get the names of all uniforms

#### Function get_program_attributes_names()

~~~ .cpp
vector<pair<string, int>> get_program_attributes_names(const gl_program& prog);
~~~

Get the names of all attributes

#### Function set_program_uniform()

~~~ .cpp
bool set_program_uniform(
    gl_program& prog, int pos, const int* val, int ncomp, int count);
~~~

Set uniform integer values val for program pid and variable loc.
The values have nc number of components (1-4) and count elements
(for arrays).

#### Function set_program_uniform()

~~~ .cpp
bool set_program_uniform(
    gl_program& prog, int pos, const float* val, int ncomp, int count);
~~~

Set uniform float values val for program pid and variable var.
The values have nc number of components (1-4) and count elements
(for arrays).

#### Function set_program_uniform()

~~~ .cpp
inline bool set_program_uniform(gl_program& prog, int var, bool val);
~~~

Set uniform float values val for program pid and variable var.

#### Function set_program_uniform()

~~~ .cpp
inline bool set_program_uniform(gl_program& prog, int var, int val);
~~~

Set uniform float values val for program pid and variable var.

#### Function set_program_uniform()

~~~ .cpp
inline bool set_program_uniform(gl_program& prog, int var, float val);
~~~

Set uniform float values val for program pid and variable var.

#### Function set_program_uniform()

~~~ .cpp
inline bool set_program_uniform(gl_program& prog, int var, const vec2f& val);
~~~

Set uniform float values val for program pid and variable var.

#### Function set_program_uniform()

~~~ .cpp
inline bool set_program_uniform(gl_program& prog, int var, const vec3f& val);
~~~

Set uniform float values val for program pid and variable var.

#### Function set_program_uniform()

~~~ .cpp
inline bool set_program_uniform(gl_program& prog, int var, const vec4f& val);
~~~

Set uniform float values val for program pid and variable var.

#### Function set_program_uniform()

~~~ .cpp
inline bool set_program_uniform(gl_program& prog, int var, const vec2i& val);
~~~

Set uniform float values val for program pid and variable var.

#### Function set_program_uniform()

~~~ .cpp
inline bool set_program_uniform(gl_program& prog, int var, const vec3i& val);
~~~

Set uniform float values val for program pid and variable var.

#### Function set_program_uniform()

~~~ .cpp
inline bool set_program_uniform(gl_program& prog, int var, const vec4i& val);
~~~

Set uniform float values val for program pid and variable var.

#### Function set_program_uniform()

~~~ .cpp
inline bool set_program_uniform(gl_program& prog, int var, const mat4f& val);
~~~

Set uniform float values val for program pid and variable var.

#### Function set_program_uniform()

~~~ .cpp
inline bool set_program_uniform(gl_program& prog, int var, const frame3f& val);
~~~

Set uniform float values val for program pid and variable var.

#### Function set_program_uniform()

~~~ .cpp
inline bool set_program_uniform(
    gl_program& prog, int var, const int* val, int num);
~~~

Set uniform float values val for program pid and variable var.

#### Function set_program_uniform()

~~~ .cpp
inline bool set_program_uniform(
    gl_program& prog, int var, const float* val, int num);
~~~

Set uniform float values val for program pid and variable var.

#### Function set_program_uniform()

~~~ .cpp
inline bool set_program_uniform(
    gl_program& prog, int var, const vec2f* val, int num);
~~~

Set uniform float values val for program pid and variable var.

#### Function set_program_uniform()

~~~ .cpp
inline bool set_program_uniform(
    gl_program& prog, int var, const vec3f* val, int num);
~~~

Set uniform float values val for program pid and variable var.

#### Function set_program_uniform()

~~~ .cpp
inline bool set_program_uniform(
    gl_program& prog, int var, const vec4f* val, int num);
~~~

Set uniform float values val for program pid and variable var.

#### Function set_program_uniform()

~~~ .cpp
inline bool set_program_uniform(
    gl_program& prog, int var, const vec2i* val, int num);
~~~

Set uniform float values val for program pid and variable var.

#### Function set_program_uniform()

~~~ .cpp
inline bool set_program_uniform(
    gl_program& prog, int var, const vec3i* val, int num);
~~~

Set uniform float values val for program pid and variable var.

#### Function set_program_uniform()

~~~ .cpp
inline bool set_program_uniform(
    gl_program& prog, int var, const vec4i* val, int num);
~~~

Set uniform float values val for program pid and variable var.

#### Function set_program_uniform()

~~~ .cpp
inline bool set_program_uniform(
    gl_program& prog, int var, const mat4f* val, int num);
~~~

Set uniform float values val for program pid and variable var.

#### Function set_program_uniform()

~~~ .cpp
template <typename T>
inline bool set_program_uniform(
    gl_program& prog, const string& var, const T& val);
~~~

Set uniform float values val for program pid and variable var.

#### Function set_program_uniform()

~~~ .cpp
template <typename T>
inline bool set_program_uniform(
    gl_program& prog, const string& var, const T* val, int num);
~~~

Set uniform float values val for program pid and variable var.

#### Function set_program_uniform_texture()

~~~ .cpp
bool set_program_uniform_texture(
    gl_program& prog, int pos, const gl_texture_info& tinfo, uint tunit);
~~~

Set uniform texture id tid and unit tunit for program pid and variable
var.

#### Function set_program_uniform_texture()

~~~ .cpp
inline bool set_program_uniform_texture(gl_program& prog, int var, int varon,
    const gl_texture_info& tinfo, uint tunit);
~~~

Set uniform texture id tid and unit tunit for program pid and variable
var. Optionally sets the int variable varon to 0/1 whether the texture
is enable on not.

#### Function set_program_uniform_texture()

~~~ .cpp
inline bool set_program_uniform_texture(gl_program& prog, const string& var,
    const gl_texture_info& tinfo, uint tunit);
~~~

Set uniform texture id tid and unit tunit for program pid and variable
var.

#### Function set_program_uniform_texture()

~~~ .cpp
inline bool set_program_uniform_texture(gl_program& prog, const string& var,
    const string& varon, const gl_texture_info& tinfo, uint tunit);
~~~

Set uniform texture id tid and unit tunit for program pid and variable
var. Optionally sets the int variable varon to 0/1 whether the texture
is enable on not.

#### Function set_program_vertattr()

~~~ .cpp
bool set_program_vertattr(
    gl_program& prog, int pos, const float* value, int nc);
~~~

Sets a constant value for a vertex attribute for program pid and
variable var. The attribute has nc components.

#### Function set_program_vertattr()

~~~ .cpp
bool set_program_vertattr(gl_program& prog, int pos, const int* value, int nc);
~~~

Sets a constant value for a vertex attribute for program pid and
variable var. The attribute has nc components.

#### Function set_program_vertattr()

~~~ .cpp
bool set_program_vertattr(
    gl_program& prog, const string& var, const gl_vertex_buffer& buf);
~~~

Sets a vartex attribute for program pid and variable var to the buffer
bid. The attribute has nc components and per-vertex values values.

#### Function set_program_vertattr()

~~~ .cpp
bool set_program_vertattr(gl_program& prog, int pos,
    const gl_vertex_buffer& buf, int nc, const float* def);
~~~

Sets a vartex attribute for program pid and variable var. The attribute
has nc components and either buffer bid or a single value def
(if bid is zero). Convenience wrapper to above functions.

#### Function set_program_vertattr()

~~~ .cpp
inline bool set_program_vertattr(
    gl_program& prog, int var, const gl_vertex_buffer& buf, const vec2f& def);
~~~

Sets a vartex attribute for program pid and variable var. The attribute
is either a buffer bid or a single value def
(if bid is zero). Convenience wrapper to above functions.

#### Function set_program_vertattr()

~~~ .cpp
inline bool set_program_vertattr(
    gl_program& prog, int var, const gl_vertex_buffer& buf, const vec3f& def);
~~~

Sets a vartex attribute for program pid and variable var. The attribute
is either a buffer bid or a single value def
(if bid is zero). Convenience wrapper to above functions.

#### Function set_program_vertattr()

~~~ .cpp
inline bool set_program_vertattr(
    gl_program& prog, int var, const gl_vertex_buffer& buf, const vec4f& def);
~~~

Sets a vartex attribute for program pid and variable var. The attribute
is either a buffer bid or a single value def
(if bid is zero). Convenience wrapper to above functions.

#### Function set_program_vertattr()

~~~ .cpp
inline bool set_program_vertattr(gl_program& prog, const string& var,
    const gl_vertex_buffer& buf, const vec2f& def);
~~~

Sets a vartex attribute for program pid and variable var. The attribute
is either a buffer bid or a single value def
(if bid is zero). Convenience wrapper to above functions.

#### Function set_program_vertattr()

~~~ .cpp
inline bool set_program_vertattr(gl_program& prog, const string& var,
    const gl_vertex_buffer& buf, const vec3f& def);
~~~

Sets a vartex attribute for program pid and variable var. The attribute
is either a buffer bid or a single value def
(if bid is zero). Convenience wrapper to above functions.

#### Function set_program_vertattr()

~~~ .cpp
inline bool set_program_vertattr(gl_program& prog, const string& var,
    const gl_vertex_buffer& buf, const vec4f& def);
~~~

Sets a vartex attribute for program pid and variable var. The attribute
is either a buffer bid or a single value def
(if bid is zero). Convenience wrapper to above functions.

#### Function is_program_valid()

~~~ .cpp
inline bool is_program_valid(const gl_program& prog);
~~~

Check whether it is valid

#### Function bind_program()

~~~ .cpp
void bind_program(const gl_program& prog);
~~~

Binds a program

#### Function unbind_program()

~~~ .cpp
void unbind_program(const gl_program& prog);
~~~

Unbind a program

#### Struct gl_stdimage_program

~~~ .cpp
struct gl_stdimage_program {
~~~

A shader for displaying images

#### Function make_stdimage_program()

~~~ .cpp
gl_stdimage_program make_stdimage_program();
~~~

Initialize the program. Call with true only after the GL is initialized.

#### Function draw_image()

~~~ .cpp
inline void draw_image(gl_stdimage_program& prog, const gl_texture& txt,
    const vec2i& win_size, const vec2f& offset, float zoom, float exposure,
    float gamma, bool filmic);
~~~

As above but includes an exposure/gamma correction.

#### Function draw_image()

~~~ .cpp
inline void draw_image(gl_stdimage_program& prog, const gl_texture& txt,
    const vec2i& win_size, const vec2f& offset, float zoom);
~~~

Draw an texture tid of size img_w, img_h on a window of size win_w,
win_h with top-left corner at ox, oy with a zoom zoom.

#### Struct gl_stdsurface_program

~~~ .cpp
struct gl_stdsurface_program {
~~~

Shade with a physically-based standard shader based on Phong/GGX.
Filmic tone mapping from
https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/

#### Function make_stdsurface_program()

~~~ .cpp
gl_stdsurface_program make_stdsurface_program();
~~~

Initialize a standard shader. Call with true only after the gl has
been initialized

#### Function is_program_valid()

~~~ .cpp
inline bool is_program_valid(const gl_stdsurface_program& prog);
~~~

Check if the program is valid

#### Function begin_stdsurface_frame()

~~~ .cpp
inline void begin_stdsurface_frame(gl_stdsurface_program& prog,
    bool shade_eyelight, float tonemap_exposure, float tonemap_gamma,
    bool tonemap_filmic, const mat4f& camera_xform,
    const mat4f& camera_xform_inv, const mat4f& camera_proj);
~~~

Starts a frame by setting exposure/gamma values, camera transforms and
projection. Sets also whether to use full shading or a quick eyelight
preview.

#### Function end_stdsurface_frame()

~~~ .cpp
inline void end_stdsurface_frame(gl_stdsurface_program& prog);
~~~

Ends a frame.

#### Function set_stdsurface_lights()

~~~ .cpp
inline void set_stdsurface_lights(gl_stdsurface_program& prog, const vec3f& amb,
    int num, vec3f* pos, vec3f* ke, gl_ltype* type);
~~~

Set num lights with position pos, color ke, type ltype. Also set the
ambient illumination amb.

#### Function begin_stdsurface_shape()

~~~ .cpp
inline void begin_stdsurface_shape(
    gl_stdsurface_program& prog, const mat4f& xform);
~~~

Begins drawing a shape with transform xform.

#### Function end_stdsurface_shape()

~~~ .cpp
inline void end_stdsurface_shape(gl_stdsurface_program& prog);
~~~

End shade drawing.

#### Function set_stdsurface_highlight()

~~~ .cpp
inline void set_stdsurface_highlight(
    gl_stdsurface_program& prog, const vec4f& highlight);
~~~

Set the object as highlighted.

#### Function set_stdsurface_material()

~~~ .cpp
inline void set_stdsurface_material(gl_stdsurface_program& prog,
    material_type mtype, gl_etype etype, const vec3f& ke, const vec3f& kd,
    const vec3f& ks, float rs, float op, const gl_texture_info& ke_txt,
    const gl_texture_info& kd_txt, const gl_texture_info& ks_txt,
    const gl_texture_info& rs_txt, const gl_texture_info& norm_txt,
    const gl_texture_info& occ_txt, bool use_phong, bool double_sided,
    bool alpha_cutout);
~~~

Set material values with emission ke, diffuse kd, specular ks and
specular roughness rs, opacity op. Indicates textures ids with the
correspoinding XXX_txt variables. Sets also normal and occlusion
maps. Works for points/lines/triangles (diffuse for points,
Kajiya-Kay for lines, GGX/Phong for triangles).
Material type matches the scene material type.

#### Function set_stdsurface_vert()

~~~ .cpp
inline void set_stdsurface_vert(gl_stdsurface_program& prog,
    const gl_vertex_buffer& pos, const gl_vertex_buffer& norm,
    const gl_vertex_buffer& texcoord, const gl_vertex_buffer& color,
    const gl_vertex_buffer& tangsp);
~~~

Set vertex data with buffers for position pos, normals norm, texture
coordinates texcoord, per-vertex color color and tangent space tangsp.

#### Function set_stdsurface_vert_skinning()

~~~ .cpp
inline void set_stdsurface_vert_skinning(gl_stdsurface_program& prog,
    const gl_vertex_buffer& weights, const gl_vertex_buffer& joints,
    int nxforms, const mat4f* xforms);
~~~

Set vertex data with buffers for skinning.

#### Function set_stdsurface_vert_gltf_skinning()

~~~ .cpp
inline void set_stdsurface_vert_gltf_skinning(gl_stdsurface_program& prog,
    const gl_vertex_buffer& weights, const gl_vertex_buffer& joints,
    int nxforms, const mat4f* xforms);
~~~

Set vertex data with buffers for skinning.

#### Function set_stdsurface_vert_skinning_off()

~~~ .cpp
inline void set_stdsurface_vert_skinning_off(gl_stdsurface_program& prog);
~~~

Disables vertex skinning.

#### Struct gl_stdsurface_state

~~~ .cpp
struct gl_stdsurface_state {
~~~

State object for gl_stdsurface_program drawing. Members are not part of the
public API.

#### Struct gl_stdsurface_params

~~~ .cpp
struct gl_stdsurface_params {
    int camera_id = 0;
    int width = 360;
    int height = 360;
    float exposure = 0;
    float gamma = 2.2f;
    bool filmic = false;
    bool wireframe = false;
    bool edges = false;
    bool cutout = false;
    bool camera_lights = false;
    vec4f background = {0, 0, 0, 0};
    vec3f ambient = {0, 0, 0};
    void* hilighted = nullptr;
}
~~~

Params for  gl_stdsurface_program drawing

- Members:
    - camera_id:      camera id
    - width:      image width
    - height:      image height
    - exposure:      image exposure
    - gamma:      image gamma
    - filmic:      image filmic tonemapping
    - wireframe:      draw as wireframe
    - edges:      draw with overlaid edges
    - cutout:      draw with an alpha cutout for binary transparency
    - camera_lights:      camera light mode
    - background:      window background
    - ambient:      ambient illumination
    - hilighted:      highlighted object


#### Function make_stdsurface_state()

~~~ .cpp
gl_stdsurface_state* make_stdsurface_state();
~~~

Initialize gl_stdsurface_program draw state

#### Function update_stdsurface_state()

~~~ .cpp
void update_stdsurface_state(gl_stdsurface_state* st, const scene* scn,
    const gl_stdsurface_params& params);
~~~

Update gl_stdsurface_program draw state. This updates stdsurface meshes
and textures on the GPU.

#### Function clear_stdsurface_state()

~~~ .cpp
void clear_stdsurface_state(gl_stdsurface_state* st);
~~~

Clear gl_stdsurface_program draw state

#### Function draw_stdsurface_scene()

~~~ .cpp
void draw_stdsurface_scene(gl_stdsurface_state* st, const scene* scn,
    const gl_stdsurface_params& params);
~~~

Draw whole scene

#### Function void()

~~~ .cpp
typedef void (*gl_text_callback)(gl_window*, unsigned int);
~~~

Text callback

#### Function void()

~~~ .cpp
typedef void (*gl_mouse_callback)(gl_window*, int button, bool press, int mods);
~~~

Mouse callback

#### Function void()

~~~ .cpp
typedef void (*gl_refresh_callback)(gl_window*);
~~~

Window refresh callback

#### Struct gl_window

~~~ .cpp
struct gl_window {
~~~

Window

#### Function make_window()

~~~ .cpp
gl_window* make_window(
    int width, int height, const string& title, void* user_pointer = nullptr);
~~~

Initialize gl_window

#### Function set_window_callbacks()

~~~ .cpp
void set_window_callbacks(gl_window* win, gl_text_callback text_cb,
    gl_mouse_callback mouse_cb, gl_refresh_callback refresh_cb);
~~~

Set gl_window callbacks

#### Function clear_window()

~~~ .cpp
void clear_window(gl_window* win);
~~~

Clear gl_window

#### Function get_user_pointer()

~~~ .cpp
inline void* get_user_pointer(gl_window* win);
~~~

Gets the user poiner

#### Function set_window_title()

~~~ .cpp
void set_window_title(gl_window* win, const string& title);
~~~

Set gl_window title

#### Function wait_events()

~~~ .cpp
void wait_events(gl_window* win);
~~~

Wait events

#### Function poll_events()

~~~ .cpp
void poll_events(gl_window* win);
~~~

Poll events

#### Function swap_buffers()

~~~ .cpp
void swap_buffers(gl_window* win);
~~~

Swap buffers

#### Function should_close()

~~~ .cpp
bool should_close(gl_window* win);
~~~

Should close

#### Function get_mouse_button()

~~~ .cpp
int get_mouse_button(gl_window* win);
~~~

Mouse button

#### Function get_mouse_pos()

~~~ .cpp
vec2i get_mouse_pos(gl_window* win);
~~~

Mouse position

#### Function get_mouse_posf()

~~~ .cpp
vec2f get_mouse_posf(gl_window* win);
~~~

Mouse position

#### Function get_window_size()

~~~ .cpp
vec2i get_window_size(gl_window* win);
~~~

Window size

#### Function get_key()

~~~ .cpp
bool get_key(gl_window* win, int key);
~~~

Check if a key is pressed (not all keys are supported)

#### Function get_framebuffer_size()

~~~ .cpp
vec2i get_framebuffer_size(gl_window* win);
~~~

Framebuffer size

#### Function get_widget_size()

~~~ .cpp
inline int get_widget_size(gl_window* win);
~~~

Widgets

#### Function get_screenshot()

~~~ .cpp
vector<vec4b> get_screenshot(
    gl_window* win, vec2i& wh, bool flipy = true, bool back = false);
~~~

Read pixels

#### Function save_screenshot()

~~~ .cpp
inline void save_screenshot(gl_window* win, const string& imfilename);
~~~

Save a screenshot to disk

#### Function handle_camera_navigation()

~~~ .cpp
bool handle_camera_navigation(gl_window* win, camera* cam, bool navigation_fps);
~~~

Handle camera navigation.

#### Function init_widgets()

~~~ .cpp
void init_widgets(gl_window* win);
~~~

Initialize widgets

#### Function begin_widgets()

~~~ .cpp
bool begin_widgets(gl_window* win, const string& title);
~~~

Begin draw widget

#### Function end_widgets()

~~~ .cpp
void end_widgets(gl_window* win);
~~~

End draw widget

#### Function get_widget_active()

~~~ .cpp
bool get_widget_active(gl_window* win);
~~~

Whether widget are active

#### Function draw_separator_widget()

~~~ .cpp
void draw_separator_widget(gl_window* win);
~~~

Horizontal separator

#### Function draw_indent_widget_begin()

~~~ .cpp
void draw_indent_widget_begin(gl_window* win);
~~~

Indent widget

#### Function draw_indent_widget_end()

~~~ .cpp
void draw_indent_widget_end(gl_window* win);
~~~

Indent widget

#### Function draw_continue_widget()

~~~ .cpp
void draw_continue_widget(gl_window* win);
~~~

Continue line with next widget

#### Function draw_label_widget()

~~~ .cpp
void draw_label_widget(gl_window* win, const string& lbl, const string& msg);
~~~

Label widget

#### Function draw_label_widget()

~~~ .cpp
template <typename... Args>
inline void draw_label_widget(
    gl_window* win, const string& lbl, const string& fmt, const Args&... args);
~~~

Label widget

#### Function draw_label_widget()

~~~ .cpp
template <typename T>
inline void draw_label_widget(gl_window* win, const string& lbl, const T& val);
~~~

Label widget

#### Function draw_value_widget()

~~~ .cpp
bool draw_value_widget(gl_window* win, const string& lbl, string& str);
~~~

Value widget

#### Function draw_value_widget()

~~~ .cpp
bool draw_value_widget(gl_window* win, const string& lbl, int* val, int ncomp,
    int min = 0, int max = 1, int incr = 1);
~~~

Value widget

#### Function draw_value_widget()

~~~ .cpp
bool draw_value_widget(gl_window* win, const string& lbl, float* val, int ncomp,
    float min = 0, float max = 1, float incr = 1);
~~~

Value widget

#### Function draw_value_widget()

~~~ .cpp
inline bool draw_value_widget(gl_window* win, const string& lbl, int& val,
    int min = 0, int max = 1, int incr = 1);
~~~

Value widget

#### Function draw_value_widget()

~~~ .cpp
inline bool draw_value_widget(gl_window* win, const string& lbl, vec2i& val,
    int min = 0, int max = 1, int incr = 1);
~~~

Value widget

#### Function draw_value_widget()

~~~ .cpp
inline bool draw_value_widget(gl_window* win, const string& lbl, vec3i& val,
    int min = 0, int max = 1, int incr = 1);
~~~

Value widget

#### Function draw_value_widget()

~~~ .cpp
inline bool draw_value_widget(gl_window* win, const string& lbl, vec4i& val,
    int min = 0, int max = 1, int incr = 1);
~~~

Value widget

#### Function draw_value_widget()

~~~ .cpp
inline bool draw_value_widget(gl_window* win, const string& lbl, float& val,
    float min = 0, float max = 1, float incr = 1);
~~~

Value widget

#### Function draw_value_widget()

~~~ .cpp
inline bool draw_value_widget(gl_window* win, const string& lbl, vec2f& val,
    float min = 0, float max = 1, float incr = 1);
~~~

Value widget

#### Function draw_value_widget()

~~~ .cpp
inline bool draw_value_widget(gl_window* win, const string& lbl, vec3f& val,
    float min = 0, float max = 1, float incr = 1);
~~~

Value widget

#### Function draw_value_widget()

~~~ .cpp
inline bool draw_value_widget(gl_window* win, const string& lbl, vec4f& val,
    float min = 0, float max = 1, float incr = 1);
~~~

Value widget

#### Function draw_value_widget()

~~~ .cpp
inline bool draw_value_widget(gl_window* win, const string& lbl, mat4f& val,
    float min = 0, float max = 1, float incr = 1);
~~~

Slider widget

#### Function draw_value_widget()

~~~ .cpp
inline bool draw_value_widget(gl_window* win, const string& lbl, frame3f& val,
    float min = -1, float max = 1, float incr = 1);
~~~

Slider widget

#### Function draw_value_widget()

~~~ .cpp
inline bool draw_value_widget(
    gl_window* win, const string& lbl, quat4f& val, float incr = 1);
~~~

Slider widget

#### Function draw_color_widget()

~~~ .cpp
bool draw_color_widget(gl_window* win, const string& lbl, vec4f& val);
~~~

Color widget

#### Function draw_color_widget()

~~~ .cpp
bool draw_color_widget(gl_window* win, const string& lbl, vec4b& val);
~~~

Color widget

#### Function draw_color_widget()

~~~ .cpp
bool draw_color_widget(gl_window* win, const string& lbl, vec3f& val);
~~~

Color widget

#### Function draw_value_widget()

~~~ .cpp
bool draw_value_widget(gl_window* win, const string& lbl, int& val,
    const vector<pair<string, int>>& labels);
~~~

Enum widget

#### Function draw_value_widget()

~~~ .cpp
template <typename T>
inline bool draw_value_widget(gl_window* win, const string& lbl, T& val,
    const vector<pair<string, T>>& labels);
~~~

Enum widget

#### Function draw_value_widget()

~~~ .cpp
bool draw_value_widget(gl_window* win, const string& lbl, bool& val);
~~~

Bool widget

#### Function draw_button_widget()

~~~ .cpp
bool draw_button_widget(gl_window* win, const string& lbl);
~~~

Button widget

#### Function draw_header_widget()

~~~ .cpp
bool draw_header_widget(gl_window* win, const string& lbl);
~~~

Collapsible header

#### Function draw_tree_widget_begin()

~~~ .cpp
bool draw_tree_widget_begin(gl_window* win, const string& lbl);
~~~

Start tree node

#### Function draw_tree_widget_end()

~~~ .cpp
void draw_tree_widget_end(gl_window* win);
~~~

Collapsible header

#### Function draw_tree_widget_begin()

~~~ .cpp
bool draw_tree_widget_begin(
    gl_window* win, const string& lbl, void*& selection, void* content);
~~~

Start selectable tree node

#### Function draw_tree_widget_begin()

~~~ .cpp
bool draw_tree_widget_begin(gl_window* win, const string& lbl, void*& selection,
    void* content, const vec4f& col);
~~~

Start selectable tree node

#### Function draw_tree_widget_end()

~~~ .cpp
void draw_tree_widget_end(gl_window* win, void* content);
~~~

End selectable tree node

#### Function draw_tree_widget_leaf()

~~~ .cpp
void draw_tree_widget_leaf(
    gl_window* win, const string& lbl, void*& selection, void* content);
~~~

Selectable tree leaf node

#### Function draw_tree_widget_leaf()

~~~ .cpp
void draw_tree_widget_leaf(gl_window* win, const string& lbl, void*& selection,
    void* content, const vec4f& col);
~~~

Selectable tree leaf node

#### Function draw_image_widget()

~~~ .cpp
void draw_image_widget(
    gl_window* win, int tid, const vec2i& size, const vec2i& imsize);
~~~

Image widget

#### Function draw_scroll_widget_begin()

~~~ .cpp
void draw_scroll_widget_begin(
    gl_window* win, const string& lbl, int height, bool border);
~~~

Scroll region

#### Function draw_scroll_widget_end()

~~~ .cpp
void draw_scroll_widget_end(gl_window* win);
~~~

Scroll region

#### Function draw_scroll_widget_here()

~~~ .cpp
void draw_scroll_widget_here(gl_window* win);
~~~

Scroll region

#### Function draw_groupid_widget_begin()

~~~ .cpp
void draw_groupid_widget_begin(gl_window* win, int gid);
~~~

Group ids

#### Function draw_groupid_widget_begin()

~~~ .cpp
void draw_groupid_widget_begin(gl_window* win, void* gid);
~~~

Group ids

#### Function draw_groupid_widget_end()

~~~ .cpp
void draw_groupid_widget_end(gl_window* win);
~~~

Group ids

#### Function draw_tree_widget_color_begin()

~~~ .cpp
void draw_tree_widget_color_begin(gl_window* win, const vec4f& color);
~~~

Text color

#### Function draw_tree_widget_color_end()

~~~ .cpp
void draw_tree_widget_color_end(gl_window* win);
~~~

Text color

#### Function draw_tonemap_widgets()

~~~ .cpp
inline void draw_tonemap_widgets(gl_window* win, const string& lbl,
    float& exposure, float& gamma, bool& filmic);
~~~

Tonemapping widgets

#### Function draw_camera_widget()

~~~ .cpp
inline bool draw_camera_widget(
    gl_window* win, const string& lbl, scene* scn, int& cam_idx);
~~~

Draws a widget that can selected the camera

#### Function draw_scene_widgets()

~~~ .cpp
bool draw_scene_widgets(gl_window* win, const string& lbl, scene* scn,
    void*& selection, const unordered_map<texture*, gl_texture>& gl_txt);
~~~

Draws widgets for a whole scene. Used for quickly making demos.

