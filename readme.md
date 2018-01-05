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

