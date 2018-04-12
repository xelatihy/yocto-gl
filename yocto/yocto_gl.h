///
// # Yocto/GL: Tiny C++ Library for Physically-based Graphics
///
// Yocto/GL is a collection utilities for building physically-based graphics
// algorithms implemented as a two-file library (`yocto_gl.h`, `yocto_gl.cpp`),
// and released under the MIT license. Features include:
///
// - convenience math functions for graphics
// - static length vectors for 2, 3, 4 length of arbitrary type
// - static length matrices for 2x2, 3x3, 4x4 of arbitrary type
// - static length rigid transforms (frames), specialized for 2d and 3d space
// - linear algebra operations and transforms
// - axis aligned bounding boxes
// - rays and ray-primitive intersection
// - point-primitive distance and overlap tests
// - normal and tangent computation for meshes and lines
// - generation of tesselated meshes
// - mesh refinement with linear tesselation and Catmull-Cark subdivision
// - keyframed animation, skinning and morphing
// - random number generation via PCG32
// - simple image data structure and a few image operations
// - simple scene format
// - generation of image examples
// - generation of scene examples
// - procedural sun and sky HDR
// - procedural Perlin noise
// - BVH for intersection and closest point query
// - Python-like string and path operations
// - utilities to load and save entire text and binary files
// - immediate mode command line parser
// - simple logger
// - path tracer supporting surfaces and hairs, GGX and MIS
// - support for loading and saving Wavefront OBJ and Khronos glTF
// - OpenGL utilities to manage textures, buffers and prograrms
// - OpenGL shader for image viewing and GGX microfacet and hair rendering
///
// The current version is 0.4.0.
///
// ## Credits
///
// This library includes code from the PCG random number generator,
// boost hash_combine,
// code from "Real-Time Collision Detection" by Christer Ericson, base64
// encode/decode by René Nyffenegger and public domain code from
// github.com/sgorsten/linalg, gist.github.com/badboy/6267743 and
// github.com/nothings/stb_perlin.h.
///
// This library imports many symbols from std for three reasons: avoid
// verbosity , ensuring better conventions when calling math functions and
// allowing easy overriding of std containers if desired. Just do not
// flatten this namespace into yours if this is a concern.
///
// For most components of the library, the use should be relatively easy to
// understand if you are familiar with 3d computer graphics. For more complex
// components, we follow the usage below.
///
///
// ## Design Considerations
///
// Yocto/GL tries to follow a simple programming model inspired by C but with
// heavy use of operator overloading for math readability. We attempt to make
// the code easy to use use rather than as performant as possible.
// We adopt a functional style and only rarely use classes and methods.
// Using a function style makes the code easier to extend, more explicit in
// the function requirements, and easier to write parallel-friendly APIs.
// I guess you could call this "data-driven programming".
///
// The use of templates in Yocto was the reason for many refactorings, going
// from no template to heavy templates use. At this time, templates are used
// in the basic types to make the codebase shorter and reduce bugs,
// at the price of accessibility for beginners. The truth is that modern C++,
// a tenant of Yocto, is heavily templated anyway, so being able to read
// template code is necessary no matter how Yocto does things.
///
// We make use of exception for error reporting. This makes the code
// much cleaner and more in line with the expectation of most other programming
// languages.
///
// Finally, we often import symbols from the standard library rather than
// using the `std::name` pattern. We found that this improves consistency
// when using math functions, and is more readable with templates. We realize
// this is not standard, but the imports are hidden within the ygl namespace,
// so library users do not have to be concern about it.
///
///
// ## Compilation
///
// Yocto/GL is written in C++14 and compiles on OSX (clang from Xcode 9+),
// Linux (gcc 6+, clang 4+) and Windows (MSVC 2015, MSVC 2017).
///
// For image loading and saving, Yocto/GL depends on `stb_image.h`,
// `stb_image_write.h`, `stb_image_resize.h` and `tinyexr.h`. These features
// can be disabled by defining YGL_IMAGEIO to 0 before including this file.
// If these features are useful, then the implementation files need to
// included in the manner described by the respective libraries. To simplify
// builds, we provide a file that builds these libraries, `stb_image.cpp`.
///
// To support Khronos glTF, Yocto/GL depends on `json.hpp`. This feature can
// be disabled by defining YGL_GLTF to 0 before including this file.
///
// OpenGL utilities include the OpenGL libraries, use GLEW on Windows/Linux,
// GLFW for windows handling and Dear ImGui for UI support.
// Since OpenGL is quite onerous and hard to link, its support can be disabled
// by defining YGL_OPENGL to 1 before including this file. If you use any of
// the OpenGL calls, make sure to properly link to the OpenGL libraries on
// your system. For ImGUI, build with the libraries `imgui.cpp`,
// `imgui_draw.cpp`, `imgui_impl_glfw_gl3.cpp`.
///
///
// ## Example Applications
///
// You can see Yocto/GL in action in the following applications written to
// test the library:
///
// - `yview.cpp`: simple OpenGL viewer for OBJ and glTF scenes
// - `ytrace.cpp`: offline path-tracer
// - `yitrace.cpp.cpp`: interactive path-tracer
// - `yscnproc.cpp`: scene manipulation and conversion to/from OBJ and glTF
// - `ytestgen.cpp`: creates test cases for the path tracer and GL viewer
// - `yimview.cpp`: HDR/PNG/JPG image viewer with exposure/gamma tone mapping
// - `yimproc.cpp`: offline image manipulation.
///
// You can build the example applications using CMake with
//     `mkdir build; cd build; cmake ..; cmake --build`
///
// Here are two images rendered with the builtin path tracer, where the
// scenes are crated with the test generator.
///
// ![Yocto/GL](images/shapes.png)
///
// ![Yocto/GL](images/lines.png)
///
///
// ## Usage
///
// To use the library simply include this file and setup the compilation
// option as described above.
// All library features are documented at their definition and should be
// relatively easy to use if you are familiar with writing graphics code.
// You can find the extracted documentation at `yocto_gl.md`.
// Here we give an overview of some of the main features.
///
///
// ### Small Vectors and Matrices, Frames, Bounding Boxes and Transforms
///
// We provide common operations for small vectors and matrices typically used
// in graphics. In particular, we support 2-4 dimensional vectors of arbitrary
// `vec2f`, `vec3f`, `vec<T, 4>` with specializarion for float
// (`vec2f`, `vec3f`, `vec4f`), int (`vec2i`, `vec3i`, `vec4i`) and bytes
// (`vec4b`). Vector operations are templated so they work on every type, but
// many of them are well-defined only for float types.
///
// We support 2-4 dimensional generic matrices `mat2f`, `mat3f`,
// `mat4f`, with matrix-matrix and matrix-vector products, transposes and
// inverses. Matrices are stored in column-major ordered and are accessed and
// constructed by column.
///
// To represent transformations, most of the library facilities prefer the use
// coordinate frames, aka rigid transforms, represented as `frame3f`.
// The structure store three coordinate axis and the frame origin. This is
// equivalent to a rigid transform written as a column-major affine
// matrix. Transform operations are better behaved with this representation.
///
// We represent coordinate bounds with axis-aligned bounding boxes in 1-4
// dimensions: `bbox<T, 1>`, `bbox<T, 2>`, `bbox3f`, `bbox<T, 4>`. These
// types support expansion operation, union and containment. We provide
// operations to compute bounds for points, lines, triangles and quads.
///
// For all basic types we support iteration with `begin()`/`end()` pairs,
// data access with `data()`, `empty()` and `size()` and stream inout and
// output.
///
// For both matrices and frames we support transform operations for points,
// vectors and directions (`trasform_point()`, `trasform_vector()`,
// `trasform_direction()`). For frames we also the support inverse operations
// (`transform_xxx_inverse()`). Transform matrices and frames can be
// constructed from basic translation, rotation and scaling, e.g. with
// `translation_mat4f()` or `translation_frame3f()` respectively, etc. For
// rotation we support axis-angle and quaternions, with slerp.
///
///
// ### Random Number Generation, Noise, and Monte Carlo support
///
// This library supports many facilities helpful in writing sampling
// functions targeting path tracing and shape generations.
///
// 1. Random number generation with PCG32:
//     1. initialize the random number generator with `make_rng()`
//     2. advance the random number state with `advance_rng()`
//     3. if necessary, you can reseed the rng with `seed_rng()`
//     4. generate random integers in an interval with `next_rand1i()`
//     5. generate random floats and double in the [0,1) range with
//        `next_rand1f()`, `next_rand2f()`, `next_rand3f()`, `next_rand1d()`
//     6. you can skip random numbers with `advance_rng()` and get the skipped
//        length with `rng_distance()`
//     7. generate random shuffled sequences with `rng_shuffle()`
// 2. Perlin noise: `perlin_noise()` to generate Perlin noise with optional
//    wrapping, with fractal variations `perlin_ridge_noise()`,
//    `perlin_fbm_noise()`, `perlin_turbulence_noise()`
// 3. Monte Carlo support: warp functions from [0,1)^k domains to domains
//    commonly used in path tracing. In particular, use `sample_hemisphere()`,
//    `sample_sphere()`, `sample_hemisphere_cosine()`,
//    `sample_hemisphere_cospower()`. `sample_disk()`. `sample_cylinder()`.
//    `sample_triangle()`. For each warp, you can compute the PDF with
//    `sample_xxx_pdf()`.
///
///
// ### Shape Utilities
///
// The library contains a few function to help with typically geometry
// manipulation useful to support scene viewing and path tracing.
///
// 1. compute line tangents, and triangle and quad areas and normals
// 2. interpolate values over primitives with `eval_line()`,
//    `eval_triangle()` and `eval_quad()`
// 3. evaluate Bezier curves and derivatives with `eval_bezier()` and
//    `eval_bezier_derivative()`
// 4. compute smooth normals and tangents with `compute_normals()`
///   `compute_tangents()`
// 5. compute tangent frames from texture coordinates with
//    `compute_tangent_space()`
// 6. compute skinning with `compute_skinning()` and
//    `compute_matrix_skinning()`
// 6. create shapes with `make_cube()`, `make_sphere()`, `make_quad()`,
//    `make_fvcube()`, `make_hair()`, `make_suzanne()`, `make_lines()`,
//    `make_points()`, `make_sphere_cube()`, `make_cube_rounded()`,
//    `make_sphere_flipcap()`, `make_cylinder()`, `make_cylinder_rounded()`,
//    `make_disk()`, `make_cylinder_side()`, `make_disk_quad()`
// 7. merge element with `marge_lines()`, `marge_triangles()`, `marge_quads()`
// 8. shape sampling with `sample_points()`, `sample_lines()`,
//    `sample_triangles()`; initialize the sampling CDFs with
//    `sample_points_cdf()`, `sample_lines_cdf()`, `sample_triangles_cdf()`
// 9.  samnple a could of point over a surface with `sample_triangles_points()`
// 10. get edges and boundaries with `get_edges()`
// 11. convert quads to triangles with `convert_quads_to_triangles()`
// 12. convert face varying to vertex shared representations with
//     `convert_face_varying()`
// 13. subdivide elements by edge splits with `subdivide_lines()`,
//     `subdivide_triangles()`, `subdivide_quads()`, `subdivide_beziers()`
// 14. Catmull-Clark subdivision surface with `subdivide_catmullclark()`
///
///
// ### Animation utilities
///
// The library contains a few function to help with typical animation
// manipulation useful to support scene viewing.
///
// 1. evaluate keyframed values with step, linear and bezier interpolation with
//    `eval_keyframed_step()`, `eval_keyframed_linear()`,
//    `eval_keyframed_bezier()`
// 2. mesh skinning with `compute_matrix_skinning()`
///
///
// ### Image and color
///
// Images are stored with the `image` templated structure. The two most used
// image types are 4-byte per pixel sRGB images `image4b`, or 4-float per
// pixel HDR images `image4f`.
///
// 1. convert between byte and float images with `srgb_to_linear()` and
//    `linear_to_srgb()`
// 2. color conversion with `hsv_to_rgb()`, `xyz_to_rgb()` and `rgb_to_xyz()`
// 3. exposure-gamma tonemapping, with optional filmic curve, with
//    `tonemap_image()`
// 4. compositing support with `image_over()`
// 5. example image generation with `m,ake_grid_image()`,
//    `make_checker_image()`, `make_bumpdimple_image()`, `make_ramp_image()`,
//    `make_gammaramp_image()`, `make_gammaramp_imagef()`, `make_uv_image()`,
//    `make_uvgrid_image()`, `make_recuvgrid_image()`
// 6. bump to normal mapping with `bump_to_normal_map()`
// 7. HDR sun-sky with `m ake_sunsky_image()`
// 8. various noise images with `make_noise_image()`, `make_fbm_image()`,
//    `make_ridge_image()`, `make_turbulence_image()`
// 9. image loading and saving with `load_image4b()`, `load_image4f()`,
//    `save_image4b()`, `save_image4f()`
// 10. image resizing with `resize_image()`
///
///
// ### Ray Intersection and Point Overlap Queries
///
// We support ray-scene intersection for points, lines and triangles
// accelerated by a simple BVH data structure.  Our BVH is written for minimal
// code and not maximum speed, but still gives reasonable results. We suggest
// the use of Intel's Embree as a fast alternative.
///
// 1. use `ray3f` to represent rays
// 2. build the BVH with `build_points_bvh()`, `build_points_bvh()` or
//   `build_points_bvh()`
// 3. perform ray-element intersection with `intersect_points_bvh()`,
//   `intersect_lines_bvh()` and `intersect_triangles_bvh()`
// 4. perform point overlap queries with `overlap_points_bvh()`,
//   `overlap_lines_bvh()` and `overlap_triangles_bvh()`
// 5. to support custom elements, use `buid_bvh()`, `intersect_bvh()` and
//   `overlap_bvh()` and provide them with proper callbacks
// 6. we also experimentally support quads with the `xxx_quads_xxx()` functions
///
///
// ### Simple scene
///
// We support a simple scene model used to quickly write demos that lets you
// load/save Wavefront OBJ and Khronos glTF and perform several simple scene
// manipulation including ray-scene intersection and closest point queries.
///
// The geometry model is comprised of a set of shapes, which are indexed
// collections of points, lines, triangles and quads. Each shape may contain
// only one element type. Shapes are organized into a scene by creating shape
// instances, each its own transform. Materials are specified like in glTF and
// include emission, base-metallic and diffuse-specular parametrization,
// normal, occlusion and displacement mapping. Finally, the scene containers
// cameras and environment maps. Quad support in shapes is experimental and
// mostly supported for loading and saving.
///
// For low-level access to OBJ/glTF formats, you are best accessing the formats
// directly with Yocto/Obj and Yocto/glTF. This components provides a
// simplified high-level access to each format which is sufficient for most
// applications and tuned for quick creating viewers, renderers and simulators.
///
// 1. load a scene with `load_scene()` and save it with `save_scene()`.
// 2. add missing data with `add_normals()`, `add_names()`, `add_hierarchy()`
//    `add_tangent_space()`
// 3. use `compute_bbox()` to compute element bounds
// 4. can merge scene together with `merge_into()`
// 5. make example scenes with `make_test_scene()`
// 6. validate scene with `validate_scene()`
///
// Ray-intersection and closet-point routines supporting points,
// lines and triangles accelerated by a two-level bounding volume
// hierarchy (BVH). Quad support is experimental.
///
// 1. build the bvh with `make_bvh()`
// 2. perform ray-interseciton tests with `intersect_ray()`
//     - use early_exit=false if you want to know the closest hit point
//     - use early_exit=false if you only need to know whether there is a hit
//     - for points and lines, a radius is required
//     - for triangles, the radius is ignored
// 2. perform point overlap tests with `overlap_point()` to check whether
//    a point overlaps with an element within a maximum distance
//     - use early_exit as above
//     - for all primitives, a radius is used if defined, but should
//       be very small compared to the size of the primitive since the radius
//       overlap is approximate
// 3. perform instance overlap queries with `overlap_instance_bounds()`
// 4. use `refit_bvh()` to recompute the bvh bounds if transforms or vertices
//    are changed (you should rebuild the bvh for large changes)
///
// Notes: Quads are internally handled as a pair of two triangles v0,v1,v3 and
// v2,v3,v1, with the u/v coordinates of the second triangle corrected as 1-u
// and 1-v to produce a quad parametrization where u and v go from 0 to 1. This
// is equivalent to Intel's Embree.
///
///
// ### Pathtracing
///
// We supply a path tracer implementation with support for textured mesh
// lights, GGX/Phong materials, environment mapping. The interface supports
// progressive parallel execution. The path tracer takes as input a scene
// and update pixels in image with traced samples. We use a straightfoward
// path tracer with MIS and also a few simpler shaders for debugging or
// quick image generation.
///
// Materials are represented as sums of an emission term, a diffuse term and
// a specular microfacet term (GGX or Phong). Only opaque for now. We pick
// a proper material type for each shape element type (points, lines,
// triangles).
///
// Lights are defined as any shape with a material emission term. Additionally
// one can also add environment maps. But even if you can, you might want to
// add a large triangle mesh with inward normals instead. The latter is more
// general (you can even more an arbitrary shape sun). For now only the first
// env is used.
///
// 1. build the ray-tracing acceleration structure with `make_bvh()`
// 2. prepare lights for rendering `update_lights()`
// 3. define rendering params with the `trace_params` structure
// 4. render blocks of samples with `trace_block()`
///
// The code can also run in fully asynchronous mode to preview images in a
// window.
///
// 1. build the ray-tracing acceleration structure with `make_bvh()`
// 2. prepare lights for rendering `update_lights()`
// 3. define rendering params with the `trace_params` structure
// 4. initialize the progressive rendering buffers
// 5. start the progressive renderer with `trace_async_start()`
// 7. stop the progressive renderer with `trace_async_stop()`
///
///
// ### Wavefront OBJ
///
// Wavefront OBJ/MTL loader and writer with support for points,
// lines, triangles and general polygons and all materials properties.
// Contains also a few extensions to easily create demos such as per-vertex
// color and radius, cameras, environment maps and instances.
// Can use either a low-level OBJ representation, from this files,
// or a high level flattened representation included in Yocto/Scn.
///
// Both in reading and writing, OBJ has no clear convention on the orientation
// of textures Y axis. So in many cases textures appears flipped. To handle
// that, use the option to flip textures coordinates on either saving or
// loading. By default texture coordinates are flipped since this seems
// the convention found on test cases collected on the web. The value Tr
// has similar problems, since its relation to opacity is software specific.
// Again we let the user chose the conversion and set the default to the
// one found on the web.
///
// In the high level interface, shapes are indexed meshes and are described
// by arrays of vertex indices for points/lines/triangles and arrays for vertex
// positions, normals, texcoords, color and radius. The latter two as
// extensions. Since OBJ is a complex formats that does not match well with
// current GPU rendering / path tracing algorithms, we adopt a simplification
// similar to other single file libraries:
// 1. vertex indices are unique, as in OpenGL and al standard indexed triangle
//   meshes data structures, and not OBJ triplets; YOCTO_OBJ ensures that no
//   vertex duplication happens thought for same triplets
// 2. we split shapes on changes to groups and materials, instead of keeping
//   per-face group/material data; this makes the data usable right away in
//   a GPU viewer; this is not a major limitation if we accept the previous
//   point that already changes shapes topology.
///
// 1. load a obj data with `load_obj()`; can load also textues
// 2. look at the `obj_XXX` data structures for access to individual elements
// 3. use obj back to disk with `save_obj()`; can also save textures
// 4. use get_shape() to get a flattened shape version that contains only
//    triangles, lines or points
///
///
// ### Khronos glTF
///
// Khronos GLTF loader and writer for Khronos glTF format. Supports
// all the glTF spec and the Khronos extensions. All parsing and writing code
// is autogenerated form the schema. Supports glTF version 2.0 and the
// following extensions: `KHR_binary_glTF` and `KHR_specular_glossiness`.
///
// This component depends on `json.hpp` and, for image loading and saving,
// it depends on `stb_image.h`, `stb_image_write.h`, `stb_image_resize.h` and
// `tinyexr.h`. This feature can be disabled as before.
///
// The library provides a low  level interface that is a direct
// C++ translation of the glTF schemas and should be used if one wants
// complete control over the format or an application wants to have their
// own scene code added. A higher-level interface is provided by the scene
// or by `yocto_gltf.h`.
///
// glTF is a very complex file format and was designed mainly with untyped
// languages in mind. We attempt to match the glTF low-level interface
// to C++ as best as it can. Since the code is generated from the schema, we
// follow glTF naming conventions and typing quite well. To simplify adoption
// and keep the API relatively simple we use vector as arrays and use
// pointers to reference to all glTF objects. While this makes it less
// efficient than it might have been, glTF heavy use of optional values makes
// this necessary. At the same time, we do not keep track of set/unset values
// for basic types (int, float, bool) as a compromise for efficiency.
///
// glTF uses integer indices to access objects.
// While writing code ourselves we found that we add significant problems
// since we would use an index to access the wrong type of scene objects.
// For this reasons, we use an explicit index `glTFid<T>` that can only access
// an object of type T. Internally this is just the same old glTF index. But
// this can used to access the scene data with `glTF::get<T>(index)`.
///
// 1. load a glTF model with `load_gltf()`
// 2. look at the `glTFXXX` data structures for access to individual elements
// 3. save glTF back to disk with `save_gltf()`
///
///
// ### OpenGL support
///
// We include a set of utilities to draw on screen with OpenGL 3.3, manage
// windows with GLFW and draw immediate-mode widgets with ImGui.
///
// 1. texture and buffer objects with `gltexture` and `gl_buffer`
//     - create textures/buffers with appropriate constructors
//     - check validity with `is_valid()`
//     - update textures/buffers with `update()` functions
//     - delete textures/buffers with `clear()`
//     - bind/unbind textures/buffers with `bind()`/`unbind()`
//     - draw elements with `gl_buffer::draw_elems()`
// 2. program objects with `glprogram`
//     - program creation with constructor
//     - check validity with `is_valid()`
//     - delete with `clear()`
//     - uniforms with `set_gluniform()`
//     - vertex attrib with `set_glattribute()`
//     - draw elements with `gl_buffer::draw_elems()`
// 3. image viewing with `glimage_program`, with support for tone mapping.
// 4. draw surfaces and hair with GGX/Kayjia-Kay with `glsurface_program`
//     - initialize the program with constructor
//     - check validity with `is_valid()`
//     - start/end each frame with `begin_frame()`, `end_frame()`
//     - define lights with `set_lights()`
//     - start/end each shape with `begin_shape()`, `end_shape()`
//     - define material Parameters with `set_material()`
//     - define vertices with `set_vert()`
//     - draw elements with `draw_elems()`
// 5. draw yocto scenes using the above shader
//     - initialize the rendering state with `init_stdsurface_state()`
//     - load/update meshes and textures with `update_stdsurface_state()`
//     - setup draw params using a `glsurface_params` struct
//     - draw scene with `draw_glsurface_scene()`
// 6. also includes other utlities for quick OpenGL hacking
// 7. GLFW window with `glwindow`
//     - create with constructor
//     - delete with `clear()`
//     - set callbacks with `set_callbacks()`
//     - includes carious utilities to query window, mouse and keyboard
// 8. immediate mode widgets using ImGui
//     - init with `init_widget()`
//     - use the various widget calls to draw the widget and handle events
///
///
// ### Other Utilities
///
// We include additional utilities for writing command line applications and
// manipulating files.
///
// 1. Python-like string operations: `startswith()`, `endswith()`,
// `contains()`,
//    `splitlines()`, `partition()`, `split()`, `splitlines()`, `strip()`,
//    `rstrip()`, `lstrip()`, `join()`, `lower()`, `upper()`, `isspace()`,
//    `replace()`
// 2. Path-like path operations: `path_dirname()`, `path_extension()`,
//    `path_basename()`, `path_filename()`, `replace_path_extension()`,
//    `prepend_path_extension()`, `split_path()`
// 3. Python-like format strings (only support for position arguments and no
//    formatting commands): `format()`, `print()`
// 5. load/save entire files: `load_binary()`, `load_text()`,
//    `save_text()` and `save_binary()`
// 4. simple logger with support for console and file streams:
//     1. create a `logger`
//     2. add more streams with `addconsole_stream()` or `add_file_stream()`
//     3. write log messages with `log_msg()` and its variants
//     4. you can also use a global default logger with the free functions
//        `log_XXX()`
// 5. timer for simple access to `std::chrono`:
//     1. create a `timer`
//     2. start and stop the clock with `start()` and `stop()`
//     3. get time with `elapsed_time()`
///
///
// ### Command Line Parsing
///
// The library includes a simple command line parser that parses commands in
// immediate mode, i.e. when an option is declared. The parser supports options
// and unnamed arguments with generic types parsed using C++ stream. The
// parser autogenerates its own documentation. This allows to write complex
// command lines with a tiny amount of implementation code on both the library
// and user end.
///
// 1. create a `cmdline` parser object by passing `argc, argv, name, help`
//     - an option for printing help is automatically added
// 2. for each option, parse it calling the functions `parse_opt()`
//     - options are parsed on the fly and a comprehensive help is
//       automatically generated
//     - supports bool (flags), int, float, double, string, enums
//     - options names are "--longname" for longname and "-s" for short
//     - command line format is "--longname value", "-s v" for all but flags
//     - values are parsed with `iostream <<` operators
//     - for general use `opt = parse_opt<type>()`
//     - for boolean flags is `parse_flag()`
//     - for enums use `parse_opte()`
// 3. for each unnamed argument, parse it calling the functions parse_arg()
//     - names are only used for help
//     - supports types as above
//     - for general use `arg = parse_arg<type>()`
//     - to parse all remaining values use `args = parse_arga<type>(...)`
// 4. end cmdline parsing with `check_parsing()` to check for unused values,
//    missing arguments
// 5. to check for error use `should_exit()` and to print the message use
//    `get_message()`
// 6. since arguments are parsed immediately, one can easily implement
//    subcommands by just branching the command line code based on a read
//    argument without any need for complex syntax
///
namespace ygl {}

//
// LICENSE:
//
// Copyright (c) 2016 -- 2017 Fabio Pellacini
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
//
// LICENSE OF INCLUDED SOFTWARE for Pcg random number generator
//
// This code also includes a small exerpt from http://www.pcg-random.org/
// licensed as follows
// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)
//
//
// LICENSE OF INCLUDED CODE FOR BASE64 (base64.h, base64.cpp)
//
// Copyright (C) 2004-2008 René Nyffenegger
//
// This source code is provided 'as-is', without any express or implied
// warranty. In no event will the author be held liable for any damages
// arising from the use of this software.
//
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
//
// 1. The origin of this source code must not be misrepresented; you must not
// claim that you wrote the original source code. If you use this source code
// in a product, an acknowledgment in the product documentation would be
// appreciated but is not required.
//
// 2. Altered source versions must be plainly marked as such, and must not be
// misrepresented as being the original source code.
//
// 3. This notice may not be removed or altered from any source distribution.
//
// René Nyffenegger rene.nyffenegger@adp-gmbh.ch
//
//
//
// LICENSE OF INCLUDED CODE FOR PERLIN NOISE (stb_perlin.h)
//
// Copyright (c) 2017 Sean Barrett
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions: The above copyright
// notice and this permission notice shall be included in all copies or
// substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS",
// WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
// TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
// FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
// TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
// THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//

#ifndef _YGL_H_
#define _YGL_H_

// -----------------------------------------------------------------------------
// COMPILATION OPTIONS
// -----------------------------------------------------------------------------

// enable image io
#ifndef YGL_IMAGEIO
#define YGL_IMAGEIO 1
#endif

// enable glTF
#ifndef YGL_GLTF
#define YGL_GLTF 1
#endif

// enable OpenGL
#ifndef YGL_OPENGL
#define YGL_OPENGL 1
#endif

// enable fast OBJ parsing
#ifndef YGL_FASTOBJ
#define YGL_FASTOBJ 1
#endif

// enable explicit json objects in glTF
#ifndef YGL_GLTFJSON
#define YGL_GLTFJSON 0
#endif

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <map>
#include <string>
#include <thread>
#include <tuple>
#include <typeinfo>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// -----------------------------------------------------------------------------
// MATH CONSTANTS AND FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

using std::acos;
using std::asin;
using std::atan;
using std::atan2;
using std::cos;
using std::exp;
using std::floor;
using std::isfinite;
using std::log;
using std::pow;
using std::round;
using std::sin;
using std::sqrt;
using std::tan;
using namespace std::string_literals;
using namespace std::literals;

using byte = unsigned char;
using uint = unsigned int;

// const auto pi = 3.14159265f;
const auto pi = (float)M_PI;
const auto flt_max = FLT_MAX;
const auto flt_min = -FLT_MAX;
const auto flt_eps = FLT_EPSILON;

inline int abs(int x) { return (x < 0) ? -x : x; }
inline float abs(float x) { return (x < 0) ? -x : x; }
inline int min(int x, int y) { return (x < y) ? x : y; }
inline float min(float x, float y) { return (x < y) ? x : y; }
inline int max(int x, int y) { return (x > y) ? x : y; }
inline float max(float x, float y) { return (x > y) ? x : y; }
inline int clamp(int x, int min_, int max_) { return min(max(x, min_), max_); }
inline float clamp(float x, float min_, float max_) {
    return min(max(x, min_), max_);
}
inline float lerp(float a, float b, float u) { return a * (1 - u) + b * u; }

}  // namespace ygl

// -----------------------------------------------------------------------------
// VECTORS
// -----------------------------------------------------------------------------
namespace ygl {

// Small size vectors.
struct vec2f {
    float x = 0, y = 0;
};
struct vec3f {
    float x = 0, y = 0, z = 0;
};
struct vec4f {
    float x = 0, y = 0, z = 0, w = 0;
};

// Integer small-size vectors.
struct vec2i {
    int x = 0, y = 0;
};
struct vec3i {
    int x = 0, y = 0, z = 0;
};
struct vec4i {
    int x = 0, y = 0, z = 0, w = 0;
};

// Byte small-sized vector for color.
struct vec4b {
    byte x = 0, y = 0, z = 0, w = 0;
};

// Zero vector constants.
const auto zero2f = vec2f();
const auto zero3f = vec3f();
const auto zero4f = vec4f();
const auto zero2i = vec2i();
const auto zero3i = vec3i();
const auto zero4i = vec4i();
const auto zero4b = vec4b();

// Vector comparison operations.
inline bool operator==(const vec2f& a, const vec2f& b) {
    return a.x == b.x && a.y == b.y;
}
inline bool operator!=(const vec2f& a, const vec2f& b) {
    return a.x != b.x || a.y != b.y;
}
inline bool operator==(const vec2i& a, const vec2i& b) {
    return a.x == b.x && a.y == b.y;
}
inline bool operator!=(const vec2i& a, const vec2i& b) {
    return a.x != b.x || a.y != b.y;
}

inline bool operator==(const vec3f& a, const vec3f& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}
inline bool operator!=(const vec3f& a, const vec3f& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z;
}
inline bool operator==(const vec3i& a, const vec3i& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}
inline bool operator!=(const vec3i& a, const vec3i& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z;
}

inline bool operator==(const vec4f& a, const vec4f& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}
inline bool operator!=(const vec4f& a, const vec4f& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w;
}
inline bool operator==(const vec4i a, const vec4i& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}
inline bool operator!=(const vec4i& a, const vec4i& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w;
}

// Vector operations.
inline vec2f operator-(const vec2f& a) { return {-a.x, -a.y}; }
inline vec2f operator+(const vec2f& a, const vec2f& b) {
    return {a.x + b.x, a.y + b.y};
}
inline vec2f operator-(const vec2f& a, const vec2f& b) {
    return {a.x - b.x, a.y - b.y};
}
inline vec2f operator*(const vec2f& a, const vec2f& b) {
    return {a.x * b.x, a.y * b.y};
}
inline vec2f operator*(const vec2f& a, float b) { return {a.x * b, a.y * b}; }
inline vec2f operator*(float a, const vec2f& b) { return {a * b.x, a * b.y}; }
inline vec2f operator/(const vec2f& a, const vec2f& b) {
    return {a.x / b.x, a.y / b.y};
}
inline vec2f operator/(const vec2f& a, float b) { return {a.x / b, a.y / b}; }
inline vec2f operator/(float a, const vec2f& b) { return {a / b.x, a / b.y}; }

// Vector operations.
inline vec3f operator+(const vec3f& a) { return a; }
inline vec3f operator-(const vec3f& a) { return {-a.x, -a.y, -a.z}; }
inline vec3f operator+(const vec3f& a, const vec3f& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
inline vec3f operator-(const vec3f& a, const vec3f& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
inline vec3f operator*(const vec3f& a, const vec3f& b) {
    return {a.x * b.x, a.y * b.y, a.z * b.z};
}
inline vec3f operator*(const vec3f& a, float b) {
    return {a.x * b, a.y * b, a.z * b};
}
inline vec3f operator*(float a, const vec3f& b) {
    return {a * b.x, a * b.y, a * b.z};
}
inline vec3f operator/(const vec3f& a, const vec3f& b) {
    return {a.x / b.x, a.y / b.y, a.z / b.z};
}
inline vec3f operator/(const vec3f& a, float b) {
    return {a.x / b, a.y / b, a.z / b};
}
inline vec3f operator/(float a, const vec3f& b) {
    return {a / b.x, a / b.y, a / b.z};
}

// Vector operations.
inline vec4f operator-(const vec4f& a) { return {-a.x, -a.y, -a.z, -a.w}; }
inline vec4f operator+(const vec4f& a, const vec4f& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}
inline vec4f operator-(const vec4f& a, const vec4f& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w};
}
inline vec4f operator*(const vec4f& a, const vec4f& b) {
    return {a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w};
}
inline vec4f operator*(const vec4f& a, float b) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
}
inline vec4f operator*(float a, const vec4f& b) {
    return {a * b.x, a * b.y, a * b.z, a * b.w};
}
inline vec4f operator/(const vec4f& a, const vec4f& b) {
    return {a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w};
}
inline vec4f operator/(const vec4f& a, float b) {
    return {a.x / b, a.y / b, a.z / b, a.w / b};
}
inline vec4f operator/(float a, const vec4f& b) {
    return {a / b.x, a / b.y, a / b.z, a / b.w};
}

// Vector assignments
inline vec2f& operator+=(vec2f& a, const vec2f& b) { return a = a + b; }
inline vec2f& operator-=(vec2f& a, const vec2f& b) { return a = a - b; }
inline vec2f& operator*=(vec2f& a, const vec2f& b) { return a = a * b; }
inline vec2f& operator*=(vec2f& a, float b) { return a = a * b; }
inline vec2f& operator/=(vec2f& a, const vec2f& b) { return a = a / b; }
inline vec2f& operator/=(vec2f& a, float b) { return a = a / b; }

inline vec3f& operator+=(vec3f& a, const vec3f& b) { return a = a + b; }
inline vec3f& operator-=(vec3f& a, const vec3f& b) { return a = a - b; }
inline vec3f& operator*=(vec3f& a, const vec3f& b) { return a = a * b; }
inline vec3f& operator*=(vec3f& a, float b) { return a = a * b; }
inline vec3f& operator/=(vec3f& a, const vec3f& b) { return a = a / b; }
inline vec3f& operator/=(vec3f& a, float b) { return a = a / b; }

inline vec4f& operator+=(vec4f& a, const vec4f& b) { return a = a + b; }
inline vec4f& operator-=(vec4f& a, const vec4f& b) { return a = a - b; }
inline vec4f& operator*=(vec4f& a, const vec4f& b) { return a = a * b; }
inline vec4f& operator*=(vec4f& a, float b) { return a = a * b; }
inline vec4f& operator/=(vec4f& a, const vec4f& b) { return a = a / b; }
inline vec4f& operator/=(vec4f& a, float b) { return a = a / b; }

// Vector products and lengths.
inline float dot(const vec2f& a, const vec2f& b) {
    return a.x * b.x + a.y * b.y;
}
inline float dot(const vec3f& a, const vec3f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline float dot(const vec4f& a, const vec4f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}
inline float cross(const vec2f& a, const vec2f& b) {
    return a.x * b.y - a.y * b.x;
}
inline vec3f cross(const vec3f& a, const vec3f& b) {
    return {
        a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline float length(const vec2f& a) { return sqrt(dot(a, a)); }
inline float length(const vec3f& a) { return sqrt(dot(a, a)); }
inline float length(const vec4f& a) { return sqrt(dot(a, a)); }
inline float length_sqr(const vec2f& a) { return dot(a, a); }
inline float length_sqr(const vec3f& a) { return dot(a, a); }
inline float length_sqr(const vec4f& a) { return dot(a, a); }
inline vec2f normalize(const vec2f& a) { return length(a) ? a / length(a) : a; }
inline vec3f normalize(const vec3f& a) { return length(a) ? a / length(a) : a; }
inline vec4f normalize(const vec4f& a) { return length(a) ? a / length(a) : a; }

// Vecror angles and slerps.
inline float angle(const vec3f& a, const vec3f& b) {
    return acos(clamp(dot(normalize(a), normalize(b)), -1.0f, 1.0f));
}
inline vec4f slerp(const vec4f& a, const vec4f& b, float u) {
    // https://en.wikipedia.org/wiki/Slerp
    auto an = normalize(a), bn = normalize(b);
    auto d = dot(an, bn);
    if (d < 0) {
        bn = -bn;
        d = -d;
    }
    if (d > 0.9995f) return normalize(an + u * (bn - an));
    auto th = acos(clamp(d, -1.0f, 1.0f));
    if (!th) return an;
    return an * (sin(th * (1 - u)) / sin(th)) + bn * (sin(th * u) / sin(th));
}

// Orthogonal vectors.
inline vec3f orthogonal(const vec3f& v) {
    // http://lolengine.net/blog/2013/09/21/picking-orthogonal-vector-combing-coconuts)
    return fabs(v.x) > fabs(v.z) ? vec3f{-v.y, v.x, 0} : vec3f{0, -v.z, v.y};
}
inline vec3f orthonormalize(const vec3f& a, const vec3f& b) {
    return normalize(a - b * dot(a, b));
}

// Reflected and refracted vector.
inline vec3f reflect(const vec3f& w, const vec3f& n) {
    return -w + 2 * dot(n, w) * n;
}
inline vec3f refract(const vec3f& w, const vec3f& n, float eta) {
    // auto k = 1.0 - eta * eta * (1.0 - dot(n, w) * dot(n, w));
    auto k = 1 - eta * eta * max(0.0f, 1 - dot(n, w) * dot(n, w));
    if (k < 0) return zero3f;  // tir
    return -w * eta + (eta * dot(n, w) - sqrt(k)) * n;
}

// Max element and clamp.
inline vec2f clamp(const vec2f& x, float min, float max) {
    return {clamp(x.x, min, max), clamp(x.y, min, max)};
}
inline vec3f clamp(const vec3f& x, float min, float max) {
    return {clamp(x.x, min, max), clamp(x.y, min, max), clamp(x.z, min, max)};
}
inline vec4f clamp(const vec4f& x, float min, float max) {
    return {clamp(x.x, min, max), clamp(x.y, min, max), clamp(x.z, min, max),
        clamp(x.w, min, max)};
}
inline float max(const vec2f& a) { return max(a.x, a.y); }
inline float max(const vec3f& a) { return max(max(a.x, a.y), a.z); }
inline float max(const vec4f& a) { return max(max(max(a.x, a.y), a.z), a.w); }
inline float min(const vec2f& a) { return min(a.x, a.y); }
inline float min(const vec3f& a) { return min(min(a.x, a.y), a.z); }
inline float min(const vec4f& a) { return min(min(min(a.x, a.y), a.z), a.w); }

// Quaternion operatons represented as xi + yj + zk + w
const auto identity_quat4f = vec4f{0, 0, 0, 1};
inline vec4f quat_mul(const vec4f& a, float b) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
}
inline vec4f quat_mul(const vec4f& a, const vec4f& b) {
    return {a.x * b.w + a.w * b.x + a.y * b.w - a.z * b.y,
        a.y * b.w + a.w * b.y + a.z * b.x - a.x * b.z,
        a.z * b.w + a.w * b.z + a.x * b.y - a.y * b.x,
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z};
}
inline vec4f quat_conjugate(const vec4f& a) { return {-a.x, -a.y, -a.z, a.w}; }
inline vec4f quat_inverse(const vec4f& a) {
    return quat_conjugate(a) / length_sqr(a);
}

// Element-wise float to byte conversion.
inline vec4b float_to_byte(const vec4f& a) {
    return {(byte)max(0, min(int(a.x * 256), 255)),
        (byte)max(0, min(int(a.y * 256), 255)),
        (byte)max(0, min(int(a.z * 256), 255)),
        (byte)max(0, min(int(a.w * 256), 255))};
}
inline vec4f byte_to_float(const vec4b& a) {
    return {a.x / 255.0f, a.y / 255.0f, a.z / 255.0f, a.w / 255.0f};
}

// Coordinate conversions. Spherical (phi,theta,r). Cylindrical (phi,r,z).
inline vec3f cartesian_to_spherical(const vec3f& p) {
    auto r = length(p);
    if (!r) return {0, 0, 0};
    auto r2 = length(vec2f{p.x, p.y});
    return {atan2(p.y / r2, p.x / r2), acos(clamp(p.z / r, -1.0f, 1.0f)), r};
}
inline vec3f spherical_to_cartesian(const vec3f& s) {
    return {
        cos(s.x) * sin(s.y) * s.z, sin(s.x) * sin(s.y) * s.z, cos(s.y) * s.z};
}
inline vec3f cartesian_to_sphericaly(const vec3f& p) {
    auto r = length(p);
    if (!r) return {0, 0, 0};
    auto r2 = length(vec2f{p.x, p.z});
    return {atan2(p.z / r2, p.x / r2), acos(clamp(p.y / r, -1.0f, 1.0f)), r};
}
inline vec3f sphericaly_to_cartesian(const vec3f& s) {
    return {
        cos(s.x) * sin(s.y) * s.z, cos(s.y) * s.z, sin(s.x) * sin(s.y) * s.z};
}
inline vec3f cartesian_to_cylindrical(const vec3f& p) {
    auto r = length(vec2f{p.x, p.y});
    if (!r) return {0, 0, p.z};
    return {atan2(p.y / r, p.x / r), r, p.z};
}
inline vec3f cylindrical_to_cartesian(const vec3f& c) {
    return {cos(c.x) * c.y, sin(c.x) * c.y, c.z};
}

// Cartesian to elliptical mapping
// Cartesian (x,y) in [-1,1]^2. Elliptical |(u,v)|<=1.
inline vec2f cartesian_to_elliptical(const vec2f& xy) {
    // Analytical Methods for Squaring the Disc, by C. Fong
    // https://arxiv.org/abs/1509.06344
    auto x = xy.x, y = xy.y;
    auto u = x * sqrt(1 - y * y / 2);
    auto v = y * sqrt(1 - x * x / 2);
    return {u, v};
}
inline vec2f elliptical_to_cartesian(const vec2f& uv) {
    // Analytical Methods for Squaring the Disc, by C. Fong
    // https://arxiv.org/abs/1509.06344
    auto u = uv.x, v = uv.y;
    auto x = sqrt(2 + u * u - v * v + 2 * sqrt(2.0f) * u) / 2.0f -
             sqrt(2 + u * u - v * v - 2 * sqrt(2.0f) * u) / 2.0f;
    auto y = sqrt(2 + v * v - u * u + 2 * sqrt(2.0f) * v) / 2.0f -
             sqrt(2 + v * v - u * u - 2 * sqrt(2.0f) * v) / 2.0f;
    return {x, y};
}

}  // namespace ygl

namespace std {
// Hash functor for vector for use with unordered_map
template <typename T, int N>
inline size_t array_hash(const T* v) {
    auto vh = hash<T>();
    auto h = (size_t)0;
    for (auto i = 0; i < N; i++)
        h ^= vh(v[i]) + 0x9e3779b9 + (h << 6) + (h >> 2);
    return h;
}
// Hash functor for vector for use with unordered_map
template <>
struct hash<ygl::vec2i> {
    size_t operator()(const ygl::vec2i& v) const {
        return array_hash<int, 2>(&v.x);
    }
};
template <>
struct hash<ygl::vec3i> {
    size_t operator()(const ygl::vec3i& v) const {
        return array_hash<int, 3>(&v.x);
    }
};
template <>
struct hash<ygl::vec4i> {
    size_t operator()(const ygl::vec4i& v) const {
        return array_hash<int, 4>(&v.x);
    }
};
}  // namespace std

// -----------------------------------------------------------------------------
// MATRICES
// -----------------------------------------------------------------------------
namespace ygl {

// Small Fixed-size square matrices stored in column major format.
struct mat2f {
    vec2f x = {1, 0}, y = {0, 1};
};
struct mat3f {
    vec3f x = {1, 0, 0}, y = {0, 1, 0}, z = {0, 0, 1};
};
struct mat4f {
    vec4f x = {1, 0, 0, 0}, y = {0, 1, 0, 0}, z = {0, 0, 1, 0},
          w = {0, 0, 0, 1};
};

// Identity matrices constants.
const auto identity_mat2f = mat2f();
const auto identity_mat3f = mat3f();
const auto identity_mat4f = mat4f();

// Matrix comparisons.
inline bool operator==(const mat2f& a, const mat2f& b) {
    return a.x == b.x && a.y == b.y;
}
inline bool operator!=(const mat2f& a, const mat2f& b) { return !(a == b); }
inline bool operator==(const mat3f& a, const mat3f& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}
inline bool operator!=(const mat3f& a, const mat3f& b) { return !(a == b); }
inline bool operator==(const mat4f& a, const mat4f& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}
inline bool operator!=(const mat4f& a, const mat4f& b) { return !(a == b); }

// Matrix operations.
inline mat2f operator+(const mat2f& a, const mat2f& b) {
    return {a.x + b.x, a.y + b.y};
}
inline mat2f operator*(const mat2f& a, float b) { return {a.x * b, a.y * b}; }
inline mat2f operator/(const mat2f& a, float b) { return {a.x / b, a.y / b}; }
inline vec2f operator*(const mat2f& a, const vec2f& b) {
    return a.x * b.x + a.y * b.y;
}
inline vec2f operator*(const vec2f& a, const mat2f& b) {
    return {dot(a, b.x), dot(a, b.y)};
}
inline mat2f operator*(const mat2f& a, const mat2f& b) {
    return {a * b.x, a * b.y};
}

// Matrix operations.
inline mat3f operator+(const mat3f& a, const mat3f& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
inline mat3f operator*(const mat3f& a, float b) {
    return {a.x * b, a.y * b, a.z * b};
}
inline mat3f operator/(const mat3f& a, float b) {
    return {a.x / b, a.y / b, a.z / b};
}
inline vec3f operator*(const mat3f& a, const vec3f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline vec3f operator*(const vec3f& a, const mat3f& b) {
    return {dot(a, b.x), dot(a, b.y), dot(a, b.z)};
}
inline mat3f operator*(const mat3f& a, const mat3f& b) {
    return {a * b.x, a * b.y, a * b.z};
}

// Matrix operations.
inline mat4f operator+(const mat4f& a, const mat4f& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}
inline mat4f operator*(const mat4f& a, float b) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
}
inline vec4f operator*(const mat4f& a, const vec4f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}
inline vec4f operator*(const vec4f& a, const mat4f& b) {
    return {dot(a, b.x), dot(a, b.y), dot(a, b.z), dot(a, b.w)};
}
inline mat4f operator*(const mat4f& a, const mat4f& b) {
    return {a * b.x, a * b.y, a * b.z, a * b.w};
}

// Matrix assignments.
inline mat2f& operator+=(mat2f& a, const mat2f& b) { return a = a + b; }
inline mat2f& operator*=(mat2f& a, const mat2f& b) { return a = a * b; }
inline mat2f& operator*=(mat2f& a, float b) { return a = a * b; }
inline mat3f& operator+=(mat3f& a, const mat3f& b) { return a = a + b; }
inline mat3f& operator*=(mat3f& a, const mat3f& b) { return a = a * b; }
inline mat3f& operator*=(mat3f& a, float b) { return a = a * b; }
inline mat4f& operator+=(mat4f& a, const mat4f& b) { return a = a + b; }
inline mat4f& operator*=(mat4f& a, const mat4f& b) { return a = a * b; }
inline mat4f& operator*=(mat4f& a, float b) { return a = a * b; }

// Matrix diagonals and transposes.
inline vec2f diagonal(const mat2f& a) { return {a.x.x, a.y.y}; }
inline vec3f diagonal(const mat3f& a) { return {a.x.x, a.y.y, a.z.z}; }
inline vec4f diagonal(const mat4f& a) { return {a.x.x, a.y.y, a.z.z, a.w.w}; }
inline mat2f transpose(const mat2f& a) {
    return {{a.x.x, a.y.x}, {a.x.y, a.y.y}};
}
inline mat3f transpose(const mat3f& a) {
    return {
        {a.x.x, a.y.x, a.z.x}, {a.x.y, a.y.y, a.z.y}, {a.x.z, a.y.z, a.z.z}};
}
inline mat4f transpose(const mat4f& a) {
    return {{a.x.x, a.y.x, a.z.x, a.w.x}, {a.x.y, a.y.y, a.z.y, a.w.y},
        {a.x.z, a.y.z, a.z.z, a.w.z}, {a.x.w, a.y.w, a.z.w, a.w.w}};
}

// Matrix adjugates, determinant and inverses.
inline mat2f adjugate(const mat2f& a);
inline mat3f adjugate(const mat3f& a);
inline mat4f adjugate(const mat4f& a);
inline float determinant(const mat2f& a);
inline float determinant(const mat3f& a);
inline float determinant(const mat4f& a);
inline mat2f inverse(const mat2f& a) {
    return adjugate(a) * (1 / determinant(a));
}
inline mat3f inverse(const mat3f& a) {
    return adjugate(a) * (1 / determinant(a));
}
inline mat4f inverse(const mat4f& a) {
    return adjugate(a) * (1 / determinant(a));
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// RIGID BODY TRANSFORMS/FRAMES
// -----------------------------------------------------------------------------
namespace ygl {

// Rigid frames stored as a column-major affine transform matrix.
struct frame2f {
    vec2f x = {1, 0}, y = {0, 1}, o = {0, 0};
};
struct frame3f {
    vec3f x = {1, 0, 0}, y = {0, 1, 0}, z = {0, 0, 1}, o = {0, 0, 0};
};

// Indentity frames.
const auto identity_frame2f = frame2f{{1, 0}, {0, 1}, {0, 0}};
const auto identity_frame3f =
    frame3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};

// Frame construction from axis.
inline frame3f make_frame_fromz(const vec3f& o, const vec3f& z_) {
    auto z = normalize(z_);
    auto x = normalize(orthogonal(z));
    auto y = normalize(cross(z, x));
    return {x, y, z, o};
}
inline frame3f make_frame_fromzx(
    const vec3f& o, const vec3f& z_, const vec3f& x_) {
    auto z = normalize(z_);
    auto x = orthonormalize(x_, z);
    auto y = normalize(cross(z, x));
    return {x, y, z, o};
}

// Frame to matrix conversion.
inline mat4f frame_to_mat(const frame3f& a) {
    return {{a.x.x, a.x.y, a.x.z, 0}, {a.y.x, a.y.y, a.y.z, 0},
        {a.z.x, a.z.y, a.z.z, 0}, {a.o.x, a.o.y, a.o.z, 1}};
}
inline frame3f mat_to_frame(const mat4f& a) {
    return {{a.x.x, a.x.y, a.x.z}, {a.y.x, a.y.y, a.y.z}, {a.z.x, a.z.y, a.z.z},
        {a.w.x, a.w.y, a.w.z}};
}

// Frame comparisons.
inline bool operator==(const frame2f& a, const frame2f& b) {
    return a.x == b.x && a.y == b.y && a.o == b.o;
}
inline bool operator!=(const frame2f& a, const frame2f& b) { return !(a == b); }
inline bool operator==(const frame3f& a, const frame3f& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.o == b.o;
}
inline bool operator!=(const frame3f& a, const frame3f& b) { return !(a == b); }

// Frame composition, equivalent to affine matrix product.
inline frame3f operator*(const frame3f& a, const frame3f& b) {
    auto rot = mat3f{a.x, a.y, a.z} * mat3f{b.x, b.y, b.z};
    auto pos = mat3f{a.x, a.y, a.z} * b.o + a.o;
    return {rot.x, rot.y, rot.z, pos};
}
// Frame inverse, equivalent to rigid affine inverse.
inline frame3f inverse(const frame3f& a) {
    auto minv = transpose(mat3f{a.x, a.y, a.z});
    return {minv.x, minv.y, minv.z, -(minv * a.o)};
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// AXIS ALIGNED BOUNDING BOXES
// -----------------------------------------------------------------------------
namespace ygl {

// Axis aligned bounding box represented as a min/max vector pairs.
struct bbox3f {
    vec3f min = {flt_max, flt_max, flt_max}, max = {flt_min, flt_min, flt_min};
};

// Empty bbox constant.
const auto invalid_bbox3f = bbox3f();

// Bounding box comparisons.
inline bool operator==(const bbox3f& a, const bbox3f& b) {
    return a.min == b.min && a.max == b.max;
}
inline bool operator!=(const bbox3f& a, const bbox3f& b) {
    return a.min != b.min || a.max != b.max;
}

// Bounding box expansions with points and other boxes.
inline bbox3f& operator+=(bbox3f& a, const vec3f& b) {
    a.min = {min(a.min.x, b.x), min(a.min.y, b.y), min(a.min.z, b.z)};
    a.max = {max(a.max.x, b.x), max(a.max.y, b.y), max(a.max.z, b.z)};
    return a;
}
inline bbox3f& operator+=(bbox3f& a, const bbox3f& b) {
    a.min = {
        min(a.min.x, b.min.x), min(a.min.y, b.min.y), min(a.min.z, b.min.z)};
    a.max = {
        max(a.max.x, b.max.x), max(a.max.y, b.max.y), max(a.max.z, b.max.z)};
    return a;
}

// Primitive bounds.
inline bbox3f point_bbox(const vec3f& p, float r = 0) {
    auto bbox = invalid_bbox3f;
    bbox += p - vec3f{r, r, r};
    bbox += p + vec3f{r, r, r};
    return bbox;
}
inline bbox3f line_bbox(
    const vec3f& v0, const vec3f& v1, float r0 = 0, float r1 = 0) {
    auto bbox = invalid_bbox3f;
    bbox += v0 - vec3f{r0, r0, r0};
    bbox += v0 + vec3f{r0, r0, r0};
    bbox += v1 - vec3f{r1, r1, r1};
    bbox += v1 + vec3f{r1, r1, r1};
    return bbox;
}
inline bbox3f triangle_bbox(const vec3f& v0, const vec3f& v1, const vec3f& v2) {
    auto bbox = invalid_bbox3f;
    bbox += v0;
    bbox += v1;
    bbox += v2;
    return bbox;
}
inline bbox3f quad_bbox(
    const vec3f& v0, const vec3f& v1, const vec3f& v2, const vec3f& v3) {
    auto bbox = invalid_bbox3f;
    bbox += v0;
    bbox += v1;
    bbox += v2;
    bbox += v3;
    return bbox;
}
inline bbox3f tetrahedron_bbox(
    const vec3f& v0, const vec3f& v1, const vec3f& v2, const vec3f& v3) {
    auto bbox = invalid_bbox3f;
    bbox += v0;
    bbox += v1;
    bbox += v2;
    bbox += v3;
    return bbox;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// RAYS
// -----------------------------------------------------------------------------
namespace ygl {

// Rays with origin, direction and min/max t value.
struct ray3f {
    vec3f o = {0, 0, 0}, d = {0, 0, 1};
    float tmin = 0, tmax = flt_max;
};

// Construct a ray from dirction or segments using a default epsilon.
inline ray3f make_ray(const vec3f& o, const vec3f& d, float eps = 1e-4f) {
    return ray3f{o, d, eps, flt_max};
}
inline ray3f make_segment(const vec3f& p1, const vec3f& p2, float eps = 1e-4f) {
    return ray3f{p1, normalize(p2 - p1), eps, length(p2 - p1) - 2 * eps};
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// TRANSFORMS
// -----------------------------------------------------------------------------
namespace ygl {

// Transforms points, vectors and directions by matrices.
inline vec2f transform_point(const mat3f& a, const vec2f& b) {
    auto tvb = a * vec3f{b.x, b.y, 1};
    return vec2f{tvb.x, tvb.y} / tvb.z;
}
inline vec3f transform_point(const mat4f& a, const vec3f& b) {
    auto tvb = a * vec4f{b.x, b.y, b.z, 1};
    return vec3f{tvb.x, tvb.y, tvb.z} / tvb.w;
}
inline vec2f transform_vector(const mat3f& a, const vec2f& b) {
    auto tvb = a * vec3f{b.x, b.y, 0};
    return vec2f{tvb.x, tvb.y} / tvb.z;
}
inline vec3f transform_vector(const mat3f& a, const vec3f& b) { return a * b; }
inline vec3f transform_vector(const mat4f& a, const vec3f& b) {
    auto tvb = a * vec4f{b.x, b.y, b.z, 0};
    return vec3f{tvb.x, tvb.y, tvb.z};
}
inline vec3f transform_direction(const mat4f& a, const vec3f& b) {
    return normalize(transform_vector(a, b));
}

// Transforms points, vectors and directions by frames.
inline vec2f transform_point(const frame2f& a, const vec2f& b) {
    return a.x * b.x + a.y * b.y + a.o;
}
inline vec3f transform_point(const frame3f& a, const vec3f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.o;
}
inline vec2f transform_vector(const frame2f& a, const vec2f& b) {
    return a.x * b.x + a.y * b.y;
}
inline vec3f transform_vector(const frame3f& a, const vec3f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline vec3f transform_direction(const frame3f& a, const vec3f& b) {
    return normalize(transform_vector(a, b));
}

// Transforms rays and bounding boxes by matrices.
inline ray3f transform_ray(const frame3f& a, const ray3f& b) {
    return {transform_point(a, b.o), transform_vector(a, b.d), b.tmin, b.tmax};
}
inline ray3f transform_ray(const mat4f& a, const ray3f& b) {
    return {transform_point(a, b.o), transform_vector(a, b.d), b.tmin, b.tmax};
}
inline bbox3f transform_bbox(const frame3f& a, const bbox3f& b) {
    auto corners = {vec3f{b.min.x, b.min.y, b.min.z},
        vec3f{b.min.x, b.min.y, b.max.z}, vec3f{b.min.x, b.max.y, b.min.z},
        vec3f{b.min.x, b.max.y, b.max.z}, vec3f{b.max.x, b.min.y, b.min.z},
        vec3f{b.max.x, b.min.y, b.max.z}, vec3f{b.max.x, b.max.y, b.min.z},
        vec3f{b.max.x, b.max.y, b.max.z}};
    auto xformed = bbox3f();
    for (auto& corner : corners) xformed += transform_point(a, corner);
    return xformed;
}
inline bbox3f transform_bbox(const mat4f& a, const bbox3f& b) {
    auto corners = {vec3f{b.min.x, b.min.y, b.min.z},
        vec3f{b.min.x, b.min.y, b.max.z}, vec3f{b.min.x, b.max.y, b.min.z},
        vec3f{b.min.x, b.max.y, b.max.z}, vec3f{b.max.x, b.min.y, b.min.z},
        vec3f{b.max.x, b.min.y, b.max.z}, vec3f{b.max.x, b.max.y, b.min.z},
        vec3f{b.max.x, b.max.y, b.max.z}};
    auto xformed = bbox3f();
    for (auto& corner : corners) xformed += transform_point(a, corner);
    return xformed;
}

// Inverse transforms by frames, assuming they are rigid transforms.
inline vec2f transform_point_inverse(const frame2f& a, const vec2f& b) {
    return {dot(b - a.o, a.x), dot(b - a.o, a.y)};
}
inline vec3f transform_point_inverse(const frame3f& a, const vec3f& b) {
    return {dot(b - a.o, a.x), dot(b - a.o, a.y), dot(b - a.o, a.z)};
}
inline vec2f transform_vector_inverse(const frame2f& a, const vec2f& b) {
    return {dot(b, a.x), dot(b, a.y)};
}
inline vec3f transform_vector_inverse(const frame3f& a, const vec3f& b) {
    return {dot(b, a.x), dot(b, a.y), dot(b, a.z)};
}
inline vec3f transform_direction_inverse(const frame3f& a, const vec3f& b) {
    return normalize(transform_vector_inverse(a, b));
}
inline ray3f transform_ray_inverse(const frame3f& a, const ray3f& b) {
    return {transform_point_inverse(a, b.o),
        transform_direction_inverse(a, b.d), b.tmin, b.tmax};
}
inline bbox3f transform_bbox_inverse(const frame3f& a, const bbox3f& b) {
    return transform_bbox(inverse(a), b);
}

// Translation, scaling and rotations transforms.
inline frame3f translation_frame(const vec3f& a) {
    return {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, a};
}
inline frame3f scaling_frame(const vec3f& a) {
    return {{a.x, 0, 0}, {0, a.y, 0}, {0, 0, a.z}, {0, 0, 0}};
}
inline frame3f rotation_frame(const vec3f& axis, float angle) {
    auto s = sin(angle), c = cos(angle);
    auto vv = normalize(axis);
    return {{c + (1 - c) * vv.x * vv.x, (1 - c) * vv.x * vv.y + s * vv.z,
                (1 - c) * vv.x * vv.z - s * vv.y},
        {(1 - c) * vv.x * vv.y - s * vv.z, c + (1 - c) * vv.y * vv.y,
            (1 - c) * vv.y * vv.z + s * vv.x},
        {(1 - c) * vv.x * vv.z + s * vv.y, (1 - c) * vv.y * vv.z - s * vv.x,
            c + (1 - c) * vv.z * vv.z},
        {0, 0, 0}};
}
inline frame3f rotation_frame(const vec4f& quat) {
    auto v = quat;
    return {{v.w * v.w + v.x * v.x - v.y * v.y - v.z * v.z,
                (v.x * v.y + v.z * v.w) * 2, (v.z * v.x - v.y * v.w) * 2},
        {(v.x * v.y - v.z * v.w) * 2,
            v.w * v.w - v.x * v.x + v.y * v.y - v.z * v.z,
            (v.y * v.z + v.x * v.w) * 2},
        {(v.z * v.x + v.y * v.w) * 2, (v.y * v.z - v.x * v.w) * 2,
            v.w * v.w - v.x * v.x - v.y * v.y + v.z * v.z},
        {0, 0, 0}};
}
inline frame3f rotation_frame(const mat3f& rot) {
    return {rot.x, rot.y, rot.z, {0, 0, 0}};
}

// Lookat frame. Z-axis can be inverted with inv_xz.
inline frame3f lookat_frame(const vec3f& eye, const vec3f& center,
    const vec3f& up, bool inv_xz = false) {
    auto w = normalize(eye - center);
    auto u = normalize(cross(up, w));
    auto v = normalize(cross(w, u));
    if (inv_xz) {
        w = -w;
        u = -u;
    }
    return {u, v, w, eye};
}

// OpenGL frustum, ortho and perspecgive matrices.
inline mat4f frustum_mat(float l, float r, float b, float t, float n, float f) {
    return {{2 * n / (r - l), 0, 0, 0}, {0, 2 * n / (t - b), 0, 0},
        {(r + l) / (r - l), (t + b) / (t - b), -(f + n) / (f - n), -1},
        {0, 0, -2 * f * n / (f - n), 0}};
}
inline mat4f ortho_mat(float l, float r, float b, float t, float n, float f) {
    return {{2 / (r - l), 0, 0, 0}, {0, 2 / (t - b), 0, 0},
        {0, 0, -2 / (f - n), 0},
        {-(r + l) / (r - l), -(t + b) / (t - b), -(f + n) / (f - n), 1}};
}
inline mat4f ortho2d_mat(float left, float right, float bottom, float top) {
    return ortho_mat(left, right, bottom, top, -1, 1);
}
inline mat4f ortho_mat(float xmag, float ymag, float near, float far) {
    return {{1 / xmag, 0, 0, 0}, {0, 1 / ymag, 0, 0},
        {0, 0, 2 / (near - far), 0}, {0, 0, (far + near) / (near - far), 1}};
}
inline mat4f perspective_mat(float fovy, float aspect, float near, float far) {
    auto tg = tan(fovy / 2);
    return {{1 / (aspect * tg), 0, 0, 0}, {0, 1 / tg, 0, 0},
        {0, 0, (far + near) / (near - far), -1},
        {0, 0, 2 * far * near / (near - far), 0}};
}
inline mat4f perspective_mat(float fovy, float aspect, float near) {
    auto tg = tan(fovy / 2);
    return {{1 / (aspect * tg), 0, 0, 0}, {0, 1 / tg, 0, 0}, {0, 0, -1, -1},
        {0, 0, 2 * near, 0}};
}

// Rotation conversions.
inline std::pair<vec3f, float> rotation_axisangle(const vec4f& quat) {
    return {normalize(vec3f{quat.x, quat.y, quat.z}), 2 * acos(quat.w)};
}
inline vec4f rotation_quat(const vec3f& axis, float angle) {
    auto len = length(axis);
    if (!len) return {0, 0, 0, 1};
    return vec4f{sin(angle / 2) * axis.x / len, sin(angle / 2) * axis.y / len,
        sin(angle / 2) * axis.z / len, cos(angle / 2)};
}

// Turntable and FPS Camera navigation.
void camera_turntable(vec3f& from, vec3f& to, vec3f& up, const vec2f& rotate,
    float dolly, const vec2f& pan);
void camera_turntable(frame3f& frame, float& focus, const vec2f& rotate,
    float dolly, const vec2f& pan);
void camera_fps(frame3f& frame, const vec3f& transl, const vec2f& rotate);

}  // namespace ygl

// -----------------------------------------------------------------------------
// RANDOM NUMBER GENERATION
// -----------------------------------------------------------------------------
namespace ygl {

// PCG random numbers from http://www.pcg-random.org/
struct rng_state {
    uint64_t state = 0x853c49e6748fea9bULL;
    uint64_t inc = 0xda3e39cb94b95bdbULL;
};

// Next random number.
inline uint32_t advance_rng(rng_state& rng) {
    uint64_t oldstate = rng.state;
    rng.state = oldstate * 6364136223846793005ULL + rng.inc;
    uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
    uint32_t rot = (uint32_t)(oldstate >> 59u);
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

// Init a random number generator with a state state from the sequence seq.
inline rng_state make_rng(uint64_t state, uint64_t seq = 1) {
    auto rng = rng_state();
    rng.state = 0U;
    rng.inc = (seq << 1u) | 1u;
    advance_rng(rng);
    rng.state += state;
    advance_rng(rng);
    return rng;
}

// Next random numbers: floats in [0,1), ints in [0,n).
inline uint32_t next_rand1i(rng_state& rng, int n) {
    return advance_rng(rng) % n;
}
inline float next_rand1f(rng_state& rng) {
    union {
        uint32_t u;
        float f;
    } x;
    x.u = (advance_rng(rng) >> 9) | 0x3f800000u;
    return x.f - 1.0f;
    // alternate implementation
    // const static auto scale = (float)(1.0 / numeric_limits<uint32_t>::max());
    // return advance_rng(rng) * scale;
}
inline vec2f next_rand2f(rng_state& rng) {
    return {next_rand1f(rng), next_rand1f(rng)};
}
inline vec3f next_rand3f(rng_state& rng) {
    return {next_rand1f(rng), next_rand1f(rng), next_rand1f(rng)};
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// MONETACARLO SAMPLING FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Sample an hemispherical direction with uniform distribution.
inline vec3f sample_hemisphere(const vec2f& ruv) {
    auto z = ruv.y;
    auto r = sqrt(1 - z * z);
    auto phi = 2 * pi * ruv.x;
    return {r * cos(phi), r * sin(phi), z};
}
inline float sample_hemisphere_pdf(const vec3f& w) {
    return (w.z <= 0) ? 0 : 1 / (2 * pi);
}

// Sample a spherical direction with uniform distribution.
inline vec3f sample_sphere(const vec2f& ruv) {
    auto z = 2 * ruv.y - 1;
    auto r = sqrt(1 - z * z);
    auto phi = 2 * pi * ruv.x;
    return {r * cos(phi), r * sin(phi), z};
}
inline float sample_sphere_pdf(const vec3f& w) { return 1 / (4 * pi); }

// Sample spherical coordinates uniformly.
inline vec2f sample_spherical(const vec2f& ruv) {
    // BUG: FIXME this is not uniform at all!!!!
    return {ruv.x, ruv.y};
}
inline float sample_spherical_pdf(const vec2f& w) { return 1 / (4 * pi); }

// Sample an hemispherical direction with cosine distribution.
inline vec3f sample_hemisphere_cosine(const vec2f& ruv) {
    auto z = sqrt(ruv.y);
    auto r = sqrt(1 - z * z);
    auto phi = 2 * pi * ruv.x;
    return {r * cos(phi), r * sin(phi), z};
}
inline float sample_hemisphere_cosine_pdf(const vec3f& w) {
    return (w.z <= 0) ? 0 : w.z / pi;
}

// Sample an hemispherical direction with cosine power distribution.
inline vec3f sample_hemisphere_cospower(float n, const vec2f& ruv) {
    auto z = pow(ruv.y, 1 / (n + 1));
    auto r = sqrt(1 - z * z);
    auto phi = 2 * pi * ruv.x;
    return {r * cos(phi), r * sin(phi), z};
}
inline float sample_hemisphere_cospower_pdf(float n, const vec3f& w) {
    return (w.z <= 0) ? 0 : pow(w.z, n) * (n + 1) / (2 * pi);
}

// Sample a point uniformly on a disk.
inline vec3f sample_disk(const vec2f& ruv) {
    auto r = sqrt(ruv.y);
    auto phi = 2 * pi * ruv.x;
    return {cos(phi) * r, sin(phi) * r, 0};
}
inline float sample_disk_pdf() { return 1 / pi; }

// Sample a point uniformly on a cylinder, without caps.
inline vec3f sample_cylinder(const vec2f& ruv) {
    auto phi = 2 * pi * ruv.x;
    return {sin(phi), cos(phi), ruv.y * 2 - 1};
}
inline float sample_cylinder_pdf() { return 1 / pi; }

// Sample a point uniformly on a triangle.
inline vec2f sample_triangle(const vec2f& ruv) {
    return {1 - sqrt(ruv.x), ruv.y * sqrt(ruv.x)};
}
inline vec3f sample_triangle(
    const vec3f& v0, const vec3f& v1, const vec3f& v2, const vec2f& ruv) {
    auto uv = sample_triangle(ruv);
    return v0 * (1 - uv.x - uv.y) + v1 * uv.x + v2 * uv.y;
}
// Pdf for uniform triangle sampling, i.e. triangle area.
inline float sample_triangle_pdf(
    const vec3f& v0, const vec3f& v1, const vec3f& v2) {
    return 2 / length(cross(v1 - v0, v2 - v0));
}

// Sample an index with uniform distribution.
inline int sample_index(int size, float r) {
    return clamp((int)(r * size), 0, size - 1);
}
inline float sample_index_pdf(int size) { return 1.0f / size; }

// Sample a discrete distribution represented by its cdf.
inline int sample_discrete(const std::vector<float>& cdf, float r) {
    // todo: implement binary search better
    r = clamp(r * cdf.back(), 0.0f, cdf.back() - 0.00001f);
    for (auto i = 0; i < cdf.size(); i++) {
        if (cdf[i] > r) return i;
    }
    return (int)cdf.size() - 1;
}
// Pdf for uniform discrete distribution sampling.
inline float sample_discrete_pdf(const std::vector<float>& cdf, int idx) {
    if (idx == 0) return cdf.at(0);
    return cdf.at(idx) - cdf.at(idx - 1);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// PERLIN NOISE FUNCTION
// -----------------------------------------------------------------------------
namespace ygl {

// Compute the revised Pelin noise function. Wrap provides a wrapping noise
// but must be power of two (wraps at 256 anyway). For octave based noise,
// good values are obtained with octaves=6 (numerber of noise calls),
// lacunarity=~2.0 (spacing between successive octaves: 2.0 for warpping
// output), gain=0.5 (relative weighting applied to each successive octave),
// offset=1.0 (used to invert the ridges).
float perlin_noise(const vec3f& p, const vec3i& wrap = zero3i);
float perlin_ridge_noise(const vec3f& p, float lacunarity = 2.0f,
    float gain = 0.5f, float offset = 1.0f, int octaves = 6,
    const vec3i& wrap = zero3i);
float perlin_fbm_noise(const vec3f& p, float lacunarity = 2.0f,
    float gain = 0.5f, int octaves = 6, const vec3i& wrap = zero3i);
float perlin_turbulence_noise(const vec3f& p, float lacunarity = 2.0f,
    float gain = 0.5f, int octaves = 6, const vec3i& wrap = zero3i);

}  // namespace ygl

// -----------------------------------------------------------------------------
// GEOMETRY UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Line properties.
inline vec3f line_tangent(const vec3f& v0, const vec3f& v1) {
    return normalize(v1 - v0);
}
inline float line_length(const vec3f& v0, const vec3f& v1) {
    return length(v1 - v0);
}

// Triangle properties.
inline vec3f triangle_normal(
    const vec3f& v0, const vec3f& v1, const vec3f& v2) {
    return normalize(cross(v1 - v0, v2 - v0));
}
inline float triangle_area(const vec3f& v0, const vec3f& v1, const vec3f& v2) {
    return length(cross(v1 - v0, v2 - v0)) / 2;
}

// Quad propeties.
inline vec3f quad_normal(
    const vec3f& v0, const vec3f& v1, const vec3f& v2, const vec3f& v3) {
    return normalize(triangle_normal(v0, v1, v3) + triangle_normal(v2, v3, v1));
}
inline float quad_area(
    const vec3f& v0, const vec3f& v1, const vec3f& v2, const vec3f& v3) {
    return triangle_area(v0, v1, v3) + triangle_area(v2, v3, v1);
}

// Triangle tangent and bitangent from uv
inline std::pair<vec3f, vec3f> triangle_tangents_fromuv(const vec3f& v0,
    const vec3f& v1, const vec3f& v2, const vec2f& uv0, const vec2f& uv1,
    const vec2f& uv2) {
    // Follows the definition in http://www.terathon.com/code/tangent.html and
    // https://gist.github.com/aras-p/2843984
    // normal points up from texture space
    auto p = v1 - v0;
    auto q = v2 - v0;
    auto s = vec2f{uv1.x - uv0.x, uv2.x - uv0.x};
    auto t = vec2f{uv1.y - uv0.y, uv2.y - uv0.y};
    auto div = s.x * t.y - s.y * t.x;

    if (div != 0) {
        auto tu = vec3f{t.y * p.x - t.x * q.x, t.y * p.y - t.x * q.y,
                      t.y * p.z - t.x * q.z} /
                  div;
        auto tv = vec3f{s.x * q.x - s.y * p.x, s.x * q.y - s.y * p.y,
                      s.x * q.z - s.y * p.z} /
                  div;
        return {tu, tv};
    } else {
        return {{1, 0, 0}, {0, 1, 0}};
    }
}

// Copies of point value. Here only for completeness.
template <typename T>
inline T interpolate_point(const std::vector<T>& vals, int p) {
    if (vals.empty()) return T();
    return vals[p];
}

// Interpolates values over a line parametrized from a to b by u. Same as lerp.
template <typename T, typename T1>
inline T interpolate_line(const T& v0, const T& v1, const T1 u) {
    return v0 * (1 - u) + v1 * u;
}
template <typename T, typename T1>
inline T interpolate_line(const std::vector<T>& vals, const vec2i& l, T1 u) {
    if (vals.empty()) return T();
    return vals[l.x] * (1 - u) + vals[l.y] * u;
}

// Interpolates values over a triangle parametrized by u and v along the
// (v1-v0) and (v2-v0) directions. Same as barycentric interpolation.
template <typename T>
inline T interpolate_triangle(
    const T& v0, const T& v1, const T& v2, const vec2f& uv) {
    return v0 * (1 - uv.x - uv.y) + v1 * uv.x + v2 * uv.y;
}
template <typename T>
inline T interpolate_triangle(
    const std::vector<T>& vals, const vec3i& t, const vec2f& uv) {
    if (vals.empty()) return T();
    return vals[t.x] * (1 - uv.x - uv.y) + vals[t.y] * uv.x + vals[t.z] * uv.y;
}
// Interpolates values over a quad parametrized by u and v along the
// (v1-v0) and (v2-v1) directions. Same as bilear interpolation.
template <typename T>
inline T interpolate_quad(
    const T& v0, const T& v1, const T& v2, const T& v3, const vec2f& uv) {
    return v0 * (1 - uv.x) * (1 - uv.y) + v1 * uv.x * (1 - uv.y) +
           v2 * uv.x * uv.y + v3 * (1 - uv.x) * uv.y;
}
template <typename T>
inline T interpolate_quad(
    const std::vector<T>& vals, const vec4i& t, const vec2f& uv) {
    if (vals.empty()) return T();
    return vals[t.x] * (1 - uv.x) * (1 - uv.y) + vals[t.y] * uv.x * (1 - uv.y) +
           vals[t.z] * uv.x * uv.y + vals[t.w] * (1 - uv.x) * uv.y;
}

// Evaluates the i-th Bernstein polynomial of degree degree at u.
template <typename T>
inline T eval_bernstein(T u, int i, int degree) {
    if (i < 0 or i > degree) return 0;
    if (degree == 0)
        return 1;
    else if (degree == 1) {
        if (i == 0)
            return 1 - u;
        else if (i == 1)
            return u;
    } else if (degree == 2) {
        if (i == 0)
            return (1 - u) * (1 - u);
        else if (i == 1)
            return 2 * u * (1 - u);
        else if (i == 2)
            return u * u;
    } else if (degree == 3) {
        if (i == 0)
            return (1 - u) * (1 - u) * (1 - u);
        else if (i == 1)
            return 3 * u * (1 - u) * (1 - u);
        else if (i == 2)
            return 3 * u * u * (1 - u);
        else if (i == 3)
            return u * u * u;
    } else {
        return (1 - u) * eval_bernstein(u, i, degree - 1) +
               u * eval_bernstein(u, i - 1, degree - 1);
    }
    return 0;
}

// Evaluates the derivative of i-th Bernstein polynomial of degree degree at u.
template <typename T>
inline T eval_bernstein_derivative(T u, int i, int degree) {
    return degree * (eval_bernstein(u, i - 1, degree - 1) -
                        eval_bernstein(u, i, degree - 1));
}

// Interpolates values along a cubic Bezier segment parametrized by u.
template <typename T>
inline T interpolate_bezier(
    const T& v0, const T& v1, const T& v2, const T& v3, float u) {
    return v0 * (1 - u) * (1 - u) * (1 - u) + v1 * 3 * u * (1 - u) * (1 - u) +
           v2 * 3 * u * u * (1 - u) + v3 * u * u * u;
}
template <typename T>
inline T interpolate_bezier(
    const std::vector<T>& vals, const vec4i& b, float u) {
    if (vals.empty()) return T();
    return interpolate_bezier(vals[b.x], vals[b.y], vals[b.z], vals[b.w], u);
}
// Computes the derivative of a cubic Bezier segment parametrized by u.
template <typename T>
inline T interpolate_bezier_derivative(
    const T& v0, const T& v1, const T& v2, const T& v3, float u) {
    return (v1 - v0) * 3 * (1 - u) * (1 - u) + (v2 - v1) * 6 * u * (1 - u) +
           (v3 - v2) * 3 * u * u;
}
template <typename T>
inline T interpolate_bezier_derivative(
    const std::vector<T>& vals, const vec4i& b, float u) {
    if (vals.empty()) return T();
    return interpolate_bezier_derivative(
        vals[b.x], vals[b.y], vals[b.z], vals[b.w], u);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// ANIMATION UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Find the first keyframe value that is greater than the argument.
inline int eval_keyframed_index(
    const std::vector<float>& times, const float& time) {
    for (auto i = 0; i < times.size(); i++)
        if (times[i] > time) return i;
    return (int)times.size();
}

// Evalautes a keyframed value using step interpolation.
template <typename T>
inline T eval_keyframed_step(
    const std::vector<float>& times, const std::vector<T>& vals, float time) {
    if (time <= times.front()) return vals.front();
    if (time >= times.back()) return vals.back();
    time = clamp(time, times.front(), times.back() - 0.001f);
    auto idx = eval_keyframed_index(times, time);
    return vals.at(idx - 1);
}

// Evalautes a keyframed value using linear interpolation.
inline vec4f eval_keyframed_slerp(const std::vector<float>& times,
    const std::vector<vec4f>& vals, float time) {
    if (time <= times.front()) return vals.front();
    if (time >= times.back()) return vals.back();
    time = clamp(time, times.front(), times.back() - 0.001f);
    auto idx = eval_keyframed_index(times, time);
    auto t = (time - times.at(idx - 1)) / (times.at(idx) - times.at(idx - 1));
    return slerp(vals.at(idx - 1), vals.at(idx), t);
}

// Evalautes a keyframed value using linear interpolation.
template <typename T>
inline T eval_keyframed_linear(
    const std::vector<float>& times, const std::vector<T>& vals, float time) {
    if (time <= times.front()) return vals.front();
    if (time >= times.back()) return vals.back();
    time = clamp(time, times.front(), times.back() - 0.001f);
    auto idx = eval_keyframed_index(times, time);
    auto t = (time - times.at(idx - 1)) / (times.at(idx) - times.at(idx - 1));
    return vals.at(idx - 1) * (1 - t) + vals.at(idx) * t;
}

// Evalautes a keyframed value using Bezier interpolation.
template <typename T>
inline T eval_keyframed_bezier(
    const std::vector<float>& times, const std::vector<T>& vals, float time) {
    if (time <= times.front()) return vals.front();
    if (time >= times.back()) return vals.back();
    time = clamp(time, times.front(), times.back() - 0.001f);
    auto idx = eval_keyframed_index(times, time);
    auto t = (time - times.at(idx - 1)) / (times.at(idx) - times.at(idx - 1));
    return interpolate_bezier(
        vals.at(idx - 3), vals.at(idx - 2), vals.at(idx - 1), vals.at(idx), t);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// SHAPE UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Compute per-vertex normals/tangents for lines/triangles/quads.
void compute_tangents(const std::vector<vec2i>& lines,
    const std::vector<vec3f>& pos, std::vector<vec3f>& tang,
    bool weighted = true);
void compute_normals(const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, std::vector<vec3f>& norm,
    bool weighted = true);
void compute_normals(const std::vector<vec4i>& quads,
    const std::vector<vec3f>& pos, std::vector<vec3f>& norm,
    bool weighted = true);

// Compute per-vertex tangent frames for triangle meshes.
// Tangent space is defined by a four component vector.
// The first three components are the tangent with respect to the u texcoord.
// The fourth component is the sign of the tangent wrt the v texcoord.
// Tangent frame is useful in normal mapping.
void compute_tangent_frames(const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, const std::vector<vec3f>& norm,
    const std::vector<vec2f>& texcoord, std::vector<vec4f>& tangsp,
    bool weighted = true);

// Apply skinning to vertex position and normals.
void compute_skinning(const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec4f>& weights,
    const std::vector<vec4i>& joints, const std::vector<mat4f>& xforms,
    std::vector<vec3f>& skinned_pos, std::vector<vec3f>& skinned_norm);
void compute_skinning(const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec4f>& weights,
    const std::vector<vec4i>& joints, const std::vector<frame3f>& xforms,
    std::vector<vec3f>& skinned_pos, std::vector<vec3f>& skinned_norm);
// Apply skinning as specified in Khronos glTF.
void compute_matrix_skinning(const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec4f>& weights,
    const std::vector<vec4i>& joints, const std::vector<mat4f>& xforms,
    std::vector<vec3f>& skinned_pos, std::vector<vec3f>& skinned_norm);

// Create an array of edges.
std::vector<vec2i> get_edges(const std::vector<vec2i>& lines,
    const std::vector<vec3i>& triangles, const std::vector<vec4i>& quads);

// Convert quads to triangles
std::vector<vec3i> convert_quads_to_triangles(const std::vector<vec4i>& quads);
// Convert quads to triangles with a diamond-like topology.
// Quads have to be consecutive one row after another.
std::vector<vec3i> convert_quads_to_triangles(
    const std::vector<vec4i>& quads, int row_length);

// Convert beziers to lines using 3 lines for each bezier.
std::vector<vec2i> convert_bezier_to_lines(const std::vector<vec4i>& beziers);

// Convert face-varying data to single primitives. Returns the quads indices
// and face ids and filled vectors for pos, norm and texcoord.
std::tuple<std::vector<vec4i>, std::vector<vec3f>, std::vector<vec3f>,
    std::vector<vec2f>>
convert_face_varying(const std::vector<vec4i>& quads_pos,
    const std::vector<vec4i>& quads_norm,
    const std::vector<vec4i>& quads_texcoord, const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec2f>& texcoord);

// Subdivide lines by splitting each line in half.
template <typename T>
void subdivide_lines(
    std::vector<vec2i>& lines, std::vector<T>& vert, bool update_lines = true);
void subdivide_lines(std::vector<vec2i>& lines, std::vector<vec3f>& pos,
    std::vector<vec3f>& tang, std::vector<vec2f>& texcoord,
    std::vector<vec4f>& color, std::vector<float>& radius);
// Subdivide triangle by splitting each triangle in four, creating new
// vertices for each edge.
template <typename T>
void subdivide_triangles(std::vector<vec3i>& triangles, std::vector<T>& vert,
    bool update_triangles = true);
void subdivide_triangles(std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<vec4f>& color, std::vector<float>& radius);
// Subdivide quads by splitting each quads in four, creating new
// vertices for each edge and for each face.
template <typename T>
void subdivide_quads(
    std::vector<vec4i>& quads, std::vector<T>& vert, bool update_quads = true);
void subdivide_quads(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<vec4f>& color, std::vector<float>& radius);
// Subdivide beziers by splitting each segment in two.
template <typename T>
void subdivide_beziers(std::vector<vec4i>& beziers, std::vector<T>& vert,
    bool update_beziers = true);
void subdivide_beziers(std::vector<vec4i>& beziers, std::vector<vec3f>& pos,
    std::vector<vec3f>& tang, std::vector<vec2f>& texcoord,
    std::vector<vec4f>& color, std::vector<float>& radius);
// Subdivide quads using Carmull-Clark subdivision rules.
template <typename T>
void subdivide_catmullclark(std::vector<vec4i>& beziers, std::vector<T>& vert,
    bool update_quads = true);
void subdivide_catmullclark(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<vec4f>& color, std::vector<float>& radius);

// Merge lines between shapes.
void merge_lines(std::vector<vec2i>& lines, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    const std::vector<vec2i>& lines1, const std::vector<vec3f>& pos1,
    const std::vector<vec3f>& norm1, const std::vector<vec2f>& texcoord1);
// Merge triangles between shapes.
void merge_triangles(std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    const std::vector<vec3i>& triangles1, const std::vector<vec3f>& pos1,
    const std::vector<vec3f>& norm1, const std::vector<vec2f>& texcoord1);
// Merge quads between shapes.
void merge_quads(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    const std::vector<vec4i>& quads1, const std::vector<vec3f>& pos1,
    const std::vector<vec3f>& norm1, const std::vector<vec2f>& texcoord1);

// Pick a point in a point set uniformly.
inline int sample_points(int npoints, float re) {
    return sample_index(npoints, re);
}
inline std::vector<float> sample_points_cdf(int npoints) {
    auto cdf = std::vector<float>(npoints);
    for (auto i = 0; i < cdf.size(); i++) cdf[i] = 1 + (i ? cdf[i - 1] : 0);
    return cdf;
}
inline int sample_points(const std::vector<float>& cdf, float re) {
    return sample_discrete(cdf, re);
}

// Pick a point on lines uniformly.
inline std::vector<float> sample_lines_cdf(
    const std::vector<vec2i>& lines, const std::vector<vec3f>& pos) {
    auto cdf = std::vector<float>(lines.size());
    for (auto i = 0; i < cdf.size(); i++) {
        auto l = lines[i];
        auto w = line_length(pos[l.x], pos[l.y]);
        cdf[i] = w + (i ? cdf[i - 1] : 0);
    }
    return cdf;
}
inline std::pair<int, float> sample_lines(
    const std::vector<float>& cdf, float re, float ru) {
    return {sample_discrete(cdf, re), ru};
}

// Pick a point on a triangle mesh uniformly.
inline std::vector<float> sample_triangles_cdf(
    const std::vector<vec3i>& triangles, const std::vector<vec3f>& pos) {
    auto cdf = std::vector<float>(triangles.size());
    for (auto i = 0; i < cdf.size(); i++) {
        auto t = triangles[i];
        auto w = triangle_area(pos[t.x], pos[t.y], pos[t.z]);
        cdf[i] = w + (i ? cdf[i - 1] : 0);
    }
    return cdf;
}
inline std::pair<int, vec2f> sample_triangles(
    const std::vector<float>& cdf, float re, const vec2f& ruv) {
    return {sample_discrete(cdf, re), sample_triangle(ruv)};
}

// Pick a point on a quad mesh uniformly.
inline std::vector<float> sample_quads_cdf(
    const std::vector<vec4i>& quads, const std::vector<vec3f>& pos) {
    auto cdf = std::vector<float>(quads.size());
    for (auto i = 0; i < cdf.size(); i++) {
        auto q = quads[i];
        auto w = quad_area(pos[q.x], pos[q.y], pos[q.z], pos[q.w]);
        cdf[i] = w + (i ? cdf[i - 1] : 0);
    }
    return cdf;
}
inline std::pair<int, vec2f> sample_quads(
    const std::vector<float>& cdf, float re, const vec2f& ruv) {
    return {sample_discrete(cdf, re), ruv};
}

// Samples a set of points over a triangle mesh uniformly. Returns pos, norm
// and tecoord of the sampled points.
std::tuple<std::vector<vec3f>, std::vector<vec3f>, std::vector<vec2f>>
sample_triangles_points(const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, const std::vector<vec3f>& norm,
    const std::vector<vec2f>& texcoord, int npoints, uint64_t seed = 0);

// Make examples shapes.
void make_quad(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec2i& steps,
    const vec2f& size, const vec2f& uvsize);
void make_cube(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec3i& steps,
    const vec3f& size, const vec3f& uvsize);
void make_cube_rounded(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec3i& steps,
    const vec3f& size, const vec3f& uvsize, float radius);
void make_sphere(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec2i& steps,
    float size, const vec2f& uvsize);
void make_sphere_cube(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, int steps,
    float size, float uvsize);
void make_sphere_flipcap(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec2i& steps,
    float size, const vec2f& uvsize, const vec2f& zflip);
void make_disk(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec2i& steps,
    float size, const vec2f& uvsize);
void make_disk_quad(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, int steps,
    float size, float uvsize);
void make_cylinder_side(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec2i& steps,
    const vec2f& size, const vec2f& uvsize, bool capped);
void make_cylinder(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec3i& steps,
    const vec2f& size, const vec3f& uvsize);
void make_cylinder_rounded(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec3i& steps,
    const vec2f& size, const vec3f& uvsize, float radius);
void make_sphere(
    std::vector<vec4i>& quads, std::vector<vec3f>& pos, int tesselation);
void make_geodesic_sphere(
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos, int tesselation);
void make_cube(
    std::vector<vec4i>& quads, std::vector<vec3f>& pos, int tesselation);
void make_fvcube(std::vector<vec4i>& quads_pos, std::vector<vec3f>& pos,
    std::vector<vec4i>& quads_norm, std::vector<vec3f>& norm,
    std::vector<vec4i>& quads_texcoord, std::vector<vec2f>& texcoord,
    int tesselation);
void make_suzanne(
    std::vector<vec4i>& quads, std::vector<vec3f>& pos, int tesselation);

// Generate lines set along a quad.
void make_lines(std::vector<vec2i>& lines, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<float>& radius, const vec2i& steps, const vec2f& size,
    const vec2f& uvsize, const vec2f& line_radius = {0.001f, 0.001f});

// Make point primitives
void make_point(std::vector<int>& points, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<float>& radius, float point_radius = 0.001f);
void make_points(std::vector<int>& points, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<float>& radius, int num, float uvsize,
    float point_radius = 0.001f);
void make_random_points(std::vector<int>& points, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<float>& radius, int num, const vec3f& size, float uvsize,
    float point_radius = 0.001f, uint64_t seed = 0);

// Make a bezier circle. Returns bezier, pos.
void make_bezier_circle(std::vector<vec4i>& beziers, std::vector<vec3f>& pos);

// Parameters for the make hair function
// TODO: remove parameters from the struct -> move to function
struct make_hair_params {
    // minimum and maximum length
    vec2f length = {0.1f, 0.1f};
    // minimum and maximum radius from base to tip
    vec2f radius = {0.001f, 0.0001f};
    // noise added to hair (strength/scale)
    vec2f noise = zero2f;
    // clump added to hair (number/strength)
    vec2f clump = zero2f;
    // rotation
    vec2f rotation = zero2f;
    // random seed
    int seed = 0;
};

// Make a hair ball around a shape.
void make_hair(std::vector<vec2i>& lines, std::vector<vec3f>& pos,
    std::vector<vec3f>& tang, std::vector<vec2f>& texcoord,
    std::vector<float>& radius, const vec2i& steps,
    const std::vector<vec3i>& striangles, const std::vector<vec4i>& squads,
    const std::vector<vec3f>& spos, const std::vector<vec3f>& snorm,
    const std::vector<vec2f>& stexcoord, const make_hair_params& params);

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMAGE TYPE AND UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Ldr image of (r,g,b,a) pixels. Access pixels with at().
struct image4b {
    std::vector<vec4b> pixels;  // image pixels
    int width = 0;              // image width
    int height = 0;             // image height

    // pixel access
    vec4b& at(int i, int j) { return pixels.at(j * width + i); }
    const vec4b& at(int i, int j) const { return pixels.at(j * width + i); }
};

// Hdr image container. Access pixels with at() or operator [].
// TODO: remove container
struct image4f {
    std::vector<vec4f> pixels;  // image pixels
    int width = 0;              // image width
    int height = 0;             // image height

    // pixel access
    vec4f& at(int i, int j) { return pixels.at(j * width + i); }
    const vec4f& at(int i, int j) const { return pixels.at(j * width + i); }
};

// Initializes empty images
inline image4b make_image4b(
    int width, int height, const vec4b& val = {0, 0, 0, 0}) {
    auto img = image4b();
    img.width = width;
    img.height = height;
    img.pixels.assign(width * height, val);
    return img;
}
inline image4f make_image4f(
    int width, int height, const vec4f& val = {0, 0, 0, 0}) {
    auto img = image4f();
    img.width = width;
    img.height = height;
    img.pixels.assign(width * height, val);
    return img;
}

// Create an image with values stored in an array in scanline order.
inline image4b make_image4b(int width, int height, const vec4b* vals) {
    auto img = image4b();
    img.width = width;
    img.height = height;
    img.pixels.assign(vals, vals + width * height);
    return img;
}
inline image4f make_image4f(int width, int height, const vec4f* vals) {
    auto img = image4f();
    img.width = width;
    img.height = height;
    img.pixels.assign(vals, vals + width * height);
    return img;
}

// Create a 4 channel image with the given number of channels
inline image4f make_image4f(
    int w, int h, int nc, const float* vals, const vec4f& def) {
    auto img = make_image4f(w, h);
    for (auto j = 0; j < h; j++) {
        for (auto i = 0; i < w; i++) {
            auto pixel = vals + (j * w + i) * nc;
            img.at(i, j) = def;
            auto img_pixel = &img.at(i, j).x;
            for (auto c = 0; c < nc; c++) img_pixel[c] = pixel[c];
        }
    }
    return img;
}
inline image4b make_image4b(
    int w, int h, int nc, const byte* vals, const vec4b& def) {
    auto img = make_image4b(w, h);
    for (auto j = 0; j < h; j++) {
        for (auto i = 0; i < w; i++) {
            auto pixel = vals + (j * w + i) * nc;
            img.at(i, j) = def;
            auto img_pixel = &img.at(i, j).x;
            for (auto c = 0; c < nc; c++) img_pixel[c] = pixel[c];
        }
    }
    return img;
}

// Check if a pixel is inside an image.
inline bool contains(const image4b& img, int i, int j) {
    return i >= 0 && i < img.width && j >= 0 && j < img.height;
}
inline bool contains(const image4f& img, int i, int j) {
    return i >= 0 && i < img.width && j >= 0 && j < img.height;
}

// Approximate conversion from srgb.
inline vec4f srgb_to_linear(const vec4b& srgb) {
    return {pow(srgb.x / 255.0f, 2.2f), pow(srgb.y / 255.0f, 2.2f),
        pow(srgb.z / 255.0f, 2.2f), srgb.w / 255.0f};
}
// Approximate conversion to srgb.
inline vec4b linear_to_srgb(const vec4f& lin) {
    return float_to_byte({pow(lin.x, 1 / 2.2f), pow(lin.y, 1 / 2.2f),
        pow(lin.z, 1 / 2.2f), lin.w});
}

// Approximate conversion from srgb.
inline image4f srgb_to_linear(const image4b& srgb) {
    auto lin = make_image4f(srgb.width, srgb.height);
    for (auto j = 0; j < srgb.height; j++)
        for (auto i = 0; i < srgb.width; i++)
            lin.at(i, j) = srgb_to_linear(srgb.at(i, j));
    return lin;
}
// Approximate conversion to srgb.
inline image4b linear_to_srgb(const image4f& lin) {
    auto srgb = make_image4b(lin.width, lin.height);
    for (auto j = 0; j < srgb.height; j++)
        for (auto i = 0; i < srgb.width; i++)
            srgb.at(i, j) = linear_to_srgb(lin.at(i, j));
    return srgb;
}

// Convert between CIE XYZ and xyY
inline vec3f xyz_to_xyY(const vec3f& xyz) {
    if (xyz == zero3f) return zero3f;
    return {xyz.x / (xyz.x + xyz.y + xyz.z), xyz.y / (xyz.x + xyz.y + xyz.z),
        xyz.y};
}
// Convert between CIE XYZ and xyY
inline vec3f xyY_to_xyz(const vec3f& xyY) {
    if (xyY.y == 0) return zero3f;
    return {xyY.x * xyY.z / xyY.y, xyY.z, (1 - xyY.x - xyY.y) * xyY.z / xyY.y};
}
// Convert between CIE XYZ and RGB
inline vec3f xyz_to_rgb(const vec3f& xyz) {
    // from http://www.brucelindbloom.com/index.html?Eqn_RGB_to_XYZ.html
    if (xyz == zero3f) return zero3f;
    return {+3.2404542f * xyz.x - 1.5371385f * xyz.y - 0.4985314f * xyz.z,
        -0.9692660f * xyz.x + 1.8760108f * xyz.y + 0.0415560f * xyz.z,
        +0.0556434f * xyz.x - 0.2040259f * xyz.y + 1.0572252f * xyz.z};
}
// Convert between CIE XYZ and RGB
inline vec3f rgb_to_xyz(const vec3f& rgb) {
    // from http://www.brucelindbloom.com/index.html?Eqn_RGB_to_XYZ.html
    if (rgb == zero3f) return zero3f;
    return {0.4124564f * rgb.x + 0.3575761f * rgb.y + 0.1804375f * rgb.z,
        0.2126729f * rgb.x + 0.7151522f * rgb.y + 0.0721750f * rgb.z,
        0.0193339f * rgb.x + 0.1191920f * rgb.y + 0.9503041f * rgb.z};
}

// Image over operator.
void image_over(vec4f* img, int width, int height, int nlayers, vec4f** layers);

// Image over operator.
void image_over(vec4b* img, int width, int height, int nlayers, vec4b** layers);

// Converts HSV to RGB.
vec4b hsv_to_rgb(const vec4b& hsv);

// Tone mapping type.
enum struct tonemap_type {
    linear = 0,
    gamma = 1,
    srgb = 2,
    filmic1 = 3,
    filmic2 = 4,
    filmic3 = 5
};

// Tone mapping HDR to LDR images.
image4b tonemap_image(const image4f& hdr, tonemap_type type, float exposure);

// Names of enum values.
inline const std::map<tonemap_type, std::string>& tonemap_type_names() {
    static auto names = std::map<tonemap_type, std::string>{
        {tonemap_type::linear, "linear"},
        {tonemap_type::gamma, "gamma"},
        {tonemap_type::srgb, "srgb"},
        {tonemap_type::filmic1, "filmic1"},
        {tonemap_type::filmic2, "filmic2"},
        {tonemap_type::filmic3, "filmic3"},
    };
    return names;
}

// Make example images.
image4b make_grid_image(int width, int height, int tile = 64,
    const vec4b& c0 = {64, 64, 64, 255},
    const vec4b& c1 = {128, 128, 128, 255});
image4b make_checker_image(int width, int height, int tile = 64,
    const vec4b& c0 = {64, 64, 64, 255},
    const vec4b& c1 = {128, 128, 128, 255});
image4b make_bumpdimple_image(int width, int height, int tile = 64);
image4b make_ramp_image(
    int width, int height, const vec4b& c0, const vec4b& c1, bool srgb = false);
image4b make_gammaramp_image(int width, int height);
image4f make_gammaramp_imagef(int width, int height);
image4b make_uv_image(int width, int height);
image4b make_uvgrid_image(
    int width, int height, int tile = 64, bool colored = true);
image4b make_recuvgrid_image(
    int width, int height, int tile = 64, bool colored = true);

// Comvert a bump map to a normal map.
image4b bump_to_normal_map(const image4b& img, float scale = 1);

// Make a sunsky HDR model with sun at theta elevation in [0,pi/2], turbidity
// in [1.7,10] with or without sun.
image4f make_sunsky_image(int res, float thetaSun, float turbidity = 3,
    bool has_sun = false, bool has_ground = true);

// Make a noise image. Wrap works only if both resx and resy are powers of two.
image4b make_noise_image(int resx, int resy, float scale = 1, bool wrap = true);
image4b make_fbm_image(int resx, int resy, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, int octaves = 6, bool wrap = true);
image4b make_ridge_image(int resx, int resy, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, float offset = 1.0f,
    int octaves = 6, bool wrap = true);
image4b make_turbulence_image(int resx, int resy, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, int octaves = 6, bool wrap = true);

#if YGL_IMAGEIO

// Check if an image is HDR based on filename.
bool is_hdr_filename(const std::string& filename);

// Loads/saves a 4 channel ldr/hdr image.
image4b load_image4b(const std::string& filename);
image4f load_image4f(const std::string& filename);
bool save_image4b(const std::string& filename, const image4b& img);
bool save_image4f(const std::string& filename, const image4f& img);
// Save a 4 channel HDR or LDR image with tonemapping based on filename.
bool save_image(const std::string& filename, const image4f& hdr,
    tonemap_type tonemapper, float exposure);

// Loads.saves an image with variable number of channels.
std::vector<float> load_imagef(
    const std::string& filename, int& width, int& height, int& ncomp);
std::vector<byte> load_imageb(
    const std::string& filename, int& width, int& height, int& ncomp);
std::vector<float> load_imagef_from_memory(const std::string& filename,
    const byte* data, int length, int& width, int& height, int& ncomp);
std::vector<byte> load_imageb_from_memory(const std::string& filename,
    const byte* data, int length, int& width, int& height, int& ncomp);
bool save_imagef(const std::string& filename, int width, int height, int ncomp,
    const float* hdr);
bool save_imageb(const std::string& filename, int width, int height, int ncomp,
    const byte* ldr);

// Filter type and edge mode for resizing.
enum struct resize_filter {
    def,
    box,
    triangle,
    cubic_spline,
    catmull_rom,
    mitchell
};
enum struct resize_edge { def, clamp, reflect, wrap, zero };

// Resize an image.
void resize_image(const image4f& img, image4f& res_img,
    resize_filter filter = resize_filter::def,
    resize_edge edge = resize_edge::def, bool premultiplied_alpha = false);
void resize_image(const image4b& img, image4b& res_img,
    resize_filter filter = resize_filter::def,
    resize_edge edge = resize_edge::def, bool premultiplied_alpha = false);

#endif

}  // namespace ygl

// -----------------------------------------------------------------------------
// RAY INTERSECTION AND CLOSEST POINT FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Intersect a ray with a point (approximate).
// Based on http://geomalgorithms.com/a02-lines.html.
bool intersect_point(const ray3f& ray, const vec3f& p, float r, float& ray_t);

// Intersect a ray with a line (approximate).
// Based on http://geomalgorithms.com/a05-intersect-1.html and
// http://geomalgorithms.com/a07-distance.html#
//     dist3D_Segment_to_Segment
bool intersect_line(const ray3f& ray, const vec3f& v0, const vec3f& v1,
    float r0, float r1, float& ray_t, vec2f& euv);

// Intersect a ray with a triangle.
bool intersect_triangle(const ray3f& ray, const vec3f& v0, const vec3f& v1,
    const vec3f& v2, float& ray_t, vec2f& euv);

// Intersect a ray with a quad represented as two triangles (0,1,3) and
// (2,3,1), with the uv coordinates of the second triangle corrected by u =
// 1-u' and v = 1-v' to produce a quad parametrization where u and v go from 0
// to 1. This is equivalent to Intel's Embree.
bool intersect_quad(const ray3f& ray, const vec3f& v0, const vec3f& v1,
    const vec3f& v2, const vec3f& v3, float& ray_t, vec2f& euv);

// Intersect a ray with a axis-aligned bounding box.
bool intersect_bbox(const ray3f& ray, const bbox3f& bbox);

// Intersect a ray with a axis-aligned bounding box, implemented as
// "Robust BVH Ray Traversal" by T. Ize published at
// http://jcgt.org/published/0002/02/02/paper.pdf
bool intersect_bbox(const ray3f& ray, const vec3f& ray_dinv,
    const vec3i& ray_dsign, const bbox3f& bbox);

// Check if a point overlaps a position within a max distance.
bool overlap_point(
    const vec3f& pos, float dist_max, const vec3f& v0, float r0, float& dist);

// Find closest line point to a position.
float closestuv_line(const vec3f& pos, const vec3f& v0, const vec3f& v1);

// Check if a line overlaps a position within a max distance.
bool overlap_line(const vec3f& pos, float dist_max, const vec3f& v0,
    const vec3f& v1, float r0, float r1, float& dist, vec2f& euv);

// Find closest triangle point to a position.
vec2f closestuv_triangle(
    const vec3f& pos, const vec3f& v0, const vec3f& v1, const vec3f& v2);

// Check if a triangle overlaps a position within a max distance.
bool overlap_triangle(const vec3f& pos, float dist_max, const vec3f& v0,
    const vec3f& v1, const vec3f& v2, float r0, float r1, float r2, float& dist,
    vec2f& euv);

// Check if a quad overlaps a position within a max distance.
bool overlap_quad(const vec3f& pos, float dist_max, const vec3f& v0,
    const vec3f& v1, const vec3f& v2, const vec3f& v3, float r0, float r1,
    float r2, float r3, float& dist, vec2f& euv);

// Check if a bouning box overlaps a position within a max distance.
bool overlap_bbox(const vec3f& pos, float dist_max, const bbox3f& bbox);

// Check if two bouning boxes overlap.
bool overlap_bbox(const bbox3f& bbox1, const bbox3f& bbox2);

}  // namespace ygl

// -----------------------------------------------------------------------------
// BVH FOR RAY INTERSECTION AND CLOSEST ELEMENT
// -----------------------------------------------------------------------------
namespace ygl {

// Type of BVH node.
enum struct bvh_node_type : uint8_t {
    internal = 0,
    point,
    line,
    triangle,
    quad,
    vertex,
    instance,
};

// Maximum number of primitives per BVH node.
const int bvh_max_prims = 4;

// BVH tree node containing its bounds, indices to the BVH arrays of either
// primitives or internal nodes, the node element type,
// and the split axis. Leaf and internal nodes are identical, except that
// indices refer to primitives for leaf nodes or other nodes for internal
// nodes. See bvh_tree for more details.
// This is an internal data structure.
struct bvh_node {
    bbox3f bbox;                    // bouds
    uint32_t prims[bvh_max_prims];  // primitives
    uint16_t count;                 // number of prims
    bvh_node_type type;             // node type
    uint8_t split_axis;             // split axis
};

// Shape instance for two-level BVH.
// This is an internal data structure.
struct bvh_instance {
    frame3f frame = identity_frame3f;      // frame
    frame3f frame_inv = identity_frame3f;  // frame inverse
    int shape_id = 0;                      // shape id
};

// BVH tree, stored as a node array. The tree structure is encoded using array
// indices instead of pointers, both for speed but also to simplify code.
// BVH nodes indices refer to either the node array, for internal nodes,
// or the primitive arrays, for leaf nodes. BVH trees may contain only one type
// of geometric primitive, like points, lines, triangle or instances of other
// BVHs. To handle multiple primitive types and transformed primitives, build
// a two-level hierarchy with the outer BVH, the scene BVH, containing inner
// BVHs, shape BVHs, each of which of a uniform primitive type.
// This is an internal data structure.
struct bvh_tree {
    std::vector<bvh_node> nodes;                   // Internal nodes.
    bvh_node_type type = bvh_node_type::internal;  // Bvh leaF type

    // data for shape BVHs
    std::vector<vec3f> pos;        // Positions for shape BVHs.
    std::vector<float> radius;     // Radius for shape BVHs.
    std::vector<int> points;       // Points for shape BVHs.
    std::vector<vec2i> lines;      // Lines for shape BVHs.
    std::vector<vec3i> triangles;  // Triangles for shape BVHs.
    std::vector<vec4i> quads;      // Quads for shape BVHs.

    // data for instance BVHs
    std::vector<bvh_instance>
        instances;  // Instance ids (iid, sid, shape bvh index).
    std::vector<bvh_tree*> shape_bvhs;  // Shape BVHs.
    bool own_shape_bvhs = false;        // Owns shape BVHs.

    // Cleanup
    ~bvh_tree() {
        if (own_shape_bvhs)
            for (auto v : shape_bvhs) delete v;
    }
};

// Build a shape BVH from a set of primitives.
bvh_tree* make_bvh(const std::vector<int>& points,
    const std::vector<vec2i>& lines, const std::vector<vec3i>& triangles,
    const std::vector<vec4i>& quads, const std::vector<vec3f>& pos,
    const std::vector<float>& radius, float def_radius, bool equalsize);
// Build a scene BVH from a set of shape instances.
bvh_tree* make_bvh(const std::vector<bvh_instance>& instances,
    const std::vector<bvh_tree*>& shape_bvhs, bool equal_size,
    bool owns_shape_bvhs);

// Grab the shape BVHs
inline const std::vector<bvh_tree*>& get_shape_bvhs(const bvh_tree* bvh) {
    return bvh->shape_bvhs;
}

// Update the node bounds for a shape bvh.
void refit_bvh(bvh_tree* bvh, const std::vector<vec3f>& pos,
    const std::vector<float>& radius, float def_radius);
// Update the node bounds for a scene bvh
void refit_bvh(bvh_tree* bvh, const std::vector<frame3f>& frames,
    const std::vector<frame3f>& frames_inv);

// Intersect ray with a bvh returning either the first or any intersection
// depending on `find_any`. Returns the ray distance `ray_t`, the instance
// id `iid`, the shape id `sid`, the shape element index `eid` and the
// shape barycentric coordinates `euv`.
bool intersect_bvh(const bvh_tree* bvh, const ray3f& ray, bool find_any,
    float& ray_t, int& iid, int& eid, vec2f& euv);

// Find a shape element that overlaps a point within a given distance
// `max_dist`, returning either the closest or any overlap depending on
// `find_any`. Returns the point distance `dist`, the instance id `iid`, the
// shape id `sid`, the shape element index `eid` and the shape barycentric
// coordinates `euv`.
bool overlap_bvh(const bvh_tree* bvh, const vec3f& pos, float max_dist,
    bool find_any, float& dist, int& iid, int& eid, vec2f& euv);

// Intersection point.
struct intersection_point {
    float dist = 0;      // Distance of the hit.
    int iid = -1;        // Instance index.
    int eid = -1;        // Shape element index.
    vec2f euv = zero2f;  // Shape barycentric coordinates.
    operator bool() const { return eid >= 0; }  // Check if valid
};

// Intersect a ray with a bvh (convenience wrapper).
intersection_point intersect_bvh(
    const bvh_tree* bvh, const ray3f& ray, bool early_exit);
// Finds the closest element with a bvh (convenience wrapper).
intersection_point overlap_bvh(
    const bvh_tree* bvh, const vec3f& pos, float max_dist, bool early_exit);

}  // namespace ygl

// -----------------------------------------------------------------------------
// SIMPLE SCENE SUPPORT
// -----------------------------------------------------------------------------
namespace ygl {

// #codegen begin refl-scene

// Camera.
struct camera {
    std::string name = "";             // name
    frame3f frame = identity_frame3f;  // transform frame
    bool ortho = false;                // orthographic
    float yfov = 2;                    // vertical field of view.
    float aspect = 16.0f / 9.0f;       // aspect ratio
    float focus = 1;                   // focu distance
    float aperture = 0;                // lens aperture
    float near = 0.01f;                // near plane distance
    float far = 10000;                 // far plane distance
};

// Texture containing either an LDR or HDR image.
struct texture {
    std::string name = "";  // name
    std::string path = "";  // file path
    image4b ldr = {};       // ldr image
    image4f hdr = {};       // hdr image
};

// Texture information to use for lookup.
struct texture_info {
    texture* txt = nullptr;  // texture
    bool wrap_s = true;      // wrap in s coordinate
    bool wrap_t = true;      // wrop in t coordinate
    bool linear = true;      // linear interpolation
    bool mipmap = true;      // mipmapping
    float scale = 1;         // scale for occ, normal, bumps
};

// Material type.
enum struct material_type {
    specular_roughness = 0,   // microfacet model (OBJ)
    metallic_roughness = 1,   // base/metallic model (glTF)
    specular_glossiness = 2,  // sepcular/glossiness (glTF)
};

// Material for surfaces, lines and triangles.
// For surfaces, uses a microfacet model with thin sheet transmission.
// The model is based on glTF for compatibility and adapted to OBJ.
// For lines, uses Kajija-Kay model. For points, a hacked up shading.
struct material {
    std::string name = "";      // name
    bool double_sided = false;  // double-sided rendering
    material_type type = material_type::specular_roughness;  // type

    // base values
    vec3f ke = {0, 0, 0};  // emission color
    vec3f kd = {0, 0, 0};  // diffuse/base color
    vec3f ks = {0, 0, 0};  // specular color / metallic factor
    vec3f kr = {0, 0, 0};  // clear coat reflection
    vec3f kt = {0, 0, 0};  // transmission color
    float rs = 0.0001;     // roughness mapped as glTF
    float op = 1;          // opacity

    // textures
    texture_info ke_txt;    // emission texture
    texture_info kd_txt;    // diffuse texture
    texture_info ks_txt;    // specular texture
    texture_info kr_txt;    // clear coat reflection texture
    texture_info kt_txt;    // transmission texture
    texture_info rs_txt;    // roughness texture
    texture_info op_txt;    // opacity texture
    texture_info bump_txt;  // bump map texture (heighfield)
    texture_info disp_txt;  // displacement map texture (heighfield)
    texture_info norm_txt;  // normal texture
    texture_info occ_txt;   // occlusion texture
};

// Shape data represented as an indexed meshes of elements.
// May contain only element type (points/lines/triangles/quads/beziers).
struct shape {
    std::string name = "";  // name
    std::string path = "";  // path for glTF buffers

    // primitives
    std::vector<int> points;       // points
    std::vector<vec2i> lines;      // lines
    std::vector<vec3i> triangles;  // triangles
    std::vector<vec4i> quads;      // quads
    std::vector<vec4i> beziers;    // beziers

    // face-varying quad primitives
    std::vector<vec4i> quads_pos;       // pos indices
    std::vector<vec4i> quads_norm;      // norm indices
    std::vector<vec4i> quads_texcoord;  // texcoord indices

    // vertex data
    std::vector<vec3f> pos;        // positions
    std::vector<vec3f> norm;       // normals/tangents
    std::vector<vec2f> texcoord;   // texcoord coordinates
    std::vector<vec2f> texcoord1;  // second set of texture coordinates
    std::vector<vec4f> color;      // colors
    std::vector<float> radius;     // radia for lines/points
    std::vector<vec4f> tangsp;     // tangent space for triangles

    // tesselation data
    int subdivision = 0;        // subdivision [deprecated]
    bool catmullclark = false;  // catmull-clark [deprecated]
};

// Shape instance.
struct instance {
    std::string name;                  // name
    frame3f frame = identity_frame3f;  // transform frame
    shape* shp = nullptr;              // shape
    material* mat = nullptr;           // material
};

// Distance at which we set environment map positions.
const auto environment_distance = 1000000.0f;

// Envinonment map.
struct environment {
    std::string name = "";             // name
    frame3f frame = identity_frame3f;  // transform frame
    vec3f ke = {0, 0, 0};              // emission color
    texture_info ke_txt = {};          // emission texture
};

// Node in a transform hierarchy.
struct node {
    std::string name = "";             // name
    node* parent = nullptr;            // parent
    frame3f frame = identity_frame3f;  // transform frame
    vec3f translation = {0, 0, 0};     // translation
    vec4f rotation = {0, 0, 0, 1};     // rotation
    vec3f scale = {1, 1, 1};           // scale
    std::vector<float> weights = {};   // morph weights
    camera* cam = nullptr;             // camera
    instance* ist = nullptr;           // instance
    environment* env = nullptr;        // environment

    // compute properties
    std::vector<node*> children_ = {};  // child nodes
};

// Keyframe type.
enum struct animation_type { linear, step, bezier };

// Keyframe data.
struct animation {
    std::string name;                              // name
    std::string path = "";                         // path for glTF buffer
    std::string group;                             // group
    animation_type type = animation_type::linear;  // type
    std::vector<float> times;                      // keyframe times
    std::vector<vec3f> translation;                // translation keyframes
    std::vector<vec4f> rotation;                   // rotation keyframes
    std::vector<vec3f> scale;                      // scale keyframes
    std::vector<std::vector<float>> weights;       // mprph weight keyframes
    std::vector<node*> targets;                    // target nodes
};

// Scene comprised an array of objects whose memory is owened by the scene.
// All members are optional,Scene objects (camera, instances, environments)
// have transforms defined internally. A scene can optionally contain a
// node hierarchy where each node might point to a camera, instance or
// environment. In that case, the element transforms are computed from
// the hierarchy. Animation is also optional, with keyframe data that
// updates node transformations only if defined.
struct scene {
    std::vector<camera*> cameras = {};            // cameras
    std::vector<shape*> shapes = {};              // shapes
    std::vector<instance*> instances = {};        // instances
    std::vector<material*> materials = {};        // materials
    std::vector<texture*> textures = {};          // textures
    std::vector<environment*> environments = {};  // environments

    std::vector<node*> nodes = {};            // node hierarchy [optional]
    std::vector<animation*> animations = {};  // animations [optional]

    // Cleanup
    ~scene() {
        for (auto v : shapes) delete v;
        for (auto v : instances) delete v;
        for (auto v : materials) delete v;
        for (auto v : textures) delete v;
        for (auto v : cameras) delete v;
        for (auto v : environments) delete v;
        for (auto v : nodes) delete v;
        for (auto v : animations) delete v;
    }
};

// Shape elements type.
enum struct shape_elem_type {
    none,
    points,
    lines,
    triangles,
    quads,
    beziers,
    vertices,
    facevarying
};

// Get shape element type.
shape_elem_type get_shape_type(const shape* shp);
// Shape element normal.
vec3f eval_elem_norm(const shape* shp, int eid);

// Shape values interpolated using barycentric coordinates.
vec3f eval_pos(const shape* shp, int eid, const vec2f& euv);
vec3f eval_norm(const shape* shp, int eid, const vec2f& euv);
vec2f eval_texcoord(const shape* shp, int eid, const vec2f& euv);
vec4f eval_color(const shape* shp, int eid, const vec2f& euv);
float eval_radius(const shape* shp, int eid, const vec2f& euv);
vec4f eval_tangsp(const shape* shp, int eid, const vec2f& euv);

// Environment values interpolated using uv parametrization.
vec3f eval_pos(
    const environment* env, const vec2f& uv, bool transformed = false);
vec3f eval_norm(
    const environment* env, const vec2f& uv, bool transformed = false);
// Environment texture coordinates from uv parametrization.
vec2f eval_texcoord(const environment* env, const vec2f& uv);
// Evaluate uv parameters for environment.
vec2f eval_uv(const environment* env, const vec3f& w, bool transformed = false);

// Evaluate a texture.
vec4f eval_texture(const texture_info& info, const vec2f& texcoord,
    bool srgb = true, const vec4f& def = {1, 1, 1, 1});
inline vec4f eval_texture(const texture* txt, const vec2f& texcoord,
    bool srgb = true, const vec4f& def = {1, 1, 1, 1}) {
    auto info = texture_info();
    info.txt = (texture*)txt;
    return eval_texture(info, texcoord, srgb, def);
}
// Generates a ray from a camera image coordinate `uv` and lens coordinates
// `luv`.
ray3f eval_camera_ray(const camera* cam, const vec2f& uv, const vec2f& luv);
// Generates a ray from a camera for pixel coordinates `ij`, the resolution
// `res`, the sub-pixel coordinates `puv` and the lens coordinates `luv` and
// the image resolution `res`.
ray3f eval_camera_ray(const camera* cam, const vec2i& ij, int res,
    const vec2f& puv, const vec2f& luv);
// Synchronizes a camera aspect with image width and height. Set image
// values any one is 0 or less. Set camera aspect otherwise.
void sync_camera_aspect(const camera* cam, int& width, int& height);

// Generate a distribution for sampling a shape uniformly based on area/length.
std::vector<float> sample_shape_cdf(const shape* shp);
// Sample a shape based on a distribution.
std::pair<int, vec2f> sample_shape(const shape* shp,
    const std::vector<float>& cdf, float re, const vec2f& ruv);

// Sample an environment uniformly.
vec2f sample_environment(const environment* env, const vec2f& ruv);

// Update the normals of a shape.  Supports only non-facevarying shapes.
void compute_normals(shape* shp);
// Subdivides shape elements. Apply subdivision surface rules if subdivide
// is true.
void subdivide_shape_once(shape* shp, bool subdiv = false);
// Tesselate a shape into basic primitives.
void tesselate_shape(shape* shp, bool subdivide,
    bool facevarying_to_sharedvertex, bool quads_to_triangles,
    bool bezier_to_lines);
// Tesselate scene shapes.
void tesselate_shapes(scene* scn, bool subdivide,
    bool facevarying_to_sharedvertex, bool quads_to_triangles,
    bool bezier_to_lines);

// Update node transforms.
void update_transforms(
    scene* scn, float time = 0, const std::string& anim_group = "");
// Compute animation range.
vec2f compute_animation_range(
    const scene* scn, const std::string& anim_group = "");

// Make a view camera either copying a given one or building a default one.
camera* make_view_camera(const scene* scn, int camera_id);

// Computes shape/scene approximate bounds.
bbox3f compute_bbox(const shape* shp);
bbox3f compute_bbox(const scene* scn, bool skip_emitting = false);

// Build a shape/scene BVH.
bvh_tree* make_bvh(
    const shape* shp, float def_radius = 0.001f, bool equalsize = true);
bvh_tree* make_bvh(
    const scene* scn, float def_radius = 0.001f, bool equalsize = true);

// Refits shape/scene BVH.
void refit_bvh(bvh_tree* bvh, const shape* shp, float def_radius = 0.001f);
void refit_bvh(
    bvh_tree* bvh, const scene* scn, bool do_shapes, float def_radius = 0.001f);

// Add missing names and resolve duplicated names.
void add_names(scene* scn);

// Add missing normals.
void add_normals(scene* scn);

// Add missing tangent space if needed.
void add_tangent_space(scene* scn);

// Add hierarchy.
void add_hierarchy(scene* scn);

// Checks for validity of the scene.
std::vector<std::string> validate(
    const scene* scn, bool skip_missing = false, bool log_as_warning = false);

// Merge scene into one another. Note that the objects are _moved_ from
// merge_from to merged_into, so merge_from will be empty after this function.
void merge_into(scene* merge_into, scene* merge_from);

// #codegen begin refl-scene

// Loading options.
struct load_options {
    bool load_textures = true;              // load textures
    bool skip_missing = true;               // whether to skip errors
    bool obj_flip_texcoord = true;          // flip texcoord in obj
    bool obj_flip_tr = true;                // flip tr in obj
    bool obj_preserve_quads = false;        // preserve quads
    bool obj_preserve_facevarying = false;  // preserve facevarying
    bool obj_split_shapes = false;          // split shapes in obj
};

// Save options.
struct save_options {
    bool save_textures = true;           // save textures
    bool skip_missing = true;            // whether to skip errors
    bool obj_flip_texcoord = true;       // flip obj texcoords
    bool obj_flip_tr = true;             // flip tr in obj
    bool obj_save_instances = false;     // preserve instances in obj
    bool gltf_separate_buffers = false;  // use separate buffers in glTF
};

// Loads/sa ves a scene in OBJ and glTF formats.
scene* load_scene(const std::string& filename, const load_options& opts = {});
void save_scene(
    const std::string& filename, const scene* scn, const save_options& opts);

// Scene selection.
struct scene_selection {
    // initialize selection
    scene_selection() : ptr(nullptr), tinfo(nullptr) {}
    template <typename T>
    scene_selection(T* val) : ptr(val), tinfo(&typeid(T)) {}

    // Checking whether it is empty.
    operator bool() const { return (bool)ptr; }
    bool empty() const { return (bool)ptr; }

    // Check if points to type T.
    template <typename T>
    bool is() const {
        return &typeid(T) == tinfo;
    }
    // Gets a pointer cast to the specific type.
    template <typename T>
    T* get() {
        return (is<T>()) ? (T*)ptr : nullptr;
    }
    // Get the raw untyped pointer.
    void* get_raw() { return ptr; }
    // Get untyped.
    void* get_untyped() { return ptr; }

   private:
    void* ptr = nullptr;                    // selected pointer
    const std::type_info* tinfo = nullptr;  // type info
};

// Print scene statistics.
void print_stats(const scene* scn);

// Names of enum values.
inline const std::map<material_type, std::string>& material_type_names() {
    static auto names = std::map<material_type, std::string>{
        {material_type::specular_roughness, "specular_roughness"},
        {material_type::metallic_roughness, "metallic_roughness"},
        {material_type::specular_glossiness, "specular_glossiness"},
    };
    return names;
}
inline const std::map<animation_type, std::string>& animation_type_names() {
    static auto names = std::map<animation_type, std::string>{
        {animation_type::linear, "linear"},
        {animation_type::step, "step"},
        {animation_type::bezier, "bezier"},
    };
    return names;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// EXAMPLE SCENES
// -----------------------------------------------------------------------------
namespace ygl {

// make scene elements
camera* make_camera(const std::string& name, const vec3f& from, const vec3f& to,
    float yfov, float aspect);
texture* make_texture(const std::string& name, const std::string& path = "",
    const image4b& ldr = {}, const image4f& hdr = {});
material* make_material(const std::string& name,
    const vec3f& kd = {0.2, 0.2, 0.2}, const vec3f& ks = {0, 0, 0},
    float rs = 1);
instance* make_instance(const std::string& name, shape* shp = nullptr,
    material* mat = nullptr, const frame3f& frame = identity_frame3f);
node* make_node(const std::string& name, camera* cam = nullptr,
    instance* ist = nullptr, environment* env = nullptr,
    const frame3f& frame = identity_frame3f);
environment* make_environment(const std::string& name,
    const vec3f& ke = {1, 1, 1}, texture* ke_txt = nullptr,
    const frame3f& frame = identity_frame3f);

// make procedural scene elements
inline std::vector<std::string>& make_camera_types() {
    static auto names = std::vector<std::string>{"cam1", "cam2", "cam3"};
    return names;
}
camera* make_proc_camera(const std::string& name, const std::string& type);
inline std::vector<std::string>& proc_texture_types() {
    static auto names = std::vector<std::string>{"grid", "checker", "colored",
        "rcolored", "bump", "uv", "gamma", "noise", "ridge", "fbm",
        "turbulence", "grid_norm", "bump_norm", "gammaf", "sky"};
    return names;
}
texture* make_proc_texture(const std::string& name, const std::string& type,
    int resolution = 512, float scale = 8.0f, float sky_sunangle = pi / 4,
    float bump_scale = 4);
inline std::vector<std::string>& proc_material_types() {
    static auto names = std::vector<std::string>{"emission", "matte", "plastic",
        "metal", "glass", "transparent", "carpaint"};
    return names;
}
material* make_proc_material(const std::string& name, const std::string& type,
    const vec3f& color = {1, 1, 1}, float roughness = 1,
    const std::string& txt = "", const std::string& norm = "");
inline std::vector<std::string>& proc_shape_types() {
    static auto names = std::vector<std::string>{"floor", "quad", "cube",
        "cube_rounded", "sphere", "sphere_cube", "geodesic_sphere",
        "sphere_flipcap", "disk", "disk_quad", "disk_bulged", "cylinder",
        "cylinder_rounded", "cylindery", "cylindery_rounded", "suzanne",
        "cube_subdiv", "suzanne_subdiv", "fvcube_subdiv", "matball", "point",
        "pointscube", "hairball", "beziercircle"};
    return names;
}
shape* make_proc_shape(const std::string& name, const std::string& type,
    const vec3i& tesselation = zero3i, const vec3f& size = zero3f,
    const vec3f& uvsize = zero3f, float rounded = 0.75f, float radius = 0.001f,
    const make_hair_params& hair_params = {});
instance* make_proc_instance(const std::string& name, const std::string& stype,
    const std::string& mtype, const frame3f& frame = identity_frame3f);

// Makes the Cornell Box scene.
scene* make_cornell_box_scene();
scene* make_simple_scene(const std::vector<std::string>& shapes,
    const std::vector<std::string>& mats, const std::string& lights,
    bool nodes = false, const std::vector<std::string>& animations = {},
    const std::string& floor_mat = "matte_grid");
scene* make_random_scene(
    const vec2i& num, const bbox3f& bbox, uint64_t seed = 13);

}  // namespace ygl

// -----------------------------------------------------------------------------
// PATH TRACING SUPPORT FUNCTION
// -----------------------------------------------------------------------------
namespace ygl {

// Phong exponent to roughness.
float specular_exponent_to_roughness(float n);

// Specular to fresnel eta.
void specular_fresnel_from_ks(const vec3f& ks, vec3f& es, vec3f& esk);
// Compute the fresnel term for dielectrics.
vec3f fresnel_dielectric(float cosw, const vec3f& eta_);
// Compute the fresnel term for metals.
vec3f fresnel_metal(float cosw, const vec3f& eta, const vec3f& etak);
// Schlick approximation of Fresnel term.
vec3f fresnel_schlick(const vec3f& ks, float cosw);
// Schlick approximation of Fresnel term weighted by roughness.
vec3f fresnel_schlick(const vec3f& ks, float cosw, float rs);

// Evaluates the GGX distribution and geometric term.
float eval_ggx(float rs, float ndh, float ndi, float ndo);
// Sample the GGX distribution.
vec3f sample_ggx(float rs, const vec2f& rn);
// Evaluates the GGX pdf.
float sample_ggx_pdf(float rs, float ndh);

}  // namespace ygl

// -----------------------------------------------------------------------------
// PATH TRACING
// -----------------------------------------------------------------------------
namespace ygl {

// #codegen begin refl-trace

// Type of rendering algorithm.
enum struct trace_type {
    pathtrace = 0,
    eyelight,
    direct,
    pathtrace_nomis,
    debug_normal,
    debug_albedo,
    debug_texcoord,
    debug_frontfacing
};

// Rendering params.
struct trace_params {
    int resolution = 512;                       // image vertical resolution
    int nsamples = 256;                         // number of samples
    trace_type tracer = trace_type::pathtrace;  // trace type
    bool notransmission = false;    // whether to test transmission in shadows
    bool double_sided = false;      // force double sided rendering
    vec3f ambient = {0, 0, 0};      // ambient lighting
    bool envmap_invisible = false;  // view environment map
    int min_depth = 3;              // minimum ray depth
    int max_depth = 8;              // maximum ray depth
    float pixel_clamp = 100;        // final pixel clamping
    float ray_eps = 1e-4f;          // ray intersection epsilon
    bool parallel = true;           // parallel execution
    int seed = 0;                   // seed for the random number generators
    int preview_resolution = 64;    // preview resolution for async rendering
    int batch_size = 16;            // sample batch size
};

// #codegen end refl-trace

// Trace pixel state. Handles image accumulation and random number generation
// for uniform and stratified sequences. The members are not part of the
// the public API.
struct trace_pixel {
    // Accumulated radiance and coverage
    vec4f acc = zero4f;
    // Random number state
    rng_state rng = rng_state();
    // Pixel coordinates
    int i = 0, j = 0;
    // Number of samples computed
    int sample = 0;
};

// Trace light as either an instance or an environment.
struct trace_light {
    const instance* ist = nullptr;     // instance for the light
    const environment* env = nullptr;  // environment for the light
};

// Trace lights. Handles sampling of illumination.
struct trace_lights {
    std::vector<trace_light> lights;  // lights
    std::unordered_map<const shape*, std::vector<float>>
        shape_distribs;                // shape dist
    std::vector<float> light_distrib;  // light distribution
};

// Initialize trace pixels.
std::vector<trace_pixel> make_trace_pixels(
    const image4f& img, const trace_params& params);
// Initialize trace lights.
trace_lights make_trace_lights(const scene* scn);

// Trace the next `nsamples` samples.
void trace_samples(const scene* scn, const camera* cam, const bvh_tree* bvh,
    const trace_lights& lights, image4f& img, std::vector<trace_pixel>& pixels,
    int nsamples, const trace_params& params);

// Trace the next `nsamples` samples with image filtering.
void trace_samples_filtered(const scene* scn, const camera* cam,
    const bvh_tree* bvh, const trace_lights& lights, image4f& img,
    std::vector<trace_pixel>& pixels, int nsamples, const trace_params& params);

// Trace the whole image.
inline void trace_image(const scene* scn, const camera* cam,
    const bvh_tree* bvh, const trace_lights& lights, image4f& img,
    std::vector<trace_pixel>& pixels, const trace_params& params,
    const std::function<void(int)>& callback) {
    for (auto& p : img.pixels) p = zero4f;
    for (auto cur_sample = 0; cur_sample < params.nsamples;
         cur_sample += params.batch_size) {
        if (callback) callback(cur_sample);
        trace_samples(scn, cam, bvh, lights, img, pixels,
            std::min(params.batch_size, params.nsamples - cur_sample), params);
    }
}

// Trace the whole image.
inline image4f trace_image(const scene* scn, const camera* cam,
    const bvh_tree* bvh, const trace_params& params,
    const std::function<void(int)>& callback) {
    auto img = make_image4f(
        (int)std::round(cam->aspect * params.resolution), params.resolution);
    auto pixels = make_trace_pixels(img, params);
    auto lights = make_trace_lights(scn);
    trace_image(scn, cam, bvh, lights, img, pixels, params, callback);
    return img;
}

// Starts an anyncrhounous renderer.
void trace_async_start(const scene* scn, const camera* cam, const bvh_tree* bvh,
    const trace_lights& lights, image4f& img, std::vector<trace_pixel>& pixels,
    std::vector<std::thread>& threads, bool& stop_flag,
    const trace_params& params, const std::function<void(int, int)>& callback);
// Stop the asynchronous renderer.
void trace_async_stop(std::vector<std::thread>& threads, bool& stop_flag);

// Names of enum values.
inline const std::map<trace_type, std::string>& trace_type_names() {
    static auto names = std::map<trace_type, std::string>{
        {trace_type::pathtrace, "pathtrace"},
        {trace_type::eyelight, "eyelight"},
        {trace_type::direct, "direct"},
        {trace_type::pathtrace_nomis, "pathtrace_nomis"},
        {trace_type::debug_normal, "debug_normal"},
        {trace_type::debug_albedo, "debug_albedo"},
        {trace_type::debug_texcoord, "debug_texcoord"},
        {trace_type::debug_frontfacing, "debug_frontfacing"},
    };
    return names;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// WAVEFRONT OBJ SUPPORT
// -----------------------------------------------------------------------------
namespace ygl {

// Obj face vertex.
struct obj_vertex {
    int pos = -1;       // position index
    int texcoord = -1;  // texcoord index
    int norm = -1;      // normal index
    int color = -1;     // color index [extension]
    int radius = -1;    // radius index [extension]
};

// Comparison for unordred_map.
inline bool operator==(const obj_vertex& a, const obj_vertex& b) {
    return a.pos == b.pos && a.texcoord == b.texcoord && a.norm == b.norm &&
           a.color == b.color && a.radius == b.radius;
}

// Obj element type.
enum struct obj_element_type : uint8_t {
    point = 1,
    line = 2,
    face = 3,
    bezier = 4
};

// Obj element (point/polyline/polygon)
struct obj_element {
    uint32_t start;         // starting vertex index
    uint8_t size;           // number of vertices
    obj_element_type type;  // element type
    uint16_t groupid = 0;   // group id
};

// Obj group properties.
struct obj_group_props {
    std::string name = "";     // group name
    std::string matname = "";  // material name
    bool faceted = false;      // faceted or smooth
};

// Obj object.
struct obj_object {
    std::string name;                     // name
    std::vector<obj_group_props> groups;  // groups
    std::vector<obj_vertex> verts;        // vertices
    std::vector<obj_element> elems;       // faces
    // Properties not explicitly handled [extension].
    std::unordered_map<std::string, std::vector<std::string>> props;
};

// Obj texture information.
struct obj_texture_info {
    std::string path = "";  // file path
    bool clamp = false;     // clamp to edge
    float scale = 1;        // scale for bump/displacement
    // Properties not explicitly handled.
    std::unordered_map<std::string, std::vector<std::string>> props;
};

// Comparison for texture info.
inline bool operator==(const obj_texture_info& a, const obj_texture_info& b) {
    if (a.path.empty() && b.path.empty()) return true;
    if (a.path != b.path) return false;
    return a.clamp == b.clamp && a.scale == b.scale && a.props == b.props;
}

// Obj texture data.
struct obj_texture {
    std::string path;            // file path
    int width = 0;               // width
    int height = 0;              // height
    int ncomp = 0;               // number of components [1-4]
    std::vector<uint8_t> datab;  // data for LDRs
    std::vector<float> dataf;    // data for HDRs
};

// Obj material.
struct obj_material {
    std::string name;  // name
    int illum = 0;     // MTL illum mode

    vec3f ke = {0, 0, 0};  // emission color
    vec3f ka = {0, 0, 0};  // ambient color
    vec3f kd = {0, 0, 0};  // diffuse color
    vec3f ks = {0, 0, 0};  // specular color
    vec3f kr = {0, 0, 0};  // reflection color
    vec3f kt = {0, 0, 0};  // transmission color
    float ns = 0;          // Phong exponent color
    float ior = 1;         // index of refraction
    float op = 1;          // opacity

    obj_texture_info ke_txt;    // emission texture
    obj_texture_info ka_txt;    // ambient texture
    obj_texture_info kd_txt;    // diffuse texture
    obj_texture_info ks_txt;    // specular texture
    obj_texture_info kr_txt;    // reflection texture
    obj_texture_info kt_txt;    // transmission texture
    obj_texture_info ns_txt;    // Phong exponent texture
    obj_texture_info op_txt;    // opacity texture
    obj_texture_info ior_txt;   // ior texture
    obj_texture_info bump_txt;  // bump map
    obj_texture_info disp_txt;  // displacement map
    obj_texture_info norm_txt;  // normal map

    // Properties not explicitly handled.
    std::unordered_map<std::string, std::vector<std::string>> props;
};

// Obj camera [extension].
struct obj_camera {
    std::string name;                  // name
    frame3f frame = identity_frame3f;  // transform
    bool ortho = false;                // orthographic
    float yfov = 2 * atan(0.5f);       // vertical field of view
    float aspect = 16.0f / 9.0f;       // aspect ratio
    float aperture = 0;                // lens aperture
    float focus = 1;                   // lens focus
};

// Obj environment [extension].
struct obj_environment {
    std::string name;                  // name
    frame3f frame = identity_frame3f;  // transform
    vec3f ke = zero3f;                 // emission color
    obj_texture_info ke_txt;           // emission texture
};

// Obj node [extension].
struct obj_node {
    std::string name;                  // name
    std::string parent;                // parent node
    std::string camname;               // camera name
    std::string objname;               // object name
    std::string envname;               // environment name
    frame3f frame = identity_frame3f;  // transform
    vec3f translation = zero3f;        // translation
    vec4f rotation = {0, 0, 0, 1};     // rotation
    vec3f scale = {1, 1, 1};           // scale
};

// Obj scene.
struct obj_scene {
    std::vector<vec3f> pos;       // vertex positions
    std::vector<vec3f> norm;      // vertex normals
    std::vector<vec2f> texcoord;  // vertex texcoords
    std::vector<vec4f> color;     // vertex colors [extension]
    std::vector<float> radius;    // vertex radia [extension]

    std::vector<obj_object*> objects;            // objects
    std::vector<obj_material*> materials;        // materials
    std::vector<obj_texture*> textures;          // textures
    std::vector<obj_camera*> cameras;            // cameras [extension]
    std::vector<obj_environment*> environments;  // environments [extension]
    std::vector<obj_node*> nodes;                // nodes [extension]

    // Cleanup.
    ~obj_scene() {
        for (auto v : objects) delete v;
        for (auto v : materials) delete v;
        for (auto v : textures) delete v;
        for (auto v : cameras) delete v;
        for (auto v : environments) delete v;
        for (auto v : nodes) delete v;
    }
};

// Load an OBJ from file `filename`. Split shapes at material and group
// boundaries, if `split_shapes` is true.
// Load textures if `load_textures` is true, and report errors only if
// `skip_missing` is false. Texture coordinates and material Tr are flipped
// if `flip_texcoord` and `flip_tp` are respectively true.
obj_scene* load_obj(const std::string& filename, bool split_shapes,
    bool load_textures = false, bool skip_missing = false,
    bool flip_texcoord = true, bool flip_tr = true);

// Save an OBJ to file `filename`. Save textures if `save_textures` is true,
// and report errors only if `skip_missing` is false.
// Texture coordinates and material Tr are flipped if `flip_texcoord` and
// `flip_tp` are respectively true.
void save_obj(const std::string& filename, const obj_scene* model,
    bool save_textures = false, bool skip_missing = false,
    bool flip_texcoord = true, bool flip_tr = true);

}  // namespace ygl

#if YGL_GLTF

// include json for glTF
#if YGL_GLTFJSON
#include "ext/json.hpp"
#endif

// -----------------------------------------------------------------------------
// KHRONOS GLTF SUPPORT
// -----------------------------------------------------------------------------
namespace ygl {

// Generic buffer data.
using buffer_data = std::vector<unsigned char>;

// Generic image data.
struct image_data {
    // Width.
    int width = 0;
    // Height.
    int height = 0;
    // Number of Channels.
    int ncomp = 0;
    // Buffer data for 8-bit images.
    std::vector<uint8_t> datab;
    // Buffer data for float images.
    std::vector<float> dataf;
};

// Id for glTF references.
template <typename T>
struct glTFid {
    // Defaoult constructor to an invalid id.
    glTFid() : _id(-1) {}
    // Explicit conversion from integer.
    explicit glTFid(int id) : _id(id) {}
    // Explicit convcersion to integer.
    explicit operator int() const { return _id; }
    // Check if it is valid.
    bool is_valid() const { return _id >= 0; }
    // Check if it is valid.
    explicit operator bool() const { return _id >= 0; }

   private:
    // id
    int _id = -1;
};

// Generic glTF object.
struct glTFProperty {
#if YGL_GLTFJSON
    // Extensions.
    map<string, nlohmann::json> extensions = {};
    // Extra data.
    nlohmann::json extras = {};
#endif
};

// #codegen begin gltf-type

// forward declaration
struct glTFChildOfRootProperty;
struct glTFAccessorSparseIndices;
struct glTFAccessorSparseValues;
struct glTFAccessorSparse;
struct glTFAccessor;
struct glTFAnimationChannelTarget;
struct glTFAnimationChannel;
struct glTFAnimationSampler;
struct glTFAnimation;
struct glTFAsset;
struct glTFBuffer;
struct glTFBufferView;
struct glTFCameraOrthographic;
struct glTFCameraPerspective;
struct glTFCamera;
struct glTFImage;
struct glTFTextureInfo;
struct glTFTexture;
struct glTFMaterialNormalTextureInfo;
struct glTFMaterialOcclusionTextureInfo;
struct glTFMaterialPbrMetallicRoughness;
struct glTFMaterialPbrSpecularGlossiness;
struct glTFMaterial;
struct glTFMeshPrimitive;
struct glTFMesh;
struct glTFNode;
struct glTFSampler;
struct glTFScene;
struct glTFSkin;
struct glTF;

// Generic glTF named object
struct glTFChildOfRootProperty : glTFProperty {
    // The user-defined name of this object.
    std::string name = "";
};

// Values for glTFAccessorSparseIndices::componentType
enum class glTFAccessorSparseIndicesComponentType {
    // Not set
    NotSet = -1,
    // UnsignedByte
    UnsignedByte = 5121,
    // UnsignedShort
    UnsignedShort = 5123,
    // UnsignedInt
    UnsignedInt = 5125,
};

// Indices of those attributes that deviate from their initialization value.
struct glTFAccessorSparseIndices : glTFProperty {
    // The index of the bufferView with sparse indices. Referenced bufferView
    // can't have ARRAY_BUFFER or ELEMENT_ARRAY_BUFFER target. [required]
    glTFid<glTFBufferView> bufferView = {};
    // The offset relative to the start of the bufferView in bytes. Must be
    // aligned.
    int byteOffset = 0;
    // The indices data type. [required]
    glTFAccessorSparseIndicesComponentType componentType =
        glTFAccessorSparseIndicesComponentType::NotSet;
};

// Array of size `accessor.sparse.count` times number of components storing the
// displaced accessor attributes pointed by `accessor.sparse.indices`.
struct glTFAccessorSparseValues : glTFProperty {
    // The index of the bufferView with sparse values. Referenced bufferView
    // can't have ARRAY_BUFFER or ELEMENT_ARRAY_BUFFER target. [required]
    glTFid<glTFBufferView> bufferView = {};
    // The offset relative to the start of the bufferView in bytes. Must be
    // aligned.
    int byteOffset = 0;
};

// Sparse storage of attributes that deviate from their initialization value.
struct glTFAccessorSparse : glTFProperty {
    // Number of entries stored in the sparse array. [required]
    int count = 0;
    // Index array of size `count` that points to those accessor attributes
    // that deviate from their initialization value. Indices must strictly
    // increase. [required]
    glTFAccessorSparseIndices* indices = nullptr;
    // Array of size `count` times number of components, storing the displaced
    // accessor attributes pointed by `indices`. Substituted values must have
    // the same `componentType` and number of components as the base accessor.
    // [required]
    glTFAccessorSparseValues* values = nullptr;

    ~glTFAccessorSparse() {
        if (indices) delete indices;
        if (values) delete values;
    }
};

// Values for glTFAccessor::componentType
enum class glTFAccessorComponentType {
    // Not set
    NotSet = -1,
    // Byte
    Byte = 5120,
    // UnsignedByte
    UnsignedByte = 5121,
    // Short
    Short = 5122,
    // UnsignedShort
    UnsignedShort = 5123,
    // UnsignedInt
    UnsignedInt = 5125,
    // Float
    Float = 5126,
};

// Values for glTFAccessor::type
enum class glTFAccessorType {
    // Not set
    NotSet = -1,
    // Scalar
    Scalar = 0,
    // Vec2
    Vec2 = 1,
    // Vec3
    Vec3 = 2,
    // Vec4
    Vec4 = 3,
    // Mat2
    Mat2 = 4,
    // Mat3
    Mat3 = 5,
    // Mat4
    Mat4 = 6,
};

// A typed view into a bufferView.  A bufferView contains raw binary data.  An
// accessor provides a typed view into a bufferView or a subset of a bufferView
// similar to how WebGL's `vertexAttribPointer()` defines an attribute in a
// buffer.
struct glTFAccessor : glTFChildOfRootProperty {
    // The index of the bufferView.
    glTFid<glTFBufferView> bufferView = {};
    // The offset relative to the start of the bufferView in bytes.
    int byteOffset = 0;
    // The datatype of components in the attribute. [required]
    glTFAccessorComponentType componentType = glTFAccessorComponentType::NotSet;
    // Specifies whether integer data values should be normalized.
    bool normalized = false;
    // The number of attributes referenced by this accessor. [required]
    int count = 0;
    // Specifies if the attribute is a scalar, vector, or matrix. [required]
    glTFAccessorType type = glTFAccessorType::NotSet;
    // Maximum value of each component in this attribute.
    std::vector<float> max = {};
    // Minimum value of each component in this attribute.
    std::vector<float> min = {};
    // Sparse storage of attributes that deviate from their initialization
    // value.
    glTFAccessorSparse* sparse = nullptr;

    ~glTFAccessor() {
        if (sparse) delete sparse;
    }
};

// Values for glTFAnimationChannelTarget::path
enum class glTFAnimationChannelTargetPath {
    // Not set
    NotSet = -1,
    // Translation
    Translation = 0,
    // Rotation
    Rotation = 1,
    // Scale
    Scale = 2,
    // Weights
    Weights = 3,
};

// The index of the node and TRS property that an animation channel targets.
struct glTFAnimationChannelTarget : glTFProperty {
    // The index of the node to target. [required]
    glTFid<glTFNode> node = {};
    // The name of the node's TRS property to modify, or the "weights" of the
    // Morph Targets it instantiates. [required]
    glTFAnimationChannelTargetPath path =
        glTFAnimationChannelTargetPath::NotSet;
};

// Targets an animation's sampler at a node's property.
struct glTFAnimationChannel : glTFProperty {
    // The index of a sampler in this animation used to compute the value for
    // the target. [required]
    glTFid<glTFAnimationSampler> sampler = {};
    // The index of the node and TRS property to target. [required]
    glTFAnimationChannelTarget* target = nullptr;

    ~glTFAnimationChannel() {
        if (target) delete target;
    }
};

// Values for glTFAnimationSampler::interpolation
enum class glTFAnimationSamplerInterpolation {
    // Not set
    NotSet = -1,
    // The animated values are linearly interpolated between keyframes. When
    // targeting a rotation, spherical linear interpolation (slerp) should be
    // used to interpolate quaternions. The number output of elements must equal
    // the number of input elements.
    Linear = 0,
    // The animated values remain constant to the output of the first keyframe,
    // until the next keyframe. The number of output elements must equal the
    // number of input elements.
    Step = 1,
    // The animation's interpolation is computed using a cubic spline with
    // specified tangents. The number of output elements must equal three times
    // the number of input elements. For each input element, the output stores
    // three elements, an in-tangent, a spline vertex, and an out-tangent. There
    // must be at least two keyframes when using this interpolation.
    CubicSpline = 3,
};

// Combines input and output accessors with an interpolation algorithm to
// define a keyframe graph (but not its target).
struct glTFAnimationSampler : glTFProperty {
    // The index of an accessor containing keyframe input values, e.g., time.
    // [required]
    glTFid<glTFAccessor> input = {};
    // Interpolation algorithm.
    glTFAnimationSamplerInterpolation interpolation =
        glTFAnimationSamplerInterpolation::Linear;
    // The index of an accessor, containing keyframe output values. [required]
    glTFid<glTFAccessor> output = {};
};

// A keyframe animation.
struct glTFAnimation : glTFChildOfRootProperty {
    // An array of channels, each of which targets an animation's sampler at a
    // node's property. Different channels of the same animation can't have
    // equal targets. [required]
    std::vector<glTFAnimationChannel*> channels = {};
    // An array of samplers that combines input and output accessors with an
    // interpolation algorithm to define a keyframe graph (but not its target).
    // [required]
    std::vector<glTFAnimationSampler*> samplers = {};

    // typed access for nodes
    glTFAnimationChannel* get(const glTFid<glTFAnimationChannel>& id) const {
        if (!id) return nullptr;
        return channels.at((int)id);
    }
    // typed access for nodes
    glTFAnimationSampler* get(const glTFid<glTFAnimationSampler>& id) const {
        if (!id) return nullptr;
        return samplers.at((int)id);
    }

    ~glTFAnimation() {
        for (auto v : channels) delete v;
        for (auto v : samplers) delete v;
    }
};

// Metadata about the glTF asset.
struct glTFAsset : glTFProperty {
    // A copyright message suitable for display to credit the content creator.
    std::string copyright = "";
    // Tool that generated this glTF model.  Useful for debugging.
    std::string generator = "";
    // The glTF version that this asset targets. [required]
    std::string version = "";
    // The minimum glTF version that this asset targets.
    std::string minVersion = "";
};

// A buffer points to binary geometry, animation, or skins.
struct glTFBuffer : glTFChildOfRootProperty {
    // The uri of the buffer.
    std::string uri = "";
    // The length of the buffer in bytes. [required]
    int byteLength = 0;
    // Stores buffer content after loading. [required]
    buffer_data data = {};
};

// Values for glTFBufferView::target
enum class glTFBufferViewTarget {
    // Not set
    NotSet = -1,
    // ArrayBuffer
    ArrayBuffer = 34962,
    // ElementArrayBuffer
    ElementArrayBuffer = 34963,
};

// A view into a buffer generally representing a subset of the buffer.
struct glTFBufferView : glTFChildOfRootProperty {
    // The index of the buffer. [required]
    glTFid<glTFBuffer> buffer = {};
    // The offset into the buffer in bytes.
    int byteOffset = 0;
    // The length of the bufferView in bytes. [required]
    int byteLength = 0;
    // The stride, in bytes.
    int byteStride = 0;
    // The target that the GPU buffer should be bound to.
    glTFBufferViewTarget target = glTFBufferViewTarget::NotSet;
};

// An orthographic camera containing properties to create an orthographic
// projection matrix.
struct glTFCameraOrthographic : glTFProperty {
    // The floating-point horizontal magnification of the view. [required]
    float xmag = 0;
    // The floating-point vertical magnification of the view. [required]
    float ymag = 0;
    // The floating-point distance to the far clipping plane. `zfar` must be
    // greater than `znear`. [required]
    float zfar = 0;
    // The floating-point distance to the near clipping plane. [required]
    float znear = 0;
};

// A perspective camera containing properties to create a perspective
// projection matrix.
struct glTFCameraPerspective : glTFProperty {
    // The floating-point aspect ratio of the field of view.
    float aspectRatio = 0;
    // The floating-point vertical field of view in radians. [required]
    float yfov = 0;
    // The floating-point distance to the far clipping plane.
    float zfar = 0;
    // The floating-point distance to the near clipping plane. [required]
    float znear = 0;
};

// Values for glTFCamera::type
enum class glTFCameraType {
    // Not set
    NotSet = -1,
    // Perspective
    Perspective = 0,
    // Orthographic
    Orthographic = 1,
};

// A camera's projection.  A node can reference a camera to apply a transform
// to place the camera in the scene.
struct glTFCamera : glTFChildOfRootProperty {
    // An orthographic camera containing properties to create an orthographic
    // projection matrix.
    glTFCameraOrthographic* orthographic = nullptr;
    // A perspective camera containing properties to create a perspective
    // projection matrix.
    glTFCameraPerspective* perspective = nullptr;
    // Specifies if the camera uses a perspective or orthographic projection.
    // [required]
    glTFCameraType type = glTFCameraType::NotSet;

    ~glTFCamera() {
        if (orthographic) delete orthographic;
        if (perspective) delete perspective;
    }
};

// Values for glTFImage::mimeType
enum class glTFImageMimeType {
    // Not set
    NotSet = -1,
    // ImageJpeg
    ImageJpeg = 0,
    // ImagePng
    ImagePng = 1,
};

// Image data used to create a texture. Image can be referenced by URI or
// `bufferView` index. `mimeType` is required in the latter case.
struct glTFImage : glTFChildOfRootProperty {
    // The uri of the image.
    std::string uri = "";
    // The image's MIME type.
    glTFImageMimeType mimeType = glTFImageMimeType::NotSet;
    // The index of the bufferView that contains the image. Use this instead of
    // the image's uri property.
    glTFid<glTFBufferView> bufferView = {};
    // Stores image content after loading.
    image_data data = {};
};

// Reference to a texture.
struct glTFTextureInfo : glTFProperty {
    // The index of the texture. [required]
    glTFid<glTFTexture> index = {};
    // The set index of texture's TEXCOORD attribute used for texture
    // coordinate mapping.
    int texCoord = 0;
};

// A texture and its sampler.
struct glTFTexture : glTFChildOfRootProperty {
    // The index of the sampler used by this texture. When undefined, a sampler
    // with repeat wrapping and auto filtering should be used.
    glTFid<glTFSampler> sampler = {};
    // The index of the image used by this texture.
    glTFid<glTFImage> source = {};
};

// Normal texture information.
struct glTFMaterialNormalTextureInfo : glTFTextureInfo {
    // The scalar multiplier applied to each normal vector of the normal
    // texture.
    float scale = 1;
};

// Occlusion texture information.
struct glTFMaterialOcclusionTextureInfo : glTFTextureInfo {
    // A scalar multiplier controlling the amount of occlusion applied.
    float strength = 1;
};

// A set of parameter values that are used to define the metallic-roughness
// material model from Physically-Based Rendering (PBR) methodology.
struct glTFMaterialPbrMetallicRoughness : glTFProperty {
    // The material's base color factor.
    vec4f baseColorFactor = {1, 1, 1, 1};
    // The base color texture.
    glTFTextureInfo* baseColorTexture = nullptr;
    // The metalness of the material.
    float metallicFactor = 1;
    // The roughness of the material.
    float roughnessFactor = 1;
    // The metallic-roughness texture.
    glTFTextureInfo* metallicRoughnessTexture = nullptr;

    ~glTFMaterialPbrMetallicRoughness() {
        if (baseColorTexture) delete baseColorTexture;
        if (metallicRoughnessTexture) delete metallicRoughnessTexture;
    }
};

// glTF extension that defines the specular-glossiness material model from
// Physically-Based Rendering (PBR) methodology.
struct glTFMaterialPbrSpecularGlossiness : glTFProperty {
    // The reflected diffuse factor of the material.
    vec4f diffuseFactor = {1, 1, 1, 1};
    // The diffuse texture.
    glTFTextureInfo* diffuseTexture = nullptr;
    // The specular RGB color of the material.
    vec3f specularFactor = {1, 1, 1};
    // The glossiness or smoothness of the material.
    float glossinessFactor = 1;
    // The specular-glossiness texture.
    glTFTextureInfo* specularGlossinessTexture = nullptr;

    ~glTFMaterialPbrSpecularGlossiness() {
        if (diffuseTexture) delete diffuseTexture;
        if (specularGlossinessTexture) delete specularGlossinessTexture;
    }
};

// Values for glTFMaterial::alphaMode
enum class glTFMaterialAlphaMode {
    // Not set
    NotSet = -1,
    // The alpha value is ignored and the rendered output is fully opaque.
    Opaque = 0,
    // The rendered output is either fully opaque or fully transparent depending
    // on the alpha value and the specified alpha cutoff value.
    Mask = 1,
    // The alpha value is used to composite the source and destination areas.
    // The rendered output is combined with the background using the normal
    // painting operation (i.e. the Porter and Duff over operator).
    Blend = 2,
};

// The material appearance of a primitive.
struct glTFMaterial : glTFChildOfRootProperty {
    // A set of parameter values that are used to define the metallic-roughness
    // material model from Physically-Based Rendering (PBR) methodology. When
    // not specified, all the default values of `pbrMetallicRoughness` apply.
    glTFMaterialPbrMetallicRoughness* pbrMetallicRoughness = nullptr;
    // A set of parameter values that are used to define the
    // specular-glossiness material model from Physically-Based Rendering (PBR)
    // methodology. When not specified, all the default values of
    // `pbrMetallicRoughness` apply.
    glTFMaterialPbrSpecularGlossiness* pbrSpecularGlossiness = nullptr;
    // The normal map texture.
    glTFMaterialNormalTextureInfo* normalTexture = nullptr;
    // The occlusion map texture.
    glTFMaterialOcclusionTextureInfo* occlusionTexture = nullptr;
    // The emissive map texture.
    glTFTextureInfo* emissiveTexture = nullptr;
    // The emissive color of the material.
    vec3f emissiveFactor = {0, 0, 0};
    // The alpha rendering mode of the material.
    glTFMaterialAlphaMode alphaMode = glTFMaterialAlphaMode::Opaque;
    // The alpha cutoff value of the material.
    float alphaCutoff = 0.5;
    // Specifies whether the material is double sided.
    bool doubleSided = false;

    ~glTFMaterial() {
        if (pbrMetallicRoughness) delete pbrMetallicRoughness;
        if (pbrSpecularGlossiness) delete pbrSpecularGlossiness;
        if (normalTexture) delete normalTexture;
        if (occlusionTexture) delete occlusionTexture;
        if (emissiveTexture) delete emissiveTexture;
    }
};

// Values for glTFMeshPrimitive::mode
enum class glTFMeshPrimitiveMode {
    // Not set
    NotSet = -1,
    // Points
    Points = 0,
    // Lines
    Lines = 1,
    // LineLoop
    LineLoop = 2,
    // LineStrip
    LineStrip = 3,
    // Triangles
    Triangles = 4,
    // TriangleStrip
    TriangleStrip = 5,
    // TriangleFan
    TriangleFan = 6,
};

// Geometry to be rendered with the given material.
struct glTFMeshPrimitive : glTFProperty {
    // A dictionary object, where each key corresponds to mesh attribute
    // semantic and each value is the index of the accessor containing
    // attribute's data. [required]
    std::map<std::string, glTFid<glTFAccessor>> attributes = {};
    // The index of the accessor that contains the indices.
    glTFid<glTFAccessor> indices = {};
    // The index of the material to apply to this primitive when rendering.
    glTFid<glTFMaterial> material = {};
    // The type of primitives to render.
    glTFMeshPrimitiveMode mode = glTFMeshPrimitiveMode::Triangles;
    // An array of Morph Targets, each  Morph Target is a dictionary mapping
    // attributes (only `POSITION`, `NORMAL`, and `TANGENT` supported) to their
    // deviations in the Morph Target.
    std::vector<std::map<std::string, glTFid<glTFAccessor>>> targets = {};
};

// A set of primitives to be rendered.  A node can contain one mesh.  A node's
// transform places the mesh in the scene.
struct glTFMesh : glTFChildOfRootProperty {
    // An array of primitives, each defining geometry to be rendered with a
    // material. [required]
    std::vector<glTFMeshPrimitive*> primitives = {};
    // Array of weights to be applied to the Morph Targets.
    std::vector<float> weights = {};

    ~glTFMesh() {
        for (auto v : primitives) delete v;
    }
};

// A node in the node hierarchy.  When the node contains `skin`, all
// `mesh.primitives` must contain `JOINTS_0` and `WEIGHTS_0` attributes.  A
// node can have either a `matrix` or any combination of
// `translation`/`rotation`/`scale` (TRS) properties. TRS properties are
// converted to matrices and postmultiplied in the `T * R * S` order to compose
// the transformation matrix; first the scale is applied to the vertices, then
// the rotation, and then the translation. If none are provided, the transform
// is the identity. When a node is targeted for animation (referenced by an
// animation.channel.target), only TRS properties may be present; `matrix` will
// not be present.
struct glTFNode : glTFChildOfRootProperty {
    // The index of the camera referenced by this node.
    glTFid<glTFCamera> camera = {};
    // The indices of this node's children.
    std::vector<glTFid<glTFNode>> children = {};
    // The index of the skin referenced by this node.
    glTFid<glTFSkin> skin = {};
    // A floating-point 4x4 transformation matrix stored in column-major order.
    mat4f matrix = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
    // The index of the mesh in this node.
    glTFid<glTFMesh> mesh = {};
    // The node's unit quaternion rotation in the order (x, y, z, w), where w
    // is the scalar.
    vec4f rotation = {0, 0, 0, 1};
    // The node's non-uniform scale.
    vec3f scale = {1, 1, 1};
    // The node's translation.
    vec3f translation = {0, 0, 0};
    // The weights of the instantiated Morph Target. Number of elements must
    // match number of Morph Targets of used mesh.
    std::vector<float> weights = {};
};

// Values for glTFSampler::magFilter
enum class glTFSamplerMagFilter {
    // Not set
    NotSet = -1,
    // Nearest
    Nearest = 9728,
    // Linear
    Linear = 9729,
};

// Values for glTFSampler::minFilter
enum class glTFSamplerMinFilter {
    // Not set
    NotSet = -1,
    // Nearest
    Nearest = 9728,
    // Linear
    Linear = 9729,
    // NearestMipmapNearest
    NearestMipmapNearest = 9984,
    // LinearMipmapNearest
    LinearMipmapNearest = 9985,
    // NearestMipmapLinear
    NearestMipmapLinear = 9986,
    // LinearMipmapLinear
    LinearMipmapLinear = 9987,
};

// glTFSampler::wrapS
enum class glTFSamplerWrapS {
    // Not set
    NotSet = -1,
    // ClampToEdge
    ClampToEdge = 33071,
    // MirroredRepeat
    MirroredRepeat = 33648,
    // Repeat
    Repeat = 10497,
};

// glTFSampler::wrapT
enum class glTFSamplerWrapT {
    // Not set
    NotSet = -1,
    // ClampToEdge
    ClampToEdge = 33071,
    // MirroredRepeat
    MirroredRepeat = 33648,
    // Repeat
    Repeat = 10497,
};

// Texture sampler properties for filtering and wrapping modes.
struct glTFSampler : glTFChildOfRootProperty {
    // Magnification filter.
    glTFSamplerMagFilter magFilter = glTFSamplerMagFilter::NotSet;
    // Minification filter.
    glTFSamplerMinFilter minFilter = glTFSamplerMinFilter::NotSet;
    // s wrapping mode.
    glTFSamplerWrapS wrapS = glTFSamplerWrapS::Repeat;
    // t wrapping mode.
    glTFSamplerWrapT wrapT = glTFSamplerWrapT::Repeat;
};

// The root nodes of a scene.
struct glTFScene : glTFChildOfRootProperty {
    // The indices of each root node.
    std::vector<glTFid<glTFNode>> nodes = {};
};

// Joints and matrices defining a skin.
struct glTFSkin : glTFChildOfRootProperty {
    // The index of the accessor containing the floating-point 4x4 inverse-bind
    // matrices.  The default is that each matrix is a 4x4 identity matrix,
    // which implies that inverse-bind matrices were pre-applied.
    glTFid<glTFAccessor> inverseBindMatrices = {};
    // The index of the node used as a skeleton root. When undefined, joints
    // transforms resolve to scene root.
    glTFid<glTFNode> skeleton = {};
    // Indices of skeleton nodes, used as joints in this skin. [required]
    std::vector<glTFid<glTFNode>> joints = {};
};

// The root object for a glTF asset.
struct glTF : glTFProperty {
    // Names of glTF extensions used somewhere in this asset.
    std::vector<std::string> extensionsUsed = {};
    // Names of glTF extensions required to properly load this asset.
    std::vector<std::string> extensionsRequired = {};
    // An array of accessors.
    std::vector<glTFAccessor*> accessors = {};
    // An array of keyframe animations.
    std::vector<glTFAnimation*> animations = {};
    // Metadata about the glTF asset. [required]
    glTFAsset* asset = nullptr;
    // An array of buffers.
    std::vector<glTFBuffer*> buffers = {};
    // An array of bufferViews.
    std::vector<glTFBufferView*> bufferViews = {};
    // An array of cameras.
    std::vector<glTFCamera*> cameras = {};
    // An array of images.
    std::vector<glTFImage*> images = {};
    // An array of materials.
    std::vector<glTFMaterial*> materials = {};
    // An array of meshes.
    std::vector<glTFMesh*> meshes = {};
    // An array of nodes.
    std::vector<glTFNode*> nodes = {};
    // An array of samplers.
    std::vector<glTFSampler*> samplers = {};
    // The index of the default scene.
    glTFid<glTFScene> scene = {};
    // An array of scenes.
    std::vector<glTFScene*> scenes = {};
    // An array of skins.
    std::vector<glTFSkin*> skins = {};
    // An array of textures.
    std::vector<glTFTexture*> textures = {};

    // typed access for nodes
    glTFAccessor* get(const glTFid<glTFAccessor>& id) const {
        if (!id) return nullptr;
        return accessors.at((int)id);
    }
    // typed access for nodes
    glTFAnimation* get(const glTFid<glTFAnimation>& id) const {
        if (!id) return nullptr;
        return animations.at((int)id);
    }
    // typed access for nodes
    glTFBuffer* get(const glTFid<glTFBuffer>& id) const {
        if (!id) return nullptr;
        return buffers.at((int)id);
    }
    // typed access for nodes
    glTFBufferView* get(const glTFid<glTFBufferView>& id) const {
        if (!id) return nullptr;
        return bufferViews.at((int)id);
    }
    // typed access for nodes
    glTFCamera* get(const glTFid<glTFCamera>& id) const {
        if (!id) return nullptr;
        return cameras.at((int)id);
    }
    // typed access for nodes
    glTFImage* get(const glTFid<glTFImage>& id) const {
        if (!id) return nullptr;
        return images.at((int)id);
    }
    // typed access for nodes
    glTFMaterial* get(const glTFid<glTFMaterial>& id) const {
        if (!id) return nullptr;
        return materials.at((int)id);
    }
    // typed access for nodes
    glTFMesh* get(const glTFid<glTFMesh>& id) const {
        if (!id) return nullptr;
        return meshes.at((int)id);
    }
    // typed access for nodes
    glTFNode* get(const glTFid<glTFNode>& id) const {
        if (!id) return nullptr;
        return nodes.at((int)id);
    }
    // typed access for nodes
    glTFSampler* get(const glTFid<glTFSampler>& id) const {
        if (!id) return nullptr;
        return samplers.at((int)id);
    }
    // typed access for nodes
    glTFScene* get(const glTFid<glTFScene>& id) const {
        if (!id) return nullptr;
        return scenes.at((int)id);
    }
    // typed access for nodes
    glTFSkin* get(const glTFid<glTFSkin>& id) const {
        if (!id) return nullptr;
        return skins.at((int)id);
    }
    // typed access for nodes
    glTFTexture* get(const glTFid<glTFTexture>& id) const {
        if (!id) return nullptr;
        return textures.at((int)id);
    }

    ~glTF() {
        for (auto v : accessors) delete v;
        for (auto v : animations) delete v;
        if (asset) delete asset;
        for (auto v : buffers) delete v;
        for (auto v : bufferViews) delete v;
        for (auto v : cameras) delete v;
        for (auto v : images) delete v;
        for (auto v : materials) delete v;
        for (auto v : meshes) delete v;
        for (auto v : nodes) delete v;
        for (auto v : samplers) delete v;
        for (auto v : scenes) delete v;
        for (auto v : skins) delete v;
        for (auto v : textures) delete v;
    }
};
// #codegen end gltf-type

// Load a gltf file `filename` from disk. Load binaries and images only if
// `load_bin` and `load_img` are true, reporting errors only if `skip_missing`
// is false.
glTF* load_gltf(const std::string& filename, bool load_bin = true,
    bool load_img = false, bool skip_missing = false);

// Load a binary gltf file `filename` from disk. Load binaries and images only
// if `load_bin` and `load_img` are true, reporting errors only if
// `skip_missing` is false.
glTF* load_binary_gltf(const std::string& filename, bool load_bin = true,
    bool load_img = false, bool skip_missing = false);

// Save a gltf file `filename` to disk. Save binaries and images only if
// `save_bin` and `save_img` are true.
void save_gltf(const std::string& filename, const glTF* gltf,
    bool save_bin = true, bool save_img = false);

// Save a gltf file `filename` to disk. Save binaries and images only if
// `save_bin` and `save_img` are true.
void save_binary_gltf(const std::string& filename, const glTF* gltf,
    bool save_bin = true, bool save_img = false);

// Computes the local node transform and its inverse.
inline mat4f node_transform(const glTFNode* node) {
    return frame_to_mat(translation_frame(node->translation) *
                        rotation_frame(node->rotation) *
                        scaling_frame(node->scale)) *
           node->matrix;
}

// A view for gltf array buffers that allows for typed access.
struct accessor_view {
    // Construct a view from an accessor.
    accessor_view(const glTF* gltf, const glTFAccessor* accessor);

    // Number of elements in the view.
    int size() const { return _size; }
    // Number of elements in the view
    int count() const { return _size; }
    // Number of components per element
    int ncomp() const { return _ncomp; }
    // Check whether the view is valid.
    bool valid() const { return _valid; }

    // Get the idx-th element of fixed length width default values.
    vec2f getv2f(int idx, const vec2f& def = {0, 0}) const {
        auto v = def;
        for (auto i = 0; i < min(_ncomp, 2); i++) (&v.x)[i] = get(idx, i);
        return v;
    }
    // Get the idx-th element of fixed length width default values.
    vec3f getv3f(int idx, const vec3f& def = {0, 0, 0}) const {
        auto v = def;
        for (auto i = 0; i < min(_ncomp, 3); i++) (&v.x)[i] = get(idx, i);
        return v;
    }
    // Get the idx-th element of fixed length width default values.
    vec4f getv4f(int idx, const vec4f& def = {0, 0, 0, 0}) const {
        auto v = def;
        for (auto i = 0; i < min(_ncomp, 4); i++) (&v.x)[i] = get(idx, i);
        return v;
    }

    // Get the idx-th element of fixed length as a matrix.
    mat4f getm4f(int idx) const {
        auto v = mat4f();
        assert(_ncomp == 16);
        auto vm = &v.x.x;
        for (auto i = 0; i < 16; i++) vm[i] = get(idx, i);
        return v;
    }

    // Get the c-th component of the idx-th element.
    float get(int idx, int c = 0) const;

    // Get the idx-th element as integer with fixed length.
    vec2i getv2i(int idx, const vec2i& def = {0, 0}) const {
        auto v = def;
        for (auto i = 0; i < min(_ncomp, 2); i++) { (&v.x)[i] = geti(idx, i); }
        return v;
    }
    // Get the idx-th element as integer with fixed length.
    vec3i getv3i(int idx, const vec3i& def = {0, 0, 0}) const {
        auto v = def;
        for (auto i = 0; i < min(_ncomp, 3); i++) { (&v.x)[i] = geti(idx, i); }
        return v;
    }
    // Get the idx-th element as integer with fixed length.
    vec4i getv4i(int idx, const vec4i& def = {0, 0, 0, 0}) const {
        auto v = def;
        for (auto i = 0; i < min(_ncomp, 4); i++) { (&v.x)[i] = geti(idx, i); }
        return v;
    }

    // Get the c-th component of the idx-th element as integer.
    int geti(int idx, int c = 0) const;

   private:
    const unsigned char* _data = nullptr;
    int _size = 0;
    int _stride = 0;
    int _ncomp = 0;
    glTFAccessorComponentType _ctype;
    bool _normalize = false;
    bool _valid = false;

    static int _num_components(glTFAccessorType type);
    static int _ctype_size(glTFAccessorComponentType componentType);
};

}  // namespace ygl

#endif

// -----------------------------------------------------------------------------
// STRING, PATH AND FILE UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Get directory name (including '/').
inline std::string path_dirname(const std::string& filename) {
    auto pos = filename.rfind('/');
    if (pos == std::string::npos) pos = filename.rfind('\\');
    if (pos == std::string::npos) return "";
    return filename.substr(0, pos + 1);
}

// Get extension (including '.').
inline std::string path_extension(const std::string& filename) {
    auto pos = filename.rfind('.');
    if (pos == std::string::npos) return "";
    return filename.substr(pos);
}

// Get file basename.
inline std::string path_basename(const std::string& filename) {
    auto dirname = path_dirname(filename);
    auto extension = path_extension(filename);
    return filename.substr(
        dirname.size(), filename.size() - dirname.size() - extension.size());
}

// Get filename without directory (equiv to get_basename() +
// get_extension()).
inline std::string path_filename(const std::string& filename) {
    return path_basename(filename) + path_extension(filename);
}

// Replace extension.
inline std::string replace_path_extension(
    const std::string& filename, const std::string& ext) {
    return path_dirname(filename) + path_basename(filename) + ext;
}

// Prepend a string to the extension.
inline std::string prepend_path_extension(
    const std::string& filename, const std::string& prep) {
    return path_dirname(filename) + path_basename(filename) + prep +
           path_extension(filename);
}

// Really-minimal Python like string format. The implementation is not fast
// nor memory efficient. But it is good enough for some needs.
inline std::string format(
    const std::string& fmt, const std::vector<std::string>& args) {
    auto open = false;
    auto cur = 0;
    auto str = std::string();
    for (auto c : fmt) {
        if (c == '{') {
            str += args[cur++];
            open = true;
        } else if (c == '}') {
            if (!open) throw std::runtime_error("bad format");
            open = false;
        } else {
            str += c;
        }
    }
    return str;
}

// format value
inline std::string format_value(const std::string& val) { return val; }
inline std::string format_value(const int& val) { return std::to_string(val); }
inline std::string format_value(const uint64_t& val) {
    return std::to_string(val);
}
inline std::string format_value(const vec2i& val) {
    return std::to_string(val.x) + " " + std::to_string(val.y);
}
inline std::string format_value(const vec3i& val) {
    return std::to_string(val.x) + " " + std::to_string(val.y) + " " +
           std::to_string(val.z);
}
inline std::string format_value(const vec4i& val) {
    return std::to_string(val.x) + " " + std::to_string(val.y) + " " +
           std::to_string(val.z) + " " + std::to_string(val.w);
}
inline std::string format_value(const float& val) {
    return std::to_string(val);
}
inline std::string format_value(const vec2f& val) {
    return std::to_string(val.x) + " " + std::to_string(val.y);
}
inline std::string format_value(const vec3f& val) {
    return std::to_string(val.x) + " " + std::to_string(val.y) + " " +
           std::to_string(val.z);
}
inline std::string format_value(const vec4f& val) {
    return std::to_string(val.x) + " " + std::to_string(val.y) + " " +
           std::to_string(val.z) + " " + std::to_string(val.w);
}
inline std::string format_value(const bool& val) {
    return (val) ? "true" : "false";
}

// Implementation of the function below.
inline void _format_one(std::vector<std::string>& vals) {}
template <typename Arg, typename... Args>
inline void _format_one(
    std::vector<std::string>& vals, const Arg& arg, const Args&... args) {
    vals.push_back(format_value(arg));
    _format_one(vals, args...);
}

// Really-minimal Python like string format. Internally uses streams for
// generality and supports for now only the '{}' operator. The implementation
// is not fast nor memory efficient. But it is good enough for some needs.
template <typename... Args>
inline std::string format(const std::string& fmt, const Args&... args) {
    auto vals = std::vector<std::string>();
    _format_one(vals, args...);
    return format(fmt, vals);
}

// Wrapper for the above function that prints to stdout.
template <typename... Args>
inline void print(const std::string& fmt, const Args&... args) {
    printf("%s", format(fmt, args...).c_str());
}

// Wrapper for the above function that prints to stdout with endline.
template <typename... Args>
inline void println(const std::string& fmt, const Args&... args) {
    printf("%s\n", format(fmt, args...).c_str());
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMMEDIATE MODE COMMAND LINE PARSER
// -----------------------------------------------------------------------------
namespace ygl {

// Immediate mode command line parser. Members are not part of the public API.
struct cmdline_parser {
    std::vector<std::string> _to_parse;    // args left to parse
    std::vector<std::string> _used_names;  // used names for check
    std::string _usage_prog;               // usage prog line
    std::string _usage_help;               // usage help line
    std::string _usage_opts;               // usage option lines
    std::string _usage_args;               // usage argument lines
    bool _usage = false;                   // help option triggered
    std::string _error;                    // parse error
};

// Initialize a command line parser.
inline cmdline_parser make_parser(
    int argc, char** argv, const std::string& prog, const std::string& help);

// Check unused arguments.
inline bool should_exit(cmdline_parser& parser);

// Returns the usage string.
inline std::string get_usage(const cmdline_parser& parser);

// Pase a flag from the command line.
inline bool parse_flag(cmdline_parser& parser, const std::string& name,
    const std::string& flag, const std::string& help, bool def = false,
    bool req = false);

// Pase an option from the command line.
template <typename T>
inline T parse_opt(cmdline_parser& parser, const std::string& name,
    const std::string& flag, const std::string& help, const T& def = {},
    bool req = false, const std::vector<T>& choices = {});

// Parse an enum option from the command line.
template <typename T>
inline T parse_opt(cmdline_parser& parser, const std::string& name,
    const std::string& flag, const std::string& help,
    const std::map<T, std::string>& key_values, const T& def, bool req = false,
    const std::vector<T>& choices = {});

// Parse positional argument from the command line.
template <typename T>
inline T parse_arg(cmdline_parser& parser, const std::string& name,
    const std::string& help, const T& def = {}, bool req = true,
    const std::vector<T>& choices = {});

// Parse all remaining positional argument from the command line.
template <typename T>
inline std::vector<T> parse_args(cmdline_parser& parser,
    const std::string& name, const std::string& help,
    const std::vector<T>& def = {}, bool req = true,
    const std::vector<T>& choices = {});

}  // namespace ygl

// -----------------------------------------------------------------------------
// SIMPLE LOGGER
// -----------------------------------------------------------------------------
namespace ygl {

// Logger object. A logger can output messages to console an a file.
struct logger {
    // whether to output verbose
    bool verbose = true;
    // whether to output to console
    bool console = true;
    // file stream for stream output
    FILE* file = nullptr;

    // cleanup
    ~logger() {
        if (file) fclose(file);
    }
};

// Make a logger with an optional console stream, an optional file stram
// and the specified verbosity level.
inline logger* make_logger(const std::string& filename = "",
    bool console = true, bool verbose = true, bool file_append = true) {
    auto lgr = new logger();
    lgr->verbose = verbose;
    lgr->console = console;
    if (filename.empty()) {
        lgr->file = nullptr;
    } else {
        lgr->file = fopen(filename.c_str(), (file_append) ? "at" : "wt");
        if (!lgr->file)
            throw std::runtime_error("could not open file " + filename);
    }
    return lgr;
}

// Get the default logger.
inline logger* get_default_logger() {
    static auto default_logger = new logger();
    return default_logger;
}

// Log a message. Used internally.
inline void _log_msg(
    const logger* lgr, const std::string& msg, const char* type) {
    char time_buf[1024];
    auto tm = time(nullptr);
    auto ttm = localtime(&tm);  // TODO: use thread safe version

    // short message for console
    if (lgr->console) {
        strftime(time_buf, 1024, "%H:%M:%S", ttm);
        printf("%s %s %s\n", time_buf, type, msg.c_str());
        fflush(stdout);
    }

    // long message for file
    if (lgr->file) {
        strftime(time_buf, 1024, "%Y-%m-%d %H:%M:%S", ttm);
        fprintf(lgr->file, "%s %s %s\n", time_buf, type, msg.c_str());
    }
}

// Log an info message.
template <typename... Args>
inline void log_info(
    const logger* lgr, const std::string& msg, const Args&... args) {
    if (!lgr->verbose) return;
    _log_msg(lgr, format(msg, args...), "INFO ");
}

// Log an info message.
template <typename... Args>
inline void log_warning(
    const logger* lgr, const std::string& msg, const Args&... args) {
    if (!lgr->verbose) return;
    _log_msg(lgr, format(msg, args...), "WARN ");
}

// Log an error message.
template <typename... Args>
inline void log_error(
    const logger* lgr, const std::string& msg, const Args&... args) {
    _log_msg(lgr, format(msg, args...), "ERROR");
}

// Log a fatal message and exit.
template <typename... Args>
inline void log_fatal(
    const logger* lgr, const std::string& msg, const Args&... args) {
    _log_msg(lgr, format(msg, args...), "FATAL");
    exit(1);
}

// Logs a message to the default loggers.
template <typename... Args>
inline void log_info(const std::string& msg, const Args&... args) {
    log_info(get_default_logger(), msg, args...);
}

// Logs a message to the default loggers.
template <typename... Args>
inline void log_warning(const std::string& msg, const Args&... args) {
    log_warning(get_default_logger(), msg, args...);
}

// Logs a message to the default loggers.
template <typename... Args>
inline void log_error(const std::string& msg, const Args&... args) {
    log_error(get_default_logger(), msg, args...);
}

// Logs a message to the default loggers.
template <typename... Args>
inline void log_fatal(const std::string& msg, const Args&... args) {
    log_fatal(get_default_logger(), msg, args...);
}

}  // namespace ygl

#if YGL_OPENGL

// -----------------------------------------------------------------------------
// OPENGL OBJECTS AND FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// OpenGL shape element types.
enum struct glelem_type : int { point = 1, line = 2, triangle = 3 };

// OpenGL light types.
enum struct gllight_type : int {
    point = 0,
    directional = 1,
};

// OpenGL lights
struct gllights {
    std::vector<vec3f> pos;          // light positions
    std::vector<vec3f> ke;           // light intensities
    std::vector<gllight_type> type;  // light types
};

// Checks for GL error and then prints.
bool check_glerror(bool print = true);

// Clear window.
void clear_glbuffers(const vec4f& background = {0, 0, 0, 0});

// Enable/disable depth test, culling, wireframe and blending.
void enable_gldepth_test(bool enabled);
void enable_glculling(bool enabled, bool front = false, bool back = true);
void enable_glwireframe(bool enabled);
void enable_glblending(bool enabled);
void set_glblend_over();

// Set viewport.
void set_glviewport(const vec4i& v);
void set_glviewport(const vec2i& v);

// Reads an image from the the framebuffer.
void read_glimagef(float* pixels, int w, int h, int nc);

// OpenGL texture object. Members are not part of the public API.
struct gltexture {
    uint tid = 0;           // texture id
    int width = 0;          // width
    int height = 0;         // height
    int ncomp = 0;          // number of components
    bool as_float = false;  // stored as float
    bool as_srgb = true;    // stored as sRGB
    bool mipmap = true;     // store with mipmaps
    bool linear = true;     // use linear interpolation
};

// Implementation of update_texture.
void update_gltexture(gltexture& txt, int w, int h, int nc, const void* pixels,
    bool floats, bool linear, bool mipmap, bool as_float, bool as_srgb);

// Updates a texture with pixels values of size w, h with nc number of
// components (1-4). Internally use bytes/floats (as_float), linear/sRGB
// (as_srgb) nearest/linear filtering (linear) and mipmmapping (mipmap).
inline void update_gltexture(gltexture& txt, int w, int h, int nc,
    const float* pixels, bool linear, bool mipmap, bool as_float) {
    update_gltexture(
        txt, w, h, nc, pixels, true, linear, mipmap, as_float, false);
}
inline void update_gltexture(gltexture& txt, int w, int h, int nc,
    const unsigned char* pixels, bool linear, bool mipmap, bool as_srgb) {
    update_gltexture(
        txt, w, h, nc, pixels, false, linear, mipmap, false, as_srgb);
}

// Updates a texture with pixels values from an image.
inline void update_gltexture(gltexture& txt, const image4f& img, bool linear,
    bool mipmap, bool as_float) {
    update_gltexture(txt, img.width, img.height, 4, (float*)img.pixels.data(),
        linear, mipmap, as_float);
}
inline void update_gltexture(gltexture& txt, const image4b& img, bool linear,
    bool mipmap, bool as_srgb) {
    update_gltexture(txt, img.width, img.height, 4,
        (unsigned char*)img.pixels.data(), linear, mipmap, as_srgb);
}

// Binds/unbinds a texture to a texture unit.
void bind_gltexture(const gltexture& txt, uint unit);
void unbind_gltexture(const gltexture& txt, uint unit);

// Clears the texture.
void clear_gltexture(gltexture& txt);

// Get texture id and check if defined.
inline uint get_gltexture_id(const gltexture& txt) { return txt.tid; }
inline bool is_gltexture_valid(const gltexture& txt) { return (bool)txt.tid; }

// Wrap values for OpenGL texture.
enum struct gltexture_wrap { not_set, repeat, clamp, mirror };

// Filter values for OpenGL texture.
enum struct gltexture_filter {
    not_set,
    linear,
    nearest,
    linear_mipmap_linear,
    nearest_mipmap_nearest,
    linear_mipmap_nearest,
    nearest_mipmap_linear
};

// OpenGL texture parameters.
struct gltexture_info {
    gltexture txt = {};                               // texture
    int texcoord = 0;                                 // texture coordinate set
    float scale = 1;                                  // texture scale
    gltexture_wrap wrap_s = gltexture_wrap::not_set;  // wrap mode s
    gltexture_wrap wrap_t = gltexture_wrap::not_set;  // wrap mode s
    gltexture_filter filter_mag = gltexture_filter::not_set;  // mag filter
    gltexture_filter filter_min = gltexture_filter::not_set;  // min filter
};

// OpenGL vertex/element buffer. Members are not part of the public API.
struct glvertex_buffer {
    uint bid = 0;        // buffer id
    int num = 0;         // number of elements
    int ncomp = 0;       // number of components
    bool elems = false;  // element buffer
};

// Updates vertex/element buffers of floats/ints respectively.
void update_glbuffer(glvertex_buffer& buf, bool elems, int num, int ncomp,
    const void* values, bool dynamic);

// Updates the buffer with new data.
inline void update_glbuffer(glvertex_buffer& buf, bool elems,
    const std::vector<float>& values, bool dynamic = false) {
    update_glbuffer(buf, elems, values.size(), 1, values.data(), dynamic);
}
inline void update_glbuffer(glvertex_buffer& buf, bool elems,
    const std::vector<vec2f>& values, bool dynamic = false) {
    update_glbuffer(
        buf, elems, values.size(), 2, (const float*)values.data(), dynamic);
}
inline void update_glbuffer(glvertex_buffer& buf, bool elems,
    const std::vector<vec3f>& values, bool dynamic = false) {
    update_glbuffer(
        buf, elems, values.size(), 3, (const float*)values.data(), dynamic);
}
inline void update_glbuffer(glvertex_buffer& buf, bool elems,
    const std::vector<vec4f>& values, bool dynamic = false) {
    update_glbuffer(
        buf, elems, values.size(), 4, (const float*)values.data(), dynamic);
}
inline void update_glbuffer(glvertex_buffer& buf, bool elems,
    const std::vector<int>& values, bool dynamic = false) {
    update_glbuffer(buf, elems, values.size(), 1, values.data(), dynamic);
}
inline void update_glbuffer(glvertex_buffer& buf, bool elems,
    const std::vector<vec2i>& values, bool dynamic = false) {
    update_glbuffer(
        buf, elems, values.size(), 2, (const int*)values.data(), dynamic);
}
inline void update_glbuffer(glvertex_buffer& buf, bool elems,
    const std::vector<vec3i>& values, bool dynamic = false) {
    update_glbuffer(
        buf, elems, values.size(), 3, (const int*)values.data(), dynamic);
}
inline void update_glbuffer(glvertex_buffer& buf, bool elems,
    const std::vector<vec4i>& values, bool dynamic = false) {
    update_glbuffer(
        buf, elems, values.size(), 4, (const int*)values.data(), dynamic);
}

// Binds/unbinds the buffer at a particular attribute location.
void bind_glbuffer(const glvertex_buffer& buf, uint vattr);
void unbind_glbuffer(const glvertex_buffer& buf, uint vattr);
void unbind_glbuffer(uint vattr);

// Get buffer id and if valid.
inline uint get_glbuffer_id(const glvertex_buffer& buf) { return buf.bid; }
inline bool is_glbuffer_valid(const glvertex_buffer& buf) {
    return (bool)buf.bid;
}
inline bool is_glbuffer_empty(const glvertex_buffer& buf) {
    return !buf.bid || !buf.num;
}

// Clears OpenGL state.
void clear_glbuffer(glvertex_buffer& buf);

// Draws elements.
void draw_glelems(const glvertex_buffer& buf);

// OpenGL program. Members are not part of the public API.
struct glprogram {
    uint pid = 0;  // program id
    uint vid = 0;  // vertex shader is
    uint fid = 0;  // fragment shader id
    uint vao = 0;  // vertex array object id
};

// Creates an OpenGL program from vertex and fragment code.
glprogram make_glprogram(
    const std::string& vertex, const std::string& fragment);

// Get uniform and attribute locations.
int get_gluniform_location(const glprogram& prog, const std::string& name);
int get_glattrib_location(const glprogram& prog, const std::string& name);

// Set uniform values.
void set_gluniform(const glprogram& prog, int var, bool val);
void set_gluniform(const glprogram& prog, int var, int val);
void set_gluniform(const glprogram& prog, int var, float val);
void set_gluniform(const glprogram& prog, int var, const vec2f& val);
void set_gluniform(const glprogram& prog, int var, const vec3f& val);
void set_gluniform(const glprogram& prog, int var, const vec4f& val);
void set_gluniform(const glprogram& prog, int var, const mat4f& val);
void set_gluniform(const glprogram& prog, int var, const frame3f& val);
template <typename T>
inline void set_gluniform(
    const glprogram& prog, const std::string& var, const T& val) {
    set_gluniform(prog, get_gluniform_location(prog, var), val);
}

// Set uniform texture.
void set_gluniform_texture(
    const glprogram& prog, int pos, const gltexture_info& tinfo, uint tunit);
// Set uniform texture with an additionasl texture enable flags.
inline void set_gluniform_texture(const glprogram& prog, int var, int varon,
    const gltexture_info& tinfo, uint tunit) {
    set_gluniform_texture(prog, var, tinfo, tunit);
    set_gluniform(prog, varon, is_gltexture_valid(tinfo.txt));
}

// Set uniform texture.
inline void set_gluniform_texture(const glprogram& prog, const std::string& var,
    const gltexture_info& tinfo, uint tunit) {
    auto loc = get_gluniform_location(prog, var);
    if (loc < 0) throw std::runtime_error("bad OpenGL id");
    return set_gluniform_texture(prog, loc, tinfo, tunit);
}
// Set uniform texture with an additionasl texture enable flags.
inline void set_gluniform_texture(const glprogram& prog, const std::string& var,
    const std::string& varon, const gltexture_info& tinfo, uint tunit) {
    auto loc = get_gluniform_location(prog, var);
    if (loc < 0) throw std::runtime_error("bad OpenGL id");
    auto locon = get_gluniform_location(prog, varon);
    if (locon < 0) throw std::runtime_error("bad OpenGL id");
    return set_gluniform_texture(prog, loc, locon, tinfo, tunit);
}

// Binds a buffer to a vertex attribute, or a constant if the buffer is empty.
void set_glattribute(
    const glprogram& prog, int var, const glvertex_buffer& buf, float def);
void set_glattribute(const glprogram& prog, int var, const glvertex_buffer& buf,
    const vec2f& def);
void set_glattribute(const glprogram& prog, int var, const glvertex_buffer& buf,
    const vec3f& def);
void set_glattribute(const glprogram& prog, int var, const glvertex_buffer& buf,
    const vec4f& def);

// Binds a buffer or constant to a vertex attribute.
template <typename T>
inline void set_glattribute(const glprogram& prog, const std::string& var,
    const glvertex_buffer& buf, const T& def) {
    auto loc = get_glattrib_location(prog, var);
    if (loc < 0) throw std::runtime_error("bad OpenGL id");
    set_glattribute(prog, loc, buf, def);
}

// Check whether the program is valid.
inline bool is_glprogram_valid(const glprogram& prog) { return (bool)prog.pid; }

// Binds/unbinds a program.
void bind_glprogram(const glprogram& prog);
void unbind_glprogram(const glprogram& prog);

// Clears OpenGL state.
void clear_program(glprogram& prog);

}  // namespace ygl

// -----------------------------------------------------------------------------
// OPENGL SCENE SHADER FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Vertex buffers for scene drawing. Members are not part of the public API.
struct glshape {
    glvertex_buffer pos;        // position
    glvertex_buffer norm;       // normals
    glvertex_buffer texcoord;   // texcoord
    glvertex_buffer texcoord1;  // texcoord
    glvertex_buffer tangsp;     // tangent space
    glvertex_buffer color;      // color
    glvertex_buffer points;     // point elements
    glvertex_buffer lines;      // line elements
    glvertex_buffer triangles;  // triangle elements
    glvertex_buffer quads;      // quad elements as tris
    glvertex_buffer beziers;    // bezier elements as l.
    glvertex_buffer edges;      // edge elements
};

// Initialize gl lights.
gllights make_gllights(const scene* scn);

// Update scene textures on the GPU.
void update_gltexture(const texture* txt, gltexture& gtxt);

// Update scene textures on the GPU.
inline std::unordered_map<texture*, gltexture> make_gltextures(
    const scene* scn) {
    auto gtextures = std::unordered_map<texture*, gltexture>();
    for (auto txt : scn->textures) update_gltexture(txt, gtextures[txt]);
    return gtextures;
}

// Clear OpenGL state.
inline void clear_gltextures(std::unordered_map<texture*, gltexture>& txts) {
    for (auto& kv : txts) clear_gltexture(kv.second);
    txts.clear();
}

// Update scene shapes on the GPU.
void update_glshape(const shape* shp, glshape& gshp);

// Clear OpenGL state.
void clear_glshape(glshape& gshp);

// Update scene shapes on the GPU.
inline std::unordered_map<shape*, glshape> make_glshapes(const scene* scn) {
    auto gshapes = std::unordered_map<shape*, glshape>();
    for (auto shp : scn->shapes) update_glshape(shp, gshapes[shp]);
    return gshapes;
}

// Clear OpenGL state.
inline void clear_glshapes(std::unordered_map<shape*, glshape>& shps) {
    for (auto& kv : shps) clear_glshape(kv.second);
    shps.clear();
}

// A shader for displaying images.  Members are not part of the public API.
struct glimage_program {
    glprogram prog;       // program
    glvertex_buffer vbo;  // vertex array
    glvertex_buffer ebo;  // element array
};

// Initialize a stdimage program.
glimage_program make_glimage_program();

// Draws an image texture the stdimage program.
void draw_glimage(const glimage_program& prog, const gltexture& txt,
    const vec2i& win_size, const vec2f& offset, float zoom,
    tonemap_type tonemapper, float exposure);

// Draws an image texture the stdimage program.
inline void draw_glimage(const glimage_program& prog, const gltexture& txt,
    const vec2i& win_size, const vec2f& offset, float zoom) {
    draw_glimage(prog, txt, win_size, offset, zoom, tonemap_type::linear, 0);
}

// Params for stdimage drawing.
struct glimage_params {
    vec2f offset = {0, 0};                          // iomage offset
    float zoom = 1;                                 // image oom
    vec4f background = zero4f;                      // background
    tonemap_type tonemapper = tonemap_type::gamma;  // hdr tonemapping type
    float exposure = 0;                             // hdr exposure
};

// Draws an image texture the stdimage program.
inline void draw_glimage(const glimage_program& prog, const gltexture& txt,
    const vec2i& win_size, const glimage_params& params,
    bool clear_background = true) {
    if (clear_background) clear_glbuffers(params.background);
    draw_glimage(prog, txt, win_size, params.offset, params.zoom,
        params.tonemapper, params.exposure);
}

// Computes the image uv coordinates corresponding to the view parameters.
vec2i get_glimage_coords(const vec2f& mouse_pos, const glimage_params& params);

// Program to shade surfaces with a physically-based standard shader based on
// Phong/GGX. Members are not part of public API.
struct glsurface_program {
    glprogram prog;  // program
    // uniform variable location
    int eyelight_id, exposure_id, tonemap_id, cam_xform_id, cam_xform_inv_id,
        cam_proj_id, lamb_id, lnum_id, lpos_id[16], lke_id[16], ltype_id[16],
        shp_xform_id, shp_normal_offset_id, highlight_id, mtype_id, ke_id,
        kd_id, ks_id, rs_id, op_id, ke_txt_id, ke_txt_on_id, kd_txt_id,
        kd_txt_on_id, ks_txt_id, ks_txt_on_id, rs_txt_id, rs_txt_on_id,
        norm_txt_id, norm_txt_on_id, occ_txt_id, occ_txt_on_id, norm_scale_id,
        occ_scale_id, use_phong_id, double_sided_id, alpha_cutout_id, etype_id,
        efaceted_id;
    // vertex attribute locations
    int pos_id, norm_id, texcoord_id, color_id, tangsp_id;
};

// Initialize a stdsurface shader.
glsurface_program make_glsurface_program();

// Check if the program is valid.
inline bool is_glprogram_valid(const glsurface_program& prog) {
    return is_glprogram_valid(prog.prog);
}

// Starts a frame by setting exposure/gamma values, camera transforms and
// projection. Sets also whether to use full shading or a quick eye light
// preview.
void begin_glsurface_frame(const glsurface_program& prog, bool shade_eyelight,
    tonemap_type tonemap, float exposure, const mat4f& camera_xform,
    const mat4f& camera_xform_inv, const mat4f& camera_proj);

// Ends a frame.
void end_glsurface_frame(const glsurface_program& prog);

// Set shading lights and ambient.
void set_glsurface_lights(
    const glsurface_program& prog, const vec3f& amb, const gllights& lights);

// Begins drawing a shape with transform `xform`.
void begin_glsurface_shape(
    const glsurface_program& prog, const mat4f& xform, float normal_offset = 0);

// End shade drawing.
void end_glsurface_shape(const glsurface_program& prog);

// Sets normal offset.
void set_glsurface_normaloffset(
    const glsurface_program& prog, float normal_offset);

// Set the object as highlighted.
void set_glsurface_highlight(
    const glsurface_program& prog, const vec4f& highlight);

// Set material values with emission `ke`, diffuse `kd`, specular `ks` and
// specular roughness `rs`, opacity `op`. Indicates textures ids with the
// correspoinding `XXX_txt` variables. Sets also normal and occlusion
// maps. Works for points/lines/triangles indicated by `etype`, (diffuse for
// points, Kajiya-Kay for lines, GGX/Phong for triangles). Material `type`
// matches the scene material type.
void set_glsurface_material(const glsurface_program& prog, material_type type,
    const vec3f& ke, const vec3f& kd, const vec3f& ks, float rs, float op,
    const gltexture_info& ke_txt, const gltexture_info& kd_txt,
    const gltexture_info& ks_txt, const gltexture_info& rs_txt,
    const gltexture_info& norm_txt, const gltexture_info& occ_txt,
    bool use_phong, bool double_sided, bool alpha_cutout);

// Set constant material with emission `ke` and opacity `op`.
void set_glsurface_constmaterial(
    const glsurface_program& prog, const vec3f& ke, float op);

// Set element properties.
void set_glsurface_elems(
    const glsurface_program& prog, glelem_type etype, bool faceted);

// Set vertex data with buffers for position pos, normals norm, texture
// coordinates texcoord, per-vertex color color and tangent space tangsp.
void set_glsurface_vert(const glsurface_program& prog,
    const glvertex_buffer& pos, const glvertex_buffer& norm,
    const glvertex_buffer& texcoord, const glvertex_buffer& color,
    const glvertex_buffer& tangsp);

// Params for stdsurface drawing.
struct glsurface_params {
    int resolution = 512;       // image resolution
    bool wireframe = false;     // wireframe drawing
    bool edges = false;         // draw edges
    float edge_offset = 0.01f;  // offset for edges
    bool cutout = false;        // draw with binary transparency
    bool eyelight = false;      // camera light mode
    tonemap_type tonemapper = tonemap_type::gamma;  // tonemapper
    float exposure = 0;                             // exposure
    vec4f background = {0, 0, 0, 0};                // background color
    vec3f ambient = {0, 0, 0};                      // ambient lighting
    vec3f highlight_color = {1, 1, 0};              // highlight color
    vec3f edge_color = {0, 0, 0};                   // edge color
    bool double_sided = false;                      // double sided rendering
    bool cull_backface = false;                     // culling back face
};

// Draw scene with stdsurface program.
void draw_glsurface_scene(const scene* scn, const camera* cam,
    glsurface_program& prog, std::unordered_map<shape*, glshape>& shapes,
    std::unordered_map<texture*, gltexture>& textures, const gllights& lights,
    const vec2i& viewport_size, const void* highlighted,
    const glsurface_params& params);

}  // namespace ygl

// Forward declaration
struct GLFWwindow;

// -----------------------------------------------------------------------------
// OPENGL WINDOWS
// -----------------------------------------------------------------------------
namespace ygl {

// Forward declaration
struct glwindow;

// Callbacks.
using text_glcallback = std::function<void(unsigned int key)>;
using mouse_glcallback = std::function<void(int button, bool press, int mods)>;
using refresh_glcallback = std::function<void()>;

// OpenGL window. Members are not part of the public API.
struct glwindow {
    GLFWwindow* gwin = nullptr;               // GLFW window
    bool widget_enabled = false;              // whether we have widgets
    text_glcallback text_cb = nullptr;        // text callback
    mouse_glcallback mouse_cb = nullptr;      // mouse callback
    refresh_glcallback refresh_cb = nullptr;  // refresh callback

    ~glwindow();  // cleaup
};

// Initialize a window.
glwindow* make_glwindow(
    int width, int height, const std::string& title, bool opengl4 = true);

// Set window callbacks.
void set_glwindow_callbacks(glwindow* win, text_glcallback text_cb,
    mouse_glcallback mouse_cb, refresh_glcallback refresh_cb);

// Set window title.
void set_glwindow_title(glwindow* win, const std::string& title);

// Event processing.
void wait_glwindow_events(glwindow* win);
void poll_glwindow_events(glwindow* win);
void swap_glwindow_buffers(glwindow* win);
#ifdef __APPLE__
void wait_glwindow_events_timeout(glwindow* win, double timeout_sec);
void post_glwindow_event(glwindow* win);
#endif
// Whether the window should exit the event processing loop.
bool should_glwindow_close(glwindow* win);

// Window/framebuffer size.
vec2i get_glwindow_size(glwindow* win);
vec2i get_glframebuffer_size(glwindow* win);

// Mouse/keyboard state queries.
int get_glmouse_button(glwindow* win);
vec2i get_glmouse_pos(glwindow* win);
vec2f get_glmouse_posf(glwindow* win);
bool get_glkey(glwindow* win, int key);
bool get_glalt_key(glwindow* win);
bool get_glctrl_key(glwindow* win);
bool get_glshift_key(glwindow* win);

// Read pixels.
image4b take_glscreenshot4b(
    glwindow* win, bool flipy = true, bool back = false);

// Handle camera navigation and scene selection
bool handle_glcamera_navigation(
    glwindow* win, camera* cam, bool navigation_fps);
bool handle_glscene_selection(glwindow* win, const scene* scn,
    const camera* cam, const bvh_tree* bvh, int res,
    const glimage_params& params, scene_selection& sel);

}  // namespace ygl

// -----------------------------------------------------------------------------
// OPENGL WIDGETS
// -----------------------------------------------------------------------------
namespace ygl {

// Initialize widgets.
void init_imgui(
    glwindow* win, bool light_style = false, bool extra_font = true);

// Begin/end draw widgets.
bool begin_imgui_frame(glwindow* win, const std::string& title);
void end_imgui_frame(glwindow* win);

// Whether widgets are active.
bool get_imgui_active(glwindow* win);

// Horizontal separator.
void draw_imgui_separator(glwindow* win);

// Indent and line continuation widget.
void begin_imgui_indent(glwindow* win);
void end_imgui_indent(glwindow* win);
void continue_imgui_line(glwindow* win);

// Label widget.
void draw_imgui_label(
    glwindow* win, const std::string& lbl, const std::string& msg);
template <typename... Args>
inline void draw_imgui_label(glwindow* win, const std::string& lbl,
    const std::string& fmt, const Args&... args) {
    draw_imgui_label(win, lbl, format(fmt, args...));
}
template <typename T>
inline void draw_imgui_label(
    glwindow* win, const std::string& lbl, const T& val) {
    draw_imgui_label(win, lbl, "{}", val);
}

// Checkbox widget
bool draw_imgui_checkbox(glwindow* win, const std::string& lbl, bool& val);
// Text widget.
bool draw_imgui_text(glwindow* win, const std::string& lbl, std::string& str);
bool draw_imgui_multiline_text(
    glwindow* win, const std::string& lbl, std::string& str);

// Drag widget scale (defaults to 1/100).
void draw_drag_speedscale(float scale);
// Drag widget.
bool draw_imgui_dragbox(
    glwindow* win, const std::string& lbl, int& val, int min = 0, int max = 1);
bool draw_imgui_dragbox(glwindow* win, const std::string& lbl, vec2i& val,
    int min = 0, int max = 1);
bool draw_imgui_dragbox(glwindow* win, const std::string& lbl, vec3i& val,
    int min = 0, int max = 1);
bool draw_imgui_dragbox(glwindow* win, const std::string& lbl, vec4i& val,
    int min = 0, int max = 1);
bool draw_imgui_dragbox(glwindow* win, const std::string& lbl, float& val,
    float min = 0, float max = 1);
bool draw_imgui_dragbox(glwindow* win, const std::string& lbl, vec2f& val,
    float min = 0, float max = 1);
bool draw_imgui_dragbox(glwindow* win, const std::string& lbl, vec3f& val,
    float min = 0, float max = 1);
bool draw_imgui_dragbox(glwindow* win, const std::string& lbl, vec4f& val,
    float min = 0, float max = 1);
bool draw_imgui_dragbox(glwindow* win, const std::string& lbl, mat4f& val,
    float min = -1, float max = 1);
bool draw_imgui_dragbox(glwindow* win, const std::string& lbl, frame3f& val,
    float min = -10, float max = 10);

// Color widget.
bool draw_imgui_colorbox(glwindow* win, const std::string& lbl, vec4b& val);
bool draw_imgui_colorbox(glwindow* win, const std::string& lbl, vec4f& val);
bool draw_imgui_colorbox(glwindow* win, const std::string& lbl, vec3f& val);
bool draw_hdr_color_widget(
    glwindow* win, const std::string& lbl, vec3f& val, float max = 10);

// Combo widget.
bool begin_imgui_combobox(
    glwindow* win, const std::string& lbl, const std::string& label);
bool draw_imgui_item(
    glwindow* win, const std::string& label, int idx, bool selected);
void end_imgui_combobox(glwindow* win);
template <typename T>
bool draw_imgui_item(
    glwindow* win, const std::string& label, int idx, T& val, const T& item) {
    auto selected = draw_imgui_item(win, label, idx, val == item);
    if (selected) val = item;
    return selected;
}
inline bool draw_imgui_combobox(glwindow* win, const std::string& lbl,
    std::string& val, const std::vector<std::string>& labels);
template <typename T>
inline bool draw_imgui_combobox(glwindow* win, const std::string& lbl, T& val,
    const std::map<T, std::string>& labels);
template <typename T>
inline bool draw_imgui_combobox(glwindow* win, const std::string& lbl, T*& val,
    const std::vector<T*>& vals, bool extra = true, T* extra_val = nullptr);

// Button widget.
bool draw_imgui_button(glwindow* win, const std::string& lbl);

// Collapsible header widget.
bool draw_imgui_header(glwindow* win, const std::string& lbl);

// Tree widget.
bool begin_imgui_tree(glwindow* win, const std::string& lbl);
void end_imgui_tree(glwindow* win);
bool begin_imgui_tree(
    glwindow* win, const std::string& lbl, void*& selection, void* content);
template <typename T>
inline bool begin_imgui_tree(
    glwindow* win, const std::string& lbl, T*& selection, T* content) {
    auto sel = selection;
    auto open = begin_imgui_tree(win, lbl, (void*&)sel, (void*)content);
    if (sel == content) selection = content;
    return open;
}
void end_imgui_tree(glwindow* win, void* content);
void draw_imgui_tree_leaf(
    glwindow* win, const std::string& lbl, void*& selection, void* content);
template <typename T>
inline void draw_imgui_tree_leaf(
    glwindow* win, const std::string& lbl, void*& selection, T* content) {
    auto sel = selection;
    draw_imgui_tree_leaf(win, lbl, sel, content);
    if (sel == content) selection = content;
}
void draw_imgui_tree_leaf(glwindow* win, const std::string& lbl,
    void*& selection, void* content, const vec4f& col);

// Image widget.
void draw_imgui_imagebox(
    glwindow* win, int tid, const vec2i& size, const vec2i& imsize);
void draw_imgui_imagebox(glwindow* win, gltexture& txt, const vec2i& size);

// Scroll region widget.
void begin_imgui_scrollarea(
    glwindow* win, const std::string& lbl, int height, bool border);
void end_imgui_scrollarea(glwindow* win);
void move_imgui_scrollarea(glwindow* win);

// Group ids widget.
void push_imgui_groupid(glwindow* win, int gid);
void push_imgui_groupid(glwindow* win, void* gid);
void push_imgui_groupid(glwindow* win, const void* gid);
void push_imgui_groupid(glwindow* win, const char* gid);
void pop_imgui_groupid(glwindow* win);

// Widget style.
void push_imgui_style(glwindow* win, const vec4f& color);
void pop_imgui_style(glwindow* win);

// Image inspection widgets.
void draw_imgui_image_inspector(glwindow* win, const std::string& lbl,
    const image4f& hdr, const image4b& ldr, const vec2f& mouse_pos,
    const glimage_params& params);

// Draws widgets for params.
bool draw_imgui_stdimage_inspector(
    glwindow* win, const std::string& lbl, glimage_params& params);
bool draw_imgui_stdsurface_inspector(
    glwindow* win, const std::string& lbl, glsurface_params& params);
bool draw_imgui_trace_inspector(
    glwindow* win, const std::string& lbl, trace_params& params);

// Draws a widget that can selected the camera.
inline bool draw_imgui_camera_selector(glwindow* win, const std::string& lbl,
    camera*& cam, const scene* scn, camera* view) {
    return draw_imgui_combobox(win, lbl, cam, scn->cameras, true, view);
}

// Draws widgets for a camera. Used for quickly making demos.
bool draw_imgui_camera_inspector(
    glwindow* win, const std::string& lbl, camera* cam);

// Draws widgets for a whole scene. Used for quickly making demos.
bool draw_imgui_scene_tree(glwindow* win, const std::string& lbl, scene* scn,
    scene_selection& sel, std::vector<ygl::scene_selection>& update_list,
    const std::unordered_map<std::string, std::string>& inspector_highlights =
        {});

// Draws widgets for a whole scene. Used for quickly making demos.
bool draw_imgui_scene_inspector(glwindow* win, const std::string& lbl,
    scene* scn, scene_selection& sel,
    std::vector<ygl::scene_selection>& update_list,
    const std::unordered_map<std::string, std::string>& inspector_highlights =
        {});

}  // namespace ygl

#endif

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR MATRICES
// -----------------------------------------------------------------------------
namespace ygl {

// Matrix adjugates, determinant and inverses.
inline mat2f adjugate(const mat2f& a) {
    return {{a.y.y, -a.x.y}, {-a.y.x, a.x.x}};
}
inline mat3f adjugate(const mat3f& a) {
    return {{a.y.y * a.z.z - a.z.y * a.y.z, a.z.y * a.x.z - a.x.y * a.z.z,
                a.x.y * a.y.z - a.y.y * a.x.z},
        {a.y.z * a.z.x - a.z.z * a.y.x, a.z.z * a.x.x - a.x.z * a.z.x,
            a.x.z * a.y.x - a.y.z * a.x.x},
        {a.y.x * a.z.y - a.z.x * a.y.y, a.z.x * a.x.y - a.x.x * a.z.y,
            a.x.x * a.y.y - a.y.x * a.x.y}};
}
inline mat4f adjugate(const mat4f& a) {
    return {{a.y.y * a.z.z * a.w.w + a.w.y * a.y.z * a.z.w +
                    a.z.y * a.w.z * a.y.w - a.y.y * a.w.z * a.z.w -
                    a.z.y * a.y.z * a.w.w - a.w.y * a.z.z * a.y.w,
                a.x.y * a.w.z * a.z.w + a.z.y * a.x.z * a.w.w +
                    a.w.y * a.z.z * a.x.w - a.w.y * a.x.z * a.z.w -
                    a.z.y * a.w.z * a.x.w - a.x.y * a.z.z * a.w.w,
                a.x.y * a.y.z * a.w.w + a.w.y * a.x.z * a.y.w +
                    a.y.y * a.w.z * a.x.w - a.x.y * a.w.z * a.y.w -
                    a.y.y * a.x.z * a.w.w - a.w.y * a.y.z * a.x.w,
                a.x.y * a.z.z * a.y.w + a.y.y * a.x.z * a.z.w +
                    a.z.y * a.y.z * a.x.w - a.x.y * a.y.z * a.z.w -
                    a.z.y * a.x.z * a.y.w - a.y.y * a.z.z * a.x.w},
        {a.y.z * a.w.w * a.z.x + a.z.z * a.y.w * a.w.x + a.w.z * a.z.w * a.y.x -
                a.y.z * a.z.w * a.w.x - a.w.z * a.y.w * a.z.x -
                a.z.z * a.w.w * a.y.x,
            a.x.z * a.z.w * a.w.x + a.w.z * a.x.w * a.z.x +
                a.z.z * a.w.w * a.x.x - a.x.z * a.w.w * a.z.x -
                a.z.z * a.x.w * a.w.x - a.w.z * a.z.w * a.x.x,
            a.x.z * a.w.w * a.y.x + a.y.z * a.x.w * a.w.x +
                a.w.z * a.y.w * a.x.x - a.x.z * a.y.w * a.w.x -
                a.w.z * a.x.w * a.y.x - a.y.z * a.w.w * a.x.x,
            a.x.z * a.y.w * a.z.x + a.z.z * a.x.w * a.y.x +
                a.y.z * a.z.w * a.x.x - a.x.z * a.z.w * a.y.x -
                a.y.z * a.x.w * a.z.x - a.z.z * a.y.w * a.x.x},
        {a.y.w * a.z.x * a.w.y + a.w.w * a.y.x * a.z.y + a.z.w * a.w.x * a.y.y -
                a.y.w * a.w.x * a.z.y - a.z.w * a.y.x * a.w.y -
                a.w.w * a.z.x * a.y.y,
            a.x.w * a.w.x * a.z.y + a.z.w * a.x.x * a.w.y +
                a.w.w * a.z.x * a.x.y - a.x.w * a.z.x * a.w.y -
                a.w.w * a.x.x * a.z.y - a.z.w * a.w.x * a.x.y,
            a.x.w * a.y.x * a.w.y + a.w.w * a.x.x * a.y.y +
                a.y.w * a.w.x * a.x.y - a.x.w * a.w.x * a.y.y -
                a.y.w * a.x.x * a.w.y - a.w.w * a.y.x * a.x.y,
            a.x.w * a.z.x * a.y.y + a.y.w * a.x.x * a.z.y +
                a.z.w * a.y.x * a.x.y - a.x.w * a.y.x * a.z.y -
                a.z.w * a.x.x * a.y.y - a.y.w * a.z.x * a.x.y},
        {a.y.x * a.w.y * a.z.z + a.z.x * a.y.y * a.w.z + a.w.x * a.z.y * a.y.z -
                a.y.x * a.z.y * a.w.z - a.w.x * a.y.y * a.z.z -
                a.z.x * a.w.y * a.y.z,
            a.x.x * a.z.y * a.w.z + a.w.x * a.x.y * a.z.z +
                a.z.x * a.w.y * a.x.z - a.x.x * a.w.y * a.z.z -
                a.z.x * a.x.y * a.w.z - a.w.x * a.z.y * a.x.z,
            a.x.x * a.w.y * a.y.z + a.y.x * a.x.y * a.w.z +
                a.w.x * a.y.y * a.x.z - a.x.x * a.y.y * a.w.z -
                a.w.x * a.x.y * a.y.z - a.y.x * a.w.y * a.x.z,
            a.x.x * a.y.y * a.z.z + a.z.x * a.x.y * a.y.z +
                a.y.x * a.z.y * a.x.z - a.x.x * a.z.y * a.y.z -
                a.y.x * a.x.y * a.z.z - a.z.x * a.y.y * a.x.z}};
}
inline float determinant(const mat2f& a) {
    return a.x.x * a.y.y - a.x.y * a.y.x;
}
inline float determinant(const mat3f& a) {
    return a.x.x * (a.y.y * a.z.z - a.z.y * a.y.z) +
           a.x.y * (a.y.z * a.z.x - a.z.z * a.y.x) +
           a.x.z * (a.y.x * a.z.y - a.z.x * a.y.y);
}
inline float determinant(const mat4f& a) {
    return a.x.x * (a.y.y * a.z.z * a.w.w + a.w.y * a.y.z * a.z.w +
                       a.z.y * a.w.z * a.y.w - a.y.y * a.w.z * a.z.w -
                       a.z.y * a.y.z * a.w.w - a.w.y * a.z.z * a.y.w) +
           a.x.y * (a.y.z * a.w.w * a.z.x + a.z.z * a.y.w * a.w.x +
                       a.w.z * a.z.w * a.y.x - a.y.z * a.z.w * a.w.x -
                       a.w.z * a.y.w * a.z.x - a.z.z * a.w.w * a.y.x) +
           a.x.z * (a.y.w * a.z.x * a.w.y + a.w.w * a.y.x * a.z.y +
                       a.z.w * a.w.x * a.y.y - a.y.w * a.w.x * a.z.y -
                       a.z.w * a.y.x * a.w.y - a.w.w * a.z.x * a.y.y) +
           a.x.w * (a.y.x * a.w.y * a.z.z + a.z.x * a.y.y * a.w.z +
                       a.w.x * a.z.y * a.y.z - a.y.x * a.z.y * a.w.z -
                       a.w.x * a.y.y * a.z.z - a.z.x * a.w.y * a.y.z);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR IMMEDIATE MODE COMMAND LINE PARSER
// -----------------------------------------------------------------------------
namespace ygl {

// cmdline implementation
inline void _check_name(cmdline_parser& parser, const std::string& name,
    const std::string& flag, bool opt) {
    if (opt) {
        if (name.size() < 3 || name[0] != '-' || name[1] != '-' ||
            name[2] == '-')
            throw std::runtime_error("bad name " + name);
    } else {
        if (name.size() < 1 || name[0] == '-')
            throw std::runtime_error("bad name " + name);
    }
    if (find(parser._used_names.begin(), parser._used_names.end(), name) !=
        parser._used_names.end())
        throw std::runtime_error("already used " + name);
    parser._used_names.push_back(name);
    if (flag.empty()) return;
    if (flag.size() < 2 || flag[0] != '-' || flag[1] == '-')
        throw std::runtime_error("bad name " + flag);
    if (find(parser._used_names.begin(), parser._used_names.end(), flag) !=
        parser._used_names.end())
        throw std::runtime_error("already used " + flag);
    parser._used_names.push_back(flag);
}

// cmdline implementation
template <typename T>
inline void _add_usage_str(cmdline_parser& parser, const std::string& name,
    const std::string& flag, bool opt, const std::string& metavar,
    const std::string& help, const std::string& def, bool req,
    const std::vector<T>& choices) {
    auto str = ""s;
    str += "  " + name;
    if (!flag.empty()) str += "/" + flag;
    if (!metavar.empty()) str += " " + metavar;
    while (str.length() < 32) str += " ";
    str += help + " ";
    if (!req && !def.empty()) str += "[" + def + "]";
    if (req) str += "(required)";
    str += "\n";
    if (!choices.empty()) {
        for (auto i = 0; i < 32; i++) str += " ";
        str += "(";
        auto first = true;
        for (auto&& c : choices) {
            if (!first) str += ",";
            str += format_value(c);
            first = false;
        }
        str += ")";
        str += "\n";
    }
    if (opt)
        parser._usage_opts += str;
    else
        parser._usage_args += str;
}

// cmdline implementation
inline void _set_error(cmdline_parser& parser, const std::string& err) {
    if (parser._error.empty()) parser._error = err;
}

// Initialize a command line parser.
inline cmdline_parser make_parser(
    int argc, char** argv, const std::string& prog, const std::string& help) {
    auto parser = cmdline_parser();
    parser._to_parse = std::vector<std::string>(argv + 1, argv + argc);
    parser._usage_prog = (prog.empty()) ? std::string(argv[0]) : prog;
    parser._usage_help = help;
    parser._usage =
        parse_flag(parser, "--help", "-?", "prints and help message");
    return parser;
}

// Check unused arguments.
inline bool should_exit(cmdline_parser& parser) {
    for (auto&& v : parser._to_parse) {
        if (v[0] == '-')
            _set_error(parser, "unknown option " + v);
        else
            _set_error(parser, "unknown argument " + v);
    }
    return parser._usage || !parser._error.empty();
}

// Returns the usage string.
inline std::string get_usage(const cmdline_parser& parser) {
    auto str = std::string();
    if (!parser._error.empty()) str += "error: " + parser._error + "\n\n";
    str += parser._usage_prog;
    if (!parser._usage_opts.empty()) str += " [options]";
    if (!parser._usage_args.empty()) str += " <arguments>";
    str += "\n";
    // while (str.size() < 32) str += " ";
    str += parser._usage_help + "\n\n";
    if (!parser._usage_opts.empty())
        str += "options:\n" + parser._usage_opts + "\n";
    if (!parser._usage_args.empty())
        str += "arguments:\n" + parser._usage_args + "\n";
    return str;
}

// Pase a flag from the command line.
inline bool parse_flag(cmdline_parser& parser, const std::string& name,
    const std::string& flag, const std::string& help, bool def, bool req) {
    // check names
    _check_name(parser, name, flag, true);
    // update usage
    _add_usage_str(
        parser, name, flag, true, "", help, "", req, std::vector<bool>{});
    // skip if error
    if (!parser._error.empty()) return def;
    // find location of option
    auto pos = find(parser._to_parse.begin(), parser._to_parse.end(), name);
    if (pos == parser._to_parse.end())
        pos = find(parser._to_parse.begin(), parser._to_parse.end(), flag);
    if (pos == parser._to_parse.end()) {
        if (req) _set_error(parser, "missing required flag " + name);
        return def;
    }
    // remove parsed arg
    parser._to_parse.erase(pos, pos + 1);
    // done
    return !def;
}

// Parse a value
inline bool parse_value(const std::string& str, int& val) {
    return sscanf(str.c_str(), "%d", &val) == 1;
}
inline bool parse_value(const std::string& str, float& val) {
    return sscanf(str.c_str(), "%f", &val) == 1;
}
inline bool parse_value(const std::string& str, vec2f& val) {
    return sscanf(str.c_str(), "%f%f", &val.x, &val.y) == 2;
}
inline bool parse_value(const std::string& str, vec3f& val) {
    return sscanf(str.c_str(), "%f%f%f", &val.x, &val.y, &val.z) == 3;
}
inline bool parse_value(const std::string& str, vec4f& val) {
    return sscanf(str.c_str(), "%f%f%f%f", &val.x, &val.y, &val.z, &val.w) == 4;
}
inline bool parse_value(const std::string& str, std::string& val) {
    val = str;
    return true;
}
inline bool parse_value(const std::string& str, bool& val) {
    if (str == "true") {
        val = true;
        return true;
    }
    if (str == "false") {
        val = false;
        return true;
    }
    return false;
}

// Pase an option from the command line.
template <typename T>
inline T parse_opt(cmdline_parser& parser, const std::string& name,
    const std::string& flag, const std::string& help, const T& def, bool req,
    const std::vector<T>& choices) {
    // check names
    _check_name(parser, name, flag, true);
    // update usage
    _add_usage_str(parser, name, flag, true, "<val>", help, format_value(def),
        req, choices);
    // skip if error
    if (!parser._error.empty()) return def;
    // find location of option
    auto pos = find(parser._to_parse.begin(), parser._to_parse.end(), name);
    if (pos == parser._to_parse.end())
        pos = find(parser._to_parse.begin(), parser._to_parse.end(), flag);
    if (pos == parser._to_parse.end()) {
        if (req) _set_error(parser, "missing option " + name);
        return def;
    }
    // check if value exists
    if (pos == parser._to_parse.end() - 1) {
        _set_error(parser, "no value for parameter " + name);
        return def;
    }
    // get value
    auto val = def;
    const auto& arg = *(pos + 1);
    // parse
    if (!parse_value(arg, val)) {
        _set_error(
            parser, "incorrect value \"" + arg + "\" for option " + name);
        return def;
    }
    // validate if necessary
    if (!choices.empty()) {
        if (find(choices.begin(), choices.end(), val) == choices.end())
            _set_error(
                parser, "incorrect value \"" + arg + "\" for option " + name);
    }
    // remove parsed arg
    parser._to_parse.erase(pos, pos + 2);
    // done
    return val;
}

// Parse an enum option from the command line.
template <typename T>
inline T parse_opt(cmdline_parser& parser, const std::string& name,
    const std::string& flag, const std::string& help,
    const std::map<T, std::string>& key_values, const T& def, bool req,
    const std::vector<T>& choices) {
    auto keys = std::vector<std::string>{};
    auto key_def = key_values.at(def);
    for (auto&& kv : key_values) keys.push_back(kv.second);
    auto key =
        parse_opt<std::string>(parser, name, flag, help, key_def, req, keys);
    if (!parser._error.empty()) return def;
    auto val = def;
    for (auto&& kv : key_values) {
        if (kv.second == key) val = kv.first;
    }
    return val;
}

// Parse positional argument from the command line.
template <typename T>
inline T parse_arg(cmdline_parser& parser, const std::string& name,
    const std::string& help, const T& def, bool req,
    const std::vector<T>& choices) {
    // check names
    _check_name(parser, name, "", false);
    // update usage
    _add_usage_str(
        parser, name, "", false, "", help, format_value(def), req, choices);
    // skip if error
    if (!parser._error.empty()) return def;
    // find location of argument
    auto pos = std::find_if(parser._to_parse.begin(), parser._to_parse.end(),
        [](const auto& s) { return s.size() > 0 && s[0] != '-'; });
    if (pos == parser._to_parse.end()) {
        if (req) _set_error(parser, "missing argument " + name);
        return def;
    }
    // get value
    auto val = def;
    const auto& arg = *(pos);
    // parse
    if (!parse_value(arg, val)) {
        _set_error(
            parser, "incorrect value \"" + arg + "\" for argument " + name);
    }
    // validate if necessary
    if (!choices.empty()) {
        if (find(choices.begin(), choices.end(), val) == choices.end())
            _set_error(
                parser, "incorrect value \"" + arg + "\" for argument " + name);
    }
    // remove parsed arg
    parser._to_parse.erase(pos, pos + 1);
    // done
    return val;
}

// Parse all remaining positional argument from the command line.
template <typename T>
inline std::vector<T> parse_args(cmdline_parser& parser,
    const std::string& name, const std::string& help, const std::vector<T>& def,
    bool req, const std::vector<T>& choices) {
    // check names
    _check_name(parser, name, "", false);
    // update usage
    _add_usage_str(parser, name, "", false, "", help, "", req, choices);
    // skip if error
    if (!parser._error.empty()) return def;
    // search for all params
    auto vals = std::vector<T>();
    while (true) {
        // find location of argument
        auto pos =
            std::find_if(parser._to_parse.begin(), parser._to_parse.end(),
                [](const auto& s) { return s.size() > 0 && s[0] != '-'; });
        if (pos == parser._to_parse.end()) break;
        // get value
        auto val = T{};
        const auto& arg = *(pos);
        // parse
        if (!parse_value(arg, val)) {
            _set_error(
                parser, "incorrect value \"" + arg + "\" for argument " + name);
        }
        // validate if necessary
        if (!choices.empty()) {
            if (find(choices.begin(), choices.end(), val) == choices.end())
                _set_error(parser,
                    "incorrect value \"" + arg + "\" for argument " + name);
        }
        // remove parsed arg
        parser._to_parse.erase(pos, pos + 1);
        // append value
        vals.push_back(val);
    }
    // check missing
    if (vals.empty()) {
        if (req) _set_error(parser, "missing argument " + name);
        return def;
    }
    // done
    return vals;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR OPENGL WIDGETS
// -----------------------------------------------------------------------------
namespace ygl {

// Combo widget.
inline bool draw_imgui_combobox(glwindow* win, const std::string& lbl,
    std::string& val, const std::vector<std::string>& labels) {
    if (!begin_imgui_combobox(win, lbl, val)) return false;
    auto old_val = val;
    for (auto i = 0; i < labels.size(); i++) {
        draw_imgui_item(win, labels[i], i, val, labels[i]);
    }
    end_imgui_combobox(win);
    return val != old_val;
}

// Combo widget.
template <typename T>
inline bool draw_imgui_combobox(glwindow* win, const std::string& lbl, T& val,
    const std::map<T, std::string>& labels) {
    if (!begin_imgui_combobox(win, lbl, labels.at(val))) return false;
    auto old_val = val;
    auto lid = 0;
    for (auto& kv : labels)
        draw_imgui_item(win, kv.second, lid++, val, kv.first);
    end_imgui_combobox(win);
    return val != old_val;
}

// Combo widget
template <typename T>
inline bool draw_imgui_combobox(glwindow* win, const std::string& lbl, T*& val,
    const std::vector<T*>& vals, bool extra, T* extra_val) {
    if (!begin_imgui_combobox(win, lbl, (val) ? val->name : "<none>"))
        return false;
    auto old_val = val;
    if (extra)
        draw_imgui_item(
            win, (extra_val) ? extra_val->name : "<none>", -1, val, extra_val);
    for (auto i = 0; i < vals.size(); i++) {
        draw_imgui_item(win, vals[i]->name, i, val, vals[i]);
    }
    end_imgui_combobox(win);
    return val != old_val;
}

}  // namespace ygl

#endif
