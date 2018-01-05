///
/// # Yocto/GL: Tiny C++ Library for Physically-based Graphics
///
/// Yocto/GL is a collection utiliies for building physically-based graphics
/// algorithms implemented as a two-file library (`yocto_gl.h`, `yocto_gl.cpp`),
/// and released under the MIT license. Features include:
///
/// - convenience math functions for graphics
/// - static length vectors for 2, 3, 4 length and int and float type
/// - static length matrices for 2x2, 3x3, 4x4 and float type
/// - static length rigid transforms (frames), specialized for 2d and 3d space
/// - linear algebra operations and transforms
/// - axis aligned bounding boxes
/// - rays and ray-primitive intersection
/// - point-primitive distance and overlap tests
/// - normal and tangent computation for meshes and lines
/// - generation of tesselated meshes
/// - mesh refinement with linear tesselation and Catmull-Cark subdivision
/// - random number generation via PCG32
/// - simple image data structure and a few image operations
/// - simple scene format
/// - generation of image examples
/// - generation of scene examples
/// - procedural sun and sky HDR
/// - procedural Perlin noise
/// - BVH for intersection and closest point query
/// - Python-like iterators, string, path and container operations
/// - utilities to load and save entire text and binary files
/// - immediate mode command line parser
/// - simple logger and thread pool
/// - path tracer supporting surfaces and hairs, GGX and MIS
/// - support for loading and saving Wavefront OBJ and Khronos glTF
/// - OpenGL utilities to manage textures, buffers and prograrms
/// - OpenGL shader for image viewing and GGX microfacet and hair rendering
///
/// The current version is 0.1.0. You can access the previous multi-file version
/// with tag "v0.0.1" in this repository.
///
/// ## Credits
///
/// This library includes code from the PCG random number generator,
/// the LLVM thread pool, boost hash_combine, Pixar multijittered sampling,
/// code from "Real-Time Collision Detection" by Christer Ericson, base64
/// encode/decode by René Nyffenegger and public domain code from
/// github.com/sgorsten/linalg, gist.github.com/badboy/6267743 and
/// github.com/nothings/stb_perlin.h.
///
/// This library imports many symbols from std for three reasons: avoid
/// verbosity , esnuring better conventions when calling math functions and
/// allowing easy overriding of std containers if desired. Just do not
/// flatten this namespace into yours if this is a concern.
///
/// For most components of the library, the use should be relatively easy to
/// understand if you are familiar with 3d computer graphics. For more complex
/// components, we follow the usage below.
///
///
/// ## Design Considerations
///
/// Yocto/GL tries to follow a simple programming model inspired by C but with
/// heavy use of operator overloading for math readability. We attempt tp make
/// the code weasy to use rather than as performant as possible. The APIs
/// attempt to make using the code as little error prone as possible, sometimes
/// at the price of some slowdown. We adopt a functional style and only rarely
/// use classes and methods. Using a function style makes the code easier to
/// extend, more explicit in the requirements, and easier to write
/// parallel-friendly APIs. I guess you could call this "data-driven
/// programming". We use templates very little now, after a major refactoring,
/// to improve error reporting, reduce compilation times and make the codebase
/// more accessible to beginners. This lead to a small increase in copied code
/// that we deem ok at this time. Finally, we often import symbols from the
/// standard library rather than using the `std::name` pattern. We found that
/// this improves consistency, especially when using math functions, is
/// significantly more readable when using templates and allows to to more
/// easily switch STL implementation if desired.
///
///
/// ## Compilation
///
/// Yocto/GL is written in C++14, with compilation supported on C++11, and
/// compiles on OSX (clang from Xcode 9+), Linux (gcc 6+, clang 4+)
/// and Windows (MSVC 2017).
///
/// For image loading and saving, Yocto/GL depends on `stb_image.h`,
/// `stb_image_write.h`, `stb_image_resize.h` and `tinyexr.h`. These features
/// can be disabled by defining YGL_IMAGEIO to 0 before including this file.
/// If these features are useful, then the implementation files need to
/// included in the manner described by the respective libraries. To simplify
/// builds, we provice a file that builds these libraries, `stb_image.cpp`.
///
/// To support Khronos glTF, Yocto/GL depends on `json.hpp`. These feature can
/// be disabled by defining YGL_GLTF to 0 before including this file.
///
/// OpenGL utilities include the OpenGL libaries, use GLEW on Windows/Linux,
/// GLFW for windows handling and Dear ImGui for UI support.
/// Since OpenGL is quite onerous and hard to link, its support is disabled by
/// default. You can enable it by defining YGL_OPENGL to 1 before including
/// this file. If you use any of the OpenGL calls, make sure to properly link to
/// the OpenGL libraries on your system. For ImGUI, build with the libraries
/// `imgui.cpp`, `imgui_draw.cpp`, `imgui_impl_glfw_gl3.cpp`.
///
///
/// ## Example Applications
///
/// You can see Yocto/GL in action in the following applications written to
/// test the library:
///
/// - `yview.cpp`: simple OpenGL viewer for OBJ and glTF scenes
/// - `ytrace.cpp`: offline path-tracer
/// - `yitrace.cpp.cpp`: interactive path-tracer
/// - `yscnproc.cpp`: scene manipulation and conversion to/from OBJ and glTF
/// - `ytestgen.cpp`: creates test cases for the path tracer and GL viewer
/// - `yimview.cpp`: HDR/PNG/JPG image viewer with exposure/gamma tone mapping
/// - `yimproc.cpp`: offline image manipulation.
///
/// You can build the example applications using CMake with
///     `mkdir build; cd build; cmake ..; cmake --build`
///
/// Here are two images rendered with the buildin path tracer, where the
/// scenes are crated with the test generator.
///
/// ![Yocto/GL](images/shapes.png)
///
/// ![Yocto/GL](images/lines.png)
///
///
/// ## Usage
///
/// To use the library simply include this file and setup the compilation
/// option as described above.
/// All library features are documented at the definition and should be
/// relatively easy to use if you are familiar with writing graphics code.
/// You can find the extracted documentation at `yocto_gl.md`.
/// Here we give an overview of some of the main features.
///
///
/// ### Small Vectors and Matrices, Frames, Bounding Boxes and Transforms
///
/// We provide common operations for small vectors and matrices typically used
/// in graphics. In particular, we support 2-4 dimensional float vectors
/// `vec2f`, `vec3f`, `vec4f`, 2-4 dimensional int vectors `vec2i`, `vec3i`,
/// `vec4i` and a 4 dimensional byte vector `vec4b`. The float vectors
/// support most arithmetic and vector operations.
///
/// We support 2-4 dimensional float matrices `mat2f`, `mat3f`, `mat4f`, with
/// matrix-matrix and matrix-vector products, trasposes and inverses. Matrices
/// are stored in column-major ordered and are accessed and constructed by
/// column.
///
/// To represent transformations, most of the library facilities prefer the use
/// cooordinate frames, aka rigid transforms, represented as `frame3f`.
/// The structure store three coodinate axis and the frame origin. This is
/// equivalenent to a rigid transform written as a column-major affine
/// matrix. Transform operations are better behaved with this representation.
///
/// We represent coordinate bounds with axis-aligned bounding boxes in 1-4
/// dimensions: `bbox1f`, `bbox2f`, `bbox3f`, `bbox4f`. These types support
/// expansion operation, union and containment. We provide operations to
/// compute bounds for points, lines, triangles and quads.
///
/// For all basic types we support iteration with `begin()`/`end()` pairs
/// and stream inout and output.
///
/// For both matrices and frames we support transform operations for points,
/// vectors and directions (`trasform_point()`, `trasform_vector()`,
/// `trasform_direction()`). For frames we also the support inverse operations
/// (`transform_xxx_inverse()`). Transform matrices and frames can be
/// constructed from basic translation, rotation and scaling, e.g. with
/// `translation_mat4f()` or `translation_frame3f()` repsectively, etc. For
/// rotation we support axis-angle and quaternions, with slerp.
///
///
/// ### Random Number Generation, Noise, Hashing and Monte Carlo support
///
/// This library supportds many facitlities helpful in writing sampling
/// functions targeting path tracing and shape generations.
///
/// 1. Random number generation with PCG32:
///     1. initialize the random number generator with `init_rng()`
///     2. advance the random number state with `advance_rng()`
///     3. if necessary, you can reseed the rng with `seed_rng()`
///     4. generate random integers in an interval with `next_rand1i()`
///     5. generate random floats and double in the [0,1) range with
///        `next_rand1f()`, `next_rand2f()`, `next_rand3f()`, `next_rand1d()`
///     6. you can skip random numbers with `advance_rng()` and get the skipped
///        length with `rng_distance()`
///     7. generate random shaffled sequences with `rng_shuffle()`
/// 2. Perlin noise: `perlin_noise()` to generate Perlin noise with optional
///    wrapping, with fractal variations `perlin_ridge_noise()`,
///    `perlin_fbm_noise()`, `perlin_turbulence_noise()`
/// 3. Integer hashing: public domain hash functions for integer values as
///    `hash_permute()`, `hash_uint32()`, `hash_uint64()`, `hash_uint64_32()`
///    and `hash_combine()`.
/// 4. Monte Carlo support: warp functions from [0,1)^k domains to domains
///    commonly used in path tracing. In particular, use `sample_hemisphere()`,
///    `sample_sphere()`, `sample_hemisphere_cosine()`,
///    `sample_hemisphere_cospower()`. `sample_disk()`. `sample_cylinder()`.
///    `sample_triangle()`. For each warp, you can compute the PDF with
///    `sample_xxx_pdf()`.
///
///
/// ### Python-like container operations and iterators
///
/// To make the code more readable, we adopt Python-like iterations and
/// container operations extensively throughout Yocto/GL. These operations
/// are mostly for internal use but could also be used externally.
///
/// 1. Python iterators with `range()` and `enumerate()`
/// 2. Python operators for containers: support for + and += for `std::vector`
/// 3. Check for containment with `contains`  similarly to `in` in Python
///
///
/// ### Shape Utilities
///
/// The library contains a few function to help with typically geometry
/// manipulation useful to support scene viewing and path tracing.
///
/// 1. compute line tangents, and triangle and quad areas and normals
/// 2. compute barycentric interpolation with `eval_barycentric_line()`,
///    `eval_barycentric_triangle()` and `eval_barycentric_quad()`
/// 3. evaluate Bezier curve and derivatives with `eval_bezier_cubic()` and
///    `eval_bezier_cubic_derivative()`
/// 4. compute smooth normals and tangents with `compute_normals()`
/// 5. compute tangent frames from texture coordinates with
///    `compute_tangent_space()`
/// 6. compute skinning with `compute_skinning()` and
///    `compute_matrix_skinning()`
/// 6. shape creation with `make_points()`, `make_lines()`, `make_uvgrid()`
/// 7. element merging with `marge_elems()`
/// 8. facet elements with `facet_elems()`
/// 9. shape sampling with `sample_points()`, `sample_lines()`,
///    `sample_triangles()`; initialize the sampling CDFs with
///    `sample_points_cdf()`, `sample_lines_cdf()`, `sample_triangles_cdf()`
/// 10. samnple a could of point over a surface with `sample_triangles_points()`
/// 11. get edges and boundaries with `get_edges()` and `get_boundary_edges()`
/// 12. convert quads to triangles with `convert_quads_to_triangles()`
/// 13. convert face varying to vertex shared representations with
///     `convert_face_varying()`
/// 14. subdivide elements by edge splits with `subdivide_elems()` and
///     `subdivide_vert()`
/// 15. Catmull-Clark subdivision surface with `subdivide_catmullclark()` with
///     support for edge and vertex creasing
/// 16. example shapes: `make_cube()`, `make_uvsphere()`, `make_uvhemisphere()`,
///     `make_uvquad()`, `make_uvcube()`, `make_fvcube()`, `make_hair()`,
///     `make_suzanne()`
///
///
/// ### Image and color
///
/// We support simple containers for either 4-byte per pixel sRGB images
/// `image4b`, or 4-float per pixel HDR images `image4f`.
///
/// 1. convert between byte and float images with `srgb_to_linear()` and
///    `linear_to_srgb()`
/// 2. color conversion with `hsv_to_rgb()`, `xyz_to_rgb()` and `rgb_to_xyz()`
/// 3. exposure-gamma tonemapping, with optional filmic curve, with
///    `tonemap_image()`
/// 4. compositing support with `image_over()`
/// 5. example image generation with `m,ake_grid_image()`,
///    `make_checker_image()`, `make_bumpdimple_image()`, `make_ramp_image()`,
///    `make_gammaramp_image()`, `make_gammaramp_imagef()`, `make_uv_image()`,
///    `make_uvgrid_image()`, `make_recuvgrid_image()`
/// 6. bump to normal mapping with `bump_to_normal_map()`
/// 7. HDR sun-sky with `m ake_sunsky_image()`
/// 8. various noise images with `make_noise_image()`, `make_fbm_image()`,
///    `make_ridge_image()`, `make_turbulence_image()`
/// 9. image loading and saving with `load_image4b()`, `load_image4f()`,
///    `save_image4b()`, `save_image4f()`
/// 10. image resizing with `resize_image()`
///
///
/// ### Ray Intersection and Point Overlap Queries
///
/// We support ray-scene intersection for points, lines and triangles
/// accelerated by a simple BVH data structure.  Our BVH is written for minimal
/// code and not maximum speed, but still gives reasonable results. We suggest
/// the use of Intel's Embree as a fast alternative.
///
/// 1. use `ray3f` to represent rays
/// 2. build the BVH with `build_points_bvh()`, `build_points_bvh()` or
///   `build_points_bvh()`
/// 3. perform ray-element intersection with `intersect_points_bvh()`,
///   `intersect_lines_bvh()` and `intersect_triangles_bvh()`
/// 4. perform point overlap queries with `overlap_points_bvh()`,
///   `overlap_lines_bvh()` and `overlap_triangles_bvh()`
/// 5. to support custom elements, use `buid_bvh()`, `intersect_bvh()` and
///   `overlap_bvh()` and provide them with proper callbacks
/// 6. we also experimentally support quads with the `xxx_quads_xxx()` functions
///
///
/// ### Simple scene
///
/// We support a simple scene model used to quickly write demos that lets you
/// load/save Wavefront OBJ and Khronos glTF and perform several simple scene
/// manipulation including ray-scene intersection and closest point queries.
///
/// The geometry model is comprised of a set of shapes, which are indexed
/// collections of points, lines, triangles and quads. Each shape may contain
/// only one element type. Shapes are organized into a scene by creating shape
/// instances, each its own transform. Materials are specified like in glTF and
/// include emission, base-metallic and diffuse-specular parametrization,
/// normal, occlusion and displacement mapping. Finally, the scene containes
/// caemras and environement maps. Quad support in shapes is experimental and
/// mostly supported for loading and saving.
///
/// For low-level access to OBJ/glTF formats, you are best accssing the formats
/// directly with Yocto/Obj and Yocto/glTF. This components provides a
/// simplified high-level access to each format which is sufficient for most
/// applications and tuned for quick creating viewers, renderers and simulators.
///
/// 1. load a scene with `load_scene()` and save it with `save_scene()`.
/// 2. add missing data with `add_elements()`
/// 3. use `compute_bounds()` to compute element bounds
/// 4. can merge scene together with `merge_into()`
/// 5. make example scenes with `make_test_scene()`
///
/// Ray-intersection and closet-point routines supporting points,
/// lines and triangles accelerated by a two-level bounding volume
/// hierarchy (BVH). Quad support is experimental.
///
/// 1. build the bvh with `build_bvh()`
/// 2. perform ray-interseciton tests with `intersect_ray()`
///     - use early_exit=false if you want to know the closest hit point
///     - use early_exit=false if you only need to know whether there is a hit
///     - for points and lines, a radius is required
///     - for triangles, the radius is ignored
/// 2. perform point overlap tests with `overlap_point()` to check whether
///    a point overlaps with an element within a maximum distance
///     - use early_exit as above
///     - for all primitives, a radius is used if defined, but should
///       be very small compared to the size of the primitive since the radius
///       overlap is approximate
/// 3. perform instance overlap queries with `overlap_instance_bounds()`
/// 4. use `refit_bvh()` to recompute the bvh bounds if transforms or vertices
///    are changed (you should rebuild the bvh for large changes)
///
/// Notes: Quads are internally handled as a pair of two triangles v0,v1,v3 and
/// v2,v3,v1, with the u/v coordinates of the second triangle corrected as 1-u
/// and 1-v to produce a quad parametrization where u and v go from 0 to 1. This
/// is equivalent to Intel's Embree.
///
///
/// ### Pathtracing
///
/// We supply a path tracer implementation with support for textured mesh
/// lights, GGX/Phong materials, environment mapping. The interface supports
/// progressive parallel execution. The path tracer takes as input a scene
/// and update pixels in image with traced samples. We use a straightfoward
/// path tracer with MIS and also a few simpler shaders for debugging or
/// quick image generation.
///
/// Materials are represented as sums of an emission term, a diffuse term and
/// a specular microfacet term (GGX or Phong). Only opaque for now. We pick
/// a proper material type for each shape element type (points, lines,
/// triangles).
///
/// Lights are defined as any shape with a material emission term. Additionally
/// one can also add environment maps. But even if you can, you might want to
/// add a large triangle mesh with inward normals instead. The latter is more
/// general (you can even more an arbitrary shape sun). For now only the first
/// env is used.
///
/// 1. build the ray-tracing acceleration structure with `build_bvh()`
/// 2. prepare lights for rendering `update_lights()`
/// 3. define rendering params with the `trace_params` structure
/// 4. render blocks of samples with `trace_block()`
///
/// The code can also run in fully asynchronous mode to preview images in a
/// window.
///
/// 1. build the ray-tracing acceleration structure with `build_bvh()`
/// 2. prepare lights for rendering `update_lights()`
/// 3. define rendering params with the `trace_params` structure
/// 4. initialize the prograssive rendering buffers
/// 5. start the progressive renderer with `trace_async_start()`
/// 7. stop the progressive renderer with `trace_async_stop()`
///
///
/// ### Wavefront OBJ
///
/// Wavefront OBJ/MTL loader and writer with support for points,
/// lines, triangles and general polygons and all materials properties.
/// Contains also a few extensions to easily create demos such as per-vertex
/// color and radius, cameras, environment maps and instances.
/// Can use either a low-level OBJ representation, from this files,
/// or a high level flattened representation included in Yocto/Scn.
///
/// Both in reading and writing, OBJ has no clear convention on the orientation
/// of textures Y axis. So in many cases textures appears flipped. To handle
/// that, use the option to flip textures coordinates on either saving or
/// loading. By default texture coordinates are flipped since this seems
/// the convention found on test cases collected on the web. The value Tr
/// has similar problems, since its relation to opacity is software specific.
/// Again we let the user chose the convension and set the default to the
/// one found on the web.
///
/// In the high level interface, shapes are indexed meshes and are described
/// by arrays of vertex indices for points/lines/triangles and arrays for vertex
/// positions, normals, texcoords, color and radius. The latter two as
/// extensions. Since OBJ is a complex formats that does not match well with
/// current GPU rendering / path tracing algorithms, we adopt a simplification
/// similar to other single file libraries:
/// 1. vertex indices are unique, as in OpenGL and al standard indexed triangle
///   meshes data structures, and not OBJ triplets; YOCTO_OBJ ensures that no
///   vertex dusplication happens thought for same triplets
/// 2. we split shapes on changes to groups and materials, instead of keeping
///   per-face group/material data; this makes the data usable right away in
///   a GPU viewer; this is not a major limitation if we accept the previous
///   point that already changes shapes topology.
///
/// 1. load a obj data with `load_obj()`; can load also textues
/// 2. look at the `obj_XXX` data structures for access to individual elements
/// 3. use obj back to disk with `save_obj()`; can also save textures
/// 4. use get_shape() to get a flattened shape version that contains only
///    triangles, lines or points
///
///
/// ### Khronos glTF
///
/// Khronos GLTF loader and writer for Khronos glTF format. Supports
/// all the glTF spec and the Khronos extensions. All parsing and writing code
/// is autogenerated form the schema. Supports glTF version 2.0 and the
/// following extensions: `KHR_binary_glTF` and `KHR_specular_glossiness`.
///
/// This component depends on `json.hpp` and, for image loading and saving,
/// it depends on `stb_image.h`, `stb_image_write.h`, `stb_image_resize.h` and
/// `tinyexr.h`. This feature can be disabled as before.
///
/// The library provides a low  level interface that is a direct
/// C++ translation of the glTF schemas and should be used if one wants
/// complete control over the fromat or an application wants to have their
/// own scene code added. A higher-level interface is provided by the scene
/// or by `yocto_gltf.h`.
///
/// glTF is a very complex file format and was designed mainly with untyped
/// languages in mind. We attempt to match the glTF low-level interface
/// to C++ as best as it can. Since the code is generated from the schema, we
/// follow glTF naming conventions and typing quite well. To simplify adoption
/// and keep the API relatively simple we use vector as arrays and use
/// pointers to reference to all glTF objects. While this makes it less effcient
/// than it might have been, glTF heavy use of optional values makes this
/// necessary. At the same time, we do not keep track of set/unset values
/// for basic types (int, float, bool) as a compromise for efficieny.
///
/// glTF uses integer indices to access objects.
/// While writing code ourselves we found that we add signiicant problems
/// since we would use an index to access the wriong type of scene objects.
/// For this reasons, we use an explit index `glTFid<T>` that can only access
/// an object of type T. Internally this is just the same old glTF index. But
/// this can used to access the scene data with `glTF::get<T>(index)`.
///
/// 1. load a glTF model with `load_gltf()`
/// 2. look at the `glTFXXX` data structures for access to individual elements
/// 3. save glTF back to disk with `save_gltf()`
///
///
/// ### OpenGL support
///
/// We include a set of utilities to draw on screen with OpenGL 3.3, manage
/// windows with GLFW and draw immediate-mode widgets with ImGui.
///
/// 1. texture and buffer objects with `gl_texture` and `gl_buffer`
///     - create textures/buffers with appropriate constructors
///     - check validity wiht `is_valid()`
///     - update textures/buffers with `update()` functions
///     - delete textures/buffers with `clear()`
///     - bind/unbind textures/buffers with `bind()`/`unbind()`
///     - draw elements with `gl_buffer::draw_elems()`
/// 2. program objects with `gl_program`
///     - program creation with constructor
///     - check validity wiht `is_valid()`
///     - delete with `clear()`
///     - uniforms with `set_program_uniform()`
///     - vertex attrib with `set_program_vertattr()`
///     - draw elements with `gl_buffer::draw_elems()`
/// 3. image viewing with `gl_stdimage_program`, with support for tone mapping.
/// 4. draw surfaces and hair with GGX/Kayjia-Kay with `gl_stdsurface_program`
///     - initialize the program with constructor
///     - check validity wiht `is_valid()`
///     - start/end each frame with `begin_frame()`, `end_frame()`
///     - define lights with `set_lights()`
///     - start/end each shape with `begin_shape()`, `end_shape()`
///     - define material Parameters with `set_material()`
///     - define vertices with `set_vert()`
///     - draw elements with `draw_elems()`
/// 5. also includes other utlities for quick OpenGL
/// 6. GLFW window with `gl_window`
///     - create with constructor
///     - delete with `clear()`
///     - set callbacks with `set_callbacks()`
///     - includes carious utiliies to query window, mouse and keyboard
/// 7. immediate mode widgets
///     - init with `init_widget()`
///     - use the various widget calls to draw the widget and handle events
///
///
/// ### Other Utilities
///
/// We include additional utilities for writing command line applications and
/// manipulating files.
///
/// 1. Python-like string opeations: `startswith()`, `endswith()`, `contains()`,
///    `splitlines()`, `partition()`, `split()`, `splitlines()`, `strip()`,
///    `rstrip()`, `lstrip()`, `join()`, `lower()`, `upper()`, `isspace()`,
///    `replace()`
/// 2. Path-like path operations: `path_dirname()`, `path_extension()`,
///    `path_basename()`, `path_filename()`, `replace_path_extension()`,
///    `prepend_path_extension()`, `split_path()`
/// 3. Python-like format strings (only support for position arguments and no
///    formatting commands): `format()`, `print()`
/// 5. load/save entire files: `load_binfile()`, `load_txtfile()`,
///    `save_binfile()` and `save_binfile()`
/// 4. simple logger with support for console and file streams:
///     1. create a `logger`
///     2. add more streams with `add_console_stream()` or `add_file_stream()`
///     3. write log messages with `log_msg()` and its variants
///     4. you can also use a global default logger with the free functions
///        `log_XXX()`
/// 5. thead pool for concurrent execution (waiting the standard to catch up):
///     1. either create a `thread_pool` or use the global one
///     2. run tasks in parallel `parallel_for()`
///     3. run tasks asynchronously `async()`
/// 6. timer for simple access to `std::chrono`:
///     1. create a `timer`
///     2. start and stop the clock with `start()` and `stop()`
///     3. get time with `elapsed_time()`
///
///
/// ### Command Line Parsing
///
/// The library includes a simple command line parser that parses commands in
/// immediate mode, i.e. when an option is declared. The parser supports options
/// and unnamed arguments with generic types parsed using C++ stream. The
/// parser autogenerates its own documentation. This allows to write complex
/// command lines with a tiny amount of implementation code on both the library
/// and user end.
///
/// 1. create a `cmdline` parser object by passing `argc, argv, name, help`
///     - an option for printing help is automatically added
/// 2. for each option, parse it calling the functions `parse_opt()`
///     - options are parsed on the fly and a comprehensive help is
///       automatically generated
///     - supports bool (flags), int, float, double, string, enums
///     - options names are "--longname" for longname and "-s" for short
///     - command line format is "--longname value", "-s v" for all but flags
///     - values are parsed with `iostream <<` operators
///     - for general use `opt = parse_opt<type>()`
///     - for boolean flags is `parse_flag()`
///     - for enums use `parse_opte()`
/// 3. for each unnamed argument, parse it calling the functions parse_arg()
///     - names are only used for help
///     - supports types as above
///     - for general use `arg = parse_arg<type>()`
///     - to parse all remaining values use `args = parse_arga<type>(...)`
/// 4. end cmdline parsing with `check_parsing()` to check for unsued values,
///    missing arguments
/// 5. to check for error use `should_exit()` and to print the message use
///    `get_message()`
/// 6. since arguments are parsed immediately, one can easily implement
///    subcommands by just branching the command line code based on a read
///    argument without any need for complex syntax
///
///
/// ## History
///
/// Here we mark only major features added to the library. Small refactorings
/// and bug fixes are reported here.
///
/// - v 0.1.0: initial release after refactoring
///

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
// LICENSE OF INCLUDED SOFTWARE for ThreadPool code from LLVM code base
//
// Copyright (c) 2003-2016 University of Illinois at Urbana-Champaign.
// All rights reserved.
//
// Developed by:
//
//     LLVM Team
//
//     University of Illinois at Urbana-Champaign
//
//     http://llvm.org
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// with the Software without restriction, including without limitation the
// rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
// sell copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
//     * Redistributions of source code must retain the above copyright notice,
//       this list of conditions and the following disclaimers.
//
//     * Redistributions in binary form must reproduce the above copyright
//     notice,
//       this list of conditions and the following disclaimers in the
//       documentation and/or other materials provided with the distribution.
//
//     * Neither the names of the LLVM Team, University of Illinois at
//       Urbana-Champaign, nor the names of its contributors may be used to
//       endorse or promote products derived from this Software without specific
//       prior written permission.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
// CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH
// THE SOFTWARE.
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
#define YGL_OPENGL 0
#endif

// enable explicit json objects in glTF
#ifndef YGL_GLTFJSON
#define YGL_GLTFJSON 0
#endif

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <algorithm>
#include <array>
#include <cassert>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <fstream>
#include <functional>
#include <future>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <map>
#include <random>
#include <sstream>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// Compilation option
#define YGL_FAST_RANDFLOAT 1

// -----------------------------------------------------------------------------
// IMPORTED MATH FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

/// sqrt
using std::sqrt;
/// pow
using std::pow;
/// pow
using std::exp;
/// log
using std::log;
/// log10
using std::log10;
/// sin
using std::sin;
/// cos
using std::cos;
/// tan
using std::tan;
/// asin
using std::asin;
/// acos
using std::acos;
/// atan
using std::atan;
/// atan2
using std::atan2;
/// absolute value
using std::abs;
/// floating point absolute value
using std::fabs;
/// floor
using std::floor;
/// ceil
using std::ceil;
/// round
using std::round;
/// isfinate
using std::isfinite;

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPORTED CONTAINERS AND RELATED FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

/// string
using std::string;
/// vector
using std::vector;
/// array
using std::array;
/// map
using std::map;
/// unordered map
using std::unordered_map;
/// unordered set
using std::unordered_set;
/// pair
using std::pair;
/// tuple
using std::tuple;
/// unique pointer
using std::unique_ptr;
/// function
using std::function;
/// string literals
using namespace std::string_literals;
/// numeric limits
using std::numeric_limits;
/// initializer list
using std::initializer_list;
/// output stream
using std::ostream;
/// input stream
using std::istream;
/// string stream
using std::stringstream;
/// file stream
using std::fstream;
/// runtime error
using std::runtime_error;
/// exception
using std::exception;
/// ios base
using std::ios_base;
/// find algorithms
using std::find;
/// swap algorithms
using std::swap;
/// get line from streams
using std::getline;
/// convert to string
using std::to_string;
/// cout object for printing
using std::cout;

// makes literals available
using namespace std::literals;

}  // namespace ygl

// -----------------------------------------------------------------------------
// BASIC TYPEDEFS, MATH CONSTANTS AND FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

/// convenient typedef for bytes
using byte = unsigned char;

/// convenient typedef for bytes
using uint = unsigned int;

/// pi (float)
const auto pif = 3.14159265f;
/// pi (double)
const auto pi = 3.1415926535897932384626433832795;

/// shortcat for float max value
const auto flt_max = numeric_limits<float>::max();
/// shortcat for float min value
const auto flt_min = numeric_limits<float>::lowest();
/// shortcat for float epsilon
const auto flt_eps = numeric_limits<float>::epsilon();
/// shortcat for int max value
const auto int_max = numeric_limits<int>::max();
/// shortcat for int min value
const auto int_min = numeric_limits<int>::min();

/// Safe minimum value.
inline int min(int x, int y) { return (x < y) ? x : y; }
/// Safe minimum value.
inline float min(float x, float y) { return (x < y) ? x : y; }
/// Safe minimum value.
inline int min(initializer_list<int> vs) {
    auto m = int_max;
    for (auto v : vs) m = min(m, v);
    return m;
}
/// Safe minimum value.
inline float min(initializer_list<float> vs) {
    auto m = flt_max;
    for (auto v : vs) m = min(m, v);
    return m;
}

/// Safe maximum value.
inline int max(int x, int y) { return (x > y) ? x : y; }
/// Safe maximum value.
inline float max(float x, float y) { return (x > y) ? x : y; }
/// Safe maximum value.
inline int max(initializer_list<int> vs) {
    auto m = int_min;
    for (auto v : vs) m = max(m, v);
    return m;
}
/// Safe maximum value.
inline float max(initializer_list<float> vs) {
    auto m = flt_min;
    for (auto v : vs) m = max(m, v);
    return m;
}

/// Clamp a value between a minimum and a maximum.
inline int clamp(int x, int min_, int max_) { return min(max(x, min_), max_); }
/// Clamp a value between a minimum and a maximum.
inline float clamp(float x, float min_, float max_) {
    return min(max(x, min_), max_);
}

/// Linear interpolation.
inline float lerp(float a, float b, float t) { return a + (b - a) * t; }
/// bilinear interpolation
inline float bilerp(float aa, float ba, float ab, float bb, float s, float t) {
    return aa * (1 - s) * (1 - t) + ba * s * (1 - t) + ab * (1 - s) * t +
           bb * s * t;
}

/// Integer power of two
inline int pow2(int x) { return 1 << x; }

/// Fast floor
inline int fastfloor(float x) {
    auto xi = (int)x;
    return (x < xi) ? xi - 1 : xi;
}
/// Safe float to byte conversion
inline byte float_to_byte(float x) {
    return (byte)max(0, min(int(x * 256), 255));
}

/// Safe byte to float conversion
inline float byte_to_float(byte x) { return (float)x / 255.0f; }

}  // namespace ygl

// -----------------------------------------------------------------------------
// VECTORS
// -----------------------------------------------------------------------------
namespace ygl {

/// Vector of 2 float elements.
struct vec2f {
    /// default constructor
    vec2f() : x{0}, y{0} {}
    /// element constructor
    explicit vec2f(float vv) : x(vv), y(vv) {}
    /// element constructor
    vec2f(float x, float y) : x{x}, y{y} {}

    /// element access
    float& operator[](int i) { return (&x)[i]; }
    /// element access
    const float& operator[](int i) const { return (&x)[i]; }

    /// data access
    float* data() { return &x; }
    /// data access
    const float* data() const { return &x; }

    /// element data
    float x;
    /// element data
    float y;
};

/// Vector of 3 float elements.
struct vec3f {
    /// default constructor
    vec3f() : x{0}, y{0}, z{0} {}
    /// element constructor
    explicit vec3f(float vv) : x(vv), y(vv), z(vv) {}
    /// element constructor
    vec3f(float x, float y, float z) : x{x}, y{y}, z{z} {}

    /// element access
    float& operator[](int i) { return (&x)[i]; }
    /// element access
    const float& operator[](int i) const { return (&x)[i]; }

    /// data access
    float* data() { return &x; }
    /// data access
    const float* data() const { return &x; }

    /// element data
    float x;
    /// element data
    float y;
    /// element data
    float z;
};

/// Vector of 4 float elements.
struct vec4f {
    /// default constructor
    vec4f() : x{0}, y{0}, z{0}, w{0} {}
    /// element constructor
    explicit vec4f(float vv) : x(vv), y(vv), z(vv), w(vv) {}
    /// element constructor
    vec4f(float x, float y, float z, float w) : x{x}, y{y}, z{z}, w{w} {}
    /// constructor from smaller vector
    vec4f(const vec3f& xyz, float w) : x{xyz.x}, y{xyz.y}, z{xyz.z}, w{w} {}

    /// element access
    float& operator[](int i) { return (&x)[i]; }
    /// element access
    const float& operator[](int i) const { return (&x)[i]; }

    /// data access
    float* data() { return &x; }
    /// data access
    const float* data() const { return &x; }

    /// access xyz components
    vec3f& xyz() { return *(vec3f*)&x; }
    /// access xyz components
    const vec3f& xyz() const { return *(vec3f*)&x; }

    /// element data
    float x;
    /// element data
    float y;
    /// element data
    float z;
    /// element data
    float w;
};

/// Vector of 2 int elements.
struct vec2i {
    /// default constructor
    vec2i() : x{0}, y{0} {}
    /// element constructor
    explicit vec2i(int vv) : x(vv), y(vv) {}
    /// element constructor
    vec2i(int x, int y) : x{x}, y{y} {}

    /// element access
    int& operator[](int i) { return (&x)[i]; }
    /// element access
    const int& operator[](int i) const { return (&x)[i]; }

    /// data access
    int* data() { return &x; }
    /// data access
    const int* data() const { return &x; }

    /// element data
    int x;
    /// element data
    int y;
};

/// Vector of 3 int elements.
struct vec3i {
    /// default constructor
    vec3i() : x{0}, y{0}, z{0} {}
    /// element constructor
    explicit vec3i(int vv) : x(vv), y(vv), z(vv) {}
    /// element constructor
    vec3i(int x, int y, int z) : x{x}, y{y}, z{z} {}

    /// element access
    int& operator[](int i) { return (&x)[i]; }
    /// element access
    const int& operator[](int i) const { return (&x)[i]; }

    /// data access
    int* data() { return &x; }
    /// data access
    const int* data() const { return &x; }

    /// element data
    int x;
    /// element data
    int y;
    /// element data
    int z;
};

/// Vector of 4 int elements.
struct vec4i {
    /// default constructor
    vec4i() : x{0}, y{0}, z{0}, w{0} {}
    /// element constructor
    explicit vec4i(int vv) : x(vv), y(vv), z(vv), w(vv) {}
    /// element constructor
    vec4i(int x, int y, int z, int w) : x{x}, y{y}, z{z}, w{w} {}
    /// constructor from smaller vector
    vec4i(const vec3i& xyz, int w) : x{xyz.x}, y{xyz.y}, z{xyz.z}, w{w} {}

    /// element access
    int& operator[](int i) { return (&x)[i]; }
    /// element access
    const int& operator[](int i) const { return (&x)[i]; }

    /// data access
    int* data() { return &x; }
    /// data access
    const int* data() const { return &x; }

    /// access xyz components
    vec3i& xyz() { return *(vec3i*)&x; }
    /// access xyz components
    const vec3i& xyz() const { return *(vec3i*)&x; }

    /// element data
    int x;
    /// element data
    int y;
    /// element data
    int z;
    /// element data
    int w;
};

/// Vector of 3 byte elements.
struct vec3b {
    /// default constructor
    vec3b() : x{0}, y{0}, z{0} {}
    /// element constructor
    explicit vec3b(int vv) : x(vv), y(vv), z(vv) {}
    /// element constructor
    vec3b(byte x, byte y, byte z) : x{x}, y{y}, z{z} {}

    /// element access
    byte& operator[](int i) { return (&x)[i]; }
    /// element access
    const byte& operator[](int i) const { return (&x)[i]; }

    /// data access
    byte* data() { return &x; }
    /// data access
    const byte* data() const { return &x; }

    /// element data
    byte x;
    /// element data
    byte y;
    /// element data
    byte z;
};

/// Vector of 4 byte elements.
struct vec4b {
    /// default constructor
    vec4b() : x{0}, y{0}, z{0}, w{0} {}
    /// element constructor
    explicit vec4b(byte vv) : x(vv), y(vv), z(vv), w(vv) {}
    /// element constructor
    vec4b(byte x, byte y, byte z, byte w) : x{x}, y{y}, z{z}, w{w} {}
    /// constructor from smaller vector
    vec4b(const vec3b& xyz, byte w) : x{xyz.x}, y{xyz.y}, z{xyz.z}, w{w} {}

    /// element access
    byte& operator[](int i) { return (&x)[i]; }
    /// element access
    const byte& operator[](int i) const { return (&x)[i]; }

    /// data access
    byte* data() { return &x; }
    /// data access
    const byte* data() const { return &x; }

    /// access xyz components
    vec3b& xyz() { return *(vec3b*)&x; }
    /// access xyz components
    const vec3b& xyz() const { return *(vec3b*)&x; }

    /// element data
    byte x;
    /// element data
    byte y;
    /// element data
    byte z;
    /// element data
    byte w;
};

/// 2-dimensional float zero vector
const auto zero2f = vec2f();
/// 3-dimensional float zero vector
const auto zero3f = vec3f();
/// 4-dimensional float zero vector
const auto zero4f = vec4f();

/// 2-dimensional int zero vector
const auto zero2i = vec2i();
/// 3-dimensional int zero vector
const auto zero3i = vec3i();
/// 4-dimensional int zero vector
const auto zero4i = vec4i();

/// 4-dimensional byte zero vector
const auto zero4b = vec4b();

/// iteration support
inline int* begin(vec2i& a) { return &a.x; }
/// iteration support
inline const int* begin(const vec2i& a) { return &a.x; }
/// iteration support
inline int* end(vec2i& a) { return &a.x + 2; }
/// iteration support
inline const int* end(const vec2i& a) { return &a.x + 2; }

/// iteration support
inline int* begin(vec3i& a) { return &a.x; }
/// iteration support
inline const int* begin(const vec3i& a) { return &a.x; }
/// iteration support
inline int* end(vec3i& a) { return &a.x + 3; }
/// iteration support
inline const int* end(const vec3i& a) { return &a.x + 3; }

/// iteration support
inline int* begin(vec4i& a) { return &a.x; }
/// iteration support
inline const int* begin(const vec4i& a) { return &a.x; }
/// iteration support
inline int* end(vec4i& a) { return &a.x + 4; }
/// iteration support
inline const int* end(const vec4i& a) { return &a.x + 4; }

/// vector operator ==
inline bool operator==(const vec2f& a, const vec2f& b) {
    return a.x == b.x && a.y == b.y;
}
/// vector operator !=
inline bool operator!=(const vec2f& a, const vec2f& b) {
    return a.x != b.x || a.y != b.y;
}

/// vector operator ==
inline bool operator==(const vec3f& a, const vec3f& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}
/// vector operator !=
inline bool operator!=(const vec3f& a, const vec3f& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z;
}

/// vector operator ==
inline bool operator==(const vec4f& a, const vec4f& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}
/// vector operator !=
inline bool operator!=(const vec4f& a, const vec4f& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w;
}

/// vector operator ==
inline bool operator==(const vec2i& a, const vec2i& b) {
    return a.x == b.x && a.y == b.y;
}
/// vector operator !=
inline bool operator!=(const vec2i& a, const vec2i& b) {
    return a.x != b.x || a.y != b.y;
}

/// vector operator ==
inline bool operator==(const vec3i& a, const vec3i& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}
/// vector operator !=
inline bool operator!=(const vec3i& a, const vec3i& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z;
}

/// vector operator ==
inline bool operator==(const vec4i& a, const vec4i& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}
/// vector operator !=
inline bool operator!=(const vec4i& a, const vec4i& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w;
}

/// vector operator < (lexicographic order - useful for map)
inline bool operator<(const vec2i& a, const vec2i& b) {
    for (auto i = 0; i < 2; i++) {
        if (a[i] < b[i]) return true;
        if (a[i] > b[i]) return false;
    }
    return false;
}
/// vector operator < (lexicographic order - useful for map)
inline bool operator<(const vec3i& a, const vec3i& b) {
    for (auto i = 0; i < 3; i++) {
        if (a[i] < b[i]) return true;
        if (a[i] > b[i]) return false;
    }
    return false;
}
/// vector operator < (lexicographic order - useful for map)
inline bool operator<(const vec4i& a, const vec4i& b) {
    for (auto i = 0; i < 4; i++) {
        if (a[i] < b[i]) return true;
        if (a[i] > b[i]) return false;
    }
    return false;
}

/// vector operator +
inline vec2f operator+(const vec2f& a) { return a; }
/// vector operator -
inline vec2f operator-(const vec2f& a) { return {-a.x, -a.y}; }
/// vector operator +
inline vec2f operator+(const vec2f& a, const vec2f& b) {
    return {a.x + b.x, a.y + b.y};
}
/// vector operator -
inline vec2f operator-(const vec2f& a, const vec2f& b) {
    return {a.x - b.x, a.y - b.y};
}
/// vector operator *
inline vec2f operator*(const vec2f& a, const vec2f& b) {
    return {a.x * b.x, a.y * b.y};
}
/// vector operator *
inline vec2f operator*(const vec2f& a, float b) { return {a.x * b, a.y * b}; }
/// vector operator *
inline vec2f operator*(float a, const vec2f& b) { return {a * b.x, a * b.y}; }
/// vector operator /
inline vec2f operator/(const vec2f& a, const vec2f& b) {
    return {a.x / b.x, a.y / b.y};
}
/// vector operator /
inline vec2f operator/(const vec2f& a, float b) { return {a.x / b, a.y / b}; }
/// vector operator /
inline vec2f operator/(float a, const vec2f& b) { return {a / b.x, a / b.y}; }

/// vector operator +
inline vec3f operator+(const vec3f& a) { return a; }
/// vector operator -
inline vec3f operator-(const vec3f& a) { return {-a.x, -a.y, -a.z}; }
/// vector operator +
inline vec3f operator+(const vec3f& a, const vec3f& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
/// vector operator -
inline vec3f operator-(const vec3f& a, const vec3f& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
/// vector operator *
inline vec3f operator*(const vec3f& a, const vec3f& b) {
    return {a.x * b.x, a.y * b.y, a.z * b.z};
}
/// vector operator *
inline vec3f operator*(const vec3f& a, float b) {
    return {a.x * b, a.y * b, a.z * b};
}
/// vector operator *
inline vec3f operator*(float a, const vec3f& b) {
    return {a * b.x, a * b.y, a * b.z};
}
/// vector operator /
inline vec3f operator/(const vec3f& a, const vec3f& b) {
    return {a.x / b.x, a.y / b.y, a.z / b.z};
}
/// vector operator /
inline vec3f operator/(const vec3f& a, float b) {
    return {a.x / b, a.y / b, a.z / b};
}
/// vector operator /
inline vec3f operator/(float a, const vec3f& b) {
    return {a / b.x, a / b.y, a / b.z};
}

/// vector operator +
inline vec4f operator+(const vec4f& a) { return a; }
/// vector operator -
inline vec4f operator-(const vec4f& a) { return {-a.x, -a.y, -a.z, -a.w}; }
/// vector operator +
inline vec4f operator+(const vec4f& a, const vec4f& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}
/// vector operator -
inline vec4f operator-(const vec4f& a, const vec4f& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w};
}
/// vector operator *
inline vec4f operator*(const vec4f& a, const vec4f& b) {
    return {a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w};
}
/// vector operator *
inline vec4f operator*(const vec4f& a, float b) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
}
/// vector operator *
inline vec4f operator*(float a, const vec4f& b) {
    return {a * b.x, a * b.y, a * b.z, a * b.w};
}
/// vector operator /
inline vec4f operator/(const vec4f& a, const vec4f& b) {
    return {a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w};
}
/// vector operator /
inline vec4f operator/(const vec4f& a, float b) {
    return {a.x / b, a.y / b, a.z / b, a.w / b};
}
/// vector operator /
inline vec4f operator/(float a, const vec4f& b) {
    return {a / b.x, a / b.y, a / b.z, a / b.w};
}

/// vector operator +=
inline vec2f& operator+=(vec2f& a, const vec2f& b) { return a = a + b; }
/// vector operator -=
inline vec2f& operator-=(vec2f& a, const vec2f& b) { return a = a - b; }
/// vector operator *=
inline vec2f& operator*=(vec2f& a, const vec2f& b) { return a = a * b; }
/// vector operator *=
inline vec2f& operator*=(vec2f& a, float b) { return a = a * b; }
/// vector operator /=
inline vec2f& operator/=(vec2f& a, const vec2f& b) { return a = a / b; }
/// vector operator /=
inline vec2f& operator/=(vec2f& a, float b) { return a = a / b; }

/// vector operator +=
inline vec3f& operator+=(vec3f& a, const vec3f& b) { return a = a + b; }
/// vector operator -=
inline vec3f& operator-=(vec3f& a, const vec3f& b) { return a = a - b; }
/// vector operator *=
inline vec3f& operator*=(vec3f& a, const vec3f& b) { return a = a * b; }
/// vector operator *=
inline vec3f& operator*=(vec3f& a, float b) { return a = a * b; }
/// vector operator /=
inline vec3f& operator/=(vec3f& a, const vec3f& b) { return a = a / b; }
/// vector operator /=
inline vec3f& operator/=(vec3f& a, float b) { return a = a / b; }

/// vector operator +=
inline vec4f& operator+=(vec4f& a, const vec4f& b) { return a = a + b; }
/// vector operator -=
inline vec4f& operator-=(vec4f& a, const vec4f& b) { return a = a - b; }
/// vector operator *=
inline vec4f& operator*=(vec4f& a, const vec4f& b) { return a = a * b; }
/// vector operator *=
inline vec4f& operator*=(vec4f& a, float b) { return a = a * b; }
/// vector operator /=
inline vec4f& operator/=(vec4f& a, const vec4f& b) { return a = a / b; }
/// vector operator /=
inline vec4f& operator/=(vec4f& a, float b) { return a = a / b; }

/// vector operator +
inline vec2i operator+(const vec2i& a) { return a; }
/// vector operator -
inline vec2i operator-(const vec2i& a) { return {-a.x, -a.y}; }
/// vector operator +
inline vec2i operator+(const vec2i& a, const vec2i& b) {
    return {a.x + b.x, a.y + b.y};
}
/// vector operator -
inline vec2i operator-(const vec2i& a, const vec2i& b) {
    return {a.x - b.x, a.y - b.y};
}

/// vector operator +
inline vec3i operator+(const vec3i& a) { return a; }
/// vector operator -
inline vec3i operator-(const vec3i& a) { return {-a.x, -a.y, -a.z}; }
/// vector operator +
inline vec3i operator+(const vec3i& a, const vec3i& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
/// vector operator -
inline vec3i operator-(const vec3i& a, const vec3i& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

/// vector operator +
inline vec4i operator+(const vec4i& a) { return a; }
/// vector operator -
inline vec4i operator-(const vec4i& a) { return {-a.x, -a.y, -a.z, -a.w}; }
/// vector operator +
inline vec4i operator+(const vec4i& a, const vec4i& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}
/// vector operator -
inline vec4i operator-(const vec4i& a, const vec4i& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w};
}

/// vector operator +=
inline vec2i& operator+=(vec2i& a, const vec2i& b) { return a = a + b; }
/// vector operator -=
inline vec2i& operator-=(vec2i& a, const vec2i& b) { return a = a - b; }

/// vector operator +=
inline vec3i& operator+=(vec3i& a, const vec3i& b) { return a = a + b; }
/// vector operator -=
inline vec3i& operator-=(vec3i& a, const vec3i& b) { return a = a - b; }

/// vector operator +=
inline vec4i& operator+=(vec4i& a, const vec4i& b) { return a = a + b; }
/// vector operator -=
inline vec4i& operator-=(vec4i& a, const vec4i& b) { return a = a - b; }

/// vector dot product
inline float dot(const vec2f& a, const vec2f& b) {
    return a.x * b.x + a.y * b.y;
}
/// vector dot product
inline float dot(const vec3f& a, const vec3f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
/// vector dot product
inline float dot(const vec4f& a, const vec4f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

/// vector cross product
inline float cross(const vec2f& a, const vec2f& b) {
    return a.x * b.y - a.y * b.x;
}
/// vector cross product
inline vec3f cross(const vec3f& a, const vec3f& b) {
    return {
        a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

/// vector length
inline float length(const vec2f& a) { return sqrt(dot(a, a)); }
/// vector length
inline float length(const vec3f& a) { return sqrt(dot(a, a)); }
/// vector length
inline float length(const vec4f& a) { return sqrt(dot(a, a)); }

/// vector normalization
inline vec2f normalize(const vec2f& a) {
    auto l = length(a);
    if (l == 0) return a;
    return a * (1 / l);
}
/// vector normalization
inline vec3f normalize(const vec3f& a) {
    auto l = length(a);
    if (l == 0) return a;
    return a * (1 / l);
}
/// vector normalization
inline vec4f normalize(const vec4f& a) {
    auto l = length(a);
    if (l == 0) return a;
    return a * (1 / l);
}

/// angle between normalized vectors
inline float uangle(const vec3f& a, const vec3f& b) {
    auto d = dot(a, b);
    return d > 1 ? 0 : acos(d < -1 ? -1 : d);
}

/// angle between normalized vectors
inline float uangle(const vec4f& a, const vec4f& b) {
    auto d = dot(a, b);
    return d > 1 ? 0 : acos(d < -1 ? -1 : d);
}

/// angle between vectors
inline float angle(const vec3f& a, const vec3f& b) {
    return uangle(normalize(a), normalize(b));
}

/// vector linear interpolation
inline vec2f lerp(const vec2f& a, const vec2f& b, float t) {
    return a * (1 - t) + b * t;
}
/// vector linear interpolation
inline vec3f lerp(const vec3f& a, const vec3f& b, float t) {
    return a * (1 - t) + b * t;
}
/// vector linear interpolation
inline vec4f lerp(const vec4f& a, const vec4f& b, float t) {
    return a * (1 - t) + b * t;
}

/// vector bilinear interpolation
inline vec3f bilerp(const vec3f& aa, const vec3f& ba, const vec3f& ab,
    const vec3f& bb, float s, float t) {
    return aa * (1 - s) * (1 - t) + ba * s * (1 - t) + ab * (1 - s) * t +
           bb * s * t;
}

/// vector normalized linear interpolation
inline vec3f nlerp(const vec3f& a, const vec3f& b, float t) {
    return normalize(lerp(a, b, t));
}

/// vector spherical linear interpolation (vectors have to be normalized)
inline vec3f slerp(const vec3f& a, const vec3f& b, float t) {
    auto th = uangle(a, b);
    return th == 0 ?
               a :
               a * (sin(th * (1 - t)) / sin(th)) + b * (sin(th * t) / sin(th));
}

/// vector normalized linear interpolation
inline vec4f nlerp(const vec4f& a, const vec4f& b, float t) {
    return normalize(lerp(a, b, t));
}

/// vector spherical linear interpolation (vectors have to be normalized)
inline vec4f slerp(const vec4f& a, const vec4f& b, float t) {
    auto th = uangle(a, b);
    return th == 0 ?
               a :
               a * (sin(th * (1 - t)) / sin(th)) + b * (sin(th * t) / sin(th));
}

/// orthogonal vector
inline vec3f orthogonal(const vec3f& v) {
    // http://lolengine.net/blog/2013/09/21/picking-orthogonal-vector-combing-coconuts)
    return abs(v.x) > abs(v.z) ? vec3f{-v.y, v.x, 0} : vec3f{0, -v.z, v.y};
}

/// orthonormalize two vectors
inline vec3f orthonormalize(const vec3f& a, const vec3f& b) {
    return normalize(a - b * dot(a, b));
}

/// vector component-wise clamp
inline vec2f clamp(const vec2f& x, float min, float max) {
    return {clamp(x.x, min, max), clamp(x.y, min, max)};
}
/// vector component-wise clamp
inline vec3f clamp(const vec3f& x, float min, float max) {
    return {clamp(x.x, min, max), clamp(x.y, min, max), clamp(x.z, min, max)};
}
/// vector component-wise clamp
inline vec4f clamp(const vec4f& x, float min, float max) {
    return {clamp(x.x, min, max), clamp(x.y, min, max), clamp(x.z, min, max),
        clamp(x.w, min, max)};
}

/// clamp the length of a vector
inline vec2f clamplen(const vec2f& x, float max) {
    auto l = length(x);
    return (l > max) ? x * max / l : x;
}
/// clamp the length of a vector
inline vec3f clamplen(const vec3f& x, float max) {
    auto l = length(x);
    return (l > max) ? x * max / l : x;
}
/// clamp the length of a vector
inline vec4f clamplen(const vec4f& x, float max) {
    auto l = length(x);
    return (l > max) ? x * max / l : x;
}

// implementation
inline pair<int, float> _min_element(int N, const float* a) {
    auto v = flt_max;
    auto pos = -1;
    for (auto i = 0; i < N; i++) {
        if (v > a[i]) {
            v = a[i];
            pos = i;
        }
    }
    return {pos, v};
}

/// min vector element
inline pair<int, float> min_element(const vec2f& a) {
    return _min_element(2, &a.x);
}
/// min vector element
inline pair<int, float> min_element(const vec3f& a) {
    return _min_element(3, &a.x);
}
/// min vector element
inline pair<int, float> min_element(const vec4f& a) {
    return _min_element(4, &a.x);
}

/// index of the max vector element
inline pair<int, float> _max_element(int N, const float* a) {
    auto v = flt_min;
    auto pos = -1;
    for (auto i = 0; i < N; i++) {
        if (v < a[i]) {
            v = a[i];
            pos = i;
        }
    }
    return {pos, v};
}

/// index of the min vector element
inline pair<int, float> max_element(const vec2f& a) {
    return _max_element(2, &a.x);
}
/// index of the min vector element
inline pair<int, float> max_element(const vec3f& a) {
    return _max_element(3, &a.x);
}
/// index of the min vector element
inline pair<int, float> max_element(const vec4f& a) {
    return _max_element(4, &a.x);
}

/// Element-wise conversion
inline vec3b float_to_byte(const vec3f& a) {
    return {float_to_byte(a.x), float_to_byte(a.y), float_to_byte(a.z)};
}
/// Element-wise conversion
inline vec3f byte_to_float(const vec3b& a) {
    return {byte_to_float(a.x), byte_to_float(a.y), byte_to_float(a.z)};
}
/// Element-wise conversion
inline vec4b float_to_byte(const vec4f& a) {
    return {float_to_byte(a.x), float_to_byte(a.y), float_to_byte(a.z),
        float_to_byte(a.w)};
}
/// Element-wise conversion
inline vec4f byte_to_float(const vec4b& a) {
    return {byte_to_float(a.x), byte_to_float(a.y), byte_to_float(a.z),
        byte_to_float(a.w)};
}

/// stream write
inline ostream& operator<<(ostream& os, const vec2f& a) {
    return os << a.x << ' ' << a.y;
}
/// stream write
inline ostream& operator<<(ostream& os, const vec3f& a) {
    return os << a.x << ' ' << a.y << ' ' << a.z;
}
/// stream write
inline ostream& operator<<(ostream& os, const vec4f& a) {
    return os << a.x << ' ' << a.y << ' ' << a.z << ' ' << a.w;
}
/// stream write
inline ostream& operator<<(ostream& os, const vec2i& a) {
    return os << a.x << ' ' << a.y;
}
/// stream write
inline ostream& operator<<(ostream& os, const vec3i& a) {
    return os << a.x << ' ' << a.y << ' ' << a.z;
}
/// stream write
inline ostream& operator<<(ostream& os, const vec4i& a) {
    return os << a.x << ' ' << a.y << ' ' << a.z << ' ' << a.w;
}

/// stream read
inline istream& operator>>(istream& is, vec2f& a) { return is >> a.x >> a.y; }
/// stream read
inline istream& operator>>(istream& is, vec3f& a) {
    return is >> a.x >> a.y >> a.z;
}
/// stream read
inline istream& operator>>(istream& is, vec4f& a) {
    return is >> a.x >> a.y >> a.z >> a.w;
}
/// stream read
inline istream& operator>>(istream& is, vec2i& a) { return is >> a.x >> a.y; }
/// stream read
inline istream& operator>>(istream& is, vec3i& a) {
    return is >> a.x >> a.y >> a.z;
}
/// stream read
inline istream& operator>>(istream& is, vec4i& a) {
    return is >> a.x >> a.y >> a.z >> a.w;
}

}  // namespace ygl

namespace std {
/// Hash functor for vector for use with unordered_map
template <>
struct hash<ygl::vec2i> {
    // from boost::hash_combine
    static size_t hash_combine(size_t h, size_t h1) {
        h ^= h1 + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
    size_t operator()(const ygl::vec2i& v) const {
        auto vh = hash<int>();
        auto h = (size_t)0;
        for (auto i = 0; i < 2; i++) h = hash_combine(h, vh(v[i]));
        return h;
    }
};
/// Hash functor for vector for use with unordered_map
template <>
struct hash<ygl::vec3i> {
    // from boost::hash_combine
    static size_t hash_combine(size_t h, size_t h1) {
        h ^= h1 + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
    size_t operator()(const ygl::vec3i& v) const {
        auto vh = hash<int>();
        auto h = (size_t)0;
        for (auto i = 0; i < 3; i++) h = hash_combine(h, vh(v[i]));
        return h;
    }
};
/// Hash functor for vector for use with unordered_map
template <>
struct hash<ygl::vec4i> {
    // from boost::hash_combine
    static size_t hash_combine(size_t h, size_t h1) {
        h ^= h1 + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
    size_t operator()(const ygl::vec4i& v) const {
        auto vh = hash<int>();
        auto h = (size_t)0;
        for (auto i = 0; i < 4; i++) h = hash_combine(h, vh(v[i]));
        return h;
    }
};
}  // namespace std

// -----------------------------------------------------------------------------
// MATRICES
// -----------------------------------------------------------------------------
namespace ygl {

/// Matrix of 2x2 elements stored in column major format.
/// Colums access via operator[].
struct mat2f {
    /// default constructor
    mat2f() : x{1, 0}, y{0, 1} {}
    /// diagonal constructor
    explicit mat2f(float vv) : x{vv, 0}, y{0, vv} {}
    /// list constructor
    mat2f(const vec2f& x, const vec2f& y) : x(x), y(y) {}

    /// element access
    vec2f& operator[](int i) { return (&x)[i]; }
    /// element access
    const vec2f& operator[](int i) const { return (&x)[i]; }

    /// data access
    vec2f* data() { return &x; }
    /// data access
    const vec2f* data() const { return &x; }

    /// element data
    vec2f x;
    /// element data
    vec2f y;
};

/// Matrix of 3x3 elements stored in column major format.
/// Colums access via operator[].
struct mat3f {
    /// default constructor
    mat3f() : x{1, 0, 0}, y{0, 1, 0}, z{0, 0, 1} {}
    /// diagonal constructor
    explicit mat3f(float vv) : x{vv, 0, 0}, y{0, vv, 0}, z{0, 0, vv} {}
    /// list constructor
    mat3f(const vec3f& x, const vec3f& y, const vec3f& z) : x(x), y(y), z(z) {}

    /// element access
    vec3f& operator[](int i) { return (&x)[i]; }
    /// element access
    const vec3f& operator[](int i) const { return (&x)[i]; }

    /// data access
    vec3f* data() { return &x; }
    /// data access
    const vec3f* data() const { return &x; }

    /// element data
    vec3f x;
    /// element data
    vec3f y;
    /// element data
    vec3f z;
};

/// Matrix of 4x4 elements stored in column major format.
/// Colums access via operator[].
struct mat4f {
    /// default constructor
    mat4f() : x{1, 0, 0, 0}, y{0, 1, 0, 0}, z{0, 0, 1, 0}, w{0, 0, 0, 1} {}
    /// diagonal constructor
    explicit mat4f(float vv)
        : x{vv, 0, 0, 0}, y{0, vv, 0, 0}, z{0, 0, vv, 0}, w{0, 0, 0, vv} {}
    /// list constructor
    mat4f(const vec4f& x, const vec4f& y, const vec4f& z, const vec4f& w)
        : x(x), y(y), z(z), w(w) {}

    /// element access
    vec4f& operator[](int i) { return (&x)[i]; }
    /// element access
    const vec4f& operator[](int i) const { return (&x)[i]; }

    /// data access
    vec4f* data() { return &x; }
    /// data access
    const vec4f* data() const { return &x; }

    /// element data
    vec4f x;
    /// element data
    vec4f y;
    /// element data
    vec4f z;
    /// element data
    vec4f w;
};

/// 2-dimensional float identity matrix
const auto identity_mat2f = mat2f();
/// 3-dimensional float identity matrix
const auto identity_mat3f = mat3f();
/// 4-dimensional float identity matrix
const auto identity_mat4f = mat4f();

/// matrix operator ==
inline bool operator==(const mat2f& a, const mat2f& b) {
    return a.x == b.x && a.y == b.y;
}
/// matrix operator !=
inline bool operator!=(const mat2f& a, const mat2f& b) { return !(a == b); }
/// matrix operator ==
inline bool operator==(const mat3f& a, const mat3f& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}
/// matrix operator !=
inline bool operator!=(const mat3f& a, const mat3f& b) { return !(a == b); }
/// matrix operator ==
inline bool operator==(const mat4f& a, const mat4f& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}
/// matrix operator !=
inline bool operator!=(const mat4f& a, const mat4f& b) { return !(a == b); }

/// matrix operator +
inline mat2f operator+(const mat2f& a, const mat2f& b) {
    return {a.x + b.x, a.y + b.y};
}
/// matrix operator +
inline mat3f operator+(const mat3f& a, const mat3f& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
/// matrix operator +
inline mat4f operator+(const mat4f& a, const mat4f& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}

/// matrix scalar multiply
inline mat2f operator*(const mat2f& a, float b) { return {a.x * b, a.y * b}; }

/// matrix scalar division
inline mat2f operator/(const mat2f& a, float b) { return {a.x / b, a.y / b}; }

/// matrix-vector right multiply
inline vec2f operator*(const mat2f& a, const vec2f& b) {
    return a.x * b.x + a.y * b.y;
}

/// matrix-vector left multiply
inline vec2f operator*(const vec2f& a, const mat2f& b) {
    return {dot(a, b.x), dot(a, b.y)};
}

/// matrix-matrix multiply
inline mat2f operator*(const mat2f& a, const mat2f& b) {
    return {a * b.x, a * b.y};
}

/// matrix scalar multiply
inline mat3f operator*(const mat3f& a, float b) {
    return {a.x * b, a.y * b, a.z * b};
}

/// matrix scalar division
inline mat3f operator/(const mat3f& a, float b) {
    return {a.x / b, a.y / b, a.z / b};
}

/// matrix scalar multiply
inline mat4f operator*(const mat4f& a, float b) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
}

/// matrix scalar division
inline mat4f operator/(const mat4f& a, float b) {
    return {a.x / b, a.y / b, a.z / b, a.w / b};
}

/// matrix-vector right multiply
inline vec3f operator*(const mat3f& a, const vec3f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/// matrix-vector left multiply
inline vec3f operator*(const vec3f& a, const mat3f& b) {
    return {dot(a, b.x), dot(a, b.y), dot(a, b.z)};
}

/// matrix-matrix multiply
inline mat3f operator*(const mat3f& a, const mat3f& b) {
    return {a * b.x, a * b.y, a * b.z};
}

/// matrix-vector right multiply
inline vec4f operator*(const mat4f& a, const vec4f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

/// matrix-vector left multiply
inline vec4f operator*(const vec4f& a, const mat4f& b) {
    return {dot(a, b.x), dot(a, b.y), dot(a, b.z), dot(a, b.w)};
}

/// matrix-matrix multiply
inline mat4f operator*(const mat4f& a, const mat4f& b) {
    return {a * b.x, a * b.y, a * b.z, a * b.w};
}

/// matrix sum assignment
inline mat2f& operator+=(mat2f& a, const mat2f& b) { return a = a + b; }

/// matrix-matrix multiply assignment
inline mat2f& operator*=(mat2f& a, const mat2f& b) { return a = a * b; }

/// matrix scaling assignment
inline mat2f& operator*=(mat2f& a, float b) { return a = a * b; }

/// matrix scaling assignment
inline mat2f& operator/=(mat2f& a, float b) { return a = a / b; }

/// matrix sum assignment
inline mat3f& operator+=(mat3f& a, const mat3f& b) { return a = a + b; }

/// matrix-matrix multiply assignment
inline mat3f& operator*=(mat3f& a, const mat3f& b) { return a = a * b; }

/// matrix scaling assignment
inline mat3f& operator*=(mat3f& a, float b) { return a = a * b; }

/// matrix scaling assignment
inline mat3f& operator/=(mat3f& a, float b) { return a = a / b; }

/// matrix sum assignment
inline mat4f& operator+=(mat4f& a, const mat4f& b) { return a = a + b; }

/// matrix-matrix multiply assignment
inline mat4f& operator*=(mat4f& a, const mat4f& b) { return a = a * b; }

/// matrix scaling assignment
inline mat4f& operator*=(mat4f& a, float b) { return a = a * b; }

/// matrix scaling assignment
inline mat4f& operator/=(mat4f& a, float b) { return a = a / b; }

/// matrix diagonal
inline vec2f mat_diagonal(const mat2f& a) { return {a.x.x, a.y.y}; }
/// matrix diagonal
inline vec3f mat_diagonal(const mat3f& a) { return {a.x.x, a.y.y, a.z.z}; }
/// matrix diagonal
inline vec4f mat_diagonal(const mat4f& a) {
    return {a.x.x, a.y.y, a.z.z, a.w.w};
}

/// matrix transpose
inline mat2f transpose(const mat2f& a) {
    return {{a.x.x, a.y.x}, {a.x.y, a.y.y}};
}
/// matrix transpose
inline mat3f transpose(const mat3f& a) {
    return {
        {a.x.x, a.y.x, a.z.x}, {a.x.y, a.y.y, a.z.y}, {a.x.z, a.y.z, a.z.z}};
}
/// matrix transpose
inline mat4f transpose(const mat4f& a) {
    return {{a.x.x, a.y.x, a.z.x, a.w.x}, {a.x.y, a.y.y, a.z.y, a.w.y},
        {a.x.z, a.y.z, a.z.z, a.w.z}, {a.x.w, a.y.w, a.z.w, a.w.w}};
}

/// matrix adjugate
inline mat2f adjugate(const mat2f& a) {
    return {{a.y.y, -a.x.y}, {-a.y.x, a.x.x}};
}
/// matrix adjugate
inline mat3f adjugate(const mat3f& a) {
    return {{a.y.y * a.z.z - a.z.y * a.y.z, a.z.y * a.x.z - a.x.y * a.z.z,
                a.x.y * a.y.z - a.y.y * a.x.z},
        {a.y.z * a.z.x - a.z.z * a.y.x, a.z.z * a.x.x - a.x.z * a.z.x,
            a.x.z * a.y.x - a.y.z * a.x.x},
        {a.y.x * a.z.y - a.z.x * a.y.y, a.z.x * a.x.y - a.x.x * a.z.y,
            a.x.x * a.y.y - a.y.x * a.x.y}};
}
/// matrix adjugate
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

/// matrix determinant
inline float determinant(const mat2f& a) {
    return a.x.x * a.y.y - a.x.y * a.y.x;
}
/// matrix determinant
inline float determinant(const mat3f& a) {
    return a.x.x * (a.y.y * a.z.z - a.z.y * a.y.z) +
           a.x.y * (a.y.z * a.z.x - a.z.z * a.y.x) +
           a.x.z * (a.y.x * a.z.y - a.z.x * a.y.y);
}
/// matrix determinant
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

/// matrix inverse (uses adjugate and determinant)
inline mat2f inverse(const mat2f& a) { return adjugate(a) / determinant(a); }
/// matrix inverse (uses adjugate and determinant)
inline mat3f inverse(const mat3f& a) { return adjugate(a) / determinant(a); }
/// matrix inverse (uses adjugate and determinant)
inline mat4f inverse(const mat4f& a) { return adjugate(a) / determinant(a); }

/// stream write
inline ostream& operator<<(ostream& os, const mat2f& a) {
    return os << a.x << ' ' << a.y;
}
/// stream read
inline istream& operator>>(istream& is, mat2f& a) { return is >> a.x >> a.y; }
/// stream write
inline ostream& operator<<(ostream& os, const mat3f& a) {
    return os << a.x << ' ' << a.y << ' ' << a.z;
}
/// stream read
inline istream& operator>>(istream& is, mat3f& a) {
    return is >> a.x >> a.y >> a.z;
}
/// stream write
inline ostream& operator<<(ostream& os, const mat4f& a) {
    return os << a.x << ' ' << a.y << ' ' << a.z << ' ' << a.w;
}
/// stream read
inline istream& operator>>(istream& is, mat4f& a) {
    return is >> a.x >> a.y >> a.z >> a.w;
}
}  // namespace ygl

// -----------------------------------------------------------------------------
// RIGID BODY TRANSFORMS/FRAMES
// -----------------------------------------------------------------------------
namespace ygl {

/// Rigid transforms stored as a column-major affine matrix.
/// In memory, this representation is equivalent to storing an 3x3 rotation
/// followed by a 3x1 translation. Viewed this way, the representation allows
/// also to retrive the axis of the coordinate frame as the first 3 columns and
/// the translation as the 4th column. Colums access via operator[].
/// Access rotation and position with pos() and rot().
struct frame3f {
    /// size
    static const int N = 3;

    /// default constructor
    frame3f() : x{1, 0, 0}, y{0, 1, 0}, z{0, 0, 1}, o{0, 0, 0} {}

    /// element constructor
    frame3f(const vec3f& x, const vec3f& y, const vec3f& z, const vec3f& o)
        : x(x), y(y), z(z), o(o) {}

    /// element constructor
    frame3f(const mat3f& m, const vec3f& t) : x(m.x), y(m.y), z(m.z), o(t) {}

    /// conversion from matrix (assumes the matrix is a frame, so dangerous!)
    frame3f(const mat4f& m)
        : x(m.x.x, m.x.y, m.x.z)
        , y(m.y.x, m.y.y, m.y.z)
        , z(m.z.x, m.z.y, m.z.z)
        , o(m.w.x, m.w.y, m.w.z) {}

    /// element access
    vec3f& operator[](int i) { return (&x)[i]; }
    /// element access
    const vec3f& operator[](int i) const { return (&x)[i]; }

    /// data access
    vec3f* data() { return &x; }
    /// data access
    const vec3f* data() const { return &x; }

    /// access position
    vec3f& pos() { return o; }
    /// access position
    const vec3f& pos() const { return o; }

    /// access rotation
    mat3f& rot() { return *(mat3f*)(&x); }
    /// access rotation
    const mat3f& rot() const { return *(mat3f*)(&x); }

    /// element data
    vec3f x;
    /// element data
    vec3f y;
    /// element data
    vec3f z;
    /// element data
    vec3f o;
};

/// indentity frame
const auto identity_frame3f =
    frame3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};

// initializes a frame3 from origin and z.
inline frame3f make_frame_fromz(const vec3f& o, const vec3f& z_) {
    auto z = normalize(z_);
    auto x = normalize(orthogonal(z));
    auto y = normalize(cross(z, x));
    return {x, y, z, o};
}

// initializes a frame3 from origin, z and x.
inline frame3f make_frame3_fromzx(
    const vec3f& o, const vec3f& z_, const vec3f& x_) {
    auto z = normalize(z_);
    auto x = orthonormalize(x_, z);
    auto y = normalize(cross(z, x));
    return {x, y, z, o};
}

/// iteration support
inline vec3f* begin(frame3f& a) { return &a.x; }
/// iteration support
inline const vec3f* begin(const frame3f& a) { return &a.x; }
/// iteration support
inline vec3f* end(frame3f& a) { return &a.x + 4; }
/// iteration support
inline const vec3f* end(const frame3f& a) { return &a.x + 4; }

/// frame to matrix conversion
inline mat4f to_mat4f(const frame3f& a) {
    return {{a.x.x, a.x.y, a.x.z, 0}, {a.y.x, a.y.y, a.y.z, 0},
        {a.z.x, a.z.y, a.z.z, 0}, {a.o.x, a.o.y, a.o.z, 1}};
}

/// matrix to frame conversion
inline frame3f to_frame3f(const mat4f& a) {
    return {{a.x.x, a.x.y, a.x.z}, {a.y.x, a.y.y, a.y.z}, {a.z.x, a.z.y, a.z.z},
        {a.w.x, a.w.y, a.w.z}};
}

/// vector operator ==
inline bool operator==(const frame3f& a, const frame3f& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.o == b.o;
}
/// vector operator !=
inline bool operator!=(const frame3f& a, const frame3f& b) { return !(a == b); }

/// frame composition (equivalent to affine matrix multiply)
inline frame3f operator*(const frame3f& a, const frame3f& b) {
    return {a.rot() * b.rot(), a.rot() * b.pos() + a.pos()};
}

/// frame inverse (equivalent to rigid affine inverse)
inline frame3f inverse(const frame3f& a) {
    auto minv = transpose(a.rot());
    return {minv, -(minv * a.pos())};
}

/// stream write
inline ostream& operator<<(ostream& os, const frame3f& a) {
    return os << a.x << ' ' << a.y << ' ' << a.z << ' ' << a.o;
}

/// stream read
inline istream& operator>>(istream& is, frame3f& a) {
    return is >> a.x >> a.y >> a.z >> a.o;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// QUATERNIONS
// -----------------------------------------------------------------------------
namespace ygl {

/// Quaternions implemented as a vec<T,4>. Data access via operator[].
/// Quaterions are xi + yj + zk + w.
struct quat4f {
    /// default constructor
    quat4f() : x{0}, y{0}, z{0}, w{1} {}

    // list constructor
    quat4f(float x, float y, float z, float w) : x{x}, y{y}, z{z}, w{w} {}

    /// conversion from vec
    explicit quat4f(const vec4f& vv) : x{vv.x}, y{vv.y}, z{vv.z}, w{vv.w} {}
    /// conversion to vec
    explicit operator vec4f() const { return {x, y, z, w}; }

    /// element access
    float& operator[](int i) { return (&x)[i]; }
    /// element access
    const float& operator[](int i) const { return (&x)[i]; }

    /// data access
    float* data() { return &x; }
    /// data access
    const float* data() const { return &x; }

    /// data
    float x;
    /// data
    float y;
    /// data
    float z;
    /// data
    float w;
};

/// float identity quaterion
const auto identity_quat4f = quat4f{0, 0, 0, 1};

/// vector operator ==
inline bool operator==(const quat4f& a, const quat4f& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}

/// vector operator !=
inline bool operator!=(const quat4f& a, const quat4f& b) { return !(a == b); }

/// quaterion multiply
inline quat4f operator*(const quat4f& a, const quat4f& b) {
    return {a.x * b.w + a.w * b.x + a.y * b.w - a.z * b.y,
        a.y * b.w + a.w * b.y + a.z * b.x - a.x * b.z,
        a.z * b.w + a.w * b.z + a.x * b.y - a.y * b.x,
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z};
}

/// quaterion multiply
inline quat4f operator*(const quat4f& a, float b) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
}

/// quaterion division
inline quat4f operator/(const quat4f& a, float b) {
    return {a.x / b, a.y / b, a.z / b, a.w / b};
}

/// quaterion conjugate
inline quat4f conjugate(const quat4f& v) { return {-v.x, -v.y, -v.z, v.w}; }

/// quaterion inverse
inline quat4f inverse(const quat4f& v) {
    return conjugate(v) / dot(vec4f(v), vec4f(v));
}

/// quaterion inverse
inline quat4f normalize(const quat4f& v) {
    auto l = length(vec4f{v.x, v.y, v.z, v.w});
    if (!l) return {0, 0, 0, 1};
    return {v.x / l, v.y / l, v.z / l, v.w / l};
}

/// quaterion normalized linear interpolation
inline quat4f nlerp(const quat4f& a, const quat4f& b, float t) {
    return (quat4f)nlerp(
        vec4f(a), dot(vec4f(a), vec4f(b)) < 0 ? -vec4f(b) : vec4f(b), t);
}

/// quaterion spherical linear interpolation
inline quat4f slerp(const quat4f& a, const quat4f& b, float t) {
    return (quat4f)slerp(
        vec4f(a), dot(vec4f(a), vec4f(b)) < 0 ? -vec4f(b) : vec4f(b), t);
}

/// stream write
inline ostream& operator<<(ostream& os, const quat4f& a) {
    return os << a.x << ' ' << a.y << ' ' << a.z << ' ' << a.w;
}

/// stream read
inline istream& operator>>(istream& is, quat4f& a) {
    return is >> a.x >> a.y >> a.z >> a.w;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// AXIS ALIGNED BOUNDING BOXES
// -----------------------------------------------------------------------------
namespace ygl {

/// Axis aligned bounding box represented as a min/max vector pair.
struct bbox1f {
    /// initializes an invalid bbox
    bbox1f() : min{flt_max}, max{flt_min} {}
    /// list constructor
    bbox1f(float m, float M) : min{m}, max{M} {}

    /// element access
    float& operator[](int i) { return (&min)[i]; }
    /// element access
    const float& operator[](int i) const { return (&min)[i]; }

    /// element data
    float min;
    /// element data
    float max;
};

/// Axis aligned bounding box represented as a min/max vector pair.
struct bbox2f {
    /// initializes an invalid bbox
    bbox2f() : min{flt_max, flt_max}, max{flt_min, flt_min} {}
    /// list constructor
    bbox2f(const vec2f& m, const vec2f& M) : min{m}, max{M} {}

    /// element access
    vec2f& operator[](int i) { return (&min)[i]; }
    /// element access
    const vec2f& operator[](int i) const { return (&min)[i]; }

    /// element data
    vec2f min;
    /// element data
    vec2f max;
};

/// Axis aligned bounding box represented as a min/max vector pair.
struct bbox3f {
    /// initializes an invalid bbox
    bbox3f() : min{flt_max, flt_max, flt_max}, max{flt_min, flt_min, flt_min} {}
    /// list constructor
    bbox3f(const vec3f& m, const vec3f& M) : min{m}, max{M} {}

    /// element access
    vec3f& operator[](int i) { return (&min)[i]; }
    /// element access
    const vec3f& operator[](int i) const { return (&min)[i]; }

    /// element data
    vec3f min;
    /// element data
    vec3f max;
};

/// Axis aligned bounding box represented as a min/max vector pair.
struct bbox4f {
    /// initializes an invalid bbox
    bbox4f()
        : min{flt_max, flt_max, flt_max, flt_max}
        , max{flt_min, flt_min, flt_min, flt_min} {}
    /// list constructor
    bbox4f(const vec4f& m, const vec4f& M) : min{m}, max{M} {}

    /// element access
    vec4f& operator[](int i) { return (&min)[i]; }
    /// element access
    const vec4f& operator[](int i) const { return (&min)[i]; }

    /// element data
    vec4f min;
    /// element data
    vec4f max;
};

/// 1-dimensional float empty bbox
const auto invalid_bbox1f = bbox1f();
/// 2-dimensional float empty bbox
const auto invalid_bbox2f = bbox2f();
/// 3-dimensional float empty bbox
const auto invalid_bbox3f = bbox3f();
/// 4-dimensional float empty bbox
const auto invalid_bbox4f = bbox4f();

/// bbox operator ==
inline bool operator==(const bbox1f& a, const bbox1f& b) {
    return a.min == b.min && a.max == b.max;
}
/// bbox operator !=
inline bool operator!=(const bbox1f& a, const bbox1f& b) {
    return a.min != b.min || a.max != b.max;
}

/// bbox operator ==
inline bool operator==(const bbox2f& a, const bbox2f& b) {
    return a.min == b.min && a.max == b.max;
}
/// bbox operator !=
inline bool operator!=(const bbox2f& a, const bbox2f& b) {
    return a.min != b.min || a.max != b.max;
}

/// bbox operator ==
inline bool operator==(const bbox3f& a, const bbox3f& b) {
    return a.min == b.min && a.max == b.max;
}
/// bbox operator !=
inline bool operator!=(const bbox3f& a, const bbox3f& b) {
    return a.min != b.min || a.max != b.max;
}

/// bbox operator ==
inline bool operator==(const bbox4f& a, const bbox4f& b) {
    return a.min == b.min && a.max == b.max;
}
/// bbox operator !=
inline bool operator!=(const bbox4f& a, const bbox4f& b) {
    return a.min != b.min || a.max != b.max;
}

/// computes the center of a bbox
inline float bbox_center(const bbox1f& a) { return (a.min + a.max) / 2; }
/// computes the diagonal of a bbox
inline float bbox_diagonal(const bbox1f& a) { return a.max - a.min; }

/// computes the center of a bbox
inline vec2f bbox_center(const bbox2f& a) { return (a.min + a.max) / 2.0f; }
/// computes the diagonal of a bbox
inline vec2f bbox_diagonal(const bbox2f& a) { return a.max - a.min; }

/// computes the center of a bbox
inline vec3f bbox_center(const bbox3f& a) { return (a.min + a.max) / 2.0f; }
/// computes the diagonal of a bbox
inline vec3f bbox_diagonal(const bbox3f& a) { return a.max - a.min; }

/// computes the center of a bbox
inline vec4f bbox_center(const bbox4f& a) { return (a.min + a.max) / 2.0f; }
/// computes the diagonal of a bbox
inline vec4f bbox_diagonal(const bbox4f& a) { return a.max - a.min; }

/// expands a bounding box with a point
inline bbox1f expand(const bbox1f& a, float b) {
    return {min(a.min, b), max(a.max, b)};
}

/// expands a bounding box with a point
inline bbox2f expand(const bbox2f& a, const vec2f& b) {
    return {{min(a.min.x, b.x), min(a.min.y, b.y)},
        {max(a.max.x, b.x), max(a.max.y, b.y)}};
}

/// expands a bounding box with a point
inline bbox3f expand(const bbox3f& a, const vec3f& b) {
    return {{min(a.min.x, b.x), min(a.min.y, b.y), min(a.min.z, b.z)},
        {max(a.max.x, b.x), max(a.max.y, b.y), max(a.max.z, b.z)}};
}

/// expands a bounding box with a point
inline bbox4f expand(const bbox4f& a, const vec4f& b) {
    return {{min(a.min.x, b.x), min(a.min.y, b.y), min(a.min.z, b.z),
                min(a.min.w, b.w)},
        {max(a.max.x, b.x), max(a.max.y, b.y), max(a.max.z, b.z),
            max(a.max.w, b.w)}};
}

/// expands a bounding box with a bounding box
inline bbox1f expand(const bbox1f& a, const bbox1f& b) {
    return {min(a.min, b.min), max(a.max, b.max)};
}

/// expands a bounding box with a bounding box
inline bbox2f expand(const bbox2f& a, const bbox2f& b) {
    return {{min(a.min.x, b.min.x), min(a.min.y, b.min.y)},
        {max(a.max.x, b.max.x), max(a.max.y, b.max.y)}};
}

/// expands a bounding box with a bounding box
inline bbox3f expand(const bbox3f& a, const bbox3f& b) {
    return {
        {min(a.min.x, b.min.x), min(a.min.y, b.min.y), min(a.min.z, b.min.z)},
        {max(a.max.x, b.max.x), max(a.max.y, b.max.y), max(a.max.z, b.max.z)}};
}

/// expands a bounding box with a bounding box
inline bbox4f expand(const bbox4f& a, const bbox4f& b) {
    return {{min(a.min.x, b.min.x), min(a.min.y, b.min.y),
                min(a.min.z, b.min.z), min(a.min.w, b.min.w)},
        {max(a.max.x, b.max.x), max(a.max.y, b.max.y), max(a.max.z, b.max.z),
            max(a.max.w, b.max.w)}};
}

/// check if a bounding box contains a point
inline bool contains(const bbox3f& a, const vec3f& b) {
    if (a.min.x > b.x || a.max.x < b.x) return false;
    if (a.min.y > b.y || a.max.y < b.y) return false;
    if (a.min.z > b.z || a.max.z < b.z) return false;
    return true;
}

/// check if a bounding box contains a bounding box
inline bool contains(const bbox3f& a, const bbox3f& b) {
    if (a.min.x > b.max.x || a.max.x < b.min.x) return false;
    if (a.min.y > b.max.y || a.max.y < b.min.y) return false;
    if (a.min.z > b.max.z || a.max.z < b.min.z) return false;
    return true;
}

/// assign to expand()
inline bbox1f& operator+=(bbox1f& a, float b) { return a = expand(a, b); }
/// assign to expand()
inline bbox1f& operator+=(bbox1f& a, const bbox1f& b) {
    return a = expand(a, b);
}

/// assign to expand()
inline bbox2f& operator+=(bbox2f& a, const vec2f& b) {
    return a = expand(a, b);
}
/// assign to expand()
inline bbox2f& operator+=(bbox2f& a, const bbox2f& b) {
    return a = expand(a, b);
}

/// assign to expand()
inline bbox3f& operator+=(bbox3f& a, const vec3f& b) {
    return a = expand(a, b);
}
/// assign to expand()
inline bbox3f& operator+=(bbox3f& a, const bbox3f& b) {
    return a = expand(a, b);
}

/// assign to expand()
inline bbox4f& operator+=(bbox4f& a, const vec4f& b) {
    return a = expand(a, b);
}
/// assign to expand()
inline bbox4f& operator+=(bbox4f& a, const bbox4f& b) {
    return a = expand(a, b);
}

/// initialize a bonding box from a list of points
inline bbox1f make_bbox(int count, const float* v) {
    auto a = invalid_bbox1f;
    for (auto j = 0; j < count; j++) a += v[j];
    return a;
}

/// initialize a bonding box from a list of points
inline bbox2f make_bbox(int count, const vec2f* v) {
    auto a = invalid_bbox2f;
    for (auto j = 0; j < count; j++) a += v[j];
    return a;
}

/// initialize a bonding box from a list of points
inline bbox3f make_bbox(int count, const vec3f* v) {
    auto a = invalid_bbox3f;
    for (auto j = 0; j < count; j++) a += v[j];
    return a;
}

/// initialize a bonding box from a list of points
inline bbox4f make_bbox(int count, const vec4f* v) {
    auto a = invalid_bbox4f;
    for (auto j = 0; j < count; j++) a += v[j];
    return a;
}

/// initialize a bonding box from a list of points
inline bbox3f make_bbox(const initializer_list<vec3f>& v) {
    auto a = invalid_bbox3f;
    for (auto&& vv : v) a += vv;
    return a;
}

/// stream write
inline ostream& operator<<(ostream& os, const bbox1f& a) {
    return os << a.min << ' ' << a.max;
}
/// stream read
inline istream& operator>>(istream& is, bbox1f& a) {
    return is >> a.min >> a.max;
}

/// stream write
inline ostream& operator<<(ostream& os, const bbox2f& a) {
    return os << a.min << ' ' << a.max;
}
/// stream read
inline istream& operator>>(istream& is, bbox2f& a) {
    return is >> a.min >> a.max;
}

/// stream write
inline ostream& operator<<(ostream& os, const bbox3f& a) {
    return os << a.min << ' ' << a.max;
}
/// stream read
inline istream& operator>>(istream& is, bbox3f& a) {
    return is >> a.min >> a.max;
}

/// stream write
inline ostream& operator<<(ostream& os, const bbox4f& a) {
    return os << a.min << ' ' << a.max;
}
/// stream read
inline istream& operator>>(istream& is, bbox4f& a) {
    return is >> a.min >> a.max;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// PRIMITIVE BBOX FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

/// Point bounds
inline bbox3f point_bbox(const vec3f& p, float r = 0) {
    return bbox3f{p - vec3f{r, r, r}, p + vec3f{r, r, r}};
}

/// Line bounds
inline bbox3f line_bbox(
    const vec3f& v0, const vec3f& v1, float r0 = 0, float r1 = 0) {
    return make_bbox({v0 - vec3f{r0, r0, r0}, v0 + vec3f{r0, r0, r0},
        v1 - vec3f{r1, r1, r1}, v1 + vec3f{r1, r1, r1}});
}

/// Triangle bounds
inline bbox3f triangle_bbox(const vec3f& v0, const vec3f& v1, const vec3f& v2) {
    return make_bbox({v0, v1, v2});
}

/// Quad bounds
inline bbox3f quad_bbox(
    const vec3f& v0, const vec3f& v1, const vec3f& v2, const vec3f& v3) {
    return make_bbox({v0, v1, v2, v3});
}

/// Tetrahedron bounds
inline bbox3f tetrahedron_bbox(
    const vec3f& v0, const vec3f& v1, const vec3f& v2, const vec3f& v3) {
    return make_bbox({v0, v1, v2, v3});
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// RAYS
// -----------------------------------------------------------------------------
namespace ygl {

/// Rays with origin, direction and min/max t value.
struct ray3f {
    /// size
    static const int N = 3;
    /// type
    using T = float;

    /// origin
    vec3f o;
    /// direction
    vec3f d;
    /// minimum distance
    float tmin;
    /// maximum distance
    float tmax;

    /// default constructor
    ray3f() : o{0, 0, 0}, d{0, 0, 1}, tmin{0}, tmax{flt_max} {}
    /// initializes a ray from its elements
    ray3f(const vec3f& o, const vec3f& d, float tmin = 0, float tmax = flt_max)
        : o(o), d(d), tmin(tmin), tmax(tmax) {}
};

/// stream write
inline ostream& operator<<(ostream& os, const ray3f& a) {
    os << a.o << ' ' << a.d << ' ' << a.tmin << ' ' << a.tmax;
    return os;
}

/// stream read
inline istream& operator>>(istream& is, ray3f& a) {
    is >> a.o >> a.d >> a.tmin >> a.tmax;
    return is;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// TRANSFORMS
// -----------------------------------------------------------------------------
namespace ygl {

/// transforms a point by a matrix
inline vec3f transform_point(const mat4f& a, const vec3f& b) {
    auto vb = vec4f{b.x, b.y, b.z, 1};
    auto tvb = a * vb;
    return vec3f{tvb.x, tvb.y, tvb.z} / tvb.w;
}

/// transforms a vector by a matrix
inline vec3f transform_vector(const mat4f& a, const vec3f& b) {
    auto vb = vec4f{b.x, b.y, b.z, 0};
    auto tvb = a * vb;
    return vec3f{tvb.x, tvb.y, tvb.z};
}

/// transforms a direction by a matrix
inline vec3f transform_direction(const mat4f& a, const vec3f& b) {
    return normalize(transform_vector(a, b));
}

/// transforms a point by a frame (rigid affine transform)
inline vec3f transform_point(const frame3f& a, const vec3f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.o;
}

/// transforms a vector by a frame (rigid affine transform)
inline vec3f transform_vector(const frame3f& a, const vec3f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/// transforms a direction by a frame (rigid affine transform)
inline vec3f transform_direction(const frame3f& a, const vec3f& b) {
    return normalize(transform_vector(a, b));
}

/// transforms a frame by a frame (rigid affine transform)
inline frame3f transform_frame(const frame3f& a, const frame3f& b) {
    return {a.rot() * b.rot(), a.rot() * b.pos() + a.pos()};
}

/// inverse transforms a point by a frame (rigid affine transform)
inline vec3f transform_point_inverse(const frame3f& a, const vec3f& b) {
    return {dot(b - a.o, a.x), dot(b - a.o, a.y), dot(b - a.o, a.z)};
}

/// inverse transforms a vector by a frame (rigid affine transform)
inline vec3f transform_vector_inverse(const frame3f& a, const vec3f& b) {
    return {dot(b, a.x), dot(b, a.y), dot(b, a.z)};
}

/// inverse transforms a direction by a frame (rigid affine transform)
inline vec3f transform_direction_inverse(const frame3f& a, const vec3f& b) {
    return normalize(transform_vector_inverse(a, b));
}

/// transforms a ray by a matrix (direction is not normalized after)
inline ray3f transform_ray(const mat4f& a, const ray3f& b) {
    return {transform_point(a, b.o), transform_vector(a, b.d), b.tmin, b.tmax};
}

/// transforms a bbox by a matrix
inline bbox3f transform_bbox(const mat4f& a, const bbox3f& b) {
    vec3f corners[8] = {
        {b.min.x, b.min.y, b.min.z},
        {b.min.x, b.min.y, b.max.z},
        {b.min.x, b.max.y, b.min.z},
        {b.min.x, b.max.y, b.max.z},
        {b.max.x, b.min.y, b.min.z},
        {b.max.x, b.min.y, b.max.z},
        {b.max.x, b.max.y, b.min.z},
        {b.max.x, b.max.y, b.max.z},
    };
    auto xformed = bbox3f();
    for (auto j = 0; j < 8; j++) xformed += transform_point(a, corners[j]);
    return xformed;
}

/// transforms a ray by a frame (rigid affine transform)
inline ray3f transform_ray(const frame3f& a, const ray3f& b) {
    return {
        transform_point(a, b.o), transform_direction(a, b.d), b.tmin, b.tmax};
}

/// transforms a bbox by a frame (rigid affine transform)
inline bbox3f transform_bbox(const frame3f& a, const bbox3f& b) {
#if 0
    vec3f corners[8] = {
        {b.min.x, b.min.y, b.min.z}, {b.min.x, b.min.y, b.max.z},
        {b.min.x, b.max.y, b.min.z}, {b.min.x, b.max.y, b.max.z},
        {b.max.x, b.min.y, b.min.z}, {b.max.x, b.min.y, b.max.z},
        {b.max.x, b.max.y, b.min.z}, {b.max.x, b.max.y, b.max.z},
    };
    auto xformed = bbox<T, 3>();
    for (auto j = 0; j < 8; j++) xformed += transform_point(a, corners[j]);
    return xformed;
#else
    // Code from Real-time Collision Detection by Christer Ericson Sect. 4.2.6
    // Transform AABB a by the matrix m and translation t,
    // find maximum extents, and store result into AABB b.
    // start by adding in translation
    auto c = bbox3f{a.pos(), a.pos()};
    // for all three axes
    for (auto i = 0; i < 3; i++) {
        // form extent by summing smaller and larger terms respectively
        for (auto j = 0; j < 3; j++) {
            auto e = a.rot()[j][i] * b.min[j];
            auto f = a.rot()[j][i] * b.max[j];
            if (e < f) {
                c.min[i] += e;
                c.max[i] += f;
            } else {
                c.min[i] += f;
                c.max[i] += e;
            }
        }
    }
    return c;
#endif
}

/// inverse transforms a ray by a frame (rigid affine transform)
inline ray3f transform_ray_inverse(const frame3f& a, const ray3f& b) {
    return {transform_point_inverse(a, b.o),
        transform_direction_inverse(a, b.d), b.tmin, b.tmax};
}

/// inverse transforms a bbox by a frame (rigid affine transform)
inline bbox3f transform_bbox_inverse(const frame3f& a, const bbox3f& b) {
    return transform_bbox(inverse(a), b);
}

/// rotation matrix from axis-angle
inline mat3f rotation_mat3f(const vec3f& axis, float angle) {
    auto s = sin(angle), c = cos(angle);
    auto vv = normalize(axis);
    return {{c + (1 - c) * vv.x * vv.x, (1 - c) * vv.x * vv.y + s * vv.z,
                (1 - c) * vv.x * vv.z - s * vv.y},
        {(1 - c) * vv.x * vv.y - s * vv.z, c + (1 - c) * vv.y * vv.y,
            (1 - c) * vv.y * vv.z + s * vv.x},
        {(1 - c) * vv.x * vv.z + s * vv.y, (1 - c) * vv.y * vv.z - s * vv.x,
            c + (1 - c) * vv.z * vv.z}};
}

/// translation frame
inline frame3f translation_frame3f(const vec3f& a) {
    return {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, a};
}

/// translation matrix
inline mat4f translation_mat4f(const vec3f& a) {
    return to_mat4f(translation_frame3f(a));
}

/// scaling frame (this is not rigid and here for symmatry of API)
inline frame3f scaling_frame3f(const vec3f& a) {
    return {{a.x, 0, 0}, {0, a.y, 0}, {0, 0, a.z}, {0, 0, 0}};
}

/// scaling matrix
inline mat4f scaling_mat4f(const vec3f& a) {
    return to_mat4f(scaling_frame3f(a));
}

/// rotation frame
inline frame3f rotation_frame3f(const vec3f& axis, float angle) {
    return {rotation_mat3f(axis, angle), {0, 0, 0}};
}

/// rotation matrix
inline mat4f rotation_mat4f(const mat3f& rot) {
    return mat4f{{rot.x.x, rot.x.y, rot.x.z, 0}, {rot.y.x, rot.y.y, rot.y.z, 0},
        {rot.z.x, rot.z.y, rot.z.z, 0}, {0, 0, 0, 1}};
}

/// rotation matrix
inline mat4f rotation_mat4f(const vec3f& axis, float angle) {
    return rotation_mat4f(rotation_frame3f(axis, angle).rot());
}

/// quaternion axis-angle conversion
inline vec4f rotation_axisangle4(const quat4f& a) {
    auto axis = normalize(vec3f{a.x, a.y, a.z});
    auto angle = acos(a.w) * 2;
    return {axis.x, axis.y, axis.z, angle};
}

/// axis-angle to quaternion
inline quat4f rotation_quat4f(const vec4f& axis_angle) {
    auto axis = vec3f{axis_angle.x, axis_angle.y, axis_angle.z};
    auto len = length(axis);
    auto angle = atan2(len, axis_angle.w);
    if (len)
        axis /= len;
    else
        axis = {0, 0, 1};
    return {axis.x, axis.y, axis.z, angle};
}

/// quaterion to matrix conversion
inline mat3f rotation_mat3f(const quat4f& v) {
    return {{v.w * v.w + v.x * v.x - v.y * v.y - v.z * v.z,
                (v.x * v.y + v.z * v.w) * 2, (v.z * v.x - v.y * v.w) * 2},
        {(v.x * v.y - v.z * v.w) * 2,
            v.w * v.w - v.x * v.x + v.y * v.y - v.z * v.z,
            (v.y * v.z + v.x * v.w) * 2},
        {(v.z * v.x + v.y * v.w) * 2, (v.y * v.z - v.x * v.w) * 2,
            v.w * v.w - v.x * v.x - v.y * v.y + v.z * v.z}};
}

/// rotation matrix
inline mat4f rotation_mat4f(const quat4f& v) {
    return rotation_mat4f(rotation_mat3f(v));
}

/// matrix to quaternion
inline quat4f rotation_quat4f(const mat3f& m_) {
    // http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm
    auto q = quat4f();
    auto m = transpose(m_);
#if 1
    auto trace = m.x.x + m.y.y + m.z.z;
    if (trace > 0) {
        float s = 0.5f / sqrt(trace + 1);
        q.w = 0.25f / s;
        q.x = (m.z.y - m.y.z) * s;
        q.y = (m.x.z - m.z.x) * s;
        q.z = (m.y.x - m.x.y) * s;
    } else {
        if (m.x.x > m.y.y && m.x.x > m.z.z) {
            float s = 2 * sqrt(max(0.0f, 1 + m.x.x - m.y.y - m.z.z));
            q.w = (m.z.y - m.y.z) / s;
            q.x = 0.25f * s;
            q.y = (m.x.y + m.y.x) / s;
            q.z = (m.x.z + m.z.x) / s;
        } else if (m.y.y > m.z.z) {
            float s = 2 * sqrt(max(0.0f, 1 + m.y.y - m.x.x - m.z.z));
            q.w = (m.x.z - m.z.x) / s;
            q.x = (m.x.y + m.y.x) / s;
            q.y = 0.25f * s;
            q.z = (m.y.z + m.z.y) / s;
        } else {
            float s = 2 * sqrt(max(0.0f, 1 + m.z.z - m.x.x - m.y.y));
            q.w = (m.y.x - m.x.y) / s;
            q.x = (m.x.z + m.z.x) / s;
            q.y = (m.y.z + m.z.y) / s;
            q.z = 0.25f * s;
        }
    }

#else
    q.w = sqrt(max(0, 1 + m.x.x + m.y.y + m.z.z)) / 2;
    q.x = sqrt(max(0, 1 + m.x.x - m.y.y - m.z.z)) / 2;
    q.y = sqrt(max(0, 1 - m.x.x + m.y.y - m.z.z)) / 2;
    q.z = sqrt(max(0, 1 - m.x.x - m.y.y + m.z.z)) / 2;
    Q.x = copysign(q.x, m.z.y - m.y.z);
    Q.y = copysign(q.y, m.x.z - m.z.x);
    Q.z = copysign(q.z, m.y.x - m.x.y);
#endif

    return q;
}

/// OpenGL lookat frame
inline frame3f lookat_frame3f(const vec3f& eye, const vec3f& center,
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

/// OpenGL lookat matrix
inline mat4f lookat_mat4f(
    const vec3f& eye, const vec3f& center, const vec3f& up) {
    return to_mat4f(lookat_frame3f(eye, center, up));
}

/// OpenGL frustum matrix
inline mat4f frustum_mat4f(
    float l, float r, float b, float t, float n, float f) {
    return {{2 * n / (r - l), 0, 0, 0}, {0, 2 * n / (t - b), 0, 0},
        {(r + l) / (r - l), (t + b) / (t - b), -(f + n) / (f - n), -1},
        {0, 0, -2 * f * n / (f - n), 0}};
}

/// OpenGL orthographic matrix
inline mat4f ortho_mat4f(float l, float r, float b, float t, float n, float f) {
    return {{2 / (r - l), 0, 0, 0}, {0, 2 / (t - b), 0, 0},
        {0, 0, -2 / (f - n), 0},
        {-(r + l) / (r - l), -(t + b) / (t - b), -(f + n) / (f - n), 1}};
}

/// OpenGL orthographic 2D matrix
inline mat4f ortho2d_mat4f(float left, float right, float bottom, float top) {
    return ortho_mat4f(left, right, bottom, top, -1, 1);
}

/// OpenGL/GLTF orthographic matrix
inline mat4f ortho_mat4f(float xmag, float ymag, float near, float far) {
    return {{1 / xmag, 0, 0, 0}, {0, 1 / ymag, 0, 0},
        {0, 0, 2 / (near - far), 0}, {0, 0, (far + near) / (near - far), 1}};
}

/// OpenGL/GLTF perspective matrix
inline mat4f perspective_mat4f(
    float fovy, float aspect, float near, float far) {
    auto tg = tan(fovy / 2);
    return {{1 / (aspect * tg), 0, 0, 0}, {0, 1 / tg, 0, 0},
        {0, 0, (far + near) / (near - far), -1},
        {0, 0, 2 * far * near / (near - far), 0}};
}

/// OpenGL/GLTF infinite perspective matrix
inline mat4f perspective_mat4f(float fovy, float aspect, float near) {
    auto tg = tan(fovy / 2);
    return {{1 / (aspect * tg), 0, 0, 0}, {0, 1 / tg, 0, 0}, {0, 0, -1, -1},
        {0, 0, 2 * near, 0}};
}

/// Decompose an affine matrix into translation, rotation, scale.
/// Assumes there is no shear and the matrix is affine.
inline void decompose_mat4f(
    const mat4f& m, vec3f& translation, mat3f& rotation, vec3f& scale) {
    translation = {m.w.x, m.w.y, m.w.z};
    rotation.x = {m.x.x, m.x.y, m.x.z};
    rotation.y = {m.y.x, m.y.y, m.y.z};
    rotation.z = {m.z.x, m.z.y, m.z.z};
    scale = {length(rotation.x), length(rotation.y), length(rotation.z)};
    rotation = {
        normalize(rotation.x), normalize(rotation.y), normalize(rotation.z)};
}

/// Convert a rotation matrix to a quaternion
inline quat4f to_quat4f(const mat3f& a) {
    auto q = quat4f();
    // from
    // http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/
    float trace = a[0][0] + a[1][1] + a[2][2];
    if (trace > 0) {  // I changed M_EPSILON to 0
        float s = 0.5f / sqrtf(trace + 1.0f);
        q.w = 0.25f / s;
        q.x = (a[1][2] - a[2][1]) * s;
        q.y = (a[2][0] - a[0][2]) * s;
        q.z = (a[0][1] - a[1][0]) * s;
    } else {
        if (a[0][0] > a[1][1] && a[0][0] > a[2][2]) {
            float s = 2.0f * sqrt(1.0f + a[0][0] - a[1][1] - a[2][2]);
            q.w = (a[1][2] - a[2][1]) / s;
            q.x = 0.25f * s;
            q.y = (a[1][0] + a[0][1]) / s;
            q.z = (a[2][0] + a[0][2]) / s;
        } else if (a[1][1] > a[2][2]) {
            float s = 2.0f * sqrtf(1.0f + a[1][1] - a[0][0] - a[2][2]);
            q.w = (a[2][0] - a[0][2]) / s;
            q.x = (a[1][0] + a[0][1]) / s;
            q.y = 0.25f * s;
            q.z = (a[2][1] + a[1][2]) / s;
        } else {
            float s = 2.0f * sqrtf(1.0f + a[2][2] - a[0][0] - a[1][1]);
            q.w = (a[0][1] - a[1][0]) / s;
            q.x = (a[2][0] + a[0][2]) / s;
            q.y = (a[2][1] + a[1][2]) / s;
            q.z = 0.25f * s;
        }
    }
    return q;
}

/// Decompose an affine matrix into translation, rotation, scale.
/// Assumes there is no shear and the matrix is affine.
inline void decompose_mat4f(
    const mat4f& m, vec3f& translation, quat4f& rotation, vec3f& scale) {
    auto rot_matrix = mat3f();
    decompose_mat4f(m, translation, rot_matrix, scale);
    rotation = to_quat4f(rot_matrix);
}

/// Decompose an affine matrix into translation, rotation, scale.
/// Assumes there is no shear and the matrix is affine.
inline mat4f compose_mat4f(
    const vec3f& translation, const mat3f& rotation, const vec3f& scale) {
    return translation_mat4f(translation) * scaling_mat4f(scale) *
           rotation_mat4f(rotation);
}

/// Decompose an affine matrix into translation, rotation, scale.
/// Assumes there is no shear and the matrix is affine.
inline mat4f compose_mat4f(
    const vec3f& translation, const quat4f& rotation, const vec3f& scale) {
    return translation_mat4f(translation) * scaling_mat4f(scale) *
           rotation_mat4f(rotation);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// UI UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

/// Turntable for UI navigation from a from/to/up parametrization of the
/// camera.
inline void camera_turntable(vec3f& from, vec3f& to, vec3f& up,
    const vec3f& rotate, float dolly, const vec3f& pan) {
    // rotate if necessary
    if (rotate.x || rotate.y) {
        auto z = normalize(to - from);
        auto lz = length(to - from);
        auto phi = atan2(z.z, z.x) + rotate.x;
        auto theta = acos(z.y) + rotate.y;
        theta = clamp(theta, 0.001f, pif - 0.001f);
        auto nz = vec3f{sin(theta) * cos(phi) * lz, cos(theta) * lz,
            sin(theta) * sin(phi) * lz};
        from = to - nz;
    }

    // dolly if necessary
    if (dolly) {
        auto z = normalize(to - from);
        auto lz = max(0.001f, length(to - from) * (1 + dolly));
        z *= lz;
        from = to - z;
    }

    // pan if necessary
    if (pan.x || pan.y) {
        auto z = normalize(to - from);
        auto x = normalize(cross(up, z));
        auto y = normalize(cross(z, x));
        auto t = vec3f{pan.x * x.x + pan.y * y.x, pan.x * x.y + pan.y * y.y,
            pan.x * x.z + pan.y * y.z};
        from += t;
        to += t;
    }
}

/// Turntable for UI navigation for a frame/distance parametrization of the
/// camera.
inline void camera_turntable(frame3f& frame, float& focus, const vec2f& rotate,
    float dolly, const vec2f& pan) {
    // rotate if necessary
    if (rotate.x || rotate.y) {
        auto phi = atan2(frame.z.z, frame.z.x) + rotate.x;
        auto theta = acos(frame.z.y) + rotate.y;
        theta = clamp(theta, 0.001f, pif - 0.001f);
        auto new_z =
            vec3f{sin(theta) * cos(phi), cos(theta), sin(theta) * sin(phi)};
        auto new_center = frame.o - frame.z * focus;
        auto new_o = new_center + new_z * focus;
        frame = lookat_frame3f(new_o, new_center, {0, 1, 0});
        focus = length(new_o - new_center);
    }

    // pan if necessary
    if (dolly) {
        auto c = frame.o - frame.z * focus;
        focus = max(focus + dolly, 0.001f);
        frame.o = c + frame.z * focus;
    }

    // pan if necessary
    if (pan.x || pan.y) { frame.o += frame.x * pan.x + frame.y * pan.y; }
}

/// FPS camera for UI navigation for a frame parametrization.
/// https://gamedev.stackexchange.com/questions/30644/how-to-keep-my-quaternion-using-fps-camera-from-tilting-and-messing-up
inline void camera_fps(
    frame3f& frame, const vec3f& transl, const vec2f& rotate) {
    auto y = vec3f{0, 1, 0};
    auto z = orthonormalize(frame.z, y);
    auto x = cross(y, z);

    frame.rot() = rotation_mat3f({1, 0, 0}, rotate.y) * frame.rot() *
                  rotation_mat3f({0, 1, 0}, rotate.x);
    frame.pos() += transl.x * x + transl.y * y + transl.z * z;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// RANDOM NUMBER GENERATION
// -----------------------------------------------------------------------------
namespace ygl {

/// PCG random numbers. A family of random number generators that supports
/// multiple sequences. In our code, we allocate one sequence for each sample.
/// PCG32 from http://www.pcg-random.org/
struct rng_pcg32 {
    /// RNG state.
    uint64_t state = 0x853c49e6748fea9bULL;
    /// RNG sequence. Must be odd.
    uint64_t inc = 0xda3e39cb94b95bdbULL;
};

/// Next random number
inline uint32_t advance_rng(rng_pcg32& rng) {
    uint64_t oldstate = rng.state;
    rng.state = oldstate * 6364136223846793005ULL + rng.inc;
    uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
    uint32_t rot = (uint32_t)(oldstate >> 59u);
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
    // return (xorshifted >> rot) | (xorshifted << ((~rot + 1u) & 31));
}

/// Multi-step advance function (jump-ahead, jump-back).
inline void advance_rng(rng_pcg32& rng, uint64_t delta) {
    // The method used here is based on Brown, "Random Number Generation
    // with Arbitrary Stride", Transactions of the American Nuclear
    // Society (Nov. 1994). The algorithm is very similar to fast
    // exponentiation.
    uint64_t cur_mult = 6364136223846793005ULL, cur_plus = rng.inc,
             acc_mult = 1u, acc_plus = 0u;

    while (delta > 0) {
        if (delta & 1) {
            acc_mult *= cur_mult;
            acc_plus = acc_plus * cur_mult + cur_plus;
        }
        cur_plus = (cur_mult + 1) * cur_plus;
        cur_mult *= cur_mult;
        delta /= 2;
    }
    rng.state = acc_mult * rng.state + acc_plus;
}

/// Multi-step advance function (jump-ahead, jump-back).
inline void advance_rng(rng_pcg32& rng, int64_t delta) {
    // Even though delta is an unsigned integer, we can pass a signed
    // integer to go backwards, it just goes "the long way round".
    advance_rng(rng, (uint64_t)delta);
}

/// Seeds a random number generator with a state state from the sequence seq.
inline void seed_rng(rng_pcg32& rng, uint64_t state, uint64_t seq = 1) {
    rng.state = 0U;
    rng.inc = (seq << 1u) | 1u;
    advance_rng(rng);
    rng.state += state;
    advance_rng(rng);
}

/// Init a random number generator with a state state from the sequence seq.
inline rng_pcg32 init_rng(uint64_t state, uint64_t seq = 1) {
    auto rng = rng_pcg32();
    seed_rng(rng, state, seq);
    return rng;
}

/// Next random uint in [0,n) range with proper weighting
inline uint32_t next_rand1i(rng_pcg32& rng, uint32_t n) {
#if YGL_RNG_FASTUINT
    return advance_rng(rng) % n;
#else
    // To avoid bias, we need to make the range of the RNG a multiple of
    // bound, which we do by dropping output less than a threshold.
    // A naive scheme to calculate the threshold would be to do
    //
    //     uint32_t threshold = 0x100000000ull % bound;
    //
    // but 64-bit div/mod is slower than 32-bit div/mod (especially on
    // 32-bit platforms).  In essence, we do
    //
    //     uint32_t threshold = (0x100000000ull-bound) % bound;
    //
    // because this version will calculate the same modulus, but the LHS
    // value is less than 2^32.
    uint32_t threshold = (~n + 1u) % n;

    // Uniformity guarantees that this loop will terminate.  In practice, it
    // should usually terminate quickly; on average (assuming all bounds are
    // equally likely), 82.25% of the time, we can expect it to require just
    // one iteration.  In the worst case, someone passes a bound of 2^31 + 1
    // (i.e., 2147483649), which invalidates almost 50% of the range.  In
    // practice, bounds are typically small and only a tiny amount of the
    // range is eliminated.
    while (true) {
        uint32_t r = advance_rng(rng);
        if (r >= threshold) return r % n;
    }
#endif
}

/// Next random float in [0,1).
inline float next_rand1f(rng_pcg32& rng) {
#if 1
    // Trick from MTGP: generate an uniformly distributed
    // single precision number in [1,2) and subtract 1.
    union {
        uint32_t u;
        float f;
    } x;
    x.u = (advance_rng(rng) >> 9) | 0x3f800000u;
    return x.f - 1.0f;
#else
    const static auto scale = (float)(1.0 / numeric_limits<uint32_t>::max());
    return advance_rng(rng) * scale;
#endif
}

/// Next random float in [a,b).
inline float next_rand1f(rng_pcg32& rng, float a, float b) {
    return a + (b - a) * next_rand1f(rng);
}

/// Next random float2 in [0,1)x[0,1).
inline vec2f next_rand2f(rng_pcg32& rng) {
    return {next_rand1f(rng), next_rand1f(rng)};
}

/// Next random float in [a.x,b.x)x[a.y,b.y).
inline vec2f next_rand2f(rng_pcg32& rng, const vec2f& a, const vec2f& b) {
    return {next_rand1f(rng, a.x, b.x), next_rand1f(rng, a.y, b.y)};
}

/// Next random float3 in [0,1)x[0,1)x[0,1).
inline vec3f next_rand3f(rng_pcg32& rng) {
    return {next_rand1f(rng), next_rand1f(rng), next_rand1f(rng)};
}

/// Next random float in [a.x,b.x)x[a.y,b.y)x[a.z,b.z).
inline vec3f next_rand2f(rng_pcg32& rng, const vec3f& a, const vec3f& b) {
    return {next_rand1f(rng, a.x, b.x), next_rand1f(rng, a.y, b.y),
        next_rand1f(rng, a.z, b.z)};
}

/// Next random double in [0, 1). Only 32 mantissa bits are filled, but still
/// better than float that uses 23.
inline double next_rand1d(rng_pcg32& rng) {
#if 1
    // Trick from MTGP: generate an uniformly distributed
    // double precision number in [1,2) and subtract 1.
    union {
        uint64_t u;
        double d;
    } x;
    x.u = ((uint64_t)advance_rng(rng) << 20) | 0x3ff0000000000000ull;
    return x.d - 1.0;
#else
    const static auto scale = (double)(1.0 / numeric_limits<uint32_t>::max());
    return advance_rng(rng) * scale;
#endif
}

/// Distance between random number generators
inline int64_t rng_distance(const rng_pcg32& a, const rng_pcg32& b) {
    assert(a.inc == b.inc);

    uint64_t cur_mult = 6364136223846793005ULL, cur_plus = a.inc,
             cur_state = b.state, the_bit = 1u, distance = 0u;

    while (a.state != cur_state) {
        if ((a.state & the_bit) != (cur_state & the_bit)) {
            cur_state = cur_state * cur_mult + cur_plus;
            distance |= the_bit;
        }
        assert((a.state & the_bit) == (cur_state & the_bit));
        the_bit <<= 1;
        cur_plus = (cur_mult + 1ULL) * cur_plus;
        cur_mult *= cur_mult;
    }

    return (int64_t)distance;
}

/// Random shuffle of a sequence.
template <typename Iterator>
inline void rng_shuffle(rng_pcg32& rng, Iterator begin, Iterator end) {
    // Draw uniformly distributed permutation and permute the
    // given STL container. From Knuth, TAoCP Vol.2(3rd 3d),
    // Section 3.4.2
    for (Iterator it = end - 1; it > begin; --it)
        std::iter_swap(
            it, begin + next_rand1i(rng, (uint32_t)(it - begin + 1)));
}

/// Random shuffle of a sequence.
template <typename T>
inline void rng_shuffle(rng_pcg32& rng, T* vals, int num) {
    // Draw uniformly distributed permutation and permute the
    // given STL container
    for (auto i = num - 1; i > 0; --i)
        swap(vals[i], vals[next_rand1i(rng, (uint32_t)(i - 1))]);
}

/// Random shuffle of a sequence.
template <typename T>
inline void rng_shuffle(rng_pcg32& rng, vector<T>& vals) {
    shuffle(rng, vals.data(), vals.size());
}

/// Equality operator
inline bool operator==(const rng_pcg32& a, const rng_pcg32& b) {
    return a.state == b.state && a.inc == b.inc;
}

/// Inequality operator
inline bool operator!=(const rng_pcg32& a, const rng_pcg32& b) {
    return a.state != b.state || a.inc != b.inc;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// MONETACARLO SAMPLING FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

/// sample hemispherical direction with uniform distribution
inline vec3f sample_hemisphere(const vec2f& ruv) {
    auto z = ruv.y;
    auto r = sqrt(1 - z * z);
    auto phi = 2 * pif * ruv.x;
    return vec3f(r * cos(phi), r * sin(phi), z);
}

/// pdf for hemispherical direction with uniform distribution
inline float sample_hemisphere_pdf(const vec3f& w) {
    return (w.z <= 0) ? 0 : 1 / (2 * pif);
}

/// spherical direction with uniform distribution
inline vec3f sample_sphere(const vec2f ruv) {
    auto z = 2 * ruv.y - 1;
    auto r = sqrt(1 - z * z);
    auto phi = 2 * pif * ruv.x;
    return vec3f(r * cos(phi), r * sin(phi), z);
}

/// pdf for spherical direction with uniform distribution
inline float sample_sphere_pdf(const vec3f& w) { return 1 / (4 * pif); }

/// hemispherical direction with cosine distribution
inline vec3f sample_hemisphere_cosine(const vec2f& ruv) {
    auto z = sqrt(ruv.y);
    auto r = sqrt(1 - z * z);
    auto phi = 2 * pif * ruv.x;
    return vec3f(r * cos(phi), r * sin(phi), z);
}

/// pdf for hemispherical direction with cosine distribution
inline float sample_hemisphere_cosine_pdf(const vec3f& w) {
    return (w.z <= 0) ? 0 : w.z / pif;
}

/// hemispherical direction with cosine power distribution
inline vec3f sample_hemisphere_cospower(const vec2f& ruv, float n) {
    auto z = pow(ruv.y, 1 / (n + 1));
    auto r = sqrt(1 - z * z);
    auto phi = 2 * pif * ruv.x;
    return vec3f(r * cos(phi), r * sin(phi), z);
}

/// pdf for hemispherical direction with cosine power distribution
inline float sample_hemisphere_cospower_pdf(const vec3f& w, float n) {
    return (w.z <= 0) ? 0 : pow(w.z, n) * (n + 1) / (2 * pif);
}

/// uniform disk
inline vec3f sample_disk(const vec2f& ruv) {
    auto r = sqrt(ruv.y);
    auto phi = 2 * pif * ruv.x;
    return vec3f(cos(phi) * r, sin(phi) * r, 0);
}

/// pdf for uniform disk
inline float sample_disk_pdf() { return 1 / pif; }

/// uniform cylinder
inline vec3f sample_cylinder(const vec2f& ruv) {
    auto phi = 2 * pif * ruv.x;
    return vec3f(sin(phi), cos(phi), ruv.y * 2 - 1);
}

/// pdf for uniform cylinder
inline float sample_cylinder_pdf() { return 1 / pif; }

/// uniform triangle
inline vec2f sample_triangle(const vec2f& ruv) {
    return {1 - sqrt(ruv.x), ruv.y * sqrt(ruv.x)};
}

/// uniform triangle
inline vec3f sample_triangle(
    const vec2f& ruv, const vec3f& v0, const vec3f& v1, const vec3f& v2) {
    auto uv = sample_triangle(ruv);
    return v0 * (1 - uv.x - uv.y) + v1 * uv.x + v2 * uv.y;
}

/// pdf for uniform triangle (triangle area)
inline float sample_triangle_pdf(
    const vec3f& v0, const vec3f& v1, const vec3f& v2) {
    return 2 / length(cross(v1 - v0, v2 - v0));
}

/// index with uniform distribution
inline int sample_index(float r, int size) {
    return clamp((int)(r * size), 0, size - 1);
}

/// pdf for index with uniform distribution
inline float sample_index_pdf(int size) { return 1.0f / size; }

}  // namespace ygl

// -----------------------------------------------------------------------------
// HASHING
// -----------------------------------------------------------------------------
namespace ygl {

/// Computes the i-th term of a permutation of l values keyed by p.
/// From Correlated Multi-Jittered Sampling by Kensler @ Pixar
inline uint32_t hash_permute(uint32_t i, uint32_t n, uint32_t key) {
    uint32_t w = n - 1;
    w |= w >> 1;
    w |= w >> 2;
    w |= w >> 4;
    w |= w >> 8;
    w |= w >> 16;
    do {
        i ^= key;
        i *= 0xe170893du;
        i ^= key >> 16;
        i ^= (i & w) >> 4;
        i ^= key >> 8;
        i *= 0x0929eb3f;
        i ^= key >> 23;
        i ^= (i & w) >> 1;
        i *= 1 | key >> 27;
        i *= 0x6935fa69;
        i ^= (i & w) >> 11;
        i *= 0x74dcb303;
        i ^= (i & w) >> 2;
        i *= 0x9e501cc3;
        i ^= (i & w) >> 2;
        i *= 0xc860a3df;
        i &= w;
        i ^= i >> 5;
    } while (i >= n);
    return (i + key) % n;
}

/// Computes a float value by hashing i with a key p.
/// From Correlated Multi-Jittered Sampling by Kensler @ Pixar
inline float hash_randfloat(uint32_t i, uint32_t key) {
    i ^= key;
    i ^= i >> 17;
    i ^= i >> 10;
    i *= 0xb36534e5;
    i ^= i >> 12;
    i ^= i >> 21;
    i *= 0x93fc4795;
    i ^= 0xdf6e307f;
    i ^= i >> 17;
    i *= 1 | key >> 18;
    return i * (1.0f / 4294967808.0f);
}

/// 32 bit integer hash. Public domain code.
inline uint32_t hash_uint32(uint64_t a) {
    a -= (a << 6);
    a ^= (a >> 17);
    a -= (a << 9);
    a ^= (a << 4);
    a -= (a << 3);
    a ^= (a << 10);
    a ^= (a >> 15);
    return a;
}

/// 64 bit integer hash. Public domain code.
inline uint64_t hash_uint64(uint64_t a) {
    a = (~a) + (a << 21);  // a = (a << 21) - a - 1;
    a ^= (a >> 24);
    a += (a << 3) + (a << 8);  // a * 265
    a ^= (a >> 14);
    a += (a << 2) + (a << 4);  // a * 21
    a ^= (a >> 28);
    a += (a << 31);
    return a;
}

/// 64-to-32 bit integer hash. Public domain code.
inline uint32_t hash_uint64_32(uint64_t a) {
    a = (~a) + (a << 18);  // a = (a << 18) - a - 1;
    a ^= (a >> 31);
    a *= 21;  // a = (a + (a << 2)) + (a << 4);
    a ^= (a >> 11);
    a += (a << 6);
    a ^= (a >> 22);
    return (uint32_t)a;
}

/// Combines two 64 bit hashes as in boost::hash_combine
inline size_t hash_combine(size_t a, size_t b) {
    return a ^ (b + 0x9e3779b9 + (a << 6) + (a >> 2));
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// PERLIN NOISE FUNCTION
// -----------------------------------------------------------------------------
namespace ygl {
// noise functions from stb_perlin.h

/// Compute the revised Pelin noise function. Wrap provides a wrapping noise
/// but must be power of two (wraps at 256 anyway). For octave based noise,
/// good values are obtained with octaves=6 (numerber of noise calls),
/// lacunarity=~2.0 (spacing between successive octaves: 2.0 for warpping
/// output), gain=0.5 (relative weighting applied to each successive octave),
/// offset=1.0 (used to invert the ridges).
float perlin_noise(const vec3f& p, const vec3i& wrap = zero3i);
/// Ridge noise function
float perlin_ridge_noise(const vec3f& p, float lacunarity = 2.0f,
    float gain = 0.5f, float offset = 1.0f, int octaves = 6,
    const vec3i& wrap = zero3i);
/// Fractal brownian motion noise - see perlin_noise() for params.
float perlin_fbm_noise(const vec3f& p, float lacunarity = 2.0f,
    float gain = 0.5f, int octaves = 6, const vec3i& wrap = zero3i);
/// Fractal turbulence noise - see perlin_noise() for params.
float perlin_turbulence_noise(const vec3f& p, float lacunarity = 2.0f,
    float gain = 0.5f, int octaves = 6, const vec3i& wrap = zero3i);

}  // namespace ygl

// -----------------------------------------------------------------------------
// PYTHON-LIKE ITERATORS AND CONTAINER OPERATIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Implementation of Python-like range generator. Create it with the the
// `range()` functions to use argument deduction.
// clang-format off
struct range_generator {
    struct iterator {
        iterator(int pos, int step) : _pos(pos), _step(step) {}
        bool operator!=(const iterator& a) const { return _pos < a._pos; }
        iterator& operator++() { _pos += _step; return *this; }
        int operator*() const { return _pos; }
        private: int _pos, _step;
    };
    range_generator(int min, int max, int step)
        : _min(min), _max(max), _step(step) {}
    iterator begin() const { return {_min, _step}; }
    iterator end() const { return {_max, _step}; }
    private: int _min, _max, _step;
};
// clang-format on

/// Python-like range
inline range_generator range(int max) { return {0, max, 1}; }

/// Python-like range
inline range_generator range(int min, int max, int step = 1) {
    return {min, max, step};
}

// Implemenetation of Python-like enumerate. Create it with the function
// `enumerate()` to use argument deduction.
// clang-format off
template <typename T>
struct enumerate_generator {
    struct iterator {
        iterator(int pos, T* data) : _pos(pos), _data(data) {}
        bool operator!=(const iterator& a) const { return _pos < a._pos; }
        iterator& operator++() { _pos += 1; _data += 1; return *this; }
        pair<int, T&> operator*() const { return {_pos, *_data}; }
        private: int _pos; T* _data;
    };
    enumerate_generator(int max, T* data) : _max(max), _data(data) {}
    iterator begin() const { return {0, _data}; }
    iterator end() const { return {_max, _data + _max}; }
    private: int _max; T* _data;
};
// clang-format on

/// Python-like range
template <typename T>
inline enumerate_generator<const T> enumerate(const vector<T>& vv) {
    return {(int)vv.size(), vv.data()};
}

/// Python-like range
template <typename T>
inline enumerate_generator<T> enumerate(vector<T>& vv) {
    return {(int)vv.size(), vv.data()};
}

/// Append an element to a vector
template <typename T>
inline vector<T> operator+(const vector<T>& v, const T& vv) {
    auto vc = vector<T>();
    vc.reserve(v.size() + 1);
    vc.insert(vc.end(), v.begin(), v.end());
    vc.push_back(vv);
    return vc;
}

/// Append an element to a vector
template <typename T>
inline vector<T>& operator+=(vector<T>& v, const T& vv) {
    v.push_back(vv);
    return v;
}

/// Append an element to a vector
template <typename T, typename ET>
inline vector<T> operator+(const vector<T>& v, const ET& vv) {
    auto vc = vector<T>();
    vc.reserve(v.size() + 1);
    vc.insert(vc.end(), v.begin(), v.end());
    vc.push_back(vv);
    return vc;
}

/// Append an element to a vector
template <typename T, typename ET>
inline vector<T>& operator+=(vector<T>& v, const ET& vv) {
    v.push_back(vv);
    return v;
}

/// Append a vector to a vector
template <typename T>
inline vector<T> operator+(const vector<T>& v, const vector<T>& vv) {
    auto vc = vector<T>();
    vc.reserve(v.size() + vv.size());
    vc.insert(vc.end(), v.begin(), v.end());
    vc.insert(vc.end(), vv.begin(), vv.end());
    return vc;
}

/// Append a vector to a vector
template <typename T>
inline vector<T>& operator+=(vector<T>& v, const vector<T>& vv) {
    v.insert(v.end(), vv.begin(), vv.end());
    return v;
}

/// Get a key
template <typename K, typename V>
inline K get_key(const std::vector<pair<K, V>>& kvs, const V& v) {
    for (auto& kv : kvs)
        if (kv.second == v) return kv.first;
    throw runtime_error("key not found");
}

/// Get a value
template <typename K, typename V>
inline V get_value(const std::vector<pair<K, V>>& kvs, const K& k) {
    for (auto& kv : kvs)
        if (kv.first == k) return kv.second;
    throw runtime_error("key not found");
}

/// Find a value in an array
template <typename T>
inline int find_idx(const vector<T>& v, const T& vv) {
    auto pos = find(v.begin(), v.end(), vv);
    if (pos == v.end()) return -1;
    return pos = v.begin();
}

/// Checks if a containers contains a value
template <typename T>
inline bool contains(const vector<T>& v, const T& vv) {
    return find(v.begin(), v.end(), vv) != v.end();
}

/// Checks if a containers contains a value
template <typename K, typename V>
inline bool contains(const map<K, V>& v, const K& vv) {
    return v.find(vv) != v.end();
}

/// Checks if a containers contains a value
template <typename K, typename V>
inline bool contains(const unordered_map<K, V>& v, const K& vv) {
    return v.find(vv) != v.end();
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// GEOMETRY UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

/// line tangent
inline vec3f line_tangent(const vec3f& v0, const vec3f& v1) {
    return normalize(v1 - v0);
}

/// line length
inline float line_length(const vec3f& v0, const vec3f& v1) {
    return length(v1 - v0);
}

/// triangle normal
inline vec3f triangle_normal(
    const vec3f& v0, const vec3f& v1, const vec3f& v2) {
    return normalize(cross(v1 - v0, v2 - v0));
}

/// triangle area
inline float triangle_area(const vec3f& v0, const vec3f& v1, const vec3f& v2) {
    return length(cross(v1 - v0, v2 - v0)) / 2;
}

/// quad area
inline float quad_area(
    const vec3f& v0, const vec3f& v1, const vec3f& v2, const vec3f& v3) {
    return triangle_area(v0, v1, v3) + triangle_area(v3, v2, v1);
}

/// tetrahedron volume
inline float tetrahedron_volume(
    const vec3f& v0, const vec3f& v1, const vec3f& v2, const vec3f& v3) {
    return dot(cross(v1 - v0, v2 - v0), v3 - v0) / 6;
}

// Triangle tangent and bitangent from uv (not othornormalized with themselfves
// not the normal). Follows the definition in
// http://www.terathon.com/code/tangent.html and
// https://gist.github.com/aras-p/2843984
inline pair<vec3f, vec3f> triangle_tangents_fromuv(const vec3f& v0,
    const vec3f& v1, const vec3f& v2, const vec2f& uv0, const vec2f& uv1,
    const vec2f& uv2) {
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

/// line barycentric interpolation
template <typename T>
inline T eval_barycentric_point(const vector<T>& vals, const int& p, float w) {
    if (vals.empty()) return T();
    return vals[p] * w;
}

/// line barycentric interpolation
template <typename T>
inline T eval_barycentric_line(
    const vector<T>& vals, const vec2i& l, const vec2f& w) {
    if (vals.empty()) return T();
    return vals[l.x] * w.x + vals[l.y] * w.y;
}

/// triangle barycentric interpolation
template <typename T>
inline T eval_barycentric_triangle(
    const vector<T>& vals, const vec3i& t, const vec3f& w) {
    if (vals.empty()) return T();
    return vals[t.x] * w.x + vals[t.y] * w.y + vals[t.z] * w.z;
}

/// tetrahedron barycentric interpolation
template <typename T>
inline T eval_barycentric_tetra(
    const vector<T>& vals, const vec4i& t, const vec4f& w) {
    if (vals.empty()) return T();
    return vals[t.x] * w.x + vals[t.y] * w.y + vals[t.z] * w.z +
           vals[t.w] * w.w;
}

/// quad interpolation based on the two-triangle representation
template <typename T>
inline T eval_barycentric_quad(
    const vector<T>& vals, const vec4i& t, const vec4f& w) {
    if (vals.empty()) return T();
    return vals[t.x] * w.x + vals[t.y] * w.y + vals[t.z] * w.z +
           vals[t.w] * w.w;
}

/// bernstein polynomials (for Bezier)
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

/// bernstein polynomials (for Bezier)
template <typename T>
inline T eval_bernstein_derivative(T u, int i, int degree) {
    return degree * (eval_bernstein(u, i - 1, degree - 1) -
                        eval_bernstein(u, i, degree - 1));
}

/// eval bezier
template <typename T, typename T1>
inline T eval_bezier_cubic(
    const T& v0, const T& v1, const T& v2, const T& v3, T1 t) {
    return v0 * eval_bernstein(t, 0, 3) + v1 * eval_bernstein(t, 1, 3) +
           v2 * eval_bernstein(t, 2, 3) + v3 * eval_bernstein(t, 3, 3);
}

/// eval bezier
template <typename T, typename T1>
inline T eval_bezier_cubic(const vector<T>& vals, const vec4i& b, T1 t) {
    if (vals.empty()) return T();
    return eval_bezier_cubic(vals[b.x], vals[b.y], vals[b.z], vals[b.w], t);
}

/// eval bezier derivative
template <typename T, typename T1>
inline T eval_bezier_cubic_derivative(
    const T& v0, const T& v1, const T& v2, const T& v3, T1 t) {
    return v0 * eval_bernstein_derivative(t, 0, 3) +
           v1 * eval_bernstein_derivative(t, 1, 3) +
           v2 * eval_bernstein_derivative(t, 2, 3) +
           v3 * eval_bernstein_derivative(t, 3, 3);
}

/// eval bezier derivative
template <typename T, typename T1>
inline T eval_bezier_cubic_derivative(
    const vector<T>& vals, const vec4i& b, T1 t) {
    if (vals.empty()) return T();
    return eval_bezier_cubic_derivative(
        vals[b.x], vals[b.y], vals[b.z], vals[b.w], t);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// SHAPE UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

/// Compute per-vertex normals/tangents for lines, triangles and quads with
/// positions pos. Weighted indicated whether the normals/tangents are
/// weighted by line length.
inline vector<vec3f> compute_normals(const vector<vec2i>& lines,
    const vector<vec3i>& triangles, const vector<vec4i>& quads,
    const vector<vec3f>& pos, bool weighted = true) {
    auto norm = vector<vec3f>(pos.size(), zero3f);
    for (auto& l : lines) {
        auto n = pos[l.y] - pos[l.x];
        if (!weighted) n = normalize(n);
        for (auto vid : l) norm[vid] += n;
    }
    for (auto& t : triangles) {
        auto n = cross(pos[t.y] - pos[t.x], pos[t.z] - pos[t.x]);
        if (!weighted) n = normalize(n);
        for (auto vid : t) norm[vid] += n;
    }
    for (auto& q : quads) {
        auto n = cross(pos[q.y] - pos[q.x], pos[q.w] - pos[q.x]) +
                 cross(pos[q.w] - pos[q.z], pos[q.x] - pos[q.z]);
        if (!weighted) n = normalize(n);
        for (auto vid : q) norm[vid] += n;
    }
    for (auto& n : norm) n = normalize(n);
    return norm;
}

/// Compute per-vertex tangent frame for triangle meshes.
/// Tangent space is defined by a four component vector.
/// The first three components are the tangent with respect to the U texcoord.
/// The fourth component is the sign of the tangent wrt the V texcoord.
/// Tangent frame is useful in normal mapping.
inline vector<vec4f> compute_tangent_frames(const vector<vec3i>& triangles,
    const vector<vec3f>& pos, const vector<vec3f>& norm,
    const vector<vec2f>& texcoord, bool weighted = true) {
    auto tangu = vector<vec3f>(pos.size(), zero3f);
    auto tangv = vector<vec3f>(pos.size(), zero3f);
    for (auto& t : triangles) {
        auto tutv = triangle_tangents_fromuv(pos[t.x], pos[t.y], pos[t.z],
            texcoord[t.x], texcoord[t.y], texcoord[t.z]);
        if (!weighted) tutv = {normalize(tutv.first), normalize(tutv.second)};
        for (auto vid : t) tangu[vid] += tutv.first;
        for (auto vid : t) tangv[vid] += tutv.second;
    }
    for (auto& t : tangu) t = normalize(t);
    for (auto& t : tangv) t = normalize(t);
    auto tangsp = vector<vec4f>(pos.size(), zero4f);
    for (auto i = 0; i < pos.size(); i++) {
        tangu[i] = orthonormalize(tangu[i], norm[i]);
        auto s = (dot(cross(norm[i], tangu[i]), tangv[i]) < 0) ? -1.0f : 1.0f;
        tangsp[i] = {tangu[i].x, tangu[i].y, tangu[i].z, s};
    }
    return tangsp;
}

/// Apply skinning
inline void compute_skinning(const vector<vec3f>& pos,
    const vector<vec3f>& norm, const vector<vec4f>& weights,
    const vector<vec4i>& joints, const vector<mat4f>& xforms,
    vector<vec3f>& skinned_pos, vector<vec3f>& skinned_norm) {
    skinned_pos.resize(pos.size());
    skinned_norm.resize(norm.size());
    for (auto i = 0; i < pos.size(); i++) {
        skinned_pos[i] =
            transform_point(xforms[joints[i].x], pos[i]) * weights[i].x +
            transform_point(xforms[joints[i].y], pos[i]) * weights[i].y +
            transform_point(xforms[joints[i].z], pos[i]) * weights[i].z +
            transform_point(xforms[joints[i].w], pos[i]) * weights[i].w;
    }
    for (auto i = 0; i < pos.size(); i++) {
        skinned_norm[i] = normalize(
            transform_direction(xforms[joints[i].x], norm[i]) * weights[i].x +
            transform_direction(xforms[joints[i].y], norm[i]) * weights[i].y +
            transform_direction(xforms[joints[i].z], norm[i]) * weights[i].z +
            transform_direction(xforms[joints[i].w], norm[i]) * weights[i].w);
    }
}

/// Apply skinning
inline void compute_skinning(const vector<vec3f>& pos,
    const vector<vec3f>& norm, const vector<vec4f>& weights,
    const vector<vec4i>& joints, const vector<frame3f>& xforms,
    vector<vec3f>& skinned_pos, vector<vec3f>& skinned_norm) {
    skinned_pos.resize(pos.size());
    skinned_norm.resize(norm.size());
    for (auto i = 0; i < pos.size(); i++) {
        skinned_pos[i] =
            transform_point(xforms[joints[i].x], pos[i]) * weights[i].x +
            transform_point(xforms[joints[i].y], pos[i]) * weights[i].y +
            transform_point(xforms[joints[i].z], pos[i]) * weights[i].z +
            transform_point(xforms[joints[i].w], pos[i]) * weights[i].w;
    }
    for (auto i = 0; i < pos.size(); i++) {
        skinned_norm[i] = normalize(
            transform_direction(xforms[joints[i].x], norm[i]) * weights[i].x +
            transform_direction(xforms[joints[i].y], norm[i]) * weights[i].y +
            transform_direction(xforms[joints[i].z], norm[i]) * weights[i].z +
            transform_direction(xforms[joints[i].w], norm[i]) * weights[i].w);
    }
}

/// Apply skinning as specified in Khronos glTF
inline void compute_matrix_skinning(const vector<vec3f>& pos,
    const vector<vec3f>& norm, const vector<vec4f>& weights,
    const vector<vec4i>& joints, const vector<mat4f>& xforms,
    vector<vec3f>& skinned_pos, vector<vec3f>& skinned_norm) {
    skinned_pos.resize(pos.size());
    skinned_norm.resize(norm.size());
    for (auto i = 0; i < pos.size(); i++) {
        auto xform = xforms[joints[i].x] * weights[i].x +
                     xforms[joints[i].y] * weights[i].y +
                     xforms[joints[i].z] * weights[i].z +
                     xforms[joints[i].w] * weights[i].w;
        skinned_pos[i] = transform_point(xform, pos[i]);
        skinned_norm[i] = normalize(transform_direction(xform, norm[i]));
    }
}

/// Create an array of edges.
inline vector<vec2i> get_edges(const vector<vec2i>& lines,
    const vector<vec3i>& triangles, const vector<vec4i>& quads) {
    auto edges = vector<vec2i>();
    auto eset = unordered_set<vec2i>();
    for (auto e : lines) {
        e = {min(e.x, e.y), max(e.x, e.y)};
        if (!eset.insert(e).second) continue;
        eset.insert({e.y, e.x});
        edges += e;
    }
    for (auto& t : triangles) {
        for (auto e : {vec2i{t.x, t.y}, vec2i{t.y, t.z}, vec2i{t.z, t.x}}) {
            e = {min(e.x, e.y), max(e.x, e.y)};
            if (!eset.insert(e).second) continue;
            eset.insert({e.y, e.x});
            edges += e;
        }
    }
    for (auto& q : quads) {
        for (auto e : {vec2i{q.x, q.y}, vec2i{q.y, q.z}, vec2i{q.z, q.w},
                 vec2i{q.w, q.x}}) {
            if (e.x == e.y) continue;
            e = {min(e.x, e.y), max(e.x, e.y)};
            if (!eset.insert(e).second) continue;
            eset.insert({e.y, e.x});
            edges += e;
        }
    }

    return edges;
}

/// Create an array of boundary edges. Lines are always considered boundaries.
inline vector<vec2i> get_boundary_edges(const vector<vec2i>& lines,
    const vector<vec3i>& triangles, const vector<vec4i>& quads) {
    auto ecount = unordered_map<vec2i, int>();

    // lines are added manually later
    for (auto l : lines) { ecount.insert({l, 2}); }
    for (auto t : triangles) {
        for (auto e : {vec2i{t.x, t.y}, vec2i{t.y, t.z}, vec2i{t.z, t.x}}) {
            e = {min(e.x, e.y), max(e.x, e.y)};
            auto ins = ecount.insert({e, 1});
            if (!ins.second) ins.first->second += 1;
        }
    }
    for (auto q : quads) {
        for (auto e : {vec2i{q.x, q.y}, vec2i{q.y, q.z}, vec2i{q.z, q.w},
                 vec2i{q.w, q.x}}) {
            if (e.x == e.y) continue;
            e = {min(e.x, e.y), max(e.x, e.y)};
            auto ins = ecount.insert({e, 1});
            if (!ins.second) ins.first->second += 1;
        }
    }

    auto boundary = lines;
    for (auto ec : ecount) {
        if (ec.second > 1) continue;
        boundary += ec.first;
    }

    return boundary;
}

/// Get a list of all unique vertices.
inline vector<int> get_verts(const vector<vec2i>& lines,
    const vector<vec3i>& triangles, const vector<vec4i>& quads) {
    auto verts = vector<int>();
    auto vset = unordered_set<int>();
    for (auto l : lines)
        for (auto vid : l)
            if (vset.insert(vid).second) verts += vid;
    for (auto t : triangles)
        for (auto vid : t)
            if (vset.insert(vid).second) verts += vid;
    for (auto q : quads)
        for (auto vid : q)
            if (vset.insert(vid).second) verts += vid;
    return verts;
}

/// Create an array of boundary vertices. Lines are always considered
/// boundaries.
inline vector<int> get_boundary_verts(const vector<vec2i>& lines,
    const vector<vec3i>& triangles, const vector<vec4i>& quads) {
    return get_verts(get_boundary_edges(lines, triangles, quads), {}, {});
}

/// Convert quads to triangles
inline vector<vec3i> convert_quads_to_triangles(const vector<vec4i>& quads) {
    auto triangles = vector<vec3i>();
    triangles.reserve(quads.size() * 2);
    for (auto& q : quads) {
        triangles += {q.x, q.y, q.w};
        if (q.z != q.w) triangles += {q.z, q.w, q.y};
    }
    return triangles;
}

/// Convert quads to triangles with a diamond-like topology.
/// Quads have to be consecutive one row after another.
inline vector<vec3i> convert_quads_to_triangles(
    const vector<vec4i>& quads, int row_length) {
    auto triangles = vector<vec3i>();
    triangles.reserve(quads.size() * 2);
    for (auto& q : quads) {
        triangles += {q.x, q.y, q.w};
        if (q.z != q.w) triangles += {q.z, q.w, q.y};
    }
    return triangles;
#if 0
        triangles.resize(usteps * vsteps * 2);
        for (auto j = 0; j < vsteps; j++) {
            for (auto i = 0; i < usteps; i++) {
                auto& f1 = triangles[(j * usteps + i) * 2 + 0];
                auto& f2 = triangles[(j * usteps + i) * 2 + 1];
                if ((i + j) % 2) {
                    f1 = {vid(i, j), vid(i + 1, j), vid(i + 1, j + 1)};
                    f2 = {vid(i + 1, j + 1), vid(i, j + 1), vid(i, j)};
                } else {
                    f1 = {vid(i, j), vid(i + 1, j), vid(i, j + 1)};
                    f2 = {vid(i + 1, j + 1), vid(i, j + 1), vid(i + 1, j)};
                }
            }
        }
#endif
    return triangles;
}

/// Convert face varying data to single primitives. Returns the quads indices
/// and filled vectors for pos, norm and texcoord.
inline tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>
convert_face_varying(const vector<vec4i>& quads_pos,
    const vector<vec4i>& quads_norm, const vector<vec4i>& quads_texcoord,
    const vector<vec3f>& pos, const vector<vec3f>& norm,
    const vector<vec2f>& texcoord) {
    // make faces unique
    unordered_map<vec3i, int> vert_map;
    auto quads = vector<vec4i>(quads_pos.size());
    for (auto fid = 0; fid < quads_pos.size(); fid++) {
        for (auto c = 0; c < 4; c++) {
            auto v = vec3i{
                quads_pos[fid][c],
                (!quads_norm.empty()) ? quads_norm[fid][c] : -1,
                (!quads_texcoord.empty()) ? quads_texcoord[fid][c] : -1,
            };
            if (vert_map.find(v) == vert_map.end()) {
                auto s = (int)vert_map.size();
                vert_map[v] = s;
            }
            quads[fid][c] = vert_map.at(v);
        }
    }

    // fill vert data
    auto qpos = vector<vec3f>();
    if (!pos.empty()) {
        qpos.resize(vert_map.size());
        for (auto& kv : vert_map) { qpos[kv.second] = pos[kv.first.x]; }
    }
    auto qnorm = vector<vec3f>();
    if (!norm.empty()) {
        qnorm.resize(vert_map.size());
        for (auto& kv : vert_map) { qnorm[kv.second] = norm[kv.first.y]; }
    }
    auto qtexcoord = vector<vec2f>();
    if (!texcoord.empty()) {
        qtexcoord.resize(vert_map.size());
        for (auto& kv : vert_map) {
            qtexcoord[kv.second] = texcoord[kv.first.z];
        }
    }

    // done
    return {quads, qpos, qnorm, qtexcoord};
}

// wrapper for implementation below
inline float _subdivide_normalize(float x) { return x; }
inline vec2f _subdivide_normalize(const vec2f& x) { return normalize(x); }
inline vec3f _subdivide_normalize(const vec3f& x) { return normalize(x); }
inline vec4f _subdivide_normalize(const vec4f& x) { return normalize(x); }

/// Tesselate lines, triangles and quads by spolitting edges.
/// Returns the tesselated elements and dictionaries for vertex calculations.
inline tuple<vector<vec2i>, vector<vec3i>, vector<vec4i>, vector<vec2i>,
    vector<vec4i>>
subdivide_elems(const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads, int nverts) {
    if (!nverts) return {};
    auto emap = unordered_map<vec2i, int>();
    auto edges = vector<vec2i>();
    for (auto e : lines) {
        if (contains(emap, e)) continue;
        emap[{e.x, e.y}] = nverts + (int)edges.size();
        emap[{e.y, e.x}] = nverts + (int)edges.size();
        edges += e;
    }
    for (auto& t : triangles) {
        for (auto e : {vec2i{t.x, t.y}, vec2i{t.y, t.z}, vec2i{t.z, t.x}}) {
            if (contains(emap, e)) continue;
            emap[{e.x, e.y}] = nverts + (int)edges.size();
            emap[{e.y, e.x}] = nverts + (int)edges.size();
            edges += e;
        }
    }
    for (auto& q : quads) {
        for (auto e : {vec2i{q.x, q.y}, vec2i{q.y, q.z}, vec2i{q.z, q.w},
                 vec2i{q.w, q.x}}) {
            if (e.x == e.y) continue;
            if (contains(emap, e)) continue;
            emap[{e.x, e.y}] = nverts + (int)edges.size();
            emap[{e.y, e.x}] = nverts + (int)edges.size();
            edges += e;
        }
    }

    auto tlines = vector<vec2i>();
    tlines.reserve(lines.size() * 2);
    for (auto& l : lines) {
        tlines += {l.x, emap.at(l)};
        tlines += {emap.at(l), l.y};
    }

    auto ttriangles = vector<vec3i>();
    ttriangles.reserve(triangles.size() * 4);
    for (auto& t : triangles) {
        ttriangles.push_back({t.x, emap.at({t.x, t.y}), emap.at({t.z, t.x})});
        ttriangles.push_back({t.y, emap.at({t.y, t.z}), emap.at({t.x, t.y})});
        ttriangles.push_back({t.z, emap.at({t.z, t.x}), emap.at({t.y, t.z})});
        ttriangles.push_back(
            {emap.at({t.x, t.y}), emap.at({t.y, t.z}), emap.at({t.z, t.x})});
    }

    auto tquads = vector<vec4i>();
    tquads.reserve(quads.size() * 4);
    for (auto fkv : enumerate(quads)) {
        auto f = fkv.second;
        auto fvert = nverts + (int)edges.size() + fkv.first;
        if (f.z != f.w) {
            tquads += {f.x, emap.at({f.x, f.y}), fvert, emap.at({f.w, f.x})};
            tquads += {f.y, emap.at({f.y, f.z}), fvert, emap.at({f.x, f.y})};
            tquads += {f.z, emap.at({f.z, f.w}), fvert, emap.at({f.y, f.z})};
            tquads += {f.w, emap.at({f.w, f.x}), fvert, emap.at({f.z, f.w})};
        } else {
            tquads += {f.x, emap.at({f.x, f.y}), fvert, emap.at({f.z, f.x})};
            tquads += {f.y, emap.at({f.y, f.z}), fvert, emap.at({f.x, f.y})};
            tquads += {f.z, emap.at({f.z, f.x}), fvert, emap.at({f.y, f.z})};
        }
    }
    tquads.shrink_to_fit();

    return {tlines, ttriangles, tquads, edges, quads};
}

/// Subdivide vertex properties given the maps
template <typename T>
inline vector<T> subdivide_vert(const vector<T>& vert,
    const vector<vec2i>& edges, const vector<vec4i>& faces,
    bool normalized = false) {
    if (vert.empty()) return {};

    auto tvert = vector<T>();
    tvert.reserve(vert.size() + edges.size() + faces.size());

    tvert += vert;
    for (auto e : edges) tvert += (vert[e.x] + vert[e.y]) / 2;
    for (auto f : faces) {
        if (f.z != f.w)
            tvert += (vert[f.x] + vert[f.y] + vert[f.z] + vert[f.w]) / 4;
        else
            tvert += (vert[f.x] + vert[f.y] + vert[f.z]) / 3;
    }

    if (normalized) {
        for (auto& n : tvert) n = _subdivide_normalize(n);
    }

    return tvert;
}

/// Performs the smoothing step of Catmull-Clark. Start with a tesselate quad
/// mesh obtained with subdivide_elems() and subdivide_vert(). To handle open
/// meshes with boundary, get the boundary from make_boundary_edge() and pass it
/// as crease_lines. To fix the boundary entirely, just get the boundary
/// vertices and pass it as creases.
template <typename T>
inline vector<T> subdivide_catmullclark(const vector<vec4i>& quads,
    const vector<T>& vert, const vector<vec2i>& crease_tlines,
    const vector<int>& crease_tpoints, bool normalized = false) {
    if (quads.empty() || vert.empty()) return vert;

    // define vertex valence ---------------------------
    auto val = vector<int>(vert.size(), 2);
    for (auto e : crease_tlines)
        for (auto vid : e) val[vid] = 1;
    for (auto vid : crease_tpoints) val[vid] = 0;

    // averaging pass ----------------------------------
    auto tvert = vector<T>(vert.size(), T());
    auto count = vector<int>(vert.size(), 0);
    for (auto p : crease_tpoints) {
        auto c = vert[p];
        if (val[p] == 0) tvert[p] += c;
        if (val[p] == 0) count[p] += 1;
    }
    for (auto e : crease_tlines) {
        auto c = (vert[e.x] + vert[e.y]) / 2.0f;
        for (auto vid : e) {
            if (val[vid] == 1) tvert[vid] += c;
            if (val[vid] == 1) count[vid] += 1;
        }
    }
    for (auto& f : quads) {
        auto c = (vert[f.x] + vert[f.y] + vert[f.z] + vert[f.w]) / 4.0f;
        for (auto vid : f) {
            if (val[vid] == 2) tvert[vid] += c;
            if (val[vid] == 2) count[vid] += 1;
        }
    }
    for (auto i = 0; i < vert.size(); i++) { tvert[i] /= (float)count[i]; }

    // correction pass ----------------------------------
    // p = p + (avg_p - p) * (4/avg_count)
    for (auto i = 0; i < vert.size(); i++) {
        if (val[i] != 2) continue;
        tvert[i] = vert[i] + (tvert[i] - vert[i]) * (4.0f / count[i]);
    }

    if (normalized) {
        for (auto& v : tvert) v = _subdivide_normalize(v);
    }

    return tvert;
}

/// Generate a rectangular grid of usteps x vsteps uv values for parametric
/// surface generation.
inline tuple<vector<vec4i>, vector<vec2f>> make_uvquads(int usteps, int vsteps,
    bool uwrap = false, bool vwrap = false, bool vpole0 = false,
    bool vpole1 = false) {
    auto uvert = (uwrap) ? usteps : usteps + 1;
    auto vvert = (vwrap) ? vsteps : vsteps + 1;
    auto vid = [=](int i, int j) {
        if (uwrap) i = i % usteps;
        if (vwrap) j = j % vsteps;
        return j * uvert + i;
    };

    auto uv = vector<vec2f>(uvert * vvert);
    for (auto j = 0; j < vvert; j++) {
        for (auto i = 0; i < uvert; i++) {
            uv[vid(i, j)] = {i / (float)usteps, j / (float)vsteps};
        }
    }

    auto quads = vector<vec4i>(usteps * vsteps);
    for (auto j = 0; j < vsteps; j++) {
        for (auto i = 0; i < usteps; i++) {
            quads[j * usteps + i] = {
                vid(i, j), vid(i + 1, j), vid(i + 1, j + 1), vid(i, j + 1)};
        }
    }

    if (vpole0) {
        if (vwrap) throw runtime_error("cannot have a pole with wrapping");
        uv = vector<vec2f>(uv.begin() + uvert, uv.end());
        uv.insert(uv.begin(), {0, 0});
        for (auto& q : quads) {
            for (auto& vid : q) { vid = (vid < usteps) ? 0 : vid - uvert + 1; }
            if (q.x == 0 && q.y == 0) q = {q.z, q.w, q.x, q.y};
        }
    }

    if (vpole1) {
        if (vwrap) throw runtime_error("cannot have a pole with wrapping");
        auto pid = (int)uv.size() - uvert;
        uv = vector<vec2f>(uv.begin(), uv.end() - uvert);
        uv.insert(uv.end(), {0, 1});
        for (auto& q : quads) {
            for (auto& vid : q) { vid = (vid < pid) ? vid : pid; }
        }
    }

    return {quads, uv};
}

/// Generate parametric num lines of usteps segments.
inline tuple<vector<vec2i>, vector<vec2f>> make_uvlines(int num, int usteps) {
    auto vid = [usteps](int i, int j) { return j * (usteps + 1) + i; };
    auto uv = vector<vec2f>((usteps + 1) * num);
    for (auto j = 0; j < num; j++) {
        for (auto i = 0; i <= usteps; i++) {
            uv[vid(i, j)] = {i / (float)usteps, j / (float)num};
        }
    }

    auto lines = vector<vec2i>(usteps * num);
    for (int j = 0; j < num; j++) {
        for (int i = 0; i < usteps; i++) {
            lines[j * usteps + i] = {vid(i, j), vid(i + 1, j)};
        }
    }

    return {lines, uv};
}

/// Generate a parametric point set. Mostly here for completeness.
inline tuple<vector<int>, vector<vec2f>> make_uvpoints(int num) {
    auto uv = vector<vec2f>(num);
    for (auto i = 0; i < num; i++) { uv[i] = {i / (float)num, 0}; }

    auto points = vector<int>(num);
    for (auto i = 0; i < num; i++) points[i] = i;

    return {points, uv};
}

/// Merge elements between shapes. The elements are merged by increasing the
/// array size of the second array by the number of vertices of the first.
/// Vertex data can then be concatenated successfully.
inline tuple<vector<vec2i>, vector<vec3i>, vector<vec4i>> merge_elems(
    int nverts, const vector<vec2i>& lines1, const vector<vec3i>& triangles1,
    const vector<vec4i>& quads1, const vector<vec2i>& lines2,
    const vector<vec3i>& triangles2, const vector<vec4i>& quads2) {
    auto lines = lines1 + lines2;
    auto triangles = triangles1 + triangles2;
    auto quads = quads1 + quads2;
    for (auto i = lines1.size(); i < lines.size(); i++)
        lines[i] += {nverts, nverts};
    for (auto i = triangles1.size(); i < triangles.size(); i++)
        triangles[i] += {nverts, nverts, nverts};
    for (auto i = quads1.size(); i < quads.size(); i++)
        quads[i] += {nverts, nverts, nverts, nverts};
    return {lines, triangles, quads};
}

/// Unshare shape data by duplicating all vertex data for each element,
/// giving a faceted look. Note that faceted tangents are not computed.
inline tuple<vector<vec2i>, vector<vec3i>, vector<vec4i>, vector<int>>
facet_elems(const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads) {
    auto verts = vector<int>();
    auto nlines = vector<vec2i>();
    for (auto l : lines) {
        nlines.push_back({(int)verts.size(), (int)verts.size() + 1});
        for (auto v : l) verts += v;
    }

    auto ntriangles = vector<vec3i>();
    for (auto t : triangles) {
        ntriangles.push_back(
            {(int)verts.size(), (int)verts.size() + 1, (int)verts.size() + 2});
        for (auto v : t) verts += v;
    }

    auto nquads = vector<vec4i>();
    for (auto q : quads) {
        if (q.z != q.w) {
            nquads.push_back({(int)verts.size(), (int)verts.size() + 1,
                (int)verts.size() + 2, (int)verts.size() + 3});
            for (auto v : q) verts += v;
        } else {
            nquads.push_back({(int)verts.size(), (int)verts.size() + 1,
                (int)verts.size() + 2, (int)verts.size() + 2});
            for (auto v : q.xyz()) verts += v;
        }
    }

    return {nlines, ntriangles, nquads, verts};
}

/// Unshare vertices for faceting
template <typename T>
inline vector<T> facet_vert(const vector<T>& vert, const vector<int>& vmap) {
    if (vert.empty()) return vert;
    auto tvert = vector<T>(vmap.size());
    for (auto vkv : enumerate(vmap)) tvert[vkv.first] = vert[vkv.second];
    return tvert;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// SHAPE SAMPLING
// -----------------------------------------------------------------------------
namespace ygl {

/// Pick a point
inline int sample_points(int npoints, float re) {
    return clamp(0, npoints - 1, (int)(re * npoints));
}

/// Compute a distribution for sampling points uniformly
inline vector<float> sample_points_cdf(int npoints) {
    auto cdf = vector<float>(npoints);
    for (auto i = 0; i < npoints; i++) cdf[i] = i + 1;
    return cdf;
}

/// Pick a point
inline int sample_points(const vector<float>& cdf, float re) {
    re = clamp(re * cdf.back(), 0.0f, cdf.back() - 0.00001f);
    return (int)(std::upper_bound(cdf.begin(), cdf.end(), re) - cdf.begin());
}

/// Compute a distribution for sampling lines uniformly
inline vector<float> sample_lines_cdf(
    const vector<vec2i>& lines, const vector<vec3f>& pos) {
    auto cdf = vector<float>(lines.size());
    for (auto i = 0; i < lines.size(); i++)
        cdf[i] = length(pos[lines[i].x] - pos[lines[i].y]);
    for (auto i = 1; i < lines.size(); i++) cdf[i] += cdf[i - 1];
    return cdf;
}

/// Pick a point on lines
inline pair<int, vec2f> sample_lines(
    const vector<float>& cdf, float re, float ruv) {
    re = clamp(re * cdf.back(), 0.0f, cdf.back() - 0.00001f);
    auto eid =
        (int)(std::upper_bound(cdf.begin(), cdf.end(), re) - cdf.begin());
    return {eid, {1 - ruv, ruv}};
}

/// Compute a distribution for sampling triangle meshes uniformly
inline vector<float> sample_triangles_cdf(
    const vector<vec3i>& triangles, const vector<vec3f>& pos) {
    auto cdf = vector<float>(triangles.size());
    for (auto i = 0; i < triangles.size(); i++)
        cdf[i] = triangle_area(
            pos[triangles[i].x], pos[triangles[i].y], pos[triangles[i].z]);
    for (auto i = 1; i < triangles.size(); i++) cdf[i] += cdf[i - 1];
    return cdf;
}

/// Pick a point on a triangle mesh
inline pair<int, vec3f> sample_triangles(
    const vector<float>& cdf, float re, const vec2f& ruv) {
    re = clamp(re * cdf.back(), 0.0f, cdf.back() - 0.00001f);
    auto eid =
        (int)(std::upper_bound(cdf.begin(), cdf.end(), re) - cdf.begin());
    return {
        eid, {sqrt(ruv.x) * (1 - ruv.y), 1 - sqrt(ruv.x), ruv.y * sqrt(ruv.x)}};
}

/// Compute a distribution for sampling quad meshes uniformly
inline vector<float> sample_quads_cdf(
    const vector<vec4i>& quads, const vector<vec3f>& pos) {
    auto cdf = vector<float>(quads.size());
    for (auto i = 0; i < quads.size(); i++)
        cdf[i] = quad_area(
            pos[quads[i].x], pos[quads[i].y], pos[quads[i].z], pos[quads[i].w]);
    for (auto i = 1; i < quads.size(); i++) cdf[i] += cdf[i - 1];
    return cdf;
}

/// Pick a point on a quad mesh
inline pair<int, vec4f> sample_quads(
    const vector<float>& cdf, float re, const vec2f& ruv) {
    if (ruv.x < 0.5f) {
        auto eid = 0;
        auto euv = zero3f;
        std::tie(eid, euv) = sample_triangles(cdf, re, {ruv.x * 2, ruv.y});
        return {eid, {euv.x, euv.y, 0, euv.z}};
    } else {
        auto eid = 0;
        auto euv = zero3f;
        std::tie(eid, euv) =
            sample_triangles(cdf, re, {(ruv.x - 0.5f) * 2, ruv.y});
        return {eid, {0, euv.z, euv.x, euv.y}};
    }
}

/// Samples a set of points over a triangle mesh uniformly. The rng function
/// takes the point index and returns vec3f numbers uniform directibuted in
/// [0,1]^3. unorm and texcoord are optional.
inline tuple<vector<vec3f>, vector<vec3f>, vector<vec2f>>
sample_triangles_points(const vector<vec3i>& triangles,
    const vector<vec3f>& pos, const vector<vec3f>& norm,
    const vector<vec2f>& texcoord, int npoints, uint64_t seed = 0) {
    auto sampled_pos = vector<vec3f>(npoints);
    auto sampled_norm = vector<vec3f>(norm.empty() ? 0 : npoints);
    auto sampled_texcoord = vector<vec2f>(texcoord.empty() ? 0 : npoints);
    auto cdf = sample_triangles_cdf(triangles, pos);
    auto rng = init_rng(seed);
    for (auto i = 0; i < npoints; i++) {
        auto eid = 0;
        auto euv = zero3f;
        std::tie(eid, euv) = sample_triangles(
            cdf, next_rand1f(rng), {next_rand1f(rng), next_rand1f(rng)});
        auto t = triangles[eid];
        sampled_pos[i] = pos[t.x] * euv.x + pos[t.y] * euv.y + pos[t.z] * euv.z;
        if (!sampled_norm.empty())
            sampled_norm[i] = normalize(
                norm[t.x] * euv.x + norm[t.y] * euv.y + norm[t.z] * euv.z);
        if (!sampled_texcoord.empty())
            sampled_texcoord[i] = texcoord[t.x] * euv.x +
                                  texcoord[t.y] * euv.y + texcoord[t.z] * euv.z;
    }

    return {sampled_pos, sampled_norm, sampled_texcoord};
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// EXAMPLE SHAPES
// -----------------------------------------------------------------------------
namespace ygl {

/// Make a sphere. Returns quads, pos.
tuple<vector<vec3i>, vector<vec3f>> make_sphere(int level);

/// Make a geodesic sphere. Returns quads, pos.
tuple<vector<vec3i>, vector<vec3f>> make_geodesicsphere(int level);

/// Make a cube with unique vertices. This is watertight but has no
/// texture coordinates or normals. Returns quads, pos.
tuple<vector<vec4i>, vector<vec3f>> make_cube();

/// Make a sphere. This is not watertight. Returns quads, pos, norm, texcoord.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>> make_uvsphere(
    int level, bool flipped = false);

/// Make a sphere. This is not watertight. Returns quads, pos, norm, texcoord.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>
make_uvhemisphere(int level, bool flipped = false);

/// Make a quad. Returns quads, pos, norm, texcoord.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>> make_uvquad(
    int level);

/// Make a facevarying sphere with unique vertices but different texture
/// coordinates. Returns (quads, pos), (quads, norm), (quads, texcoord).
tuple<vector<vec4i>, vector<vec3f>, vector<vec4i>, vector<vec3f>, vector<vec4i>,
    vector<vec2f>>
make_fvsphere();

/// Make a facevarying cube with unique vertices but different texture
/// coordinates. Returns (quads, pos), (quads, norm), (quads, texcoord).
tuple<vector<vec4i>, vector<vec3f>, vector<vec4i>, vector<vec3f>, vector<vec4i>,
    vector<vec2f>>
make_fvcube();

/// Make a suzanne monkey model for testing. Note that some quads are
/// degenerate. Returns quads, pos.
tuple<vector<vec4i>, vector<vec3f>> make_suzanne();

/// Make a cube with uv. This is not watertight. Returns quads, pos, norm,
/// texcoord.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>> make_uvcube(
    int level);

/// Make a sphere from a cube. This is not watertight. Returns quads, pos, norm,
/// texcoord.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>
make_uvspherecube(int level);

/// Make a cube than stretch it towards a sphere. This is not watertight.
/// Returns quads, pos, norm, texcoord.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>
make_uvspherizedcube(int level, float radius);

/// Make a flipped sphere. This is not watertight. Returns quads, pos, norm,
/// texcoord.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>
make_uvflipcapsphere(int level, float z, bool flipped = false);

/// Make a cutout sphere. This is not watertight. Returns quads, pos, norm,
/// texcoord.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>
make_uvcutsphere(int level, float z, bool flipped = false);

/// Make a hair ball around a shape. Returns lines, pos, norm, texcoord, radius.
tuple<vector<vec2i>, vector<vec3f>, vector<vec3f>, vector<vec2f>, vector<float>>
make_hair(int num, int level, const vec2f& len, const vec2f& rad,
    const vector<vec3i>& striangles, const vector<vec4i>& squads,
    const vector<vec3f>& spos, const vector<vec3f>& snorm,
    const vector<vec2f>& stexcoord, const vec2f& noise = zero2f,
    const vec2f& clump = zero2f, const vec2f& rotation = zero2f,
    uint32_t seed = 0);

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMAGE CONTAINERS
// -----------------------------------------------------------------------------
namespace ygl {

/// HDR image
struct image4f {
    /// empty image constructor
    image4f() : _w{0}, _h{0}, _d{} {}
    /// image constructor
    image4f(int w, int h, const vec4f& v = zero4f)
        : _w{w}, _h{h}, _d(size_t(w * h), v) {}
    /// image constructor
    image4f(int w, int h, const vec4f* v) : _w{w}, _h{h}, _d(v, v + w * h) {}

    /// width
    int width() const { return _w; }
    /// height
    int height() const { return _h; }
    /// size
    vec2i size() const { return {_w, _h}; }
    /// check for empty
    bool empty() const { return _w == 0 || _h == 0; }
    /// check for empty
    explicit operator bool() const { return _w != 0 && _h != 0; }

    /// reallocate memory
    void resize(int w, int h, const vec4f& v = zero4f) {
        _w = w;
        _h = h;
        _d.resize(_w * _h);
    }
    /// reallocate memory
    void assign(int w, int h, const vec4f& v) {
        _w = w;
        _h = h;
        _d.assign(_w * _h, v);
    }

    /// set values
    void set(const vec4f& v) { _d.assign(_w * _h, v); }

    /// element access
    vec4f& operator[](const vec2i& ij) { return _d[ij.y * _w + ij.x]; }
    /// element access
    const vec4f& operator[](const vec2i& ij) const {
        return _d[ij.y * _w + ij.x];
    }
    /// element access
    vec4f& at(const vec2i& ij) { return _d.at(ij.y * _w + ij.x); }
    /// element access
    const vec4f& at(const vec2i& ij) const { return _d.at(ij.y * _w + ij.x); }
    /// element access
    vec4f& at(int i, int j) { return _d.at(j * _w + i); }
    /// element access
    const vec4f& at(int i, int j) const { return _d.at(j * _w + i); }

    /// data access
    vec4f* data() { return _d.data(); }
    /// data access
    const vec4f* data() const { return _d.data(); }

   private:
    int _w, _h;
    vector<vec4f> _d;
};

/// LDR image
struct image4b {
    /// empty image constructor
    image4b() : _w{0}, _h{0}, _d{} {}
    /// image constructor
    image4b(int w, int h, const vec4b& v = zero4b)
        : _w{w}, _h{h}, _d(size_t(w * h), v) {}
    /// image constructor
    image4b(int w, int h, const vec4b* v) : _w{w}, _h{h}, _d(v, v + w * h) {}

    /// width
    int width() const { return _w; }
    /// height
    int height() const { return _h; }
    /// size
    vec2i size() const { return {_w, _h}; }
    /// check for empty
    bool empty() const { return _w == 0 || _h == 0; }
    /// check for empty
    explicit operator bool() const { return _w != 0 && _h != 0; }

    /// reallocate memory
    void resize(int w, int h, const vec4b& v = zero4b) {
        _w = w;
        _h = h;
        _d.resize(_w * _h);
    }
    /// reallocate memory
    void assign(int w, int h, const vec4b& v) {
        _w = w;
        _h = h;
        _d.assign(_w * _h, v);
    }

    /// set values
    void set(const vec4b& v) { _d.assign(_w * _h, v); }

    /// element access
    vec4b& operator[](const vec2i& ij) { return _d[ij.y * _w + ij.x]; }
    /// element access
    const vec4b& operator[](const vec2i& ij) const {
        return _d[ij.y * _w + ij.x];
    }
    /// element access
    vec4b& at(const vec2i& ij) { return _d.at(ij.y * _w + ij.x); }
    /// element access
    const vec4b& at(const vec2i& ij) const { return _d.at(ij.y * _w + ij.x); }
    /// element access
    vec4b& at(int i, int j) { return _d.at(j * _w + i); }
    /// element access
    const vec4b& at(int i, int j) const { return _d.at(j * _w + i); }

    /// data access
    vec4b* data() { return _d.data(); }
    /// data access
    const vec4b* data() const { return _d.data(); }

   private:
    int _w, _h;
    vector<vec4b> _d;
};

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMAGE OPERATIONS
// -----------------------------------------------------------------------------
namespace ygl {

/// Approximate conversion from srgb.
inline vec3f srgb_to_linear(const vec3b& srgb) {
    return {pow(byte_to_float(srgb.x), 2.2f), pow(byte_to_float(srgb.y), 2.2f),
        pow(byte_to_float(srgb.z), 2.2f)};
}
/// Approximate conversion from srgb.
inline vec4f srgb_to_linear(const vec4b& srgb) {
    return {pow(byte_to_float(srgb.x), 2.2f), pow(byte_to_float(srgb.y), 2.2f),
        pow(byte_to_float(srgb.z), 2.2f), byte_to_float(srgb.w)};
}
/// Approximate conversion to srgb.
inline vec3b linear_to_srgb(const vec3f& lin) {
    return {float_to_byte(pow(lin.x, 1 / 2.2f)),
        float_to_byte(pow(lin.y, 1 / 2.2f)),
        float_to_byte(pow(lin.z, 1 / 2.2f))};
}
/// Approximate conversion to srgb.
inline vec4b linear_to_srgb(const vec4f& lin) {
    return {float_to_byte(pow(lin.x, 1 / 2.2f)),
        float_to_byte(pow(lin.y, 1 / 2.2f)),
        float_to_byte(pow(lin.z, 1 / 2.2f)), float_to_byte(lin.w)};
}

/// Convert between CIE XYZ and xyY
inline vec3f xyz_to_xyY(const vec3f& xyz) {
    if (xyz == zero3f) return zero3f;
    return {xyz.x / (xyz.x + xyz.y + xyz.z), xyz.y / (xyz.x + xyz.y + xyz.z),
        xyz.y};
}
/// Convert between CIE XYZ and xyY
inline vec3f xyY_to_xyz(const vec3f& xyY) {
    if (xyY.y == 0) return zero3f;
    return {xyY.x * xyY.z / xyY.y, xyY.z, (1 - xyY.x - xyY.y) * xyY.z / xyY.y};
}
/// Convert between CIE XYZ and RGB
inline vec3f xyz_to_rgb(const vec3f& xyz) {
    // from http://www.brucelindbloom.com/index.html?Eqn_RGB_to_XYZ.html
    if (xyz == zero3f) return zero3f;
    return {+3.2404542f * xyz.x - 1.5371385f * xyz.y - 0.4985314f * xyz.z,
        -0.9692660f * xyz.x + 1.8760108f * xyz.y + 0.0415560f * xyz.z,
        +0.0556434f * xyz.x - 0.2040259f * xyz.y + 1.0572252f * xyz.z};
}
/// Convert between CIE XYZ and RGB
inline vec3f rgb_to_xyz(const vec3f& rgb) {
    // from http://www.brucelindbloom.com/index.html?Eqn_RGB_to_XYZ.html
    if (rgb == zero3f) return zero3f;
    return {0.4124564f * rgb.x + 0.3575761f * rgb.y + 0.1804375f * rgb.z,
        0.2126729f * rgb.x + 0.7151522f * rgb.y + 0.0721750f * rgb.z,
        0.0193339f * rgb.x + 0.1191920f * rgb.y + 0.9503041f * rgb.z};
}

#if 1
/// Tone map with a fitted filmic curve.
///
/// Implementation from
/// https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
inline float tonemap_filmic(float hdr) {
    // rescale
    auto x = hdr * 2.05f;
    // fitted values
    float a = 2.51f, b = 0.03f, c = 2.43f, d = 0.59f, e = 0.14f;
    auto y = ((x * (a * x + b)) / (x * (c * x + d) + e));
    return pow(clamp(y, 0.0f, 1.0f), 1 / 2.2f);
}
#else
inline float tonemap_filmic(float x) {
    auto y =
        (x * (x * (x * (x * 2708.7142 + 6801.1525) + 1079.5474) + 1.1614649) -
            0.00004139375) /
        (x * (x * (x * (x * 983.38937 + 4132.0662) + 2881.6522) + 128.35911) +
            1.0);
    return (float)std::max(y, 0.0);
}
#endif

/// Tone mapping HDR to LDR images.
inline image4b tonemap_image(
    const image4f& hdr, float exposure, float gamma, bool filmic = false) {
    auto ldr = image4b(hdr.width(), hdr.height());
    auto scale = pow(2.0f, exposure);
    for (auto j = 0; j < hdr.height(); j++) {
        for (auto i = 0; i < hdr.width(); i++) {
            auto h = hdr[{i, j}];
            h.xyz() *= scale;
            if (filmic) {
                h.xyz() = {tonemap_filmic(h.x), tonemap_filmic(h.y),
                    tonemap_filmic(h.z)};
            } else {
                h.xyz() = {pow(h.x, 1 / gamma), pow(h.y, 1 / gamma),
                    pow(h.z, 1 / gamma)};
            }
            ldr[{i, j}] = float_to_byte(h);
        }
    }
    return ldr;
}

/// Image over operator
inline void image_over(
    vec4f* img, int width, int height, int nlayers, vec4f** layers) {
    for (auto i = 0; i < width * height; i++) {
        img[i] = {0, 0, 0, 0};
        auto weight = 1.0f;
        for (auto l = 0; l < nlayers; l++) {
            img[i].x += layers[l][i].x * layers[l][i].w * weight;
            img[i].y += layers[l][i].y * layers[l][i].w * weight;
            img[i].z += layers[l][i].z * layers[l][i].w * weight;
            img[i].w += layers[l][i].w * weight;
            weight *= (1 - layers[l][i].w);
        }
        if (img[i].w) {
            img[i].x /= img[i].w;
            img[i].y /= img[i].w;
            img[i].z /= img[i].w;
        }
    }
}

/// Image over operator
inline void image_over(
    vec4b* img, int width, int height, int nlayers, vec4b** layers) {
    for (auto i = 0; i < width * height; i++) {
        auto comp = zero4f;
        auto weight = 1.0f;
        for (auto l = 0; l < nlayers && weight > 0; l++) {
            auto w = byte_to_float(layers[l][i].w);
            comp.x += byte_to_float(layers[l][i].x) * w * weight;
            comp.y += byte_to_float(layers[l][i].y) * w * weight;
            comp.z += byte_to_float(layers[l][i].z) * w * weight;
            comp.w += w * weight;
            weight *= (1 - w);
        }
        if (comp.w) {
            img[i].x = float_to_byte(comp.x / comp.w);
            img[i].y = float_to_byte(comp.y / comp.w);
            img[i].z = float_to_byte(comp.z / comp.w);
            img[i].w = float_to_byte(comp.w);
        } else {
            img[i] = {0, 0, 0, 0};
        }
    }
}

/// Convert HSV to RGB
///
/// Implementatkion from
/// http://stackoverflow.com/questions/3018313/algorithm-to-convert-rgb-to-hsv-and-hsv-to-rgb-in-range-0-255-for-both
inline vec4b hsv_to_rgb(const vec4b& hsv) {
    vec4b rgb = {0, 0, 0, hsv.w};
    byte region, remainder, p, q, t;

    byte h = hsv.x, s = hsv.y, v = hsv.z;

    if (s == 0) {
        rgb.x = v;
        rgb.y = v;
        rgb.z = v;
        return rgb;
    }

    region = h / 43;
    remainder = (h - (region * 43)) * 6;

    p = (v * (255 - s)) >> 8;
    q = (v * (255 - ((s * remainder) >> 8))) >> 8;
    t = (v * (255 - ((s * (255 - remainder)) >> 8))) >> 8;

    switch (region) {
        case 0:
            rgb.x = v;
            rgb.y = t;
            rgb.z = p;
            break;
        case 1:
            rgb.x = q;
            rgb.y = v;
            rgb.z = p;
            break;
        case 2:
            rgb.x = p;
            rgb.y = v;
            rgb.z = t;
            break;
        case 3:
            rgb.x = p;
            rgb.y = q;
            rgb.z = v;
            break;
        case 4:
            rgb.x = t;
            rgb.y = p;
            rgb.z = v;
            break;
        default:
            rgb.x = v;
            rgb.y = p;
            rgb.z = q;
            break;
    }

    return rgb;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// EXAMPLE IMAGES
// -----------------------------------------------------------------------------
namespace ygl {

/// Make a grid image
inline image4b make_grid_image(int width, int height, int tile = 64,
    const vec4b& c0 = {90, 90, 90, 255},
    const vec4b& c1 = {128, 128, 128, 255}) {
    image4b pixels(width, height);
    for (int j = 0; j < width; j++) {
        for (int i = 0; i < height; i++) {
            auto c = i % tile == 0 || i % tile == tile - 1 || j % tile == 0 ||
                     j % tile == tile - 1;
            pixels.at(i, j) = (c) ? c0 : c1;
        }
    }
    return pixels;
}

/// Make a checkerboard image
inline image4b make_checker_image(int width, int height, int tile = 64,
    const vec4b& c0 = {90, 90, 90, 255},
    const vec4b& c1 = {128, 128, 128, 255}) {
    image4b pixels(width, height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto c = (i / tile + j / tile) % 2 == 0;
            pixels.at(i, j) = (c) ? c0 : c1;
        }
    }
    return pixels;
}

/// Make an image with bumps and dimples.
inline image4b make_bumpdimple_image(int width, int height, int tile = 64) {
    image4b pixels(width, height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto c = (i / tile + j / tile) % 2 == 0;
            auto ii = i % tile - tile / 2, jj = j % tile - tile / 2;
            auto r =
                sqrt(float(ii * ii + jj * jj)) / sqrt(float(tile * tile) / 4);
            auto h = 0.5f;
            if (r < 0.5f) { h += (c) ? (0.5f - r) : -(0.5f - r); }
            auto g = float_to_byte(h);
            pixels.at(i, j) = vec4b{g, g, g, 255};
        }
    }
    return pixels;
}

/// Make a uv colored grid
inline image4b make_ramp_image(int width, int height, const vec4b& c0,
    const vec4b& c1, bool srgb = false) {
    image4b pixels(width, height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto u = (float)i / (float)width;
            if (srgb) {
                pixels.at(i, j) = linear_to_srgb(
                    srgb_to_linear(c0) * (1 - u) + srgb_to_linear(c1) * u);
            } else {
                pixels.at(i, j) = float_to_byte(
                    byte_to_float(c0) * (1 - u) + byte_to_float(c1) * u);
            }
        }
    }
    return pixels;
}

/// Make a gamma ramp image
inline image4b make_gammaramp_image(int width, int height) {
    image4b pixels(width, height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto u = j / float(height - 1);
            if (i < width / 3) u = pow(u, 2.2f);
            if (i > (width * 2) / 3) u = pow(u, 1 / 2.2f);
            auto c = (unsigned char)(u * 255);
            pixels.at(i, j) = {c, c, c, 255};
        }
    }
    return pixels;
}

/// Make a gamma ramp image
inline image4f make_gammaramp_imagef(int width, int height) {
    image4f pixels(width, height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto u = j / float(height - 1);
            if (i < width / 3) u = pow(u, 2.2f);
            if (i > (width * 2) / 3) u = pow(u, 1 / 2.2f);
            pixels.at(i, j) = {u, u, u, 1};
        }
    }
    return pixels;
}

/// Make an image color with red/green in the [0,1] range. Helpful to visualize
/// uv texture coordinate application.
inline image4b make_uv_image(int width, int height) {
    image4b pixels(width, height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto r = float_to_byte(i / (float)(width - 1));
            auto g = float_to_byte(j / (float)(height - 1));
            pixels.at(i, j) = vec4b{r, g, 0, 255};
        }
    }
    return pixels;
}

/// Make a uv colored grid
inline image4b make_uvgrid_image(
    int width, int height, int tile = 64, bool colored = true) {
    image4b pixels(width, height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            byte ph = 32 * (i / (height / 8));
            byte pv = 128;
            byte ps = 64 + 16 * (7 - j / (height / 8));
            if (i % (tile / 2) && j % (tile / 2)) {
                if ((i / tile + j / tile) % 2)
                    pv += 16;
                else
                    pv -= 16;
            } else {
                pv = 196;
                ps = 32;
            }
            pixels.at(i, j) = (colored) ? hsv_to_rgb({ph, ps, pv, 255}) :
                                          vec4b{pv, pv, pv, 255};
        }
    }
    return pixels;
}

/// Make a uv recusive colored grid
inline image4b make_recuvgrid_image(
    int width, int height, int tile = 64, bool colored = true) {
    image4b pixels(width, height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            byte ph = 32 * (i / (height / 8));
            byte pv = 128;
            byte ps = 64 + 16 * (7 - j / (height / 8));
            if (i % (tile / 2) && j % (tile / 2)) {
                if ((i / tile + j / tile) % 2)
                    pv += 16;
                else
                    pv -= 16;
                if ((i / (tile / 4) + j / (tile / 4)) % 2)
                    pv += 4;
                else
                    pv -= 4;
                if ((i / (tile / 8) + j / (tile / 8)) % 2)
                    pv += 1;
                else
                    pv -= 1;
            } else {
                pv = 196;
                ps = 32;
            }
            pixels.at(i, j) = (colored) ? hsv_to_rgb({ph, ps, pv, 255}) :
                                          vec4b{pv, pv, pv, 255};
        }
    }
    return pixels;
}

/// Comvert a bump map to a normal map.
inline image4b bump_to_normal_map(const image4b& img, float scale = 1) {
    image4b norm(img.width(), img.height());
    for (int j = 0; j < img.height(); j++) {
        for (int i = 0; i < img.width(); i++) {
            auto i1 = (i + 1) % img.width(), j1 = (j + 1) % img.height();
            auto p00 = img.at(i, j), p10 = img.at(i1, j), p01 = img.at(i, j1);
            auto g00 = (float(p00.x) + float(p00.y) + float(p00.z)) / (3 * 255);
            auto g01 = (float(p01.x) + float(p01.y) + float(p01.z)) / (3 * 255);
            auto g10 = (float(p10.x) + float(p10.y) + float(p10.z)) / (3 * 255);
            auto n = vec3f{scale * (g00 - g10), scale * (g00 - g01), 1.0f};
            n = normalize(n) * 0.5f + vec3f{0.5f, 0.5f, 0.5f};
            auto c =
                vec4b{byte(n.x * 255), byte(n.y * 255), byte(n.z * 255), 255};
            norm.at(i, j) = c;
        }
    }
    return norm;
}

/// Make a sunsky HDR model with sun at theta elevation in [0,pi/2], turbidity
/// in [1.7,10] with or without sun.
image4f make_sunsky_image(int res, float thetaSun, float turbidity = 3,
    bool has_sun = false, bool has_ground = true);

/// Compute the revised Pelin noise function. Wrap provides a wrapping noise
/// but must be power of two (wraps at 256 anyway). For octave based noise,
/// good values are obtained with octaves=6 (numerber of noise calls),
/// lacunarity=~2.0 (spacing between successive octaves: 2.0 for warpping
/// output), gain=0.5 (relative weighting applied to each successive octave),
/// offset=1.0 (used to invert the ridges).

/// Make a noise image. Wrap works only if both resx and resy are powers of two.
image4b make_noise_image(int resx, int resy, float scale = 1, bool wrap = true);
/// Make a noise image. Wrap works only if both resx and resy are powers of two.
image4b make_fbm_image(int resx, int resy, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, int octaves = 6, bool wrap = true);
/// Make a noise image. Wrap works only if both resx and resy are powers of two.
image4b make_ridge_image(int resx, int resy, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, float offset = 1.0f,
    int octaves = 6, bool wrap = true);
/// Make a noise image. Wrap works only if both resx and resy are powers of two.
image4b make_turbulence_image(int resx, int resy, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, int octaves = 6, bool wrap = true);

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMAGE LOADING/SAVING
// -----------------------------------------------------------------------------
namespace ygl {

#if YGL_IMAGEIO

/// Check if an image is HDR based on filename
bool is_hdr_filename(const string& filename);

/// Loads an ldr image.
image4b load_image4b(const string& filename);

/// Loads an hdr image.
image4f load_image4f(const string& filename);

/// Saves an ldr image.
bool save_image4b(const string& filename, const image4b& img);

/// Saves an hdr image.
bool save_image4f(const string& filename, const image4f& img);

/// Loads an image
vector<float> load_imagef(
    const string& filename, int& width, int& height, int& ncomp);

/// Loads an image
vector<byte> load_image(
    const string& filename, int& width, int& height, int& ncomp);

/// Loads an image from memory.
vector<float> load_imagef_from_memory(const string& filename, const byte* data,
    int length, int& width, int& height, int& ncomp);

/// Loads an image from memory.
vector<byte> load_image_from_memory(const string& filename, const byte* data,
    int length, int& width, int& height, int& ncomp);

/// Saves an image
bool save_imagef(
    const string& filename, int width, int height, int ncomp, const float* hdr);

/// Saves an image
bool save_image(
    const string& filename, int width, int height, int ncomp, const byte* ldr);

/// Save an HDR or LDR image with tonemapping based on filename
inline bool save_image(const string& filename, const image4f& hdr,
    float exposure, float gamma, bool filmic = false) {
    if (is_hdr_filename(filename)) {
        return save_image4f(filename, hdr);
    } else {
        auto ldr = tonemap_image(hdr, exposure, gamma, filmic);
        return save_image4b(filename, ldr);
    }
}

/// Filter for resizing
enum struct resize_filter {
    /// default
    def = 0,
    /// box filter
    box = 1,
    /// triangle filter
    triangle = 2,
    /// cubic spline
    cubic_spline = 3,
    /// Catmull-Rom interpolating sline
    catmull_rom = 4,
    /// Mitchel-Netrevalli filter with B=1/3, C=1/3
    mitchell = 5
};

/// Edge mode for resizing
enum struct resize_edge {
    /// default
    def = 0,
    /// clamp
    clamp = 1,
    /// reflect
    reflect = 2,
    /// wrap
    wrap = 3,
    /// zero
    zero = 4
};

/// Resize image.
void resize_image(const image4f& img, image4f& res_img,
    resize_filter filter = resize_filter::def,
    resize_edge edge = resize_edge::def, bool premultiplied_alpha = false);

/// Resize image.
void resize_image(const image4b& img, image4b& res_img,
    resize_filter filter = resize_filter::def,
    resize_edge edge = resize_edge::def, bool premultiplied_alpha = false);

#endif

}  // namespace ygl

// -----------------------------------------------------------------------------
// RAY-PRIMITIVE INTERSECTION FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

/// Intersect a ray with a point (approximate)
///
/// Parameters:
/// - ray: ray origin and direction, parameter min, max range
/// - p: point position
/// - r: point radius
///
/// Out Parameters:
/// - ray_t: ray parameter at the intersection point
/// - euv: primitive uv ( {0,0} for points )
///
/// Returns:
/// - whether the intersection occurred
///
/// Iplementation Notes:
/// - out Parameters and only writtent o if an intersection occurs
/// - algorithm finds the closest point on the ray segment to the point and
///    test their distance with the point radius
/// - based on http://geomalgorithms.com/a02-lines.html.
inline bool intersect_point(
    const ray3f& ray, const vec3f& p, float r, float& ray_t) {
    // find parameter for line-point minimum distance
    auto w = p - ray.o;
    auto t = dot(w, ray.d) / dot(ray.d, ray.d);

    // exit if not within bounds
    if (t < ray.tmin || t > ray.tmax) return false;

    // test for line-point distance vs point radius
    auto rp = ray.o + ray.d * t;
    auto prp = p - rp;
    if (dot(prp, prp) > r * r) return false;

    // intersection occurred: set params and exit
    ray_t = t;

    return true;
}

/// Intersect a ray with a line
///
/// Parameters:
/// - ray: ray origin and direction, parameter min, max range
/// - v0, v1: line segment points
/// - r0, r1: line segment radia
///
/// Out Parameters:
/// - ray_t: ray parameter at the intersection point
/// - euv: euv.x is the line parameter at the intersection ( euv.y is zero )
///
/// Returns:
/// - whether the intersection occurred
///
/// Notes:
/// - out Parameters and only writtent o if an intersection occurs
/// - algorithm find the closest points on line and ray segment and test
///   their distance with the line radius at that location
/// - based on http://geomalgorithms.com/a05-intersect-1.html
/// - based on http://geomalgorithms.com/a07-distance.html#
///     dist3D_Segment_to_Segment
inline bool intersect_line(const ray3f& ray, const vec3f& v0, const vec3f& v1,
    float r0, float r1, float& ray_t, vec2f& euv) {
    // setup intersection params
    auto u = ray.d;
    auto v = v1 - v0;
    auto w = ray.o - v0;

    // compute values to solve a linear system
    auto a = dot(u, u);
    auto b = dot(u, v);
    auto c = dot(v, v);
    auto d = dot(u, w);
    auto e = dot(v, w);
    auto det = a * c - b * b;

    // check determinant and exit if lines are parallel
    // (could use EPSILONS if desired)
    if (det == 0) return false;

    // compute Parameters on both ray and segment
    auto t = (b * e - c * d) / det;
    auto s = (a * e - b * d) / det;

    // exit if not within bounds
    if (t < ray.tmin || t > ray.tmax) return false;

    // clamp segment param to segment corners
    s = clamp(s, (float)0, (float)1);

    // compute segment-segment distance on the closest points
    auto p0 = ray.o + ray.d * t;
    auto p1 = v0 + (v1 - v0) * s;
    auto p01 = p0 - p1;

    // check with the line radius at the same point
    auto r = r0 * (1 - s) + r1 * s;
    if (dot(p01, p01) > r * r) return false;

    // intersection occurred: set params and exit
    ray_t = t;
    euv = {1 - s, s};

    return true;
}

/// Intersect a ray with a triangle
///
/// Parameters:
/// - ray: ray origin and direction, parameter min, max range
/// - v0, v1, v2: triangle vertices
///
/// Out Parameters:
/// - ray_t: ray parameter at the intersection point
/// - euv: baricentric coordinates of the intersection
///
/// Returns:
/// - whether the intersection occurred
///
/// Notes:
/// - out Parameters and only writtent o if an intersection occurs
/// - algorithm based on Muller-Trombone intersection test
inline bool intersect_triangle(const ray3f& ray, const vec3f& v0,
    const vec3f& v1, const vec3f& v2, float& ray_t, vec3f& euv) {
    // compute triangle edges
    auto edge1 = v1 - v0;
    auto edge2 = v2 - v0;

    // compute determinant to solve a linear system
    auto pvec = cross(ray.d, edge2);
    auto det = dot(edge1, pvec);

    // check determinant and exit if triangle and ray are parallel
    // (could use EPSILONS if desired)
    if (det == 0) return false;
    auto inv_det = 1.0f / det;

    // compute and check first bricentric coordinated
    auto tvec = ray.o - v0;
    auto u = dot(tvec, pvec) * inv_det;
    if (u < 0 || u > 1) return false;

    // compute and check second bricentric coordinated
    auto qvec = cross(tvec, edge1);
    auto v = dot(ray.d, qvec) * inv_det;
    if (v < 0 || u + v > 1) return false;

    // compute and check ray parameter
    auto t = dot(edge2, qvec) * inv_det;
    if (t < ray.tmin || t > ray.tmax) return false;

    // intersection occurred: set params and exit
    ray_t = t;
    euv = {1 - u - v, u, v};

    return true;
}

/// Intersect a ray with a quad represented as two triangles (0,1,3) and
/// (2,3,1), with the uv coordinates of the second triangle corrected by u =
/// 1-u' and v = 1-v' to produce a quad parametrization where u and v go from 0
/// to 1. This is equivalent to Intel's Embree. The external user does not have
/// to be concerned about the parametrization and can just use the euv as
/// specified.
///
/// Parameters:
/// - ray: ray origin and direction, parameter min, max range
/// - v0, v1, v2, v3: quad vertices
///
/// Out Parameters:
/// - ray_t: ray parameter at the intersection point
/// - euv: baricentric coordinates of the intersection
///
/// Returns:
/// - whether the intersection occurred
inline bool intersect_quad(const ray3f& ray, const vec3f& v0, const vec3f& v1,
    const vec3f& v2, const vec3f& v3, float& ray_t, vec4f& euv) {
    auto hit = false;
    auto tray = ray;
    if (intersect_triangle(tray, v0, v1, v3, ray_t, (vec3f&)euv)) {
        euv = {euv.x, euv.y, 0, euv.z};
        tray.tmax = ray_t;
        hit = true;
    }
    if (intersect_triangle(tray, v2, v3, v1, ray_t, (vec3f&)euv)) {
        euv = {0, 1 - euv.y, euv.y + euv.z - 1, 1 - euv.z};
        tray.tmax = ray_t;
        hit = true;
    }
    return hit;
}

/// Intersect a ray with a tetrahedron. Note that we consider only
/// intersection wiht the tetrahedra surface and discount intersction with
/// the interior.
///
/// Parameters:
/// - ray: ray to intersect with
/// - v0, v1, v2: triangle vertices
///
/// Out Parameters:
/// - ray_t: ray parameter at the intersection point
/// - euv: baricentric coordinates of the intersection
///
/// Returns:
/// - whether the intersection occurred
///
/// TODO: check order
/// TODO: uv
inline bool intersect_tetrahedron(const ray3f& ray_, const vec3f& v0,
    const vec3f& v1, const vec3f& v2, const vec3f& v3, float& ray_t,
    vec4f& euv) {
    // check intersction for each face
    auto hit = false;
    auto ray = ray_;
    auto tuv = zero3f;
    if (intersect_triangle(ray, v0, v1, v2, ray_t, tuv)) {
        hit = true;
        ray.tmax = ray_t;
    }
    if (intersect_triangle(ray, v0, v1, v3, ray_t, tuv)) {
        hit = true;
        ray.tmax = ray_t;
    }
    if (intersect_triangle(ray, v0, v2, v3, ray_t, tuv)) {
        hit = true;
        ray.tmax = ray_t;
    }
    if (intersect_triangle(ray, v1, v2, v3, ray_t, tuv)) {
        hit = true;
        ray.tmax = ray_t;
    }

    return hit;
}

/// Intersect a ray with a axis-aligned bounding box
///
/// Parameters:
/// - ray: ray to intersect with
/// - bbox: bounding box min/max bounds
///
/// Returns:
/// - whether the intersection occurred
inline bool intersect_check_bbox(const ray3f& ray, const bbox3f& bbox) {
    // set up convenient pointers for looping over axes
    auto tmin = ray.tmin, tmax = ray.tmax;

    // for each axis, clip intersection against the bounding planes
    for (int i = 0; i < 3; i++) {
        // determine intersection ranges
        auto invd = 1.0f / ray.d[i];
        auto t0 = (bbox.min[i] - ray.o[i]) * invd;
        auto t1 = (bbox.max[i] - ray.o[i]) * invd;
        // flip based on range directions
        if (invd < 0.0f) {
            float a = t0;
            t0 = t1;
            t1 = a;
        }
        // clip intersection
        tmin = t0 > tmin ? t0 : tmin;
        tmax = t1 < tmax ? t1 : tmax;
        // if intersection is empty, exit
        if (tmin > tmax) return false;
    }

    // passed all planes, then intersection occurred
    return true;
}

/// Min/max used in BVH traversal. Copied here since the traversal code
/// relies on the specific behaviour wrt NaNs.
static inline const float& _safemin(const float& a, const float& b) {
    return (a < b) ? a : b;
}
/// Min/max used in BVH traversal. Copied here since the traversal code
/// relies on the specific behaviour wrt NaNs.
static inline const float& _safemax(const float& a, const float& b) {
    return (a > b) ? a : b;
}

/// Intersect a ray with a axis-aligned bounding box
///
/// Parameters:
/// - ray_o, ray_d: ray origin and direction
/// - ray_tmin, ray_tmax: ray parameter min, max range
/// - ray_dinv: ray inverse direction
/// - ray_dsign: ray direction sign
/// - bbox_min, bbox_max: bounding box min/max bounds
///
/// Returns:
/// - whether the intersection occurred
///
/// Implementation Notes:
/// - based on "Robust BVH Ray Traversal" by T. Ize published at
/// http://jcgt.org/published/0002/02/02/paper.pdf
inline bool intersect_check_bbox(const ray3f& ray, const vec3f& ray_dinv,
    const vec3i& ray_dsign, const bbox3f& bbox_) {
    auto bbox = &bbox_.min;
    auto txmin = (bbox[ray_dsign.x].x - ray.o.x) * ray_dinv.x;
    auto txmax = (bbox[1 - ray_dsign.x].x - ray.o.x) * ray_dinv.x;
    auto tymin = (bbox[ray_dsign.y].y - ray.o.y) * ray_dinv.y;
    auto tymax = (bbox[1 - ray_dsign.y].y - ray.o.y) * ray_dinv.y;
    auto tzmin = (bbox[ray_dsign.z].z - ray.o.z) * ray_dinv.z;
    auto tzmax = (bbox[1 - ray_dsign.z].z - ray.o.z) * ray_dinv.z;
    auto tmin = _safemax(tzmin, _safemax(tymin, _safemax(txmin, ray.tmin)));
    auto tmax = _safemin(tzmax, _safemin(tymax, _safemin(txmax, ray.tmax)));
    tmax *= 1.00000024f;  // for double: 1.0000000000000004
    return tmin <= tmax;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// POINT-PRIMITIVE DISTANCE FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// TODO: documentation
inline bool overlap_point(
    const vec3f& pos, float dist_max, const vec3f& p, float r, float& dist) {
    auto d2 = dot(pos - p, pos - p);
    if (d2 > (dist_max + r) * (dist_max + r)) return false;
    dist = sqrt(d2);
    return true;
}

// TODO: documentation
inline vec2f closestuv_line(
    const vec3f& pos, const vec3f& v0, const vec3f& v1) {
    auto ab = v1 - v0;
    auto d = dot(ab, ab);
    // Project c onto ab, computing parameterized position d(t) = a + t*(b –
    // a)
    auto u = dot(pos - v0, ab) / d;
    u = clamp(u, (float)0, (float)1);
    return {1 - u, u};
}

// TODO: documentation
inline bool overlap_line(const vec3f& pos, float dist_max, const vec3f& v0,
    const vec3f& v1, float r0, float r1, float& dist, vec2f& euv) {
    auto uv = closestuv_line(pos, v0, v1);
    // Compute projected position from the clamped t d = a + t * ab;
    auto p = lerp(v0, v1, uv.y);
    auto r = lerp(r0, r1, uv.y);
    auto d2 = dot(pos - p, pos - p);
    // check distance
    if (d2 > (dist_max + r) * (dist_max + r)) return false;
    // done
    dist = sqrt(d2);
    euv = uv;
    return true;
}

// TODO: documentation
// this is a complicated test -> I probably prefer to use a sequence of test
// (triangle body, and 3 edges)
inline vec3f closestuv_triangle(
    const vec3f& pos, const vec3f& v0, const vec3f& v1, const vec3f& v2) {
    auto ab = v1 - v0;
    auto ac = v2 - v0;
    auto ap = pos - v0;

    auto d1 = dot(ab, ap);
    auto d2 = dot(ac, ap);

    // corner and edge cases
    if (d1 <= 0 && d2 <= 0) return vec3f{1, 0, 0};

    auto bp = pos - v1;
    auto d3 = dot(ab, bp);
    auto d4 = dot(ac, bp);
    if (d3 >= 0 && d4 <= d3) return vec3f{0, 1, 0};

    auto vc = d1 * d4 - d3 * d2;
    if ((vc <= 0) && (d1 >= 0) && (d3 <= 0))
        return vec3f{1 - d1 / (d1 - d3), d1 / (d1 - d3), 0};

    auto cp = pos - v2;
    auto d5 = dot(ab, cp);
    auto d6 = dot(ac, cp);
    if (d6 >= 0 && d5 <= d6) return vec3f{0, 0, 1};

    auto vb = d5 * d2 - d1 * d6;
    if ((vb <= 0) && (d2 >= 0) && (d6 <= 0))
        return vec3f{1 - d2 / (d2 - d6), 0, d2 / (d2 - d6)};

    auto va = d3 * d6 - d5 * d4;
    if ((va <= 0) && (d4 - d3 >= 0) && (d5 - d6 >= 0)) {
        auto w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return vec3f{0, 1 - w, w};
    }

    // face case
    auto denom = 1 / (va + vb + vc);
    auto v = vb * denom;
    auto w = vc * denom;
    return vec3f{1 - v - w, v, w};
}

// TODO: documentation
inline bool overlap_triangle(const vec3f& pos, float dist_max, const vec3f& v0,
    const vec3f& v1, const vec3f& v2, float r0, float r1, float r2, float& dist,
    vec3f& euv) {
    auto uv = closestuv_triangle(pos, v0, v1, v2);
    auto p = v0 * uv.x + v1 * uv.y + v2 * uv.z;
    auto r = r0 * uv.x + r1 * uv.y + r2 * uv.z;
    auto dd = dot(p - pos, p - pos);
    if (dd > (dist_max + r) * (dist_max + r)) return false;
    dist = sqrt(dd);
    euv = uv;
    return true;
}

// TODO: documentation
inline bool overlap_quad(const vec3f& pos, float dist_max, const vec3f& v0,
    const vec3f& v1, const vec3f& v2, const vec3f& v3, float r0, float r1,
    float r2, float r3, float& dist, vec4f& euv) {
    auto hit = false;
    if (overlap_triangle(
            pos, dist_max, v0, v1, v3, r0, r1, r3, dist, (vec3f&)euv)) {
        euv = {euv.x, euv.y, 0, euv.z};
        dist_max = dist;
        hit = true;
    }
    if (overlap_triangle(
            pos, dist_max, v2, v3, v1, r2, r3, r1, dist, (vec3f&)euv)) {
        // dist_max = dist;
        euv = {0, 1 - euv.y, euv.y + euv.z - 1, 1 - euv.z};
        hit = true;
    }
    return hit;
}

// TODO: documentation
inline bool overlap_tetrahedron(const vec3f& pos, const vec3f& v0,
    const vec3f& v1, const vec3f& v2, const vec3f& v3, vec4f& euv) {
    auto vol = dot(v3 - v0, cross(v3 - v1, v3 - v0));
    if (vol == 0) return false;
    auto u = dot(v3 - v0, cross(v3 - v1, v3 - v0)) / vol;
    if (u < 0 || u > 1) return false;
    auto v = dot(v3 - v0, cross(v3 - v1, v3 - v0)) / vol;
    if (v < 0 || v > 1 || u + v > 1) return false;
    auto w = dot(v3 - v0, cross(v3 - v1, v3 - v0)) / vol;
    if (w < 0 || w > 1 || u + v + w > 1) return false;
    euv = {u, v, w, 1 - u - v - w};
    return true;
}

// TODO: documentation
inline bool overlap_tetrahedron(const vec3f& pos, float dist_max,
    const vec3f& v0, const vec3f& v1, const vec3f& v2, const vec3f& v3,
    float r0, float r1, float r2, float r3, float& dist, vec4f& euv) {
    // check interior
    if (overlap_tetrahedron(pos, v0, v1, v2, v3, euv)) {
        dist = 0;
        return true;
    }

    // check faces
    auto hit = false;
    auto tuv = zero3f;
    if (overlap_triangle(pos, dist_max, v0, v1, v2, r0, r1, r2, dist, tuv)) {
        hit = true;
        dist_max = dist;
    }
    if (overlap_triangle(pos, dist_max, v0, v1, v3, r0, r1, r3, dist, tuv)) {
        hit = true;
        dist_max = dist;
    }
    if (overlap_triangle(pos, dist_max, v0, v2, v3, r0, r2, r3, dist, tuv)) {
        hit = true;
        dist_max = dist;
    }
    if (overlap_triangle(pos, dist_max, v1, v2, v3, r1, r2, r3, dist, tuv)) {
        hit = true;
        // dist_max = dist;
    }

    return hit;
}

// TODO: documentation
inline bool distance_check_bbox(
    const vec3f& pos, float dist_max, const bbox3f& bbox) {
    // computing distance
    auto dd = 0.0f;

    // For each axis count any excess distance outside box extents
    for (int i = 0; i < 3; i++) {
        auto v = pos[i];
        if (v < bbox.min[i]) dd += (bbox.min[i] - v) * (bbox.min[i] - v);
        if (v > bbox.max[i]) dd += (v - bbox.max[i]) * (v - bbox.max[i]);
    }

    // check distance
    return dd < dist_max * dist_max;
}

// TODO: doc
inline bool overlap_bbox(const bbox3f& bbox1, const bbox3f& bbox2) {
    if (bbox1.max.x < bbox2.min.x || bbox1.min.x > bbox2.max.x) return false;
    if (bbox1.max.y < bbox2.min.y || bbox1.min.y > bbox2.max.y) return false;
    if (bbox1.max.z < bbox2.min.z || bbox1.min.z > bbox2.max.z) return false;
    return true;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// BVH FOR RAY INTERSECTION AND CLOSEST ELEMENT
// -----------------------------------------------------------------------------
namespace ygl {

// number of primitives to avoid splitting on
const int bvh_minprims = 4;

/// BVH tree node containing its bounds, indices to the BVH arrays of either
/// sorted primitives or internal nodes, whether its a leaf or an internal node,
/// and the split axis. Leaf and internal nodes are identical, except that
/// indices refer to primitives for leaf nodes or other nodes for internal
/// nodes. See bvh_tree for more details.
///
/// This is an internal data structure.
struct bvh_node {
    /// bounding box
    bbox3f bbox;
    /// index to the first sorted primitive/node
    uint32_t start;
    /// number of primitives/nodes
    uint16_t count;
    /// whether it is a leaf
    uint8_t isleaf;
    /// slit axis
    uint8_t axis;
};

/// BVH tree, stored as a node array. The tree structure is encoded using array
/// indices instead of pointers, both for speed but also to simplify code.
/// BVH nodes indices refer to either the node array, for internal nodes,
/// or a primitive array, for leaf nodes. BVH trees may contain only one type
/// of geometric primitive, like points, lines, triangle or shape other BVHs.
/// To handle multiple primitive types and transformed primitices, build
/// a two-level hierarchy with the outer BVH, the scene BVH, containing inner
/// BVHs, shape BVHs, each of which of a uniform primitive type.
///
/// This is an internal data structure.
struct bvh_tree {
    /// sorted array of internal nodes
    vector<bvh_node> nodes;
    /// sorted elements
    vector<int> sorted_prim;
};

// Struct that pack a bounding box, its associate primitive index, and other
// data for faster hierarchy build.
// This is internal only and should not be used externally.
struct bvh_bound_prim {
    bbox3f bbox;   // bounding box
    vec3f center;  // bounding box center (for faster sort)
    int pid;       // primitive id
};

// Comparison function for each axis
struct bvh_bound_prim_comp {
    int axis;
    float middle;

    bvh_bound_prim_comp(int a, float m = 0) : axis(a), middle(m) {}

    bool operator()(const bvh_bound_prim& a, const bvh_bound_prim& b) const {
        return a.center[axis] < b.center[axis];
    }

    bool operator()(const bvh_bound_prim& a) const {
        return a.center[axis] < middle;
    }
};

// Initializes the BVH node node that contains the primitives sorted_prims
// from start to end, by either splitting it into two other nodes,
// or initializing it as a leaf. When splitting, the heuristic heuristic is
// used and nodes added sequentially in the preallocated nodes array and
// the number of nodes nnodes is updated.
inline void make_bvh_node(bvh_node* node, vector<bvh_node>& nodes,
    bvh_bound_prim* sorted_prims, int start, int end, bool equalsize) {
    // compute node bounds
    node->bbox = invalid_bbox3f;
    for (auto i = start; i < end; i++) node->bbox += sorted_prims[i].bbox;

    // decide whether to create a leaf
    if (end - start <= bvh_minprims) {
        // makes a leaf node
        node->isleaf = true;
        node->start = start;
        node->count = end - start;
    } else {
        // choose the split axis and position
        // init to default values
        auto axis = 0;
        auto mid = (start + end) / 2;

        // compute primintive bounds and size
        auto centroid_bbox = invalid_bbox3f;
        for (auto i = start; i < end; i++)
            centroid_bbox += sorted_prims[i].center;
        auto centroid_size = bbox_diagonal(centroid_bbox);

        // check if it is not possible to split
        if (centroid_size == zero3f) {
            // we failed to split for some reasons
            node->isleaf = true;
            node->start = start;
            node->count = end - start;
        } else {
            // split along largest
            auto largest_axis = max_element(centroid_size).first;

            // check heuristic
            if (equalsize) {
                // split the space in the middle along the largest axis
                axis = largest_axis;
                mid = (int)(std::partition(sorted_prims + start,
                                sorted_prims + end,
                                bvh_bound_prim_comp(largest_axis,
                                    bbox_center(centroid_bbox)[largest_axis])) -
                            sorted_prims);
            } else {
                // balanced tree split: find the largest axis of the bounding
                // box and split along this one right in the middle
                axis = largest_axis;
                mid = (start + end) / 2;
                std::nth_element(sorted_prims + start, sorted_prims + mid,
                    sorted_prims + end, bvh_bound_prim_comp(largest_axis));
            }

            // check correctness
            assert(axis >= 0 && mid > 0);
            assert(mid > start && mid < end);

            // makes an internal node
            node->isleaf = false;
            // perform the splits by preallocating the child nodes and recurring
            node->axis = axis;
            node->start = (int)nodes.size();
            node->count = 2;
            nodes.emplace_back();
            nodes.emplace_back();
            // build child nodes
            make_bvh_node(&nodes[node->start], nodes, sorted_prims, start, mid,
                equalsize);
            make_bvh_node(&nodes[node->start + 1], nodes, sorted_prims, mid,
                end, equalsize);
        }
    }
}

/// Build a BVH from a set of primitives.
inline bvh_tree* build_bvh(
    int nprims, bool equalsize, const function<bbox3f(int)>& elem_bbox) {
    // allocate if needed
    auto bvh = new bvh_tree();

    // prepare prims
    auto bound_prims = vector<bvh_bound_prim>(nprims);
    for (auto i = 0; i < nprims; i++) {
        bound_prims[i].pid = i;
        bound_prims[i].bbox = elem_bbox(i);
        bound_prims[i].center = bbox_center(bound_prims[i].bbox);
    }

    // clear bvh
    bvh->nodes.clear();
    bvh->sorted_prim.clear();

    // allocate nodes (over-allocate now then shrink)
    bvh->nodes.reserve(nprims * 2);

    // start recursive splitting
    bvh->nodes.emplace_back();
    make_bvh_node(
        &bvh->nodes[0], bvh->nodes, bound_prims.data(), 0, nprims, equalsize);

    // shrink back
    bvh->nodes.shrink_to_fit();

    // init sorted element arrays
    // for shared memory, stored pointer to the external data
    // store the sorted primitive order for BVH walk
    bvh->sorted_prim.resize(nprims);
    for (int i = 0; i < nprims; i++) {
        bvh->sorted_prim[i] = bound_prims[i].pid;
    }

    // done
    return bvh;
}

/// Build a triangles BVH.
inline bvh_tree* build_triangles_bvh(const vector<vec3i>& triangles,
    const vector<vec3f>& pos, bool equal_size = true) {
    return build_bvh(
        (int)triangles.size(), equal_size, [&triangles, &pos](int eid) {
            auto f = triangles[eid];
            return triangle_bbox(pos[f.x], pos[f.y], pos[f.z]);
        });
}

/// Build a quads BVH.
inline bvh_tree* build_quads_bvh(const vector<vec4i>& quads,
    const vector<vec3f>& pos, bool equal_size = true) {
    return build_bvh((int)quads.size(), equal_size, [&quads, &pos](int eid) {
        auto f = quads[eid];
        return quad_bbox(pos[f.x], pos[f.y], pos[f.z], pos[f.w]);
    });
}

/// Build a lines BVH.
inline bvh_tree* build_lines_bvh(const vector<vec2i>& lines,
    const vector<vec3f>& pos, const vector<float>& radius,
    bool equal_size = true) {
    return build_bvh(
        (int)lines.size(), equal_size, [&lines, &pos, &radius](int eid) {
            auto f = lines[eid];
            return line_bbox(pos[f.x], pos[f.y], radius[f.x], radius[f.y]);
        });
}

/// Build a points BVH.
inline bvh_tree* build_points_bvh(const vector<int>& points,
    const vector<vec3f>& pos, const vector<float>& radius,
    bool equal_size = true) {
    return build_bvh(
        (int)points.size(), equal_size, [&points, &pos, &radius](int eid) {
            auto f = points[eid];
            return point_bbox(pos[f], radius[f]);
        });
}

/// Build a points BVH.
inline bvh_tree* build_points_bvh(const vector<vec3f>& pos,
    const vector<float>& radius, bool equal_size = true) {
    return build_bvh((int)pos.size(), equal_size, [&pos, &radius](int eid) {
        auto r = (radius.empty()) ? 0.00001f : radius[eid];
        return point_bbox(pos[eid], r);
    });
}

/// Recursively recomputes the node bounds for a shape bvh
inline void refit_bvh(
    bvh_tree* bvh, int nodeid, const function<bbox3f(int)>& elem_bbox) {
    // refit
    auto node = &bvh->nodes[nodeid];
    node->bbox = invalid_bbox3f;
    if (node->isleaf) {
        for (auto i = 0; i < node->count; i++) {
            auto idx = bvh->sorted_prim[node->start + i];
            node->bbox += elem_bbox(idx);
        }
    } else {
        for (auto i = 0; i < node->count; i++) {
            auto idx = node->start + i;
            refit_bvh(bvh, idx, elem_bbox);
            node->bbox += bvh->nodes[idx].bbox;
        }
    }
}

/// Refit triangles bvh
inline void refit_triangles_bvh(
    bvh_tree* bvh, const vec3i* triangles, const vec3f* pos) {
    refit_bvh(bvh, 0, [triangles, pos](int eid) {
        auto f = triangles[eid];
        return triangle_bbox(pos[f.x], pos[f.y], pos[f.z]);
    });
}

/// Refit triangles bvh
inline void refit_triangles_bvh(
    bvh_tree* bvh, const vector<vec3i>& triangles, const vector<vec3f>& pos) {
    refit_triangles_bvh(bvh, triangles.data(), pos.data());
}

/// Refit quads bvh
inline void refit_quads_bvh(
    bvh_tree* bvh, const vec4i* quads, const vec3f* pos) {
    refit_bvh(bvh, 0, [quads, pos](int eid) {
        auto f = quads[eid];
        return quad_bbox(pos[f.x], pos[f.y], pos[f.z], pos[f.w]);
    });
}

/// Refit quads bvh
inline void refit_quads_bvh(
    bvh_tree* bvh, const vector<vec4i>& quads, const vector<vec3f>& pos) {
    refit_quads_bvh(bvh, quads.data(), pos.data());
}

/// Refit lines bvh
inline void refit_lines_bvh(
    bvh_tree* bvh, const vec2i* lines, const vec3f* pos, const float* radius) {
    refit_bvh(bvh, 0, [lines, pos, radius](int eid) {
        auto f = lines[eid];
        return line_bbox(pos[f.x], pos[f.y], radius[f.x], radius[f.y]);
    });
}

/// Refit lines bvh
inline void refit_lines_bvh(bvh_tree* bvh, const vector<vec2i>& lines,
    const vector<vec3f>& pos, const vector<float>& radius) {
    refit_lines_bvh(bvh, lines.data(), pos.data(), radius.data());
}

/// Refit points bvh
inline void refit_points_bvh(
    bvh_tree* bvh, const int* points, const vec3f* pos, const float* radius) {
    refit_bvh(bvh, 0, [points, pos, radius](int eid) {
        auto f = points[eid];
        return point_bbox(pos[f], (radius) ? radius[f] : 0);
    });
}

/// Refit points bvh
inline void refit_points_bvh(bvh_tree* bvh, const vector<int>& points,
    const vector<vec3f>& pos, const vector<float>& radius) {
    refit_points_bvh(bvh, points.data(), pos.data(), radius.data());
}
/// Refit points bvh
inline void refit_points_bvh(
    bvh_tree* bvh, const vec3f* pos, const float* radius) {
    refit_bvh(bvh, 0,
        [pos, radius](int eid) { return point_bbox(pos[eid], radius[eid]); });
}

/// Refit lines bvh
inline void refit_points_bvh(
    bvh_tree* bvh, const vector<vec3f>& pos, const vector<float>& radius) {
    refit_points_bvh(bvh, pos.data(), radius.data());
}

/// Intersect ray with a bvh.
inline bool intersect_bvh(const bvh_tree* bvh, const ray3f& ray_,
    bool early_exit, float& ray_t, int& eid,
    const function<bool(int, const ray3f&, float&)>& intersect_elem) {
    // node stack
    int node_stack[64];
    auto node_cur = 0;
    node_stack[node_cur++] = 0;

    // shared variables
    auto hit = false;

    // copy ray to modify it
    auto ray = ray_;

    // prepare ray for fast queries
    auto ray_dinv = vec3f{1, 1, 1} / ray.d;
    auto ray_dsign = vec3i{(ray_dinv.x < 0) ? 1 : 0, (ray_dinv.y < 0) ? 1 : 0,
        (ray_dinv.z < 0) ? 1 : 0};
    auto ray_reverse = array<bool, 4>{
        {(bool)ray_dsign.x, (bool)ray_dsign.y, (bool)ray_dsign.z, false}};

    // walking stack
    while (node_cur) {
        // grab node
        auto node = bvh->nodes[node_stack[--node_cur]];

        // intersect bbox
        if (!intersect_check_bbox(ray, ray_dinv, ray_dsign, node.bbox))
            continue;

        // intersect node, switching based on node type
        // for each type, iterate over the the primitive list
        if (!node.isleaf) {
            // for internal nodes, attempts to proceed along the
            // split axis from smallest to largest nodes
            if (ray_reverse[node.axis]) {
                for (auto i = 0; i < node.count; i++) {
                    auto idx = node.start + i;
                    node_stack[node_cur++] = idx;
                    assert(node_cur < 64);
                }
            } else {
                for (auto i = node.count - 1; i >= 0; i--) {
                    auto idx = node.start + i;
                    node_stack[node_cur++] = idx;
                    assert(node_cur < 64);
                }
            }
        } else {
            for (auto i = 0; i < node.count; i++) {
                auto idx = bvh->sorted_prim[node.start + i];
                if (intersect_elem(idx, ray, ray_t)) {
                    hit = true;
                    ray.tmax = ray_t;
                    eid = idx;
                    if (early_exit) return true;
                }
            }
        }
    }

    return hit;
}

/// Finds the closest element with a bvh.
inline bool overlap_bvh(const bvh_tree* bvh, const vec3f& pos, float max_dist,
    bool early_exit, float& dist, int& eid,
    const function<bool(int, const vec3f&, float, float&)>& overlap_elem) {
    // node stack
    int node_stack[64];
    auto node_cur = 0;
    node_stack[node_cur++] = 0;

    // hit
    auto hit = false;

    // walking stack
    while (node_cur) {
        // grab node
        auto node = bvh->nodes[node_stack[--node_cur]];

        // intersect bbox
        if (!distance_check_bbox(pos, max_dist, node.bbox)) continue;

        // intersect node, switching based on node type
        // for each type, iterate over the the primitive list
        if (!node.isleaf) {
            // internal node
            for (auto idx = node.start; idx < node.start + node.count; idx++) {
                node_stack[node_cur++] = idx;
                assert(node_cur < 64);
            }
        } else {
            for (auto i = 0; i < node.count; i++) {
                auto idx = bvh->sorted_prim[node.start + i];
                if (overlap_elem(idx, pos, max_dist, dist)) {
                    hit = true;
                    max_dist = dist;
                    eid = idx;
                    if (early_exit) return true;
                }
            }
        }
    }

    return hit;
}

/// Intersect a triangle BVH
inline bool intersect_triangles_bvh(const bvh_tree* bvh, const vec3i* triangles,
    const vec3f* pos, const ray3f& ray, bool early_exit, float& ray_t, int& eid,
    vec3f& euv) {
    return intersect_bvh(bvh, ray, early_exit, ray_t, eid,
        [&triangles, &pos, &euv](int eid, const ray3f& ray, float& ray_t) {
            const auto& f = triangles[eid];
            return intersect_triangle(
                ray, pos[f.x], pos[f.y], pos[f.z], ray_t, euv);
        });
}

/// Intersect a triangle BVH
inline bool intersect_triangles_bvh(const bvh_tree* bvh,
    const vector<vec3i>& triangles, const vector<vec3f>& pos, const ray3f& ray,
    bool early_exit, float& ray_t, int& eid, vec3f& euv) {
    return intersect_triangles_bvh(
        bvh, triangles.data(), pos.data(), ray, early_exit, ray_t, eid, euv);
}

/// Intersect a quad BVH
inline bool intersect_quads_bvh(const bvh_tree* bvh, const vec4i* quads,
    const vec3f* pos, const ray3f& ray, bool early_exit, float& ray_t, int& eid,
    vec4f& euv) {
    return intersect_bvh(bvh, ray, early_exit, ray_t, eid,
        [&quads, &pos, &euv](int eid, const ray3f& ray, float& ray_t) {
            const auto& f = quads[eid];
            return intersect_quad(
                ray, pos[f.x], pos[f.y], pos[f.z], pos[f.w], ray_t, euv);
        });
}

/// Intersect a quad BVH
inline bool intersect_quads_bvh(const bvh_tree* bvh, const vector<vec4i>& quads,
    const vector<vec3f>& pos, const ray3f& ray, bool early_exit, float& ray_t,
    int& eid, vec4f& euv) {
    return intersect_quads_bvh(
        bvh, quads.data(), pos.data(), ray, early_exit, ray_t, eid, euv);
}

/// Intersect a line BVH
inline bool intersect_lines_bvh(const bvh_tree* bvh, const vec2i* lines,
    const vec3f* pos, const float* radius, const ray3f& ray, bool early_exit,
    float& ray_t, int& eid, vec2f& euv) {
    return intersect_bvh(bvh, ray, early_exit, ray_t, eid,
        [&lines, &pos, &radius, &euv](int eid, const ray3f& ray, float& ray_t) {
            auto f = lines[eid];
            return intersect_line(
                ray, pos[f.x], pos[f.y], radius[f.x], radius[f.y], ray_t, euv);
        });
}

/// Intersect a line BVH
inline bool intersect_lines_bvh(const bvh_tree* bvh, const vector<vec2i>& lines,
    const vector<vec3f>& pos, const vector<float>& radius, const ray3f& ray,
    bool early_exit, float& ray_t, int& eid, vec2f& euv) {
    return intersect_lines_bvh(bvh, lines.data(), pos.data(), radius.data(),
        ray, early_exit, ray_t, eid, euv);
}

/// Intersect a point BVH
inline bool intersect_points_bvh(const bvh_tree* bvh, const int* points,
    const vec3f* pos, const float* radius, const ray3f& ray, bool early_exit,
    float& ray_t, int& eid) {
    return intersect_bvh(bvh, ray, early_exit, ray_t, eid,
        [&points, &pos, &radius](int eid, const ray3f& ray, float& ray_t) {
            auto f = points[eid];
            return intersect_point(ray, pos[f], radius[f], ray_t);
        });
}

/// Intersect a point BVH
inline bool intersect_points_bvh(const bvh_tree* bvh, const vector<int>& points,
    const vector<vec3f>& pos, const vector<float>& radius, const ray3f& ray,
    bool early_exit, float& ray_t, int& eid) {
    return intersect_points_bvh(bvh, points.data(), pos.data(), radius.data(),
        ray, early_exit, ray_t, eid);
}

/// Intersect a point BVH
inline bool intersect_points_bvh(const bvh_tree* bvh, const vec3f* pos,
    const float* radius, const ray3f& ray, bool early_exit, float& ray_t,
    int& eid) {
    return intersect_bvh(bvh, ray, early_exit, ray_t, eid,
        [&pos, &radius](int eid, const ray3f& ray, float& ray_t) {
            return intersect_point(ray, pos[eid], radius[eid], ray_t);
        });
}

/// Intersect a point BVH
inline bool intersect_points_bvh(const bvh_tree* bvh, const vector<vec3f>& pos,
    const vector<float>& radius, const ray3f& ray, bool early_exit,
    float& ray_t, int& eid) {
    return intersect_points_bvh(
        bvh, pos.data(), radius.data(), ray, early_exit, ray_t, eid);
}

/// Intersect a triangle BVH
inline bool overlap_triangles_bvh(const bvh_tree* bvh, const vec3i* triangles,
    const vec3f* pos, const float* radius, const vec3f& pt, float max_dist,
    bool early_exit, float& dist, int& eid, vec3f& euv) {
    return overlap_bvh(bvh, pt, max_dist, early_exit, dist, eid,
        [&triangles, &pos, &radius, &euv](
            int eid, const vec3f& pt, float max_dist, float& dist) {
            auto f = triangles[eid];
            return overlap_triangle(pt, max_dist, pos[f.x], pos[f.y], pos[f.z],
                (radius) ? radius[f.x] : 0, (radius) ? radius[f.y] : 0,
                (radius) ? radius[f.z] : 0, dist, euv);
        });
}

/// Intersect a triangle BVH
inline bool overlap_triangles_bvh(const bvh_tree* bvh,
    const vector<vec3i>& triangles, const vector<vec3f>& pos,
    const vector<float>& radius, const vec3f& pt, float max_dist,
    bool early_exit, float& dist, int& eid, vec3f& euv) {
    return overlap_triangles_bvh(bvh, triangles.data(), pos.data(),
        radius.data(), pt, max_dist, early_exit, dist, eid, euv);
}

/// Intersect a quad BVH
inline bool overlap_quads_bvh(const bvh_tree* bvh, const vec4i* quads,
    const vec3f* pos, const float* radius, const vec3f& pt, float max_dist,
    bool early_exit, float& dist, int& eid, vec4f& euv) {
    return overlap_bvh(bvh, pt, max_dist, early_exit, dist, eid,
        [&quads, &pos, &radius, &euv](
            int eid, const vec3f& pt, float max_dist, float& dist) {
            auto f = quads[eid];
            return overlap_quad(pt, max_dist, pos[f.x], pos[f.y], pos[f.z],
                pos[f.w], (radius) ? radius[f.x] : 0,
                (radius) ? radius[f.y] : 0, (radius) ? radius[f.z] : 0,
                (radius) ? radius[f.w] : 0, dist, euv);
        });
}

/// Intersect a quad BVH
inline bool overlap_quads_bvh(const bvh_tree* bvh, const vector<vec4i>& quads,
    const vector<vec3f>& pos, const vector<float>& radius, const vec3f& pt,
    float max_dist, bool early_exit, float& dist, int& eid, vec4f& euv) {
    return overlap_quads_bvh(bvh, quads.data(), pos.data(), radius.data(), pt,
        max_dist, early_exit, dist, eid, euv);
}

/// Intersect a line BVH
inline bool overlap_lines_bvh(const bvh_tree* bvh, const vec2i* lines,
    const vec3f* pos, const float* radius, const vec3f& pt, float max_dist,
    bool early_exit, float& dist, int& eid, vec2f& euv) {
    return overlap_bvh(bvh, pt, max_dist, early_exit, dist, eid,
        [&lines, &pos, &radius, &euv](
            int eid, const vec3f& pt, float max_dist, float& dist) {
            auto f = lines[eid];
            return overlap_line(pt, max_dist, pos[f.x], pos[f.y],
                (radius) ? radius[f.x] : 0, (radius) ? radius[f.y] : 0, dist,
                euv);
        });
}

/// Intersect a line BVH
inline bool overlap_lines_bvh(const bvh_tree* bvh, const vector<vec2i>& lines,
    const vector<vec3f>& pos, const vector<float>& radius, const vec3f& pt,
    float max_dist, bool early_exit, float& dist, int& eid, vec2f& euv) {
    return overlap_lines_bvh(bvh, lines.data(), pos.data(), radius.data(), pt,
        max_dist, early_exit, dist, eid, euv);
}

/// Intersect a point BVH
inline bool overlap_points_bvh(const bvh_tree* bvh, const int* points,
    const vec3f* pos, const float* radius, const vec3f& pt, float max_dist,
    bool early_exit, float& dist, int& eid) {
    return overlap_bvh(bvh, pt, max_dist, early_exit, dist, eid,
        [&points, &pos, &radius](
            int eid, const vec3f& pt, float max_dist, float& dist) {
            auto f = points[eid];
            return overlap_point(
                pt, max_dist, pos[f], (radius) ? radius[f] : 0, dist);
        });
}

/// Intersect a point BVH
inline bool overlap_points_bvh(const bvh_tree* bvh, const vector<int>& points,
    const vector<vec3f>& pos, const vector<float>& radius, const vec3f& pt,
    float max_dist, bool early_exit, float& dist, int& eid) {
    return overlap_points_bvh(bvh, points.data(), pos.data(), radius.data(), pt,
        max_dist, early_exit, dist, eid);
}

/// Intersect a point BVH
inline bool overlap_points_bvh(const bvh_tree* bvh, const vec3f* pos,
    const float* radius, const vec3f& pt, float max_dist, bool early_exit,
    float& dist, int& eid) {
    return overlap_bvh(bvh, pt, max_dist, early_exit, dist, eid,
        [&pos, &radius](int eid, const vec3f& pt, float max_dist, float& dist) {
            return overlap_point(pt, max_dist, pos[eid], radius[eid], dist);
        });
}

/// Intersect a point BVH
inline bool overlap_points_bvh(const bvh_tree* bvh, const vector<vec3f>& pos,
    const vector<float>& radius, const vec3f& pt, float max_dist,
    bool early_exit, float& dist, int& eid) {
    return overlap_points_bvh(
        bvh, pos.data(), radius.data(), pt, max_dist, early_exit, dist, eid);
}

/// Finds the overlap between BVH leaf nodes.
template <typename OverlapElem>
void overlap_bvh_elems(const bvh_tree* bvh1, const bvh_tree* bvh2,
    bool skip_duplicates, bool skip_self, vector<vec2i>& overlaps,
    const OverlapElem& overlap_elems) {
    // node stack
    vec2i node_stack[128];
    auto node_cur = 0;
    node_stack[node_cur++] = {0, 0};

    // walking stack
    while (node_cur) {
        // grab node
        auto node_idx = node_stack[--node_cur];
        const auto node1 = bvh1->nodes[node_idx.x];
        const auto node2 = bvh2->nodes[node_idx.y];

        // intersect bbox
        if (!overlap_bbox(node1.bbox, node2.bbox)) continue;

        // check for leaves
        if (node1.isleaf && node2.isleaf) {
            // collide primitives
            for (auto i1 = node1.start; i1 < node1.start + node1.count; i1++) {
                for (auto i2 = node2.start; i2 < node2.start + node2.count;
                     i2++) {
                    auto idx1 = bvh1->sorted_prim[i1];
                    auto idx2 = bvh2->sorted_prim[i2];
                    if (skip_duplicates && idx1 > idx2) continue;
                    if (skip_self && idx1 == idx2) continue;
                    if (overlap_elems(idx1, idx2))
                        overlaps.push_back({idx1, idx2});
                }
            }
        } else {
            // descend
            if (node1.isleaf) {
                for (auto idx2 = node2.start; idx2 < node2.start + node2.count;
                     idx2++) {
                    node_stack[node_cur++] = {node_idx.x, (int)idx2};
                    assert(node_cur < 128);
                }
            } else if (node2.isleaf) {
                for (auto idx1 = node1.start; idx1 < node1.start + node1.count;
                     idx1++) {
                    node_stack[node_cur++] = {(int)idx1, node_idx.y};
                    assert(node_cur < 128);
                }
            } else {
                for (auto idx2 = node2.start; idx2 < node2.start + node2.count;
                     idx2++) {
                    for (auto idx1 = node1.start;
                         idx1 < node1.start + node1.count; idx1++) {
                        node_stack[node_cur++] = {(int)idx1, (int)idx2};
                        assert(node_cur < 128);
                    }
                }
            }
        }
    }
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// SIMPLE SCENE SUPPORT
// -----------------------------------------------------------------------------
namespace ygl {

/// Scene Texture
struct texture {
    /// name
    string name;
    /// path
    string path;
    /// if loaded, ldr image
    image4b ldr;
    /// if loaded, hdr image
    image4f hdr;

    /// get texture width
    int width() const {
        if (ldr) return ldr.width();
        if (hdr) return hdr.width();
        return 0;
    }
    /// get texture height
    int height() const {
        if (ldr) return ldr.height();
        if (hdr) return hdr.height();
        return 0;
    }
};

/// Scene Texture Additional Information
struct texture_info {
    /// texture pointer
    texture* txt = nullptr;
    /// wrap s coordinate
    bool wrap_s = true;
    /// wrap t coordinate
    bool wrap_t = true;
    /// linear interpolation
    bool linear = true;
    /// mipmaping
    bool mipmap = true;
    /// texture strength (occlusion and normal)
    float scale = 1;

    /// check whether the texture if present
    operator bool() const { return (bool)txt; }
};

/// Material type
enum struct material_type {
    /// Microfacet material type (OBJ)
    specular_roughness = 0,
    /// Base and metallic material (metallic-roughness in glTF)
    metallic_roughness = 1,
    /// Diffuse and specular material (specular-glossness in glTF)
    specular_glossiness = 2,
};

/// Scene Material
struct material {
    // whole material data -------------------
    /// material name
    string name;
    /// double-sided rendering
    bool double_sided = false;
    /// material type
    material_type mtype = material_type::specular_roughness;

    // color information ---------------------
    /// emission color
    vec3f ke = {0, 0, 0};
    /// diffuse color / base color
    vec3f kd = {0, 0, 0};
    /// specular color / metallic factor
    vec3f ks = {0, 0, 0};
    /// clear coat reflection
    vec3f kr = {0, 0, 0};
    /// transmission color
    vec3f kt = {0, 0, 0};
    /// roughness
    float rs = 0.0001;
    /// opacity
    float op = 1;

    // textures -------------------------------
    /// emission texture
    texture_info ke_txt = {};
    /// diffuse texture
    texture_info kd_txt = {};
    /// specular texture
    texture_info ks_txt = {};
    /// reflection texture
    texture_info kr_txt = {};
    /// transmission texture
    texture_info kt_txt = {};
    /// roughness texture
    texture_info rs_txt = {};
    /// bump map texture (heighfield)
    texture_info bump_txt = {};
    /// displacement map texture (heighfield)
    texture_info disp_txt = {};
    /// normal texture
    texture_info norm_txt = {};
    /// occlusion texture
    texture_info occ_txt = {};
};

/// Shape data represented as an indexed array.
/// May contain only one of the points/lines/triangles/quads.
struct shape {
    /// shape name
    string name = "";
    /// path (used for saving in glTF)
    string path = "";
    /// shape material
    material* mat = nullptr;

    // shape elements -------------------------
    /// points
    vector<int> points;
    /// lines
    vector<vec2i> lines;
    /// triangles
    vector<vec3i> triangles;
    /// quads
    vector<vec4i> quads;
    /// face-varying indices for position
    vector<vec4i> quads_pos;
    /// face-varying indices for normal
    vector<vec4i> quads_norm;
    /// face-varying indices for texcoord
    vector<vec4i> quads_texcoord;

    // vertex data ----------------------------
    /// per-vertex position (3 float)
    vector<vec3f> pos;
    /// per-vertex normals (3 float)
    vector<vec3f> norm;
    /// per-vertex texcoord (2 float)
    vector<vec2f> texcoord;
    /// per-vertex second texcoord (2 float)
    vector<vec2f> texcoord1;
    /// per-vertex color (4 float)
    vector<vec4f> color;
    /// per-vertex radius (1 float)
    vector<float> radius;
    /// per-vertex tangent space (4 float)
    vector<vec4f> tangsp;

    // computed data --------------------------
    /// element CDF for sampling
    vector<float> elem_cdf;
    /// BVH
    bvh_tree* bvh = nullptr;
    /// bounding box (needs to be updated explicitly)
    bbox3f bbox = invalid_bbox3f;

    // clean
    ~shape() {
        if (bvh) delete bvh;
    }
};

/// Shape instance.
struct instance {
    // name
    string name;
    /// transform frame
    frame3f frame = identity_frame3f;
    /// shape instance
    shape* shp = nullptr;

    // computed data --------------------------
    /// bounding box (needs to be updated explicitly)
    bbox3f bbox = invalid_bbox3f;

    /// instance transform as matrix
    mat4f xform() const { return to_mat4f(frame); }
};

/// Scene Camera
struct camera {
    /// name
    string name;
    /// transform frame
    frame3f frame = identity_frame3f;
    /// ortho cam
    bool ortho = false;
    /// vertical field of view
    float yfov = 2;
    /// aspect ratio
    float aspect = 16.0f / 9.0f;
    /// focus distance
    float focus = 1;
    /// lens aperture
    float aperture = 0;
    /// near plane distance
    float near = 0.01f;
    /// far plane distance
    float far = 10000;
};

/// Envinonment map
struct environment {
    /// name
    string name;
    /// transform frame
    frame3f frame = identity_frame3f;
    /// emission coefficient
    vec3f ke = {0, 0, 0};
    /// emission texture
    texture_info ke_txt = {};
};

/// Light, either an instance or an environment.
/// This is only used internally to avoid looping over all objects every time.
struct light {
    /// instance
    instance* ist = nullptr;
    /// environment
    environment* env = nullptr;
};

/// Scene
struct scene {
    /// shape array
    vector<shape*> shapes;
    /// instance array
    vector<instance*> instances;
    /// material array
    vector<material*> materials;
    /// texture array
    vector<texture*> textures;
    /// camera array
    vector<camera*> cameras;
    /// environment array
    vector<environment*> environments;

    /// light array
    vector<light*> lights;

    // computed data --------------------------
    /// BVH
    bvh_tree* bvh = nullptr;
    /// bounding box (needs to be updated explicitly)
    bbox3f bbox = invalid_bbox3f;

    /// cleanup
    ~scene() {
        for (auto v : shapes)
            if (v) delete v;
        for (auto v : instances)
            if (v) delete v;
        for (auto v : materials)
            if (v) delete v;
        for (auto v : textures)
            if (v) delete v;
        for (auto v : cameras)
            if (v) delete v;
        for (auto v : environments)
            if (v) delete v;
        for (auto light : lights)
            if (light) delete light;
        if (bvh) delete bvh;
    }
};

/// Shape value interpolated using barycentric coordinates
template <typename T>
inline T eval_barycentric(
    const shape* shp, const vector<T>& vals, int eid, const vec4f& euv) {
    if (vals.empty()) return T();
    if (!shp->triangles.empty()) {
        return eval_barycentric_triangle(
            vals, shp->triangles[eid], {euv.x, euv.y, euv.z});
    } else if (!shp->lines.empty()) {
        return eval_barycentric_line(vals, shp->lines[eid], {euv.x, euv.y});
    } else if (!shp->points.empty()) {
        return eval_barycentric_point(vals, shp->points[eid], euv.x);
    } else if (!shp->quads.empty()) {
        return eval_barycentric_quad(vals, shp->quads[eid], euv);
    } else {
        return vals[eid];  // points
    }
}

/// Shape position interpolated using barycentric coordinates
inline vec3f eval_pos(const shape* shp, int eid, const vec4f& euv) {
    return eval_barycentric(shp, shp->pos, eid, euv);
}

/// Shape normal interpolated using barycentric coordinates
inline vec3f eval_norm(const shape* shp, int eid, const vec4f& euv) {
    return normalize(eval_barycentric(shp, shp->norm, eid, euv));
}

/// Shape texcoord interpolated using barycentric coordinates
inline vec2f eval_texcoord(const shape* shp, int eid, const vec4f& euv) {
    return eval_barycentric(shp, shp->texcoord, eid, euv);
}

/// Shape texcoord interpolated using barycentric coordinates
inline vec4f eval_color(const shape* shp, int eid, const vec4f& euv) {
    return eval_barycentric(shp, shp->color, eid, euv);
}

/// Shape tangent space interpolated using barycentric coordinates
inline vec4f eval_tangsp(const shape* shp, int eid, const vec4f& euv) {
    return eval_barycentric(shp, shp->tangsp, eid, euv);
}

/// Instance position interpolated using barycentric coordinates
inline vec3f eval_pos(const instance* ist, int eid, const vec4f& euv) {
    return transform_point(
        ist->frame, eval_barycentric(ist->shp, ist->shp->pos, eid, euv));
}

/// Instance normal interpolated using barycentric coordinates
inline vec3f eval_norm(const instance* ist, int eid, const vec4f& euv) {
    return transform_direction(ist->frame,
        normalize(eval_barycentric(ist->shp, ist->shp->norm, eid, euv)));
}

/// Evaluate a texture
inline vec4f eval_texture(const texture_info& info, const vec2f& texcoord,
    bool srgb = true, const vec4f& def = {1, 1, 1, 1}) {
    if (!info.txt) return def;

    // get texture
    auto txt = info.txt;
    assert(txt->hdr || txt->ldr);

    auto lookup = [&def, &txt, &srgb](int i, int j) {
        if (txt->ldr)
            return (srgb) ? srgb_to_linear(txt->ldr[{i, j}]) :
                            byte_to_float(txt->ldr[{i, j}]);
        else if (txt->hdr)
            return txt->hdr[{i, j}];
        else
            return def;
    };

    // get image width/height
    auto w = txt->width(), h = txt->height();

    // get coordinates normalized for tiling
    auto s = 0.0f, t = 0.0f;
    if (!info.wrap_s) {
        s = clamp(texcoord.x, 0.0f, 1.0f) * w;
    } else {
        s = std::fmod(texcoord.x, 1.0f) * w;
        if (s < 0) s += w;
    }
    if (!info.wrap_t) {
        t = clamp(texcoord.y, 0.0f, 1.0f) * h;
    } else {
        t = std::fmod(texcoord.y, 1.0f) * h;
        if (t < 0) t += h;
    }

    // get image coordinates and residuals
    auto i = clamp((int)s, 0, w - 1), j = clamp((int)t, 0, h - 1);
    auto ii = (i + 1) % w, jj = (j + 1) % h;
    auto u = s - i, v = t - j;

    // nearest lookup
    if (!info.linear) return lookup(i, j);

    // handle interpolation
    return lookup(i, j) * (1 - u) * (1 - v) + lookup(i, jj) * (1 - u) * v +
           lookup(ii, j) * u * (1 - v) + lookup(ii, jj) * u * v;
}

/// Subdivides shape elements. Apply subdivision surface rules if subdivide
/// is true.
inline void subdivide_shape(shape* shp, bool subdiv = false) {
    if (!shp->lines.empty() || !shp->triangles.empty() || !shp->quads.empty()) {
        vector<vec2i> edges;
        vector<vec4i> faces;
        tie(shp->lines, shp->triangles, shp->quads, edges, faces) =
            subdivide_elems(
                shp->lines, shp->triangles, shp->quads, (int)shp->pos.size());
        shp->pos = subdivide_vert(shp->pos, edges, faces);
        shp->norm = subdivide_vert(shp->norm, edges, faces);
        shp->texcoord = subdivide_vert(shp->texcoord, edges, faces);
        shp->color = subdivide_vert(shp->color, edges, faces);
        shp->radius = subdivide_vert(shp->radius, edges, faces);
        if (subdiv && !shp->quads.empty()) {
            auto boundary = get_boundary_edges({}, {}, shp->quads);
            shp->pos =
                subdivide_catmullclark(shp->quads, shp->pos, boundary, {});
            shp->norm =
                subdivide_catmullclark(shp->quads, shp->norm, boundary, {});
            shp->texcoord =
                subdivide_catmullclark(shp->quads, shp->texcoord, boundary, {});
            shp->color =
                subdivide_catmullclark(shp->quads, shp->color, boundary, {});
            shp->radius =
                subdivide_catmullclark(shp->quads, shp->radius, boundary, {});
            shp->norm = compute_normals({}, {}, shp->quads, shp->pos);
        }
    } else if (!shp->quads_pos.empty()) {
        vector<vec2i> _lines;
        vector<vec3i> _triangles;
        vector<vec2i> edges;
        vector<vec4i> faces;
        tie(_lines, _triangles, shp->quads_pos, edges, faces) =
            subdivide_elems({}, {}, shp->quads_pos, shp->pos.size());
        shp->pos = subdivide_vert(shp->pos, edges, faces);
        tie(_lines, _triangles, shp->quads_norm, edges, faces) =
            subdivide_elems({}, {}, shp->quads_norm, shp->norm.size());
        shp->norm = subdivide_vert(shp->norm, edges, faces);
        tie(_lines, _triangles, shp->quads_texcoord, edges, faces) =
            subdivide_elems({}, {}, shp->quads_texcoord, shp->texcoord.size());
        shp->texcoord = subdivide_vert(shp->texcoord, edges, faces);
        if (subdiv) {
            shp->pos = subdivide_catmullclark(shp->quads_pos, shp->pos,
                get_boundary_edges({}, {}, shp->quads_pos), {});
            shp->norm = subdivide_catmullclark(shp->quads_norm, shp->norm,
                get_boundary_edges({}, {}, shp->quads_norm), {});
            shp->texcoord =
                subdivide_catmullclark(shp->quads_texcoord, shp->texcoord, {},
                    get_boundary_verts({}, {}, shp->quads_texcoord));
        }
    }
}

/// Facet a shape. Supports only non-face0varying shapes
inline void facet_shape(shape* shp) {
    if (!shp->lines.empty() || !shp->triangles.empty() || !shp->quads.empty()) {
        vector<int> verts;
        tie(shp->lines, shp->triangles, shp->quads, verts) =
            facet_elems(shp->lines, shp->triangles, shp->quads);
        shp->pos = facet_vert(shp->pos, verts);
        shp->norm = facet_vert(shp->norm, verts);
        shp->texcoord = facet_vert(shp->texcoord, verts);
        shp->color = facet_vert(shp->color, verts);
        shp->radius = facet_vert(shp->radius, verts);
    }
}

/// Tesselate a shape into basic primitives
inline void tesselate_shape(shape* shp) {
    if (shp->quads_pos.empty()) return;
    std::tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
        convert_face_varying(shp->quads_pos, shp->quads_norm,
            shp->quads_texcoord, shp->pos, shp->norm, shp->texcoord);
}

/// Tesselate scene shapes and update pointers
inline void tesselate_shapes(scene* scn) {
    for (auto shp : scn->shapes) tesselate_shape(shp);
}

/// Loading options
struct load_options {
    /// Whether to load textures
    bool load_textures = true;
    /// Skip missing files without giving and error
    bool skip_missing = true;
    /// Whether to flip the v coordinate in OBJ
    bool obj_flip_texcoord = true;
    /// Duplicate vertices if smoothing off in OBJ
    bool obj_facet_non_smooth = false;
    /// Whether to flip tr in OBJ
    bool obj_flip_tr = true;
    /// whether to preserve quads
    bool preserve_quads = false;
    /// whether to preserve face-varying faces
    bool preserve_facevarying = false;
};

/// Loads a scene. For now OBJ or glTF are supported.
/// Throws an exception if an error occurs.
scene* load_scene(const string& filename, const load_options& opts = {});

/// Save options
struct save_options {
    /// Whether to save textures
    bool save_textures = true;
    /// Skip missing files without giving and error
    bool skip_missing = true;
    /// Whether to flip the v coordinate in OBJ
    bool obj_flip_texcoord = true;
    /// Whether to flip tr in OBJ
    bool obj_flip_tr = true;
    /// Whether to use separate buffers in gltf
    bool gltf_separate_buffers = false;
};

/// Saves a scene. For now OBJ and glTF are supported.
/// Throws an exception if an error occurs.
void save_scene(
    const string& filename, const scene* scn, const save_options& opts);

/// Add elements options
struct add_elements_options {
    /// Add missing normal
    bool smooth_normals = true;
    /// Add missing radius for points and lines (<=0 for no adding)
    float pointline_radius = 0;
    /// Add missing trangent space
    bool tangent_space = true;
    /// texture data
    bool texture_data = true;
    /// Add instances
    bool shape_instances = true;
    /// Add default camera
    bool default_camera = true;
    /// Add an empty default environment
    bool default_environment = false;
    /// Add default names
    bool default_names = true;
    /// Add default paths
    bool default_paths = true;

    /// initialize to no element
    static add_elements_options none() {
        auto opts = add_elements_options();
        memset(&opts, 0, sizeof(opts));
        return opts;
    }
};

/// Add elements
void add_elements(scene* scn, const add_elements_options& opts = {});

/// Merge scene into one another. Note that the objects are _moved_ from
/// merge_from to merged_into, so merge_from will be empty after this function.
void merge_into(scene* merge_into, scene* merge_from);

/// Computes a shape bounding box (quick computation that ignores radius)
inline void update_bounds(shape* shp) {
    shp->bbox = invalid_bbox3f;
    for (auto p : shp->pos) shp->bbox += vec3f(p);
}

/// Updates the instance bounding box
inline void update_bounds(instance* ist, bool do_shape = true) {
    if (do_shape) update_bounds(ist->shp);
    ist->bbox = transform_bbox(ist->frame, ist->shp->bbox);
}

/// Updates the scene and scene's instances bounding boxes
inline void update_bounds(scene* scn, bool do_shapes = true) {
    if (do_shapes) {
        for (auto shp : scn->shapes) update_bounds(shp);
    }
    scn->bbox = invalid_bbox3f;
    if (!scn->instances.empty()) {
        for (auto ist : scn->instances) {
            update_bounds(ist, false);
            scn->bbox += ist->bbox;
        }
    } else {
        for (auto shp : scn->shapes) { scn->bbox += shp->bbox; }
    }
}

/// Flatten scene instances into separate meshes.
inline void flatten_instances(scene* scn) {
    if (scn->instances.empty()) return;
    auto shapes = scn->shapes;
    scn->shapes.clear();
    auto instances = scn->instances;
    scn->instances.clear();
    for (auto ist : instances) {
        if (!ist->shp) continue;
        auto xf = ist->xform();
        auto nshp = new shape(*ist->shp);
        for (auto& p : nshp->pos) p = transform_point(xf, p);
        for (auto& n : nshp->norm) n = transform_direction(xf, n);
        scn->shapes.push_back(nshp);
    }
    for (auto e : shapes) delete e;
    for (auto e : instances) delete e;
}

/// Initialize the lights
void update_lights(scene* scn, bool point_only);

/// Print scene information (call update bounds bes before)
void print_info(const scene* scn);

/// Build a shape BVH
inline void build_bvh(shape* shp, bool equalsize = true) {
    if (!shp->points.empty()) {
        shp->bvh =
            build_points_bvh(shp->points, shp->pos, shp->radius, equalsize);
    } else if (!shp->lines.empty()) {
        shp->bvh =
            build_lines_bvh(shp->lines, shp->pos, shp->radius, equalsize);
    } else if (!shp->triangles.empty()) {
        shp->bvh = build_triangles_bvh(shp->triangles, shp->pos, equalsize);
    } else if (!shp->quads.empty()) {
        shp->bvh = build_quads_bvh(shp->quads, shp->pos, equalsize);
    } else {
        shp->bvh = build_points_bvh(shp->pos, shp->radius, equalsize);
    }
    shp->bbox = shp->bvh->nodes[0].bbox;
}

/// Build a scene BVH
inline void build_bvh(
    scene* scn, bool equalsize = true, bool do_shapes = true) {
    // do shapes
    if (do_shapes) {
        for (auto shp : scn->shapes) build_bvh(shp, equalsize);
    }

    // update instance bbox
    for (auto ist : scn->instances)
        ist->bbox = transform_bbox(ist->frame, ist->shp->bbox);

    // tree bvh
    scn->bvh = build_bvh((int)scn->instances.size(), equalsize,
        [scn](int eid) { return scn->instances[eid]->bbox; });
}

/// Refits a scene BVH
inline void refit_bvh(shape* shp) {
    if (!shp->points.empty()) {
        refit_points_bvh(shp->bvh, shp->points, shp->pos, shp->radius);
    } else if (!shp->lines.empty()) {
        refit_lines_bvh(shp->bvh, shp->lines, shp->pos, shp->radius);
    } else if (!shp->triangles.empty()) {
        refit_triangles_bvh(shp->bvh, shp->triangles, shp->pos);
    } else if (!shp->quads.empty()) {
        refit_quads_bvh(shp->bvh, shp->quads, shp->pos);
    } else {
        refit_points_bvh(shp->bvh, shp->pos, shp->radius);
    }
    shp->bbox = shp->bvh->nodes[0].bbox;
}

/// Refits a scene BVH
inline void refit_bvh(scene* scn, bool do_shapes = true) {
    if (do_shapes) {
        for (auto shp : scn->shapes) refit_bvh(shp);
    }

    // update instance bbox
    for (auto ist : scn->instances)
        ist->bbox = transform_bbox(ist->frame, ist->shp->bbox);

    // recompute bvh bounds
    refit_bvh(
        scn->bvh, 0, [scn](int eid) { return scn->instances[eid]->bbox; });
}

/// Intersect the shape with a ray. Find any interstion if early_exit,
/// otherwise find first intersection.
///
/// - Parameters:
///     - scn: scene to intersect
///     - ray: ray to be intersected
///     - early_exit: whether to stop at the first found hit
///     - ray_t: ray distance at intersection
///     - eid: shape element index
///     - euv: element barycentric coordinates
/// - Returns:
///     - whether it intersected
inline bool intersect_ray(const shape* shp, const ray3f& ray, bool early_exit,
    float& ray_t, int& eid, vec4f& euv) {
    // switch over shape type
    if (!shp->triangles.empty()) {
        if (intersect_triangles_bvh(shp->bvh, shp->triangles, shp->pos, ray,
                early_exit, ray_t, eid, (vec3f&)euv)) {
            euv = {euv.x, euv.y, euv.z, 0};
            return true;
        }
    } else if (!shp->quads.empty()) {
        if (intersect_quads_bvh(shp->bvh, shp->quads, shp->pos, ray, early_exit,
                ray_t, eid, euv)) {
            return true;
        }
    } else if (!shp->lines.empty()) {
        if (intersect_lines_bvh(shp->bvh, shp->lines, shp->pos, shp->radius,
                ray, early_exit, ray_t, eid, (vec2f&)euv)) {
            euv = {euv.x, euv.y, 0, 0};
            return true;
        }
    } else if (!shp->points.empty()) {
        if (intersect_points_bvh(shp->bvh, shp->points, shp->pos, shp->radius,
                ray, early_exit, ray_t, eid)) {
            euv = {1, 0, 0, 0};
            return true;
        }
    } else {
        if (intersect_points_bvh(
                shp->bvh, shp->pos, shp->radius, ray, early_exit, ray_t, eid)) {
            euv = {1, 0, 0, 0};
            return true;
        }
    }

    return false;
}

/// Intersect the instance with a ray. Find any interstion if early_exit,
/// otherwise find first intersection.
///
/// - Parameters:
///     - scn: scene to intersect
///     - ray: ray to be intersected
///     - early_exit: whether to stop at the first found hit
///     - ray_t: ray distance at intersection
///     - eid: shape element index
///     - euv: element barycentric coordinates
/// - Returns:
///     - whether it intersected
inline bool intersect_ray(const instance* ist, const ray3f& ray,
    bool early_exit, float& ray_t, int& eid, vec4f& euv) {
    return intersect_ray(ist->shp, transform_ray_inverse(ist->frame, ray),
        early_exit, ray_t, eid, euv);
}

/// Intersect the scene with a ray. Find any interstion if early_exit,
/// otherwise find first intersection.
///
/// - Parameters:
///     - scn: scene to intersect
///     - ray: ray to be intersected
///     - early_exit: whether to stop at the first found hit
///     - ray_t: ray distance at intersection
///     - iid: instance index
///     - eid: shape element index
///     - euv: element barycentric coordinates
/// - Returns:
///     - whether it intersected
inline bool intersect_ray(const scene* scn, const ray3f& ray, bool early_exit,
    float& ray_t, int& iid, int& eid, vec4f& euv) {
    return intersect_bvh(scn->bvh, ray, early_exit, ray_t, iid,
        [&eid, &euv, early_exit, scn](int iid, const ray3f& ray, float& ray_t) {
            return intersect_ray(
                scn->instances[iid], ray, early_exit, ray_t, eid, euv);
        });
}

/// Surface point.
struct intersection_point {
    /// distance of the hit along the ray or from the point
    float dist = 0;
    /// instance index
    int iid = -1;
    /// shape element index
    int eid = -1;
    /// shape barycentric coordinates
    vec4f euv = zero4f;

    /// check if intersection is valid
    operator bool() const { return eid >= 0; }
};

/// Intersect the scene with a ray. Find any interstion if early_exit,
/// otherwise find first intersection.
///
/// - Parameters:
///     - scn: scene to intersect
///     - ray: ray to be intersected
///     - early_exit: whether to stop at the first found hit
/// - Returns:
///     - intersection record
inline intersection_point intersect_ray(
    const scene* scn, const ray3f& ray, bool early_exit) {
    auto isec = intersection_point();
    if (!intersect_ray(
            scn, ray, early_exit, isec.dist, isec.iid, isec.eid, isec.euv))
        return {};
    return isec;
}

/// Finds the closest element that overlaps a point within a given distance.
///
/// - Parameters:
///     - scn: scene to intersect
///     - pos: point position
///     - max_dist: maximu valid distance
///     - early_exit: whether to stop at the first found hit
///     - dist: distance at intersection
///     - eid: shape element index
///     - euv: element barycentric coordinates
/// - Returns:
///     - whether it intersected
inline bool overlap_point(const shape* shp, const vec3f& pos, float max_dist,
    bool early_exit, float& dist, int& eid, vec4f& euv) {
    // switch over shape type
    if (!shp->triangles.empty()) {
        if (overlap_triangles_bvh(shp->bvh, shp->triangles, shp->pos,
                shp->radius, pos, max_dist, early_exit, dist, eid,
                (vec3f&)euv)) {
            euv = {euv.x, euv.y, euv.z, 0};
            return true;
        }
    } else if (!shp->quads.empty()) {
        if (overlap_quads_bvh(shp->bvh, shp->quads, shp->pos, shp->radius, pos,
                max_dist, early_exit, dist, eid, euv)) {
            return true;
        }
    } else if (!shp->lines.empty()) {
        if (overlap_lines_bvh(shp->bvh, shp->lines, shp->pos, shp->radius, pos,
                max_dist, early_exit, dist, eid, (vec2f&)euv)) {
            euv = {euv.x, euv.y, 0, 0};
            return true;
        }
    } else if (!shp->points.empty()) {
        if (overlap_points_bvh(shp->bvh, shp->points, shp->pos, shp->radius,
                pos, max_dist, early_exit, dist, eid)) {
            euv = {1, 0, 0, 0};
            return true;
        }
    } else {
        if (overlap_points_bvh(shp->bvh, shp->pos, shp->radius, pos, max_dist,
                early_exit, dist, eid)) {
            euv = {1, 0, 0, 0};
        }
        return true;
    }

    return false;
}

/// Finds the closest element that overlaps a point within a given distance.
///
/// - Parameters:
///     - scn: scene to intersect
///     - pos: point position
///     - max_dist: maximu valid distance
///     - early_exit: whether to stop at the first found hit
///     - dist: distance at intersection
///     - eid: shape element index
///     - euv: element barycentric coordinates
/// - Returns:
///     - whether it intersected
inline bool overlap_point(const instance* ist, const vec3f& pos, float max_dist,
    bool early_exit, float& dist, int& eid, vec4f& euv) {
    return overlap_point(ist->shp, transform_point_inverse(ist->frame, pos),
        max_dist, early_exit, dist, eid, euv);
}

/// Finds the closest element that overlaps a point within a given distance.
///
/// - Parameters:
///     - scn: scene to intersect
///     - pos: point position
///     - max_dist: maximu valid distance
///     - early_exit: whether to stop at the first found hit
///     - dist: distance at intersection
///     - iid: instance index
///     - eid: shape element index
///     - euv: element barycentric coordinates
/// - Returns:
///     - whether it intersected
inline bool overlap_point(const scene* scn, const vec3f& pos, float max_dist,
    bool early_exit, float& dist, int& iid, int& eid, vec4f& euv) {
    return overlap_bvh(scn->bvh, pos, max_dist, early_exit, dist, iid,
        [&eid, &euv, early_exit, scn](
            int iid, const vec3f& pos, float max_dist, float& dist) {
            return overlap_point(
                scn->instances[iid], pos, max_dist, early_exit, dist, eid, euv);
        });
}

/// Find the list of overlaps between instance bounds.
inline void overlap_instance_bounds(const scene* scn1, const scene* scn2,
    bool skip_duplicates, bool skip_self, vector<vec2i>& overlaps) {
    overlaps.clear();
    overlap_bvh_elems(scn1->bvh, scn2->bvh, skip_duplicates, skip_self,
        overlaps, [scn1, scn2](int i1, int i2) {
            return overlap_bbox(
                scn1->instances[i1]->bbox, scn2->instances[i2]->bbox);
        });
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// EXAMPLE SCENES
// -----------------------------------------------------------------------------
namespace ygl {

/// Makes the Cornell Box scene
scene* make_cornell_box_scene();

/// Test scene enumeration.
enum struct test_scene_type {
    cornell_box,     // Cornell box
    textures,        // Scene containing all generated textures
    shapes,          // Scene contaning a few meshes
    plane_al,        // Simple scene with a plane and area lights
    nothing_el,      // Simple scene with no objects and an environment map
    basic_pl,        // Simple scene with no textures and point lights
    simple_pl,       // Textured scene with point lights
    simple_al,       // Textured scene with area lights
    simple_el,       // Textured scene with env lights
    transparent_al,  // Scene to test trasparency
    points_al,       // Point primitive
    lines_al,        // Line primitives
    subdiv_al,       // Scene to test subdivions surfaces
    plastics_al,     // Material ball with plastics - area light
    plastics_el,     // Material ball with plastics - env light
    metals_al,       // Material ball with metals - area light
    metals_el,       // Material ball with metals - env light
    tesselation_pl,  // Scene to show different tesselation
    textureuv_pl,    // Scene to show texture uvs
    normalmap_pl,    // Scene to show normal mapping
    instances_pl,    // Scene with small number of instances
    instancel_pl,    // Scene with large number of instances
};

/// Names for enumeration
inline const vector<pair<string, test_scene_type>>& test_scene_names() {
    static auto names = vector<pair<string, test_scene_type>>{
        {"cornell_box", test_scene_type::cornell_box},
        {"textures", test_scene_type::textures},
        {"shapes", test_scene_type::shapes},
        {"plane_al", test_scene_type::plane_al},
        {"nothing_el", test_scene_type::nothing_el},
        {"basic_pl", test_scene_type::basic_pl},
        {"simple_pl", test_scene_type::simple_pl},
        {"simple_al", test_scene_type::simple_al},
        {"simple_el", test_scene_type::simple_el},
        {"transparent_al", test_scene_type::transparent_al},
        {"points_al", test_scene_type::points_al},
        {"lines_al", test_scene_type::lines_al},
        {"subdiv_al", test_scene_type::subdiv_al},
        {"plastics_al", test_scene_type::plastics_al},
        {"plastics_el", test_scene_type::plastics_el},
        {"metals_al", test_scene_type::metals_al},
        {"metals_el", test_scene_type::metals_el},
        {"tesselation_pl", test_scene_type::tesselation_pl},
        {"textureuv_pl", test_scene_type::textureuv_pl},
        {"normalmap_pl", test_scene_type::normalmap_pl},
        {"instances_pl", test_scene_type::instances_pl},
        {"instancel_pl", test_scene_type::instancel_pl},
    };
    return names;
}

/// Makes a test scene
scene* make_test_scene(test_scene_type stype);

}  // namespace ygl

// -----------------------------------------------------------------------------
// PATH TRACING
// -----------------------------------------------------------------------------
namespace ygl {

// convenient typedef for bytes
using byte = unsigned char;

/// Type of rendering algorithm (shader)
enum struct trace_shader_type {
    /// pathtrace
    pathtrace = 0,
    /// eye hight for quick previews
    eyelight,
    /// direct illumination
    direct,
    /// pathtrace without MIS (usedful ony for debugging)
    pathtrace_nomis,
    /// debug normal
    debug_normal,
    /// debug albedo
    debug_albedo,
    /// debug texcoord
    debug_texcoord,
};

/// Names for enumeration
inline const vector<pair<string, trace_shader_type>>& trace_shader_names() {
    static auto names = vector<pair<string, trace_shader_type>>{
        {"path", trace_shader_type::pathtrace},
        {"eye", trace_shader_type::eyelight},
        {"direct", trace_shader_type::direct},
        {"path_nomis", trace_shader_type::pathtrace_nomis},
        {"normal", trace_shader_type::debug_normal},
        {"albedo", trace_shader_type::debug_albedo},
        {"texcoord", trace_shader_type::debug_texcoord},
    };
    return names;
}

/// Random number generator type
enum struct trace_rng_type {
    /// uniform random numbers
    uniform = 0,
    /// stratified random numbers
    stratified,
};

/// Names for enumeration
inline const vector<pair<string, trace_rng_type>>& trace_rng_names() {
    static auto names = vector<pair<string, trace_rng_type>>{
        {"uniform", trace_rng_type::uniform},
        {"stratified", trace_rng_type::stratified}};
    return names;
}

/// Filter type
enum struct trace_filter_type {
    /// box filter
    box = 1,
    /// hat filter
    triangle = 2,
    /// cubic spline
    cubic = 3,
    /// Catmull-Rom spline
    catmull_rom = 4,
    /// Mitchell-Netrevalli
    mitchell = 5
};

/// Names for enumeration
inline const vector<pair<string, trace_filter_type>>& trace_filter_names() {
    static auto names =
        vector<pair<string, trace_filter_type>>{{"box", trace_filter_type::box},
            {"triangle", trace_filter_type::triangle},
            {"cubic", trace_filter_type::cubic},
            {"catmull-rom", trace_filter_type::catmull_rom},
            {"mitchell", trace_filter_type::mitchell}};
    return names;
}

/// Rendering params
struct trace_params {
    /// camera id
    int camera_id = 0;
    /// width
    int width = 0;
    /// height
    int height = 0;
    /// number of samples
    int nsamples = 256;
    /// sampler type
    trace_shader_type stype = trace_shader_type::pathtrace;
    /// wheter to test transmission in shadows
    bool shadow_notransmission = false;
    /// random number generation type
    trace_rng_type rtype = trace_rng_type::stratified;
    /// filter type
    trace_filter_type ftype = trace_filter_type::box;
    /// ambient lighting
    vec3f amb = {0, 0, 0};
    /// view environment map
    bool envmap_invisible = false;
    /// minimum ray depth
    int min_depth = 3;
    /// maximum ray depth
    int max_depth = 8;
    /// final pixel clamping
    float pixel_clamp = 10;
    /// ray intersection epsilon
    float ray_eps = 1e-4f;
    /// parallel execution
    bool parallel = true;
    /// seed for the random number generators
    uint32_t seed = 0;
    /// block size for parallel batches (probably leave it as is)
    int block_size = 32;
};

/// Make image blocks
inline vector<pair<vec2i, vec2i>> trace_blocks(const trace_params& params) {
    vector<pair<vec2i, vec2i>> blocks;
    for (int j = 0; j < params.height; j += params.block_size) {
        for (int i = 0; i < params.width; i += params.block_size) {
            blocks.push_back(
                {{i, j}, {min(i + params.block_size, params.width),
                             min(j + params.block_size, params.height)}});
        }
    }
    return blocks;
}

/// Make a 2D array of random number generators for parallelization
inline vector<rng_pcg32> trace_rngs(const trace_params& params) {
    auto rngs = vector<rng_pcg32>(params.width * params.height);
    for (auto j : range(params.height)) {
        for (auto i : range(params.width)) {
            rngs[j * params.width + i] =
                init_rng(params.seed, (j * params.width + i) * 2 + 1);
        }
    }
    return rngs;
}

/// Renders a block of samples
///
/// Notes: It is safe to call the function in parallel on different blocks.
/// But two threads should not access the same pixels at the same time. If
/// the same block is rendered with different samples, samples have to be
/// sequential.
///
/// - Parameters:
///     - scn: trace scene
///     - img: pixel data in RGBA format (width/height in params)
///     - block: range of pixels to render
///     - samples_min, samples_max: range of samples to render
///     - params: trace params
void trace_block(const scene* scn, image4f& img, const vec2i& block_min,
    const vec2i& block_max, int samples_min, int samples_max,
    vector<rng_pcg32>& rngs, const trace_params& params);

/// Trace the next samples in [samples_min, samples_max) range.
/// Samples have to be traced consecutively.
void trace_samples(const scene* scn, image4f& img, int samples_min,
    int samples_max, vector<rng_pcg32>& rngs, const trace_params& params);

/// Renders a filtered block of samples
///
/// Notes: It is safe to call the function in parallel on different blocks.
/// But two threads should not access the same pixels at the same time. If
/// the same block is rendered with different samples, samples have to be
/// sequential.
///
/// - Parameters:
///     - scn: trace scene
///     - img: pixel data in RGBA format (width/height in params)
///     - acc: accumulation buffer in RGBA format (width/height in params)
///     - weight: weight buffer in float format (width/height in params)
///     - block: range of pixels to render
///     - samples_min, samples_max: range of samples to render
///     - image_mutex: mutex for locking
///     - params: trace params
void trace_block_filtered(const scene* scn, image4f& img, image4f& acc,
    image4f& weight, const vec2i& block_min, const vec2i& block_max,
    int samples_min, int samples_max, vector<rng_pcg32>& rngs,
    std::mutex& image_mutex, const trace_params& params);

/// Trace the next samples in [samples_min, samples_max) range.
/// Samples have to be traced consecutively.
void trace_filtered_samples(const scene* scn, image4f& img, image4f& acc,
    image4f& weight, int samples_min, int samples_max, vector<rng_pcg32>& rngs,
    const trace_params& params);

/// Trace the whole image
inline image4f trace_image(const scene* scn, const trace_params& params) {
    auto img = image4f(params.width, params.height);
    auto rngs = trace_rngs(params);
    if (params.ftype == trace_filter_type::box) {
        trace_samples(scn, img, 0, params.nsamples, rngs, params);
    } else {
        auto acc = image4f(params.width, params.height);
        auto weight = image4f(params.width, params.height);
        trace_filtered_samples(
            scn, img, acc, weight, 0, params.nsamples, rngs, params);
    }
    return img;
}

// forward declaration
struct thread_pool;

/// Starts an anyncrhounous renderer with a maximum of 256 samples.
void trace_async_start(const scene* scn, image4f& img, vector<rng_pcg32>& rngs,
    const trace_params& params, thread_pool* pool,
    const function<void(int)>& callback);

/// Stop the asynchronous renderer.
void trace_async_stop(thread_pool* pool);

}  // namespace ygl

// -----------------------------------------------------------------------------
// WAVEFRONT OBJ SUPPORT
// -----------------------------------------------------------------------------
namespace ygl {

/// Face vertex
struct obj_vertex {
    /// position
    int pos;
    /// texcoord
    int texcoord;
    /// normal
    int norm;
    /// color [extension]
    int color;
    /// radius [extension]
    int radius;

    /// Constructor (copies members initializing missing ones to -1)
    obj_vertex(int pos = -1, int texcoord = -1, int norm = -1, int color = -1,
        int radius = -1)
        : pos(pos)
        , texcoord(texcoord)
        , norm(norm)
        , color(color)
        , radius(radius) {}
};

// Comparison for unordred_map
inline bool operator==(const obj_vertex& a, const obj_vertex& b) {
    return a.pos == b.pos && a.texcoord == b.texcoord && a.norm == b.norm &&
           a.color == b.color && a.radius == b.radius;
}

/// element type
enum struct obj_element_type : uint16_t {
    /// lists of points
    point = 1,
    /// polylines
    line = 2,
    /// polygon faces
    face = 3,
    /// tetrahedrons
    tetra = 4,
};

/// Element vertex indices
struct obj_element {
    /// starting vertex index
    uint32_t start;
    /// element type
    obj_element_type type;
    /// number of vertices
    uint16_t size;
};

/// Element group
struct obj_group {
    // group data ---------------------------
    /// material name
    string matname;
    /// group name
    string groupname;
    /// smoothing
    bool smoothing = true;

    // element data -------------------------
    /// element vertices
    vector<obj_vertex> verts;
    /// element faces
    vector<obj_element> elems;
};

/// Obj object
struct obj_object {
    // object data --------------------------
    /// object name
    string name;

    // element data -------------------------
    /// element groups
    vector<obj_group> groups;
};

/// Texture information for OBJ
struct obj_texture_info {
    /// the texture path
    string path = "";
    /// whether to clamp tp th edge
    bool clamp = false;
    /// the scale for bump and displacement
    float scale = 1;
    /// the rest of the unknown properties
    unordered_map<string, vector<string>> unknown_props;
};

// comparison for texture info
inline bool operator==(const obj_texture_info& a, const obj_texture_info& b) {
    if (a.path.empty() && b.path.empty()) return true;
    if (a.path != b.path) return false;
    return a.clamp == b.clamp && a.scale == b.scale &&
           a.unknown_props == b.unknown_props;
}

/// OBJ texture. Texture data is loaded only if desired.
struct obj_texture {
    // whole texture data ------------------
    /// texture path
    string path;
    /// Width
    int width = 0;
    /// Height
    int height = 0;
    /// Number of Channels
    int ncomp = 0;
    /// Buffer data for 8-bit images
    vector<uint8_t> datab;
    /// Buffer data for float images
    vector<float> dataf;
};

/// OBJ material
struct obj_material {
    // whole material data ------------------
    /// material name
    string name;
    /// MTL illum mode
    int illum = 0;

    // color information --------------------
    /// emission color
    vec3f ke = {0, 0, 0};
    /// ambient color
    vec3f ka = {0, 0, 0};
    /// diffuse color
    vec3f kd = {0, 0, 0};
    /// specular color
    vec3f ks = {0, 0, 0};
    /// reflection color
    vec3f kr = {0, 0, 0};
    /// transmision color
    vec3f kt = {0, 0, 0};
    /// phong exponent for ks
    float ns = 1;
    /// index of refraction
    float ior = 1;
    /// opacity
    float op = 1;

    // texture names for the above properties
    /// emission texture
    obj_texture_info ke_txt;
    /// ambient texture
    obj_texture_info ka_txt;
    /// diffuse texture
    obj_texture_info kd_txt;
    /// specular texture
    obj_texture_info ks_txt;
    /// reflection texture
    obj_texture_info kr_txt;
    /// transmission texture
    obj_texture_info kt_txt;
    /// specular exponent texture
    obj_texture_info ns_txt;
    /// opacity texture
    obj_texture_info op_txt;
    /// index of refraction
    obj_texture_info ior_txt;
    /// bump map texture (heighfield)
    obj_texture_info bump_txt;
    /// displacement map texture (heighfield)
    obj_texture_info disp_txt;
    /// normal map texture
    obj_texture_info norm_txt;

    // unknown properties ---------------------
    /// unknown string props
    unordered_map<string, vector<string>> unknown_props;
};

/// Camera [extension]
struct obj_camera {
    /// camera name
    string name;
    /// transform frame (affine matrix)
    frame3f frame = identity_frame3f;
    /// orthografic camera
    bool ortho = false;
    /// vertical field of view
    float yfov = 2 * atan(0.5f);
    /// aspect ratio
    float aspect = 16.0f / 9.0f;
    /// lens aperture
    float aperture = 0;
    /// focus distance
    float focus = 1;
};

/// Environment [extension]
struct obj_environment {
    /// environment name
    string name;
    /// transform frame (affine matrix)
    frame3f frame = identity_frame3f;
    /// material name
    string matname;
};

/// Instance [extension]
struct obj_instance {
    /// instance name
    string name;
    /// transform frame (affine matrix)
    frame3f frame = identity_frame3f;
    /// object name
    string objname;
};

/// OBJ asset
struct obj_scene {
    // vertex data -------------------------
    /// vertex positions
    vector<vec3f> pos;
    /// vertex normals
    vector<vec3f> norm;
    /// vertex texcoord
    vector<vec2f> texcoord;
    /// vertex color [extension]
    vector<vec4f> color;
    /// vertex radius [extension]
    vector<float> radius;

    // scene objects -----------------------
    /// objects
    vector<obj_object*> objects;
    /// materials
    vector<obj_material*> materials;
    /// textures
    vector<obj_texture*> textures;
    /// cameras [extension]
    vector<obj_camera*> cameras;
    /// env maps [extension]
    vector<obj_environment*> environments;
    /// instances [extension]
    vector<obj_instance*> instances;

    /// cleanup
    ~obj_scene() {
        for (auto v : objects)
            if (v) delete v;
        for (auto v : materials)
            if (v) delete v;
        for (auto v : textures)
            if (v) delete v;
        for (auto v : cameras)
            if (v) delete v;
        for (auto v : environments)
            if (v) delete v;
        for (auto v : instances)
            if (v) delete v;
    }
};

/// Load OBJ
///
/// - Parameters:
///     - filename: filename
///     - load_texture: whether to load textures
///     - skip_missing: whether to skip missing files
///     - flip_texcoord: whether to flip the v coordinate
///     - flip_tr: whether to flip the Tr value
/// - Return:
///     - obj (nullptr on error)
obj_scene* load_obj(const string& filename, bool load_textures = false,
    bool skip_missing = false, bool flip_texcoord = true, bool flip_tr = true);

/// Save OBJ
///
/// - Parameters:
///     - filename: filename
///     - model: obj data to save
///     - save_textures: whether to save textures
///     - skip_missing: whether to skip missing files
///     - flip_texcoord: whether to flip the v coordinate
///     - flip_tr: whether to flip the Tr value
/// - Returns:
///     - whether an error occurred
void save_obj(const string& filename, const obj_scene* model,
    bool save_textures = false, bool skip_missing = false,
    bool flip_texcoord = true, bool flip_tr = true);

/// Shape. May contain only one of the points/lines/triangles.
struct obj_shape {
    /// name of the group that enclosed it
    string name = "";
    /// name of the material
    string matname = "";

    // shape elements -------------------------
    /// points
    vector<int> points;
    /// lines
    vector<vec2i> lines;
    /// triangles
    vector<vec3i> triangles;
    /// tetrahedrons
    vector<vec4i> tetras;

    // vertex data ----------------------------
    /// per-vertex position (3 float)
    vector<vec3f> pos;
    /// per-vertex normals (3 float)
    vector<vec3f> norm;
    /// per-vertex texcoord (2 float)
    vector<vec2f> texcoord;
    /// [extension] per-vertex color (4 float)
    vector<vec4f> color;
    /// [extension] per-vertex radius (1 float)
    vector<float> radius;
};

/// Mesh
struct obj_mesh {
    // name
    string name;
    /// primitives
    vector<obj_shape> shapes;

    /// cleanup
    ~obj_mesh();
};

/// Gets a mesh from an OBJ object.
obj_mesh* get_mesh(
    const obj_scene* model, const obj_object& oobj, bool facet_non_smooth);

}  // namespace ygl

#if YGL_GLTF

// include json for glTF
#if YGL_GLTFJSON
#include "json.hpp"
#endif

// -----------------------------------------------------------------------------
// KHRONOS GLTF SUPPORT
// -----------------------------------------------------------------------------
namespace ygl {

/// Generic buffer data.
using buffer_data = vector<unsigned char>;

/// Generic image data.
struct image_data {
    /// Width
    int width = 0;
    /// Height
    int height = 0;
    /// Number of Channels
    int ncomp = 0;
    /// Buffer data for 8-bit images
    vector<uint8_t> datab;
    /// Buffer data for float images
    vector<float> dataf;
};

/// glTFid
template <typename T>
struct glTFid {
    /// defaoult constructor to an invalid id
    glTFid() : _id(-1) {}
    /// explicit conversion from integer
    explicit glTFid(int id) : _id(id) {}
    /// explicit convcersion to integer
    explicit operator int() const { return _id; }

    /// check if it is valid
    bool is_valid() const { return _id >= 0; }
    /// check if it is valid
    explicit operator bool() const { return _id >= 0; }

   private:
    // id
    int _id = -1;
};

/// Generic glTF object
struct glTFProperty {
#if YGL_GLTFJSON
    /// Extensions.
    map<string, nlohmann::json> extensions = {};
    /// Extra data.
    nlohmann::json extras = {};
#endif
};

// #codegen begin type

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

/// Generic glTF named object
struct glTFChildOfRootProperty : glTFProperty {
    /// The user-defined name of this object.
    string name = "";
};

/// Values for glTFAccessorSparseIndices::componentType
enum class glTFAccessorSparseIndicesComponentType {
    /// Not set
    NotSet = -1,
    // UnsignedByte
    UnsignedByte = 5121,
    // UnsignedShort
    UnsignedShort = 5123,
    // UnsignedInt
    UnsignedInt = 5125,
};

/// Indices of those attributes that deviate from their initialization value.
struct glTFAccessorSparseIndices : glTFProperty {
    /// The index of the bufferView with sparse indices. Referenced bufferView
    /// can't have ARRAY_BUFFER or ELEMENT_ARRAY_BUFFER target. [required]
    glTFid<glTFBufferView> bufferView = {};
    /// The offset relative to the start of the bufferView in bytes. Must be
    /// aligned.
    int byteOffset = 0;
    /// The indices data type. [required]
    glTFAccessorSparseIndicesComponentType componentType =
        glTFAccessorSparseIndicesComponentType::NotSet;
};

/// Array of size `accessor.sparse.count` times number of components storing the
/// displaced accessor attributes pointed by `accessor.sparse.indices`.
struct glTFAccessorSparseValues : glTFProperty {
    /// The index of the bufferView with sparse values. Referenced bufferView
    /// can't have ARRAY_BUFFER or ELEMENT_ARRAY_BUFFER target. [required]
    glTFid<glTFBufferView> bufferView = {};
    /// The offset relative to the start of the bufferView in bytes. Must be
    /// aligned.
    int byteOffset = 0;
};

/// Sparse storage of attributes that deviate from their initialization value.
struct glTFAccessorSparse : glTFProperty {
    /// Number of entries stored in the sparse array. [required]
    int count = 0;
    /// Index array of size `count` that points to those accessor attributes
    /// that deviate from their initialization value. Indices must strictly
    /// increase. [required]
    glTFAccessorSparseIndices* indices = nullptr;
    /// Array of size `count` times number of components, storing the displaced
    /// accessor attributes pointed by `indices`. Substituted values must have
    /// the same `componentType` and number of components as the base accessor.
    /// [required]
    glTFAccessorSparseValues* values = nullptr;

    /// Cleanup
    ~glTFAccessorSparse() {
        if (indices) delete indices;
        if (values) delete values;
    }
};

/// Values for glTFAccessor::componentType
enum class glTFAccessorComponentType {
    /// Not set
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

/// Values for glTFAccessor::type
enum class glTFAccessorType {
    /// Not set
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

/// A typed view into a bufferView.  A bufferView contains raw binary data.  An
/// accessor provides a typed view into a bufferView or a subset of a bufferView
/// similar to how WebGL's `vertexAttribPointer()` defines an attribute in a
/// buffer.
struct glTFAccessor : glTFChildOfRootProperty {
    /// The index of the bufferView.
    glTFid<glTFBufferView> bufferView = {};
    /// The offset relative to the start of the bufferView in bytes.
    int byteOffset = 0;
    /// The datatype of components in the attribute. [required]
    glTFAccessorComponentType componentType = glTFAccessorComponentType::NotSet;
    /// Specifies whether integer data values should be normalized.
    bool normalized = false;
    /// The number of attributes referenced by this accessor. [required]
    int count = 0;
    /// Specifies if the attribute is a scalar, vector, or matrix. [required]
    glTFAccessorType type = glTFAccessorType::NotSet;
    /// Maximum value of each component in this attribute.
    vector<float> max = {};
    /// Minimum value of each component in this attribute.
    vector<float> min = {};
    /// Sparse storage of attributes that deviate from their initialization
    /// value.
    glTFAccessorSparse* sparse = nullptr;

    /// Cleanup
    ~glTFAccessor() {
        if (sparse) delete sparse;
    }
};

/// Values for glTFAnimationChannelTarget::path
enum class glTFAnimationChannelTargetPath {
    /// Not set
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

/// The index of the node and TRS property that an animation channel targets.
struct glTFAnimationChannelTarget : glTFProperty {
    /// The index of the node to target. [required]
    glTFid<glTFNode> node = {};
    /// The name of the node's TRS property to modify, or the "weights" of the
    /// Morph Targets it instantiates. [required]
    glTFAnimationChannelTargetPath path =
        glTFAnimationChannelTargetPath::NotSet;
};

/// Targets an animation's sampler at a node's property.
struct glTFAnimationChannel : glTFProperty {
    /// The index of a sampler in this animation used to compute the value for
    /// the target. [required]
    glTFid<glTFAnimationSampler> sampler = {};
    /// The index of the node and TRS property to target. [required]
    glTFAnimationChannelTarget* target = nullptr;

    /// Cleanup
    ~glTFAnimationChannel() {
        if (target) delete target;
    }
};

/// Values for glTFAnimationSampler::interpolation
enum class glTFAnimationSamplerInterpolation {
    /// Not set
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
    // The animation's interpolation is computed using a uniform Catmull-Rom
    // spline. The number of output elements must equal two more than the number
    // of input elements. The first and last output elements represent the start
    // and end tangents of the spline. There must be at least four keyframes
    // when using this interpolation.
    CatmullRomSpline = 2,
    // The animation's interpolation is computed using a cubic spline with
    // specified tangents. The number of output elements must equal three times
    // the number of input elements. For each input element, the output stores
    // three elements, an in-tangent, a spline vertex, and an out-tangent. There
    // must be at least two keyframes when using this interpolation.
    CubicSpline = 3,
};

/// Combines input and output accessors with an interpolation algorithm to
/// define a keyframe graph (but not its target).
struct glTFAnimationSampler : glTFProperty {
    /// The index of an accessor containing keyframe input values, e.g., time.
    /// [required]
    glTFid<glTFAccessor> input = {};
    /// Interpolation algorithm.
    glTFAnimationSamplerInterpolation interpolation =
        glTFAnimationSamplerInterpolation::Linear;
    /// The index of an accessor, containing keyframe output values. [required]
    glTFid<glTFAccessor> output = {};
};

/// A keyframe animation.
struct glTFAnimation : glTFChildOfRootProperty {
    /// An array of channels, each of which targets an animation's sampler at a
    /// node's property. Different channels of the same animation can't have
    /// equal targets. [required]
    vector<glTFAnimationChannel*> channels = {};
    /// An array of samplers that combines input and output accessors with an
    /// interpolation algorithm to define a keyframe graph (but not its target).
    /// [required]
    vector<glTFAnimationSampler*> samplers = {};

    /// typed access for nodes
    glTFAnimationChannel* get(const glTFid<glTFAnimationChannel>& id) const {
        if (!id) return nullptr;
        return channels.at((int)id);
    }
    /// typed access for nodes
    glTFAnimationSampler* get(const glTFid<glTFAnimationSampler>& id) const {
        if (!id) return nullptr;
        return samplers.at((int)id);
    }
    /// Cleanup
    ~glTFAnimation() {
        for (auto v : channels)
            if (v) delete v;
        for (auto v : samplers)
            if (v) delete v;
    }
};

/// Metadata about the glTF asset.
struct glTFAsset : glTFProperty {
    /// A copyright message suitable for display to credit the content creator.
    string copyright = "";
    /// Tool that generated this glTF model.  Useful for debugging.
    string generator = "";
    /// The glTF version that this asset targets. [required]
    string version = "";
    /// The minimum glTF version that this asset targets.
    string minVersion = "";
};

/// A buffer points to binary geometry, animation, or skins.
struct glTFBuffer : glTFChildOfRootProperty {
    /// The uri of the buffer.
    string uri = "";
    /// The length of the buffer in bytes. [required]
    int byteLength = 0;
    /// Stores buffer content after loading. [required]
    buffer_data data = {};
};

/// Values for glTFBufferView::target
enum class glTFBufferViewTarget {
    /// Not set
    NotSet = -1,
    // ArrayBuffer
    ArrayBuffer = 34962,
    // ElementArrayBuffer
    ElementArrayBuffer = 34963,
};

/// A view into a buffer generally representing a subset of the buffer.
struct glTFBufferView : glTFChildOfRootProperty {
    /// The index of the buffer. [required]
    glTFid<glTFBuffer> buffer = {};
    /// The offset into the buffer in bytes.
    int byteOffset = 0;
    /// The length of the bufferView in bytes. [required]
    int byteLength = 0;
    /// The stride, in bytes.
    int byteStride = 0;
    /// The target that the GPU buffer should be bound to.
    glTFBufferViewTarget target = glTFBufferViewTarget::NotSet;
};

/// An orthographic camera containing properties to create an orthographic
/// projection matrix.
struct glTFCameraOrthographic : glTFProperty {
    /// The floating-point horizontal magnification of the view. [required]
    float xmag = 0;
    /// The floating-point vertical magnification of the view. [required]
    float ymag = 0;
    /// The floating-point distance to the far clipping plane. `zfar` must be
    /// greater than `znear`. [required]
    float zfar = 0;
    /// The floating-point distance to the near clipping plane. [required]
    float znear = 0;
};

/// A perspective camera containing properties to create a perspective
/// projection matrix.
struct glTFCameraPerspective : glTFProperty {
    /// The floating-point aspect ratio of the field of view.
    float aspectRatio = 0;
    /// The floating-point vertical field of view in radians. [required]
    float yfov = 0;
    /// The floating-point distance to the far clipping plane.
    float zfar = 0;
    /// The floating-point distance to the near clipping plane. [required]
    float znear = 0;
};

/// Values for glTFCamera::type
enum class glTFCameraType {
    /// Not set
    NotSet = -1,
    // Perspective
    Perspective = 0,
    // Orthographic
    Orthographic = 1,
};

/// A camera's projection.  A node can reference a camera to apply a transform
/// to place the camera in the scene.
struct glTFCamera : glTFChildOfRootProperty {
    /// An orthographic camera containing properties to create an orthographic
    /// projection matrix.
    glTFCameraOrthographic* orthographic = nullptr;
    /// A perspective camera containing properties to create a perspective
    /// projection matrix.
    glTFCameraPerspective* perspective = nullptr;
    /// Specifies if the camera uses a perspective or orthographic projection.
    /// [required]
    glTFCameraType type = glTFCameraType::NotSet;

    /// Cleanup
    ~glTFCamera() {
        if (orthographic) delete orthographic;
        if (perspective) delete perspective;
    }
};

/// Values for glTFImage::mimeType
enum class glTFImageMimeType {
    /// Not set
    NotSet = -1,
    // ImageJpeg
    ImageJpeg = 0,
    // ImagePng
    ImagePng = 1,
};

/// Image data used to create a texture. Image can be referenced by URI or
/// `bufferView` index. `mimeType` is required in the latter case.
struct glTFImage : glTFChildOfRootProperty {
    /// The uri of the image.
    string uri = "";
    /// The image's MIME type.
    glTFImageMimeType mimeType = glTFImageMimeType::NotSet;
    /// The index of the bufferView that contains the image. Use this instead of
    /// the image's uri property.
    glTFid<glTFBufferView> bufferView = {};
    /// Stores image content after loading.
    image_data data = {};
};

/// Reference to a texture.
struct glTFTextureInfo : glTFProperty {
    /// The index of the texture. [required]
    glTFid<glTFTexture> index = {};
    /// The set index of texture's TEXCOORD attribute used for texture
    /// coordinate mapping.
    int texCoord = 0;
};

/// A texture and its sampler.
struct glTFTexture : glTFChildOfRootProperty {
    /// The index of the sampler used by this texture. When undefined, a sampler
    /// with repeat wrapping and auto filtering should be used.
    glTFid<glTFSampler> sampler = {};
    /// The index of the image used by this texture.
    glTFid<glTFImage> source = {};
};

/// Normal texture information.
struct glTFMaterialNormalTextureInfo : glTFTextureInfo {
    /// The scalar multiplier applied to each normal vector of the normal
    /// texture.
    float scale = 1;
};

/// Occlusion texture information.
struct glTFMaterialOcclusionTextureInfo : glTFTextureInfo {
    /// A scalar multiplier controlling the amount of occlusion applied.
    float strength = 1;
};

/// A set of parameter values that are used to define the metallic-roughness
/// material model from Physically-Based Rendering (PBR) methodology.
struct glTFMaterialPbrMetallicRoughness : glTFProperty {
    /// The material's base color factor.
    vec4f baseColorFactor = {1, 1, 1, 1};
    /// The base color texture.
    glTFTextureInfo* baseColorTexture = nullptr;
    /// The metalness of the material.
    float metallicFactor = 1;
    /// The roughness of the material.
    float roughnessFactor = 1;
    /// The metallic-roughness texture.
    glTFTextureInfo* metallicRoughnessTexture = nullptr;

    /// Cleanup
    ~glTFMaterialPbrMetallicRoughness() {
        if (baseColorTexture) delete baseColorTexture;
        if (metallicRoughnessTexture) delete metallicRoughnessTexture;
    }
};

/// glTF extension that defines the specular-glossiness material model from
/// Physically-Based Rendering (PBR) methodology.
struct glTFMaterialPbrSpecularGlossiness : glTFProperty {
    /// The reflected diffuse factor of the material.
    vec4f diffuseFactor = {1, 1, 1, 1};
    /// The diffuse texture.
    glTFTextureInfo* diffuseTexture = nullptr;
    /// The specular RGB color of the material.
    vec3f specularFactor = {1, 1, 1};
    /// The glossiness or smoothness of the material.
    float glossinessFactor = 1;
    /// The specular-glossiness texture.
    glTFTextureInfo* specularGlossinessTexture = nullptr;

    /// Cleanup
    ~glTFMaterialPbrSpecularGlossiness() {
        if (diffuseTexture) delete diffuseTexture;
        if (specularGlossinessTexture) delete specularGlossinessTexture;
    }
};

/// Values for glTFMaterial::alphaMode
enum class glTFMaterialAlphaMode {
    /// Not set
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

/// The material appearance of a primitive.
struct glTFMaterial : glTFChildOfRootProperty {
    /// A set of parameter values that are used to define the metallic-roughness
    /// material model from Physically-Based Rendering (PBR) methodology. When
    /// not specified, all the default values of `pbrMetallicRoughness` apply.
    glTFMaterialPbrMetallicRoughness* pbrMetallicRoughness = nullptr;
    /// A set of parameter values that are used to define the
    /// specular-glossiness material model from Physically-Based Rendering (PBR)
    /// methodology. When not specified, all the default values of
    /// `pbrMetallicRoughness` apply.
    glTFMaterialPbrSpecularGlossiness* pbrSpecularGlossiness = nullptr;
    /// The normal map texture.
    glTFMaterialNormalTextureInfo* normalTexture = nullptr;
    /// The occlusion map texture.
    glTFMaterialOcclusionTextureInfo* occlusionTexture = nullptr;
    /// The emissive map texture.
    glTFTextureInfo* emissiveTexture = nullptr;
    /// The emissive color of the material.
    vec3f emissiveFactor = {0, 0, 0};
    /// The alpha rendering mode of the material.
    glTFMaterialAlphaMode alphaMode = glTFMaterialAlphaMode::Opaque;
    /// The alpha cutoff value of the material.
    float alphaCutoff = 0.5;
    /// Specifies whether the material is double sided.
    bool doubleSided = false;

    /// Cleanup
    ~glTFMaterial() {
        if (pbrMetallicRoughness) delete pbrMetallicRoughness;
        if (pbrSpecularGlossiness) delete pbrSpecularGlossiness;
        if (normalTexture) delete normalTexture;
        if (occlusionTexture) delete occlusionTexture;
        if (emissiveTexture) delete emissiveTexture;
    }
};

/// Values for glTFMeshPrimitive::mode
enum class glTFMeshPrimitiveMode {
    /// Not set
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

/// Geometry to be rendered with the given material.
struct glTFMeshPrimitive : glTFProperty {
    /// A dictionary object, where each key corresponds to mesh attribute
    /// semantic and each value is the index of the accessor containing
    /// attribute's data. [required]
    map<string, glTFid<glTFAccessor>> attributes = {};
    /// The index of the accessor that contains the indices.
    glTFid<glTFAccessor> indices = {};
    /// The index of the material to apply to this primitive when rendering.
    glTFid<glTFMaterial> material = {};
    /// The type of primitives to render.
    glTFMeshPrimitiveMode mode = glTFMeshPrimitiveMode::Triangles;
    /// An array of Morph Targets, each  Morph Target is a dictionary mapping
    /// attributes (only `POSITION`, `NORMAL`, and `TANGENT` supported) to their
    /// deviations in the Morph Target.
    vector<map<string, glTFid<glTFAccessor>>> targets = {};
};

/// A set of primitives to be rendered.  A node can contain one mesh.  A node's
/// transform places the mesh in the scene.
struct glTFMesh : glTFChildOfRootProperty {
    /// An array of primitives, each defining geometry to be rendered with a
    /// material. [required]
    vector<glTFMeshPrimitive*> primitives = {};
    /// Array of weights to be applied to the Morph Targets.
    vector<float> weights = {};

    /// Cleanup
    ~glTFMesh() {
        for (auto v : primitives)
            if (v) delete v;
    }
};

/// A node in the node hierarchy.  When the node contains `skin`, all
/// `mesh.primitives` must contain `JOINTS_0` and `WEIGHTS_0` attributes.  A
/// node can have either a `matrix` or any combination of
/// `translation`/`rotation`/`scale` (TRS) properties. TRS properties are
/// converted to matrices and postmultiplied in the `T * R * S` order to compose
/// the transformation matrix; first the scale is applied to the vertices, then
/// the rotation, and then the translation. If none are provided, the transform
/// is the identity. When a node is targeted for animation (referenced by an
/// animation.channel.target), only TRS properties may be present; `matrix` will
/// not be present.
struct glTFNode : glTFChildOfRootProperty {
    /// The index of the camera referenced by this node.
    glTFid<glTFCamera> camera = {};
    /// The indices of this node's children.
    vector<glTFid<glTFNode>> children = {};
    /// The index of the skin referenced by this node.
    glTFid<glTFSkin> skin = {};
    /// A floating-point 4x4 transformation matrix stored in column-major order.
    mat4f matrix = {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
    /// The index of the mesh in this node.
    glTFid<glTFMesh> mesh = {};
    /// The node's unit quaternion rotation in the order (x, y, z, w), where w
    /// is the scalar.
    quat4f rotation = {0, 0, 0, 1};
    /// The node's non-uniform scale.
    vec3f scale = {1, 1, 1};
    /// The node's translation.
    vec3f translation = {0, 0, 0};
    /// The weights of the instantiated Morph Target. Number of elements must
    /// match number of Morph Targets of used mesh.
    vector<float> weights = {};
};

/// Values for glTFSampler::magFilter
enum class glTFSamplerMagFilter {
    /// Not set
    NotSet = -1,
    // Nearest
    Nearest = 9728,
    // Linear
    Linear = 9729,
};

/// Values for glTFSampler::minFilter
enum class glTFSamplerMinFilter {
    /// Not set
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

/// glTFSampler::wrapS
enum class glTFSamplerWrapS {
    /// Not set
    NotSet = -1,
    // ClampToEdge
    ClampToEdge = 33071,
    // MirroredRepeat
    MirroredRepeat = 33648,
    // Repeat
    Repeat = 10497,
};

/// glTFSampler::wrapT
enum class glTFSamplerWrapT {
    /// Not set
    NotSet = -1,
    // ClampToEdge
    ClampToEdge = 33071,
    // MirroredRepeat
    MirroredRepeat = 33648,
    // Repeat
    Repeat = 10497,
};

/// Texture sampler properties for filtering and wrapping modes.
struct glTFSampler : glTFChildOfRootProperty {
    /// Magnification filter.
    glTFSamplerMagFilter magFilter = glTFSamplerMagFilter::NotSet;
    /// Minification filter.
    glTFSamplerMinFilter minFilter = glTFSamplerMinFilter::NotSet;
    /// s wrapping mode.
    glTFSamplerWrapS wrapS = glTFSamplerWrapS::Repeat;
    /// t wrapping mode.
    glTFSamplerWrapT wrapT = glTFSamplerWrapT::Repeat;
};

/// The root nodes of a scene.
struct glTFScene : glTFChildOfRootProperty {
    /// The indices of each root node.
    vector<glTFid<glTFNode>> nodes = {};
};

/// Joints and matrices defining a skin.
struct glTFSkin : glTFChildOfRootProperty {
    /// The index of the accessor containing the floating-point 4x4 inverse-bind
    /// matrices.  The default is that each matrix is a 4x4 identity matrix,
    /// which implies that inverse-bind matrices were pre-applied.
    glTFid<glTFAccessor> inverseBindMatrices = {};
    /// The index of the node used as a skeleton root. When undefined, joints
    /// transforms resolve to scene root.
    glTFid<glTFNode> skeleton = {};
    /// Indices of skeleton nodes, used as joints in this skin. [required]
    vector<glTFid<glTFNode>> joints = {};
};

/// The root object for a glTF asset.
struct glTF : glTFProperty {
    /// Names of glTF extensions used somewhere in this asset.
    vector<string> extensionsUsed = {};
    /// Names of glTF extensions required to properly load this asset.
    vector<string> extensionsRequired = {};
    /// An array of accessors.
    vector<glTFAccessor*> accessors = {};
    /// An array of keyframe animations.
    vector<glTFAnimation*> animations = {};
    /// Metadata about the glTF asset. [required]
    glTFAsset* asset = nullptr;
    /// An array of buffers.
    vector<glTFBuffer*> buffers = {};
    /// An array of bufferViews.
    vector<glTFBufferView*> bufferViews = {};
    /// An array of cameras.
    vector<glTFCamera*> cameras = {};
    /// An array of images.
    vector<glTFImage*> images = {};
    /// An array of materials.
    vector<glTFMaterial*> materials = {};
    /// An array of meshes.
    vector<glTFMesh*> meshes = {};
    /// An array of nodes.
    vector<glTFNode*> nodes = {};
    /// An array of samplers.
    vector<glTFSampler*> samplers = {};
    /// The index of the default scene.
    glTFid<glTFScene> scene = {};
    /// An array of scenes.
    vector<glTFScene*> scenes = {};
    /// An array of skins.
    vector<glTFSkin*> skins = {};
    /// An array of textures.
    vector<glTFTexture*> textures = {};

    /// typed access for nodes
    glTFAccessor* get(const glTFid<glTFAccessor>& id) const {
        if (!id) return nullptr;
        return accessors.at((int)id);
    }
    /// typed access for nodes
    glTFAnimation* get(const glTFid<glTFAnimation>& id) const {
        if (!id) return nullptr;
        return animations.at((int)id);
    }
    /// typed access for nodes
    glTFBuffer* get(const glTFid<glTFBuffer>& id) const {
        if (!id) return nullptr;
        return buffers.at((int)id);
    }
    /// typed access for nodes
    glTFBufferView* get(const glTFid<glTFBufferView>& id) const {
        if (!id) return nullptr;
        return bufferViews.at((int)id);
    }
    /// typed access for nodes
    glTFCamera* get(const glTFid<glTFCamera>& id) const {
        if (!id) return nullptr;
        return cameras.at((int)id);
    }
    /// typed access for nodes
    glTFImage* get(const glTFid<glTFImage>& id) const {
        if (!id) return nullptr;
        return images.at((int)id);
    }
    /// typed access for nodes
    glTFMaterial* get(const glTFid<glTFMaterial>& id) const {
        if (!id) return nullptr;
        return materials.at((int)id);
    }
    /// typed access for nodes
    glTFMesh* get(const glTFid<glTFMesh>& id) const {
        if (!id) return nullptr;
        return meshes.at((int)id);
    }
    /// typed access for nodes
    glTFNode* get(const glTFid<glTFNode>& id) const {
        if (!id) return nullptr;
        return nodes.at((int)id);
    }
    /// typed access for nodes
    glTFSampler* get(const glTFid<glTFSampler>& id) const {
        if (!id) return nullptr;
        return samplers.at((int)id);
    }
    /// typed access for nodes
    glTFScene* get(const glTFid<glTFScene>& id) const {
        if (!id) return nullptr;
        return scenes.at((int)id);
    }
    /// typed access for nodes
    glTFSkin* get(const glTFid<glTFSkin>& id) const {
        if (!id) return nullptr;
        return skins.at((int)id);
    }
    /// typed access for nodes
    glTFTexture* get(const glTFid<glTFTexture>& id) const {
        if (!id) return nullptr;
        return textures.at((int)id);
    }
    /// Cleanup
    ~glTF() {
        for (auto v : accessors)
            if (v) delete v;
        for (auto v : animations)
            if (v) delete v;
        if (asset) delete asset;
        for (auto v : buffers)
            if (v) delete v;
        for (auto v : bufferViews)
            if (v) delete v;
        for (auto v : cameras)
            if (v) delete v;
        for (auto v : images)
            if (v) delete v;
        for (auto v : materials)
            if (v) delete v;
        for (auto v : meshes)
            if (v) delete v;
        for (auto v : nodes)
            if (v) delete v;
        for (auto v : samplers)
            if (v) delete v;
        for (auto v : scenes)
            if (v) delete v;
        for (auto v : skins)
            if (v) delete v;
        for (auto v : textures)
            if (v) delete v;
    }
};
// #codegen end type
// -----------------------------------------------------------

/// Loads a gltf file from disk
///
/// - Parameters:
///     - filename: scene filename
///     - load_bin/load_img: load binary data
///     - skip_missing: do not throw an exception if a file is missing
/// - Returns:
///     - gltf data loaded (nullptr on error)
glTF* load_gltf(const string& filename, bool load_bin = true,
    bool load_img = false, bool skip_missing = false);

/// Loads a binary gltf file from disk
///
/// - Parameters:
///     - filename: scene filename
///     - other params as above
/// - Returns:
///     - gltf data loaded (nullptr on error)
glTF* load_binary_gltf(const string& filename, bool load_bin = true,
    bool load_img = false, bool skip_missing = false);

/// Saves a scene to disk
///
/// - Parameters:
///     - filename: scene filename
///     - gltf: data to save
///     - save_bin/save_images: save binary data
void save_gltf(const string& filename, const glTF* gltf, bool save_bin = true,
    bool save_images = false);

/// Saves a scene to disk
///
/// - Parameters:
///     - filename: scene filename
///     - gltf: data to save
///     - save_bin/save_images: save binary data
void save_binary_gltf(const string& filename, const glTF* gltf,
    bool save_bin = true, bool save_images = false);

/// Computes the local node transform and its inverse.
inline mat4f node_transform(const glTFNode* node) {
    return translation_mat4f(node->translation) *
           rotation_mat4f(node->rotation) * scaling_mat4f(node->scale) *
           node->matrix;
}

/// A view for gltf array buffers that allows for typed access.
struct accessor_view {
    /// construct a view from an accessor
    accessor_view(const glTF* gltf, const glTFAccessor* accessor);

    /// number of elements in the view
    int size() const { return _size; }
    /// number of elements in the view
    int count() const { return _size; }
    /// number of components per element
    int ncomp() const { return _ncomp; }
    /// check whether the view is valid
    bool valid() const { return _valid; }

    /// get the idx-th element of fixed length width default values
    vec2f getv2f(int idx, const vec2f& def = {0, 0}) const {
        auto v = def;
        for (auto i = 0; i < min(_ncomp, 2); i++) v[i] = get(idx, i);
        return v;
    }
    /// get the idx-th element of fixed length width default values
    vec3f getv3f(int idx, const vec3f& def = {0, 0, 0}) const {
        auto v = def;
        for (auto i = 0; i < min(_ncomp, 3); i++) v[i] = get(idx, i);
        return v;
    }
    /// get the idx-th element of fixed length width default values
    vec4f getv4f(int idx, const vec4f& def = {0, 0, 0, 0}) const {
        auto v = def;
        for (auto i = 0; i < min(_ncomp, 4); i++) v[i] = get(idx, i);
        return v;
    }

    /// get the idx-th element of fixed length as a matrix
    mat4f getm4f(int idx) const {
        auto v = mat4f();
        assert(_ncomp == 16);
        for (auto j = 0; j < 4; j++)
            for (auto i = 0; i < 4; i++) v[j][i] = get(idx, j * 4 + i);
        return v;
    }

    /// get the c-th component of the idx-th element
    float get(int idx, int c = 0) const;

    /// get the idx-th element as integer with fixed length
    vec2i getv2i(int idx, const vec2i& def = {0, 0}) const {
        auto v = def;
        for (auto i = 0; i < min(_ncomp, 2); i++) { v[i] = geti(idx, i); }
        return v;
    }
    /// get the idx-th element as integer with fixed length
    vec3i getv3i(int idx, const vec3i& def = {0, 0, 0}) const {
        auto v = def;
        for (auto i = 0; i < min(_ncomp, 3); i++) { v[i] = geti(idx, i); }
        return v;
    }
    /// get the idx-th element as integer with fixed length
    vec4i getv4i(int idx, const vec4i& def = {0, 0, 0, 0}) const {
        auto v = def;
        for (auto i = 0; i < min(_ncomp, 4); i++) { v[i] = geti(idx, i); }
        return v;
    }

    /// get the c-th component of the idx-th element as integer
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
// PYTHON-LIKE STRING, PATH AND FILE OPERATIONS
// -----------------------------------------------------------------------------
namespace ygl {

/// Checks if a string starts with a prefix.
inline bool startswith(const string& str, const string& substr) {
    if (str.length() < substr.length()) return false;
    for (auto i = 0; i < substr.length(); i++)
        if (str[i] != substr[i]) return false;
    return true;
}

/// Checks if a string ends with a prefix.
inline bool endswith(const string& str, const string& substr) {
    if (str.length() < substr.length()) return false;
    auto offset = str.length() - substr.length();
    for (auto i = 0; i < substr.length(); i++)
        if (str[i + offset] != substr[i]) return false;
    return true;
}

/// Check is a string contains a substring.
inline bool contains(const string& str, const string& substr) {
    return str.find(substr) != str.npos;
}

/// Splits a string into lines at the '\n' character. The line
/// terminator is kept if keep_newline. This function does not work on
/// Window if keep_newline is true.
inline vector<string> splitlines(const string& str, bool keep_newline = false) {
    if (str.empty()) return {};
    auto lines = vector<string>();
    auto line = vector<char>();
    for (auto c : str) {
        if (c == '\n') {
            if (keep_newline) line.push_back(c);
            lines.push_back(string(line.begin(), line.end()));
            line.clear();
        } else {
            line.push_back(c);
        }
    }
    if (!line.empty()) lines.push_back(string(line.begin(), line.end()));
    return lines;
}

/// Partition the string.
inline vector<string> partition(const string& str, const string& split) {
    auto pos = str.find(split);
    if (pos == str.npos) return {str, "", ""};
    return {str.substr(0, pos), split, str.substr(pos + split.length())};
}

/// Splits the string.
inline vector<string> split(const string& str) {
    if (str.empty()) return {};
    auto ret = vector<string>();
    auto lpos = (size_t)0;
    while (lpos != str.npos) {
        auto pos = str.find_first_of(" \t\n\r", lpos);
        if (pos != str.npos) {
            if (pos > lpos) ret.push_back(str.substr(lpos, pos - lpos));
            lpos = pos + 1;
        } else {
            if (lpos < str.size()) ret.push_back(str.substr(lpos));
            lpos = pos;
        }
    }
    return ret;
}

/// Splits the string.
inline vector<string> split(const string& str, const string& substr) {
    if (str.empty()) return {};
    auto ret = vector<string>();
    auto lpos = (size_t)0;
    while (lpos != str.npos) {
        auto pos = str.find(substr, lpos);
        if (pos != str.npos) {
            ret.push_back(str.substr(lpos, pos - lpos));
            lpos = pos + substr.size();
        } else {
            if (lpos < str.size()) ret.push_back(str.substr(lpos));
            lpos = pos;
        }
    }
    return ret;
}

/// Splits the string.
inline vector<string> split(const string& str, char substr) {
    if (str.empty()) return {};
    auto ret = vector<string>();
    auto lpos = (size_t)0;
    while (lpos != str.npos) {
        auto pos = str.find(substr, lpos);
        if (pos != str.npos) {
            ret.push_back(str.substr(lpos, pos - lpos));
            lpos = pos + 1;
        } else {
            if (lpos < str.size()) ret.push_back(str.substr(lpos));
            lpos = pos;
        }
    }
    return ret;
}

/// Strip the string.
inline string rstrip(const string& str) {
    auto pos = str.find_last_not_of(" \t\r\n");
    if (pos == str.npos) return "";
    return str.substr(0, pos + 1);
}

/// Strip the string.
inline string lstrip(const string& str) {
    auto pos = str.find_first_not_of(" \t\r\n");
    if (pos == str.npos) return "";
    return str.substr(pos);
}

/// Strip the string.
inline string strip(const string& str) { return rstrip(lstrip(str)); }

/// Joins a list of string with a string as separator.
inline string join(const vector<string>& strs, const string& sep) {
    auto ret = string();
    auto first = true;
    for (auto& str : strs) {
        if (!first) ret += sep;
        ret += str;
        first = false;
    }
    return ret;
}

/// Converts an ASCII string to lowercase.
inline string lower(const string& str) {
    auto s = str;
    for (auto& c : s) c = tolower(c);
    return s;
}

/// Converts an ASCII string to uppercase.
inline string upper(const string& str) {
    auto s = str;
    for (auto& c : s) c = toupper(c);
    return s;
}

/// Check if a string is space.
inline bool isspace(const string& str) {
    for (auto c : str) {
        if (c != ' ' && c != '\n' && c != '\t' && c != '\r') return false;
    }
    return true;
}

/// Replace s1 with s2 in str.
inline string replace(const string& str, const string& s1, const string& s2) {
    auto s = string();
    auto last = 0;
    auto pos = (int)str.find(s1);
    while (pos != str.npos) {
        s += str.substr(last, pos - last);
        s += s2;
        last = pos + (int)s1.length();
        pos = (int)str.find(s1, last);
    }
    s += str.substr(last);
    return s;
}

/// Get directory name (including '/').
inline string path_dirname(const string& filename) {
    auto pos = filename.rfind('/');
    if (pos == string::npos) pos = filename.rfind('\\');
    if (pos == string::npos) return "";
    return filename.substr(0, pos + 1);
}

/// Get extension (including '.').
inline string path_extension(const string& filename) {
    auto pos = filename.rfind('.');
    if (pos == string::npos) return "";
    return filename.substr(pos);
}

/// Get file basename.
inline string path_basename(const string& filename) {
    auto dirname = path_dirname(filename);
    auto extension = path_extension(filename);
    return filename.substr(
        dirname.size(), filename.size() - dirname.size() - extension.size());
}

/// Get filename without directory (equiv to get_basename() +
/// get_extension()).
inline string path_filename(const string& filename) {
    return path_basename(filename) + path_extension(filename);
}

/// Replace extension.
inline string replace_path_extension(
    const string& filename, const string& ext) {
    return path_dirname(filename) + path_basename(filename) + ext;
}

/// Prepend a string to the extension.
inline string prepend_path_extension(
    const string& filename, const string& prep) {
    return path_dirname(filename) + path_basename(filename) + prep +
           path_extension(filename);
}

/// Splits a path calling the above functions.
inline void split_path(
    const string& filename, string& dirname, string& basename, string& ext) {
    dirname = path_dirname(filename);
    basename = path_basename(filename);
    ext = path_extension(filename);
}

/// Really-minimal Python like string format. The implementation is not fast
/// nor memory efficient. But it is good enough for some needs.
inline string format(const string& fmt, const vector<string>& args) {
    auto open = false;
    auto cur = 0;
    auto str = string();
    for (auto c : fmt) {
        if (c == '{') {
            str += args[cur++];
            open = true;
        } else if (c == '}') {
            if (!open) throw runtime_error("bad format");
            open = false;
        } else {
            str += c;
        }
    }
    return str;
}

// Implementation of the function below
inline void _format_one(vector<string>& vals) {}
template <typename Arg, typename... Args>
inline void _format_one(
    vector<string>& vals, const Arg& arg, const Args&... args) {
    auto stream = stringstream();
    stream << arg;
    vals.push_back(stream.str());
    _format_one(vals, args...);
}

/// Really-minimal Python like string format. Internally uses streams for
/// generality and supports for now only the '{}' operator. The implementation
/// is not fast nor memory efficient. But it is good enough for some needs.
template <typename... Args>
inline string format(const string& fmt, const Args&... args) {
    auto vals = vector<string>();
    _format_one(vals, args...);
    return format(fmt, vals);
}

/// Wrapper for the above function that prints to stdout.
template <typename... Args>
inline void print(const string& fmt, const Args&... args) {
    printf("%s", format(fmt, args...).c_str());
}

/// Wrapper for the above function that prints to stdout with endline.
template <typename... Args>
inline void println(const string& fmt, const Args&... args) {
    printf("%s\n", format(fmt, args...).c_str());
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// FILE LOADING AND SAVING
// -----------------------------------------------------------------------------
namespace ygl {

/// Loads the contents of a binary file in an in-memory array.
inline vector<unsigned char> load_binfile(const string& filename) {
    fstream fs(filename, ios_base::in | ios_base::binary);
    if (fs.fail()) throw runtime_error("cannot read file " + filename);
    fs.seekg(0, std::ios::end);
    auto buf = vector<unsigned char>(fs.tellg());
    fs.seekg(0);
    fs.read((char*)buf.data(), buf.size());
    if (fs.fail() || fs.bad())
        throw runtime_error("cannot read file " + filename);
    return buf;
}

/// Loads the contents of a text file into a string.
inline string load_txtfile(const string& filename) {
    fstream fs(filename, ios_base::in);
    if (fs.fail()) throw runtime_error("cannot read file " + filename);
    stringstream ss;
    ss << fs.rdbuf();
    if (fs.fail()) throw runtime_error("cannot read file " + filename);
    return ss.str();
}

/// Saves binary data to a file.
inline void save_binfile(
    const string& filename, const vector<unsigned char>& data) {
    fstream fs(filename, ios_base::out | ios_base::binary);
    if (fs.fail()) throw runtime_error("cannot write file " + filename);
    fs.write((const char*)data.data(), data.size());
    if (fs.fail() || fs.bad())
        throw runtime_error("cannot write file " + filename);
}

/// Saves a string to a text file.
inline void save_txtfile(const string& filename, const string& str) {
    fstream fs(filename, ios_base::out);
    if (fs.fail()) throw runtime_error("cannot write file " + filename);
    fs << str;
    if (fs.fail()) throw runtime_error("cannot write file " + filename);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMMEDIATE MODE COMMAND LINE PARSER
// -----------------------------------------------------------------------------
namespace ygl {

/// Immediate mode command line parser (opaque type)
struct cmdline_parser;

/// Immediate mode command line parser
struct cmdline_parser {
    // private implementation
    vector<string> _to_parse;    // args left to parse
    vector<string> _used_names;  // used names for check
    string _usage_prog;          // usage prog line
    string _usage_help;          // usage help line
    string _usage_opts;          // usage option lines
    string _usage_args;          // usage argument lines
    bool _usage = false;         // help option triggered
    string _error;               // parse error
};

// cmdline implementation
inline void _check_name(
    cmdline_parser& parser, const string& name, const string& flag, bool opt) {
    if (opt) {
        if (name.size() < 3 || name[0] != '-' || name[1] != '-' ||
            name[2] == '-')
            throw runtime_error("bad name " + name);
    } else {
        if (name.size() < 1 || name[0] == '-')
            throw runtime_error("bad name " + name);
    }
    if (find(parser._used_names.begin(), parser._used_names.end(), name) !=
        parser._used_names.end())
        throw runtime_error("already used " + name);
    parser._used_names.push_back(name);
    if (flag.empty()) return;
    if (flag.size() < 2 || flag[0] != '-' || flag[1] == '-')
        throw runtime_error("bad name " + flag);
    if (find(parser._used_names.begin(), parser._used_names.end(), flag) !=
        parser._used_names.end())
        throw runtime_error("already used " + flag);
    parser._used_names.push_back(flag);
}

// cmdline implementation
template <typename T>
inline void _add_usage_str(cmdline_parser& parser, const string& name,
    const string& flag, bool opt, const string& help, const string& def,
    bool req, const vector<T>& choices) {
    auto stream = stringstream();
    stream << "  " << name;
    if (!flag.empty()) stream << "/" << flag;
    while (stream.str().length() < 32) stream << " ";
    stream << help << " ";
    if (!req) stream << "[" << def << "]";
    stream << "\n";
    if (!choices.empty()) {
        for (auto i = 0; i < 32; i++) stream << " ";
        stream << "(";
        auto first = true;
        for (auto&& c : choices) {
            if (!first) stream << ",";
            stream << c;
            first = false;
        }
        stream << ")";
        stream << "\n";
    }
    if (opt)
        parser._usage_opts += stream.str();
    else
        parser._usage_args += stream.str();
}

// cmdline implementation
template <typename T>
inline void _add_usage(cmdline_parser& parser, const string& name,
    const string& flag, bool opt, const string& help, const T& def, bool req,
    const vector<T>& choices) {
    auto stream = stringstream();
    stream << def;
    _add_usage_str(parser, name, flag, opt, help, stream.str(), req, choices);
}

// cmdline implementation
template <typename T>
inline void _add_usage(cmdline_parser& parser, const string& name,
    const string& flag, bool opt, const string& help, const vector<T>& def,
    bool req, const vector<T>& choices) {
    auto stream = stringstream();
    auto first = true;
    for (auto&& v : def) {
        if (!first) stream << ",";
        stream << v;
        first = false;
    }
    _add_usage_str(parser, name, flag, opt, help, stream.str(), req, choices);
}

// cmdline implementation
inline void _set_error(cmdline_parser& parser, const string& err) {
    if (parser._error.empty()) parser._error = err;
}

/// check unused arguments
inline bool should_exit(cmdline_parser& parser) {
    for (auto&& v : parser._to_parse) {
        if (v[0] == '-')
            _set_error(parser, "unknown option " + v);
        else
            _set_error(parser, "unknown argument " + v);
    }
    return !parser._error.empty();
}

/// returns the usage string
inline string get_usage(const cmdline_parser& parser) {
    auto str = string();
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

/// parse a flag from the command line
inline bool parse_flag(cmdline_parser& parser, const string& name,
    const string& flag, const string& help, bool def = false,
    bool req = false) {
    // check names
    _check_name(parser, name, flag, true);
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

/// parse an option from the command line
template <typename T>
inline T parse_opt(cmdline_parser& parser, const string& name,
    const string& flag, const string& help, const T& def = {}, bool req = false,
    const vector<T>& choices = {}) {
    // check names
    _check_name(parser, name, flag, true);
    // update usage
    _add_usage(parser, name, flag, true, help, def, req, choices);
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
    auto stream = stringstream(arg);
    stream >> val;
    if (stream.fail()) {
        _set_error(
            parser, "incorrect value \"" + arg + "\" for option " + name);
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

/// parse an enum option from the command line
template <typename T>
inline T parse_opt(cmdline_parser& parser, const string& name,
    const string& flag, const string& help,
    const vector<pair<string, T>>& key_values, const T& def, bool req = false,
    const vector<T>& choices = {}) {
    auto keys = vector<string>{};
    auto key_def = string();
    for (auto&& kv : key_values) {
        keys.push_back(kv.first);
        if (kv.second == def) key_def = kv.first;
    }
    auto key = parse_opt<string>(parser, name, flag, help, key_def, req, keys);
    if (!parser._error.empty()) return def;
    auto val = def;
    for (auto&& kv : key_values) {
        if (kv.first == key) val = kv.second;
    }
    return val;
}

// parse positional argument from the command line
template <typename T>
inline T parse_arg(cmdline_parser& parser, const string& name,
    const string& help, const T& def = {}, bool req = true,
    const vector<T>& choices = {}) {
    // check names
    _check_name(parser, name, "", false);
    // update usage
    _add_usage(parser, name, "", false, help, def, req, choices);
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
    auto stream = stringstream(arg);
    stream >> val;
    if (stream.fail()) {
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

// parse all remaining positional argument from the command line
template <typename T>
inline vector<T> parse_args(cmdline_parser& parser, const string& name,
    const string& help, const vector<T>& def = {}, bool req = true,
    const vector<T>& choices = {}) {
    // check names
    _check_name(parser, name, "", false);
    // update usage
    _add_usage(parser, name, "", false, help, def, req, choices);
    // skip if error
    if (!parser._error.empty()) return def;
    // search for all params
    auto vals = vector<T>();
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
        auto stream = stringstream(arg);
        stream >> val;
        if (stream.fail()) {
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

/// initialize the command line
inline cmdline_parser make_parser(
    int argc, char** argv, const string& prog, const string& help) {
    auto parser = cmdline_parser();
    parser._to_parse = vector<string>(argv + 1, argv + argc);
    parser._usage_prog = (prog.empty()) ? string(argv[0]) : prog;
    parser._usage_help = help;
    parser._usage =
        parse_flag(parser, "--help", "-h", "prints and help message");
    return parser;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// SIMPLE LOGGER
// -----------------------------------------------------------------------------
namespace ygl {

/// Logger object. A logger can output messages to multiple streams.
/// Use add streams commands for it.
struct logger {
    /// whether to output verbose
    bool _verbose = true;
    /// whether to output to console
    bool _console = true;
    /// file stream for stream output
    FILE* _file = nullptr;

    // cleanup
    ~logger() {
        if (_file) fclose(_file);
    }
};

/// Make a logger with an optional console stream and a verbosity level
inline logger* make_logger(bool console = true, bool verbose = true) {
    auto lgr = new logger();
    lgr->_verbose = verbose;
    lgr->_console = console;
    lgr->_file = nullptr;
    return lgr;
}

/// Add a file stream to a logger.
///
/// - Parameters:
///     - lgr: logger
///     - filename: filename
///     - append: append or write open mode for file logger
///     - short_message: whether to use a short message version
///     - output_level: output level
///     - flush_level: output level
/// - Returns:
///     - true if ok
inline void add_file_stream(logger* lgr, const string& filename, bool append) {
    lgr->_file = fopen(filename.c_str(), (append) ? "at" : "wt");
    if (!lgr->_file) throw runtime_error("could not open file " + filename);
}

/// Get default logger.
/// By default a non-verbose stdout logger is creater.
inline logger* get_default_logger() {
    static auto default_logger = new logger();
    return default_logger;
}

// Log a message. Used internally.
inline void _log_msg(logger* lgr, const string& msg, const char* type) {
    char time_buf[1024];
    auto tm = time(nullptr);
    auto ttm = localtime(&tm);  // TODO: use thread safe version

    // short message for console
    if (lgr->_console) {
        strftime(time_buf, 1024, "%H:%M:%S", ttm);
        printf("%s %s %s\n", time_buf, type, msg.c_str());
        fflush(stdout);
    }

    // long message for file
    if (lgr->_file) {
        strftime(time_buf, 1024, "%Y-%m-%d %H:%M:%S", ttm);
        fprintf(lgr->_file, "%s %s %s\n", time_buf, type, msg.c_str());
    }
}

/// Log a info message
template <typename... Args>
inline void log_info(logger* lgr, const string& msg, const Args&... args) {
    if (!lgr->_verbose) return;
    _log_msg(lgr, format(msg, args...), "INFO ");
}

/// Log a info message
template <typename... Args>
inline void log_warning(logger* lgr, const string& msg, const Args&... args) {
    if (!lgr->_verbose) return;
    _log_msg(lgr, format(msg, args...), "WARN ");
}

/// Log an error message
template <typename... Args>
inline void log_error(logger* lgr, const string& msg, const Args&... args) {
    _log_msg(lgr, format(msg, args...), "ERROR");
}

/// Log a fatal message and exit
template <typename... Args>
inline void log_fatal(logger* lgr, const string& msg, const Args&... args) {
    _log_msg(lgr, format(msg, args...), "FATAL");
    exit(1);
}

/// Adds a file stream to the default logger
inline void add_file_stream(const string& filename, bool append) {
    add_file_stream(get_default_logger(), filename, append);
}

/// Logs a message to the default loggers
template <typename... Args>
inline void log_info(const string& msg, const Args&... args) {
    log_info(get_default_logger(), msg, args...);
}

/// Logs a message to the default loggers
template <typename... Args>
inline void log_error(const string& msg, const Args&... args) {
    log_error(get_default_logger(), msg, args...);
}

/// Logs a message to the default loggers
template <typename... Args>
inline void log_fatal(const string& msg, const Args&... args) {
    log_fatal(get_default_logger(), msg, args...);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// THREAD POOL
// -----------------------------------------------------------------------------
namespace ygl {

/// Thread pool for concurrency. This code is derived from LLVM ThreadPool
struct thread_pool {
    // initialize the thread pool
    thread_pool(int nthreads = std::thread::hardware_concurrency())
        : _working_threads(0), _stop_flag(false) {
        _threads.reserve(nthreads);
        for (auto tid = 0; tid < nthreads; tid++) {
            _threads.emplace_back([this] { _thread_proc(); });
        }
    }

    // cleanup
    ~thread_pool() {
        {
            std::unique_lock<std::mutex> lock_guard(_queue_lock);
            _stop_flag = true;
        }
        _queue_condition.notify_all();
        for (auto& Worker : _threads) Worker.join();
    }

    // empty the queue
    void _clear_pool() {
        {
            std::unique_lock<std::mutex> lock_guard(_queue_lock);
            _tasks.clear();
        }
        _queue_condition.notify_all();
    }

    // schedule an asynchronous taks
    std::shared_future<void> _run_async(std::function<void()> task) {
        // Wrap the Task in a packaged_task to return a future object.
        std::packaged_task<void()> packaged_task(std::move(task));
        auto future = packaged_task.get_future();
        {
            std::unique_lock<std::mutex> lock_guard(_queue_lock);
            assert(!_stop_flag &&
                   "Queuing a thread during ThreadPool destruction");
            _tasks.push_back(std::move(packaged_task));
        }
        _queue_condition.notify_one();
        return future.share();
    }

    // wait for all tasks to finish
    void _wait() {
        std::unique_lock<std::mutex> lock_guard(_completion_lock);
        _completion_condition.wait(
            lock_guard, [&] { return _tasks.empty() && !_working_threads; });
    }

    // parallel for
    void _parallel_for(int count, const function<void(int idx)>& task) {
        for (auto idx = 0; idx < count; idx++) {
            _run_async([&task, idx]() { task(idx); });
        }
        _wait();
    }

    // implementation -------------------------------------------------
    void _thread_proc() {
        while (true) {
            std::packaged_task<void()> task;
            {
                std::unique_lock<std::mutex> lock_guard(_queue_lock);
                _queue_condition.wait(
                    lock_guard, [&] { return _stop_flag || !_tasks.empty(); });

                if (_stop_flag && _tasks.empty()) return;

                {
                    _working_threads++;
                    std::unique_lock<std::mutex> lock_guard(_completion_lock);
                }
                task = std::move(_tasks.front());
                _tasks.pop_front();
            }

            task();

            {
                std::unique_lock<std::mutex> lock_guard(_completion_lock);
                _working_threads--;
            }

            _completion_condition.notify_all();
        }
    }

    vector<std::thread> _threads;
    std::deque<std::packaged_task<void()>> _tasks;
    std::mutex _queue_lock;
    std::condition_variable _queue_condition;
    std::mutex _completion_lock;
    std::condition_variable _completion_condition;
    std::atomic<unsigned> _working_threads;
    bool _stop_flag = false;
};

/// Makes a thread pool
inline thread_pool* make_pool(
    int nthreads = std::thread::hardware_concurrency()) {
    return new thread_pool(nthreads);
}

/// Runs a task asynchronously onto the global thread pool
inline std::shared_future<void> run_async(
    thread_pool* pool, const function<void()>& task) {
    return pool->_run_async(task);
}

/// Wait for all jobs to finish on the global thread pool
inline void wait_pool(thread_pool* pool) { pool->_wait(); }

/// Clear all jobs on the global thread pool
inline void clear_pool(thread_pool* pool) { pool->_clear_pool(); }

/// Parallel for implementation on the global thread pool
inline void parallel_for(
    thread_pool* pool, int count, const function<void(int idx)>& task) {
    pool->_parallel_for(count, task);
}

/// Global pool
inline thread_pool* get_global_pool() {
    static auto pool = (thread_pool*)nullptr;
    if (!pool) pool = new thread_pool();
    return pool;
}

/// Runs a task asynchronously onto the global thread pool
inline std::shared_future<void> run_async(const function<void()>& task) {
    return run_async(get_global_pool(), task);
}

/// Wait for all jobs to finish on the global thread pool
inline void wait_pool() { wait_pool(get_global_pool()); }

/// Clear all jobs on the global thread pool
inline void clear_pool() { clear_pool(get_global_pool()); }

/// Parallel for implementation on the global thread pool
inline void parallel_for(int count, const function<void(int idx)>& task) {
    parallel_for(get_global_pool(), count, task);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// TIMER
// -----------------------------------------------------------------------------
namespace ygl {

/// A simple wrapper for std::chrono.
struct timer {
    /// initialize a timer and start it if necessary
    timer(bool autostart = true) {
        if (autostart) start();
    }

    /// start a timer
    void start() {
        _start = std::chrono::steady_clock::now();
        _started = true;
    }

    /// stops a timer
    void stop() {
        _end = std::chrono::steady_clock::now();
        _started = false;
    }

    /// elapsed time
    double elapsed() {
        if (_started) stop();
        std::chrono::duration<double> diff = (_end - _start);
        return diff.count();
    }

   private:
    bool _started = false;
    std::chrono::time_point<std::chrono::steady_clock> _start, _end;
};

}  // namespace ygl

#if YGL_OPENGL

// -----------------------------------------------------------------------------
// OPENGL FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

/// Shape types
enum struct gl_etype : int {
    /// points
    point = 1,
    /// lines
    line = 2,
    /// triangles
    triangle = 3,
    /// quads
    quad = 4,
};

/// Light types
enum struct gl_ltype : int {
    /// point lights
    point = 0,
    /// directional lights
    directional = 1,
};

/// Checks for GL error and then prints
bool gl_check_error(bool print = true);

/// Clear window
void gl_clear_buffers(const vec4f& background = {0, 0, 0, 0});

/// Enable/disable depth test
void gl_enable_depth_test(bool enabled);

/// Enable/disable culling
void gl_enable_culling(bool enabled);

/// Enable/disable wireframe
void gl_enable_wireframe(bool enabled);

/// Enable/disable edges. Attempts to avoid z-fighting but the method is not
/// robust.
void gl_enable_edges(bool enabled, float tolerance = 0.9999f);

/// Enable/disable blending
void gl_enable_blending(bool enabled);

/// Set blending to over operator
void gl_set_blend_over();

/// Line width
void gl_line_width(float w);

/// Set viewport
void gl_set_viewport(const vec4i& v);

/// Set viewport
void gl_set_viewport(const vec2i& v);

// This is a public API. See above for documentation.
void gl_read_imagef(float* pixels, int w, int h, int nc);

// -----------------------------------------------------------------------------
// TEXTURE FUNCTIONS
// -----------------------------------------------------------------------------

/// Opengl texture object
struct gl_texture {
    // texture handle
    uint _tid = 0;
    // width
    int _width = 0;
    // height
    int _height = 0;
    // ncomp
    int _ncomp = 0;
    // floats
    bool _float = false;
    // srgb
    bool _srgb = true;
    // mipmap
    bool _mipmap = true;
};

// Implementation of make_texture.
void _init_texture(gl_texture& txt, int w, int h, int nc, const void* pixels,
    bool floats, bool linear, bool mipmap, bool as_float, bool as_srgb);

// Implementation of update_texture.
void _update_texture(
    gl_texture& txt, int w, int h, int nc, const void* pixels, bool floats);

/// Creates a texture with pixels values of size w, h with nc number of
/// components (1-4).
/// Internally use float if as_float and filtering if filter.
/// Returns the texture id.
inline gl_texture make_texture(int w, int h, int nc, const float* pixels,
    bool linear, bool mipmap, bool as_float) {
    auto txt = gl_texture();
    _init_texture(txt, w, h, nc, pixels, true, linear, mipmap, as_float, false);
    return txt;
}

/// Creates a texture with pixels values of size w, h with nc number of
/// components (1-4).
/// Internally use srgb lookup if as_srgb and filtering if filter.
/// Returns the texture id.
inline gl_texture make_texture(int w, int h, int nc,
    const unsigned char* pixels, bool linear, bool mipmap, bool as_srgb) {
    auto txt = gl_texture();
    _init_texture(txt, w, h, nc, pixels, false, linear, mipmap, false, as_srgb);
    return txt;
}

/// Creates a texture from an image.
/// Internally use float if as_float and filtering if filter.
/// Returns the texture id.
inline gl_texture make_texture(
    const image4f& img, bool linear, bool mipmap, bool as_float) {
    return make_texture(img.width(), img.height(), 4, (float*)img.data(),
        linear, mipmap, as_float);
}

/// Creates a texture from an image.
/// Internally use srgb lookup if as_srgb and filtering if filter.
/// Returns the texture id.
inline gl_texture make_texture(
    const image4b& img, bool linear, bool mipmap, bool as_srgb) {
    return make_texture(img.width(), img.height(), 4, (byte*)img.data(), linear,
        mipmap, as_srgb);
}

/// Updates the texture tid with new image data.
inline void update_texture(
    gl_texture& txt, int w, int h, int nc, const float* pixels) {
    _update_texture(txt, w, h, nc, pixels, true);
}

/// Updates the texture tid with new image data.
inline void update_texture(
    gl_texture& txt, int w, int h, int nc, const unsigned char* pixels) {
    _update_texture(txt, w, h, nc, pixels, false);
}

/// Updates the texture tid with new image data.
inline void update_texture(gl_texture& txt, const image4f& img) {
    update_texture(txt, img.width(), img.height(), 4, (const float*)img.data());
}

/// Updates the texture tid with new image data.
inline void update_texture(gl_texture& txt, const image4b& img) {
    update_texture(
        txt, img.width(), img.height(), 4, (const unsigned char*)img.data());
}

/// Binds a texture to a texture unit
void bind_texture(const gl_texture& txt, uint unit);

/// Unbinds
void unbind_texture(const gl_texture& txt, uint unit);

/// Get id
inline uint get_texture_id(const gl_texture& txt) { return txt._tid; }

/// Check if defined
inline bool is_texture_valid(const gl_texture& txt) { return (bool)txt._tid; }

/// Destroys the texture tid.
void clear_texture(gl_texture& txt);

/// Wrap values for texture
enum struct gl_texture_wrap {
    /// not set
    not_set = 0,
    /// repeat
    repeat = 1,
    /// clamp to edge
    clamp = 2,
    /// mirror
    mirror = 3,
};

/// Filter values for texture
enum struct gl_texture_filter {
    /// not set
    not_set = 0,
    /// linear
    linear = 1,
    /// nearest
    nearest = 2,
    /// mip-mapping
    linear_mipmap_linear = 3,
    /// mip-mapping
    nearest_mipmap_nearest = 4,
    /// mip-mapping
    linear_mipmap_nearest = 5,
    /// mip-mapping
    nearest_mipmap_linear = 6,
};

/// Texture information for parameter setting.
struct gl_texture_info {
    /// texture
    gl_texture txt = {};
    /// texture coordinate set
    int texcoord = 0;
    /// texture strength/scale (used by some models)
    float scale = 1;
    /// wrap mode
    gl_texture_wrap wrap_s = gl_texture_wrap::not_set;
    /// wrap mode
    gl_texture_wrap wrap_t = gl_texture_wrap::not_set;
    /// filter mode
    gl_texture_filter filter_mag = gl_texture_filter::not_set;
    /// filter mode
    gl_texture_filter filter_min = gl_texture_filter::not_set;

    /// default constructor
    gl_texture_info() {}
    /// constructor from texture id only
    gl_texture_info(const gl_texture& tid) : txt(tid) {}
};

// -----------------------------------------------------------------------------
// VERTEX ARRAY BUFFER
// -----------------------------------------------------------------------------

/// OpenGL vertex/element buffer
struct gl_vertex_buffer {
    // buffer id
    uint _bid = 0;
    // number of elements
    int _num = 0;
    // number of components
    int _ncomp = 0;
    // whether is is floats
    bool _float = true;
};

// Creates a buffer with num elements of size size stored in values, where
// content is dyanamic if dynamic.
void _init_vertex_buffer(gl_vertex_buffer& buf, int n, int nc,
    const void* values, bool as_float, bool dynamic);

// Updates the buffer bid with new data.
void _update_vertex_buffer(
    gl_vertex_buffer& buf, int n, int nc, const void* values, bool as_float);

/// Creates a buffer.
inline gl_vertex_buffer make_vertex_buffer(
    int num, int ncomp, const float* values, bool dynamic = false) {
    auto buf = gl_vertex_buffer();
    _init_vertex_buffer(buf, num, ncomp, values, true, dynamic);
    return buf;
}

/// Creates a buffer.
inline gl_vertex_buffer make_vertex_buffer(
    int num, int ncomp, const int* values, bool dynamic = false) {
    auto buf = gl_vertex_buffer();
    _init_vertex_buffer(buf, num, ncomp, values, true, dynamic);
    return buf;
}

/// Creates a buffer.
inline gl_vertex_buffer make_vertex_buffer(
    const vector<float>& values, bool dynamic = false) {
    return make_vertex_buffer(values.size(), 1, values.data(), dynamic);
}

/// Creates a buffer.
inline gl_vertex_buffer make_vertex_buffer(
    const vector<vec2f>& values, bool dynamic = false) {
    return make_vertex_buffer(
        values.size(), 2, (const float*)values.data(), dynamic);
}

/// Creates a buffer.
inline gl_vertex_buffer make_vertex_buffer(
    const vector<vec3f>& values, bool dynamic = false) {
    return make_vertex_buffer(
        values.size(), 3, (const float*)values.data(), dynamic);
}

/// Creates a buffer.
inline gl_vertex_buffer make_vertex_buffer(
    const vector<vec4f>& values, bool dynamic = false) {
    return make_vertex_buffer(
        values.size(), 4, (const float*)values.data(), dynamic);
}

/// Creates a buffer.
inline gl_vertex_buffer make_vertex_buffer(
    const vector<int>& values, bool dynamic = false) {
    return make_vertex_buffer(values.size(), 1, values.data(), dynamic);
}

/// Creates a buffer.
inline gl_vertex_buffer make_vertex_buffer(
    const vector<vec2i>& values, bool dynamic = false) {
    return make_vertex_buffer(
        values.size(), 2, (const int*)values.data(), dynamic);
}
/// Creates a buffer.
inline gl_vertex_buffer make_vertex_buffer(
    const vector<vec3i>& values, bool dynamic = false) {
    return make_vertex_buffer(
        values.size(), 3, (const int*)values.data(), dynamic);
}

/// Creates a buffer.
inline gl_vertex_buffer make_vertex_buffer(
    const vector<vec4i>& values, bool dynamic = false) {
    return make_vertex_buffer(
        values.size(), 4, (const int*)values.data(), dynamic);
}

/// Updates the buffer with new data.
inline void update_vertex_buffer(
    gl_vertex_buffer& buf, int num, int ncomp, const float* values) {
    _update_vertex_buffer(buf, num, ncomp, values, true);
}

/// Updates the buffer with new data.
inline void update_vertex_buffer(
    gl_vertex_buffer& buf, int num, int ncomp, const int* values) {
    _update_vertex_buffer(buf, num, ncomp, values, false);
}

/// Updates the buffer bid with new data.
inline void update_vertex_buffer(
    gl_vertex_buffer& buf, const vector<float>& values) {
    update_vertex_buffer(buf, values.size(), 1, values.data());
}

/// Updates the buffer bid with new data.
inline void update_vertex_buffer(
    gl_vertex_buffer& buf, const vector<vec2f>& values) {
    update_vertex_buffer(buf, values.size(), 2, (const float*)values.data());
}

/// Updates the buffer bid with new data.
inline void update_vertex_buffer(
    gl_vertex_buffer& buf, const vector<vec3f>& values) {
    update_vertex_buffer(buf, values.size(), 3, (const float*)values.data());
}

/// Updates the buffer bid with new data.
inline void update_vertex_buffer(
    gl_vertex_buffer& buf, const vector<vec4f>& values) {
    update_vertex_buffer(buf, values.size(), 4, (const float*)values.data());
}

/// Updates the buffer bid with new data.
inline void update_vertex_buffer(
    gl_vertex_buffer& buf, const vector<int>& values) {
    update_vertex_buffer(buf, values.size(), 1, values.data());
}

/// Updates the buffer bid with new data.
inline void update_vertex_buffer(
    gl_vertex_buffer& buf, const vector<vec2i>& values) {
    update_vertex_buffer(buf, values.size(), 2, (const int*)values.data());
}

/// Updates the buffer bid with new data.
inline void update_vertex_buffer(
    gl_vertex_buffer& buf, const vector<vec3i>& values) {
    update_vertex_buffer(buf, values.size(), 3, (const int*)values.data());
}

/// Updates the buffer bid with new data.
inline void update_vertex_buffer(
    gl_vertex_buffer& buf, const vector<vec4i>& values) {
    update_vertex_buffer(buf, values.size(), 4, (const int*)values.data());
}

/// Bind the buffer at a particular attribute location
void bind_vertex_buffer(const gl_vertex_buffer& buf, uint vattr);

/// Unbind the buffer
void unbind_vertex_buffer(const gl_vertex_buffer& buf, uint vattr);

/// Unbind the buffer
void unbind_vertex_buffer(uint vattr);

/// Get id
inline uint get_vertex_buffer_id(const gl_vertex_buffer& buf) {
    return buf._bid;
}

/// Check if defined
inline bool is_vertex_buffer_valid(const gl_vertex_buffer& buf) {
    return (bool)buf._bid;
}

/// Destroys the buffer
void clear_vertex_buffer(gl_vertex_buffer& buf);

// -----------------------------------------------------------------------------
// VERTEX ELEMENTS BUFFER
// -----------------------------------------------------------------------------

/// OpenGL vertex/element buffer
struct gl_element_buffer {
    // buffer id
    uint _bid = 0;
    // number of elements
    int _num = 0;
    // number of components
    int _ncomp = 0;
};

// Creates a buffer with num elements of size size stored in values, where
// content is dyanamic if dynamic.
// Returns the buffer id.
void _init_element_buffer(
    gl_element_buffer& buf, int n, int nc, const int* values, bool dynamic);

// Updates the buffer bid with new data.
void _update_element_buffer(
    gl_element_buffer& buf, int n, int nc, const int* values);

/// Creates a buffer.
inline gl_element_buffer make_element_buffer(
    int num, int ncomp, const int* values, bool dynamic = false) {
    auto buf = gl_element_buffer();
    _init_element_buffer(buf, num, ncomp, values, dynamic);
    return buf;
}

/// Creates a buffer.
inline gl_element_buffer make_element_buffer(
    const vector<int>& values, bool dynamic = false) {
    return make_element_buffer(values.size(), 1, values.data(), dynamic);
}

/// Creates a buffer.
inline gl_element_buffer make_element_buffer(
    const vector<vec2i>& values, bool dynamic = false) {
    return make_element_buffer(
        values.size(), 2, (const int*)values.data(), dynamic);
}

/// Creates a buffer.
inline gl_element_buffer make_element_buffer(
    const vector<vec3i>& values, bool dynamic = false) {
    return make_element_buffer(
        values.size(), 3, (const int*)values.data(), dynamic);
}

/// Creates a buffer.
inline gl_element_buffer make_element_buffer(
    const vector<vec4i>& values, bool dynamic = false) {
    return make_element_buffer(
        values.size(), 4, (const int*)values.data(), dynamic);
}

/// Updates the buffer with new data.
inline void update_element_buffer(
    gl_element_buffer& buf, int num, int ncomp, const int* values) {
    _update_element_buffer(buf, num, ncomp, values);
}

/// Updates the buffer bid with new data.
inline void update_element_buffer(
    gl_element_buffer& buf, const vector<int>& values) {
    update_element_buffer(buf, values.size(), 1, values.data());
}

/// Updates the buffer bid with new data.
inline void update_element_buffer(
    gl_element_buffer& buf, const vector<vec2i>& values) {
    update_element_buffer(buf, values.size(), 2, (const int*)values.data());
}

/// Updates the buffer bid with new data.
inline void update_element_buffer(
    gl_element_buffer& buf, const vector<vec3i>& values) {
    update_element_buffer(buf, values.size(), 3, (const int*)values.data());
}

/// Updates the buffer bid with new data.
inline void update_element_buffer(
    gl_element_buffer& buf, const vector<vec4i>& values) {
    update_element_buffer(buf, values.size(), 4, (const int*)values.data());
}

/// Draws elements.
void draw_elems(const gl_element_buffer& buf);

/// Get id
inline uint get_element_buffer_id(const gl_element_buffer& buf) {
    return buf._bid;
}

/// Check if defined
inline bool is_element_buffer_valid(const gl_element_buffer& buf) {
    return (bool)buf._bid;
}

/// Destroys the buffer
void clear_element_buffer(gl_element_buffer& buf);

// -----------------------------------------------------------------------------
// PROGRAM FUNCTIONS
// -----------------------------------------------------------------------------

/// OpenGL program
struct gl_program {
    // program id
    uint _pid = 0;
    // vertex shader id
    uint _vid = 0;
    // fragment shader id
    uint _fid = 0;
    // vertex array object
    uint _vao = 0;
};

/// Creates and OpenGL program from vertex and fragment code. Returns the
/// program id. Optionally return vertex and fragment shader ids. A VAO is
/// created.
gl_program make_program(const string& vertex, const string& fragment);

/// Destroys the program pid and optionally the sahders vid and fid.
void clear_program(gl_program& prog);

/// Get uniform location (simple GL wrapper that avoids GL includes)
int get_program_uniform_location(const gl_program& prog, const string& name);

/// Get uniform location (simple GL wrapper that avoids GL includes)
int get_program_attrib_location(const gl_program& prog, const string& name);

/// Get the names of all uniforms
vector<pair<string, int>> get_program_uniforms_names(const gl_program& prog);

/// Get the names of all attributes
vector<pair<string, int>> get_program_attributes_names(const gl_program& prog);

/// Set uniform integer values val for program pid and variable loc.
/// The values have nc number of components (1-4) and count elements
/// (for arrays).
bool set_program_uniform(
    gl_program& prog, int pos, const int* val, int ncomp, int count);

/// Set uniform float values val for program pid and variable var.
/// The values have nc number of components (1-4) and count elements
/// (for arrays).
bool set_program_uniform(
    gl_program& prog, int pos, const float* val, int ncomp, int count);

/// Set uniform float values val for program pid and variable var.
inline bool set_program_uniform(gl_program& prog, int var, bool val) {
    auto vali = (int)val;
    return set_program_uniform(prog, var, &vali, 1, 1);
}

/// Set uniform float values val for program pid and variable var.
inline bool set_program_uniform(gl_program& prog, int var, int val) {
    return set_program_uniform(prog, var, &val, 1, 1);
}

/// Set uniform float values val for program pid and variable var.
inline bool set_program_uniform(gl_program& prog, int var, float val) {
    return set_program_uniform(prog, var, &val, 1, 1);
}

/// Set uniform float values val for program pid and variable var.
inline bool set_program_uniform(gl_program& prog, int var, const vec2f& val) {
    return set_program_uniform(prog, var, val.data(), 2, 1);
}
/// Set uniform float values val for program pid and variable var.
inline bool set_program_uniform(gl_program& prog, int var, const vec3f& val) {
    return set_program_uniform(prog, var, val.data(), 3, 1);
}
/// Set uniform float values val for program pid and variable var.
inline bool set_program_uniform(gl_program& prog, int var, const vec4f& val) {
    return set_program_uniform(prog, var, val.data(), 4, 1);
}

/// Set uniform float values val for program pid and variable var.
inline bool set_program_uniform(gl_program& prog, int var, const vec2i& val) {
    return set_program_uniform(prog, var, val.data(), 2, 1);
}
/// Set uniform float values val for program pid and variable var.
inline bool set_program_uniform(gl_program& prog, int var, const vec3i& val) {
    return set_program_uniform(prog, var, val.data(), 3, 1);
}
/// Set uniform float values val for program pid and variable var.
inline bool set_program_uniform(gl_program& prog, int var, const vec4i& val) {
    return set_program_uniform(prog, var, val.data(), 4, 1);
}

/// Set uniform float values val for program pid and variable var.
inline bool set_program_uniform(gl_program& prog, int var, const mat4f& val) {
    return set_program_uniform(prog, var, (float*)val.data(), 16, 1);
}

/// Set uniform float values val for program pid and variable var.
inline bool set_program_uniform(gl_program& prog, int var, const frame3f& val) {
    return set_program_uniform(prog, var, (float*)val.data(), 12, 1);
}

/// Set uniform float values val for program pid and variable var.
inline bool set_program_uniform(
    gl_program& prog, int var, const int* val, int num) {
    return set_program_uniform(prog, var, val, 1, num);
}

/// Set uniform float values val for program pid and variable var.
inline bool set_program_uniform(
    gl_program& prog, int var, const float* val, int num) {
    return set_program_uniform(prog, var, val, 1, num);
}

/// Set uniform float values val for program pid and variable var.
inline bool set_program_uniform(
    gl_program& prog, int var, const vec2f* val, int num) {
    return set_program_uniform(prog, var, (const float*)val, 2, num);
}
/// Set uniform float values val for program pid and variable var.
inline bool set_program_uniform(
    gl_program& prog, int var, const vec3f* val, int num) {
    return set_program_uniform(prog, var, (const float*)val, 3, num);
}
/// Set uniform float values val for program pid and variable var.
inline bool set_program_uniform(
    gl_program& prog, int var, const vec4f* val, int num) {
    return set_program_uniform(prog, var, (const float*)val, 4, num);
}

/// Set uniform float values val for program pid and variable var.
inline bool set_program_uniform(
    gl_program& prog, int var, const vec2i* val, int num) {
    return set_program_uniform(prog, var, (const int*)val, 2, num);
}
/// Set uniform float values val for program pid and variable var.
inline bool set_program_uniform(
    gl_program& prog, int var, const vec3i* val, int num) {
    return set_program_uniform(prog, var, (const int*)val, 3, num);
}
/// Set uniform float values val for program pid and variable var.
inline bool set_program_uniform(
    gl_program& prog, int var, const vec4i* val, int num) {
    return set_program_uniform(prog, var, (const int*)val, 4, num);
}

/// Set uniform float values val for program pid and variable var.
inline bool set_program_uniform(
    gl_program& prog, int var, const mat4f* val, int num) {
    return set_program_uniform(prog, var, (float*)val, 16, num);
}

/// Set uniform float values val for program pid and variable var.
template <typename T>
inline bool set_program_uniform(
    gl_program& prog, const string& var, const T& val) {
    auto loc = get_program_uniform_location(prog, var);
    if (loc < 0) return false;
    return set_program_uniform(prog, loc, val);
}

/// Set uniform float values val for program pid and variable var.
template <typename T>
inline bool set_program_uniform(
    gl_program& prog, const string& var, const T* val, int num) {
    auto loc = get_program_uniform_location(prog, var);
    if (loc < 0) return false;
    return set_program_uniform(loc, val, num);
}

/// Set uniform texture id tid and unit tunit for program pid and variable
/// var.
bool set_program_uniform_texture(
    gl_program& prog, int pos, const gl_texture_info& tinfo, uint tunit);

/// Set uniform texture id tid and unit tunit for program pid and variable
/// var. Optionally sets the int variable varon to 0/1 whether the texture
/// is enable on not.
inline bool set_program_uniform_texture(gl_program& prog, int var, int varon,
    const gl_texture_info& tinfo, uint tunit) {
    if (!set_program_uniform_texture(prog, var, tinfo, tunit)) return false;
    if (!set_program_uniform(prog, varon, is_texture_valid(tinfo.txt)))
        return false;
    return true;
}

/// Set uniform texture id tid and unit tunit for program pid and variable
/// var.
inline bool set_program_uniform_texture(gl_program& prog, const string& var,
    const gl_texture_info& tinfo, uint tunit) {
    auto loc = get_program_uniform_location(prog, var);
    if (loc < 0) return false;
    return set_program_uniform_texture(prog, loc, tinfo, tunit);
}

/// Set uniform texture id tid and unit tunit for program pid and variable
/// var. Optionally sets the int variable varon to 0/1 whether the texture
/// is enable on not.
inline bool set_program_uniform_texture(gl_program& prog, const string& var,
    const string& varon, const gl_texture_info& tinfo, uint tunit) {
    auto loc = get_program_uniform_location(prog, var);
    if (loc < 0) return false;
    auto locon = get_program_uniform_location(prog, varon);
    if (locon < 0) return false;
    return set_program_uniform_texture(prog, loc, locon, tinfo, tunit);
}

/// Sets a constant value for a vertex attribute for program pid and
/// variable var. The attribute has nc components.
bool set_program_vertattr(
    gl_program& prog, int pos, const float* value, int nc);

/// Sets a constant value for a vertex attribute for program pid and
/// variable var. The attribute has nc components.
bool set_program_vertattr(gl_program& prog, int pos, const int* value, int nc);

/// Sets a vartex attribute for program pid and variable var to the buffer
/// bid. The attribute has nc components and per-vertex values values.
bool set_program_vertattr(
    gl_program& prog, const string& var, const gl_vertex_buffer& buf);

/// Sets a vartex attribute for program pid and variable var. The attribute
/// has nc components and either buffer bid or a single value def
/// (if bid is zero). Convenience wrapper to above functions.
bool set_program_vertattr(gl_program& prog, int pos,
    const gl_vertex_buffer& buf, int nc, const float* def);

/// Sets a vartex attribute for program pid and variable var. The attribute
/// is either a buffer bid or a single value def
/// (if bid is zero). Convenience wrapper to above functions.
inline bool set_program_vertattr(
    gl_program& prog, int var, const gl_vertex_buffer& buf, const vec2f& def) {
    return set_program_vertattr(prog, var, buf, 2, def.data());
}
/// Sets a vartex attribute for program pid and variable var. The attribute
/// is either a buffer bid or a single value def
/// (if bid is zero). Convenience wrapper to above functions.
inline bool set_program_vertattr(
    gl_program& prog, int var, const gl_vertex_buffer& buf, const vec3f& def) {
    return set_program_vertattr(prog, var, buf, 3, def.data());
}
/// Sets a vartex attribute for program pid and variable var. The attribute
/// is either a buffer bid or a single value def
/// (if bid is zero). Convenience wrapper to above functions.
inline bool set_program_vertattr(
    gl_program& prog, int var, const gl_vertex_buffer& buf, const vec4f& def) {
    return set_program_vertattr(prog, var, buf, 4, def.data());
}

/// Sets a vartex attribute for program pid and variable var. The attribute
/// is either a buffer bid or a single value def
/// (if bid is zero). Convenience wrapper to above functions.
inline bool set_program_vertattr(gl_program& prog, const string& var,
    const gl_vertex_buffer& buf, const vec2f& def) {
    auto loc = get_program_attrib_location(prog, var);
    if (loc < 0) return false;
    return set_program_vertattr(prog, loc, buf, def);
}
/// Sets a vartex attribute for program pid and variable var. The attribute
/// is either a buffer bid or a single value def
/// (if bid is zero). Convenience wrapper to above functions.
inline bool set_program_vertattr(gl_program& prog, const string& var,
    const gl_vertex_buffer& buf, const vec3f& def) {
    auto loc = get_program_attrib_location(prog, var);
    if (loc < 0) return false;
    return set_program_vertattr(prog, loc, buf, def);
}
/// Sets a vartex attribute for program pid and variable var. The attribute
/// is either a buffer bid or a single value def
/// (if bid is zero). Convenience wrapper to above functions.
inline bool set_program_vertattr(gl_program& prog, const string& var,
    const gl_vertex_buffer& buf, const vec4f& def) {
    auto loc = get_program_attrib_location(prog, var);
    if (loc < 0) return false;
    return set_program_vertattr(prog, loc, buf, def);
}

/// Check whether it is valid
inline bool is_program_valid(const gl_program& prog) { return (bool)prog._pid; }

/// Binds a program
void bind_program(const gl_program& prog);

/// Unbind a program
void unbind_program(const gl_program& prog);

// -----------------------------------------------------------------------------
// IMAGE SHADER FUNCTIONS
// -----------------------------------------------------------------------------

/// A shader for displaying images
struct gl_stdimage_program {
    // program
    gl_program _prog = {};
    // vertex array
    gl_vertex_buffer _vbo = {};
    // element array
    gl_element_buffer _ebo = {};
};

/// Initialize the program. Call with true only after the GL is initialized.
gl_stdimage_program make_stdimage_program();

/// As above but includes an exposure/gamma correction.
inline void draw_image(gl_stdimage_program& prog, const gl_texture& txt,
    const vec2i& win_size, const vec2f& offset, float zoom, float exposure,
    float gamma, bool filmic) {
    assert(is_texture_valid(txt));

    bind_program(prog._prog);

    gl_enable_blending(true);
    gl_set_blend_over();

    bind_texture(txt, 0);
    set_program_uniform(prog._prog, "zoom", zoom);
    set_program_uniform(
        prog._prog, "win_size", vec2f{(float)win_size.x, (float)win_size.y});
    set_program_uniform(prog._prog, "offset", offset);
    set_program_uniform(prog._prog, "tonemap.filmic", filmic);
    set_program_uniform(prog._prog, "tonemap.exposure", exposure);
    set_program_uniform(prog._prog, "tonemap.gamma", gamma);
    set_program_uniform_texture(prog._prog, "img", txt, 0);

    set_program_vertattr(prog._prog, "vert_texcoord", prog._vbo, vec2f{0, 0});
    draw_elems(prog._ebo);

    unbind_program(prog._prog);

    gl_enable_blending(false);

    assert(gl_check_error());
}

/// Draw an texture tid of size img_w, img_h on a window of size win_w,
/// win_h with top-left corner at ox, oy with a zoom zoom.
inline void draw_image(gl_stdimage_program& prog, const gl_texture& txt,
    const vec2i& win_size, const vec2f& offset, float zoom) {
    assert(gl_check_error());
    draw_image(prog, txt, win_size, offset, zoom, 0, 1, false);
    assert(gl_check_error());
}

// -----------------------------------------------------------------------------
// STANDARD SHADER FUNCTIONS
// -----------------------------------------------------------------------------

/// Shade with a physically-based standard shader based on Phong/GGX.
/// Filmic tone mapping from
/// https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
struct gl_stdsurface_program {
    // gl program
    gl_program _prog = {};
};

/// Initialize a standard shader. Call with true only after the gl has
/// been initialized
gl_stdsurface_program make_stdsurface_program();

/// Check if the program is valid
inline bool is_program_valid(const gl_stdsurface_program& prog) {
    return is_program_valid(prog._prog);
}

/// Starts a frame by setting exposure/gamma values, camera transforms and
/// projection. Sets also whether to use full shading or a quick eyelight
/// preview.
inline void begin_stdsurface_frame(gl_stdsurface_program& prog,
    bool shade_eyelight, float tonemap_exposure, float tonemap_gamma,
    bool tonemap_filmic, const mat4f& camera_xform,
    const mat4f& camera_xform_inv, const mat4f& camera_proj) {
    static auto eyelight_id =
        get_program_uniform_location(prog._prog, "lighting.eyelight");
    static auto exposure_id =
        get_program_uniform_location(prog._prog, "tonemap.exposure");
    static auto gamma_id =
        get_program_uniform_location(prog._prog, "tonemap.gamma");
    static auto filmic_id =
        get_program_uniform_location(prog._prog, "tonemap.filmic");
    static auto xform_id =
        get_program_uniform_location(prog._prog, "camera.xform");
    static auto xform_inv_id =
        get_program_uniform_location(prog._prog, "camera.xform_inv");
    static auto proj_id =
        get_program_uniform_location(prog._prog, "camera.proj");
    assert(gl_check_error());
    bind_program(prog._prog);
    set_program_uniform(prog._prog, eyelight_id, shade_eyelight);
    set_program_uniform(prog._prog, exposure_id, tonemap_exposure);
    set_program_uniform(prog._prog, gamma_id, tonemap_gamma);
    set_program_uniform(prog._prog, filmic_id, tonemap_filmic);
    set_program_uniform(prog._prog, xform_id, camera_xform);
    set_program_uniform(prog._prog, xform_inv_id, camera_xform_inv);
    set_program_uniform(prog._prog, proj_id, camera_proj);
    assert(gl_check_error());
}

/// Ends a frame.
inline void end_stdsurface_frame(gl_stdsurface_program& prog) {
    assert(gl_check_error());
    unbind_program(prog._prog);
    //    glBindVertexArray(0);
    //    glUseProgram(0);
    assert(gl_check_error());
}

/// Set num lights with position pos, color ke, type ltype. Also set the
/// ambient illumination amb.
inline void set_stdsurface_lights(gl_stdsurface_program& prog, const vec3f& amb,
    int num, vec3f* pos, vec3f* ke, gl_ltype* type) {
    static auto amb_id =
        get_program_uniform_location(prog._prog, "lighting.amb");
    static auto lnum_id =
        get_program_uniform_location(prog._prog, "lighting.lnum");
    static auto lpos_id =
        get_program_uniform_location(prog._prog, "lighting.lpos");
    static auto lke_id =
        get_program_uniform_location(prog._prog, "lighting.lke");
    static auto ltype_id =
        get_program_uniform_location(prog._prog, "lighting.ltype");
    assert(gl_check_error());
    set_program_uniform(prog._prog, amb_id, amb);
    set_program_uniform(prog._prog, lnum_id, num);
    set_program_uniform(prog._prog, lpos_id, pos, num);
    set_program_uniform(prog._prog, lke_id, ke, num);
    set_program_uniform(prog._prog, ltype_id, (int*)type, num);
    assert(gl_check_error());
}

/// Begins drawing a shape with transform xform.
inline void begin_stdsurface_shape(
    gl_stdsurface_program& prog, const mat4f& xform) {
    static auto xform_id =
        get_program_uniform_location(prog._prog, "shape_xform");
    assert(gl_check_error());
    set_program_uniform(prog._prog, xform_id, xform);
    assert(gl_check_error());
}

/// End shade drawing.
inline void end_stdsurface_shape(gl_stdsurface_program& prog) {
    assert(gl_check_error());
    for (int i = 0; i < 16; i++) unbind_vertex_buffer(i);
    assert(gl_check_error());
}

/// Set the object as highlighted.
inline void set_stdsurface_highlight(
    gl_stdsurface_program& prog, const vec4f& highlight) {
    static auto highlight_id =
        get_program_uniform_location(prog._prog, "highlight");
    set_program_uniform(prog._prog, highlight_id, highlight);
}

/// Set material values with emission ke, diffuse kd, specular ks and
/// specular roughness rs, opacity op. Indicates textures ids with the
/// correspoinding XXX_txt variables. Sets also normal and occlusion
/// maps. Works for points/lines/triangles (diffuse for points,
/// Kajiya-Kay for lines, GGX/Phong for triangles).
/// Material type matches the scene material type.
inline void set_stdsurface_material(gl_stdsurface_program& prog,
    material_type mtype, gl_etype etype, const vec3f& ke, const vec3f& kd,
    const vec3f& ks, float rs, float op, const gl_texture_info& ke_txt,
    const gl_texture_info& kd_txt, const gl_texture_info& ks_txt,
    const gl_texture_info& rs_txt, const gl_texture_info& norm_txt,
    const gl_texture_info& occ_txt, bool use_phong, bool double_sided,
    bool alpha_cutout) {
    static auto mtype_id =
        get_program_uniform_location(prog._prog, "material.mtype");
    static auto etype_id =
        get_program_uniform_location(prog._prog, "material.etype");
    static auto ke_id = get_program_uniform_location(prog._prog, "material.ke");
    static auto kd_id = get_program_uniform_location(prog._prog, "material.kd");
    static auto ks_id = get_program_uniform_location(prog._prog, "material.ks");
    static auto rs_id = get_program_uniform_location(prog._prog, "material.rs");
    static auto op_id = get_program_uniform_location(prog._prog, "material.op");
    static auto ke_txt_id =
        get_program_uniform_location(prog._prog, "material.txt_ke");
    static auto ke_txt_on_id =
        get_program_uniform_location(prog._prog, "material.txt_ke_on");
    static auto kd_txt_id =
        get_program_uniform_location(prog._prog, "material.txt_kd");
    static auto kd_txt_on_id =
        get_program_uniform_location(prog._prog, "material.txt_kd_on");
    static auto ks_txt_id =
        get_program_uniform_location(prog._prog, "material.txt_ks");
    static auto ks_txt_on_id =
        get_program_uniform_location(prog._prog, "material.txt_ks_on");
    static auto rs_txt_id =
        get_program_uniform_location(prog._prog, "material.txt_rs");
    static auto rs_txt_on_id =
        get_program_uniform_location(prog._prog, "material.txt_rs_on");
    static auto norm_txt_id =
        get_program_uniform_location(prog._prog, "material.txt_norm");
    static auto norm_txt_on_id =
        get_program_uniform_location(prog._prog, "material.txt_norm_on");
    static auto occ_txt_id =
        get_program_uniform_location(prog._prog, "material.txt_occ");
    static auto occ_txt_on_id =
        get_program_uniform_location(prog._prog, "material.txt_occ_on");
    static auto norm_scale_id =
        get_program_uniform_location(prog._prog, "material.norm_scale");
    static auto occ_scale_id =
        get_program_uniform_location(prog._prog, "material.occ_scale");
    static auto use_phong_id =
        get_program_uniform_location(prog._prog, "material.use_phong");
    static auto double_sided_id =
        get_program_uniform_location(prog._prog, "material.double_sided");
    static auto alpha_cutout_id =
        get_program_uniform_location(prog._prog, "material.alpha_cutout");

    static auto mtypes = unordered_map<material_type, int>{
        {material_type::specular_roughness, 1},
        {material_type::metallic_roughness, 2},
        {material_type::specular_glossiness, 3}};

    assert(gl_check_error());
    set_program_uniform(prog._prog, mtype_id, mtypes.at(mtype));
    set_program_uniform(prog._prog, etype_id, (int)etype);
    set_program_uniform(prog._prog, ke_id, ke);
    set_program_uniform(prog._prog, kd_id, kd);
    set_program_uniform(prog._prog, ks_id, ks);
    set_program_uniform(prog._prog, rs_id, rs);
    set_program_uniform(prog._prog, op_id, op);
    set_program_uniform_texture(prog._prog, ke_txt_id, ke_txt_on_id, ke_txt, 0);
    set_program_uniform_texture(prog._prog, kd_txt_id, kd_txt_on_id, kd_txt, 1);
    set_program_uniform_texture(prog._prog, ks_txt_id, ks_txt_on_id, ks_txt, 2);
    set_program_uniform_texture(prog._prog, rs_txt_id, rs_txt_on_id, rs_txt, 3);
    set_program_uniform_texture(
        prog._prog, norm_txt_id, norm_txt_on_id, norm_txt, 4);
    set_program_uniform_texture(
        prog._prog, occ_txt_id, occ_txt_on_id, occ_txt, 5);
    set_program_uniform(prog._prog, norm_scale_id, norm_txt.scale);
    set_program_uniform(prog._prog, occ_scale_id, occ_txt.scale);
    set_program_uniform(prog._prog, use_phong_id, use_phong);
    set_program_uniform(prog._prog, double_sided_id, double_sided);
    set_program_uniform(prog._prog, alpha_cutout_id, alpha_cutout);
    assert(gl_check_error());
}

/// Set vertex data with buffers for position pos, normals norm, texture
/// coordinates texcoord, per-vertex color color and tangent space tangsp.
inline void set_stdsurface_vert(gl_stdsurface_program& prog,
    const gl_vertex_buffer& pos, const gl_vertex_buffer& norm,
    const gl_vertex_buffer& texcoord, const gl_vertex_buffer& color,
    const gl_vertex_buffer& tangsp) {
    static auto pos_id = get_program_attrib_location(prog._prog, "vert_pos");
    static auto norm_id = get_program_attrib_location(prog._prog, "vert_norm");
    static auto texcoord_id =
        get_program_attrib_location(prog._prog, "vert_texcoord");
    static auto color_id =
        get_program_attrib_location(prog._prog, "vert_color");
    static auto tangsp_id =
        get_program_attrib_location(prog._prog, "vert_tangsp");
    assert(gl_check_error());
    set_program_vertattr(prog._prog, pos_id, pos, zero3f);
    set_program_vertattr(prog._prog, norm_id, norm, zero3f);
    set_program_vertattr(prog._prog, texcoord_id, texcoord, zero2f);
    set_program_vertattr(prog._prog, color_id, color, vec4f{1, 1, 1, 1});
    set_program_vertattr(prog._prog, tangsp_id, tangsp, zero4f);
    assert(gl_check_error());
}

/// Set vertex data with buffers for skinning.
inline void set_stdsurface_vert_skinning(gl_stdsurface_program& prog,
    const gl_vertex_buffer& weights, const gl_vertex_buffer& joints,
    int nxforms, const mat4f* xforms) {
    static auto type_id = get_program_uniform_location(prog._prog, "skin_type");
    static auto xforms_id =
        get_program_uniform_location(prog._prog, "skin_xforms");
    static auto weights_id =
        get_program_attrib_location(prog._prog, "vert_skin_weights");
    static auto joints_id =
        get_program_attrib_location(prog._prog, "vert_skin_joints");
    int type = 1;
    set_program_uniform(prog._prog, type_id, type);
    set_program_uniform(prog._prog, xforms_id, xforms, min(nxforms, 32));
    set_program_vertattr(prog._prog, weights_id, weights, zero4f);
    set_program_vertattr(prog._prog, joints_id, joints, zero4f);
}

/// Set vertex data with buffers for skinning.
inline void set_stdsurface_vert_gltf_skinning(gl_stdsurface_program& prog,
    const gl_vertex_buffer& weights, const gl_vertex_buffer& joints,
    int nxforms, const mat4f* xforms) {
    static auto type_id = get_program_uniform_location(prog._prog, "skin_type");
    static auto xforms_id =
        get_program_uniform_location(prog._prog, "skin_xforms");
    static auto weights_id =
        get_program_attrib_location(prog._prog, "vert_skin_weights");
    static auto joints_id =
        get_program_attrib_location(prog._prog, "vert_skin_joints");
    int type = 2;
    set_program_uniform(prog._prog, type_id, type);
    set_program_uniform(prog._prog, xforms_id, xforms, min(nxforms, 32));
    set_program_vertattr(prog._prog, weights_id, weights, zero4f);
    set_program_vertattr(prog._prog, joints_id, joints, zero4f);
}

/// Disables vertex skinning.
inline void set_stdsurface_vert_skinning_off(gl_stdsurface_program& prog) {
    static auto type_id = get_program_uniform_location(prog._prog, "skin_type");
    // static auto xforms_id = get_program_uniform_location(prog._prog,
    // "skin_xforms");
    static auto weights_id =
        get_program_attrib_location(prog._prog, "vert_skin_weights");
    static auto joints_id =
        get_program_attrib_location(prog._prog, "vert_skin_joints");
    int type = 0;
    set_program_uniform(prog._prog, type_id, type);
    set_program_vertattr(prog._prog, weights_id, {}, zero4f);
    set_program_vertattr(prog._prog, joints_id, {}, zero4f);
}

}  // namespace ygl

// Forward declaration
struct GLFWwindow;

// -----------------------------------------------------------------------------
// OPENGL WINDOWS AND WIDGETS
// -----------------------------------------------------------------------------
namespace ygl {

// Forward declaration
struct gl_window;

/// Text callback
typedef void (*gl_text_callback)(gl_window*, unsigned int);

/// Mouse callback
typedef void (*gl_mouse_callback)(gl_window*, int button, bool press, int mods);

/// Window refresh callback
typedef void (*gl_refresh_callback)(gl_window*);

/// Window
struct gl_window {
    GLFWwindow* _gwin = nullptr;
    void* _user_pointer = nullptr;
    int _widget_width = 320;
    bool _widget_enabled = false;
    gl_text_callback _text_cb = nullptr;
    gl_mouse_callback _mouse_cb = nullptr;
    gl_refresh_callback _refresh_cb = nullptr;
};

/// Initialize gl_window
gl_window* make_window(
    int width, int height, const string& title, void* user_pointer = nullptr);

/// Set gl_window callbacks
void set_window_callbacks(gl_window* win, gl_text_callback text_cb,
    gl_mouse_callback mouse_cb, gl_refresh_callback refresh_cb);

/// Clear gl_window
void clear_window(gl_window* win);

/// Gets the user poiner
inline void* get_user_pointer(gl_window* win) { return win->_user_pointer; }

/// Set gl_window title
void set_window_title(gl_window* win, const string& title);

/// Wait events
void wait_events(gl_window* win);

/// Poll events
void poll_events(gl_window* win);

/// Swap buffers
void swap_buffers(gl_window* win);

/// Should close
bool should_close(gl_window* win);

/// Mouse button
int get_mouse_button(gl_window* win);

/// Mouse position
vec2i get_mouse_pos(gl_window* win);

/// Mouse position
vec2f get_mouse_posf(gl_window* win);

/// Window size
vec2i get_window_size(gl_window* win);

/// Check if a key is pressed (not all keys are supported)
bool get_key(gl_window* win, int key);

/// Framebuffer size
vec2i get_framebuffer_size(gl_window* win);

/// Widgets
inline int get_widget_size(gl_window* win) { return win->_widget_width; }

/// Read pixels
vector<vec4b> get_screenshot(
    gl_window* win, vec2i& wh, bool flipy = true, bool back = false);

/// Save a screenshot to disk
inline void save_screenshot(gl_window* win, const string& imfilename) {
    auto wh = vec2i{0, 0};
    auto pixels = get_screenshot(win, wh);
    save_image(imfilename, wh.x, wh.y, 4, (unsigned char*)pixels.data());
}

/// Initialize widgets
void init_widgets(gl_window* win);

/// Begin draw widget
bool begin_widgets(gl_window* win, const string& title);

/// End draw widget
void end_widgets(gl_window* win);

/// Whether widget are active
bool get_widget_active(gl_window* win);

/// Horizontal separator
void draw_separator_widget(gl_window* win);

/// Indent widget
void draw_indent_widget_begin(gl_window* win);

/// Indent widget
void draw_indent_widget_end(gl_window* win);

/// Continue line with next widget
void draw_continue_widget(gl_window* win);

/// Label widget
void draw_label_widget(gl_window* win, const string& lbl, const string& msg);

/// Label widget
template <typename... Args>
inline void draw_label_widget(
    gl_window* win, const string& lbl, const string& fmt, const Args&... args);

/// Label widget
template <typename T>
inline void draw_label_widget(gl_window* win, const string& lbl, const T& val) {
    auto sst = stringstream();
    sst << val;
    return draw_label_widget(win, lbl, sst.str());
}

/// Value widget
bool draw_value_widget(gl_window* win, const string& lbl, string& str);

/// Value widget
bool draw_value_widget(gl_window* win, const string& lbl, int* val, int ncomp,
    int min = 0, int max = 1, int incr = 1);

/// Value widget
bool draw_value_widget(gl_window* win, const string& lbl, float* val, int ncomp,
    float min = 0, float max = 1, float incr = 1);

/// Value widget
inline bool draw_value_widget(gl_window* win, const string& lbl, int& val,
    int min = 0, int max = 1, int incr = 1) {
    return draw_value_widget(win, lbl, &val, 1, min, max, incr);
}

/// Value widget
inline bool draw_value_widget(gl_window* win, const string& lbl, vec2i& val,
    int min = 0, int max = 1, int incr = 1) {
    return draw_value_widget(win, lbl, val.data(), 2, min, max, incr);
}
/// Value widget
inline bool draw_value_widget(gl_window* win, const string& lbl, vec3i& val,
    int min = 0, int max = 1, int incr = 1) {
    return draw_value_widget(win, lbl, val.data(), 3, min, max, incr);
}
/// Value widget
inline bool draw_value_widget(gl_window* win, const string& lbl, vec4i& val,
    int min = 0, int max = 1, int incr = 1) {
    return draw_value_widget(win, lbl, val.data(), 4, min, max, incr);
}

/// Value widget
inline bool draw_value_widget(gl_window* win, const string& lbl, float& val,
    float min = 0, float max = 1, float incr = 1) {
    return draw_value_widget(win, lbl, &val, 1, min, max, incr);
}

/// Value widget
inline bool draw_value_widget(gl_window* win, const string& lbl, vec2f& val,
    float min = 0, float max = 1, float incr = 1) {
    return draw_value_widget(win, lbl, val.data(), 2, min, max, incr);
}
/// Value widget
inline bool draw_value_widget(gl_window* win, const string& lbl, vec3f& val,
    float min = 0, float max = 1, float incr = 1) {
    return draw_value_widget(win, lbl, val.data(), 3, min, max, incr);
}
/// Value widget
inline bool draw_value_widget(gl_window* win, const string& lbl, vec4f& val,
    float min = 0, float max = 1, float incr = 1) {
    return draw_value_widget(win, lbl, val.data(), 4, min, max, incr);
}

/// Slider widget
inline bool draw_value_widget(gl_window* win, const string& lbl, mat4f& val,
    float min = 0, float max = 1, float incr = 1) {
    auto modx = draw_value_widget(win, lbl + ".x", val.x, min, max, incr);
    auto mody = draw_value_widget(win, lbl + ".y", val.y, min, max, incr);
    auto modz = draw_value_widget(win, lbl + ".z", val.z, min, max, incr);
    auto modw = draw_value_widget(win, lbl + ".w", val.w, min, max, incr);
    return modx || mody || modz || modw;
}

/// Slider widget
inline bool draw_value_widget(gl_window* win, const string& lbl, frame3f& val,
    float min = -1, float max = 1, float incr = 1) {
    auto modx = draw_value_widget(win, lbl + ".x", val.x, -1, 1, 0.01f);
    auto mody = draw_value_widget(win, lbl + ".y", val.y, -1, 1, 0.01f);
    auto modz = draw_value_widget(win, lbl + ".z", val.z, -1, 1, 0.01f);
    auto modo = draw_value_widget(win, lbl + ".o", val.o, min, max, incr);
    // TODO: orthonormalize
    return modx || mody || modz || modo;
}

/// Slider widget
inline bool draw_value_widget(
    gl_window* win, const string& lbl, quat4f& val, float incr = 1) {
    auto mod = draw_value_widget(win, lbl, *(vec4f*)&val, -1, 1, incr);
    if (mod) val = normalize(val);
    return mod;
}

/// Color widget
bool draw_color_widget(gl_window* win, const string& lbl, vec4f& val);

/// Color widget
bool draw_color_widget(gl_window* win, const string& lbl, vec4b& val);

/// Color widget
bool draw_color_widget(gl_window* win, const string& lbl, vec3f& val);

// Support
inline bool _enum_widget_labels_ptr(void* data, int idx, const char** out) {
    auto labels = (vector<pair<string, int>>*)data;
    *out = labels->at(idx).first.c_str();
    return true;
}

// Support
inline bool _enum_widget_labels_int(void* data, int idx, const char** out) {
    auto labels = (vector<pair<string, int>>*)data;
    *out = labels->at(idx).first.c_str();
    return true;
}

/// Enum widget
bool draw_value_widget(gl_window* win, const string& lbl, int& val,
    const vector<pair<string, int>>& labels);

/// Enum widget
template <typename T>
inline bool draw_value_widget(gl_window* win, const string& lbl, T& val,
    const vector<pair<string, T>>& labels) {
    return draw_value_widget(
        win, lbl, (int&)val, (const vector<pair<string, int>>&)labels);
}

/// Bool widget
bool draw_value_widget(gl_window* win, const string& lbl, bool& val);

/// Button widget
bool draw_button_widget(gl_window* win, const string& lbl);

/// Collapsible header
bool draw_header_widget(gl_window* win, const string& lbl);

/// Start tree node
bool draw_tree_widget_begin(gl_window* win, const string& lbl);

/// Collapsible header
void draw_tree_widget_end(gl_window* win);

/// Start selectable tree node
bool draw_tree_widget_begin(
    gl_window* win, const string& lbl, void*& selection, void* content);

/// Start selectable tree node
bool draw_tree_widget_begin(gl_window* win, const string& lbl, void*& selection,
    void* content, const vec4f& col);

/// End selectable tree node
void draw_tree_widget_end(gl_window* win, void* content);

/// Selectable tree leaf node
void draw_tree_widget_leaf(
    gl_window* win, const string& lbl, void*& selection, void* content);

/// Selectable tree leaf node
void draw_tree_widget_leaf(gl_window* win, const string& lbl, void*& selection,
    void* content, const vec4f& col);

/// Image widget
void draw_image_widget(
    gl_window* win, int tid, const vec2i& size, const vec2i& imsize);

/// Scroll region
void draw_scroll_widget_begin(
    gl_window* win, const string& lbl, int height, bool border);

/// Scroll region
void draw_scroll_widget_end(gl_window* win);

/// Scroll region
void draw_scroll_widget_here(gl_window* win);

/// Group ids
void draw_groupid_widget_begin(gl_window* win, int gid);

/// Group ids
void draw_groupid_widget_begin(gl_window* win, void* gid);

/// Group ids
void draw_groupid_widget_end(gl_window* win);

/// Text color
void draw_tree_widget_color_begin(gl_window* win, const vec4f& color);

/// Text color
void draw_tree_widget_color_end(gl_window* win);

/// Tonemapping widgets
inline void draw_tonemap_widgets(gl_window* win, const string& lbl,
    float& exposure, float& gamma, bool& filmic) {
    draw_value_widget(win, lbl + "exposure", exposure, -20, 20, 1);
    draw_value_widget(win, lbl + "gamma", gamma, 0.1, 5, 0.1);
    draw_value_widget(win, lbl + "filmic", filmic);
}

/// Draws a widget that can selected the camera
inline bool draw_camera_widget(
    gl_window* win, const string& lbl, scene* scn, camera*& cam) {
    auto camera_names = vector<pair<string, camera*>>{};
    for (auto cam : scn->cameras) camera_names.push_back({cam->name, cam});
    return draw_value_widget(win, lbl, cam, camera_names);
}

/// Draws widgets for a whole scene. Used for quickly making demos.
bool draw_scene_widgets(gl_window* win, const string& lbl, scene* scn,
    void*& selection, const unordered_map<texture*, gl_texture>& gl_txt);

}  // namespace ygl

#endif

#endif
