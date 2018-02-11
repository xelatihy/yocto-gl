///
/// # Yocto/GL: Tiny C++ Library for Physically-based Graphics
///
/// Yocto/GL is a collection utilities for building physically-based graphics
/// algorithms implemented as a two-file library (`yocto_gl.h`, `yocto_gl.cpp`),
/// and released under the MIT license. Features include:
///
/// - convenience math functions for graphics
/// - static length vectors for 2, 3, 4 length of arbitrary type
/// - static length matrices for 2x2, 3x3, 4x4 of arbitrary type
/// - static length rigid transforms (frames), specialized for 2d and 3d space
/// - linear algebra operations and transforms
/// - axis aligned bounding boxes
/// - rays and ray-primitive intersection
/// - point-primitive distance and overlap tests
/// - normal and tangent computation for meshes and lines
/// - generation of tesselated meshes
/// - mesh refinement with linear tesselation and Catmull-Cark subdivision
/// - keyframed animation, skinning and morphing
/// - random number generation via PCG32
/// - simple image data structure and a few image operations
/// - simple scene format
/// - generation of image examples
/// - generation of scene examples
/// - procedural sun and sky HDR
/// - procedural Perlin noise
/// - BVH for intersection and closest point query
/// - Python-like string, path and container operations
/// - utilities to load and save entire text and binary files
/// - immediate mode command line parser
/// - simple logger
/// - path tracer supporting surfaces and hairs, GGX and MIS
/// - support for loading and saving Wavefront OBJ and Khronos glTF
/// - support for loading Bezier curves from SVG
/// - OpenGL utilities to manage textures, buffers and prograrms
/// - OpenGL shader for image viewing and GGX microfacet and hair rendering
///
/// The current version is 0.1.0. You can access the previous multi-file version
/// with tag "v0.0.1" in this repository.
///
/// ## Credits
///
/// This library includes code from the PCG random number generator,
/// boost hash_combine, Pixar multijittered sampling,
/// code from "Real-Time Collision Detection" by Christer Ericson, base64
/// encode/decode by René Nyffenegger and public domain code from
/// github.com/sgorsten/linalg, gist.github.com/badboy/6267743 and
/// github.com/nothings/stb_perlin.h.
///
/// This library imports many symbols from std for three reasons: avoid
/// verbosity , ensuring better conventions when calling math functions and
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
/// heavy use of operator overloading for math readability. We attempt to make
/// the code easy to use use rather than as performant as possible.
/// We adopt a functional style and only rarely use classes and methods.
/// Using a function style makes the code easier to extend, more explicit in
/// the function requirements, and easier to write parallel-friendly APIs.
/// I guess you could call this "data-driven programming".
///
/// The use of templates in Yocto was the reason for many refactorings, going
/// from no template to heavy templates use. At this time, templates are used
/// in the basic types to make the codebase shorter and reduce bugs,
/// at the price of accessibility for beginners. The truth is that modern C++,
/// a tenant of Yocto, is heavily templated anyway, so being able to read
/// template code is necessary no matter how Yocto does things.
///
/// We make use of exception for error reporting. This makes the code
/// much cleaner and more in line with the expectation of most other programming
/// languages.
///
/// Finally, we often import symbols from the standard library rather than
/// using the `std::name` pattern. We found that this improves consistency
/// when using math functions, and is more readable with templates. We realize
/// this is not standard, but the imports are hidden within the ygl namespace,
/// so library users do not have to be concern about it.
///
///
/// ## Compilation
///
/// Yocto/GL is written in C++14 and compiles on OSX (clang from Xcode 9+),
/// Linux (gcc 6+, clang 4+) and Windows (MSVC 2015, MSVC 2017).
///
/// For image loading and saving, Yocto/GL depends on `stb_image.h`,
/// `stb_image_write.h`, `stb_image_resize.h` and `tinyexr.h`. These features
/// can be disabled by defining YGL_IMAGEIO to 0 before including this file.
/// If these features are useful, then the implementation files need to
/// included in the manner described by the respective libraries. To simplify
/// builds, we provide a file that builds these libraries, `stb_image.cpp`.
///
/// To support Khronos glTF, Yocto/GL depends on `json.hpp`. This feature can
/// be disabled by defining YGL_GLTF to 0 before including this file.
///
/// To support SVG, Yocto/GL depends on `nanosvg.h`. This feature can
/// be disabled by defining YGL_SVG to 0 before including this file.
///
/// OpenGL utilities include the OpenGL libraries, use GLEW on Windows/Linux,
/// GLFW for windows handling and Dear ImGui for UI support.
/// Since OpenGL is quite onerous and hard to link, its support can be disabled
/// by defining YGL_OPENGL to 1 before including this file. If you use any of
/// the OpenGL calls, make sure to properly link to the OpenGL libraries on
/// your system. For ImGUI, build with the libraries `imgui.cpp`,
/// `imgui_draw.cpp`, `imgui_impl_glfw_gl3.cpp`.
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
/// Here are two images rendered with the builtin path tracer, where the
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
/// All library features are documented at their definition and should be
/// relatively easy to use if you are familiar with writing graphics code.
/// You can find the extracted documentation at `yocto_gl.md`.
/// Here we give an overview of some of the main features.
///
///
/// ### Small Vectors and Matrices, Frames, Bounding Boxes and Transforms
///
/// We provide common operations for small vectors and matrices typically used
/// in graphics. In particular, we support 2-4 dimensional vectors of arbitrary
/// `vec<T, 2>`, `vec<T, 3>`, `vec<T, 4>` with specializarion for float
/// (`vec2f`, `vec3f`, `vec4f`), int (`vec2i`, `vec3i`, `vec4i`) and bytes
/// (`vec4b`). Vector operations are templated so they work on every type, but
/// many of them are well-defined only for float types.
///
/// We support 2-4 dimensional generic matrices `mat<T, 2>`, `mat<T, 3>`,
/// `mat<T, 4>`, with matrix-matrix and matrix-vector products, transposes and
/// inverses. Matrices are stored in column-major ordered and are accessed and
/// constructed by column.
///
/// To represent transformations, most of the library facilities prefer the use
/// coordinate frames, aka rigid transforms, represented as `frame<T, 3>`.
/// The structure store three coordinate axis and the frame origin. This is
/// equivalent to a rigid transform written as a column-major affine
/// matrix. Transform operations are better behaved with this representation.
///
/// We represent coordinate bounds with axis-aligned bounding boxes in 1-4
/// dimensions: `bbox<T, 1>`, `bbox<T, 2>`, `bbox<T, 3>`, `bbox<T, 4>`. These
/// types support expansion operation, union and containment. We provide
/// operations to compute bounds for points, lines, triangles and quads.
///
/// For all basic types we support iteration with `begin()`/`end()` pairs,
/// data access with `data()`, `empty()` and `size()` and stream inout and
/// output.
///
/// For both matrices and frames we support transform operations for points,
/// vectors and directions (`trasform_point()`, `trasform_vector()`,
/// `trasform_direction()`). For frames we also the support inverse operations
/// (`transform_xxx_inverse()`). Transform matrices and frames can be
/// constructed from basic translation, rotation and scaling, e.g. with
/// `translation_mat4f()` or `translation_frame3f()` respectively, etc. For
/// rotation we support axis-angle and quaternions, with slerp.
///
///
/// ### Random Number Generation, Noise, Hashing and Monte Carlo support
///
/// This library supports many facilities helpful in writing sampling
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
///     7. generate random shuffled sequences with `rng_shuffle()`
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
/// ### Shape Utilities
///
/// The library contains a few function to help with typically geometry
/// manipulation useful to support scene viewing and path tracing.
///
/// 1. compute line tangents, and triangle and quad areas and normals
/// 2. interpolate values over primitives with `eval_line()`,
///    `eval_triangle()` and `eval_quad()`
/// 3. evaluate Bezier curves and derivatives with `eval_bezier()` and
///    `eval_bezier_derivative()`
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
/// 14. subdivide elements by edge splits with `subdivide_elems_linear()` and
///     `subdivide_vert_linear()`; for an easier interface use
///     `subdivide_lines_linear()`, `subdivide_triangles_linear()`,
///     `subdivide_quads_linear()`
/// 15. Catmull-Clark subdivision surface with `subdivide_vert_catmullclark()`
///     with support for edge and vertex creasing
/// 16. subdivide Bezier with `subdivide_bezier_recursive()` and
///     `subdivide_vert_bezier()`
/// 17. example shapes: `make_cube()`, `make_uvsphere()`, `make_uvhemisphere()`,
///     `make_uvquad()`, `make_uvcube()`, `make_fvcube()`, `make_hair()`,
///     `make_suzanne()`
///
///
/// ### Animation utilities
///
/// The library contains a few function to help with typical animation
/// manipulation useful to support scene viewing.
///
/// 1. evaluate keyframed values with step, linear and bezier interpolation with
///    `eval_keyframed_step()`, `eval_keyframed_linear()`,
///    `eval_keyframed_bezier()`
/// 2. mesh skinning with `compute_matrix_skinning()`
///
///
/// ### Image and color
///
/// Images are stored with the `image` templated structure. The two most used
/// image types are 4-byte per pixel sRGB images `image4b`, or 4-float per
/// pixel HDR images `image4f`.
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
/// normal, occlusion and displacement mapping. Finally, the scene containers
/// cameras and environment maps. Quad support in shapes is experimental and
/// mostly supported for loading and saving.
///
/// For low-level access to OBJ/glTF formats, you are best accessing the formats
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
/// 1. build the bvh with `make_bvh()`
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
/// 1. build the ray-tracing acceleration structure with `make_bvh()`
/// 2. prepare lights for rendering `update_lights()`
/// 3. define rendering params with the `trace_params` structure
/// 4. render blocks of samples with `trace_block()`
///
/// The code can also run in fully asynchronous mode to preview images in a
/// window.
///
/// 1. build the ray-tracing acceleration structure with `make_bvh()`
/// 2. prepare lights for rendering `update_lights()`
/// 3. define rendering params with the `trace_params` structure
/// 4. initialize the progressive rendering buffers
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
/// Again we let the user chose the conversion and set the default to the
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
///   vertex duplication happens thought for same triplets
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
/// complete control over the format or an application wants to have their
/// own scene code added. A higher-level interface is provided by the scene
/// or by `yocto_gltf.h`.
///
/// glTF is a very complex file format and was designed mainly with untyped
/// languages in mind. We attempt to match the glTF low-level interface
/// to C++ as best as it can. Since the code is generated from the schema, we
/// follow glTF naming conventions and typing quite well. To simplify adoption
/// and keep the API relatively simple we use vector as arrays and use
/// pointers to reference to all glTF objects. While this makes it less
/// efficient than it might have been, glTF heavy use of optional values makes
/// this necessary. At the same time, we do not keep track of set/unset values
/// for basic types (int, float, bool) as a compromise for efficiency.
///
/// glTF uses integer indices to access objects.
/// While writing code ourselves we found that we add significant problems
/// since we would use an index to access the wrong type of scene objects.
/// For this reasons, we use an explicit index `glTFid<T>` that can only access
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
///     - check validity with `is_valid()`
///     - update textures/buffers with `update()` functions
///     - delete textures/buffers with `clear()`
///     - bind/unbind textures/buffers with `bind()`/`unbind()`
///     - draw elements with `gl_buffer::draw_elems()`
/// 2. program objects with `gl_program`
///     - program creation with constructor
///     - check validity with `is_valid()`
///     - delete with `clear()`
///     - uniforms with `set_program_uniform()`
///     - vertex attrib with `set_program_vertattr()`
///     - draw elements with `gl_buffer::draw_elems()`
/// 3. image viewing with `gl_stdimage_program`, with support for tone mapping.
/// 4. draw surfaces and hair with GGX/Kayjia-Kay with `gl_stdsurface_program`
///     - initialize the program with constructor
///     - check validity with `is_valid()`
///     - start/end each frame with `begin_frame()`, `end_frame()`
///     - define lights with `set_lights()`
///     - start/end each shape with `begin_shape()`, `end_shape()`
///     - define material Parameters with `set_material()`
///     - define vertices with `set_vert()`
///     - draw elements with `draw_elems()`
/// 5. draw yocto scenes using the above shader
///     - initialize the rendering state with `init_stdsurface_state()`
///     - load/update meshes and textures with `update_stdsurface_state()`
///     - setup draw params using a `gl_stdsurface_params` struct
///     - draw scene with `draw_stdsurface_scene()`
/// 6. also includes other utlities for quick OpenGL hacking
/// 7. GLFW window with `gl_window`
///     - create with constructor
///     - delete with `clear()`
///     - set callbacks with `set_callbacks()`
///     - includes carious utilities to query window, mouse and keyboard
/// 8. immediate mode widgets using ImGui
///     - init with `init_widget()`
///     - use the various widget calls to draw the widget and handle events
///
///
/// ### Other Utilities
///
/// We include additional utilities for writing command line applications and
/// manipulating files.
///
/// 1. Python-like string operations: `startswith()`, `endswith()`,
/// `contains()`,
///    `splitlines()`, `partition()`, `split()`, `splitlines()`, `strip()`,
///    `rstrip()`, `lstrip()`, `join()`, `lower()`, `upper()`, `isspace()`,
///    `replace()`
/// 2. Path-like path operations: `path_dirname()`, `path_extension()`,
///    `path_basename()`, `path_filename()`, `replace_path_extension()`,
///    `prepend_path_extension()`, `split_path()`
/// 3. Python-like format strings (only support for position arguments and no
///    formatting commands): `format()`, `print()`
/// 5. load/save entire files: `load_binary()`, `load_text()`,
///    `save_text()` and `save_binary()`
/// 4. simple logger with support for console and file streams:
///     1. create a `logger`
///     2. add more streams with `add_console_stream()` or `add_file_stream()`
///     3. write log messages with `log_msg()` and its variants
///     4. you can also use a global default logger with the free functions
///        `log_XXX()`
/// 5. timer for simple access to `std::chrono`:
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
/// 4. end cmdline parsing with `check_parsing()` to check for unused values,
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
/// and bug fixes are not reported here.
///
/// - v 0.3.0: templated types, animation and objects in scene, api cleanups
/// - v 0.2.0: various bug fixes and improvement to OpenGL drawing and widgets
/// - v 0.1.0: initial release after refactoring
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

// enable SVG
#ifndef YGL_SVG
#define YGL_SVG 1
#endif

// enable OpenGL
#ifndef YGL_OPENGL
#define YGL_OPENGL 1
#endif

// enable explicit json objects in glTF
#ifndef YGL_GLTFJSON
#define YGL_GLTFJSON 0
#endif

// use iostream for whole file operations
#ifndef YGL_IOSTREAM
#define YGL_IOSTREAM 0
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
#include <set>
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
// BASIC MATH CONSTANTS AND FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup math Basic math constants and functions
/// @{

/// Convenient typedef for bytes.
using byte = unsigned char;

/// Convenient typedef for unsigned ints.
using uint = unsigned int;

/// Pi (float).
const auto pif = 3.14159265f;
/// Pi (double).
const auto pi = 3.1415926535897932384626433832795;

/// Shortcat for float max value.
const auto flt_max = std::numeric_limits<float>::max();
/// Shortcat for float min value.
const auto flt_min = std::numeric_limits<float>::lowest();
/// Shortcat for float epsilon.
const auto flt_eps = std::numeric_limits<float>::epsilon();
/// Shortcat for int max value.
const auto int_max = std::numeric_limits<int>::max();
/// Shortcat for int min value.
const auto int_min = std::numeric_limits<int>::min();

/// Square root.
template <typename T>
inline T sqrt(T a) {
    return std::sqrt(a);
}
/// Power.
template <typename T, typename T1>
inline auto pow(T a, T1 b) {
    return std::pow(a, b);
}
/// Exponential.
template <typename T>
inline T exp(T a) {
    return std::exp(a);
}
/// Logarithm.
template <typename T>
inline T log(T a) {
    return std::log(a);
}
/// Sine.
template <typename T>
inline T sin(T a) {
    return std::sin(a);
}
/// Cosine.
template <typename T>
inline T cos(T a) {
    return std::cos(a);
}
/// Tangent.
template <typename T>
inline T tan(T a) {
    return std::tan(a);
}
/// Arc sine.
template <typename T>
inline T asin(T a) {
    return std::asin(a);
}
/// Arc cosine.
template <typename T>
inline T acos(T a) {
    return std::acos(a);
}
/// Arc tangent.
template <typename T>
inline T atan(T a) {
    return std::atan(a);
}
/// Arc tangent.
template <typename T, typename T1>
inline auto atan2(T a, T1 b) {
    return std::atan2(a, b);
}
/// Absolute value.
template <typename T>
inline T abs(T a) {
    return (a >= 0) ? a : -a;
}
/// Floor.
template <typename T>
inline T floor(T a) {
    return std::floor(a);
}
/// Round.
template <typename T>
inline T round(T a) {
    return std::round(a);
}

/// Safe minimum value.
template <typename T>
inline T min(T x, T y) {
    return (x < y) ? x : y;
}
/// Safe minimum value.
template <typename T>
inline T min(std::initializer_list<T> vs) {
    auto m = int_max;
    for (auto v : vs) m = min(m, v);
    return m;
}
/// Safe maximum value.
template <typename T>
inline T max(T x, T y) {
    return (x > y) ? x : y;
}
/// Safe maximum value.
template <typename T>
inline T max(std::initializer_list<T> vs) {
    auto m = int_min;
    for (auto v : vs) m = max(m, v);
    return m;
}

/// Clamp a value between a minimum and a maximum.
template <typename T>
inline T clamp(T x, T min_, T max_) {
    return min(max(x, min_), max_);
}

/// Linear interpolation.
template <typename T, typename T1>
inline T lerp(const T& a, const T& b, T1 u) {
    return a * (1 - u) + b * u;
}

/// Bilinear interpolation. Order is specified like quads counter-clockwise,
/// so a,b,c,d correspond to parameters (0,0), (0,1), (1,1), (0,1).
template <typename T, typename T1>
inline float bilerp(
    const T& a, const T& b, const T& c, const T& d, T1 u, T1 v) {
    return a * (1 - u) * (1 - v) + b * u * (1 - v) + c * u * v +
           d * (1 - u) * v;
}

/// Integer power of two.
inline int pow2(int x) { return 1 << x; }

/// Fast floor.
inline int fastfloor(float x) {
    auto xi = (int)x;
    return (x < xi) ? xi - 1 : xi;
}

/// Safe float to byte conversion.
inline byte float_to_byte(float x) {
    return (byte)max(0, min(int(x * 256), 255));
}
/// Safe byte to float conversion.
inline float byte_to_float(byte x) { return (float)x / 255.0f; }

/// String literals.
using namespace std::string_literals;
/// Makes literals available
using namespace std::literals;

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// VECTORS
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup vec Fixed-size vectors
/// @{

/// Generic vector of N elements. This is used only to define template
/// specializations for small fixed sized vectors.
template <typename T, int N>
struct vec;

/// Vector of 1 element. Defined only for completeness.
template <typename T>
struct vec<T, 1> {
    /// Default constructor. Initializes to zeros.
    vec() : x{0} {}
    /// Element constructor.
    vec(T x) : x{x} {}

    /// Element access.
    T& operator[](int i) { return (&x)[i]; }
    /// Element access.
    const T& operator[](int i) const { return (&x)[i]; }

    /// Element data.
    T x;
};

/// Vector of 2 elements.
template <typename T>
struct vec<T, 2> {
    /// Default constructor. Initializes to zeros.
    vec() : x{0}, y{0} {}
    /// Element constructor.
    explicit vec(T vv) : x(vv), y(vv) {}
    /// Element constructor.
    vec(T x, T y) : x{x}, y{y} {}

    /// Element access.
    T& operator[](int i) { return (&x)[i]; }
    /// Element access.
    const T& operator[](int i) const { return (&x)[i]; }

    /// Element data.
    T x;
    /// Element data.
    T y;
};

/// Vector of 3 elements.
template <typename T>
struct vec<T, 3> {
    /// Default constructor. Initializes to zeros.
    vec() : x{0}, y{0}, z{0} {}
    /// Element constructor
    explicit vec(T vv) : x(vv), y(vv), z(vv) {}
    /// Element constructor
    vec(T x, T y, T z) : x{x}, y{y}, z{z} {}

    /// Element access
    T& operator[](int i) { return (&x)[i]; }
    /// Element access
    const T& operator[](int i) const { return (&x)[i]; }

    /// Element data
    T x;
    /// Element data
    T y;
    /// Element data
    T z;
};

/// Vector of 4 elements.
template <typename T>
struct vec<T, 4> {
    /// Default constructor.  Initializes to zeros.
    vec() : x{0}, y{0}, z{0}, w{0} {}
    /// Element constructor.
    explicit vec(T vv) : x(vv), y(vv), z(vv), w(vv) {}
    /// Element constructor.
    vec(T x, T y, T z, T w) : x{x}, y{y}, z{z}, w{w} {}

    /// Element access.
    T& operator[](int i) { return (&x)[i]; }
    /// Element access.
    const T& operator[](int i) const { return (&x)[i]; }

    /// Element data.
    T x;
    /// Element data.
    T y;
    /// Element data.
    T z;
    /// Element data.
    T w;
};

/// 1-dimensional float vector.
using vec1f = vec<float, 1>;
/// 2-dimensional float vector.
using vec2f = vec<float, 2>;
/// 3-dimensional float vector
using vec3f = vec<float, 3>;
/// 4-dimensional float vector
using vec4f = vec<float, 4>;
/// 1-dimensional int vector.
using vec1i = vec<int, 1>;
/// 2-dimensional int vector.
using vec2i = vec<int, 2>;
/// 3-dimensional int vector.
using vec3i = vec<int, 3>;
/// 4-dimensional int vector.
using vec4i = vec<int, 4>;
/// 4-dimensional byte vector.
using vec4b = vec<byte, 4>;

/// 1-dimensional float zero vector.
const auto zero1f = vec1f();
/// 2-dimensional float zero vector.
const auto zero2f = vec2f();
/// 3-dimensional float zero vector.
const auto zero3f = vec3f();
/// 4-dimensional float zero vector.
const auto zero4f = vec4f();
/// 1-dimensional int zero vector.
const auto zero1i = vec1i();
/// 2-dimensional int zero vector.
const auto zero2i = vec2i();
/// 3-dimensional int zero vector.
const auto zero3i = vec3i();
/// 4-dimensional int zero vector.
const auto zero4i = vec4i();
/// 4-dimensional byte zero vector.
const auto zero4b = vec4b();

/// Element iteration.
template <typename T, int N>
inline T* begin(vec<T, N>& a) {
    return &a.x;
}
/// Element iteration.
template <typename T, int N>
inline const T* begin(const vec<T, N>& a) {
    return &a.x;
}
/// Element iteration.
template <typename T, int N>
inline T* end(vec<T, N>& a) {
    return &a.x + N;
}
/// Element iteration.
template <typename T, int N>
inline const T* end(const vec<T, N>& a) {
    return &a.x + N;
}
/// Element access.
template <typename T, int N>
inline T* data(vec<T, N>& a) {
    return &a.x;
}
/// Element access.
template <typename T, int N>
inline const T* data(const vec<T, N>& a) {
    return &a.x;
}
/// Number of elements.
template <typename T, int N>
inline int size(vec<T, N>& a) {
    return N;
}
/// Empty check (always false for useful for templated code).
template <typename T, int N>
inline bool empty(vec<T, N>& a) {
    return false;
}

/// Vector equality.
template <typename T>
inline bool operator==(const vec<T, 1>& a, const vec<T, 1>& b) {
    return a.x == b.x;
}
/// Vector inequality.
template <typename T>
inline bool operator!=(const vec<T, 1>& a, const vec<T, 1>& b) {
    return a.x != b.x;
}

/// Vector equality.
template <typename T>
inline bool operator==(const vec<T, 2>& a, const vec<T, 2>& b) {
    return a.x == b.x && a.y == b.y;
}
/// Vector inequality.
template <typename T>
inline bool operator!=(const vec<T, 2>& a, const vec<T, 2>& b) {
    return a.x != b.x || a.y != b.y;
}

/// Vector equality.
template <typename T>
inline bool operator==(const vec<T, 3>& a, const vec<T, 3>& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}
/// Vector inequality.
template <typename T>
inline bool operator!=(const vec<T, 3>& a, const vec<T, 3>& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z;
}

/// Vector equality.
template <typename T>
inline bool operator==(const vec<T, 4>& a, const vec<T, 4>& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}
/// Vector inequality.
template <typename T>
inline bool operator!=(const vec<T, 4>& a, const vec<T, 4>& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w;
}

/// Vector comparison using lexicographic order, useful for map.
template <typename T>
inline bool operator<(const vec<T, 2>& a, const vec<T, 2>& b) {
    for (auto i = 0; i < 2; i++) {
        if (a[i] < b[i]) return true;
        if (a[i] > b[i]) return false;
    }
    return false;
}
/// Vector comparison using lexicographic order, useful for map.
template <typename T>
inline bool operator<(const vec<T, 3>& a, const vec<T, 3>& b) {
    for (auto i = 0; i < 3; i++) {
        if (a[i] < b[i]) return true;
        if (a[i] > b[i]) return false;
    }
    return false;
}
/// Vector comparison using lexicographic order, useful for map.
template <typename T>
inline bool operator<(const vec<T, 4>& a, const vec<T, 4>& b) {
    for (auto i = 0; i < 4; i++) {
        if (a[i] < b[i]) return true;
        if (a[i] > b[i]) return false;
    }
    return false;
}

/// Vector unary plus (for completeness).
template <typename T>
inline vec<T, 2> operator+(const vec<T, 2>& a) {
    return a;
}
/// Vector negation.
template <typename T>
inline vec<T, 2> operator-(const vec<T, 2>& a) {
    return {-a.x, -a.y};
}
/// Vector sum.
template <typename T>
inline vec<T, 2> operator+(const vec<T, 2>& a, const vec<T, 2>& b) {
    return {a.x + b.x, a.y + b.y};
}
/// Vector difference.
template <typename T>
inline vec<T, 2> operator-(const vec<T, 2>& a, const vec<T, 2>& b) {
    return {a.x - b.x, a.y - b.y};
}
/// Vector scalar product.
template <typename T>
inline vec<T, 2> operator*(const vec<T, 2>& a, const vec<T, 2>& b) {
    return {a.x * b.x, a.y * b.y};
}
/// Vector scalar product.
template <typename T, typename T1>
inline vec<T, 2> operator*(const vec<T, 2>& a, T1 b) {
    return {a.x * b, a.y * b};
}
/// Vector scalar product.
template <typename T>
inline vec<T, 2> operator*(float a, const vec<T, 2>& b) {
    return {a * b.x, a * b.y};
}
/// Vector scalar division.
template <typename T>
inline vec<T, 2> operator/(const vec<T, 2>& a, const vec<T, 2>& b) {
    return {a.x / b.x, a.y / b.y};
}
/// Vector scalar division.
template <typename T, typename T1>
inline vec<T, 2> operator/(const vec<T, 2>& a, T1 b) {
    return {a.x / b, a.y / b};
}
/// Vector scalar division.
template <typename T, typename T1>
inline vec<T, 2> operator/(T1 a, const vec<T, 2>& b) {
    return {a / b.x, a / b.y};
}

/// Vector unary plus (for completeness).
template <typename T>
inline vec<T, 3> operator+(const vec<T, 3>& a) {
    return a;
}
/// Vector negation.
template <typename T>
inline vec<T, 3> operator-(const vec<T, 3>& a) {
    return {-a.x, -a.y, -a.z};
}
/// Vector sum.
template <typename T>
inline vec<T, 3> operator+(const vec<T, 3>& a, const vec<T, 3>& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
/// Vector operator -.
template <typename T>
inline vec<T, 3> operator-(const vec<T, 3>& a, const vec<T, 3>& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
/// Vector scalar product.
template <typename T>
inline vec<T, 3> operator*(const vec<T, 3>& a, const vec<T, 3>& b) {
    return {a.x * b.x, a.y * b.y, a.z * b.z};
}
/// Vector scalar product.
template <typename T, typename T1>
inline vec<T, 3> operator*(const vec<T, 3>& a, T1 b) {
    return {a.x * b, a.y * b, a.z * b};
}
/// Vector scalar product.
template <typename T, typename T1>
inline vec<T, 3> operator*(T1 a, const vec<T, 3>& b) {
    return {a * b.x, a * b.y, a * b.z};
}
/// Vector scalar division.
template <typename T>
inline vec<T, 3> operator/(const vec<T, 3>& a, const vec<T, 3>& b) {
    return {a.x / b.x, a.y / b.y, a.z / b.z};
}
/// Vector scalar division.
template <typename T, typename T1>
inline vec<T, 3> operator/(const vec<T, 3>& a, T1 b) {
    return {a.x / b, a.y / b, a.z / b};
}
/// Vector scalar division.
template <typename T, typename T1>
inline vec<T, 3> operator/(T1 a, const vec<T, 3>& b) {
    return {a / b.x, a / b.y, a / b.z};
}

/// Vector unary plus (for completeness).
template <typename T>
inline vec<T, 4> operator+(const vec<T, 4>& a) {
    return a;
}
/// Vector negation.
template <typename T>
inline vec<T, 4> operator-(const vec<T, 4>& a) {
    return {-a.x, -a.y, -a.z, -a.w};
}
/// Vector sum.
template <typename T>
inline vec<T, 4> operator+(const vec<T, 4>& a, const vec<T, 4>& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}
/// Vector difference.
template <typename T>
inline vec<T, 4> operator-(const vec<T, 4>& a, const vec<T, 4>& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w};
}
/// Vector scalar product.
template <typename T>
inline vec<T, 4> operator*(const vec<T, 4>& a, const vec<T, 4>& b) {
    return {a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w};
}
/// Vector scalar product.
template <typename T>
inline vec<T, 4> operator*(const vec<T, 4>& a, float b) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
}
/// Vector scalar product.
template <typename T>
inline vec<T, 4> operator*(float a, const vec<T, 4>& b) {
    return {a * b.x, a * b.y, a * b.z, a * b.w};
}
/// Vector scalar division.
template <typename T>
inline vec<T, 4> operator/(const vec<T, 4>& a, const vec<T, 4>& b) {
    return {a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w};
}
/// Vector scalar division.
template <typename T, typename T1>
inline vec<T, 4> operator/(const vec<T, 4>& a, T1 b) {
    return {a.x / b, a.y / b, a.z / b, a.w / b};
}
/// Vector scalar division.
template <typename T, typename T1>
inline vec<T, 4> operator/(T1 a, const vec<T, 4>& b) {
    return {a / b.x, a / b.y, a / b.z, a / b.w};
}

/// Vector assignment.
template <typename T, int N>
inline vec<T, N>& operator+=(vec<T, N>& a, const vec<T, N>& b) {
    return a = a + b;
}
/// Vector assignment.
template <typename T, int N>
inline vec<T, N>& operator-=(vec<T, 2>& a, const vec<T, N>& b) {
    return a = a - b;
}
/// Vector assignment.
template <typename T, int N>
inline vec<T, N>& operator*=(vec<T, N>& a, const vec<T, N>& b) {
    return a = a * b;
}
/// Vector assignment.
template <typename T, int N, typename T1>
inline vec<T, N>& operator*=(vec<T, N>& a, T1 b) {
    return a = a * b;
}
/// Vector assignment.
template <typename T, int N>
inline vec<T, N>& operator/=(vec<T, N>& a, const vec<T, N>& b) {
    return a = a / b;
}
/// Vector assignment.
template <typename T, int N, typename T1>
inline vec<T, N>& operator/=(vec<T, N>& a, T1 b) {
    return a = a / b;
}

/// Vector dot product.
template <typename T>
inline T dot(const vec<T, 2>& a, const vec<T, 2>& b) {
    return a.x * b.x + a.y * b.y;
}
/// Vector dot product.
template <typename T>
inline T dot(const vec<T, 3>& a, const vec<T, 3>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
/// Vector dot product.
template <typename T>
inline T dot(const vec<T, 4>& a, const vec<T, 4>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

/// Vector cross product.
template <typename T>
inline T cross(const vec<T, 2>& a, const vec<T, 2>& b) {
    return a.x * b.y - a.y * b.x;
}
/// Vector cross product.
template <typename T>
inline vec<T, 3> cross(const vec<T, 3>& a, const vec<T, 3>& b) {
    return {
        a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

/// Vector length.
template <typename T, int N>
inline T length(const vec<T, N>& a) {
    return sqrt(dot(a, a));
}

/// Vector normalization.
template <typename T, int N>
inline vec<T, N> normalize(const vec<T, N>& a) {
    auto l = length(a);
    if (l == 0) return a;
    return a * (1 / l);
}

/// Angle between vectors.
template <typename T, int N>
inline T angle(const vec<T, N>& a, const vec<T, N>& b) {
    auto d = clamp(dot(normalize(a), normalize(b)), (T)-1, (T)1);
    return acos(d);
}

/// Vector spherical linear interpolation (vectors have to be normalized).
template <typename T, int N, typename T1>
inline vec<T, N> slerp(const vec<T, N>& a, const vec<T, N>& b, T1 u) {
    auto th = angle(a, b);
    if (!th) return a;
    return a * (sin(th * (1 - u)) / sin(th)) + b * (sin(th * u) / sin(th));
}

/// Orthogonal vector.
template <typename T>
inline vec<T, 3> orthogonal(const vec<T, 3>& v) {
    // http://lolengine.net/blog/2013/09/21/picking-orthogonal-vector-combing-coconuts)
    return abs(v.x) > abs(v.z) ? vec<T, 3>{-v.y, v.x, 0} :
                                 vec<T, 3>{0, -v.z, v.y};
}

/// Orthonormalize two vectors.
template <typename T>
inline vec<T, 3> orthonormalize(const vec<T, 3>& a, const vec<T, 3>& b) {
    return normalize(a - b * dot(a, b));
}

/// Reflected vector.
template <typename T>
inline vec<T, 3> reflect(const vec<T, 3>& w, const vec<T, 3>& n) {
    return -w + 2 * dot(n, w) * n;
}

/// Refracted vector.
template <typename T>
inline vec<T, 3> refract(const vec<T, 3>& w, const vec<T, 3>& n, T eta) {
    // auto k = 1.0 - eta * eta * (1.0 - dot(n, w) * dot(n, w));
    auto k = 1 - eta * eta * max((T)0, 1 - dot(n, w) * dot(n, w));
    if (k < 0) return zero3f;  // tir
    return -w * eta + (eta * dot(n, w) - sqrt(k)) * n;
}

/// Component-wise clamp.
template <typename T, typename T1>
inline vec<T, 2> clamp(const vec<T, 2>& x, T1 min, T1 max) {
    return {clamp(x.x, min, max), clamp(x.y, min, max)};
}
/// Component-wise clamp.
template <typename T, typename T1>
inline vec<T, 3> clamp(const vec<T, 3>& x, T1 min, T1 max) {
    return {clamp(x.x, min, max), clamp(x.y, min, max), clamp(x.z, min, max)};
}
/// Component-wise clamp.
template <typename T, typename T1>
inline vec<T, 4> clamp(const vec<T, 4>& x, T1 min, T1 max) {
    return {clamp(x.x, min, max), clamp(x.y, min, max), clamp(x.z, min, max),
        clamp(x.w, min, max)};
}

/// Clamp a vector to a maximum length.
template <typename T, int N, typename T1>
inline vec<T, N> clamplen(const vec<T, N>& x, T1 max) {
    auto l = length(x);
    return (l > max) ? x * max / l : x;
}

/// Index of minimum element.
template <typename T, int N>
inline int min_element(const vec<T, N>& a) {
    auto v = std::numeric_limits<T>::max();
    auto pos = -1;
    for (auto i = 0; i < N; i++) {
        if (v > a[i]) {
            v = a[i];
            pos = i;
        }
    }
    return pos;
}
/// Value of minimum element.
template <typename T, int N>
inline T min_element_value(const vec<T, N>& a) {
    return a[min_element(a)];
}

/// Index of maximum element.
template <typename T, int N>
inline int max_element(const vec<T, N>& a) {
    auto v = std::numeric_limits<T>::lowest();
    auto pos = -1;
    for (auto i = 0; i < N; i++) {
        if (v < a[i]) {
            v = a[i];
            pos = i;
        }
    }
    return pos;
}
/// Value of maximum element.
template <typename T, int N>
inline T max_element_value(const vec<T, N>& a) {
    return a[max_element(a)];
}

/// Element-wise float to byte conversion.
inline vec4b float_to_byte(const vec4f& a) {
    return {float_to_byte(a.x), float_to_byte(a.y), float_to_byte(a.z),
        float_to_byte(a.w)};
}
/// Element-wise byte to float conversion.
inline vec4f byte_to_float(const vec4b& a) {
    return {byte_to_float(a.x), byte_to_float(a.y), byte_to_float(a.z),
        byte_to_float(a.w)};
}

/// Stream write.
template <typename T, int N>
inline std::ostream& operator<<(std::ostream& os, const vec<T, N>& a) {
    for (auto i = 0; i < N; i++) {
        if (i) os << ' ';
        os << data(a)[i];
    }
    return os;
}
/// Stream read.
template <typename T, int N>
inline std::istream& operator>>(std::istream& is, vec<T, N>& a) {
    for (auto i = 0; i < N; i++) is >> data(a)[i];
    return is;
}

/// @}

}  // namespace ygl

namespace std {
/// Hash functor for vector for use with unordered_map
template <typename T, int N>
struct hash<ygl::vec<T, N>> {
    // from boost::hash_combine
    static size_t hash_combine(size_t h, size_t h1) {
        h ^= h1 + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
    size_t operator()(const ygl::vec<T, N>& v) const {
        auto vh = hash<int>();
        auto h = (size_t)0;
        for (auto i = 0; i < N; i++) h = hash_combine(h, vh(v[i]));
        return h;
    }
};
}  // namespace std

// -----------------------------------------------------------------------------
// MATRICES
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup mat Fixed-size matrices
/// @{

/// Generic matrix of NxN elements. This is used only to define template
/// specializations for small fixed-sized matrices.
template <typename T, int N>
struct mat;

/// Matrix of 2x2 elements stored in column-major format.
template <typename T>
struct mat<T, 2> {
    /// Default constructor. Initializes to identity matrix.
    mat() : x{1, 0}, y{0, 1} {}
    /// Constructs a matrix with the given diagonal.
    explicit mat(T vv) : x{vv, 0}, y{0, vv} {}
    /// Constructs a matrix from its columns.
    mat(const vec<T, 2>& x, const vec<T, 2>& y) : x(x), y(y) {}

    /// Column access.
    vec<T, 2>& operator[](int i) { return (&x)[i]; }
    /// Column access.
    const vec<T, 2>& operator[](int i) const { return (&x)[i]; }

    /// Column data.
    vec<T, 2> x;
    /// Column data.
    vec<T, 2> y;
};

/// Matrix of 3x3 elements stored in column major format.
/// Colums access via operator[].
template <typename T>
struct mat<T, 3> {
    /// Default constructor. Initializes to identity matrix.
    mat() : x{1, 0, 0}, y{0, 1, 0}, z{0, 0, 1} {}
    /// Constructs a matrix with the given diagonal.
    explicit mat(T vv) : x{vv, 0, 0}, y{0, vv, 0}, z{0, 0, vv} {}
    /// Constructs a matrix from its columns.
    mat(const vec<T, 3>& x, const vec<T, 3>& y, const vec<T, 3>& z)
        : x(x), y(y), z(z) {}

    /// Column access.
    vec<T, 3>& operator[](int i) { return (&x)[i]; }
    /// Column access.
    const vec<T, 3>& operator[](int i) const { return (&x)[i]; }

    /// Column data.
    vec<T, 3> x;
    /// Column data.
    vec<T, 3> y;
    /// Column data.
    vec<T, 3> z;
};

/// Matrix of 4x4 elements stored in column major format.
/// Colums access via operator[].
template <typename T>
struct mat<T, 4> {
    /// Default constructor. Initializes to identity matrix.
    mat() : x{1, 0, 0, 0}, y{0, 1, 0, 0}, z{0, 0, 1, 0}, w{0, 0, 0, 1} {}
    /// Constructs a matrix with the given diagonal.
    explicit mat(float vv)
        : x{vv, 0, 0, 0}, y{0, vv, 0, 0}, z{0, 0, vv, 0}, w{0, 0, 0, vv} {}
    /// Constructs a matrix from its columns.
    mat(const vec<T, 4>& x, const vec<T, 4>& y, const vec<T, 4>& z,
        const vec<T, 4>& w)
        : x(x), y(y), z(z), w(w) {}

    /// Column access.
    vec<T, 4>& operator[](int i) { return (&x)[i]; }
    /// Column access.
    const vec<T, 4>& operator[](int i) const { return (&x)[i]; }

    /// Column data.
    vec<T, 4> x;
    /// Column data.
    vec<T, 4> y;
    /// Column data.
    vec<T, 4> z;
    /// Column data.
    vec<T, 4> w;
};

/// 2-dimensional float matrix.
using mat2f = mat<float, 2>;
/// 3-dimensional float matrix.
using mat3f = mat<float, 3>;
/// 4-dimensional float matrix.
using mat4f = mat<float, 4>;

/// 2-dimensional float identity matrix.
const auto identity_mat2f = mat2f();
/// 3-dimensional float identity matrix.
const auto identity_mat3f = mat3f();
/// 4-dimensional float identity matrix.
const auto identity_mat4f = mat4f();

/// Column iteration.
template <typename T, int N>
inline vec<T, N>* begin(mat<T, N>& m) {
    return &(m.x);
}
/// Column iteration.
template <typename T, int N>
inline vec<T, N>* end(mat<T, N>& m) {
    return &(m.x) + N;
}
/// Column iteration.
template <typename T, int N>
inline const vec<T, N>* begin(const mat<T, N>& m) {
    return &(m.x);
}
/// Column iteration.
template <typename T, int N>
inline const vec<T, N>* end(const mat<T, N>& m) {
    return &(m.x) + N;
}
/// Column access.
template <typename T, int N>
inline vec<T, N>* data(mat<T, N>& m) {
    return &(m.x);
}
/// Column access.
template <typename T, int N>
inline const vec<T, N>* data(const mat<T, N>& m) {
    return &(m.x);
}
/// Number of columns.
template <typename T, int N>
inline int size(mat<T, N>& a) {
    return N;
}
/// Empty check (always false for useful for templated code).
template <typename T, int N>
inline bool empty(mat<T, N>& a) {
    return false;
}

/// Matrix equality.
template <typename T>
inline bool operator==(const mat<T, 2>& a, const mat<T, 2>& b) {
    return a.x == b.x && a.y == b.y;
}
/// Matrix inequality.
template <typename T>
inline bool operator!=(const mat<T, 2>& a, const mat<T, 2>& b) {
    return !(a == b);
}
/// Matrix equality.
template <typename T>
inline bool operator==(const mat<T, 3>& a, const mat<T, 3>& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}
/// Matrix inequality.
template <typename T>
inline bool operator!=(const mat<T, 3>& a, const mat<T, 3>& b) {
    return !(a == b);
}
/// Matrix equality.
template <typename T>
inline bool operator==(const mat<T, 4>& a, const mat<T, 4>& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}
/// Matrix inequality.
template <typename T>
inline bool operator!=(const mat<T, 4>& a, const mat<T, 4>& b) {
    return !(a == b);
}

/// Matrix sum.
template <typename T>
inline mat<T, 2> operator+(const mat<T, 2>& a, const mat<T, 2>& b) {
    return {a.x + b.x, a.y + b.y};
}
/// Matrix scalar product.
template <typename T, typename T1>
inline mat<T, 2> operator*(const mat<T, 2>& a, T1 b) {
    return {a.x * b, a.y * b};
}
/// matrix scalar division.
template <typename T, typename T1>
inline mat<T, 2> operator/(const mat<T, 2>& a, T1 b) {
    return {a.x / b, a.y / b};
}
/// Matrix-vector product.
template <typename T>
inline vec<T, 2> operator*(const mat<T, 2>& a, const vec<T, 2>& b) {
    return a.x * b.x + a.y * b.y;
}
/// Matrix-vector product.
template <typename T>
inline vec<T, 2> operator*(const vec<T, 2>& a, const mat<T, 2>& b) {
    return {dot(a, b.x), dot(a, b.y)};
}
/// Matrix-matrix product.
template <typename T>
inline mat<T, 2> operator*(const mat<T, 2>& a, const mat<T, 2>& b) {
    return {a * b.x, a * b.y};
}

/// Matrix sum.
template <typename T>
inline mat<T, 3> operator+(const mat<T, 3>& a, const mat<T, 3>& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
/// Matrix scalar product.
template <typename T, typename T1>
inline mat<T, 3> operator*(const mat<T, 3>& a, T1 b) {
    return {a.x * b, a.y * b, a.z * b};
}
/// Matrix scalar division.
template <typename T, typename T1>
inline mat<T, 3> operator/(const mat<T, 3>& a, T1 b) {
    return {a.x / b, a.y / b, a.z / b};
}
/// Matrix-vector product.
template <typename T>
inline vec<T, 3> operator*(const mat<T, 3>& a, const vec<T, 3>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
/// Matrix-vector product.
template <typename T>
inline vec<T, 3> operator*(const vec<T, 3>& a, const mat<T, 3>& b) {
    return {dot(a, b.x), dot(a, b.y), dot(a, b.z)};
}
/// Matrix-matrix product.
template <typename T>
inline mat<T, 3> operator*(const mat<T, 3>& a, const mat<T, 3>& b) {
    return {a * b.x, a * b.y, a * b.z};
}

/// Matrix sum.
template <typename T>
inline mat<T, 4> operator+(const mat<T, 4>& a, const mat<T, 4>& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}
/// Matrix scalar product.
template <typename T, typename T1>
inline mat<T, 4> operator*(const mat<T, 4>& a, T1 b) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
}
/// Matrix scalar division.
template <typename T, typename T1>
inline mat<T, 4> operator/(const mat<T, 4>& a, T1 b) {
    return {a.x / b, a.y / b, a.z / b, a.w / b};
}
/// Matrix-vector product.
template <typename T>
inline vec<T, 4> operator*(const mat<T, 4>& a, const vec<T, 4>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}
/// Matrix-vector product.
template <typename T>
inline vec<T, 4> operator*(const vec<T, 4>& a, const mat<T, 4>& b) {
    return {dot(a, b.x), dot(a, b.y), dot(a, b.z), dot(a, b.w)};
}
/// Matrix-matrix product.
template <typename T>
inline mat<T, 4> operator*(const mat<T, 4>& a, const mat<T, 4>& b) {
    return {a * b.x, a * b.y, a * b.z, a * b.w};
}

/// Matrix assignment.
template <typename T, int N>
inline mat<T, N>& operator+=(mat<T, N>& a, const mat<T, N>& b) {
    return a = a + b;
}
/// Matrix assignment.
template <typename T, int N>
inline mat<T, N>& operator*=(mat<T, N>& a, const mat<T, N>& b) {
    return a = a * b;
}
/// Matrix assignment.
template <typename T, int N, typename T1>
inline mat<T, N>& operator*=(mat<T, N>& a, T1 b) {
    return a = a * b;
}
/// Matrix assignment.
template <typename T, int N, typename T1>
inline mat<T, N>& operator/=(mat<T, N>& a, T1 b) {
    return a = a / b;
}

/// Matrix diagonal.
template <typename T>
inline vec<T, 2> mat_diagonal(const mat<T, 2>& a) {
    return {a.x.x, a.y.y};
}
/// Matrix diagonal.
template <typename T>
inline vec<T, 3> mat_diagonal(const mat<T, 3>& a) {
    return {a.x.x, a.y.y, a.z.z};
}
/// Matrix diagonal.
template <typename T>
inline vec<T, 4> mat_diagonal(const mat<T, 4>& a) {
    return {a.x.x, a.y.y, a.z.z, a.w.w};
}

/// Matrix transpose.
template <typename T>
inline mat<T, 2> transpose(const mat<T, 2>& a) {
    return {{a.x.x, a.y.x}, {a.x.y, a.y.y}};
}
/// Matrix transpose.
template <typename T>
inline mat<T, 3> transpose(const mat<T, 3>& a) {
    return {
        {a.x.x, a.y.x, a.z.x}, {a.x.y, a.y.y, a.z.y}, {a.x.z, a.y.z, a.z.z}};
}
/// Matrix transpose.
template <typename T>
inline mat<T, 4> transpose(const mat<T, 4>& a) {
    return {{a.x.x, a.y.x, a.z.x, a.w.x}, {a.x.y, a.y.y, a.z.y, a.w.y},
        {a.x.z, a.y.z, a.z.z, a.w.z}, {a.x.w, a.y.w, a.z.w, a.w.w}};
}

/// Matrix adjugate.
template <typename T>
inline mat<T, 2> adjugate(const mat<T, 2>& a) {
    return {{a.y.y, -a.x.y}, {-a.y.x, a.x.x}};
}
/// Matrix adjugate.
template <typename T>
inline mat<T, 3> adjugate(const mat<T, 3>& a) {
    return {{a.y.y * a.z.z - a.z.y * a.y.z, a.z.y * a.x.z - a.x.y * a.z.z,
                a.x.y * a.y.z - a.y.y * a.x.z},
        {a.y.z * a.z.x - a.z.z * a.y.x, a.z.z * a.x.x - a.x.z * a.z.x,
            a.x.z * a.y.x - a.y.z * a.x.x},
        {a.y.x * a.z.y - a.z.x * a.y.y, a.z.x * a.x.y - a.x.x * a.z.y,
            a.x.x * a.y.y - a.y.x * a.x.y}};
}
/// Matrix adjugate.
template <typename T>
inline mat<T, 4> adjugate(const mat<T, 4>& a) {
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

/// Matrix determinant.
template <typename T>
inline T determinant(const mat<T, 2>& a) {
    return a.x.x * a.y.y - a.x.y * a.y.x;
}
/// Matrix determinant.
template <typename T>
inline T determinant(const mat<T, 3>& a) {
    return a.x.x * (a.y.y * a.z.z - a.z.y * a.y.z) +
           a.x.y * (a.y.z * a.z.x - a.z.z * a.y.x) +
           a.x.z * (a.y.x * a.z.y - a.z.x * a.y.y);
}
/// Matrix determinant.
template <typename T>
inline T determinant(const mat<T, 4>& a) {
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

/// Matrix inverse.
template <typename T, int N>
inline mat<T, N> inverse(const mat<T, N>& a) {
    return adjugate(a) / determinant(a);
}

/// Stream write.
template <typename T, int N>
inline std::ostream& operator<<(std::ostream& os, const mat<T, N>& a) {
    for (auto i = 0; i < N; i++) {
        if (i) os << ' ';
        os << data(a)[i];
    }
    return os;
}
/// Stream read.
template <typename T, int N>
inline std::istream& operator>>(std::istream& is, mat<T, N>& a) {
    for (auto i = 0; i < N; i++) is >> data(a)[i];
    return is;
}

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// RIGID BODY TRANSFORMS/FRAMES
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup frame Rigid-body frames
/// @{

/// Generic frame of N elements. This is used only to define template
/// specializations for small fixed sized frames.
template <typename T, int N>
struct frame;

/// Rigid transforms stored as a column-major affine matrix.
/// In memory, this representation is equivalent to storing an NxN rotation
/// followed by a Nx1 translation. Viewed this way, the representation allows
/// also to retrive the axis of the coordinate frame as the first N columns and
/// the translation as the (N+1)-th column. Colums access via operator[].
/// Access rotation and position with pos() and rot().
template <typename T>
struct frame<T, 3> {
    /// Default constructor. Initializes to the identity frame.
    frame() : x{1, 0, 0}, y{0, 1, 0}, z{0, 0, 1}, o{0, 0, 0} {}
    /// Basic and origin constructor. Equavalent to columns of affine matrix.
    frame(const vec<T, 3>& x, const vec<T, 3>& y, const vec<T, 3>& z,
        const vec<T, 3>& o)
        : x(x), y(y), z(z), o(o) {}
    /// Rotation and traslation constructor.
    frame(const mat<T, 3>& m, const vec<T, 3>& t)
        : x(m.x), y(m.y), z(m.z), o(t) {}

    /// Element/column access
    vec<T, 3>& operator[](int i) { return (&x)[i]; }
    /// Element/column access
    const vec<T, 3>& operator[](int i) const { return (&x)[i]; }

    /// Axes and origin data
    vec<T, 3> x;
    /// Axes and origin data
    vec<T, 3> y;
    /// Axes and origin data
    vec<T, 3> z;
    /// Axes and origin data
    vec<T, 3> o;
};

/// 3-dimensional float frame.
using frame3f = frame<float, 3>;

/// Indentity frame.
const auto identity_frame3f =
    frame3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};

/// Element/column iteration.
template <typename T, int N>
inline vec<T, N>* begin(frame<T, N>& a) {
    return &a.x;
}
/// Element/column iteration.
template <typename T, int N>
inline const vec<T, N>* begin(const frame<T, N>& a) {
    return &a.x;
}
/// Element/column iteration.
template <typename T, int N>
inline vec<T, N>* end(frame<T, N>& a) {
    return &a.x + 4;
}
/// Element/column iteration.
template <typename T, int N>
inline const vec<T, N>* end(const frame<T, N>& a) {
    return &a.x + 4;
}
/// Element/column access.
template <typename T, int N>
inline vec<T, N>* data(frame<T, N>& a) {
    return &a.x;
}
/// Element/column access.
template <typename T, int N>
inline const vec<T, N>* data(const frame<T, N>& a) {
    return &a.x;
}
/// Number of columns in the underlying affine matrix.
template <typename T, int N>
inline int size(frame<T, N>& a) {
    return N + 1;
}
/// Empty check (always false for useful for templated code).
template <typename T, int N>
inline bool empty(frame<T, N>& a) {
    return false;
}

// Initializes a frame from origin and z.
template <typename T>
inline frame<T, 3> make_frame_fromz(const vec<T, 3>& o, const vec<T, 3>& z_) {
    auto z = normalize(z_);
    auto x = normalize(orthogonal(z));
    auto y = normalize(cross(z, x));
    return {x, y, z, o};
}
// Initializes a frame3 from origin, z and x.
template <typename T>
inline frame<T, 3> make_frame_fromzx(
    const vec<T, 3>& o, const vec<T, 3>& z_, const vec<T, 3>& x_) {
    auto z = normalize(z_);
    auto x = orthonormalize(x_, z);
    auto y = normalize(cross(z, x));
    return {x, y, z, o};
}

/// Frame to matrix conversion.
template <typename T>
inline mat<T, 4> frame_to_mat(const frame<T, 3>& a) {
    return {{a.x.x, a.x.y, a.x.z, 0}, {a.y.x, a.y.y, a.y.z, 0},
        {a.z.x, a.z.y, a.z.z, 0}, {a.o.x, a.o.y, a.o.z, 1}};
}
/// Matrix to frame conversion.
template <typename T>
inline frame<T, 3> mat_to_frame(const mat<T, 4>& a) {
    return {{a.x.x, a.x.y, a.x.z}, {a.y.x, a.y.y, a.y.z}, {a.z.x, a.z.y, a.z.z},
        {a.w.x, a.w.y, a.w.z}};
}

/// Frame origin.
template <typename T, int N>
vec<T, N>& frame_pos(frame<T, N>& a) {
    return a.o;
}
/// Frame origin.
template <typename T, int N>
const vec<T, N>& frame_pos(const frame<T, N>& a) {
    return a.o;
}
/// Frame rotation
template <typename T, int N>
mat<T, 3>& frame_rot(frame<T, N>& a) {
    return *(mat<T, 3>*)(&a.x);
}
/// Frame rotation
template <typename T, int N>
const mat<T, 3>& frame_rot(const frame<T, N>& a) {
    return *(mat<T, 3>*)(&a.x);
}

/// Frame equality.
template <typename T>
inline bool operator==(const frame<T, 3>& a, const frame<T, 3>& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.o == b.o;
}
/// Frame inequality.
template <typename T>
inline bool operator!=(const frame<T, 3>& a, const frame<T, 3>& b) {
    return !(a == b);
}

/// Frame composition, equivalent to affine matrix product.
template <typename T>
inline frame<T, 3> operator*(const frame<T, 3>& a, const frame<T, 3>& b) {
    return {mat<T, 3>{a.x, a.y, a.z} * mat<T, 3>{b.x, b.y, b.z},
        mat<T, 3>{a.x, a.y, a.z} * b.o + a.o};
}

/// Frame inverse, equivalent to rigid affine inverse.
template <typename T>
inline frame<T, 3> inverse(const frame<T, 3>& a) {
    auto minv = transpose(mat<T, 3>{a.x, a.y, a.z});
    return {minv, -(minv * a.o)};
}

/// Stream write.
template <typename T, int N>
inline std::ostream& operator<<(std::ostream& os, const frame<T, N>& a) {
    for (auto i = 0; i < N + 1; i++) {
        if (i) os << ' ';
        os << data(a)[i];
    }
    return os;
}
/// Stream read.
template <typename T, int N>
inline std::istream& operator>>(std::istream& is, frame<T, N>& a) {
    for (auto i = 0; i < N + 1; i++) is >> data(a)[i];
    return is;
}

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// QUATERNIONS
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup quat Quaternions
/// @{

/// Generic quaternion of N elements. This is used only to define template
/// specializations for small fixed sized quaternions.
template <typename T, int N>
struct quat;

/// Quaternions implemented as 4-dimensional vector as xi + yj + zk + w.
/// Element access via operator[]. The cosde here assume the use as unit
/// quaternions for rotations.
template <typename T>
struct quat<T, 4> {
    /// Default constructor. Initializes to identity rotation.
    quat() : x{0}, y{0}, z{0}, w{1} {}
    // Element consgtructor.
    quat(T x, T y, T z, T w) : x{x}, y{y}, z{z}, w{w} {}
    /// Conversion from vec.
    explicit quat(const vec<T, 4>& vv) : x{vv.x}, y{vv.y}, z{vv.z}, w{vv.w} {}
    /// Conversion to vec.
    explicit operator vec<T, 4>() const { return {x, y, z, w}; }

    /// Element access.
    T& operator[](int i) { return (&x)[i]; }
    /// Element access.
    const T& operator[](int i) const { return (&x)[i]; }

    /// Element data.
    T x;
    /// Element data.
    T y;
    /// Element data.
    T z;
    /// Element data.
    T w;
};

/// 4-dimensional float quaternion.
using quat4f = quat<float, 4>;

/// Float identity quaternion.
const auto identity_quat4f = quat4f{0, 0, 0, 1};

/// Element iteration.
template <typename T, int N>
inline T* begin(quat<T, N>& a) {
    return &a.x;
}
/// Element iteration.
template <typename T, int N>
inline const T* begin(const quat<T, N>& a) {
    return &a.x;
}
/// Element iteration.
template <typename T, int N>
inline T* end(quat<T, N>& a) {
    return &a.x + N;
}
/// Element iteration.
template <typename T, int N>
inline const T* end(const quat<T, N>& a) {
    return &a.x + N;
}
/// Element access.
template <typename T, int N>
inline T* data(quat<T, N>& a) {
    return &a.x;
}
/// Element access.
template <typename T, int N>
inline const T* data(const quat<T, N>& a) {
    return &a.x;
}
/// Number of elements.
template <typename T, int N>
inline int size(quat<T, N>& a) {
    return N;
}
/// Empty check (always false for useful for templated code).
template <typename T, int N>
inline bool empty(quat<T, N>& a) {
    return false;
}

/// Quaternion equality.
template <typename T>
inline bool operator==(const quat<T, 4>& a, const quat<T, 4>& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}
/// Quaternion inequality.
template <typename T>
inline bool operator!=(const quat<T, 4>& a, const quat<T, 4>& b) {
    return !(a == b);
}

/// Quaternion sum.
template <typename T>
inline quat<T, 4> operator+(const quat<T, 4>& a, const quat<T, 4>& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}
/// Quaternion scalar product.
template <typename T, typename T1>
inline quat<T, 4> operator*(const quat<T, 4>& a, T1 b) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
}
/// Quaternion scalar division.
template <typename T, typename T1>
inline quat<T, 4> operator/(const quat<T, 4>& a, T1 b) {
    return {a.x / b, a.y / b, a.z / b, a.w / b};
}

/// Quaternion product.
template <typename T>
inline quat<T, 4> operator*(const quat<T, 4>& a, const quat<T, 4>& b) {
    return {a.x * b.w + a.w * b.x + a.y * b.w - a.z * b.y,
        a.y * b.w + a.w * b.y + a.z * b.x - a.x * b.z,
        a.z * b.w + a.w * b.z + a.x * b.y - a.y * b.x,
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z};
}

/// Quaternion conjugate.
template <typename T>
inline quat<T, 4> conjugate(const quat<T, 4>& v) {
    return {-v.x, -v.y, -v.z, v.w};
}
/// Quaternion inverse.
template <typename T>
inline quat<T, 4> inverse(const quat<T, 4>& v) {
    return conjugate(v) / dot(vec<T, 4>(v), vec<T, 4>(v));
}

/// Quaternion normalization.
template <typename T>
inline quat<T, 4> normalize(const quat<T, 4>& v) {
    auto l = length(vec<T, 4>{v.x, v.y, v.z, v.w});
    if (!l) return {0, 0, 0, 1};
    return {v.x / l, v.y / l, v.z / l, v.w / l};
}

/// Quaternion spherical linear interpolation.
template <typename T, typename T1>
inline quat<T, 4> slerp(const quat<T, 4>& a, const quat<T, 4>& b, T1 t) {
    return (quat<T, 4>)slerp(vec<T, 4>(a),
        dot(vec<T, 4>(a), vec<T, 4>(b)) < 0 ? -vec<T, 4>(b) : vec<T, 4>(b), t);
}

/// Stream write.
template <typename T, int N>
inline std::ostream& operator<<(std::ostream& os, const quat<T, N>& a) {
    for (auto i = 0; i < N; i++) {
        if (i) os << ' ';
        os << data(a)[i];
    }
    return os;
}
/// Stream read.
template <typename T, int N>
inline std::istream& operator>>(std::istream& is, quat<T, N>& a) {
    for (auto i = 0; i < N; i++) is >> data(a)[i];
    return is;
}

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// AXIS ALIGNED BOUNDING BOXES
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup bbox Axis-aligned bounding boxes
/// @{

/// Axis aligned bounding box represented as a min/max vector pairs.
template <typename T, int N>
struct bbox {
    /// Initializes an invalid bbox.
    bbox() : min{flt_max}, max{flt_min} {}
    /// Element constructor with min/max values.
    bbox(const vec<T, N>& m, const vec<T, N>& M) : min{m}, max{M} {}

    /// Element access.
    vec<T, N>& operator[](int i) { return (&min)[i]; }
    /// Element access.
    const vec<T, N>& operator[](int i) const { return (&min)[i]; }

    /// Minimum bounds.
    vec<T, N> min;
    /// Maximum bounds.
    vec<T, N> max;
};

/// 1-dimensional float bounding box.
using bbox1f = bbox<float, 1>;
/// 2-dimensional float bounding box.
using bbox2f = bbox<float, 2>;
/// 3-dimensional float bounding box.
using bbox3f = bbox<float, 3>;
/// 4-dimensional float bounding box.
using bbox4f = bbox<float, 4>;

/// 1-dimensional float empty bbox.
const auto invalid_bbox1f = bbox1f();
/// 2-dimensional float empty bbox.
const auto invalid_bbox2f = bbox2f();
/// 3-dimensional float empty bbox.
const auto invalid_bbox3f = bbox3f();
/// 4-dimensional float empty bbox.
const auto invalid_bbox4f = bbox4f();

/// Bounding box equality.
template <typename T, int N>
inline bool operator==(const bbox<T, N>& a, const bbox<T, N>& b) {
    return a.min == b.min && a.max == b.max;
}
/// Bounding box inequality.
template <typename T, int N>
inline bool operator!=(const bbox<T, N>& a, const bbox<T, N>& b) {
    return a.min != b.min || a.max != b.max;
}

/// Bounding box center.
template <typename T, int N>
inline vec<T, N> bbox_center(const bbox<T, N>& a) {
    return (a.min + a.max) / 2;
}
/// Bounding box diagonal.
template <typename T, int N>
inline vec<T, N> bbox_diagonal(const bbox<T, N>& a) {
    return a.max - a.min;
}

/// Expands a bounding box with a point.
template <typename T>
inline bbox<T, 1> expand(const bbox<T, 1>& a, T b) {
    return {{min(a.min.x, b)}, {max(a.max.x, b)}};
}
/// Expands a bounding box with a point.
template <typename T>
inline bbox<T, 1> expand(const bbox<T, 1>& a, const vec<T, 1>& b) {
    return {{min(a.min.x, b.x)}, {max(a.max.x, b.x)}};
}
/// Expands a bounding box with a point.
template <typename T>
inline bbox<T, 2> expand(const bbox<T, 2>& a, const vec<T, 2>& b) {
    return {{min(a.min.x, b.x), min(a.min.y, b.y)},
        {max(a.max.x, b.x), max(a.max.y, b.y)}};
}
/// Expands a bounding box with a point.
template <typename T>
inline bbox<T, 3> expand(const bbox<T, 3>& a, const vec<T, 3>& b) {
    return {{min(a.min.x, b.x), min(a.min.y, b.y), min(a.min.z, b.z)},
        {max(a.max.x, b.x), max(a.max.y, b.y), max(a.max.z, b.z)}};
}
/// Expands a bounding box with a point.
template <typename T>
inline bbox<T, 4> expand(const bbox<T, 4>& a, const vec<T, 4>& b) {
    return {{min(a.min.x, b.x), min(a.min.y, b.y), min(a.min.z, b.z),
                min(a.min.w, b.w)},
        {max(a.max.x, b.x), max(a.max.y, b.y), max(a.max.z, b.z),
            max(a.max.w, b.w)}};
}

/// Expands a bounding box with a bounding box.
template <typename T>
inline bbox<T, 1> expand(const bbox<T, 1>& a, const bbox<T, 1>& b) {
    return {min(a.min.x, b.min.x), max(a.max.x, b.max.x)};
}
/// Expands a bounding box with a bounding box.
template <typename T>
inline bbox<T, 2> expand(const bbox<T, 2>& a, const bbox<T, 2>& b) {
    return {{min(a.min.x, b.min.x), min(a.min.y, b.min.y)},
        {max(a.max.x, b.max.x), max(a.max.y, b.max.y)}};
}
/// Expands a bounding box with a bounding box.
template <typename T>
inline bbox<T, 3> expand(const bbox<T, 3>& a, const bbox<T, 3>& b) {
    return {
        {min(a.min.x, b.min.x), min(a.min.y, b.min.y), min(a.min.z, b.min.z)},
        {max(a.max.x, b.max.x), max(a.max.y, b.max.y), max(a.max.z, b.max.z)}};
}
/// Expands a bounding box with a bounding box.
template <typename T>
inline bbox<T, 4> expand(const bbox<T, 4>& a, const bbox<T, 4>& b) {
    return {{min(a.min.x, b.min.x), min(a.min.y, b.min.y),
                min(a.min.z, b.min.z), min(a.min.w, b.min.w)},
        {max(a.max.x, b.max.x), max(a.max.y, b.max.y), max(a.max.z, b.max.z),
            max(a.max.w, b.max.w)}};
}

/// Check if a bounding box contains a point.
template <typename T, int N>
inline bool contains(const bbox<T, N>& a, const vec<T, N>& b) {
    for (auto i = 0; i < N; i++)
        if (a.min[i] > b[i] || a.max[i] < b[i]) return false;
    return true;
}
/// Check if a bounding box contains a bounding box.
template <typename T, int N>
inline bool contains(const bbox<T, 3>& a, const bbox<T, 3>& b) {
    for (auto i = 0; i < N; i++)
        if (a.min[i] > b.max[i] || a.max[i] < b.min[i]) return false;
    return true;
}

/// Expands a bounding box with a point.
template <typename T, int N>
inline bbox<T, N>& operator+=(bbox<T, N>& a, const vec<T, N>& b) {
    return a = expand(a, b);
}
/// Expands a bounding box with a bounding box.
template <typename T, int N>
inline bbox<T, N>& operator+=(bbox<T, N>& a, const bbox<T, N>& b) {
    return a = expand(a, b);
}

/// Initialize a bonding box from a list of points.
template <typename T, int N>
inline bbox<T, N> make_bbox(int count, const vec<T, N>* v) {
    auto a = bbox<T, N>();
    for (auto j = 0; j < count; j++) a += v[j];
    return a;
}
/// Initialize a bonding box from a list of points.
template <typename T, int N>
inline bbox<T, N> make_bbox(const std::initializer_list<vec<T, N>>& v) {
    auto a = bbox<T, N>();
    for (auto&& vv : v) a += vv;
    return a;
}

/// Computes the corners of a bounding boxes.
template <typename T>
inline std::array<vec<T, 2>, 4> bbox_corners(const bbox<T, 2>& a) {
    return {{{a.min.x, a.min.y}, {a.min.x, a.min.y}, {a.min.x, a.max.y},
        {a.min.x, a.max.y}}};
}
/// Computes the corners of a bounding boxes.
template <typename T>
inline std::array<vec<T, 3>, 8> bbox_corners(const bbox<T, 3>& a) {
    return {{{a.min.x, a.min.y, a.min.z}, {a.min.x, a.min.y, a.max.z},
        {a.min.x, a.max.y, a.min.z}, {a.min.x, a.max.y, a.max.z},
        {a.max.x, a.min.y, a.min.z}, {a.max.x, a.min.y, a.max.z},
        {a.max.x, a.max.y, a.min.z}, {a.max.x, a.max.y, a.max.z}}};
}

/// Stream write.
template <typename T, int N>
inline std::ostream& operator<<(std::ostream& os, const bbox<T, N>& a) {
    return os << a.min << ' ' << a.max;
}
/// Stream read.
template <typename T, int N>
inline std::istream& operator>>(std::istream& is, bbox<T, N>& a) {
    return is >> a.min >> a.max;
}

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// PRIMITIVE BBOX FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup prim_bbox Primitive bounding boxes
/// @{

/// Point bounds.
template <typename T, typename T1>
inline bbox<T, 3> point_bbox(const vec<T, 3>& p, T1 r = 0) {
    return bbox<T, 3>{p - vec<T, 3>{r, r, r}, p + vec<T, 3>{r, r, r}};
}

/// Line bounds.
template <typename T, typename T1>
inline bbox<T, 3> line_bbox(
    const vec<T, 3>& v0, const vec<T, 3>& v1, T1 r0 = 0, T1 r1 = 0) {
    return make_bbox({v0 - vec<T, 3>{r0, r0, r0}, v0 + vec<T, 3>{r0, r0, r0},
        v1 - vec<T, 3>{r1, r1, r1}, v1 + vec<T, 3>{r1, r1, r1}});
}

/// Triangle bounds.
template <typename T>
inline bbox<T, 3> triangle_bbox(
    const vec<T, 3>& v0, const vec<T, 3>& v1, const vec<T, 3>& v2) {
    return make_bbox({v0, v1, v2});
}

/// Quad bounds.
template <typename T>
inline bbox<T, 3> quad_bbox(const vec<T, 3>& v0, const vec<T, 3>& v1,
    const vec<T, 3>& v2, const vec<T, 3>& v3) {
    return make_bbox({v0, v1, v2, v3});
}

/// Tetrahedron bounds.
template <typename T, typename T1>
inline bbox<T, 3> tetrahedron_bbox(const vec<T, 3>& v0, const vec<T, 3>& v1,
    const vec<T, 3>& v2, const vec<T, 3>& v3) {
    return make_bbox({v0, v1, v2, v3});
}

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// RAYS
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup ray Rays
/// @{

/// Rays with origin, direction and min/max t value. Origin and directions are
/// of N-elements.
template <typename T, int N>
struct ray {
    /// Default constructor. Initializes to an invalid ray.
    ray() : o{0}, d{0}, tmin{0}, tmax{flt_max} {}
    /// Initializes a ray from its elements.
    ray(const vec<T, N>& o, const vec<T, N>& d, T tmin = 0, T tmax = flt_max)
        : o(o), d(d), tmin(tmin), tmax(tmax) {}

    /// origin
    vec<T, N> o;
    /// direction
    vec<T, N> d;
    /// minimum distance
    T tmin;
    /// maximum distance
    T tmax;
};

/// 2-dimensional float ray.
using ray2f = ray<float, 2>;
/// 3-dimensional float ray.
using ray3f = ray<float, 3>;
/// 4-dimensional float ray.
using ray4f = ray<float, 4>;

/// Construct a ray using a default epsilon.
template <typename T, int N>
inline ray<T, N> make_ray(
    const vec<T, N>& o, const vec<T, N>& d, T eps = 1e-4f) {
    return ray<T, N>{o, d, eps, flt_max};
}
/// Construct a ray segment using a default epsilon.
template <typename T, int N>
inline ray<T, N> make_segment(
    const vec<T, N>& p1, const vec<T, N>& p2, T eps = 1e-4f) {
    return ray<T, N>{p1, normalize(p2 - p1), eps, length(p2 - p1) - 2 * eps};
}

/// Stream write.
template <typename T, int N>
inline std::ostream& operator<<(std::ostream& os, const ray<T, N>& a) {
    return os << a.o << ' ' << a.d << ' ' << a.tmin << ' ' << a.tmax;
}
/// Stream read.
template <typename T, int N>
inline std::istream& operator>>(std::istream& is, ray<T, N>& a) {
    return is >> a.o >> a.d >> a.tmin >> a.tmax;
}

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// TRANSFORMS
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup transform Transforms
/// @{

/// Transforms a point by a matrix.
template <typename T>
inline vec<T, 2> transform_point(const mat<T, 3>& a, const vec<T, 2>& b) {
    auto vb = vec<T, 2>{b.x, b.y, 1};
    auto tvb = a * vb;
    return vec<T, 2>{tvb.x, tvb.y} / tvb.w;
}
/// Transforms a point by a matrix.
template <typename T>
inline vec<T, 3> transform_point(const mat<T, 4>& a, const vec<T, 3>& b) {
    auto vb = vec<T, 4>{b.x, b.y, b.z, 1};
    auto tvb = a * vb;
    return vec<T, 3>{tvb.x, tvb.y, tvb.z} / tvb.w;
}
/// Transforms a vector by a matrix.
template <typename T>
inline vec<T, 2> transform_vector(const mat<T, 3>& a, const vec<T, 2>& b) {
    auto vb = vec<T, 2>{b.x, b.y, 0};
    auto tvb = a * vb;
    return vec<T, 2>{tvb.x, tvb.y} / tvb.z;
}
/// Transforms a vector by a matrix.
template <typename T>
inline vec<T, 3> transform_vector(const mat<T, 4>& a, const vec<T, 3>& b) {
    auto vb = vec<T, 4>{b.x, b.y, b.z, 0};
    auto tvb = a * vb;
    return vec<T, 3>{tvb.x, tvb.y, tvb.z};
}
/// Transforms a direction by a matrix.
template <typename T, int N>
inline vec<T, N> transform_direction(
    const mat<T, N + 1>& a, const vec<T, N>& b) {
    return normalize(transform_vector(a, b));
}
/// Transforms a ray by a matrix, leaving the direction not normalized.
template <typename T, int N>
inline ray<T, N> transform_ray(const mat<T, N + 1>& a, const ray<T, N>& b) {
    return {transform_point(a, b.o), transform_vector(a, b.d), b.tmin, b.tmax};
}
/// transforms a bbox by a matrix
template <typename T, int N>
inline bbox<T, N> transform_bbox(const mat<T, N + 1>& a, const bbox<T, N>& b) {
    auto corners = bbox_corners(b);
    auto xformed = bbox<T, N>();
    for (auto& corner : corners) xformed += transform_point(a, corner);
    return xformed;
}

/// Transforms a point by a frame, i.e. an affine transform.
template <typename T>
inline vec<T, 2> transform_point(const frame<T, 2>& a, const vec<T, 2>& b) {
    return a.x * b.x + a.y * b.y + a.o;
}
/// Transforms a point by a frame, i.e. an affine transform.
template <typename T>
inline vec<T, 3> transform_point(const frame<T, 3>& a, const vec<T, 3>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.o;
}
/// Transforms a vector by a frame, i.e. an affine transform.
template <typename T>
inline vec<T, 2> transform_vector(const frame<T, 2>& a, const vec<T, 2>& b) {
    return a.x * b.x + a.y * b.y;
}
/// Transforms a vector by a frame, i.e. an affine transform.
template <typename T>
inline vec<T, 3> transform_vector(const frame<T, 3>& a, const vec<T, 3>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
/// Transforms a direction by a frame, i.e. an affine transform.
template <typename T, int N>
inline vec<T, N> transform_direction(const frame<T, N>& a, const vec<T, N>& b) {
    return normalize(transform_vector(a, b));
}
/// Transforms a frame by a frame, i.e. an affine transform.
template <typename T, int N>
inline frame<T, N> transform_frame(const frame<T, N>& a, const frame<T, N>& b) {
    return {a.rot() * b.rot(), a.rot() * b.pos() + a.pos()};
}
/// Transforms a ray by a frame, i.e. an affine transform.
template <typename T, int N>
inline ray<T, N> transform_ray(const frame<T, 3>& a, const ray<T, N>& b) {
    return {
        transform_point(a, b.o), transform_direction(a, b.d), b.tmin, b.tmax};
}
/// Transforms a bbox by a frame, i.e. an affine transform.
template <typename T, int N>
inline bbox<T, N> transform_bbox(const frame<T, N>& a, const bbox<T, N>& b) {
    auto corners = bbox_corners(b);
    auto xformed = bbox<T, N>();
    for (auto& corner : corners) xformed += transform_point(a, corner);
    return xformed;
}
/// Transforms a bbox by a frame, i.e. an affine transform.
template <typename T>
inline bbox<T, 3> transform_bbox(const frame<T, 3>& a, const bbox<T, 3>& b) {
    // Code from Real-time Collision Detection by Christer Ericson Sect. 4.2.6
    // Transform AABB a by the matrix m and translation t,
    // find maximum extents, and store result into AABB b.
    // start by adding in translation
    auto c = bbox<T, 3>{a.o, a.o};
    auto rot = mat<T, 3>{a.x, a.y, a.z};
    // for all three axes
    for (auto i = 0; i < 3; i++) {
        // form extent by summing smaller and larger terms respectively
        for (auto j = 0; j < 3; j++) {
            auto e = rot[j][i] * b.min[j];
            auto f = rot[j][i] * b.max[j];
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
}

/// Inverse transforms a point by a frame, assuming a rigid transform.
template <typename T>
inline vec<T, 2> transform_point_inverse(
    const frame<T, 2>& a, const vec<T, 2>& b) {
    return {dot(b - a.o, a.x), dot(b - a.o, a.y)};
}
/// Inverse transforms a point by a frame, assuming a rigid transform.
template <typename T>
inline vec<T, 3> transform_point_inverse(
    const frame<T, 3>& a, const vec<T, 3>& b) {
    return {dot(b - a.o, a.x), dot(b - a.o, a.y), dot(b - a.o, a.z)};
}
/// Inverse transforms a vector by a frame, assuming a rigid transform.
template <typename T>
inline vec<T, 2> transform_vector_inverse(
    const frame<T, 2>& a, const vec<T, 2>& b) {
    return {dot(b, a.x), dot(b, a.y)};
}
/// Inverse transforms a vector by a frame, assuming a rigid transform.
template <typename T>
inline vec<T, 3> transform_vector_inverse(
    const frame<T, 3>& a, const vec<T, 3>& b) {
    return {dot(b, a.x), dot(b, a.y), dot(b, a.z)};
}
/// Inverse transforms a direction by a frame, assuming a rigid transform.
template <typename T, int N>
inline vec<T, N> transform_direction_inverse(
    const frame<T, N>& a, const vec<T, N>& b) {
    return normalize(transform_vector_inverse(a, b));
}
/// Inverse transforms a direction by a frame, assuming a rigid transform.
template <typename T, int N>
inline ray<T, N> transform_ray_inverse(
    const frame<T, N>& a, const ray<T, N>& b) {
    return {transform_point_inverse(a, b.o),
        transform_direction_inverse(a, b.d), b.tmin, b.tmax};
}
/// Inverse transforms a bbox by a frame, assuming a rigid transform.
template <typename T, int N>
inline bbox<T, N> transform_bbox_inverse(
    const frame<T, N>& a, const bbox<T, N>& b) {
    return transform_bbox(inverse(a), b);
}

/// Translation affine transform.
template <typename T>
inline frame<T, 3> translation_frame(const vec<T, 3>& a) {
    return {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, a};
}
/// Scaling affine transform; this is not rigid and here for symmatry of API.
template <typename T>
inline frame<T, 3> scaling_frame(const vec<T, 3>& a) {
    return {{a.x, 0, 0}, {0, a.y, 0}, {0, 0, a.z}, {0, 0, 0}};
}
/// Rotation affine transform.
template <typename T, typename T1>
inline frame<T, 3> rotation_frame(const vec<T, 3>& axis, T1 angle) {
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
/// Rotation affine transform.
template <typename T>
inline frame<T, 3> rotation_frame(const quat<T, 4>& v) {
    return {{v.w * v.w + v.x * v.x - v.y * v.y - v.z * v.z,
                (v.x * v.y + v.z * v.w) * 2, (v.z * v.x - v.y * v.w) * 2},
        {(v.x * v.y - v.z * v.w) * 2,
            v.w * v.w - v.x * v.x + v.y * v.y - v.z * v.z,
            (v.y * v.z + v.x * v.w) * 2},
        {(v.z * v.x + v.y * v.w) * 2, (v.y * v.z - v.x * v.w) * 2,
            v.w * v.w - v.x * v.x - v.y * v.y + v.z * v.z},
        {0, 0, 0}};
}
/// OpenGL lookat frame. Z-axis can be inverted with inv_xz.
template <typename T>
inline frame<T, 3> lookat_frame(const vec<T, 3>& eye, const vec<T, 3>& center,
    const vec<T, 3>& up, bool inv_xz = false) {
    auto w = normalize(eye - center);
    auto u = normalize(cross(up, w));
    auto v = normalize(cross(w, u));
    if (inv_xz) {
        w = -w;
        u = -u;
    }
    return {u, v, w, eye};
}

/// OpenGL frustum matrix.
template <typename T>
inline mat<T, 4> frustum_mat(T l, T r, T b, T t, T n, T f) {
    return {{2 * n / (r - l), 0, 0, 0}, {0, 2 * n / (t - b), 0, 0},
        {(r + l) / (r - l), (t + b) / (t - b), -(f + n) / (f - n), -1},
        {0, 0, -2 * f * n / (f - n), 0}};
}
/// OpenGL orthographic matrix.
template <typename T>
inline mat<T, 4> ortho_mat(T l, T r, T b, T t, T n, T f) {
    return {{2 / (r - l), 0, 0, 0}, {0, 2 / (t - b), 0, 0},
        {0, 0, -2 / (f - n), 0},
        {-(r + l) / (r - l), -(t + b) / (t - b), -(f + n) / (f - n), 1}};
}

/// OpenGL orthographic 2D matrix.
template <typename T>
inline mat<T, 4> ortho2d_mat(T left, T right, T bottom, T top) {
    return ortho_mat(left, right, bottom, top, (T)-1, (T)1);
}
/// OpenGL orthographic matrix.
template <typename T>
inline mat<T, 4> ortho_mat(T xmag, T ymag, T near, T far) {
    return {{1 / xmag, 0, 0, 0}, {0, 1 / ymag, 0, 0},
        {0, 0, 2 / (near - far), 0}, {0, 0, (far + near) / (near - far), 1}};
}

/// OpenGL perspective matrix.
template <typename T>
inline mat<T, 4> perspective_mat(T fovy, T aspect, T near, T far) {
    auto tg = tan(fovy / 2);
    return {{1 / (aspect * tg), 0, 0, 0}, {0, 1 / tg, 0, 0},
        {0, 0, (far + near) / (near - far), -1},
        {0, 0, 2 * far * near / (near - far), 0}};
}
/// OpenGL infinite perspective matrix.
template <typename T>
inline mat<T, 4> perspective_mat(T fovy, T aspect, T near) {
    auto tg = tan(fovy / 2);
    return {{1 / (aspect * tg), 0, 0, 0}, {0, 1 / tg, 0, 0}, {0, 0, -1, -1},
        {0, 0, 2 * near, 0}};
}

/// Rotation affine transform.
template <typename T>
inline std::pair<vec<T, 4>, T> rotation_axisangle(const quat<T, 4>& a) {
    return {normalize(vec<T, 3>{a.x, a.y, a.z}), 2 * acos(a.w)};
}
/// Axis-angle to quaternion conversion.
template <typename T>
inline quat<T, 4> rotation_quat(const vec<T, 3>& axis, T angle) {
    auto len = length(axis);
    if (!len) return {0, 0, 0, 1};
    return quat<T, 4>{sin(angle / 2) * axis.x / len,
        sin(angle / 2) * axis.y / len, sin(angle / 2) * axis.z / len,
        cos(angle / 2)};
}
/// Rotation matrix to quaternion conversion.
template <typename T>
inline quat<T, 4> rotation_quat(const mat<T, 3>& m_) {
    // http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm
    auto q = quat<T, 4>();
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

/// Decompose an affine matrix into translation, rotation, scale.
/// Assumes there is no shear.
template <typename T>
inline std::tuple<vec<T, 3>, mat<T, 3>, vec<T, 3>> decompose_frame(
    const frame<T, 3>& m) {
    return {m.o, {normalize(m.x), normalize(m.y), normalize(m.z)},
        {length(m.x), length(m.y), length(m.z)}};
}

/// Decompose an affine matrix into translation, rotation, scale.
/// Assumes there is no shear and the matrix is affine.
template <typename T>
inline std::tuple<vec<T, 3>, quat<T, 4>, vec<T, 3>> decompose_frame(
    const frame<T, 3>& m) {
    auto pos = vec<T, 3>();
    auto rot = mat<T, 3>();
    auto scl = vec<T, 3>();
    tie(pos, rot, scl) = decompose_frame(m);
    return {pos, rotation_quat(rot), scl};
}
/// Decompose an affine matrix into translation, rotation, scale.
/// Assumes there is no shear and the matrix is affine.
template <typename T>
inline frame<T, 4> compose_frame(const vec<T, 3>& translation,
    const mat<T, 3>& rotation, const vec<T, 3>& scale) {
    return translation_frame(translation) * scaling_frame(scale) *
           rotation_frame(rotation);
}
/// Decompose an affine matrix into translation, rotation, scale.
/// Assumes there is no shear and the matrix is affine.
template <typename T>
inline frame<T, 4> compose_frame(const vec<T, 3>& translation,
    const quat<T, 4>& rotation, const vec<T, 3>& scale) {
    return translation_frame(translation) * scaling_frame(scale) *
           rotation_frame(rotation);
}

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// UI UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup ui User interface utilities
/// @{

/// Turntable for UI navigation.
void camera_turntable(vec3f& from, vec3f& to, vec3f& up, const vec2f& rotate,
    float dolly, const vec2f& pan);
/// Turntable for UI navigation.
void camera_turntable(frame3f& frame, float& focus, const vec2f& rotate,
    float dolly, const vec2f& pan);
/// FPS camera for UI navigation.
void camera_fps(frame3f& frame, const vec3f& transl, const vec2f& rotate);

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// RANDOM NUMBER GENERATION
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup rng Random number generation
/// @{

/// PCG random numbers. A family of random number generators that supports
/// multiple sequences. From http://www.pcg-random.org/
struct rng_pcg32 {
    /// RNG state.
    uint64_t state = 0x853c49e6748fea9bULL;
    /// RNG sequence. Must be odd.
    uint64_t inc = 0xda3e39cb94b95bdbULL;
};

/// Next random number.
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
    uint32_t threshold = (~n + 1u) % n;
    while (true) {
        uint32_t r = advance_rng(rng);
        if (r >= threshold) return r % n;
    }
#endif
}

/// Next random float in [0,1).
inline float next_rand1f(rng_pcg32& rng) {
#if 1
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

/// Distance between random number generators.
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
        swap(vals[i], vals[next_rand1i(rng, (uint32_t)i)]);
}

/// Random shuffle of a sequence.
template <typename T>
inline void rng_shuffle(rng_pcg32& rng, std::vector<T>& vals) {
    rng_shuffle(rng, vals.data(), vals.size());
}

/// Random number generator equality.
inline bool operator==(const rng_pcg32& a, const rng_pcg32& b) {
    return a.state == b.state && a.inc == b.inc;
}
/// Random number generator inequality.
inline bool operator!=(const rng_pcg32& a, const rng_pcg32& b) {
    return a.state != b.state || a.inc != b.inc;
}

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// MONETACARLO SAMPLING FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup montecarlo Monte Carlo sampling
/// @{

/// Sample an hemispherical direction with uniform distribution.
inline vec3f sample_hemisphere(const vec2f& ruv) {
    auto z = ruv.y;
    auto r = sqrt(1 - z * z);
    auto phi = 2 * pif * ruv.x;
    return vec3f(r * cos(phi), r * sin(phi), z);
}
/// Pdf for uniform hemispherical direction.
inline float sample_hemisphere_pdf(const vec3f& w) {
    return (w.z <= 0) ? 0 : 1 / (2 * pif);
}

/// Sample a spherical direction with uniform distribution.
inline vec3f sample_sphere(const vec2f& ruv) {
    auto z = 2 * ruv.y - 1;
    auto r = sqrt(1 - z * z);
    auto phi = 2 * pif * ruv.x;
    return vec3f(r * cos(phi), r * sin(phi), z);
}
/// Pdf for uniform spherical direction.
inline float sample_sphere_pdf(const vec3f& w) { return 1 / (4 * pif); }

/// Sample an hemispherical direction with cosine distribution.
inline vec3f sample_hemisphere_cosine(const vec2f& ruv) {
    auto z = sqrt(ruv.y);
    auto r = sqrt(1 - z * z);
    auto phi = 2 * pif * ruv.x;
    return vec3f(r * cos(phi), r * sin(phi), z);
}
/// Pdf for cosine hemispherical direction.
inline float sample_hemisphere_cosine_pdf(const vec3f& w) {
    return (w.z <= 0) ? 0 : w.z / pif;
}

/// Sample an hemispherical direction with cosine power distribution.
inline vec3f sample_hemisphere_cospower(float n, const vec2f& ruv) {
    auto z = pow(ruv.y, 1 / (n + 1));
    auto r = sqrt(1 - z * z);
    auto phi = 2 * pif * ruv.x;
    return vec3f(r * cos(phi), r * sin(phi), z);
}
/// Pdf for cosine power hemispherical direction.
inline float sample_hemisphere_cospower_pdf(float n, const vec3f& w) {
    return (w.z <= 0) ? 0 : pow(w.z, n) * (n + 1) / (2 * pif);
}

/// Sample a point uniformly on a disk.
inline vec3f sample_disk(const vec2f& ruv) {
    auto r = sqrt(ruv.y);
    auto phi = 2 * pif * ruv.x;
    return {cos(phi) * r, sin(phi) * r, 0};
}
/// Pdf for uniform disk sampling.
inline float sample_disk_pdf() { return 1 / pif; }

/// Sample a point uniformly on a cylinder, without caps.
inline vec3f sample_cylinder(const vec2f& ruv) {
    auto phi = 2 * pif * ruv.x;
    return {sin(phi), cos(phi), ruv.y * 2 - 1};
}
/// Pdf for uniform cylinder sampling.
inline float sample_cylinder_pdf() { return 1 / pif; }

/// Sample a point uniformly on a triangle.
inline vec2f sample_triangle(const vec2f& ruv) {
    return {1 - sqrt(ruv.x), ruv.y * sqrt(ruv.x)};
}
/// Sample a point uniformly on a triangle.
inline vec3f sample_triangle(
    const vec3f& v0, const vec3f& v1, const vec3f& v2, const vec2f& ruv) {
    auto uv = sample_triangle(ruv);
    return v0 * (1 - uv.x - uv.y) + v1 * uv.x + v2 * uv.y;
}
/// Pdf for uniform triangle sampling, i.e. triangle area.
inline float sample_triangle_pdf(
    const vec3f& v0, const vec3f& v1, const vec3f& v2) {
    return 2 / length(cross(v1 - v0, v2 - v0));
}

/// Sample an index with uniform distribution.
inline int sample_index(int size, float r) {
    return clamp((int)(r * size), 0, size - 1);
}
/// Pdf for uniform index sampling.
inline float sample_index_pdf(int size) { return 1.0f / size; }

/// Sample a discrete distribution represented by its cdf.
inline int sample_discrete(const std::vector<float>& cdf, float r) {
    // todo: implement binary search better
    r = clamp(r * cdf.back(), 0.0f, cdf.back() - 0.00001f);
    for (auto i = 1; i < cdf.size(); i++) {
        if (cdf[i] > r) return i - 1;
    }
    return (int)cdf.size() - 1;
}
/// Pdf for uniform discrete distribution sampling.
inline float sample_discrete_pdf(const std::vector<float>& cdf, int idx) {
    if (idx == 0) return cdf.at(0);
    return cdf.at(idx) - cdf.at(idx - 1);
}

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// HASHING
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup hash Hashing
/// @{

/// Computes the i-th term of a permutation of l values keyed by p.
/// From Correlated Multi-Jittered Sampling by Kensler @ Pixar
inline uint32_t cmjs_permute(uint32_t i, uint32_t n, uint32_t key) {
    // clang-format off
    uint32_t w = n - 1;
    w |= w >> 1; w |= w >> 2; w |= w >> 4; w |= w >> 8; w |= w >> 16;
    do {
        i ^= key;
        i *= 0xe170893du; i ^= key >> 16;
        i ^= (i & w) >> 4; i ^= key >> 8;
        i *= 0x0929eb3f; i ^= key >> 23;
        i ^= (i & w) >> 1; i *= 1 | key >> 27;
        i *= 0x6935fa69; i ^= (i & w) >> 11;
        i *= 0x74dcb303; i ^= (i & w) >> 2;
        i *= 0x9e501cc3; i ^= (i & w) >> 2;
        i *= 0xc860a3df; i &= w; i ^= i >> 5;
    } while (i >= n);
    return (i + key) % n;
    // clang-format on
}

/// Computes a float value by hashing i with a key p.
/// From Correlated Multi-Jittered Sampling by Kensler @ Pixar
inline float cmjs_randfloat(uint32_t i, uint32_t key) {
    // clang-format off
    i ^= key;
    i ^= i >> 17; i ^= i >> 10; i *= 0xb36534e5;
    i ^= i >> 12; i ^= i >> 21; i *= 0x93fc4795;
    i ^= 0xdf6e307f; i ^= i >> 17; i *= 1 | key >> 18;
    return i * (1.0f / 4294967808.0f);
    // clang-format on
}

/// 32 bit integer hash. Public domain code.
inline uint32_t hash_uint32(uint64_t a) {
    // clang-format off
    a -= (a << 6); a ^= (a >> 17); a -= (a << 9); a ^= (a << 4);
    a -= (a << 3); a ^= (a << 10); a ^= (a >> 15);
    return a;
    // clang-format on
}

/// 64 bit integer hash. Public domain code.
inline uint64_t hash_uint64(uint64_t a) {
    // clang-format off
    a = (~a) + (a << 21); a ^= (a >> 24);
    a += (a << 3) + (a << 8); a ^= (a >> 14);
    a += (a << 2) + (a << 4); a ^= (a >> 28); a += (a << 31);
    return a;
    // clang-format on
}

/// 64-to-32 bit integer hash. Public domain code.
inline uint32_t hash_uint64_32(uint64_t a) {
    // clang-format off
    a = (~a) + (a << 18); a ^= (a >> 31); a *= 21;
    a ^= (a >> 11); a += (a << 6); a ^= (a >> 22);
    return (uint32_t)a;
    // clang-format on
}

/// Combines two 64 bit hashes as in boost::hash_combine.
inline size_t hash_combine(size_t a, size_t b) {
    return a ^ (b + 0x9e3779b9 + (a << 6) + (a >> 2));
}

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// PERLIN NOISE FUNCTION
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup noise Perlin noise
/// @{

// noise functions from stb_perlin.h

/// Compute the revised Pelin noise function. Wrap provides a wrapping noise
/// but must be power of two (wraps at 256 anyway). For octave based noise,
/// good values are obtained with octaves=6 (numerber of noise calls),
/// lacunarity=~2.0 (spacing between successive octaves: 2.0 for warpping
/// output), gain=0.5 (relative weighting applied to each successive octave),
/// offset=1.0 (used to invert the ridges).
float perlin_noise(const vec3f& p, const vec3i& wrap = zero3i);
/// Ridge noise function. See perlin_noise() for params.
float perlin_ridge_noise(const vec3f& p, float lacunarity = 2.0f,
    float gain = 0.5f, float offset = 1.0f, int octaves = 6,
    const vec3i& wrap = zero3i);
/// Fractal brownian motion noise. See perlin_noise() for params.
float perlin_fbm_noise(const vec3f& p, float lacunarity = 2.0f,
    float gain = 0.5f, int octaves = 6, const vec3i& wrap = zero3i);
/// Fractal turbulence noise. See perlin_noise() for params.
float perlin_turbulence_noise(const vec3f& p, float lacunarity = 2.0f,
    float gain = 0.5f, int octaves = 6, const vec3i& wrap = zero3i);

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// PYTHON-LIKE ITERATORS
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup python_iter Python-like iterators
/// @{

/// Implementation of Python-like range generator. Create it with the the
/// `range()` functions to use argument deduction.
struct range_generator {
    // clang-format off
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
    // clang-format on
};

/// Python-like range ierator.
inline range_generator range(int max) { return {0, max, 1}; }
/// Python-like range ierator.
inline range_generator range(int min, int max, int step = 1) {
    return {min, max, step};
}

/// Implemenetation of Python-like enumerate. Create it with the function
/// `enumerate()` to use argument deduction.
template <typename T>
struct enumerate_generator {
    // clang-format off
    struct iterator {
        iterator(int pos, T* data) : _pos(pos), _data(data) {}
        bool operator!=(const iterator& a) const { return _pos < a._pos; }
        iterator& operator++() { _pos += 1; _data += 1; return *this; }
        std::pair<int, T&> operator*() const { return {_pos, *_data}; }
        private: int _pos; T* _data;
    };
    enumerate_generator(int max, T* data) : _max(max), _data(data) {}
    iterator begin() const { return {0, _data}; }
    iterator end() const { return {_max, _data + _max}; }
    private: int _max; T* _data;
    // clang-format on
};

/// Python-like enumerate.
template <typename T>
inline enumerate_generator<const T> enumerate(const std::vector<T>& vv) {
    return {(int)vv.size(), vv.data()};
}
/// Python-like enumerate.
template <typename T>
inline enumerate_generator<T> enumerate(std::vector<T>& vv) {
    return {(int)vv.size(), vv.data()};
}

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// CONTAINER UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup container_ops Container operations
/// @{

/// Append a vector to a vector.
template <typename T>
inline void append(std::vector<T>& v, const std::vector<T>& vv) {
    v.insert(v.end(), vv.begin(), vv.end());
}
/// Append two vectors.
template <typename T>
inline std::vector<T> join(const std::vector<T>& a, const std::vector<T>& b) {
    auto v = std::vector<T>();
    append(v, a);
    append(v, b);
    return v;
}

/// Get a key from a list of key-value pairs and its value.
template <typename K, typename V>
inline K get_key(const std::vector<std::pair<K, V>>& kvs, const V& v) {
    for (auto& kv : kvs)
        if (kv.second == v) return kv.first;
    throw std::runtime_error("key not found");
}
/// Get a value from a list of key-value pairs and its key.
template <typename K, typename V>
inline V get_value(const std::vector<std::pair<K, V>>& kvs, const K& k) {
    for (auto& kv : kvs)
        if (kv.first == k) return kv.second;
    throw std::runtime_error("key not found");
}

/// Find the position of a value in an array. Returns -1 if not found.
/// Wrapper for std::find().
template <typename T>
inline int find(const std::vector<T>& v, const T& vv) {
    auto pos = find(v.begin(), v.end(), vv);
    if (pos == v.end()) return -1;
    return pos = v.begin();
}

/// Find the first array value that is greater than the argument.
/// Assumes that the array is sorted. Wrapper for std::upper_bound().
template <typename T>
inline int upper_bound(const std::vector<T>& v, const T& vv) {
    for (auto i = 0; i < v.size(); i++)
        if (v[i] > vv) return i;
    return (int)v.size();
}

/// Find the first array value smaller that is greater or equal to the argument.
/// Assumes that the array is sorted. Wrapper for std::lower_bound().
template <typename T>
inline int lower_bound(const std::vector<T>& v, const T& vv) {
    for (auto i = 0; i < v.size(); i++)
        if (v[i] >= vv) return i;
    return (int)v.size();
}

/// Checks if a containers contains a value.
template <typename T>
inline bool contains(const std::vector<T>& v, const T& vv) {
    return find(v.begin(), v.end(), vv) != v.end();
}
/// Checks if a containers contains a value.
template <typename K, typename V>
inline bool contains(const std::map<K, V>& v, const K& vv) {
    return v.find(vv) != v.end();
}
/// Checks if a containers contains a value.
template <typename K, typename V>
inline bool contains(const std::unordered_map<K, V>& v, const K& vv) {
    return v.find(vv) != v.end();
}
/// Checks if a containers contains a value.
template <typename K, typename V>
inline bool contains(const std::set<K, V>& v, const K& vv) {
    return v.find(vv) != v.end();
}
/// Checks if a containers contains a value.
template <typename K, typename V>
inline bool contains(const std::unordered_set<K, V>& v, const K& vv) {
    return v.find(vv) != v.end();
}
/// Checks if a containers contains a value.
template <typename K, typename V, typename K1>
inline bool contains(const std::unordered_map<K, V>& v, const K1& vv) {
    return v.find(vv) != v.end();
}
/// Checks if a containers contains a value.
template <typename K, typename V, typename K1>
inline bool contains(const std::unordered_set<K, V>& v, const K1& vv) {
    return v.find(vv) != v.end();
}

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// TYPE SUPPORT
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup type Type support
/// @{

// Implementation from
// https://stackoverflow.com/questions/3926637/how-to-intentionally-cause-a-compile-time-error-on-template-instantiation
template <typename T>
struct _always_false {
    enum { value = false };
};

/// Names of enum values. Specialized by enums that support reflection.
template <typename T>
inline const std::vector<std::pair<std::string, T>>& enum_names() {
    static_assert(_always_false<T>::value, "Specialize this function.");
}

/// Names of enum values.
template <typename T>
inline const std::vector<std::pair<std::string, T>>& enum_names(T v) {
    return enum_names<T>();
}

/// Stream write.
template <typename T,
    typename std::enable_if<std::is_enum<T>::value, int>::type = 0>
inline std::ostream& operator<<(std::ostream& os, const T& a) {
    return os << get_key(enum_names(a), a);
}
/// Stream read.
template <typename T,
    typename std::enable_if<std::is_enum<T>::value, int>::type = 0>
inline std::istream& operator>>(std::istream& is, T& a) {
    auto str = std::string();
    is >> str;
    a = get_val(enum_names(a), str);
    return is;
}

/// Types of variable semantic
enum struct visit_sem_type {
    /// Generic value.
    value = 0,
    /// Name.
    name = 1,
    /// Path.
    path = 2,
    /// Object.
    object = 3,
    /// Reference.
    reference = 4,
    /// Color.
    color = 5,
};

/// Semantic for reflected values
struct visit_sem {
    /// Type.
    visit_sem_type type = visit_sem_type::value;
    /// Minimum value for numeric types.
    float min = 0;
    /// Maximum value for numeric types.
    float max = 0;
};

/// Type trait to enable visitors.
template <class T>
struct has_visitor : std::false_type {};

/// Visit struct elements. Calls `visitor(name,val.var,sem)` for each variable
/// of a structure, where `name` is the name of the variable, `var` is the
/// variable and `sem` is one a `visit_sem` value.
/// Implemented by structures that support reflection.
template <typename T, typename Visitor>
inline void visit(T& val, Visitor&& visitor) {
    static_assert(_always_false<T>::value, "Override this function.");
}

/// Visit pointer elements.
template <typename T, typename Visitor>
inline void visit(T*& val, Visitor&& visitor) {
    if (val) visit(*val, visitor);
}

/// Stream write.
template <typename T,
    typename std::enable_if<has_visitor<T>::value, int>::type = 0>
inline std::ostream& operator<<(std::ostream& os, const T& a) {
    visit(a, [&os](auto& name, auto& val, auto&) {
        os << name << ": " << val << "\n";
    });
}

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// GEOMETRY UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup geom Geometry utilities
/// @{

/// Line tangent.
inline vec3f line_tangent(const vec3f& v0, const vec3f& v1) {
    return normalize(v1 - v0);
}
/// Line length.
inline float line_length(const vec3f& v0, const vec3f& v1) {
    return length(v1 - v0);
}

/// Triangle normal.
inline vec3f triangle_normal(
    const vec3f& v0, const vec3f& v1, const vec3f& v2) {
    return normalize(cross(v1 - v0, v2 - v0));
}
/// Triangle area.
inline float triangle_area(const vec3f& v0, const vec3f& v1, const vec3f& v2) {
    return length(cross(v1 - v0, v2 - v0)) / 2;
}

/// Quad normal.
inline vec3f quad_normal(
    const vec3f& v0, const vec3f& v1, const vec3f& v2, const vec3f& v3) {
    return normalize(triangle_normal(v0, v1, v3) + triangle_normal(v3, v2, v1));
}
/// Quad area.
inline float quad_area(
    const vec3f& v0, const vec3f& v1, const vec3f& v2, const vec3f& v3) {
    return triangle_area(v0, v1, v3) + triangle_area(v3, v2, v1);
}

/// tetrahedron volume
inline float tetrahedron_volume(
    const vec3f& v0, const vec3f& v1, const vec3f& v2, const vec3f& v3) {
    return dot(cross(v1 - v0, v2 - v0), v3 - v0) / 6;
}

/// Triangle tangent and bitangent from uv (not othornormalized with themselfves
/// not the normal). Follows the definition in
/// http://www.terathon.com/code/tangent.html and
/// https://gist.github.com/aras-p/2843984
inline std::pair<vec3f, vec3f> triangle_tangents_fromuv(const vec3f& v0,
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

/// Copies of point value. Here only for completeness.
template <typename T>
inline T interpolate_point(const std::vector<T>& vals, int p) {
    if (vals.empty()) return T();
    return vals[p];
}

/// Interpolates values over a line parametrized from a to b by u. Same as lerp.
template <typename T, typename T1>
inline T interpolate_line(const T& v0, const T& v1, const T1 u) {
    return v0 * (1 - u) + v1 * u;
}
/// Interpolates values over lines parametrized from a to b by u. Same as lerp.
template <typename T, typename T1>
inline T interpolate_line(const std::vector<T>& vals, const vec2i& l, T1 u) {
    if (vals.empty()) return T();
    return vals[l.x] * (1 - u) + vals[l.y] * u;
}

/// Interpolates values over a triangle parametrized by u and v along the
/// (v1-v0) and (v2-v0) directions. Same as barycentric interpolation.
template <typename T, typename T1>
inline T interpolate_triangle(
    const T& v0, const T& v1, const T& v2, const vec<T1, 2>& uv) {
    return v0 * (1 - uv.x - uv.y) + v1 * uv.x + v2 * uv.y;
}
/// Interpolates values over triangles parametrized by u and v along the
/// (v1-v0) and (v2-v0) directions. Same as barycentric interpolation.
template <typename T, typename T1>
inline T interpolate_triangle(
    const std::vector<T>& vals, const vec3i& t, const vec<T1, 2>& uv) {
    if (vals.empty()) return T();
    return vals[t.x] * (1 - uv.x - uv.y) + vals[t.y] * uv.x + vals[t.z] * uv.y;
}

/// Interpolates values over a quad parametrized by u and v along the
/// (v1-v0) and (v2-v1) directions. Same as bilear interpolation.
template <typename T, typename T1>
inline T interpolate_triangle(
    const T& v0, const T& v1, const T& v2, const T& v3, const vec<T1, 2>& uv) {
    return v0 * (1 - uv.x) * (1 - uv.y) + v1 * uv.x * (1 - uv.y) +
           v2 * uv.x * uv.y + v3 * (1 - uv.x) * uv.y;
}
/// Interpolates values over quads parametrized by u and v along the
/// (v1-v0) and (v2-v1) direction. Same as bilear interpolation.
template <typename T, typename T1>
inline T interpolate_quad(
    const std::vector<T>& vals, const vec4i& t, const vec<T1, 2>& uv) {
    if (vals.empty()) return T();
    return vals[t.x] * (1 - uv.x) * (1 - uv.y) + vals[t.y] * uv.x * (1 - uv.y) +
           vals[t.z] * uv.x * uv.y + vals[t.w] * (1 - uv.x) * uv.y;
}

/// Evaluates the i-th Bernstein polynomial of degree degree at u.
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

/// Evaluates the derivative of i-th Bernstein polynomial of degree degree at u.
template <typename T>
inline T eval_bernstein_derivative(T u, int i, int degree) {
    return degree * (eval_bernstein(u, i - 1, degree - 1) -
                        eval_bernstein(u, i, degree - 1));
}

/// Interpolates values along a cubic Bezier segment parametrized by u.
template <typename T, typename T1>
inline T interpolate_bezier(
    const T& v0, const T& v1, const T& v2, const T& v3, T1 u) {
    return v0 * (1 - u) * (1 - u) * (1 - u) + v1 * 3 * u * (1 - u) * (1 - u) +
           v2 * 3 * u * u * (1 - u) + v3 * u * u * u;
}
/// Interpolates values along cubic Bezier segments parametrized by u.
template <typename T, typename T1>
inline T interpolate_bezier(const std::vector<T>& vals, const vec4i& b, T1 u) {
    if (vals.empty()) return T();
    return eval_bezier(vals[b.x], vals[b.y], vals[b.z], vals[b.w], u);
}
/// Computes the derivative of a cubic Bezier segment parametrized by u.
template <typename T, typename T1>
inline T interpolate_bezier_derivative(
    const T& v0, const T& v1, const T& v2, const T& v3, T1 u) {
    return (v1 - v0) * 3 * (1 - u) * (1 - u) + (v2 - v1) * 6 * u * (1 - u) +
           (v3 - v2) * 3 * u * u;
}
/// Computes the derivative of a cubic Bezier segments parametrized by u.
template <typename T, typename T1>
inline T interpolate_bezier_derivative(
    const std::vector<T>& vals, const vec4i& b, T1 u) {
    if (vals.empty()) return T();
    return interpolate_bezier_derivative(
        vals[b.x], vals[b.y], vals[b.z], vals[b.w], u);
}

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// ANIMATION UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup anim Animation utilities
/// @{

/// Evalautes a keyframed value using step interpolation.
template <typename T>
inline T eval_keyframed_step(
    const std::vector<float>& times, const std::vector<T>& vals, float time) {
    if (time <= times.front()) return vals.front();
    if (time >= times.back()) return vals.back();
    time = clamp(time, times.front(), times.back() - 0.001f);
    auto idx = upper_bound(times, time);
    return vals.at(idx - 1);
}

// Implementation detail.
template <typename T>
inline T eval_keyframed_lerp(const T& a, const T& b, float t) {
    return lerp(a, b, t);
}
template <typename T>
inline quat<T, 4> eval_keyframed_lerp(
    const quat<T, 4>& a, const quat<T, 4>& b, float t) {
    return slerp(a, b, t);
}

/// Evalautes a keyframed value using linear interpolation.
template <typename T>
inline T eval_keyframed_linear(
    const std::vector<float>& times, const std::vector<T>& vals, float time) {
    if (time <= times.front()) return vals.front();
    if (time >= times.back()) return vals.back();
    time = clamp(time, times.front(), times.back() - 0.001f);
    auto idx = upper_bound(times, time);
    auto t = (time - times.at(idx - 1)) / (times.at(idx) - times.at(idx - 1));
    return eval_keyframed_lerp(vals.at(idx - 1), vals.at(idx), t);
}

// Implementation detail.
template <typename T>
inline T eval_keyframed_cubic(
    const T& a, const T& b, const T& c, const T& d, float t) {
    return interpolate_bezier(a, b, c, d, t);
}
template <typename T>
inline quat<T, 4> eval_keyframed_cubic(const quat<T, 4>& a, const quat<T, 4>& b,
    const quat<T, 4>& c, const quat<T, 4>& d, float t) {
    return normalize(interpolate_bezier(a, b, c, d, t));
}

/// Evalautes a keyframed value using Bezier interpolation.
template <typename T>
inline T eval_keyframed_bezier(
    const std::vector<float>& times, const std::vector<T>& vals, float time) {
    if (time <= times.front()) return vals.front();
    if (time >= times.back()) return vals.back();
    time = clamp(time, times.front(), times.back() - 0.001f);
    auto idx = upper_bound(times, time);
    auto t = (time - times.at(idx - 1)) / (times.at(idx) - times.at(idx - 1));
    return eval_keyframed_cubic(
        vals.at(idx - 3), vals.at(idx - 2), vals.at(idx - 1), vals.at(idx), t);
}

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// SHAPE UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup shape Shape utilities
/// @{

/// Compute per-vertex normals/tangents for lines, triangles and quads with
/// positions pos. Weighted indicated whether the normals/tangents are
/// weighted by length/area.
std::vector<vec3f> compute_normals(const std::vector<vec2i>& lines,
    const std::vector<vec3i>& triangles, const std::vector<vec4i>& quads,
    const std::vector<vec3f>& pos, bool weighted = true);

/// Compute per-vertex tangent frames for triangle meshes.
/// Tangent space is defined by a four component vector.
/// The first three components are the tangent with respect to the u texcoord.
/// The fourth component is the sign of the tangent wrt the v texcoord.
/// Tangent frame is useful in normal mapping.
std::vector<vec4f> compute_tangent_frames(const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, const std::vector<vec3f>& norm,
    const std::vector<vec2f>& texcoord, bool weighted = true);

/// Apply skinning to vertex position and normals.
void compute_skinning(const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec4f>& weights,
    const std::vector<vec4i>& joints, const std::vector<mat4f>& xforms,
    std::vector<vec3f>& skinned_pos, std::vector<vec3f>& skinned_norm);
/// Apply skinning.
void compute_skinning(const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec4f>& weights,
    const std::vector<vec4i>& joints, const std::vector<frame3f>& xforms,
    std::vector<vec3f>& skinned_pos, std::vector<vec3f>& skinned_norm);
/// Apply skinning as specified in Khronos glTF.
void compute_matrix_skinning(const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec4f>& weights,
    const std::vector<vec4i>& joints, const std::vector<mat4f>& xforms,
    std::vector<vec3f>& skinned_pos, std::vector<vec3f>& skinned_norm);

/// Create an array of edges.
std::vector<vec2i> get_edges(const std::vector<vec2i>& lines,
    const std::vector<vec3i>& triangles, const std::vector<vec4i>& quads);
/// Create an array of boundary edges. Lines are always considered boundaries.
std::vector<vec2i> get_boundary_edges(const std::vector<vec2i>& lines,
    const std::vector<vec3i>& triangles, const std::vector<vec4i>& quads);

/// Get a list of all unique vertices.
std::vector<int> get_verts(const std::vector<vec2i>& lines,
    const std::vector<vec3i>& triangles, const std::vector<vec4i>& quads);
/// Create an array of boundary vertices. Lines are always considered
/// boundaries.
std::vector<int> get_boundary_verts(const std::vector<vec2i>& lines,
    const std::vector<vec3i>& triangles, const std::vector<vec4i>& quads);

/// Convert quads to triangles
std::vector<vec3i> convert_quads_to_triangles(const std::vector<vec4i>& quads);
/// Convert quads to triangles with a diamond-like topology.
/// Quads have to be consecutive one row after another.
std::vector<vec3i> convert_quads_to_triangles(
    const std::vector<vec4i>& quads, int row_length);

/// Convert beziers to lines using 3 lines for each bezier.
std::vector<vec2i> convert_bezier_to_lines(const std::vector<vec4i>& beziers);

/// Convert face-varying data to single primitives. Returns the quads indices
/// and filled vectors for pos, norm and texcoord.
std::tuple<std::vector<vec4i>, std::vector<vec3f>, std::vector<vec3f>,
    std::vector<vec2f>>
convert_face_varying(const std::vector<vec4i>& quads_pos,
    const std::vector<vec4i>& quads_norm,
    const std::vector<vec4i>& quads_texcoord, const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec2f>& texcoord);

/// Tesselate lines, triangles and quads by splitting edges and faces for quads.
/// Returns the tesselated elements and edge and face array for vertex
/// calculations.
std::tuple<std::vector<vec2i>, std::vector<vec3i>, std::vector<vec4i>,
    std::vector<vec2i>, std::vector<vec4i>>
subdivide_elems_linear(const std::vector<vec2i>& lines,
    const std::vector<vec3i>& triangles, const std::vector<vec4i>& quads,
    int nverts);

/// Subdivide vertex properties for edges and faces. This is implemented for
/// vecs and floats.
template <typename T>
std::vector<T> subdivide_vert_linear(const std::vector<T>& vert,
    const std::vector<vec2i>& edges, const std::vector<vec4i>& faces,
    bool normalized = false);

/// Performs the smoothing step of Catmull-Clark. Start with a tesselate quad
/// mesh obtained with subdivide_elems_linear() and subdivide_vert_linear(). To
/// handle open meshes with boundary, get the boundary from make_boundary_edge()
/// and pass it as crease_lines. To fix the boundary entirely, just get the
/// boundary vertices and pass it as creases. This is implemented for vecs
/// and floats.
template <typename T>
std::vector<T> subdivide_vert_catmullclark(const std::vector<vec4i>& quads,
    const std::vector<T>& vert, const std::vector<vec2i>& crease_tlines,
    const std::vector<int>& crease_tpoints, bool normalized = false);

/// Subdivide Bezier segments recursively by splitting each segment into two
/// in the middle. Returns the tesselated elements and dictionaries for vertex
/// calculations.
std::tuple<std::vector<vec4i>, std::vector<int>, std::vector<vec4i>>
subdivide_bezier_recursive(const std::vector<vec4i>& beziers, int nverts);
/// Subdivide vertex properties for Bezier subdivision. This is implemente for
/// vecs and floats.
template <typename T>
std::vector<T> subdivide_vert_bezier(const std::vector<T>& vert,
    const std::vector<int>& verts, const std::vector<vec4i>& segments,
    bool normalized = false);

/// Generate a rectangular grid of usteps x vsteps uv values for parametric
/// surface generation. Values cam wrap and have poles.
std::tuple<std::vector<vec4i>, std::vector<vec2f>> make_uvquads(int usteps,
    int vsteps, bool uwrap = false, bool vwrap = false, bool vpole0 = false,
    bool vpole1 = false);

/// Generate parametric num lines of usteps segments.
std::tuple<std::vector<vec2i>, std::vector<vec2f>> make_uvlines(
    int num, int usteps);

/// Generate a parametric point set. Mostly here for completeness.
std::tuple<std::vector<int>, std::vector<vec2f>> make_uvpoints(int num);

/// Merge elements between shapes. The elements are merged by increasing the
/// array size of the second array by the number of vertices of the first.
/// Vertex data can then be concatenated successfully.
std::tuple<std::vector<vec2i>, std::vector<vec3i>, std::vector<vec4i>>
merge_elems(int nverts, const std::vector<vec2i>& lines1,
    const std::vector<vec3i>& triangles1, const std::vector<vec4i>& quads1,
    const std::vector<vec2i>& lines2, const std::vector<vec3i>& triangles2,
    const std::vector<vec4i>& quads2);

/// Unshare shape data by duplicating all vertex data for each element,
/// giving a faceted look. Note that faceted tangents are not computed.
/// Returns the indices to copy from.
std::tuple<std::vector<vec2i>, std::vector<vec3i>, std::vector<vec4i>,
    std::vector<int>>
facet_elems(const std::vector<vec2i>& lines,
    const std::vector<vec3i>& triangles, const std::vector<vec4i>& quads);

/// Unshare vertices for faceting. This is implemented for vec and float types.
template <typename T>
std::vector<T> facet_vert(
    const std::vector<T>& vert, const std::vector<int>& vmap);

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// SHAPE SAMPLING
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup shape_sampling Shape sampling
/// @{

/// Pick a point.
inline int sample_points(int npoints, float re) {
    return sample_index(npoints, re);
}
/// Compute a distribution for sampling points uniformly.
std::vector<float> sample_points_cdf(int npoints);
/// Pick a point uniformly.
inline int sample_points(const std::vector<float>& cdf, float re) {
    return sample_discrete(cdf, re);
}

/// Compute a distribution for sampling lines uniformly.
std::vector<float> sample_lines_cdf(
    const std::vector<vec2i>& lines, const std::vector<vec3f>& pos);
/// Pick a point on lines uniformly.
inline std::pair<int, float> sample_lines(
    const std::vector<float>& cdf, float re, float ru) {
    return {sample_discrete(cdf, re), ru};
}

/// Compute a distribution for sampling triangle meshes uniformly.
std::vector<float> sample_triangles_cdf(
    const std::vector<vec3i>& triangles, const std::vector<vec3f>& pos);
/// Pick a point on a triangle mesh uniformly.
inline std::pair<int, vec2f> sample_triangles(
    const std::vector<float>& cdf, float re, const vec2f& ruv) {
    return {sample_discrete(cdf, re), sample_triangle(ruv)};
}

/// Compute a distribution for sampling quad meshes uniformly.
std::vector<float> sample_quads_cdf(
    const std::vector<vec4i>& quads, const std::vector<vec3f>& pos);
/// Pick a point on a quad mesh uniformly.
inline std::pair<int, vec2f> sample_quads(
    const std::vector<float>& cdf, float re, const vec2f& ruv) {
    return {sample_discrete(cdf, re), ruv};
}

/// Samples a set of points over a triangle mesh uniformly. Returns pos, norm
/// and tecoord of the sampled points.
std::tuple<std::vector<vec3f>, std::vector<vec3f>, std::vector<vec2f>>
sample_triangles_points(const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, const std::vector<vec3f>& norm,
    const std::vector<vec2f>& texcoord, int npoints, uint64_t seed = 0);

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// EXAMPLE SHAPES
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup shape_example Example shapes
/// @{

/// Make a sphere. Returns quads, pos.
std::tuple<std::vector<vec3i>, std::vector<vec3f>> make_sphere(int tesselation);

/// Make a geodesic sphere. Returns quads, pos.
std::tuple<std::vector<vec3i>, std::vector<vec3f>> make_geodesicsphere(
    int tesselation);

/// Make a cube with unique vertices. This is watertight but has no
/// texture coordinates or normals. Returns quads, pos.
std::tuple<std::vector<vec4i>, std::vector<vec3f>> make_cube(int tesselation);

/// Make a sphere. This is not watertight. Returns quads, pos, norm, texcoord.
std::tuple<std::vector<vec4i>, std::vector<vec3f>, std::vector<vec3f>,
    std::vector<vec2f>>
make_uvsphere(int tesselation, bool flipped = false);

/// Make a sphere. This is not watertight. Returns quads, pos, norm, texcoord.
std::tuple<std::vector<vec4i>, std::vector<vec3f>, std::vector<vec3f>,
    std::vector<vec2f>>
make_uvhemisphere(int tesselation, bool flipped = false);

/// Make a quad. Returns quads, pos, norm, texcoord.
std::tuple<std::vector<vec4i>, std::vector<vec3f>, std::vector<vec3f>,
    std::vector<vec2f>>
make_uvquad(int tesselation);

/// Make a facevarying sphere with unique vertices but different texture
/// coordinates. Returns (quads, pos), (quads, norm), (quads, texcoord).
std::tuple<std::vector<vec4i>, std::vector<vec3f>, std::vector<vec4i>,
    std::vector<vec3f>, std::vector<vec4i>, std::vector<vec2f>>
make_fvsphere(int tesselation);

/// Make a facevarying cube with unique vertices but different texture
/// coordinates. Returns (quads, pos), (quads, norm), (quads, texcoord).
std::tuple<std::vector<vec4i>, std::vector<vec3f>, std::vector<vec4i>,
    std::vector<vec3f>, std::vector<vec4i>, std::vector<vec2f>>
make_fvcube(int tesselation);

/// Make a suzanne monkey model for testing. Note that some quads are
/// degenerate. Returns quads, pos.
std::tuple<std::vector<vec4i>, std::vector<vec3f>> make_suzanne(
    int tesselation);

/// Make a cube. This is not watertight. Returns quads, pos, norm, texcoord.
std::tuple<std::vector<vec4i>, std::vector<vec3f>, std::vector<vec3f>,
    std::vector<vec2f>>
make_uvcube(int tesselation);

/// Make a sphere from a cube. This is not watertight. Returns quads, pos, norm,
/// texcoord.
std::tuple<std::vector<vec4i>, std::vector<vec3f>, std::vector<vec3f>,
    std::vector<vec2f>>
make_uvspherecube(int tesselation);

/// Make a cube than stretch it towards a sphere. This is not watertight.
/// Returns quads, pos, norm, texcoord.
std::tuple<std::vector<vec4i>, std::vector<vec3f>, std::vector<vec3f>,
    std::vector<vec2f>>
make_uvspherizedcube(int tesselation, float radius);

/// Make a sphere with caps flipped. This is not watertight. Returns quads, pos,
/// norm, texcoord.
std::tuple<std::vector<vec4i>, std::vector<vec3f>, std::vector<vec3f>,
    std::vector<vec2f>>
make_uvflipcapsphere(int tesselation, float z, bool flipped = false);

/// Make a cutout sphere. This is not watertight. Returns quads, pos, norm,
/// texcoord.
std::tuple<std::vector<vec4i>, std::vector<vec3f>, std::vector<vec3f>,
    std::vector<vec2f>>
make_uvcutsphere(int tesselation, float z, bool flipped = false);

/// Make seashell params
struct make_seashell_params {
    /// Spiral revolutions.
    float spiral_revolutions = 2;
    /// Spiral angle (alpha) in [0,pi].
    float spiral_angle = 83 * pif / 180;
    /// Enlarging revolutions (beta) in [0,2pi].
    float enlarging_angle = 42 * pif / 180;
    /// Spiral aperture (A) in [0,inf].
    float spiral_aperture = 0.25f;
    /// Ellipse axis (a,b) in [0,inf].
    vec2f ellipse_axis = {0.12f, 0.20f};
    /// Curve rotatation (psi, Omega, mu) in [0,2pi].
    vec3f curve_rotation = {70 * pif / 180, 30 * pif / 180, 10 * pif / 180};
    /// Number of nodules (N) in [0,ing].
    float nodules_num = 0;
    /// Length of nodules along curve and spiral (W1,W2) in [0,inf].
    vec2f nodule_length = {0, 0};
    /// Height of nodules (L) in [0,inf].
    float nodule_height = 0;
    /// Position of nodules (P) in [0,inf].
    float nodule_pos = 0;
};

/// Make a seashell. This is not watertight. Returns quads, pos, norm, texcoord.
std::tuple<std::vector<vec4i>, std::vector<vec3f>, std::vector<vec3f>,
    std::vector<vec2f>>
make_uvseashell(int tesselation, const make_seashell_params& params);

/// Make a bezier circle. Returns bezier, pos.
std::tuple<std::vector<vec4i>, std::vector<vec3f>> make_bezier_circle();

/// Parameters for the make hair function
struct make_hair_params {
    /// minimum and maximum length
    vec2f length = {0.1f, 0.1f};
    /// minimum and maximum radius from base to tip
    vec2f radius = {0.005f, 0.001f};
    /// noise added to hair (strength/scale)
    vec2f noise = zero2f;
    /// clump added to hair (number/strength)
    vec2f clump = zero2f;
    /// rotation
    vec2f rotation = zero2f;
    /// random seed
    uint32_t seed = 0;
};

/// Make a hair ball around a shape. Returns lines, pos, norm, texcoord, radius.
std::tuple<std::vector<vec2i>, std::vector<vec3f>, std::vector<vec3f>,
    std::vector<vec2f>, std::vector<float>>
make_hair(int num, int tesselation, const std::vector<vec3i>& striangles,
    const std::vector<vec4i>& squads, const std::vector<vec3f>& spos,
    const std::vector<vec3f>& snorm, const std::vector<vec2f>& stexcoord,
    const make_hair_params& params);

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMAGE CONTAINERS
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup image Image containers
/// @{

/// Generic image container. Access pixels with at() or operator [].
template <typename T>
struct image {
    /// empty image constructor
    image() : w{0}, h{0}, pixels{} {}
    /// image constructor
    image(int width, int height, const T& val = {})
        : w{width}, h{height}, pixels(size_t(width * height), val) {}

    /// width
    int width() const { return w; }
    /// height
    int height() const { return h; }
    /// check for empty
    bool empty() const { return w == 0 || h == 0; }

    /// Element access
    T& at(int i, int j) { return pixels.at(j * w + i); }
    /// Element access
    const T& at(int i, int j) const { return pixels.at(j * w + i); }

    /// Width and height [private].
    int w, h;
    /// Pixels [private].
    std::vector<T> pixels;
};

/// HDR image
using image4f = image<vec4f>;
/// LDR image
using image4b = image<vec4b>;

/// Pixel iteration.
template <typename T>
inline T* begin(image<T>& a) {
    return a.pixels.data();
}
/// Pixel iteration.
template <typename T>
inline const T* begin(const image<T>& a) {
    return a.pixels.data();
}
/// Pixel iteration.
template <typename T>
inline T* end(image<T>& a) {
    return a.pixels.data() + a.width() * a.height();
}
/// Pixel iteration.
template <typename T>
inline const T* end(const image<T>& a) {
    return a.pixels.data() + a.width() * a.height();
}
/// Pixel access.
template <typename T>
inline T* data(image<T>& a) {
    return a.pixels.data();
}
/// Pixel access.
template <typename T>
inline const T* data(const image<T>& a) {
    return a.pixels.data();
}
/// Number of pixels.
template <typename T>
inline int size(image<T>& a) {
    return a.width() * a.height();
}
/// Chech if an image is empty.
template <typename T>
inline bool empty(image<T>& a) {
    return a.width() * a.height() == 0;
}

/// Create an image with values stored in an array in scanliine order.
template <typename T>
inline image<T> make_image(int width, int height, T* vals) {
    auto img = image<T>(width, height);
    for (auto idx = 0; idx < width * height; idx++) img.pixels[idx] = vals[idx];
    return img;
}

/// Create a 4 channel image with the given number of channels
template <typename T>
inline image<vec<T, 4>> make_image(
    int w, int h, int nc, const T* vals, const vec<T, 4>& def) {
    auto img = image<vec<T, 4>>(w, h);
    for (auto j = 0; j < h; j++) {
        for (auto i = 0; i < w; i++) {
            auto pixel = vals + (j * w + i) * nc;
            img.at(i, j) = def;
            for (auto c = 0; c < nc; c++) img.at(i, j)[c] = pixel[c];
        }
    }
    return img;
}

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMAGE OPERATIONS
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup image_ops Image operations
/// @{

/// Approximate conversion from srgb.
inline vec4f srgb_to_linear(const vec4b& srgb) {
    return {pow(byte_to_float(srgb.x), 2.2f), pow(byte_to_float(srgb.y), 2.2f),
        pow(byte_to_float(srgb.z), 2.2f), byte_to_float(srgb.w)};
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

/// Tone mapping HDR to LDR images.
image4b tonemap_image(
    const image4f& hdr, float exposure, float gamma, bool filmic = false);

/// Image over operator.
void image_over(vec4f* img, int width, int height, int nlayers, vec4f** layers);

/// Image over operator.
void image_over(vec4b* img, int width, int height, int nlayers, vec4b** layers);

/// Converts HSV to RGB.
vec4b hsv_to_rgb(const vec4b& hsv);

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// EXAMPLE IMAGES
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup image_example Example images
/// @{

/// Make a grid image.
image4b make_grid_image(int width, int height, int tile = 64,
    const vec4b& c0 = {90, 90, 90, 255},
    const vec4b& c1 = {128, 128, 128, 255});

/// Make a checkerboard image.
image4b make_checker_image(int width, int height, int tile = 64,
    const vec4b& c0 = {90, 90, 90, 255},
    const vec4b& c1 = {128, 128, 128, 255});

/// Make an image with bumps and dimples.
image4b make_bumpdimple_image(int width, int height, int tile = 64);

/// Make a uv colored grid.
image4b make_ramp_image(
    int width, int height, const vec4b& c0, const vec4b& c1, bool srgb = false);

/// Make a gamma ramp image.
image4b make_gammaramp_image(int width, int height);

/// Make a gamma ramp image.
image4f make_gammaramp_imagef(int width, int height);

/// Make an image color with red/green in the [0,1] range. Helpful to visualize
/// uv texture coordinate application.
image4b make_uv_image(int width, int height);

/// Make a uv colored grid.
image4b make_uvgrid_image(
    int width, int height, int tile = 64, bool colored = true);

/// Make a uv recusive colored grid.
image4b make_recuvgrid_image(
    int width, int height, int tile = 64, bool colored = true);

/// Comvert a bump map to a normal map.
image4b bump_to_normal_map(const image4b& img, float scale = 1);

/// Make a sunsky HDR model with sun at theta elevation in [0,pi/2], turbidity
/// in [1.7,10] with or without sun.
image4f make_sunsky_image(int res, float thetaSun, float turbidity = 3,
    bool has_sun = false, bool has_ground = true);

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

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMAGE LOADING/SAVING
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup image_io Image loading and saving
/// @{

#if YGL_IMAGEIO

/// Check if an image is HDR based on filename.
bool is_hdr_filename(const std::string& filename);

/// Loads a 4 channel ldr image.
image4b load_image4b(const std::string& filename);
/// Loads a 4 channel hdr image.
image4f load_image4f(const std::string& filename);

/// Saves a 4 channel ldr image.
bool save_image4b(const std::string& filename, const image4b& img);
/// Saves a 4 channel hdr image.
bool save_image4f(const std::string& filename, const image4f& img);

/// Loads an image with variable number of channels.
std::vector<float> load_imagef(
    const std::string& filename, int& width, int& height, int& ncomp);
/// Loads an image with variable number of channels.
std::vector<byte> load_image(
    const std::string& filename, int& width, int& height, int& ncomp);

/// Loads an image from memory with variable number of channels.
std::vector<float> load_imagef_from_memory(const std::string& filename,
    const byte* data, int length, int& width, int& height, int& ncomp);
/// Loads an image from memory with variable number of channels.
std::vector<byte> load_image_from_memory(const std::string& filename,
    const byte* data, int length, int& width, int& height, int& ncomp);

/// Saves an image with variable number of channels.
bool save_imagef(const std::string& filename, int width, int height, int ncomp,
    const float* hdr);
/// Saves an image with variable number of channels.
bool save_image(const std::string& filename, int width, int height, int ncomp,
    const byte* ldr);

/// Save a 4 channel HDR or LDR image with tonemapping based on filename.
bool save_image(const std::string& filename, const image4f& hdr, float exposure,
    float gamma, bool filmic = false);

/// Filter type for resizing.
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

/// Edge mode for resizing.
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

/// Resize an image.
void resize_image(const image4f& img, image4f& res_img,
    resize_filter filter = resize_filter::def,
    resize_edge edge = resize_edge::def, bool premultiplied_alpha = false);
/// Resize an image.
void resize_image(const image4b& img, image4b& res_img,
    resize_filter filter = resize_filter::def,
    resize_edge edge = resize_edge::def, bool premultiplied_alpha = false);

#endif

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// RAY-PRIMITIVE INTERSECTION FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup intersect Ray-primitive intersection
/// @{

/// Intersect a ray with a point (approximate).
/// Based on http://geomalgorithms.com/a02-lines.html.
bool intersect_point(const ray3f& ray, const vec3f& p, float r, float& ray_t);

/// Intersect a ray with a line (approximate).
/// Based on http://geomalgorithms.com/a05-intersect-1.html and
/// http://geomalgorithms.com/a07-distance.html#
///     dist3D_Segment_to_Segment
bool intersect_line(const ray3f& ray, const vec3f& v0, const vec3f& v1,
    float r0, float r1, float& ray_t, vec2f& euv);

/// Intersect a ray with a triangle.
bool intersect_triangle(const ray3f& ray, const vec3f& v0, const vec3f& v1,
    const vec3f& v2, float& ray_t, vec2f& euv);

/// Intersect a ray with a quad represented as two triangles (0,1,3) and
/// (2,3,1), with the uv coordinates of the second triangle corrected by u =
/// 1-u' and v = 1-v' to produce a quad parametrization where u and v go from 0
/// to 1. This is equivalent to Intel's Embree.
bool intersect_quad(const ray3f& ray, const vec3f& v0, const vec3f& v1,
    const vec3f& v2, const vec3f& v3, float& ray_t, vec2f& euv);

/// Intersect a ray with a axis-aligned bounding box.
bool intersect_bbox(const ray3f& ray, const bbox3f& bbox);

/// Intersect a ray with a axis-aligned bounding box, implemented as
/// "Robust BVH Ray Traversal" by T. Ize published at
/// http://jcgt.org/published/0002/02/02/paper.pdf
bool intersect_bbox(const ray3f& ray, const vec3f& ray_dinv,
    const vec3i& ray_dsign, const bbox3f& bbox);

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// POINT-PRIMITIVE DISTANCE FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup overlap Point-primitive overlap
/// @{

/// Check if a point overlaps a position within a max distance.
bool overlap_point(
    const vec3f& pos, float dist_max, const vec3f& v0, float r0, float& dist);

/// Find closest line point to a position.
float closestuv_line(const vec3f& pos, const vec3f& v0, const vec3f& v1);

/// Check if a line overlaps a position within a max distance.
bool overlap_line(const vec3f& pos, float dist_max, const vec3f& v0,
    const vec3f& v1, float r0, float r1, float& dist, vec2f& euv);

/// Find closest triangle point to a position.
vec2f closestuv_triangle(
    const vec3f& pos, const vec3f& v0, const vec3f& v1, const vec3f& v2);

/// Check if a triangle overlaps a position within a max distance.
bool overlap_triangle(const vec3f& pos, float dist_max, const vec3f& v0,
    const vec3f& v1, const vec3f& v2, float r0, float r1, float r2, float& dist,
    vec2f& euv);

/// Check if a quad overlaps a position within a max distance.
bool overlap_quad(const vec3f& pos, float dist_max, const vec3f& v0,
    const vec3f& v1, const vec3f& v2, const vec3f& v3, float r0, float r1,
    float r2, float r3, float& dist, vec2f& euv);

/// Check if a bouning box overlaps a position within a max distance.
bool overlap_bbox(const vec3f& pos, float dist_max, const bbox3f& bbox);

/// Check if two bouning boxes overlap.
bool overlap_bbox(const bbox3f& bbox1, const bbox3f& bbox2);

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// BVH FOR RAY INTERSECTION AND CLOSEST ELEMENT
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup bvh Bounding volume hierarchy
/// @{

/// Type of BVH node.
enum struct bvh_node_type : uint32_t {
    /// Internal.
    internal = 0,
    /// Points.
    point = 1,
    /// Lines.
    line = 2,
    /// Triangles.
    triangle = 3,
    /// Quads.
    quad = 4,
    /// Vertices.
    vertex = 8,
    /// Instances.
    instance = 16,
};

/// BVH tree node containing its bounds, indices to the BVH arrays of either
/// sorted primitives or internal nodes, the node element type,
/// and the split axis. Leaf and internal nodes are identical, except that
/// indices refer to primitives for leaf nodes or other nodes for internal
/// nodes. See bvh_tree for more details.
/// This is an internal data structure.
struct bvh_node {
    /// Bounding box.
    bbox3f bbox;
    /// Index to the first sorted primitive/node.
    uint32_t start;
    /// Number of primitives/nodes.
    uint16_t count;
    /// Type of node.
    bvh_node_type type;
    /// Slit axis for internal nodes.
    uint8_t axis;
};

// forward declaration
struct bvh_tree;

/// Shape instance for two-level BVH.
/// This is an internal data structure.
struct bvh_instance {
    /// Frame.
    frame3f frame = identity_frame3f;
    /// Frame inverse.
    frame3f frame_inv = identity_frame3f;
    /// Instance id to be returned.
    int iid = 0;
    /// Shape id to be returned.
    int sid = 0;
    /// Shape bvh.
    bvh_tree* bvh = nullptr;
};

/// BVH tree, stored as a node array. The tree structure is encoded using array
/// indices instead of pointers, both for speed but also to simplify code.
/// BVH nodes indices refer to either the node array, for internal nodes,
/// or a primitive array, for leaf nodes. BVH trees may contain only one type
/// of geometric primitive, like points, lines, triangle or instances of other
/// BVHs. To handle multiple primitive types and transformed primitives, build
/// a two-level hierarchy with the outer BVH, the scene BVH, containing inner
/// BVHs, shape BVHs, each of which of a uniform primitive type.
/// This is an internal data structure.
struct bvh_tree {
    /// Sorted array of internal nodes.
    std::vector<bvh_node> nodes;
    /// Sorted array of elements.
    std::vector<int> sorted_prim;
    /// Leaf element type.
    bvh_node_type type = bvh_node_type::internal;

    /// Positions for shape BVHs.
    std::vector<vec3f> pos;
    /// Radius for shape BVHs.
    std::vector<float> radius;
    /// Points for shape BVHs.
    std::vector<int> points;
    /// Lines for shape BVHs.
    std::vector<vec2i> lines;
    /// Triangles for shape BVHs.
    std::vector<vec3i> triangles;
    /// Quads for shape BVHs.
    std::vector<vec4i> quads;

    /// Instance ids (iid, sid, shape bvh index).
    std::vector<bvh_instance> instances;
    /// Shape BVHs.
    std::vector<bvh_tree*> shape_bvhs;
    /// Whether it owns the memory of the shape BVHs.
    bool own_shape_bvhs = false;

    /// Cleanup.
    ~bvh_tree();
};

/// Build a shape BVH from a set of primitives.
bvh_tree* make_bvh(const std::vector<int>& points,
    const std::vector<vec2i>& lines, const std::vector<vec3i>& triangles,
    const std::vector<vec4i>& quads, const std::vector<vec3f>& pos,
    const std::vector<float>& radius, float def_radius, bool equalsize);
/// Build a scene BVH from a set of shape instances.
bvh_tree* make_bvh(const std::vector<bvh_instance>& instances,
    const std::vector<bvh_tree*>& shape_bvhs, bool own_shape_bvhs,
    bool equal_size);

/// Grab the shape BVHs
inline const std::vector<bvh_tree*>& get_shape_bvhs(const bvh_tree* bvh) {
    return bvh->shape_bvhs;
}

/// Update the node bounds for a shape bvh.
void refit_bvh(bvh_tree* bvh, const std::vector<vec3f>& pos,
    const std::vector<float>& radius, float def_radius);
/// Update the node bounds for a scene bvh
void refit_bvh(bvh_tree* bvh, const std::vector<frame3f>& frames,
    const std::vector<frame3f>& frames_inv);

/// Intersect ray with a bvh returning either the first or any intersection
/// depending on `find_any`. Returns the ray distance `ray_t`, the instance
/// id `iid`, the shape id `sid`, the shape element index `eid` and the
/// shape barycentric coordinates `euv`.
bool intersect_bvh(const bvh_tree* bvh, const ray3f& ray, bool find_any,
    float& ray_t, int& iid, int& sid, int& eid, vec2f& euv);

/// Find a shape element that overlaps a point within a given distance
/// `max_dist`, returning either the closest or any overlap depending on
/// `find_any`. Returns the point distance `dist`, the instance id `iid`, the
/// shape id `sid`, the shape element index `eid` and the shape barycentric
/// coordinates `euv`.
bool overlap_bvh(const bvh_tree* bvh, const vec3f& pos, float max_dist,
    bool find_any, float& dist, int& iid, int& sid, int& eid, vec2f& euv);

/// Intersection point.
struct intersection_point {
    /// Distance of the hit along the ray or from the point.
    float dist = 0;
    /// Instance index.
    int iid = -1;
    /// Shape index.
    int sid = -1;
    /// Shape element index.
    int eid = -1;
    /// Shape barycentric coordinates.
    vec2f euv = zero2f;

    /// Check if intersection is valid.
    operator bool() const { return eid >= 0; }
};

/// Intersect a ray with a bvh (convenience wrapper).
intersection_point intersect_bvh(
    const bvh_tree* bvh, const ray3f& ray, bool early_exit);

/// Finds the closest element with a bvh (convenience wrapper).
intersection_point overlap_bvh(
    const bvh_tree* bvh, const vec3f& pos, float max_dist, bool early_exit);

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// SIMPLE SCENE SUPPORT
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup scene Simple scene
/// @{

/// Texture containing either an LDR or HDR image.
struct texture {
    /// Name.
    std::string name = "";
    /// Path.
    std::string path = "";
    /// If loaded, ldr image.
    image4b ldr = {};
    /// If loaded, hdr image.
    image4f hdr = {};
};

/// Texture information to use for lookup.
struct texture_info {
    /// Texture pointer.
    texture* txt = nullptr;
    /// Wrap s coordinate.
    bool wrap_s = true;
    /// Wrap t coordinate.
    bool wrap_t = true;
    /// Linear interpolation.
    bool linear = true;
    /// Mipmaping.
    bool mipmap = true;
    /// Texture strength (occlusion and normal).
    float scale = 1;

    /// Check whether the texture is present.
    operator bool() const { return (bool)txt; }
};

/// Material type.
enum struct material_type {
    /// Microfacet material type (OBJ).
    specular_roughness = 0,
    /// Base and metallic material (metallic-roughness in glTF).
    metallic_roughness = 1,
    /// Diffuse and specular material (specular-glossness in glTF).
    specular_glossiness = 2,
};

/// Material for surfaces, lines and triangles.
struct material {
    /// Name.
    std::string name = "";
    /// Double-sided rendering.
    bool double_sided = false;
    /// Material type.
    material_type type = material_type::specular_roughness;
    /// Emission color.
    vec3f ke = {0, 0, 0};
    /// Diffuse color / base color.
    vec3f kd = {0, 0, 0};
    /// Specular color / metallic factor.
    vec3f ks = {0, 0, 0};
    /// Clear coat reflection.
    vec3f kr = {0, 0, 0};
    /// Transmission color.
    vec3f kt = {0, 0, 0};
    /// Roughness.
    float rs = 0.0001;
    /// Opacity.
    float op = 1;
    /// Emission texture.
    texture_info ke_txt = {};
    /// Diffuse texture.
    texture_info kd_txt = {};
    /// Specular texture.
    texture_info ks_txt = {};
    /// Clear coat reflection texture.
    texture_info kr_txt = {};
    /// Transmission texture.
    texture_info kt_txt = {};
    /// Roughness texture.
    texture_info rs_txt = {};
    /// Bump map texture (heighfield).
    texture_info bump_txt = {};
    /// Displacement map texture (heighfield).
    texture_info disp_txt = {};
    /// Normal texture.
    texture_info norm_txt = {};
    /// Occlusion texture.
    texture_info occ_txt = {};
};

/// Shape data represented as an indexed array.
/// May contain only one of the points/lines/triangles/quads.
struct shape {
    /// Name.
    std::string name = "";
    /// Material.
    material* mat = nullptr;

    /// Points.
    std::vector<int> points;
    /// Lines.
    std::vector<vec2i> lines;
    /// Triangles.
    std::vector<vec3i> triangles;
    /// Quads.
    std::vector<vec4i> quads;
    /// Face-varying indices for position.
    std::vector<vec4i> quads_pos;
    /// Face-varying indices for normal.
    std::vector<vec4i> quads_norm;
    /// Face-varying indices for texcoord.
    std::vector<vec4i> quads_texcoord;
    /// Bezier.
    std::vector<vec4i> beziers;

    /// Vertex position.
    std::vector<vec3f> pos;
    /// Vertex normals.
    std::vector<vec3f> norm;
    /// Vertex texcoord.
    std::vector<vec2f> texcoord;
    /// Vertex second texcoord.
    std::vector<vec2f> texcoord1;
    /// Vertex color.
    std::vector<vec4f> color;
    /// per-vertex radius.
    std::vector<float> radius;
    /// Vertex tangent space.
    std::vector<vec4f> tangsp;

    /// Number of times to subdivide.
    int subdivision = 0;
    /// Whether to use Catmull-Clark subdivision.
    bool catmullclark = false;
};

/// Group of shapes.
struct shape_group {
    /// Name.
    std::string name = "";
    /// Path used for saving in glTF.
    std::string path = "";
    /// Shapes.
    std::vector<shape*> shapes;

    /// Cleanup.
    ~shape_group();
};

/// Shape instance.
struct instance {
    // Name.
    std::string name = "";
    /// Transform frame.
    frame3f frame = identity_frame3f;
    /// Shape instance.
    shape_group* shp = nullptr;
};

/// Camera.
struct camera {
    /// Name.
    std::string name = "";
    /// Transform frame.
    frame3f frame = identity_frame3f;
    /// Orthographic camera.
    bool ortho = false;
    /// Vertical field of view.
    float yfov = 2;
    /// Aspect ratio.
    float aspect = 16.0f / 9.0f;
    /// Focus distance.
    float focus = 1;
    /// Lens aperture.
    float aperture = 0;
    /// Near plane distance.
    float near = 0.01f;
    /// Far plane distance.
    float far = 10000;
};

/// Envinonment map.
struct environment {
    /// Name.
    std::string name = "";
    /// Transform frame.
    frame3f frame = identity_frame3f;
    /// Emission coefficient.
    vec3f ke = {0, 0, 0};
    /// Emission texture.
    texture_info ke_txt = {};
};

/// Node in a transform hierarchy.
struct node {
    /// Name.
    std::string name = "";
    /// Parent node.
    node* parent = nullptr;
    /// Transform frame.
    frame3f frame = identity_frame3f;
    /// Translation.
    vec3f translation = zero3f;
    /// Rotation.
    quat4f rotation = {0, 0, 0, 1};
    /// Scaling.
    vec3f scaling = {1, 1, 1};
    /// Weights for morphing.
    std::vector<float> weights = {};
    /// Camera the node points to.
    camera* cam = nullptr;
    /// Instance the node points to.
    instance* ist = nullptr;
    /// Environment the node points to.
    environment* env = nullptr;

    /// Child nodes. This is a computed value only stored for convenience.
    std::vector<node*> children_ = {};
};

/// Keyframe type.
enum struct keyframe_type {
    /// Linear interpolation.
    linear = 0,
    /// Step function.
    step = 1,
    /// Catmull-Rom interpolation.
    catmull_rom = 2,
    /// Cubic Bezier interpolation.
    bezier = 3,
};

/// Keyframe data.
struct animation {
    /// Name.
    std::string name;
    /// Interpolation.
    keyframe_type type = keyframe_type::linear;
    /// Times.
    std::vector<float> times;
    /// Translation.
    std::vector<vec3f> translation;
    /// Rotation.
    std::vector<quat4f> rotation;
    /// Scaling.
    std::vector<vec3f> scaling;
    /// Weights for morphing.
    std::vector<std::vector<float>> weights;
};

/// Animation made of multiple keyframed values.
struct animation_group {
    /// Name.
    std::string name;
    /// Path  used when writing files on disk with glTF.
    std::string path = "";
    /// Keyframed values.
    std::vector<animation*> animations;
    /// Binds keyframe values to nodes.
    std::vector<std::pair<animation*, node*>> targets;

    // Cleanup
    ~animation_group();
};

/// Scene comprised an array of objects whose memory is owened by the scene.
/// All members are optional, but different algorithm might require different
/// data to be laoded. Scene objects (camera, instances, environments) have
/// transforms defined internally. A scene can optionally contain a
/// node hierarchy where each node might point to a camera, instance or
/// environment. In that case, the element transforms are computed from
/// the hierarchy. Animation is also optional, with keyframe data that
/// updates node transformations only if defined.
struct scene {
    /// Shapes.
    std::vector<shape_group*> shapes = {};
    /// Shape instances.
    std::vector<instance*> instances = {};
    /// Materials.
    std::vector<material*> materials = {};
    /// Textures.
    std::vector<texture*> textures = {};
    /// Cameras.
    std::vector<camera*> cameras = {};
    /// Environments.
    std::vector<environment*> environments = {};

    /// Node hierarchy.
    std::vector<node*> nodes = {};
    /// Node animations.
    std::vector<animation_group*> animations = {};

    // Cleanup.
    ~scene();
};

/// Shape position interpolated using barycentric coordinates.
vec3f eval_pos(const shape* shp, int eid, const vec2f& euv);
/// Shape normal interpolated using barycentric coordinates.
vec3f eval_norm(const shape* shp, int eid, const vec2f& euv);
/// Shape texcoord interpolated using barycentric coordinates.
vec2f eval_texcoord(const shape* shp, int eid, const vec2f& euv);
/// Shape color interpolated using barycentric coordinates.
vec4f eval_color(const shape* shp, int eid, const vec2f& euv);
/// Shape radius interpolated using barycentric coordinates.
float eval_radius(const shape* shp, int eid, const vec2f& euv);
/// Shape tangent space interpolated using barycentric coordinates.
vec4f eval_tangsp(const shape* shp, int eid, const vec2f& euv);
/// Instance position interpolated using barycentric coordinates.
vec3f eval_pos(const instance* ist, int sid, int eid, const vec2f& euv);
/// Instance normal interpolated using barycentric coordinates.
vec3f eval_norm(const instance* ist, int sid, int eid, const vec2f& euv);

/// Evaluate a texture.
vec4f eval_texture(const texture_info& info, const vec2f& texcoord,
    bool srgb = true, const vec4f& def = {1, 1, 1, 1});
/// Generates a ray from a camera for image plane coordinate uv and the
/// lens coordinates luv.
ray3f eval_camera_ray(const camera* cam, const vec2f& uv, const vec2f& luv);

/// Finds an element by name.
template <typename T>
inline T* find_named_elem(
    const std::vector<T*>& elems, const std::string& name) {
    if (name == "") return nullptr;
    for (auto elem : elems)
        if (elem->name == name) return elem;
    return nullptr;
}
// Adds a named element or return the one with that name.
template <typename T>
inline T* add_named_elem(std::vector<T*>& elems, const std::string& name) {
    for (auto elem : elems)
        if (elem->name == name) return elem;
    auto elem = new T();
    elem->name = name;
    elems.push_back(elem);
    return elem;
}

/// Subdivides shape elements. Apply subdivision surface rules if subdivide
/// is true.
void subdivide_shape_once(shape* shp, bool subdiv = false);
/// Facet a shape. Supports only non-facevarying shapes.
void facet_shape(shape* shp, bool recompute_normals = true);
/// Tesselate a shape into basic primitives.
void tesselate_shape(shape* shp, bool subdivide,
    bool facevarying_to_sharedvertex, bool quads_to_triangles,
    bool bezier_to_lines);
/// Tesselate scene shapes.
void tesselate_shapes(scene* scn, bool subdivide,
    bool facevarying_to_sharedvertex, bool quads_to_triangles,
    bool bezier_to_lines);

/// Update node transforms.
void update_transforms(scene* scn, float time = 0);
/// Compute animation range.
vec2f compute_animation_range(const scene* scn);

/// Make a view camera either copying a given one or building a default one.
camera* make_view_camera(const scene* scn, int camera_id);

/// Computes a shape bounding box using a quick computation that ignores radius.
bbox3f compute_bounds(const shape* shp);
/// Compute a scene bounding box.
bbox3f compute_bounds(const scene* scn);

/// Flatten scene instances into separate shapes.
void flatten_instances(scene* scn);

/// Print scene information.
void print_info(const scene* scn);

/// Build a shape BVH.
bvh_tree* make_bvh(
    const shape* shp, float def_radius = 0.001f, bool equalsize = true);
/// Build a scene BVH.
bvh_tree* make_bvh(
    const scene* scn, float def_radius = 0.001f, bool equalsize = true);

/// Refits a scene BVH.
void refit_bvh(bvh_tree* bvh, const shape* shp, float def_radius = 0.001f);
/// Refits a scene BVH.
void refit_bvh(
    bvh_tree* bvh, const scene* scn, bool do_shapes, float def_radius = 0.001f);

/// Add elements options.
struct add_elements_options {
    /// Add missing normal.
    bool smooth_normals = true;
    /// Add missing trangent space.
    bool tangent_space = true;
    /// Add empty texture data.
    bool texture_data = true;
    /// Add instances.
    bool shape_instances = true;
    /// Add default names.
    bool default_names = true;
    /// Add default paths.
    bool default_paths = true;

    /// Initialize to no elements.
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

/// Loading options.
struct load_options {
    /// Whether to load textures.
    bool load_textures = true;
    /// Skip missing files without giving and error.
    bool skip_missing = true;
    /// Whether to flip the v coordinate in OBJ.
    bool obj_flip_texcoord = true;
    /// Duplicate vertices if smoothing off in OBJ.
    bool obj_facet_non_smooth = false;
    /// Whether to flip tr in OBJ.
    bool obj_flip_tr = true;
    /// Whether to preserve quads.
    bool preserve_quads = false;
    /// Whether to preserve face-varying faces.
    bool preserve_facevarying = false;
    /// Whether to preserve node hierarchy.
    bool preserve_hierarchy = false;
};

/// Loads a scene. For now OBJ or glTF are supported.
/// Throws an exception if an error occurs.
scene* load_scene(const std::string& filename, const load_options& opts = {});

/// Save options.
struct save_options {
    /// Whether to save textures.
    bool save_textures = true;
    /// Skip missing files without giving and error.
    bool skip_missing = true;
    /// Whether to flip the v coordinate in OBJ.
    bool obj_flip_texcoord = true;
    /// Whether to flip tr in OBJ.
    bool obj_flip_tr = true;
    /// Whether to use separate buffers in gltf.
    bool gltf_separate_buffers = false;
};

/// Saves a scene. For now OBJ and glTF are supported.
void save_scene(
    const std::string& filename, const scene* scn, const save_options& opts);

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// SCENE TYPE SUPPORT
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup scene_type Scene type support
/// @{

/// Names of enum values.
template <>
inline const std::vector<std::pair<std::string, material_type>>&
enum_names<material_type>() {
    static auto names = std::vector<std::pair<std::string, material_type>>{
        {"specular_roughness", material_type::specular_roughness},
        {"metallic_roughness", material_type::metallic_roughness},
        {"specular_glossiness", material_type::specular_glossiness},
    };
    return names;
}

/// Names of enum values.
template <>
inline const std::vector<std::pair<std::string, keyframe_type>>&
enum_names<keyframe_type>() {
    static auto names = std::vector<std::pair<std::string, keyframe_type>>{
        {"linear", keyframe_type::linear},
        {"step", keyframe_type::step},
        {"bezier", keyframe_type::bezier},
        {"catmull_rom", keyframe_type::catmull_rom},
    };
    return names;
}

/// Visit struct elements.
template <typename Visitor>
inline void visit(texture& val, Visitor&& visitor) {
    visitor("name", val.name, visit_sem{visit_sem_type::name});
    visitor("path", val.path, visit_sem{visit_sem_type::path});
    visitor("ldr", val.ldr, visit_sem{visit_sem_type::value});
    visitor("hdr", val.hdr, visit_sem{visit_sem_type::value});
}

/// Visit struct elements.
template <typename Visitor>
inline void visit(texture_info& val, Visitor&& visitor) {
    visitor("txt", val.txt, visit_sem{visit_sem_type::reference});
    visitor("wrap_s", val.wrap_s, visit_sem{visit_sem_type::value});
    visitor("wrap_t", val.wrap_t, visit_sem{visit_sem_type::value});
    visitor("linear", val.linear, visit_sem{visit_sem_type::value});
    visitor("mipmap", val.mipmap, visit_sem{visit_sem_type::value});
    visitor("scale", val.scale, visit_sem{visit_sem_type::value, 0, 10});
}

/// Visit struct elements.
template <typename Visitor>
inline void visit(material& val, Visitor&& visitor) {
    visitor("name", val.name, visit_sem{visit_sem_type::name});
    visitor("double_sided", val.double_sided, visit_sem{visit_sem_type::value});
    visitor("type", val.type, visit_sem{visit_sem_type::value});
    visitor("ke", val.ke, visit_sem{visit_sem_type::color, 0, 1000});
    visitor("kd", val.kd, visit_sem{visit_sem_type::color});
    visitor("ks", val.ks, visit_sem{visit_sem_type::color});
    visitor("kr", val.kr, visit_sem{visit_sem_type::color});
    visitor("kt", val.kt, visit_sem{visit_sem_type::color});
    visitor("rs", val.rs, visit_sem{visit_sem_type::value});
    visitor("op", val.op, visit_sem{visit_sem_type::value});
    visitor("ke_txt", val.ke_txt, visit_sem{visit_sem_type::value});
    visitor("kd_txt", val.kd_txt, visit_sem{visit_sem_type::value});
    visitor("ks_txt", val.ks_txt, visit_sem{visit_sem_type::value});
    visitor("kr_txt", val.kr_txt, visit_sem{visit_sem_type::value});
    visitor("kt_txt", val.kt_txt, visit_sem{visit_sem_type::value});
    visitor("rs_txt", val.rs_txt, visit_sem{visit_sem_type::value});
    visitor("bump_txt", val.bump_txt, visit_sem{visit_sem_type::value});
    visitor("disp_txt", val.disp_txt, visit_sem{visit_sem_type::value});
    visitor("norm_txt", val.norm_txt, visit_sem{visit_sem_type::value});
    visitor("occ_txt", val.occ_txt, visit_sem{visit_sem_type::value});
}

/// Visit struct elements.
template <typename Visitor>
inline void visit(shape& val, Visitor&& visitor) {
    visitor("name", val.name, visit_sem{visit_sem_type::name});
    visitor("mat", val.mat, visit_sem{visit_sem_type::reference});
    visitor("points", val.points, visit_sem{visit_sem_type::value});
    visitor("lines", val.lines, visit_sem{visit_sem_type::value});
    visitor("triangles", val.triangles, visit_sem{visit_sem_type::value});
    visitor("quads", val.quads, visit_sem{visit_sem_type::value});
    visitor("quads_pos", val.quads_pos, visit_sem{visit_sem_type::value});
    visitor("quads_norm", val.quads_norm, visit_sem{visit_sem_type::value});
    visitor(
        "quads_texcoord", val.quads_texcoord, visit_sem{visit_sem_type::value});
    visitor("beziers", val.beziers, visit_sem{visit_sem_type::value});
    visitor("pos", val.pos, visit_sem{visit_sem_type::value});
    visitor("norm", val.norm, visit_sem{visit_sem_type::value});
    visitor("texcoord", val.texcoord, visit_sem{visit_sem_type::value});
    visitor("texcoord1", val.texcoord1, visit_sem{visit_sem_type::value});
    visitor("color", val.color, visit_sem{visit_sem_type::value});
    visitor("radius", val.radius, visit_sem{visit_sem_type::value});
    visitor("tangsp", val.tangsp, visit_sem{visit_sem_type::value});
    visitor("subdivision", val.subdivision, visit_sem{visit_sem_type::value});
    visitor("catmullclark", val.catmullclark, visit_sem{visit_sem_type::value});
}

/// Visit struct elements.
template <typename Visitor>
inline void visit(shape_group& val, Visitor&& visitor) {
    visitor("name", val.name, visit_sem{visit_sem_type::name});
    visitor("path", val.path, visit_sem{visit_sem_type::value});
    visitor("shapes", val.shapes, visit_sem{visit_sem_type::object});
}

/// Visit struct elements.
template <typename Visitor>
inline void visit(instance& val, Visitor&& visitor) {
    visitor("name", val.name, visit_sem{visit_sem_type::name});
    visitor("frame", val.frame, visit_sem{visit_sem_type::value});
    visitor("shp", val.shp, visit_sem{visit_sem_type::reference});
}

/// Visit struct elements.
template <typename Visitor>
inline void visit(camera& val, Visitor&& visitor) {
    visitor("name", val.name, visit_sem{visit_sem_type::name});
    visitor("frame", val.frame, visit_sem{visit_sem_type::value});
    visitor("ortho", val.ortho, visit_sem{visit_sem_type::value});
    visitor("yfov", val.yfov, visit_sem{visit_sem_type::value, 0.01f, 10});
    visitor("aspect", val.aspect, visit_sem{visit_sem_type::value, 0.3f, 3.0f});
    visitor(
        "focus", val.focus, visit_sem{visit_sem_type::value, 0.0f, 10000.0f});
    visitor("aperture", val.aperture, visit_sem{visit_sem_type::value, 0, 10});
    visitor("near", val.near, visit_sem{visit_sem_type::value, 0.001f, 10.0f});
    visitor("far", val.far, visit_sem{visit_sem_type::value, 1.0f, 10000.0f});
}

/// Visit struct elements.
template <typename Visitor>
inline void visit(environment& val, Visitor&& visitor) {
    visitor("name", val.name, visit_sem{visit_sem_type::name});
    visitor("frame", val.frame, visit_sem{visit_sem_type::value});
    visitor("ke", val.ke, visit_sem{visit_sem_type::color, 0, 1000});
    visitor("ke_txt", val.ke_txt, visit_sem{visit_sem_type::value});
}

/// Visit struct elements.
template <typename Visitor>
inline void visit(node& val, Visitor&& visitor) {
    visitor("name", val.name, visit_sem{visit_sem_type::name});
    visitor("parent", val.parent, visit_sem{visit_sem_type::reference});
    visitor("frame", val.frame, visit_sem{visit_sem_type::value, -10, 10});
    visitor("translation", val.translation,
        visit_sem{visit_sem_type::value, -10, 10});
    visitor("rotation", val.rotation, visit_sem{visit_sem_type::value});
    visitor(
        "scaling", val.scaling, visit_sem{visit_sem_type::value, .01f, 10.0f});
    visitor("weights", val.weights, visit_sem{visit_sem_type::value});
    visitor("cam", val.cam, visit_sem{visit_sem_type::reference});
    visitor("ist", val.ist, visit_sem{visit_sem_type::reference});
    visitor("env", val.env, visit_sem{visit_sem_type::reference});
}

/// Visit struct elements.
template <typename Visitor>
inline void visit(animation& val, Visitor&& visitor) {
    visitor("name", val.name, visit_sem{visit_sem_type::name});
    visitor("type", val.type, visit_sem{visit_sem_type::value});
    visitor("times", val.times, visit_sem{visit_sem_type::value});
    visitor("translation", val.translation,
        visit_sem{visit_sem_type::value, -10, 10});
    visitor("rotation", val.rotation, visit_sem{visit_sem_type::value});
    visitor(
        "scaling", val.scaling, visit_sem{visit_sem_type::value, 0.01f, 10.0f});
    visitor("weights", val.weights, visit_sem{visit_sem_type::value});
}

/// Visit struct elements.
template <typename Visitor>
inline void visit(animation_group& val, Visitor&& visitor) {
    visitor("name", val.name, visit_sem{visit_sem_type::name});
    visitor("path", val.path, visit_sem{visit_sem_type::value});
    visitor("animations", val.animations, visit_sem{visit_sem_type::object});
    visitor("targets", val.targets, visit_sem{visit_sem_type::reference});
}

/// Visit struct elements.
template <typename Visitor>
inline void visit(scene& val, Visitor&& visitor) {
    visitor("cameras", val.cameras, visit_sem{visit_sem_type::object});
    visitor("shapes", val.shapes, visit_sem{visit_sem_type::object});
    visitor("instances", val.instances, visit_sem{visit_sem_type::object});
    visitor(
        "environments", val.environments, visit_sem{visit_sem_type::object});
    visitor("materials", val.materials, visit_sem{visit_sem_type::object});
    visitor("textures", val.textures, visit_sem{visit_sem_type::object});
    visitor("nodes", val.nodes, visit_sem{visit_sem_type::object});
    visitor("animations", val.animations, visit_sem{visit_sem_type::object});
}

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// EXAMPLE SCENES
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup scene_example Example scenes
/// @{

/// Makes the Cornell Box scene.
scene* make_cornell_box_scene();

/// Test camera parameters.
struct test_camera_params {
    /// Name.
    std::string name = "";
    /// From point.
    vec3f from = {0, 0, -1};
    /// To point.
    vec3f to = zero3f;
    /// Field of view.
    float yfov = 45 * pif / 180;
    /// Aspect ratio.
    float aspect = 1;
};

/// Updates a test camera.
void update_test_elem(
    const scene* scn, camera* cam, const test_camera_params& tcam);

/// Test camera presets.
std::unordered_map<std::string, test_camera_params>& test_camera_presets();

/// Test texture type.
enum struct test_texture_type {
    /// None (empty texture).
    none,
    /// Grid image.
    grid,
    /// Checker image.
    checker,
    /// Colored checker image.
    colored,
    /// Detailed colored checker image.
    rcolored,
    /// Bump and dimple imahe.
    bump,
    /// Uv debug image.
    uv,
    /// Gamma ramp.
    gamma,
    /// Perlin noise image.
    noise,
    /// Perlin ridge image.
    ridge,
    /// Perlin fbm image.
    fbm,
    /// Perlin turbulence image.
    turbulence,
    /// Gamma ramp (HDR).
    gammaf,
    /// Sky (HDR).
    sky,
};

/// Test texture parameters.
struct test_texture_params {
    /// Name.
    std::string name = "";
    /// Type.
    test_texture_type type = test_texture_type::none;
    /// Resolution.
    int resolution = 512;
    /// Tile size for grid-like textures.
    int tile_size = 64;
    /// Noise scale for noise-like textures.
    int noise_scale = 8;
    /// Sun angle for sunsky-like textures.
    float sky_sunangle = pif / 4;
    /// Convert to normal map.
    bool bump_to_normal = false;
    /// Bump to normal scale.
    float bump_scale = 4;
};

/// Updates a test texture.
void update_test_elem(
    const scene* scn, texture* txt, const test_texture_params& ttxt);

/// Test texture presets.
std::unordered_map<std::string, test_texture_params>& test_texture_presets();

/// Test material type.
enum struct test_material_type {
    /// None (empty material).
    none,
    /// Emission.
    emission,
    /// Matte (diffuse).
    matte,
    /// Plastic.
    plastic,
    /// Matetal.
    metal,
    /// Transparent (diffuse with opacity).
    transparent,
};

/// Test material parameters.
struct test_material_params {
    /// Name.
    std::string name = "";
    /// Type.
    test_material_type type = test_material_type::matte;
    /// Emission strenght.
    float emission = 1;
    /// Base color.
    vec3f color = {0.2, 0.2, 0.2};
    /// Opacity (only for supported materials).
    float opacity = 1;
    /// Roughness.
    float roughness = 0.1;
    /// Base texture.
    std::string texture = "";
    /// Normal map.
    std::string normal = "";
};

/// Updates a test material.
void update_test_elem(
    const scene* scn, material* mat, const test_material_params& tmat);

/// Test material presets.
std::unordered_map<std::string, test_material_params>& test_material_presets();

/// Test shape type.
enum struct test_shape_type {
    /// Floor (shared vertex, 20x20 size).
    floor,
    /// Quad (shared vertex).
    quad,
    /// Cube (shared vertex, not watertight).
    cube,
    /// Sphere (shared vertex, not watertight).
    sphere,
    /// Sphere with cube uvs (shared vertex, not watertight).
    spherecube,
    /// Spherized cube (shared vertex, not watertight).
    spherizedcube,
    /// Geodesic sphere (shared vertex, watertight, no texcoord).
    geosphere,
    /// Sphere with flipped cap (shared vertex, not watertight).
    flipcapsphere,
    /// Suzanne (shared vertex, no texcoord).
    suzanne,
    /// Position-only cube (shared vertex).
    cubep,
    /// Face-varying cube (shared vertex).
    fvcube,
    /// Face-varying sphere (shared vertex).
    fvsphere,
    /// Matball (shared vertex, not watertight).
    matball,
    /// Single point.
    point,
    /// Random points in a cube.
    pointscube,
    /// Random lines on a sphere.
    hairball,
    /// Bezier circle.
    beziercircle,
};

/// Test shape parameters.
struct test_shape_params {
    /// Shape name (if not filled, assign a default based on type).
    std::string name = "";
    /// Shape type.
    test_shape_type type = test_shape_type::sphere;
    /// Material name.
    std::string material = "";
    /// Interior material name.
    std::string interior = "";
    /// Level of shape tesselatation (-1 for default).
    int tesselation = -1;
    /// Level of shape tesselation for subdivision surfaces.
    int subdivision = 0;
    /// Shape scale.
    float scale = 1;
    /// Radius for points and lines.
    float radius = -1;
    /// Faceted shape.
    bool faceted = false;
    /// Number of elements for points and lines (-1 for default).
    int num = -1;
    /// Hair generation params.
    make_hair_params hair_params = {};
};

/// Updates a test shape, adding it to the scene if missing.
void update_test_elem(
    const scene* scn, shape* shp, const test_shape_params& tshp);

/// Test shape presets.
std::unordered_map<std::string, test_shape_params>& test_shape_presets();

/// Test instance parameters.
struct test_instance_params {
    /// Name (if not filled, assign a default one).
    std::string name = "";
    /// Shape name.
    std::string shape = "";
    /// Base frame.
    frame3f frame = identity_frame3f;
    /// Rotation in Euler angles.
    vec3f rotation = zero3f;
};

/// Updates a test instance.
void update_test_elem(
    const scene* scn, instance* ist, const test_instance_params& tist);

/// Test instance presets.
std::unordered_map<std::string, test_instance_params>& test_instance_presets();

/// Test environment parameters.
struct test_environment_params {
    /// Name.
    std::string name = "";
    /// Emission strenght.
    float emission = 1;
    /// Emission color.
    vec3f color = {1, 1, 1};
    /// Emission texture.
    std::string texture = "";
    /// Frame.
    frame3f frame = identity_frame3f;
    /// Rotation around y axis.
    float rotation = 0;
};

/// Updates a test instance.
void update_test_elem(
    const scene* scn, environment* env, const test_environment_params& tenv);

/// Test environment presets.
std::unordered_map<std::string, test_environment_params>&
test_environment_presets();

/// Test node parameters.
struct test_node_params {
    /// Name.
    std::string name = "";
    /// Parent node.
    std::string parent = "";
    /// Camera.
    std::string camera = "";
    /// Instance.
    std::string instance = "";
    /// Environment.
    std::string environment = "";
    /// Frame.
    frame3f frame = identity_frame3f;
    /// Translation.
    vec3f translation = {0, 0, 0};
    /// Roation.
    quat4f rotation = {0, 0, 0, 1};
    /// Scaling.
    vec3f scaling = {1, 1, 1};
};

/// Updates a test node.
void update_test_elem(
    const scene* scn, node* nde, const test_node_params& tndr);

/// Test nodes presets.
std::unordered_map<std::string, test_node_params>& test_node_presets();

/// Test animation parameters.
struct test_animation_params {
    /// Name.
    std::string name = "";
    /// Linear or bezier.
    bool bezier = false;
    /// Animation speed.
    float speed = 1;
    /// Animation scale.
    float scale = 1;
    /// Keyframes times.
    std::vector<float> times = {};
    /// Translation keyframes.
    std::vector<vec3f> translation = {};
    /// Rotation keyframes.
    std::vector<quat4f> rotation = {};
    /// Scale keyframes.
    std::vector<vec3f> scaling = {};
    /// Environment.
    std::vector<std::string> nodes = {};
};

/// Updates a test node.
void update_test_elem(
    const scene* scn, animation_group* anm, const test_animation_params& tndr);

/// Test nodes presets.
std::unordered_map<std::string, test_node_params>& test_node_presets();

/// Test scene.
struct test_scene_params {
    /// Name.
    std::string name;
    /// Cameras.
    std::vector<test_camera_params> cameras;
    /// Textures.
    std::vector<test_texture_params> textures;
    /// Materials.
    std::vector<test_material_params> materials;
    /// Shapes.
    std::vector<test_shape_params> shapes;
    /// Instances.
    std::vector<test_instance_params> instances;
    /// Environmennts.
    std::vector<test_environment_params> environments;
    /// Nodes.
    std::vector<test_node_params> nodes;
    /// Animations.
    std::vector<test_animation_params> animations;
};

/// Updates a test scene, adding missing objects. Objects are only added and
/// never removed.
void update_test_scene(scene* scn, const test_scene_params& tscn,
    const std::unordered_set<void*>& refresh = {});

/// Makes a test scene. Convenience wrapper around `update_test_scene()`.
inline scene* make_test_scene(const test_scene_params& tscn) {
    auto scn = new scene();
    update_test_scene(scn, tscn);
    return scn;
}

/// Test scene presets.
std::unordered_map<std::string, test_scene_params>& test_scene_presets();

/// Remove duplicates based on name.
void remove_duplicates(test_scene_params& tscn);

/// Load test scene.
test_scene_params load_test_scene(const std::string& filename);

/// Save test scene.
void save_test_scene(const std::string& filename, const test_scene_params& scn);

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// EXAMPLE SCENE TYPE SUPPORT
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup scene_example_type Example scenes type support
/// @{

/// Names of enum values.
template <>
inline const std::vector<std::pair<std::string, test_texture_type>>&
enum_names<test_texture_type>() {
    static auto names = std::vector<std::pair<std::string, test_texture_type>>{
        {"none", test_texture_type::none},
        {"grid", test_texture_type::grid},
        {"colored", test_texture_type::colored},
        {"checker", test_texture_type::checker},
        {"rcolored", test_texture_type::rcolored},
        {"bump", test_texture_type::bump},
        {"uv", test_texture_type::uv},
        {"gamma", test_texture_type::gamma},
        {"noise", test_texture_type::noise},
        {"ridge", test_texture_type::ridge},
        {"fbm", test_texture_type::fbm},
        {"turbulence", test_texture_type::turbulence},
        {"gammaf", test_texture_type::gammaf},
        {"sky", test_texture_type::sky},
    };
    return names;
}

/// Names of enum values.
template <>
inline const std::vector<std::pair<std::string, test_material_type>>&
enum_names<test_material_type>() {
    static auto names = std::vector<std::pair<std::string, test_material_type>>{
        {"none", test_material_type::none},
        {"emission", test_material_type::emission},
        {"matte", test_material_type::matte},
        {"plastic", test_material_type::plastic},
        {"metal", test_material_type::metal},
        {"transparent", test_material_type::transparent},
    };
    return names;
}

/// Names of enum values.
template <>
inline const std::vector<std::pair<std::string, test_shape_type>>&
enum_names<test_shape_type>() {
    static auto names = std::vector<std::pair<std::string, test_shape_type>>{
        {"floor", test_shape_type::floor},
        {"quad", test_shape_type::quad},
        {"cube", test_shape_type::cube},
        {"sphere", test_shape_type::sphere},
        {"spherecube", test_shape_type::spherecube},
        {"spherizedcube", test_shape_type::spherizedcube},
        {"geosphere", test_shape_type::geosphere},
        {"flipcapsphere", test_shape_type::flipcapsphere},
        {"suzanne", test_shape_type::suzanne},
        {"cubep", test_shape_type::cubep},
        {"fvcube", test_shape_type::fvcube},
        {"fvsphere", test_shape_type::fvsphere},
        {"matball", test_shape_type::matball},
        {"point", test_shape_type::point},
        {"pointscube", test_shape_type::pointscube},
        {"hairball", test_shape_type::hairball},
        {"beziercircle", test_shape_type::beziercircle},
    };
    return names;
}

/// Visit struct elements.
template <typename Visitor>
inline void visit(test_camera_params& val, Visitor&& visitor) {
    visitor("name", val.name, visit_sem{visit_sem_type::name});
    visitor("from", val.from, visit_sem{visit_sem_type::value, -10, 10});
    visitor("to", val.to, visit_sem{visit_sem_type::value, -10, 10});
    visitor("yfov", val.yfov, visit_sem{visit_sem_type::value});
    visitor("aspect", val.aspect, visit_sem{visit_sem_type::value});
}

/// Visit struct elements.
template <typename Visitor>
inline void visit(test_texture_params& val, Visitor&& visitor) {
    visitor("name", val.name, visit_sem{visit_sem_type::name});
    visitor("resolution", val.resolution,
        visit_sem{visit_sem_type::value, 256, 1024});
    visitor(
        "tile_size", val.tile_size, visit_sem{visit_sem_type::value, -1, 256});
    visitor(
        "noise_scale", val.noise_scale, visit_sem{visit_sem_type::value, 0, 8});
    visitor("sky_sunangle", val.sky_sunangle,
        visit_sem{visit_sem_type::value, 0, pif});
    visitor(
        "bump_to_normal", val.bump_to_normal, visit_sem{visit_sem_type::value});
    visitor(
        "bump_scale", val.bump_scale, visit_sem{visit_sem_type::value, 1, 10});
}

/// Visit struct elements.
template <typename Visitor>
inline void visit(test_material_params& val, Visitor&& visitor) {
    visitor("name", val.name, visit_sem{visit_sem_type::name});
    visitor("type", val.type, visit_sem{visit_sem_type::value});
    visitor("emission", val.emission, visit_sem{visit_sem_type::value, 1000});
    visitor("color", val.color, visit_sem{visit_sem_type::value});
    visitor("opacity", val.opacity, visit_sem{visit_sem_type::value});
    visitor("roughness", val.roughness, visit_sem{visit_sem_type::value});
    visitor("texture", val.texture, visit_sem{visit_sem_type::value});
    visitor("normal", val.normal, visit_sem{visit_sem_type::value});
}

/// Visit struct elements.
template <typename Visitor>
inline void visit(test_shape_params& val, Visitor&& visitor) {
    visitor("name", val.name, visit_sem{visit_sem_type::name});
    visitor("type", val.type, visit_sem{visit_sem_type::value});
    visitor("material", val.material, visit_sem{visit_sem_type::value});
    visitor("interior", val.interior, visit_sem{visit_sem_type::value});
    visitor("tesselation", val.tesselation,
        visit_sem{visit_sem_type::value, -1, 10});
    visitor("subdivision", val.subdivision,
        visit_sem{visit_sem_type::value, -1, 10});
    visitor("scale", val.scale, visit_sem{visit_sem_type::value, 0.01f, 10.0f});
    visitor(
        "radius", val.radius, visit_sem{visit_sem_type::value, 0.0001f, 0.01f});
    visitor("faceted", val.faceted, visit_sem{visit_sem_type::value});
    visitor("num", val.num, visit_sem{visit_sem_type::value, -1, 10000});
}

/// Visit struct elements.
template <typename Visitor>
inline void visit(test_instance_params& val, Visitor&& visitor) {
    visitor("name", val.name, visit_sem{visit_sem_type::name});
    visitor("shape", val.shape, visit_sem{visit_sem_type::value});
    visitor("frame", val.frame, visit_sem{visit_sem_type::value});
    visitor("rotation", val.rotation, visit_sem{visit_sem_type::value});
}

/// Visit struct elements.
template <typename Visitor>
inline void visit(test_environment_params& val, Visitor&& visitor) {
    visitor("name", val.name, visit_sem{visit_sem_type::name});
    visitor(
        "emission", val.emission, visit_sem{visit_sem_type::value, 0, 1000});
    visitor("color", val.color, visit_sem{visit_sem_type::value});
    visitor("texture", val.texture, visit_sem{visit_sem_type::value});
    visitor("frame", val.frame, visit_sem{visit_sem_type::value});
    visitor("rotation", val.rotation, visit_sem{visit_sem_type::value});
}

/// Visit struct elements.
template <typename Visitor>
inline void visit(test_node_params& val, Visitor&& visitor) {
    visitor("name", val.name, visit_sem{visit_sem_type::name});
    visitor("parent", val.parent, visit_sem{visit_sem_type::value});
    visitor("camera", val.camera, visit_sem{visit_sem_type::value});
    visitor("instance", val.instance, visit_sem{visit_sem_type::value});
    visitor("environment", val.environment, visit_sem{visit_sem_type::value});
    visitor("frame", val.frame, visit_sem{visit_sem_type::value, -10, 10});
    visitor("translation", val.translation,
        visit_sem{visit_sem_type::value, -10, 10});
    visitor("rotation", val.rotation, visit_sem{visit_sem_type::value});
    visitor("scaling", val.scaling, visit_sem{visit_sem_type::value});
}

/// Visit struct elements.
template <typename Visitor>
inline void visit(test_animation_params& val, Visitor&& visitor) {
    visitor("name", val.name, visit_sem{visit_sem_type::name});
    visitor("bezier", val.bezier, visit_sem{visit_sem_type::value});
    visitor("speed", val.speed, visit_sem{visit_sem_type::value});
    visitor("scale", val.scale, visit_sem{visit_sem_type::value, 0.01f, 10});
    visitor("times", val.times, visit_sem{visit_sem_type::value});
    visitor("translation", val.translation,
        visit_sem{visit_sem_type::value, -10, 10});
    visitor("rotation", val.rotation, visit_sem{visit_sem_type::value});
    visitor(
        "scaling", val.scaling, visit_sem{visit_sem_type::value, 0.01f, 10});
    visitor("nodes", val.nodes, visit_sem{visit_sem_type::value});
}

/// Visit struct elements.
template <typename Visitor>
inline void visit(test_scene_params& val, Visitor&& visitor) {
    visitor("name", val.name, visit_sem{visit_sem_type::name});
    visitor("cameras", val.cameras, visit_sem{visit_sem_type::value});
    visitor("shapes", val.shapes, visit_sem{visit_sem_type::value});
    visitor("instances", val.instances, visit_sem{visit_sem_type::value});
    visitor("environments", val.environments, visit_sem{visit_sem_type::value});
    visitor("materials", val.materials, visit_sem{visit_sem_type::value});
    visitor("textures", val.textures, visit_sem{visit_sem_type::value});
    visitor("nodes", val.nodes, visit_sem{visit_sem_type::value});
    visitor("animations", val.animations, visit_sem{visit_sem_type::value});
}

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// PATH TRACING SUPPORT FUNCTION
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup trace_support Path-tracing support
/// @{

/// Phong exponent to roughness.
float specular_exponent_to_roughness(float n);

/// Specular to fresnel eta.
void specular_fresnel_from_ks(const vec3f& ks, vec3f& es, vec3f& esk);

/// Compute the fresnel term for dielectrics. Implementation from
/// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
vec3f fresnel_dielectric(float cosw, const vec3f& eta_);

/// Compute the fresnel term for metals. Implementation from
/// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
vec3f fresnel_metal(float cosw, const vec3f& eta, const vec3f& etak);

/// Schlick approximation of Fresnel term.
vec3f fresnel_schlick(const vec3f& ks, float cosw);

/// Schlick approximation of Fresnel term weighted by roughness.
/// This is a hack, but works better than not doing it.
vec3f fresnel_schlick(const vec3f& ks, float cosw, float rs);

/// Evaluates the GGX distribution and geometric term.
float eval_ggx(float rs, float ndh, float ndi, float ndo);

/// Sample the GGX distribution.
vec3f sample_ggx(float rs, const vec2f& rn);

/// Evaluates the GGX pdf.
float sample_ggx_pdf(float rs, float ndh);

/// Triangle filter. Ppublic domain from stb_image_resize.
inline float filter_triangle(float x) {
    x = (float)fabs(x);
    if (x <= 1.0f) return 1 - x;
    return 0;
}
/// Cubic filter. Ppublic domain from stb_image_resize.
inline float filter_cubic(float x) {
    x = (float)fabs(x);
    if (x < 1.0f) return (4 + x * x * (3 * x - 6)) / 6;
    if (x < 2.0f) return (8 + x * (-12 + x * (6 - x))) / 6;
    return 0.0f;
}
/// Catmull-rom filter. Ppublic domain from stb_image_resize.
inline float filter_catmullrom(float x) {
    x = (float)fabs(x);
    if (x < 1.0f) return 1 - x * x * (2.5f - 1.5f * x);
    if (x < 2.0f) return 2 - x * (4 + x * (0.5f * x - 2.5f));
    return 0.0f;
}
/// Mitchell filter. Ppublic domain from stb_image_resize.
inline float filter_mitchell(float x) {
    x = (float)fabs(x);
    if (x < 1.0f) return (16 + x * x * (21 * x - 36)) / 18;
    if (x < 2.0f) return (32 + x * (-60 + x * (36 - 7 * x))) / 18;
    return 0.0f;
}

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// PATH TRACING
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup trace Path tracing
/// @{

/// Type of rendering algorithm.
enum struct trace_shader_type {
    /// Pathtrace.
    pathtrace = 0,
    /// Eye light for quick previews.
    eyelight,
    /// Direct illumination.
    direct,
    /// Pathtrace without MIS, usedful ony for debugging.
    pathtrace_nomis,
    /// Debug normal.
    debug_normal,
    /// Debug albedo.
    debug_albedo,
    /// Debug texcoord.
    debug_texcoord,
};

/// Random number generator type.
enum struct trace_rng_type {
    /// Uniform random numbers.
    uniform = 0,
    /// Stratified random numbers.
    stratified,
};

/// Filter type.
enum struct trace_filter_type {
    /// Box filter.
    box = 1,
    /// Hat filter.
    triangle = 2,
    /// Cubic spline.
    cubic = 3,
    /// Catmull-Rom spline.
    catmull_rom = 4,
    /// Mitchell-Netrevalli.
    mitchell = 5
};

/// Rendering params.
struct trace_params {
    /// Image width.
    int width = 360;
    /// Image height.
    int height = 360;
    /// Number of samples.
    int nsamples = 256;
    /// Sampler type.
    trace_shader_type stype = trace_shader_type::pathtrace;
    /// Wheter to test transmission in shadows.
    bool shadow_notransmission = false;
    /// Random number generation type.
    trace_rng_type rtype = trace_rng_type::stratified;
    /// Filter type.
    trace_filter_type ftype = trace_filter_type::box;
    /// Ambient lighting.
    vec3f ambient = {0, 0, 0};
    /// View environment map.
    bool envmap_invisible = false;
    /// Minimum ray depth.
    int min_depth = 3;
    /// Maximum ray depth.
    int max_depth = 8;
    /// Final pixel clamping.
    float pixel_clamp = 10;
    /// Ray intersection epsilon.
    float ray_eps = 1e-4f;
    /// Parallel execution.
    bool parallel = true;
    /// Seed for the random number generators.
    uint32_t seed = 0;
    /// Block size for parallel batches (probably leave it as is).
    int block_size = 32;
    /// Batch size for progressive rendering.
    int batch_size = 16;
};

/// Trace pixel state. Handles image accumulation and random number generation
/// for uniform and stratified sequences. The members are not part of the
/// the public API.
struct trace_pixel {
    /// Accumulated radiance.
    vec3f col = zero3f;
    /// Accumulated coverage.
    float alpha = 1;
    /// Random number state.
    rng_pcg32 rng = rng_pcg32();
    /// Pixel coordinates.
    int i = 0, j = 0;
    /// Number of samples computed.
    int sample = 0;
    /// Current dimension.
    int dimension = 0;
    /// Pixel weight for filtering.
    float weight = 0;
};

/// Trace light as either instances or environments. The members are not part of
/// the the public API.
struct trace_light {
    /// Instance pointer for instance lights.
    const instance* ist = nullptr;
    /// Environment pointer for environment lights.
    const environment* env = nullptr;
};

/// Trace lights. Handles sampling of illumination. The members are not part of
/// the the public API.
struct trace_lights {
    /// Shape instances.
    std::vector<trace_light> lights;
    /// Shape cdfs.
    std::unordered_map<const shape*, std::vector<float>> shape_cdfs;
    /// Shape areas.
    std::unordered_map<const shape*, float> shape_areas;
    /// Check whether there are any lights.
    bool empty() const { return lights.empty(); }
    /// Number of lights.
    int size() const { return (int)lights.size(); }
};

/// Initialize trace pixels.
image<trace_pixel> make_trace_pixels(const trace_params& params);
/// Initialize trace lights.
trace_lights make_trace_lights(const scene* scn);

/// Trace the next `nsamples` samples.
void trace_samples(const scene* scn, const camera* cam, const bvh_tree* bvh,
    const trace_lights& lights, image4f& img, image<trace_pixel>& pixels,
    int nsamples, const trace_params& params);

/// Trace the next `nsamples` samples with image filtering.
void trace_samples_filtered(const scene* scn, const camera* cam,
    const bvh_tree* bvh, const trace_lights& lights, image4f& img,
    image<trace_pixel>& pixels, int nsamples, const trace_params& params);

/// Trace the whole image.
inline image4f trace_image(const scene* scn, const camera* cam,
    const bvh_tree* bvh, const trace_params& params) {
    auto img = image4f(params.width, params.height);
    auto pixels = make_trace_pixels(params);
    auto lights = make_trace_lights(scn);
    trace_samples(scn, cam, bvh, lights, img, pixels, params.nsamples, params);
    return img;
}

/// Starts an anyncrhounous renderer.
void trace_async_start(const scene* scn, const camera* cam, const bvh_tree* bvh,
    const trace_lights& lights, image4f& img, image<trace_pixel>& pixels,
    std::vector<std::thread>& threads, bool& stop_flag,
    const trace_params& params);
/// Stop the asynchronous renderer.
void trace_async_stop(std::vector<std::thread>& threads, bool& stop_flag);

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// PATH TRACING TYPE SUPPORT
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup trace_type Path tracing type support
/// @{

/// Names of enum values.
template <>
inline const std::vector<std::pair<std::string, trace_shader_type>>&
enum_names<trace_shader_type>() {
    static auto names = std::vector<std::pair<std::string, trace_shader_type>>{
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

/// Names of enum values.
template <>
inline const std::vector<std::pair<std::string, trace_rng_type>>&
enum_names<trace_rng_type>() {
    static auto names = std::vector<std::pair<std::string, trace_rng_type>>{
        {"uniform", trace_rng_type::uniform},
        {"stratified", trace_rng_type::stratified}};
    return names;
}

/// Names of enum values.
template <>
inline const std::vector<std::pair<std::string, trace_filter_type>>&
enum_names<trace_filter_type>() {
    static auto names = std::vector<std::pair<std::string, trace_filter_type>>{
        {"box", trace_filter_type::box},
        {"triangle", trace_filter_type::triangle},
        {"cubic", trace_filter_type::cubic},
        {"catmull-rom", trace_filter_type::catmull_rom},
        {"mitchell", trace_filter_type::mitchell},
    };
    return names;
}

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// WAVEFRONT OBJ SUPPORT
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup obj Wavefront OBJ
/// @{

/// Obj face vertex.
struct obj_vertex {
    /// Position.
    int pos;
    /// Texcoord.
    int texcoord;
    /// Normal.
    int norm;
    /// Color [extension].
    int color;
    /// Radius [extension].
    int radius;

    /// Element constructor. Initializes all non-specified members as -1.
    obj_vertex(int pos = -1, int texcoord = -1, int norm = -1, int color = -1,
        int radius = -1)
        : pos(pos)
        , texcoord(texcoord)
        , norm(norm)
        , color(color)
        , radius(radius) {}
};

// Comparison for unordred_map.
inline bool operator==(const obj_vertex& a, const obj_vertex& b) {
    return a.pos == b.pos && a.texcoord == b.texcoord && a.norm == b.norm &&
           a.color == b.color && a.radius == b.radius;
}

/// Obj element type.
enum struct obj_element_type : uint16_t {
    /// List of points.
    point = 1,
    /// Polyline.
    line = 2,
    /// Polygon face.
    face = 3,
    /// Bezier segments.
    bezier = 4,
};

/// Obj element vertex indices.
struct obj_element {
    /// Starting vertex index.
    uint32_t start;
    /// Element type.
    obj_element_type type;
    /// Number of vertices.
    uint16_t size;
};

/// Obj element group.
struct obj_group {
    /// Material name.
    std::string matname;
    /// Group name.
    std::string groupname;
    /// Smoothing.
    bool smoothing = true;
    /// Element vertices.
    std::vector<obj_vertex> verts;
    /// Element faces.
    std::vector<obj_element> elems;
    /// Properties not explicitly handled [extension].
    std::unordered_map<std::string, std::vector<std::string>> props;
};

/// Obj object.
struct obj_object {
    /// Name.
    std::string name;
    /// Element groups.
    std::vector<obj_group*> groups;
    /// Properties not explicitly handled [extension].
    std::unordered_map<std::string, std::vector<std::string>> props;

    /// Cleanup.
    ~obj_object();
};

/// Obj texture information.
struct obj_texture_info {
    /// File path.
    std::string path = "";
    /// Whether to clamp to the edge.
    bool clamp = false;
    /// Scale for bump and displacement.
    float scale = 1;
    /// Properties not explicitly handled.
    std::unordered_map<std::string, std::vector<std::string>> props;
};

// Comparison for texture info.
inline bool operator==(const obj_texture_info& a, const obj_texture_info& b) {
    if (a.path.empty() && b.path.empty()) return true;
    if (a.path != b.path) return false;
    return a.clamp == b.clamp && a.scale == b.scale && a.props == b.props;
}

/// Obj texture. Texture data is loaded only if desired.
struct obj_texture {
    /// File path.
    std::string path;
    /// Width.
    int width = 0;
    /// Height.
    int height = 0;
    /// Number of Channels.
    int ncomp = 0;
    /// Buffer data for LDR images.
    std::vector<uint8_t> datab;
    /// Buffer data for HDR images.
    std::vector<float> dataf;
};

/// Obj material.
struct obj_material {
    /// Name.
    std::string name;
    /// MTL illum mode.
    int illum = 0;

    /// Emission color.
    vec3f ke = {0, 0, 0};
    /// Ambient color.
    vec3f ka = {0, 0, 0};
    /// Diffuse color.
    vec3f kd = {0, 0, 0};
    /// Specular color.
    vec3f ks = {0, 0, 0};
    /// Reflection color.
    vec3f kr = {0, 0, 0};
    /// Transmision color.
    vec3f kt = {0, 0, 0};
    /// Phong exponent for ks.
    float ns = 1;
    /// Index of refraction.
    float ior = 1;
    /// Opacity.
    float op = 1;

    /// Emission texture.
    obj_texture_info ke_txt;
    /// Ambient texture.
    obj_texture_info ka_txt;
    /// Diffuse texture.
    obj_texture_info kd_txt;
    /// Specular texture.
    obj_texture_info ks_txt;
    /// Reflection texture.
    obj_texture_info kr_txt;
    /// Transmission texture.
    obj_texture_info kt_txt;
    /// Specular exponent texture.
    obj_texture_info ns_txt;
    /// Opacity texture.
    obj_texture_info op_txt;
    /// Index of refraction.
    obj_texture_info ior_txt;
    /// Bump map texture (heighfield).
    obj_texture_info bump_txt;
    /// Displacement map texture (heighfield).
    obj_texture_info disp_txt;
    /// Normal map texture.
    obj_texture_info norm_txt;

    /// Properties not explicitly handled.
    std::unordered_map<std::string, std::vector<std::string>> props;
};

/// Obj camera [extension].
struct obj_camera {
    /// Camera name.
    std::string name;
    /// Transform frame (affine matrix).
    frame3f frame = identity_frame3f;
    /// Orthografic camera.
    bool ortho = false;
    /// Vertical field of view.
    float yfov = 2 * atan(0.5f);
    /// Aspect ratio.
    float aspect = 16.0f / 9.0f;
    /// Lens aperture.
    float aperture = 0;
    /// Focus distance.
    float focus = 1;
};

/// Obj environment [extension].
struct obj_environment {
    /// Environment name.
    std::string name;
    /// Transform frame (affine matrix).
    frame3f frame = identity_frame3f;
    /// Material name.
    std::string matname;
};

/// Obj node [extension].
struct obj_node {
    /// Node name.
    std::string name;
    /// Node parent.
    std::string parent;
    /// Camera name.
    std::string camname;
    /// Instance name.
    std::string objname;
    /// Environment name.
    std::string envname;
    /// Transform frame (affine matrix).
    frame3f frame = identity_frame3f;
    /// Translation.
    vec3f translation = zero3f;
    /// Rotation.
    quat4f rotation = {0, 0, 0, 1};
    /// Scaling.
    vec3f scaling = {1, 1, 1};
};

/// Obj scene.
struct obj_scene {
    /// Vertex positions.
    std::vector<vec3f> pos;
    /// Vertex normals.
    std::vector<vec3f> norm;
    /// Vertex texcoord.
    std::vector<vec2f> texcoord;
    /// Vertex color [extension].
    std::vector<vec4f> color;
    /// Vertex radius [extension].
    std::vector<float> radius;

    /// Objects.
    std::vector<obj_object*> objects;
    /// Materials.
    std::vector<obj_material*> materials;
    /// Textures.
    std::vector<obj_texture*> textures;
    /// Cameras [extension].
    std::vector<obj_camera*> cameras;
    /// Environments [extension].
    std::vector<obj_environment*> environments;
    /// Nodes [extension].
    std::vector<obj_node*> nodes;

    /// Cleanup.
    ~obj_scene();
};

/// Load an OBJ from file `filename`. Load textures if `load_textures` is true,
/// and report errors only if `skip_missing` is false.
/// Texture coordinates and material Tr are flipped if `flip_texcoord` and
/// `flip_tp` are respectively true.
obj_scene* load_obj(const std::string& filename, bool load_textures = false,
    bool skip_missing = false, bool flip_texcoord = true, bool flip_tr = true);

/// Save an OBJ to file `filename`. Save textures if `save_textures` is true,
/// and report errors only if `skip_missing` is false.
/// Texture coordinates and material Tr are flipped if `flip_texcoord` and
/// `flip_tp` are respectively true.
void save_obj(const std::string& filename, const obj_scene* model,
    bool save_textures = false, bool skip_missing = false,
    bool flip_texcoord = true, bool flip_tr = true);

/// @}

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

/// @defgroup gltf Khronos glTF
/// @{

/// Generic buffer data.
using buffer_data = std::vector<unsigned char>;

/// Generic image data.
struct image_data {
    /// Width.
    int width = 0;
    /// Height.
    int height = 0;
    /// Number of Channels.
    int ncomp = 0;
    /// Buffer data for 8-bit images.
    std::vector<uint8_t> datab;
    /// Buffer data for float images.
    std::vector<float> dataf;
};

/// Id for glTF references.
template <typename T>
struct glTFid {
    /// Defaoult constructor to an invalid id.
    glTFid() : _id(-1) {}
    /// Explicit conversion from integer.
    explicit glTFid(int id) : _id(id) {}
    /// Explicit convcersion to integer.
    explicit operator int() const { return _id; }
    /// Check if it is valid.
    bool is_valid() const { return _id >= 0; }
    /// Check if it is valid.
    explicit operator bool() const { return _id >= 0; }

   private:
    // id
    int _id = -1;
};

/// Generic glTF object.
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
    std::string name = "";
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
    std::vector<float> max = {};
    /// Minimum value of each component in this attribute.
    std::vector<float> min = {};
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
    std::vector<glTFAnimationChannel*> channels = {};
    /// An array of samplers that combines input and output accessors with an
    /// interpolation algorithm to define a keyframe graph (but not its target).
    /// [required]
    std::vector<glTFAnimationSampler*> samplers = {};

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
    std::string copyright = "";
    /// Tool that generated this glTF model.  Useful for debugging.
    std::string generator = "";
    /// The glTF version that this asset targets. [required]
    std::string version = "";
    /// The minimum glTF version that this asset targets.
    std::string minVersion = "";
};

/// A buffer points to binary geometry, animation, or skins.
struct glTFBuffer : glTFChildOfRootProperty {
    /// The uri of the buffer.
    std::string uri = "";
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
    std::string uri = "";
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
    std::map<std::string, glTFid<glTFAccessor>> attributes = {};
    /// The index of the accessor that contains the indices.
    glTFid<glTFAccessor> indices = {};
    /// The index of the material to apply to this primitive when rendering.
    glTFid<glTFMaterial> material = {};
    /// The type of primitives to render.
    glTFMeshPrimitiveMode mode = glTFMeshPrimitiveMode::Triangles;
    /// An array of Morph Targets, each  Morph Target is a dictionary mapping
    /// attributes (only `POSITION`, `NORMAL`, and `TANGENT` supported) to their
    /// deviations in the Morph Target.
    std::vector<std::map<std::string, glTFid<glTFAccessor>>> targets = {};
};

/// A set of primitives to be rendered.  A node can contain one mesh.  A node's
/// transform places the mesh in the scene.
struct glTFMesh : glTFChildOfRootProperty {
    /// An array of primitives, each defining geometry to be rendered with a
    /// material. [required]
    std::vector<glTFMeshPrimitive*> primitives = {};
    /// Array of weights to be applied to the Morph Targets.
    std::vector<float> weights = {};

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
    std::vector<glTFid<glTFNode>> children = {};
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
    std::vector<float> weights = {};
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
    std::vector<glTFid<glTFNode>> nodes = {};
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
    std::vector<glTFid<glTFNode>> joints = {};
};

/// The root object for a glTF asset.
struct glTF : glTFProperty {
    /// Names of glTF extensions used somewhere in this asset.
    std::vector<std::string> extensionsUsed = {};
    /// Names of glTF extensions required to properly load this asset.
    std::vector<std::string> extensionsRequired = {};
    /// An array of accessors.
    std::vector<glTFAccessor*> accessors = {};
    /// An array of keyframe animations.
    std::vector<glTFAnimation*> animations = {};
    /// Metadata about the glTF asset. [required]
    glTFAsset* asset = nullptr;
    /// An array of buffers.
    std::vector<glTFBuffer*> buffers = {};
    /// An array of bufferViews.
    std::vector<glTFBufferView*> bufferViews = {};
    /// An array of cameras.
    std::vector<glTFCamera*> cameras = {};
    /// An array of images.
    std::vector<glTFImage*> images = {};
    /// An array of materials.
    std::vector<glTFMaterial*> materials = {};
    /// An array of meshes.
    std::vector<glTFMesh*> meshes = {};
    /// An array of nodes.
    std::vector<glTFNode*> nodes = {};
    /// An array of samplers.
    std::vector<glTFSampler*> samplers = {};
    /// The index of the default scene.
    glTFid<glTFScene> scene = {};
    /// An array of scenes.
    std::vector<glTFScene*> scenes = {};
    /// An array of skins.
    std::vector<glTFSkin*> skins = {};
    /// An array of textures.
    std::vector<glTFTexture*> textures = {};

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

/// Load a gltf file `filename` from disk. Load binaries and images only if
/// `load_bin` and `load_img` are true, reporting errors only if `skip_missing`
/// is false.
glTF* load_gltf(const std::string& filename, bool load_bin = true,
    bool load_img = false, bool skip_missing = false);

/// Load a binary gltf file `filename` from disk. Load binaries and images only
/// if `load_bin` and `load_img` are true, reporting errors only if
/// `skip_missing` is false.
glTF* load_binary_gltf(const std::string& filename, bool load_bin = true,
    bool load_img = false, bool skip_missing = false);

/// Save a gltf file `filename` to disk. Save binaries and images only if
/// `save_bin` and `save_img` are true.
void save_gltf(const std::string& filename, const glTF* gltf,
    bool save_bin = true, bool save_img = false);

/// Save a gltf file `filename` to disk. Save binaries and images only if
/// `save_bin` and `save_img` are true.
void save_binary_gltf(const std::string& filename, const glTF* gltf,
    bool save_bin = true, bool save_img = false);

/// Computes the local node transform and its inverse.
inline mat4f node_transform(const glTFNode* node) {
    return frame_to_mat(translation_frame(node->translation) *
                        rotation_frame(node->rotation) *
                        scaling_frame(node->scale)) *
           node->matrix;
}

/// A view for gltf array buffers that allows for typed access.
struct accessor_view {
    /// Construct a view from an accessor.
    accessor_view(const glTF* gltf, const glTFAccessor* accessor);

    /// Number of elements in the view.
    int size() const { return _size; }
    /// Number of elements in the view
    int count() const { return _size; }
    /// Number of components per element
    int ncomp() const { return _ncomp; }
    /// Check whether the view is valid.
    bool valid() const { return _valid; }

    /// Get the idx-th element of fixed length width default values.
    vec2f getv2f(int idx, const vec2f& def = {0, 0}) const {
        auto v = def;
        for (auto i = 0; i < min(_ncomp, 2); i++) v[i] = get(idx, i);
        return v;
    }
    /// Get the idx-th element of fixed length width default values.
    vec3f getv3f(int idx, const vec3f& def = {0, 0, 0}) const {
        auto v = def;
        for (auto i = 0; i < min(_ncomp, 3); i++) v[i] = get(idx, i);
        return v;
    }
    /// Get the idx-th element of fixed length width default values.
    vec4f getv4f(int idx, const vec4f& def = {0, 0, 0, 0}) const {
        auto v = def;
        for (auto i = 0; i < min(_ncomp, 4); i++) v[i] = get(idx, i);
        return v;
    }

    /// Get the idx-th element of fixed length as a matrix.
    mat4f getm4f(int idx) const {
        auto v = mat4f();
        assert(_ncomp == 16);
        for (auto j = 0; j < 4; j++)
            for (auto i = 0; i < 4; i++) v[j][i] = get(idx, j * 4 + i);
        return v;
    }

    /// Get the c-th component of the idx-th element.
    float get(int idx, int c = 0) const;

    /// Get the idx-th element as integer with fixed length.
    vec2i getv2i(int idx, const vec2i& def = {0, 0}) const {
        auto v = def;
        for (auto i = 0; i < min(_ncomp, 2); i++) { v[i] = geti(idx, i); }
        return v;
    }
    /// Get the idx-th element as integer with fixed length.
    vec3i getv3i(int idx, const vec3i& def = {0, 0, 0}) const {
        auto v = def;
        for (auto i = 0; i < min(_ncomp, 3); i++) { v[i] = geti(idx, i); }
        return v;
    }
    /// Get the idx-th element as integer with fixed length.
    vec4i getv4i(int idx, const vec4i& def = {0, 0, 0, 0}) const {
        auto v = def;
        for (auto i = 0; i < min(_ncomp, 4); i++) { v[i] = geti(idx, i); }
        return v;
    }

    /// Get the c-th component of the idx-th element as integer.
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

/// @}

}  // namespace ygl

#endif

#if YGL_SVG

// -----------------------------------------------------------------------------
// SVG SUPPORT
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup svg Svg
/// @{

/// Svg path.
struct svg_path {
    /// Path vertices.
    std::vector<vec2f> pos;
};

/// Svg shape.
struct svg_shape {
    /// Paths.
    std::vector<svg_path*> paths;

    /// Cleanup.
    ~svg_shape() {
        for (auto e : paths) delete e;
    }
};

/// Svg scene.
struct svg_scene {
    /// Shapes
    std::vector<svg_shape*> shapes;

    /// Cleanup.
    ~svg_scene() {
        for (auto e : shapes) delete e;
    }
};

/// Load an SVG.
svg_scene* load_svg(const std::string& filename);

/// Save an SVG.
void save_svg(const std::string& filename, const svg_scene* svg);

/// @}

}  // namespace ygl

#endif

// -----------------------------------------------------------------------------
// PYTHON-LIKE STRING, PATH AND FILE OPERATIONS
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup string_ops String, path and file functions
/// @{

/// Checks if a string starts with a prefix.
inline bool startswith(const std::string& str, const std::string& substr) {
    if (str.length() < substr.length()) return false;
    for (auto i = 0; i < substr.length(); i++)
        if (str[i] != substr[i]) return false;
    return true;
}

/// Checks if a string ends with a prefix.
inline bool endswith(const std::string& str, const std::string& substr) {
    if (str.length() < substr.length()) return false;
    auto offset = str.length() - substr.length();
    for (auto i = 0; i < substr.length(); i++)
        if (str[i + offset] != substr[i]) return false;
    return true;
}

/// Check is a string contains a substring.
inline bool contains(const std::string& str, const std::string& substr) {
    return str.find(substr) != str.npos;
}

/// Splits a string into lines at the '\n' character. The line
/// terminator is kept if keep_newline. This function does not work on
/// Window if keep_newline is true.
inline std::vector<std::string> splitlines(
    const std::string& str, bool keep_newline = false) {
    if (str.empty()) return {};
    auto lines = std::vector<std::string>();
    auto line = std::vector<char>();
    for (auto c : str) {
        if (c == '\n') {
            if (keep_newline) line.push_back(c);
            lines.push_back(std::string(line.begin(), line.end()));
            line.clear();
        } else {
            line.push_back(c);
        }
    }
    if (!line.empty()) lines.push_back(std::string(line.begin(), line.end()));
    return lines;
}

/// Partition the string.
inline std::vector<std::string> partition(
    const std::string& str, const std::string& split) {
    auto pos = str.find(split);
    if (pos == str.npos) return {str, "", ""};
    return {str.substr(0, pos), split, str.substr(pos + split.length())};
}

/// Splits the string.
inline std::vector<std::string> split(const std::string& str) {
    if (str.empty()) return {};
    auto ret = std::vector<std::string>();
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
inline std::vector<std::string> split(
    const std::string& str, const std::string& substr) {
    if (str.empty()) return {};
    auto ret = std::vector<std::string>();
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
inline std::vector<std::string> split(const std::string& str, char substr) {
    if (str.empty()) return {};
    auto ret = std::vector<std::string>();
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
inline std::string rstrip(const std::string& str) {
    auto pos = str.find_last_not_of(" \t\r\n");
    if (pos == str.npos) return "";
    return str.substr(0, pos + 1);
}

/// Strip the string.
inline std::string lstrip(const std::string& str) {
    auto pos = str.find_first_not_of(" \t\r\n");
    if (pos == str.npos) return "";
    return str.substr(pos);
}

/// Strip the string.
inline std::string strip(const std::string& str) { return rstrip(lstrip(str)); }

/// Joins a list of string with a string as separator.
inline std::string join(
    const std::vector<std::string>& strs, const std::string& sep) {
    auto ret = std::string();
    auto first = true;
    for (auto& str : strs) {
        if (!first) ret += sep;
        ret += str;
        first = false;
    }
    return ret;
}

/// Converts an ASCII string to lowercase.
inline std::string lower(const std::string& str) {
    auto s = str;
    for (auto& c : s) c = tolower(c);
    return s;
}

/// Converts an ASCII string to uppercase.
inline std::string upper(const std::string& str) {
    auto s = str;
    for (auto& c : s) c = toupper(c);
    return s;
}

/// Check if a string is space.
inline bool isspace(const std::string& str) {
    for (auto c : str) {
        if (c != ' ' && c != '\n' && c != '\t' && c != '\r') return false;
    }
    return true;
}

/// Replace s1 with s2 in str.
inline std::string replace(
    const std::string& str, const std::string& s1, const std::string& s2) {
    auto s = std::string();
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
inline std::string path_dirname(const std::string& filename) {
    auto pos = filename.rfind('/');
    if (pos == std::string::npos) pos = filename.rfind('\\');
    if (pos == std::string::npos) return "";
    return filename.substr(0, pos + 1);
}

/// Get extension (including '.').
inline std::string path_extension(const std::string& filename) {
    auto pos = filename.rfind('.');
    if (pos == std::string::npos) return "";
    return filename.substr(pos);
}

/// Get file basename.
inline std::string path_basename(const std::string& filename) {
    auto dirname = path_dirname(filename);
    auto extension = path_extension(filename);
    return filename.substr(
        dirname.size(), filename.size() - dirname.size() - extension.size());
}

/// Get filename without directory (equiv to get_basename() +
/// get_extension()).
inline std::string path_filename(const std::string& filename) {
    return path_basename(filename) + path_extension(filename);
}

/// Replace extension.
inline std::string replace_path_extension(
    const std::string& filename, const std::string& ext) {
    return path_dirname(filename) + path_basename(filename) + ext;
}

/// Prepend a string to the extension.
inline std::string prepend_path_extension(
    const std::string& filename, const std::string& prep) {
    return path_dirname(filename) + path_basename(filename) + prep +
           path_extension(filename);
}

/// Splits a path calling the above functions.
inline void split_path(const std::string& filename, std::string& dirname,
    std::string& basename, std::string& ext) {
    dirname = path_dirname(filename);
    basename = path_basename(filename);
    ext = path_extension(filename);
}

/// Convert from Windows to Unix/OsX path separator
inline std::string path_convert_eparator(const std::string& path_) {
    auto path = path_;
    for (auto& c : path)
        if (c == '\\') c = '/';
    return path;
}

/// Really-minimal Python like string format. The implementation is not fast
/// nor memory efficient. But it is good enough for some needs.
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

// Implementation of the function below.
inline void _format_one(std::vector<std::string>& vals) {}
template <typename Arg, typename... Args>
inline void _format_one(
    std::vector<std::string>& vals, const Arg& arg, const Args&... args) {
    auto stream = std::stringstream();
    stream << arg;
    vals.push_back(stream.str());
    _format_one(vals, args...);
}

/// Really-minimal Python like string format. Internally uses streams for
/// generality and supports for now only the '{}' operator. The implementation
/// is not fast nor memory efficient. But it is good enough for some needs.
template <typename... Args>
inline std::string format(const std::string& fmt, const Args&... args) {
    auto vals = std::vector<std::string>();
    _format_one(vals, args...);
    return format(fmt, vals);
}

/// Wrapper for the above function that prints to stdout.
template <typename... Args>
inline void print(const std::string& fmt, const Args&... args) {
    printf("%s", format(fmt, args...).c_str());
}

/// Wrapper for the above function that prints to stdout with endline.
template <typename... Args>
inline void println(const std::string& fmt, const Args&... args) {
    printf("%s\n", format(fmt, args...).c_str());
}

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// FILE LOADING AND SAVING
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup file_io File loading and saving
/// @{

/// Loads the contents of a binary file in an in-memory array.
inline std::vector<unsigned char> load_binary(const std::string& filename) {
    // http://stackoverflow.com/questions/116038/what-is-the-best-way-to-read-an-entire-file-into-a-stdstring-in-c
    std::fstream fs(
        filename, std::ios_base::in | std::ios_base::binary | std::ios::ate);
    if (fs.fail()) throw std::runtime_error("cannot read file " + filename);
    fs.seekg(0, std::ios::end);
    auto buf = std::vector<unsigned char>(fs.tellg());
    fs.seekg(0);
    fs.read((char*)buf.data(), buf.size());
    if (fs.fail() || fs.bad())
        throw std::runtime_error("cannot read file " + filename);
    return buf;
}

/// Loads the contents of a text file into a string.
inline std::string load_text(const std::string& filename) {
    std::fstream fs(filename, std::ios_base::in);
    if (fs.fail()) throw std::runtime_error("cannot read file " + filename);
    std::stringstream ss;
    ss << fs.rdbuf();
    if (fs.fail()) throw std::runtime_error("cannot read file " + filename);
    return ss.str();
}

/// Saves binary data to a file.
inline void save_binary(
    const std::string& filename, const std::vector<unsigned char>& data) {
#if YGL_IOSTREAM
    fstream fs(filename, ios_base::out | ios_base::binary);
    if (fs.fail()) throw std::runtime_error("cannot write file " + filename);
    fs.write((const char*)data.data(), data.size());
    if (fs.fail() || fs.bad())
        throw std::runtime_error("cannot write file " + filename);
#else
    auto f = fopen(filename.c_str(), "wb");
    if (!f) throw std::runtime_error("cannot write file " + filename);
    fwrite(data.data(), 1, (int)data.size(), f);
    fclose(f);
#endif
}

/// Saves a string to a text file.
inline void save_text(const std::string& filename, const std::string& str) {
#if YGL_IOSTREAM
    fstream fs(filename, ios_base::out);
    if (fs.fail()) throw std::runtime_error("cannot write file " + filename);
    fs << str;
    if (fs.fail()) throw std::runtime_error("cannot write file " + filename);
#else
    auto f = fopen(filename.c_str(), "wt");
    if (!f) throw std::runtime_error("cannot write file " + filename);
    fwrite(str.c_str(), 1, (int)str.size(), f);
    fclose(f);
#endif
}

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMMEDIATE MODE COMMAND LINE PARSER
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup cmdline Immediate-mode command line parser
/// @{

/// Immediate mode command line parser (opaque type)
struct cmdline_parser;

/// Immediate mode command line parser. Members are not part of the public API.
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

/// Check unused arguments.
inline bool should_exit(cmdline_parser& parser);

/// Returns the usage string.
inline std::string get_usage(const cmdline_parser& parser);

/// Pase a flag from the command line.
inline bool parse_flag(cmdline_parser& parser, const std::string& name,
    const std::string& flag, const std::string& help, bool def = false,
    bool req = false);

/// Pase an option from the command line.
template <typename T>
inline T parse_opt(cmdline_parser& parser, const std::string& name,
    const std::string& flag, const std::string& help, const T& def = {},
    bool req = false, const std::vector<T>& choices = {});

/// Parse an enum option from the command line.
template <typename T>
inline T parse_opt(cmdline_parser& parser, const std::string& name,
    const std::string& flag, const std::string& help,
    const std::vector<std::pair<std::string, T>>& key_values, const T& def,
    bool req = false, const std::vector<T>& choices = {});

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

/// Initialize a command line parser.
inline cmdline_parser make_parser(
    int argc, char** argv, const std::string& prog, const std::string& help);

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// SIMPLE LOGGER
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup logger Simple logging
/// @{

/// Logger object. A logger can output messages to console an a file.
/// Members are not part of the public API.
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

/// Make a logger with an optional console stream, an optional file stram
/// and the specified verbosity level.
inline logger* make_logger(const std::string& filename = "",
    bool console = true, bool verbose = true, bool file_append = true) {
    auto lgr = new logger();
    lgr->_verbose = verbose;
    lgr->_console = console;
    if (filename.empty()) {
        lgr->_file = nullptr;
    } else {
        lgr->_file = fopen(filename.c_str(), (file_append) ? "at" : "wt");
        if (!lgr->_file)
            throw std::runtime_error("could not open file " + filename);
    }
    return lgr;
}

/// Get the default logger.
inline logger* get_default_logger() {
    static auto default_logger = new logger();
    return default_logger;
}

// Log a message. Used internally.
inline void _log_msg(logger* lgr, const std::string& msg, const char* type) {
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

/// Log an info message.
template <typename... Args>
inline void log_info(logger* lgr, const std::string& msg, const Args&... args) {
    if (!lgr->_verbose) return;
    _log_msg(lgr, format(msg, args...), "INFO ");
}

/// Log an info message.
template <typename... Args>
inline void log_warning(
    logger* lgr, const std::string& msg, const Args&... args) {
    if (!lgr->_verbose) return;
    _log_msg(lgr, format(msg, args...), "WARN ");
}

/// Log an error message.
template <typename... Args>
inline void log_error(
    logger* lgr, const std::string& msg, const Args&... args) {
    _log_msg(lgr, format(msg, args...), "ERROR");
}

/// Log a fatal message and exit.
template <typename... Args>
inline void log_fatal(
    logger* lgr, const std::string& msg, const Args&... args) {
    _log_msg(lgr, format(msg, args...), "FATAL");
    exit(1);
}

/// Logs a message to the default loggers.
template <typename... Args>
inline void log_info(const std::string& msg, const Args&... args) {
    log_info(get_default_logger(), msg, args...);
}

/// Logs a message to the default loggers.
template <typename... Args>
inline void log_error(const std::string& msg, const Args&... args) {
    log_error(get_default_logger(), msg, args...);
}

/// Logs a message to the default loggers.
template <typename... Args>
inline void log_fatal(const std::string& msg, const Args&... args) {
    log_fatal(get_default_logger(), msg, args...);
}

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// TIMER
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup timer Simple timer
/// @{

/// A simple wrapper for std::chrono.
struct timer {
    /// Initialize a timer and start it if necessary.
    timer(bool autostart = true) {
        if (autostart) start();
    }

    /// Start a timer.
    void start() {
        _start = std::chrono::steady_clock::now();
        _started = true;
    }

    /// Stops a timer.
    void stop() {
        _end = std::chrono::steady_clock::now();
        _started = false;
    }

    /// Elapsed time.
    double elapsed() {
        if (_started) stop();
        std::chrono::duration<double> diff = (_end - _start);
        return diff.count();
    }

   private:
    bool _started = false;
    std::chrono::time_point<std::chrono::steady_clock> _start, _end;
};

/// @}

}  // namespace ygl

#if YGL_OPENGL

// -----------------------------------------------------------------------------
// OPENGL FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup gl_util OpenGL utilities
/// @{

/// OpenGL shape element types.
enum struct gl_elem_type : int {
    /// Points.
    point = 1,
    /// Lines.
    line = 2,
    /// Triangles.
    triangle = 3
};

/// OpenGL light types.
enum struct gl_light_type : int {
    /// Point lights.
    point = 0,
    /// Directional lights.
    directional = 1,
};

/// OpenGL lights
struct gl_lights {
    /// light positions.
    std::vector<vec3f> pos;
    /// Light intensities.
    std::vector<vec3f> ke;
    /// Light types.
    std::vector<gl_light_type> type;
};

/// Checks for GL error and then prints.
bool gl_check_error(bool print = true);

/// Clear window.
void gl_clear_buffers(const vec4f& background = {0, 0, 0, 0});

/// Enable/disable depth test.
void gl_enable_depth_test(bool enabled);
/// Enable/disable culling.
void gl_enable_culling(bool enabled, bool front = false, bool back = true);
/// Enable/disable wireframe.
void gl_enable_wireframe(bool enabled);
/// Enable/disable blending.
void gl_enable_blending(bool enabled);
/// Set blending to over operator.
void gl_set_blend_over();

/// Line width.
void gl_line_width(float w);

/// Set viewport.
void gl_set_viewport(const vec4i& v);
/// Set viewport.
void gl_set_viewport(const vec2i& v);

/// Reads an image from the the framebuffer.
void gl_read_imagef(float* pixels, int w, int h, int nc);

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// OPENGL TEXTURE FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup gl_texture OpenGL textures
/// @{

/// OpenGL texture object. Members are not part of the public API.
struct gl_texture {
    // Texture id.
    uint tid = 0;
    // Width.
    int width = 0;
    // Height.
    int height = 0;
    // Ncomp.
    int ncomp = 0;
    // Stored as floats.
    bool as_float = false;
    // Stored as sRGB.
    bool as_srgb = true;
    // Mipmap creation.
    bool mipmap = true;
    // Linear interpolation.
    bool linear = true;
};

// Implementation of update_texture.
void _update_texture(gl_texture& txt, int w, int h, int nc, const void* pixels,
    bool floats, bool linear, bool mipmap, bool as_float, bool as_srgb);

/// Updates a texture with pixels values of size w, h with nc number of
/// components (1-4). Internally use float if as_float and filtering if filter.
inline void update_texture(gl_texture& txt, int w, int h, int nc,
    const float* pixels, bool linear, bool mipmap, bool as_float) {
    _update_texture(
        txt, w, h, nc, pixels, true, linear, mipmap, as_float, false);
}

/// Updates a texture with pixels values of size w, h with nc number of
/// components (1-4). Internally use float if as_float and filtering if filter.
inline void update_texture(gl_texture& txt, int w, int h, int nc,
    const unsigned char* pixels, bool linear, bool mipmap, bool as_srgb) {
    _update_texture(
        txt, w, h, nc, pixels, false, linear, mipmap, false, as_srgb);
}

/// Updates a texture with pixels values from an image.
/// Internally use float if as_float and filtering if filter.
inline void update_texture(gl_texture& txt, const image4f& img, bool linear,
    bool mipmap, bool as_float) {
    update_texture(txt, img.width(), img.height(), 4, (const float*)data(img),
        linear, mipmap, as_float);
}

/// Updates a texture with pixels values from an image.
/// Internally use float if as_float and filtering if filter.
inline void update_texture(gl_texture& txt, const image4b& img, bool linear,
    bool mipmap, bool as_srgb) {
    update_texture(txt, img.width(), img.height(), 4,
        (const unsigned char*)data(img), linear, mipmap, as_srgb);
}

/// Updates a texture with pixels values from an image.
inline void update_texture(gl_texture& txt, const image4f& img) {
    update_texture(txt, img, txt.linear, txt.mipmap, txt.as_float);
}

/// Updates a texture with pixels values from an image.
inline void update_texture(gl_texture& txt, const image4b& img) {
    update_texture(txt, img, txt.linear, txt.mipmap, txt.as_srgb);
}

/// Creates a texture from an image. Convenience wrapper to update_texture().
inline gl_texture make_texture(
    const image4f& img, bool linear, bool mipmap, bool as_float) {
    auto txt = gl_texture();
    update_texture(txt, img, linear, mipmap, as_float);
    return txt;
}

/// Creates a texture from an image. Convenience wrapper to update_texture().
inline gl_texture make_texture(
    const image4b& img, bool linear, bool mipmap, bool as_srgb) {
    auto txt = gl_texture();
    update_texture(txt, img, linear, mipmap, as_srgb);
    return txt;
}

/// Binds a texture to a texture unit.
void bind_texture(const gl_texture& txt, uint unit);

/// Unbinds a texture.
void unbind_texture(const gl_texture& txt, uint unit);

/// Get texture id.
inline uint get_texture_id(const gl_texture& txt) { return txt.tid; }

/// Check if defined.
inline bool is_texture_valid(const gl_texture& txt) { return (bool)txt.tid; }

/// Destroys the texture tid.
void clear_texture(gl_texture& txt);

/// Wrap values for OpenGL texture.
enum struct gl_texture_wrap {
    /// Not set.
    not_set = 0,
    /// Repeat.
    repeat = 1,
    /// Clamp to edge.
    clamp = 2,
    /// Mirror.
    mirror = 3,
};

/// Filter values for OpenGL texture.
enum struct gl_texture_filter {
    /// Not set.
    not_set = 0,
    /// Linear.
    linear = 1,
    /// Nearest.
    nearest = 2,
    /// Mip-mapping.
    linear_mipmap_linear = 3,
    /// Mip-mapping.
    nearest_mipmap_nearest = 4,
    /// Mip-mapping.
    linear_mipmap_nearest = 5,
    /// Mip-mapping.
    nearest_mipmap_linear = 6,
};

/// OpenGL texture parameters.
struct gl_texture_info {
    /// Texture.
    gl_texture txt = {};
    /// Texture coordinate set.
    int texcoord = 0;
    /// Texture strength/scale (used by some models).
    float scale = 1;
    /// Wrap s mode.
    gl_texture_wrap wrap_s = gl_texture_wrap::not_set;
    /// Wrap t mode
    gl_texture_wrap wrap_t = gl_texture_wrap::not_set;
    /// Filter mag mode.
    gl_texture_filter filter_mag = gl_texture_filter::not_set;
    /// Filter min mode.
    gl_texture_filter filter_min = gl_texture_filter::not_set;

    /// Default constructor.
    gl_texture_info() {}
    /// Constructor from texture id only.
    gl_texture_info(const gl_texture& tid) : txt(tid) {}
};

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// OPENGL VERTEX ARRAY BUFFER
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup gl_vertex_buffer OpenGL vertex array buffers
/// @{

/// OpenGL vertex/element buffer. Members are not part of the public API.
struct gl_vertex_buffer {
    // Buffer id.
    uint bid = 0;
    // Number of elements.
    int num = 0;
    // Number of components.
    int ncomp = 0;
    // Whether it is floats.
    bool as_float = true;
};

// Updates the bufferwith new data.
void _update_vertex_buffer(gl_vertex_buffer& buf, int n, int nc,
    const void* values, bool as_float, bool dynamic);

/// Updates the buffer with new data.
inline void update_vertex_buffer(gl_vertex_buffer& buf, int num, int ncomp,
    const float* values, bool dynamic = false) {
    _update_vertex_buffer(buf, num, ncomp, values, true, dynamic);
}

/// Updates the buffer with new data.
inline void update_vertex_buffer(gl_vertex_buffer& buf, int num, int ncomp,
    const int* values, bool dynamic = false) {
    _update_vertex_buffer(buf, num, ncomp, values, false, dynamic);
}

/// Updates the buffer with new data.
inline void update_vertex_buffer(gl_vertex_buffer& buf,
    const std::vector<float>& values, bool dynamic = false) {
    update_vertex_buffer(buf, values.size(), 1, values.data(), dynamic);
}

/// Updates the buffer with new data.
inline void update_vertex_buffer(gl_vertex_buffer& buf,
    const std::vector<vec2f>& values, bool dynamic = false) {
    update_vertex_buffer(
        buf, values.size(), 2, (const float*)values.data(), dynamic);
}

/// Updates the buffer with new data.
inline void update_vertex_buffer(gl_vertex_buffer& buf,
    const std::vector<vec3f>& values, bool dynamic = false) {
    update_vertex_buffer(
        buf, values.size(), 3, (const float*)values.data(), dynamic);
}

/// Updates the buffer with new data.
inline void update_vertex_buffer(gl_vertex_buffer& buf,
    const std::vector<vec4f>& values, bool dynamic = false) {
    update_vertex_buffer(
        buf, values.size(), 4, (const float*)values.data(), dynamic);
}

/// Updates the buffer with new data.
inline void update_vertex_buffer(gl_vertex_buffer& buf,
    const std::vector<int>& values, bool dynamic = false) {
    update_vertex_buffer(buf, values.size(), 1, values.data(), dynamic);
}

/// Updates the buffer with new data.
inline void update_vertex_buffer(gl_vertex_buffer& buf,
    const std::vector<vec2i>& values, bool dynamic = false) {
    update_vertex_buffer(
        buf, values.size(), 2, (const int*)values.data(), dynamic);
}

/// Updates the buffer with new data.
inline void update_vertex_buffer(gl_vertex_buffer& buf,
    const std::vector<vec3i>& values, bool dynamic = false) {
    update_vertex_buffer(
        buf, values.size(), 3, (const int*)values.data(), dynamic);
}

/// Updates the buffer with new data.
inline void update_vertex_buffer(gl_vertex_buffer& buf,
    const std::vector<vec4i>& values, bool dynamic = false) {
    update_vertex_buffer(
        buf, values.size(), 4, (const int*)values.data(), dynamic);
}

/// Make a buffer with new data.
template <typename T>
inline gl_vertex_buffer make_vertex_buffer(
    const std::vector<T>& values, bool dynamic = false) {
    auto buf = gl_vertex_buffer();
    update_vertex_buffer(buf, values, dynamic);
    return buf;
}

/// Bind the buffer at a particular attribute location.
void bind_vertex_buffer(const gl_vertex_buffer& buf, uint vattr);

/// Unbind the buffer.
void unbind_vertex_buffer(const gl_vertex_buffer& buf, uint vattr);
/// Unbind the buffer.
void unbind_vertex_buffer(uint vattr);

/// Get buffer id.
inline uint get_vertex_buffer_id(const gl_vertex_buffer& buf) {
    return buf.bid;
}

/// Check if defined.
inline bool is_vertex_buffer_valid(const gl_vertex_buffer& buf) {
    return (bool)buf.bid;
}

/// Destroys the buffer.
void clear_vertex_buffer(gl_vertex_buffer& buf);

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// OPENGL VERTEX ELEMENTS BUFFER
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup gl_element_buffer OpenGL element array buffers
/// @{

/// OpenGL element array buffer. Members are not part of the public API.
struct gl_element_buffer {
    /// Buffer id.
    uint bid = 0;
    /// Number of elements.
    int num = 0;
    /// Number of components.
    int ncomp = 0;
};

// Updates the bufferwith new data.
void _update_element_buffer(
    gl_element_buffer& buf, int n, int nc, const int* values, bool dynamic);

/// Updates the buffer with new data.
inline void update_element_buffer(gl_element_buffer& buf, int num, int ncomp,
    const int* values, bool dynamic = false) {
    _update_element_buffer(buf, num, ncomp, values, dynamic);
}

/// Updates the bufferwith new data.
inline void update_element_buffer(gl_element_buffer& buf,
    const std::vector<int>& values, bool dynamic = false) {
    update_element_buffer(buf, values.size(), 1, values.data(), dynamic);
}

/// Updates the buffer with new data.
inline void update_element_buffer(gl_element_buffer& buf,
    const std::vector<vec2i>& values, bool dynamic = false) {
    update_element_buffer(
        buf, values.size(), 2, (const int*)values.data(), dynamic);
}

/// Updates the buffer with new data.
inline void update_element_buffer(gl_element_buffer& buf,
    const std::vector<vec3i>& values, bool dynamic = false) {
    update_element_buffer(
        buf, values.size(), 3, (const int*)values.data(), dynamic);
}

/// Updates the buffer with new data.
inline void update_element_buffer(gl_element_buffer& buf,
    const std::vector<vec4i>& values, bool dynamic = false) {
    update_element_buffer(
        buf, values.size(), 4, (const int*)values.data(), dynamic);
}

/// Make a buffer with new data.
template <typename T>
inline gl_element_buffer make_element_buffer(
    const std::vector<T>& values, bool dynamic = false) {
    auto buf = gl_element_buffer();
    update_element_buffer(buf, values, dynamic);
    return buf;
}

/// Draws elements.
void draw_elems(const gl_element_buffer& buf);

/// Get id
inline uint get_element_buffer_id(const gl_element_buffer& buf) {
    return buf.bid;
}

/// Check if defined
inline bool is_element_buffer_valid(const gl_element_buffer& buf) {
    return (bool)buf.bid;
}

/// Destroys the buffer
void clear_element_buffer(gl_element_buffer& buf);

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// OPENGL PROGRAM FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup gl_program OpenGL programs
/// @{

/// OpenGL program. Members are not part of the public API.
struct gl_program {
    /// Program id.
    uint pid = 0;
    /// Vertex shader id.
    uint vid = 0;
    /// Fragment shader id.
    uint fid = 0;
    /// Vertex array object.
    uint vao = 0;
};

/// Creates an OpenGL program from vertex and fragment code.
gl_program make_program(const std::string& vertex, const std::string& fragment);

/// Destroys the program.
void clear_program(gl_program& prog);

/// Get uniform location.
int get_program_uniform_location(
    const gl_program& prog, const std::string& name);
/// Get attribute location.
int get_program_attrib_location(
    const gl_program& prog, const std::string& name);

/// Get the names of all uniforms.
std::vector<std::pair<std::string, int>> get_program_uniforms_names(
    const gl_program& prog);
/// Get the names of all attributes.
std::vector<std::pair<std::string, int>> get_program_attributes_names(
    const gl_program& prog);

/// Set uniform integer values `val` for program `prog` and variable `pos`.
/// The values have `nc` number of components (1-4) and `count` elements.
bool set_program_uniform(
    gl_program& prog, int pos, const int* val, int ncomp, int count);

/// Set uniform float values `val` for program `prog` and variable `pos`.
/// The values have `nc` number of components (1-4) and `count` elements.
bool set_program_uniform(
    gl_program& prog, int pos, const float* val, int ncomp, int count);

/// Set uniform value.
inline bool set_program_uniform(gl_program& prog, int var, bool val) {
    auto vali = (int)val;
    return set_program_uniform(prog, var, &vali, 1, 1);
}
/// Set uniform value.
inline bool set_program_uniform(gl_program& prog, int var, int val) {
    return set_program_uniform(prog, var, &val, 1, 1);
}
/// Set uniform value.
inline bool set_program_uniform(gl_program& prog, int var, float val) {
    return set_program_uniform(prog, var, &val, 1, 1);
}
/// Set uniform value.
template <typename T, int N>
inline bool set_program_uniform(
    gl_program& prog, int var, const vec<T, N>& val) {
    return set_program_uniform(prog, var, data(val), N, 1);
}
/// Set uniform value.
template <typename T>
inline bool set_program_uniform(
    gl_program& prog, int var, const mat<T, 4>& val) {
    return set_program_uniform(prog, var, (T*)data(val), 16, 1);
}
/// Set uniform value.
template <typename T>
inline bool set_program_uniform(
    gl_program& prog, int var, const frame<T, 3>& val) {
    return set_program_uniform(prog, var, (T*)data(val), 12, 1);
}

/// Set uniform array.
inline bool set_program_uniform(
    gl_program& prog, int var, const int* val, int num) {
    return set_program_uniform(prog, var, val, 1, num);
}
/// Set uniform array.
inline bool set_program_uniform(
    gl_program& prog, int var, const float* val, int num) {
    return set_program_uniform(prog, var, val, 1, num);
}
/// Set uniform array.
template <typename T, int N>
inline bool set_program_uniform(
    gl_program& prog, int var, const vec<T, N>* val, int num) {
    return set_program_uniform(prog, var, (T*)val, N, num);
}
/// Set uniform array.
template <typename T>
inline bool set_program_uniform(
    gl_program& prog, int var, const mat<T, 4>* val, int num) {
    return set_program_uniform(prog, var, (T*)val, 16, num);
}
/// Set uniform array.
template <typename T>
inline bool set_program_uniform(
    gl_program& prog, int var, const frame<T, 3>* val, int num) {
    return set_program_uniform(prog, var, (T*)val, 12, num);
}

/// Set uniform value for names variable.
template <typename T>
inline bool set_program_uniform(
    gl_program& prog, const std::string& var, const T& val) {
    auto loc = get_program_uniform_location(prog, var);
    if (loc < 0) return false;
    return set_program_uniform(prog, loc, val);
}
/// Set uniform array for names variable.
template <typename T>
inline bool set_program_uniform(
    gl_program& prog, const std::string& var, const T* val, int num) {
    auto loc = get_program_uniform_location(prog, var);
    if (loc < 0) return false;
    return set_program_uniform(loc, val, num);
}

/// Set uniform texture.
bool set_program_uniform_texture(
    gl_program& prog, int pos, const gl_texture_info& tinfo, uint tunit);
/// Set uniform texture with an additionasl texture enable flags.
inline bool set_program_uniform_texture(gl_program& prog, int var, int varon,
    const gl_texture_info& tinfo, uint tunit) {
    if (!set_program_uniform_texture(prog, var, tinfo, tunit)) return false;
    if (!set_program_uniform(prog, varon, is_texture_valid(tinfo.txt)))
        return false;
    return true;
}

/// Set uniform texture.
inline bool set_program_uniform_texture(gl_program& prog,
    const std::string& var, const gl_texture_info& tinfo, uint tunit) {
    auto loc = get_program_uniform_location(prog, var);
    if (loc < 0) return false;
    return set_program_uniform_texture(prog, loc, tinfo, tunit);
}
/// Set uniform texture with an additionasl texture enable flags.
inline bool set_program_uniform_texture(gl_program& prog,
    const std::string& var, const std::string& varon,
    const gl_texture_info& tinfo, uint tunit) {
    auto loc = get_program_uniform_location(prog, var);
    if (loc < 0) return false;
    auto locon = get_program_uniform_location(prog, varon);
    if (locon < 0) return false;
    return set_program_uniform_texture(prog, loc, locon, tinfo, tunit);
}

/// Sets a constant `value` of `nc` components for the vertex attribute at
/// `pos` location.
bool set_program_vertattr(
    gl_program& prog, int pos, const float* value, int nc);
/// Sets a constant `value` of `nc` components for the vertex attribute at
/// `pos` location.
bool set_program_vertattr(gl_program& prog, int pos, const int* value, int nc);

/// Binds a buffer to a vertex attribute.
bool set_program_vertattr(
    gl_program& prog, const std::string& var, const gl_vertex_buffer& buf);

/// Binds a buffer to a vertex attribute, or a constant value if the buffer is
/// empty.
bool set_program_vertattr(gl_program& prog, int pos,
    const gl_vertex_buffer& buf, int nc, const float* def);

/// Binds a buffer or constant to a vertex attribute.
template <typename T, int N>
inline bool set_program_vertattr(gl_program& prog, int var,
    const gl_vertex_buffer& buf, const vec<T, N>& def) {
    return set_program_vertattr(prog, var, buf, N, data(def));
}
/// Binds a buffer or constant to a vertex attribute.
template <typename T, int N>
inline bool set_program_vertattr(gl_program& prog, const std::string& var,
    const gl_vertex_buffer& buf, const vec<T, N>& def) {
    auto loc = get_program_attrib_location(prog, var);
    if (loc < 0) return false;
    return set_program_vertattr(prog, loc, buf, def);
}

/// Check whether the program is valid.
inline bool is_program_valid(const gl_program& prog) { return (bool)prog.pid; }

/// Binds a program.
void bind_program(const gl_program& prog);
/// Unbind a program.
void unbind_program(const gl_program& prog);

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// OPENGL SCENE SHADER FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup gl_shape OpenGL scene shader support
/// @{

/// Vertex buffers for scene drawing. Members are not part of the public API.
struct gl_shape {
    gl_vertex_buffer pos = {};         // position
    gl_vertex_buffer norm = {};        // normals
    gl_vertex_buffer texcoord = {};    // texcoord
    gl_vertex_buffer texcoord1 = {};   // texcoord (sond version)
    gl_vertex_buffer color = {};       // color
    gl_vertex_buffer tangsp = {};      // tangent space
    gl_element_buffer points = {};     // point elements
    gl_element_buffer lines = {};      // line elements
    gl_element_buffer triangles = {};  // triangle elements
    gl_element_buffer quads = {};      // quad elements (as 2 triangles)
    gl_element_buffer beziers = {};    // bezier elements (as 3 lines)
    gl_element_buffer edges = {};      // edge elements
};

/// Clear shape.
void clear_shape(gl_shape& shp);

/// Initialize gl lights.
gl_lights make_gl_lights(const scene* scn);

/// Clear scene textures on the GPU.
void clear_textures(std::unordered_map<texture*, gl_texture>& textures);

/// Clear scene shapes on the GPU.
void clear_shapes(std::unordered_map<shape*, gl_shape>& shapes);

/// Update scene textures on the GPU.
void update_textures(const scene* scn,
    std::unordered_map<texture*, gl_texture>& textures,
    const std::unordered_set<texture*>& refresh = {}, bool clear = false);

/// Update scene shapes on the GPU.
void update_shapes(const scene* scn,
    std::unordered_map<shape*, gl_shape>& shapes,
    const std::unordered_set<shape*>& refresh = {},
    const std::unordered_set<shape_group*>& refreshg = {}, bool clear = false);

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// OPENGL IMAGE SHADER FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup gl_stdimage_program OpenGL image shader
/// @{

/// A shader for displaying images.  Members are not part of the public API.
struct gl_stdimage_program {
    /// Program.
    gl_program prog = {};
    /// Vertex array.
    gl_vertex_buffer vbo = {};
    // Element array.
    gl_element_buffer ebo = {};
};

/// Initialize a stdimage program.
gl_stdimage_program make_stdimage_program();

/// Draws an image texture the stdimage program.
void draw_image(gl_stdimage_program& prog, const gl_texture& txt,
    const vec2i& win_size, const vec2f& offset, float zoom, float exposure,
    float gamma, bool filmic);

/// Draws an image texture the stdimage program.
inline void draw_image(gl_stdimage_program& prog, const gl_texture& txt,
    const vec2i& win_size, const vec2f& offset, float zoom) {
    draw_image(prog, txt, win_size, offset, zoom, 0, 1, false);
}

/// Params for stdimage drawing.
struct gl_stdimage_params {
    /// Window size.
    vec2i win_size = {0, 0};
    /// Image offset.
    vec2f offset = {0, 0};
    /// Image zoom.
    float zoom = 1;
    /// Tonemap exposure.
    float exposure = 1;
    /// Tonemap gamma.
    float gamma = 2.2f;
    /// Tonemap filmic.
    bool filmic = false;
    /// Image background.
    vec4f background = zero4f;
};

/// Draws an image texture the stdimage program.
inline void draw_image(gl_stdimage_program& prog, const gl_texture& txt,
    const gl_stdimage_params& params, bool clear_background = true) {
    if (clear_background) gl_clear_buffers(params.background);
    draw_image(prog, txt, params.win_size, params.offset, params.zoom,
        params.exposure, params.gamma, params.filmic);
}

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// OPENGL STANDARD SURFACE SHADER FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup gl_stdsurface_program OpenGL surface shader
/// @{

/// Program to shade surfaces with a physically-based standard shader based on
/// Phong/GGX. Members are not part of public API.
struct gl_stdsurface_program {
    /// Program.
    gl_program prog = {};
};

/// Initialize a stdsurface shader.
gl_stdsurface_program make_stdsurface_program();

/// Check if the program is valid.
inline bool is_program_valid(const gl_stdsurface_program& prog) {
    return is_program_valid(prog.prog);
}

/// Starts a frame by setting exposure/gamma values, camera transforms and
/// projection. Sets also whether to use full shading or a quick eye light
/// preview.
void begin_stdsurface_frame(gl_stdsurface_program& prog, bool shade_eyelight,
    float tonemap_exposure, float tonemap_gamma, bool tonemap_filmic,
    const mat4f& camera_xform, const mat4f& camera_xform_inv,
    const mat4f& camera_proj);

/// Ends a frame.
void end_stdsurface_frame(gl_stdsurface_program& prog);

/// Set shading lights and ambient.
void set_stdsurface_lights(
    gl_stdsurface_program& prog, const vec3f& amb, const gl_lights& lights);

/// Begins drawing a shape with transform `xform`.
void begin_stdsurface_shape(
    gl_stdsurface_program& prog, const mat4f& xform, float normal_offset = 0);

/// End shade drawing.
void end_stdsurface_shape(gl_stdsurface_program& prog);

/// Sets normal offset.
void set_stdsurface_normaloffset(
    gl_stdsurface_program& prog, float normal_offset);

/// Set the object as highlighted.
void set_stdsurface_highlight(
    gl_stdsurface_program& prog, const vec4f& highlight);

/// Set material values with emission `ke`, diffuse `kd`, specular `ks` and
/// specular roughness `rs`, opacity `op`. Indicates textures ids with the
/// correspoinding `XXX_txt` variables. Sets also normal and occlusion
/// maps. Works for points/lines/triangles indicated by `etype`, (diffuse for
/// points, Kajiya-Kay for lines, GGX/Phong for triangles). Material `type`
/// matches the scene material type.
void set_stdsurface_material(gl_stdsurface_program& prog, material_type type,
    gl_elem_type etype, const vec3f& ke, const vec3f& kd, const vec3f& ks,
    float rs, float op, const gl_texture_info& ke_txt,
    const gl_texture_info& kd_txt, const gl_texture_info& ks_txt,
    const gl_texture_info& rs_txt, const gl_texture_info& norm_txt,
    const gl_texture_info& occ_txt, bool use_phong, bool double_sided,
    bool alpha_cutout);

/// Set constant material with emission `ke` and opacity `op`.
void set_stdsurface_constmaterial(
    gl_stdsurface_program& prog, const vec3f& ke, float op);

/// Set vertex data with buffers for position pos, normals norm, texture
/// coordinates texcoord, per-vertex color color and tangent space tangsp.
void set_stdsurface_vert(gl_stdsurface_program& prog,
    const gl_vertex_buffer& pos, const gl_vertex_buffer& norm,
    const gl_vertex_buffer& texcoord, const gl_vertex_buffer& color,
    const gl_vertex_buffer& tangsp);

/// Set vertex data with buffers for skinning.
void set_stdsurface_vert_skinning(gl_stdsurface_program& prog,
    const gl_vertex_buffer& weights, const gl_vertex_buffer& joints,
    int nxforms, const mat4f* xforms);

/// Set vertex data with buffers for skinning.
void set_stdsurface_vert_gltf_skinning(gl_stdsurface_program& prog,
    const gl_vertex_buffer& weights, const gl_vertex_buffer& joints,
    int nxforms, const mat4f* xforms);

/// Disables vertex skinning.
void set_stdsurface_vert_skinning_off(gl_stdsurface_program& prog);

/// Params for stdsurface drawing.
struct gl_stdsurface_params {
    /// Image width.
    int width = 360;
    /// Image height.
    int height = 360;
    /// Image exposure.
    float exposure = 0;
    /// Image gamma.
    float gamma = 2.2f;
    /// Image filmic tonemapping.
    bool filmic = false;
    /// Draw as wireframe.
    bool wireframe = false;
    /// Draw with overlaid edges
    bool edges = false;
    /// Offset for edges.
    float edge_offset = 0.01f;
    /// Draw with an alpha cutout for binary transparency.
    bool cutout = false;
    /// Camera light mode.
    bool camera_lights = false;
    /// Window background.
    vec4f background = {0, 0, 0, 0};
    /// Ambient illumination.
    vec3f ambient = {0, 0, 0};
    /// Highlighted object.
    void* highlighted = nullptr;
    /// Highlight color.
    vec3f highlight_color = {1, 1, 0};
    /// Edge color.
    vec3f edge_color = {0, 0, 0};
    /// Cull back face.
    bool cull_backface = true;
};

/// Draw scene with stdsurface program.
void draw_stdsurface_scene(const scene* scn, const camera* cam,
    gl_stdsurface_program& prog, std::unordered_map<shape*, gl_shape>& shapes,
    std::unordered_map<texture*, gl_texture>& textures, const gl_lights& lights,
    const gl_stdsurface_params& params);

/// @}

}  // namespace ygl

// Forward declaration
struct GLFWwindow;

// -----------------------------------------------------------------------------
// OPENGL WINDOWS
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup gl_window OpenGL window
/// @{

// Forward declaration
struct gl_window;

/// Text callback.
typedef void (*gl_text_callback)(gl_window*, unsigned int);
/// Mouse callback.
typedef void (*gl_mouse_callback)(gl_window*, int button, bool press, int mods);
/// Window refresh callback.
typedef void (*gl_refresh_callback)(gl_window*);

/// OpenGL window. Members are not part of the public API.
struct gl_window {
    GLFWwindow* gwin = nullptr;
    void* user_pointer = nullptr;
    bool widget_enabled = false;
    gl_text_callback text_cb = nullptr;
    gl_mouse_callback mouse_cb = nullptr;
    gl_refresh_callback refresh_cb = nullptr;
};

/// Initialize a window.
gl_window* make_window(int width, int height, const std::string& title,
    void* user_pointer = nullptr);

/// Set window callbacks.
void set_window_callbacks(gl_window* win, gl_text_callback text_cb,
    gl_mouse_callback mouse_cb, gl_refresh_callback refresh_cb);

/// Clear window.
void clear_window(gl_window* win);

/// Gets the user poiner.
inline void* get_user_pointer(gl_window* win) { return win->user_pointer; }

/// Set window title.
void set_window_title(gl_window* win, const std::string& title);

/// Wait events
void wait_events(gl_window* win);
/// Poll events
void poll_events(gl_window* win);
/// Swap buffers
void swap_buffers(gl_window* win);

/// Should close
bool should_close(gl_window* win);

/// Window size
vec2i get_window_size(gl_window* win);
/// Framebuffer size
vec2i get_framebuffer_size(gl_window* win);

/// Mouse button
int get_mouse_button(gl_window* win);
/// Mouse position
vec2i get_mouse_pos(gl_window* win);
/// Mouse position
vec2f get_mouse_posf(gl_window* win);
/// Check if a key is pressed (not all keys are supported)
bool get_key(gl_window* win, int key);

/// Read pixels
std::vector<vec4b> get_screenshot(
    gl_window* win, vec2i& wh, bool flipy = true, bool back = false);
/// Save a screenshot to disk
inline void save_screenshot(gl_window* win, const std::string& imfilename) {
    auto wh = vec2i{0, 0};
    auto pixels = get_screenshot(win, wh);
    save_image(imfilename, wh.x, wh.y, 4, (unsigned char*)pixels.data());
}

/// Handle camera navigation.
bool handle_camera_navigation(gl_window* win, camera* cam, bool navigation_fps);

/// @}

}  // namespace ygl

// -----------------------------------------------------------------------------
// OPENGL WIDGETS
// -----------------------------------------------------------------------------
namespace ygl {

/// @defgroup gl_widgets OpenGL widgets
/// @{

/// Initialize widgets.
void init_widgets(
    gl_window* win, bool light_style = false, bool extra_font = true);

/// Begin draw widgets.
bool begin_widgets(gl_window* win, const std::string& title);
/// End draw widgets.
void end_widgets(gl_window* win);

/// Whether widgets are active.
bool get_widget_active(gl_window* win);

/// Horizontal separator.
void draw_separator_widget(gl_window* win);

/// Indent widget.
void draw_indent_widget_begin(gl_window* win);
/// Indent widget.
void draw_indent_widget_end(gl_window* win);
/// Continue line with next widget.
void draw_continue_widget(gl_window* win);

/// Label widget.
void draw_label_widget(
    gl_window* win, const std::string& lbl, const std::string& msg);

/// Label widget.
template <typename... Args>
inline void draw_label_widget(gl_window* win, const std::string& lbl,
    const std::string& fmt, const Args&... args) {
    auto msg = format(fmt, args...);
    draw_label_widget(win, lbl, msg);
}
/// Label widget.
template <typename T>
inline void draw_label_widget(
    gl_window* win, const std::string& lbl, const T& val) {
    auto sst = std::stringstream();
    sst << val;
    return draw_label_widget(win, lbl, sst.str());
}
/// Label widget.
template <typename T>
inline void draw_label_widget(gl_window* win, const std::string& lbl,
    const std::vector<T>& vals, bool skip_empty = false) {
    if (skip_empty && vals.empty()) return;
    draw_label_widget(win, lbl, (int)vals.size());
};

/// Value widget
bool draw_value_widget(gl_window* win, const std::string& lbl, bool& val);
/// Value widget.
bool draw_value_widget(
    gl_window* win, const std::string& lbl, std::string& str);
/// Value widget.
bool draw_value_widget(gl_window* win, const std::string& lbl, int* val,
    int ncomp, int min = 0, int max = 1, int incr = 1);
/// Value widget.
bool draw_value_widget(gl_window* win, const std::string& lbl, float* val,
    int ncomp, float min = 0, float max = 1, float incr = 1);
/// Value widget.
inline bool draw_value_widget(gl_window* win, const std::string& lbl, int& val,
    int min = 0, int max = 1, int incr = 1) {
    return draw_value_widget(win, lbl, &val, 1, min, max, incr);
}
/// Value widget.
inline bool draw_value_widget(gl_window* win, const std::string& lbl,
    float& val, float min = 0, float max = 1, float incr = 1) {
    return draw_value_widget(win, lbl, &val, 1, min, max, incr);
}

/// Value widget.
template <int N>
inline bool draw_value_widget(gl_window* win, const std::string& lbl,
    vec<int, N>& val, int min = 0, int max = 1, int incr = 1) {
    return draw_value_widget(win, lbl, data(val), N, min, max, incr);
}
/// Value widget.
template <int N>
inline bool draw_value_widget(gl_window* win, const std::string& lbl,
    vec<float, N>& val, float min = 0, float max = 1, float incr = 0.01f) {
    return draw_value_widget(win, lbl, data(val), N, min, max, incr);
}
/// Value widget.
inline bool draw_value_widget(gl_window* win, const std::string& lbl,
    mat<float, 4>& val, float min = 0, float max = 1, float incr = 0.01f) {
    auto modx = draw_value_widget(win, lbl + ".x", val.x, min, max, incr);
    auto mody = draw_value_widget(win, lbl + ".y", val.y, min, max, incr);
    auto modz = draw_value_widget(win, lbl + ".z", val.z, min, max, incr);
    auto modw = draw_value_widget(win, lbl + ".w", val.w, min, max, incr);
    return modx || mody || modz || modw;
}
/// Value widget.
inline bool draw_value_widget(gl_window* win, const std::string& lbl,
    frame<float, 3>& val, float min = -10, float max = 10, float incr = 0.01f) {
    auto modx = draw_value_widget(win, lbl + ".x", val.x, -1, 1, 0.01f);
    auto mody = draw_value_widget(win, lbl + ".y", val.y, -1, 1, 0.01f);
    auto modz = draw_value_widget(win, lbl + ".z", val.z, -1, 1, 0.01f);
    auto modo = draw_value_widget(win, lbl + ".o", val.o, min, max, incr);
    // TODO: orthonormalize
    return modx || mody || modz || modo;
}
/// Value widget.
template <typename T, int N>
inline bool draw_value_widget(gl_window* win, const std::string& lbl,
    quat<T, N>& val, float min = -1, float max = 1, float incr = 0.01f) {
    auto mod = draw_value_widget(win, lbl, *(vec<T, N>*)&val, min, max, incr);
    if (mod) val = normalize(val);
    return mod;
}

/// Color widget.
bool draw_color_widget(gl_window* win, const std::string& lbl, vec4f& val);
/// Color widget.
bool draw_color_widget(gl_window* win, const std::string& lbl, vec4b& val);
/// Color widget.
bool draw_color_widget(gl_window* win, const std::string& lbl, vec3f& val);

/// Combo widget.
bool draw_combo_widget_begin(
    gl_window* win, const std::string& lbl, const std::string& label);
/// Combo widget.
bool draw_combo_widget_item(
    gl_window* win, const std::string& label, int idx, bool selected);
/// Combo widget.
void draw_combo_widget_end(gl_window* win);
/// Combo widget.
template <typename T>
bool draw_combo_widget_item(
    gl_window* win, const std::string& label, int idx, T& val, const T& item) {
    auto selected = draw_combo_widget_item(win, label, idx, val == item);
    if (selected) val = item;
    return selected;
}

/// Combo widget.
template <typename T, typename T1>
inline bool draw_value_widget(gl_window* win, const std::string& lbl, T& val,
    const std::vector<T1>& vals, const std::function<T(const T1&)>& value_func,
    const std::function<std::string(const T1&)>& label_func) {
    auto label = std::string();
    for (auto& v : vals)
        if (value_func(v) == val) label = label_func(v);
    if (!draw_combo_widget_begin(win, lbl, label)) return false;
    auto changed = false;
    for (auto i = 0; i < vals.size(); i++) {
        auto selected = val == value_func(vals[i]);
        if (draw_combo_widget_item(win, label_func(vals[i]), i, selected)) {
            val = value_func(vals[i]);
            changed = true;
        }
    }
    draw_combo_widget_end(win);
    return changed;
}
/// Combo widget.
inline bool draw_value_widget(gl_window* win, const std::string& lbl,
    std::string& val, const std::vector<std::string>& labels) {
    if (!draw_combo_widget_begin(win, lbl, val)) return false;
    auto old_val = val;
    for (auto i = 0; i < labels.size(); i++) {
        draw_combo_widget_item(win, labels[i], i, val, labels[i]);
    }
    draw_combo_widget_end(win);
    return val != old_val;
}
/// Combo widget.
template <typename T>
inline bool draw_value_widget(gl_window* win, const std::string& lbl, T& val,
    const std::vector<std::pair<std::string, T>>& labels) {
    auto label = std::string();
    for (auto& kv : labels)
        if (kv.second == val) label = kv.first;
    if (!draw_combo_widget_begin(win, lbl, label)) return false;
    auto old_val = val;
    for (auto i = 0; i < labels.size(); i++) {
        draw_combo_widget_item(win, labels[i].first, i, val, labels[i].second);
    }
    draw_combo_widget_end(win);
    return val != old_val;
}

/// Combo widget
template <typename T>
inline bool draw_value_widget(gl_window* win, const std::string& lbl, T*& val,
    const std::vector<T*>& vals, bool extra = true, T* extra_val = nullptr) {
    if (!draw_combo_widget_begin(win, lbl, (val) ? val->name : "<none>"))
        return false;
    auto old_val = val;
    if (extra)
        draw_combo_widget_item(
            win, (extra_val) ? extra_val->name : "<none>", -1, val, extra_val);
    for (auto i = 0; i < vals.size(); i++) {
        draw_combo_widget_item(win, vals[i]->name, i, val, vals[i]);
    }
    draw_combo_widget_end(win);
    return val != old_val;
}

/// Button widget.
bool draw_button_widget(gl_window* win, const std::string& lbl);

/// Collapsible header widget.
bool draw_header_widget(gl_window* win, const std::string& lbl);

/// Start tree widget.
bool draw_tree_widget_begin(gl_window* win, const std::string& lbl);
/// End tree widget.
void draw_tree_widget_end(gl_window* win);
/// Start selectable tree node widget.
bool draw_tree_widget_begin(
    gl_window* win, const std::string& lbl, void*& selection, void* content);
/// Start selectable tree node widget.
bool draw_tree_widget_begin(gl_window* win, const std::string& lbl,
    void*& selection, void* content, const vec4f& col);
/// End selectable tree node widget.
void draw_tree_widget_end(gl_window* win, void* content);
/// Selectable tree leaf nodewidget.
void draw_tree_widget_leaf(
    gl_window* win, const std::string& lbl, void*& selection, void* content);
/// Selectable tree leaf node widget.
void draw_tree_widget_leaf(gl_window* win, const std::string& lbl,
    void*& selection, void* content, const vec4f& col);
/// Text color widget.
void draw_tree_widget_color_begin(gl_window* win, const vec4f& color);
/// Text color widget.
void draw_tree_widget_color_end(gl_window* win);

/// Image widget.
void draw_image_widget(
    gl_window* win, int tid, const vec2i& size, const vec2i& imsize);
/// Image widget.
void draw_image_widget(gl_window* win, gl_texture& txt, const vec2i& size);

/// Scroll region widget.
void draw_scroll_widget_begin(
    gl_window* win, const std::string& lbl, int height, bool border);
/// Scroll region widget.
void draw_scroll_widget_end(gl_window* win);
/// Scroll region widget.
void draw_scroll_widget_here(gl_window* win);

/// Group ids widget.
void draw_groupid_widget_begin(gl_window* win, int gid);
/// Group ids widget.
void draw_groupid_widget_begin(gl_window* win, void* gid);
/// Group ids widget.
void draw_groupid_widget_begin(gl_window* win, const char* gid);
/// Group ids widget.
void draw_groupid_widget_end(gl_window* win);

/// Tonemapping widgets.
inline void draw_tonemap_widgets(gl_window* win, const std::string& lbl,
    float& exposure, float& gamma, bool& filmic) {
    draw_value_widget(win, lbl + "exposure", exposure, -20, 20, 1);
    draw_value_widget(win, lbl + "gamma", gamma, 0.1, 5, 0.1);
    draw_value_widget(win, lbl + "filmic", filmic);
}

/// Image view widgets.
inline void draw_imageview_widgets(gl_window* win, const std::string& lbl,
    gl_stdimage_params& params, bool show_tonemap = true) {
    draw_value_widget(win, lbl + "offset", params.offset);
    draw_value_widget(win, lbl + "zoom", params.zoom);
    draw_color_widget(win, lbl + "background", params.background);
    if (show_tonemap) {
        draw_tonemap_widgets(
            win, lbl, params.exposure, params.gamma, params.filmic);
    }
}

/// Image inspection widgets.
inline void draw_imageinspect_widgets(gl_window* win, const std::string& lbl,
    const image4f& hdr, const image4b& ldr, const vec2f& mouse_pos,
    const gl_stdimage_params& params) {
    auto xy = (mouse_pos - params.offset) / params.zoom;
    auto i = (int)round(xy.x), j = (int)round(xy.y);
    auto v4f = zero4f;
    auto v4b = zero4b;
    if (!hdr.empty()) {
        auto w = hdr.width(), h = hdr.height();
        if (i >= 0 && i < w && j >= 0 && j < h) {
            v4f = hdr.at(i, j);
            v4b = linear_to_srgb(hdr.at(i, j));
        }
    }
    if (!ldr.empty()) {
        auto w = ldr.width(), h = ldr.height();
        if (i >= 0 && i < w && j >= 0 && j < h) {
            v4f = srgb_to_linear(ldr.at(i, j));
            v4b = ldr.at(i, j);
        }
    }
    draw_label_widget(win, lbl + "mouse pos", vec2i{i, j});
    draw_label_widget(win, lbl + "hdr val", v4f);
    draw_label_widget(win, lbl + "ldr val", v4b);
}

/// Draws a widget that can selected the camera.
inline bool draw_camera_widget(gl_window* win, const std::string& lbl,
    camera*& cam, scene* scn, camera* view) {
    return draw_value_widget(win, lbl, cam, scn->cameras, true, view);
}

/// Draws widgets for a whole scene. Used for quickly making demos.
bool draw_scene_widgets(gl_window* win, const std::string& lbl, scene* scn,
    void*& selection, const std::unordered_map<texture*, gl_texture>& gl_txt,
    test_scene_params* test_scn = nullptr);

/// @}

}  // namespace ygl

#endif

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
    auto stream = std::stringstream();
    stream << "  " << name;
    if (!flag.empty()) stream << "/" << flag;
    if (!metavar.empty()) stream << " " << metavar;
    while (stream.str().length() < 32) stream << " ";
    stream << help << " ";
    if (!req) stream << "[" << def << "]";
    if (req) stream << "(required)";
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
inline void _add_usage(cmdline_parser& parser, const std::string& name,
    const std::string& flag, bool opt, bool flag_opt, const std::string& help,
    const T& def, bool req, const std::vector<T>& choices) {
    auto stream = std::stringstream();
    stream << def;
    _add_usage_str(parser, name, flag, opt, (flag_opt) ? "" : "<val>", help,
        stream.str(), req, choices);
}

// cmdline implementation
template <typename T>
inline void _add_usage(cmdline_parser& parser, const std::string& name,
    const std::string& flag, bool opt, const std::string& help,
    const std::vector<T>& def, bool req, const std::vector<T>& choices) {
    auto stream = std::stringstream();
    auto first = true;
    for (auto&& v : def) {
        if (!first) stream << ",";
        stream << v;
        first = false;
    }
    _add_usage_str(
        parser, name, flag, opt, "<val>*", help, stream.str(), req, choices);
}

// cmdline implementation
inline void _set_error(cmdline_parser& parser, const std::string& err) {
    if (parser._error.empty()) parser._error = err;
}

// Check unused arguments.
inline bool should_exit(cmdline_parser& parser) {
    for (auto&& v : parser._to_parse) {
        if (v[0] == '-')
            _set_error(parser, "unknown option " + v);
        else
            _set_error(parser, "unknown argument " + v);
    }
    return !parser._error.empty();
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
    _add_usage(parser, name, flag, true, true, help, def, req, {});
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

// Pase an option from the command line.
template <typename T>
inline T parse_opt(cmdline_parser& parser, const std::string& name,
    const std::string& flag, const std::string& help, const T& def, bool req,
    const std::vector<T>& choices) {
    // check names
    _check_name(parser, name, flag, true);
    // update usage
    _add_usage(parser, name, flag, true, false, help, def, req, choices);
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
    auto stream = std::stringstream(arg);
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

// Parse an enum option from the command line.
template <typename T>
inline T parse_opt(cmdline_parser& parser, const std::string& name,
    const std::string& flag, const std::string& help,
    const std::vector<std::pair<std::string, T>>& key_values, const T& def,
    bool req, const std::vector<T>& choices) {
    auto keys = std::vector<std::string>{};
    auto key_def = std::string();
    for (auto&& kv : key_values) {
        keys.push_back(kv.first);
        if (kv.second == def) key_def = kv.first;
    }
    auto key =
        parse_opt<std::string>(parser, name, flag, help, key_def, req, keys);
    if (!parser._error.empty()) return def;
    auto val = def;
    for (auto&& kv : key_values) {
        if (kv.first == key) val = kv.second;
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
    _add_usage(parser, name, "", false, false, help, def, req, choices);
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
    auto stream = std::stringstream(arg);
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

// Parse all remaining positional argument from the command line.
template <typename T>
inline std::vector<T> parse_args(cmdline_parser& parser,
    const std::string& name, const std::string& help, const std::vector<T>& def,
    bool req, const std::vector<T>& choices) {
    // check names
    _check_name(parser, name, "", false);
    // update usage
    _add_usage(parser, name, "", false, help, def, req, choices);
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
        auto stream = std::stringstream(arg);
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

// Initialize a command line parser.
inline cmdline_parser make_parser(
    int argc, char** argv, const std::string& prog, const std::string& help) {
    auto parser = cmdline_parser();
    parser._to_parse = std::vector<std::string>(argv + 1, argv + argc);
    parser._usage_prog = (prog.empty()) ? std::string(argv[0]) : prog;
    parser._usage_help = help;
    parser._usage =
        parse_flag(parser, "--help", "-h", "prints and help message");
    return parser;
}

}  // namespace ygl

#endif
