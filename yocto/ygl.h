//
// # Yocto/GL: Tiny C++ Library for Physically-Based Rendering
//
//
// Yocto/GL is a collection of utilities for manipukating various graphics
// quantities and creating a simple path tracer, The list of utilities is
// described below.
//
//
// ## Design Considerations
//
// Yocto/GL follows a "data-driven programming model" that makes data explicit.
// Data is stored in simple structs and access with free functions or directly.
// All data is public, so we make no attempt at encapsulation.
// All objects is Yocto?GL have value semantic and we do not use pointers
// in data structure but indices. This means that everything can be trivially
// serialized and there is no need for memory management.
//
// We do this since this makes Yocto/GL easier to extend and quicker to learn,
// which a more explicit data flow that is easier whrn writing parallel code.
// Since Yocto/GL is mainly used for research and teaching,
// explicit data is both more hackable and easier to understand.
//
// The use of templates in Yocto was the reason for many refactorings, going
// from no template to heavy template use. After many changes, we settled
// on using templates following the established convention in the C++ standard
// library.
//
// We make use of exception for error reporting. This makes the code
// cleaner and more in line with the expectation of most other programming
// languages. At the same time, exception are not as easy to use in C++
// and are disabled in many libraries. For this reasons, this will likely
// change in the future.
//
// Finally, we import math symbols from the standard library rather than
// using the `std::name` pattern into the `ygl` namespace. This makes math code
// easier to read, and allows us to override come function implementation when
// desired.
//
//
// ## Small Vectors and Matrices, Frames, Bounding Boxes and Transforms
//
// We provide common operations for small vectors and matrices typically used
// in graphics. In particular, we support 1-4 dimensional vectors `vec<T, 1>`,
// `vec<T, 2>`, `vec<T, 3>`, `vec<T, 4>`. The one dimensional version is mostly
// for completeness.
//
// We support 1-4 dimensional generic matrices `mat<T, N, 1>`, `mat<T, N, 2>`,
// `mat<T, N, 3>`, `mat<T, N, 4>`, with matrix-matrix and matrix-vector
// products, transposes and inverses. Matrices are stored in column-major
// order and are accessed and constructed by column. The one dimensional version
// is for completeness only.
//
// To represent transformations, most of the library facilities prefer the use
// coordinate frames, aka rigid transforms, represented as `frame<T, 2>`,
// `frame<T, 3>`. The structure store three coordinate axes and the origin.
// This is equivalent to a rigid transform written as a column-major affine
// matrix. Transform operations are better behaved with this representation.
//
// We represent coordinate bounds with axis-aligned bounding boxes with
// `bbox<T, 1>`, `bbox<T, 2>`, `bbox<T, 3>`, `bbox<T, 4>`, with support for
// expansion operations for points and other bounding boxes. We provide
// operations to compute bounds for points, lines, triangles and quads.
//
// For both matrices and frames we support transform operations for points,
// vectors and directions (`transform_point()`, `transform_vector()`,
// `transform_direction()`). For frames we also the support inverse operations
// (`transform_xxx_inverse()`). Transform matrices and frames can be
// constructed from basic translation, rotation and scaling, e.g. with
// `translation_mat()` or `translation_frame()` respectively, etc.
// For rotation we support axis-angle and quaternions, with slerp.
//
// All basic types support printing and parsing using `print()` and `parse()`
// methods and iostream `<<` and `>>` operators.
//
//
// ## Geometry functions
//
// The library supports basic geomtry functions such as computing
// line/triangle/quad normals and areas, picking points on triangles
// and the like. In these functions triangles are parameterized with us written
// w.r.t the (p1-p0) and (p2-p0) axis respectively. Quads are internally handled
// as pairs of two triangles p0,p1,p3 and p2,p3,p1, with the u/v coordinates
// of the second triangle corrected as 1-u and 1-v to produce a quad
// parametrization where u and v go from 0 to 1. Degenerate quads with p2==p3
// represent triangles correctly, an this convention is used throught the
// library. This is equivalent to Intel's Embree.
//
//
// ## Shape functions
//
// We provide a small number of utilities for shape manipulation for index
// triangle and quad meshes, indexed line and point sets and indexed beziers.
// The utliities collected here are written to support a global illumination
// rendering and not for generic geometry processing. We support operation for
// shape smoothing, shape subdivision (including Catmull-Clark subdivs), and
// example shape creation.
//
// 1. compute line tangents, and triangle and quad areas and normals with
//    `line_tangent()`, `triamgle_normal()`, `quad_normal()` and
//    `line_length()`, `triangle_area()` and `quad_normal()`
// 2. interpolate values over primitives with `interpolate_line()`,
//    `interpolate_triangle()` and `interpolate_quad()`
// 3. evaluate Bezier curves and derivatives with `interpolate_bezier()` and
//    `interpolate_bezier_derivative()`
// 4. compute smooth normals and tangents with `compute_vertex_normals()`
//   `compute_vertex_tangents()`
// 5. compute tangent frames from texture coordinates with
//    `compute_tangent_spaces()`
// 6. compute skinning with `compute_skinning()` and
//    `compute_matrix_skinning()`
// 6. create shapes with `make_cube_shape()`, `make_sphere_shape()`,
// `make_quad_shape()`,
//    `make_cube_fvshape()`, `make_hair_shape()`, `make_suzanne_shape()`,
//    `make_lines_shape()`, `make_points_shape()`, `make_sphere_cube_shape()`,
//    `make_cube_rounded_shape()`, `make_sphere_flipcap_shape()`,
//    `make_cylinder_shape()`, `make_cylinder_rounded_shape()`,
//    `make_disk_shape()`, `make_cylinder_side_shape()`, `make_disk_quad_shape()`
// 7. merge element with `marge_lines()`, `marge_triangles()`, `marge_quads()`
// 8. shape sampling with `sample_points_element()`, `sample_lines_element()`,
//    `sample_triangles_element()`; initialize the sampling CDFs with
//    `sample_points_element_cdf()`, `sample_lines_element_cdf()`,
//    `sample_triangles_element_cdf()`
// 9.  sample a could of point over a surface with
// `sample_triangles_element_points()`
// 10. get edges and boundaries with `get_edges()`
// 11. convert quads to triangles with `convert_quads_to_triangles()`
// 12. convert face varying to vertex shared representations with
//     `convert_face_varying()`
// 13. subdivide elements by edge splits with `subdivide_lines()`,
//     `subdivide_triangles()`, `subdivide_quads()`, `subdivide_beziers()`
// 14. Catmull-Clark subdivision surface with `subdivide_catmullclark()`
//
//
// ## Random Number Generation, Noise, and Monte Carlo support
//
// This library supports many facilities helpful in writing sampling
// functions targeting path tracing and shape generations. Implementation of
// Perlin noise is include based on stb libraries.
//
// 1. Random number generation with PCG32:
//     1. initialize the random number generator with `make_rng()`
//     2. advance the random number state with `advance_rng()`
//     3. if necessary, you can reseed the rng with `seed_rng()`
//     4. generate random integers in an interval with `get_random_int()`
//     5. generate random floats and double in the [0,1) range with
//        `get_random_float()`, `get_random_vec2f()`, `get_random_vec3f()`,
//        `next_rand1d()`
// 2. Perlin noise: `perlin_noise()` to generate Perlin noise with optional
//    wrapping, with fractal variations `perlin_ridge_noise()`,
//    `perlin_fbm_noise()`, `perlin_turbulence_noise()`
// 3. Monte Carlo support: warp functions from [0,1)^k domains to domains
//    commonly used in path tracing. In particular, use
//    `sample_hemisphere_direction()`, `sample_sphere_direction()`,
//    `sample_hemisphere_direction_cosine()`,
//    `sample_hemisphere_direction_cospower()`. `sample_disk_point()`.
//    `sample_cylinder_point()`. `sample_triangle()`,
//    `sample_discrete_distribution()`. For each warp, you can compute
//     the PDF with `sample_xxx_pdf()`.
//
//
// ## Ray-Scene and Closest-Point Queries
//
// Yocto/GL provides ray-scene intersection for points, lines, triangles and
// quads accelerated by a two-level BVH data structure. Our BVH is written for
// minimal code and not maximum speed, but still gives reasonable results. We
// suggest the use of Intel's Embree as a more efficient alternative.
//
// In Yocto/Bvh, shapes are described as collections of indexed primitives
// (points/lines/triangles/quads) like the standard triangle mesh used in
// real-time graphics. A scene if represented as transformed instances of
// shapes. The internal data structure is a two-level BVH, with a BVH for each
// shape and one top-level BVH for the whole scene. This design support
// instancing for large scenes and easy BVH refitting for interactive
// applications.
//
// 1. fill the shape or instance data
// 2. build the BVH with `build_scene_bvh()`
// 3. perform ray-element intersection with `intersect_scene_bvh()`
// 4. perform point overlap queries with `overlap_scene_bvh()`
// 5. refit the BVH with `refit_scene_bvh()` after updating internal data
//
//
// ## Image Utilities
//
// Yocto/GL supports a very small set is color and image utilities including
// color utilities, example image creation, tone mapping, image resizing, and
// sunsky procedural images. Yocto/Image is written to support the need of a
// global illumination renderer, rather than the need of generic image editing.
//
// 0. load and save image with Yocto/GLIO
// 1. create images with `image<T>` data structure
// 2. resize images with `resize_image()`
// 3. tonemap images with `tonemap_filmic()` that convert from linear HDR to
//    sRGB LDR with exposure and an optional filmic curve
// 5. make various image examples with the `make_XXX_image()` functions
// 6. create procedural sun-sky images with `make_sunsky_image()`
//
//
// # Simple scene representation
//
// Yocto/GL define a simple scene data structure useful to create quick demos
// and as the repsetnation upon which the path tracer works.
//
// In Yocto scenes, shapes are represented as indexed collections of points,
// lines, triangles, quads and bezier segments. Each shape may contain
// only one element type. Shapes are organized into a scene by creating shape
// instances, each its own transform. Materials are specified like in OBJ and
// glTF and include emission, base-metallic and diffuse-specular
// parametrization, normal, occlusion and displacement mapping. Finally, the
// scene containers cameras and environment maps. Quad support in shapes is
// experimental and mostly supported for loading and saving. Lights in
// Yocto/Scene are pointers to either instances or environments. The scene
// supports an optional node hierarchy with animation modeled on the glTF model.
//
// 1. load a scene with Yocto/GLIO,
// 2. add missing data with `add_XXX()` functions
// 3. use `compute_shape_box()/compute_scene_box()` to compute element bounds
// 4. can merge scene together with `merge_scene()`
// 6. for ray-intersection and closest point queries, a BVH can be created with
//    `make_shape_bvh()/make_scene_bvh()` and refit with
//    `refit_shape_bvh()/refit_scene_bvh()`
// 7. compute interpolated values over scene elements with `evaluate_XXX()`
//    functions
//
//
// ## Physically-based Path Tracing
//
// Yocto/GL includes a tiny, but fully featured, path tracer with support for
// textured mesh area lights, GGX materials, environment mapping. The algorithm
// makes heavy use of MIS for fast convergence.
// The interface supports progressive parallel execution both synchronously,
// for CLI applications, and asynchronously for interactive viewing.
//
// Materials are represented as sums of an emission term, a diffuse term and
// a specular microfacet term (GGX or Phong), and a transmission term for
// this sheet glass.
// Lights are defined as any shape with a material emission term. Additionally
// one can also add environment maps. But even if you can, you might want to
// add a large triangle mesh with inward normals instead. The latter is more
// general (you can even more an arbitrary shape sun). For now only the first
// environment is used.
//
// 1. prepare the ray-tracing acceleration structure with `build_scene_bvh()`
// 2. prepare lights for rendering with `make_trace_lights()`
// 3. create the image buffer and random number generators `make_trace_state()`
// 4. render blocks of samples with `trace_samples()`
// 5. you can also start an asynchronous renderer with `trace_asynch_start()`
//
//
// ## Utilities foe containers
//
// Internally we use a few container utilities to make the use of STL simpler.
// These are just thin wrapper to STL idioms.
//
// 1. check for containment with `contains()`
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2018 Fabio Pellacini
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

#ifndef _YGL_H_
#define _YGL_H_

#ifndef YGL_EMBREE
#define YGL_EMBREE 1
#endif

#ifndef YGL_QUADS_AS_TRIANGLES
#define YGL_QUADS_AS_TRIANGLES 1
#endif

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <algorithm>
#include <array>
#include <atomic>
#include <cctype>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <deque>
#include <functional>
#include <limits>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// -----------------------------------------------------------------------------
// MATH CONSTANTS AND FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

using std::abs;
using std::acos;
using std::asin;
using std::atan;
using std::atan2;
using std::cos;
using std::exp;
using std::exp2;
using std::fabs;
using std::floor;
using std::fmod;
using std::get;
using std::isfinite;
using std::log;
using std::pow;
using std::round;
using std::sin;
using std::sqrt;
using std::swap;
using std::tan;

using std::array;
using std::atomic;
using std::deque;
using std::function;
using std::ignore;
using std::lock_guard;
using std::make_unique;
using std::mutex;
using std::runtime_error;
using std::string;
using std::thread;
using std::tie;
using std::tuple;
using std::unordered_map;
using std::vector;
using namespace std::string_literals;

using byte = unsigned char;
using uint = unsigned int;

template <typename T>
inline T min(const T& x, const T& y) {
    return (x < y) ? x : y;
}
template <typename T>
inline T max(const T& x, const T& y) {
    return (x > y) ? x : y;
}
template <typename T>
inline T clamp(const T& value, const T& min_, const T& max_) {
    return min(max(value, min_), max_);
}
template <typename T, typename T1>
inline T lerp(const T& a, const T& b, T1 u) {
    return a * (1 - u) + b * u;
}

template <class T>
const T mint_ = std::numeric_limits<T>::lowest();
template <class T>
const T maxt_ = std::numeric_limits<T>::max();
template <class T>
const T epst_ = std::numeric_limits<T>::epsilon();

// type limits with simpler interface
template <typename T>
constexpr inline T mint() {
    return std::numeric_limits<T>::lowest();
}
template <typename T>
constexpr inline T maxt() {
    return std::numeric_limits<T>::max();
}
template <typename T>
constexpr inline T epst() {
    return std::numeric_limits<T>::epsilon();
}
constexpr const auto maxf = maxt<float>();
constexpr const auto epsf = epst<float>();

template <class T>
constexpr const T    pi  = (T)3.14159265358979323846;
constexpr const auto pif = 3.14159265f;

}  // namespace ygl

// -----------------------------------------------------------------------------
// PRINT/PARSE UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Formats a string `fmt` with values taken from `args`. Uses `{}` as
// placeholder.
template <typename... Args>
inline string format(const string& fmt, const Args&... args);

// Converts to string.
template <typename T>
inline string to_string(const T& value);

// Prints a formatted string to stdout or file.
template <typename... Args>
inline bool print(FILE* fs, const string& fmt, const Args&... args);
template <typename... Args>
inline bool print(const string& fmt, const Args&... args) {
    return print(stdout, fmt, args...);
}

// Format duration string from nanoseconds
inline string format_duration(int64_t duration);
// Format a large integer number in human readable form
inline string format_num(uint64_t num);

// Parse a list of space separated values.
template <typename... Args>
inline bool parse(const string& str, Args&... args);

// get time in nanoseconds - useful only to compute difference of times
inline int64_t get_time() {
    return std::chrono::high_resolution_clock::now().time_since_epoch().count();
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// LOGGING UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Log info/error/fatal/trace message
template <typename... Args>
inline void log_info(const string& fmt, const Args&... args);
template <typename... Args>
inline void log_error(const string& fmt, const Args&... args);
template <typename... Args>
inline void log_fatal(const string& fmt, const Args&... args);

// Setup logging
inline void set_log_console(bool enabled);
inline void set_log_file(const string& filename, bool append = false);

// Log traces for timing and program debugging
struct log_scope;
template <typename... Args>
inline void log_trace(const string& fmt, const Args&... args);
template <typename... Args>
inline log_scope log_trace_begin(const string& fmt, const Args&... args);
template <typename... Args>
inline void log_trace_end(log_scope& scope);
template <typename... Args>
inline log_scope log_trace_scoped(const string& fmt, const Args&... args);

}  // namespace ygl

// -----------------------------------------------------------------------------
// VECTORS
// -----------------------------------------------------------------------------
namespace ygl {

// Small size vectors.
template <typename T, int N>
struct vec;

// Small size vectors.
template <typename T>
struct vec<T, 1> {
    T x = 0;
};
template <typename T>
struct vec<T, 2> {
    T x = 0;
    T y = 0;
};
template <typename T>
struct vec<T, 3> {
    T x = 0;
    T y = 0;
    T z = 0;
};
template <typename T>
struct vec<T, 4> {
    T x = 0;
    T y = 0;
    T z = 0;
    T w = 0;
};

// Type aliases.
using vec1f = vec<float, 1>;
using vec2f = vec<float, 2>;
using vec3f = vec<float, 3>;
using vec4f = vec<float, 4>;
using vec1i = vec<int, 2>;
using vec2i = vec<int, 2>;
using vec3i = vec<int, 3>;
using vec4i = vec<int, 4>;
using vec4b = vec<byte, 4>;

// Zero vector constants.
constexpr const auto zero1f = vec1f{0};
constexpr const auto zero2f = vec2f{0, 0};
constexpr const auto zero3f = vec3f{0, 0, 0};
constexpr const auto zero4f = vec4f{0, 0, 0, 0};
constexpr const auto zero1i = vec1i{0};
constexpr const auto zero2i = vec2i{0, 0};
constexpr const auto zero3i = vec3i{0, 0, 0};
constexpr const auto zero4i = vec4i{0, 0, 0, 0};
constexpr const auto zero4b = vec4b{0, 0, 0, 0};

// Access component by index.
template <typename T, int N>
constexpr inline T& at(vec<T, N>& v, int i) {
    return *(&v.x + i);
}
template <typename T, int N>
constexpr inline const T& at(const vec<T, N>& v, int i) {
    return *(&v.x + i);
}
template <int I, typename T, int N>
constexpr inline T& get(vec<T, N>& v) {
    return *(&v.x + N);
}
template <int I, typename T, int N>
constexpr inline const T& get(const vec<T, N>& v) {
    return *(&v.x + N);
}

// Access xyz component of a vec4 typically used for color operation.
template <typename T>
constexpr inline vec<T, 3>& xyz(const vec<T, 4>& v) {
    return (vec<T, 3>&)v;
}
template <typename T>
constexpr inline vec<T, 3>& xyz(vec<T, 4>& v) {
    return (vec<T, 3>&)v;
}

// Iteration and data
template <typename T, int N>
constexpr int size(const vec<T, N>& v) {
    return N;
}
template <typename T, int N>
constexpr T* begin(vec<T, N>& v) {
    return &v.x;
}
template <typename T, int N>
constexpr const T* begin(const vec<T, N>& v) {
    return &v.x;
}
template <typename T, int N>
constexpr T* end(vec<T, N>& v) {
    return &v.x + N;
}
template <typename T, int N>
constexpr const T* end(const vec<T, N>& v) {
    return &v.x + N;
}
template <typename T, int N>
constexpr T* data(vec<T, N>& v) {
    return &v.x;
}
template <typename T, int N>
constexpr const T* data(const vec<T, N>& v) {
    return &v.x;
}

// Vector comparison operations.
template <typename T>
constexpr inline bool operator==(const vec<T, 1>& a, const vec<T, 1>& b) {
    return a.x == b.x;
}
template <typename T>
constexpr inline bool operator!=(const vec<T, 1>& a, const vec<T, 1>& b) {
    return a.x != b.x;
}
template <typename T, typename T1>
constexpr inline bool operator==(const vec<T, 1>& a, T1 b) {
    return a.x == b;
}
template <typename T, typename T1>
constexpr inline bool operator!=(const vec<T, 1>& a, T1 b) {
    return a.x != b;
}
template <typename T>
constexpr inline bool operator==(const vec<T, 2>& a, const vec<T, 2>& b) {
    return a.x == b.x && a.y == b.y;
}
template <typename T>
constexpr inline bool operator!=(const vec<T, 2>& a, const vec<T, 2>& b) {
    return a.x != b.x || a.y != b.y;
}
template <typename T, typename T1>
constexpr inline bool operator==(const vec<T, 2>& a, T1 b) {
    return a.x == b && a.y == b;
}
template <typename T, typename T1>
constexpr inline bool operator!=(const vec<T, 2>& a, T1 b) {
    return a.x != b || a.y != b;
}
template <typename T>
constexpr inline bool operator==(const vec<T, 3>& a, const vec<T, 3>& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}
template <typename T>
constexpr inline bool operator!=(const vec<T, 3>& a, const vec<T, 3>& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z;
}
template <typename T, typename T1>
constexpr inline bool operator==(const vec<T, 3>& a, T1 b) {
    return a.x == b && a.y == b && a.z == b;
}
template <typename T, typename T1>
constexpr inline bool operator!=(const vec<T, 3>& a, T1 b) {
    return a.x != b || a.y != b || a.z != b;
}
template <typename T>
constexpr inline bool operator==(const vec<T, 4>& a, const vec<T, 4>& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}
template <typename T>
constexpr inline bool operator!=(const vec<T, 4>& a, const vec<T, 4>& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w;
}
template <typename T, typename T1>
constexpr inline bool operator==(const vec<T, 4>& a, T1 b) {
    return a.x == b && a.y == b && a.z == b && a.w == b;
}
template <typename T, typename T1>
constexpr inline bool operator!=(const vec<T, 4>& a, T1 b) {
    return a.x != b || a.y != b || a.z != b || a.w != b;
}

// Vector operations.
template <typename T>
constexpr inline vec<T, 1> operator-(const vec<T, 1>& a) {
    return {-a.x};
}
template <typename T>
constexpr inline vec<T, 1> operator+(const vec<T, 1>& a, const vec<T, 1>& b) {
    return {a.x + b.x};
}
template <typename T, typename T1>
constexpr inline vec<T, 1> operator+(const vec<T, 1>& a, T1 b) {
    return {a.x + b};
}
template <typename T, typename T1>
constexpr inline vec<T, 1> operator+(T1 a, const vec<T, 1>& b) {
    return {a + b.x};
}
template <typename T>
constexpr inline vec<T, 1> operator-(const vec<T, 1>& a, const vec<T, 1>& b) {
    return {a.x - b.x};
}
template <typename T, typename T1>
constexpr inline vec<T, 1> operator-(const vec<T, 1>& a, T1 b) {
    return {a.x - b};
}
template <typename T, typename T1>
constexpr inline vec<T, 1> operator-(T1 a, const vec<T, 1>& b) {
    return {a - b.x};
}
template <typename T>
constexpr inline vec<T, 1> operator*(const vec<T, 1>& a, const vec<T, 1>& b) {
    return {a.x * b.x};
}
template <typename T, typename T1>
constexpr inline vec<T, 1> operator*(const vec<T, 1>& a, T1 b) {
    return {a.x * b};
}
template <typename T, typename T1>
constexpr inline vec<T, 1> operator*(T1 a, const vec<T, 1>& b) {
    return {a * b.x};
}
template <typename T>
constexpr inline vec<T, 1> operator/(const vec<T, 1>& a, const vec<T, 1>& b) {
    return {a.x / b.x};
}
template <typename T, typename T1>
constexpr inline vec<T, 1> operator/(const vec<T, 1>& a, T1 b) {
    return {a.x / b};
}
template <typename T, typename T1>
constexpr inline vec<T, 1> operator/(T1 a, const vec<T, 1>& b) {
    return {a / b.x};
}

// Vector operations.
template <typename T>
constexpr inline vec<T, 2> operator-(const vec<T, 2>& a) {
    return {-a.x, -a.y};
}
template <typename T>
constexpr inline vec<T, 2> operator+(const vec<T, 2>& a, const vec<T, 2>& b) {
    return {a.x + b.x, a.y + b.y};
}
template <typename T, typename T1>
constexpr inline vec<T, 2> operator+(const vec<T, 2>& a, T1 b) {
    return {a.x + b, a.y + b};
}
template <typename T, typename T1>
constexpr inline vec<T, 2> operator+(T1 a, const vec<T, 2>& b) {
    return {a + b.x, a + b.y};
}
template <typename T>
constexpr inline vec<T, 2> operator-(const vec<T, 2>& a, const vec<T, 2>& b) {
    return {a.x - b.x, a.y - b.y};
}
template <typename T, typename T1>
constexpr inline vec<T, 2> operator-(const vec<T, 2>& a, T1 b) {
    return {a.x - b, a.y - b};
}
template <typename T, typename T1>
constexpr inline vec<T, 2> operator-(T1 a, const vec<T, 2>& b) {
    return {a - b.x, a - b.y};
}
template <typename T>
constexpr inline vec<T, 2> operator*(const vec<T, 2>& a, const vec<T, 2>& b) {
    return {a.x * b.x, a.y * b.y};
}
template <typename T, typename T1>
constexpr inline vec<T, 2> operator*(const vec<T, 2>& a, T1 b) {
    return {a.x * b, a.y * b};
}
template <typename T, typename T1>
constexpr inline vec<T, 2> operator*(T1 a, const vec<T, 2>& b) {
    return {a * b.x, a * b.y};
}
template <typename T>
constexpr inline vec<T, 2> operator/(const vec<T, 2>& a, const vec<T, 2>& b) {
    return {a.x / b.x, a.y / b.y};
}
template <typename T, typename T1>
constexpr inline vec<T, 2> operator/(const vec<T, 2>& a, T1 b) {
    return {a.x / b, a.y / b};
}
template <typename T, typename T1>
constexpr inline vec<T, 2> operator/(T1 a, const vec<T, 2>& b) {
    return {a / b.x, a / b.y};
}

// Vector operations.
template <typename T>
constexpr inline vec<T, 3> operator+(const vec<T, 3>& a) {
    return a;
}
template <typename T>
constexpr inline vec<T, 3> operator-(const vec<T, 3>& a) {
    return {-a.x, -a.y, -a.z};
}
template <typename T>
constexpr inline vec<T, 3> operator+(const vec<T, 3>& a, const vec<T, 3>& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
template <typename T, typename T1>
constexpr inline vec<T, 3> operator+(const vec<T, 3>& a, T1 b) {
    return {a.x + b, a.y + b, a.z + b};
}
template <typename T, typename T1>
constexpr inline vec<T, 3> operator+(T1 a, const vec<T, 3>& b) {
    return {a + b.x, a + b.y, a + b.z};
}
template <typename T>
constexpr inline vec<T, 3> operator-(const vec<T, 3>& a, const vec<T, 3>& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
template <typename T, typename T1>
constexpr inline vec<T, 3> operator-(const vec<T, 3>& a, T1 b) {
    return {a.x - b, a.y - b, a.z - b};
}
template <typename T, typename T1>
constexpr inline vec<T, 3> operator-(T1 a, const vec<T, 3>& b) {
    return {a - b.x, a - b.y, a - b.z};
}
template <typename T>
constexpr inline vec<T, 3> operator*(const vec<T, 3>& a, const vec<T, 3>& b) {
    return {a.x * b.x, a.y * b.y, a.z * b.z};
}
template <typename T, typename T1>
constexpr inline vec<T, 3> operator*(const vec<T, 3>& a, T1 b) {
    return {a.x * b, a.y * b, a.z * b};
}
template <typename T, typename T1>
constexpr inline vec<T, 3> operator*(T1 a, const vec<T, 3>& b) {
    return {a * b.x, a * b.y, a * b.z};
}
template <typename T>
constexpr inline vec<T, 3> operator/(const vec<T, 3>& a, const vec<T, 3>& b) {
    return {a.x / b.x, a.y / b.y, a.z / b.z};
}
template <typename T, typename T1>
constexpr inline vec<T, 3> operator/(const vec<T, 3>& a, T1 b) {
    return {a.x / b, a.y / b, a.z / b};
}
template <typename T, typename T1>
constexpr inline vec<T, 3> operator/(T1 a, const vec<T, 3>& b) {
    return {a / b.x, a / b.y, a / b.z};
}

// Vector operations.
template <typename T>
constexpr inline vec<T, 4> operator-(const vec<T, 4>& a) {
    return {-a.x, -a.y, -a.z, -a.w};
}
template <typename T>
constexpr inline vec<T, 4> operator+(const vec<T, 4>& a, const vec<T, 4>& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}
template <typename T, typename T1>
constexpr inline vec<T, 4> operator+(const vec<T, 4>& a, T1 b) {
    return {a.x + b, a.y + b, a.z + b, a.w + b};
}
template <typename T, typename T1>
constexpr inline vec<T, 4> operator+(T1 a, const vec<T, 4>& b) {
    return {a + b.x, a + b.y, a + b.z, a + b.w};
}
template <typename T>
constexpr inline vec<T, 4> operator-(const vec<T, 4>& a, const vec<T, 4>& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w};
}
template <typename T, typename T1>
constexpr inline vec<T, 4> operator-(const vec<T, 4>& a, T1 b) {
    return {a.x - b, a.y - b, a.z - b, a.w - b};
}
template <typename T, typename T1>
constexpr inline vec<T, 4> operator-(T1 a, const vec<T, 4>& b) {
    return {a - b.x, a - b.y, a - b.z, a - b.w};
}
template <typename T>
constexpr inline vec<T, 4> operator*(const vec<T, 4>& a, const vec<T, 4>& b) {
    return {a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w};
}
template <typename T, typename T1>
constexpr inline vec<T, 4> operator*(const vec<T, 4>& a, T1 b) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
}
template <typename T, typename T1>
constexpr inline vec<T, 4> operator*(T1 a, const vec<T, 4>& b) {
    return {a * b.x, a * b.y, a * b.z, a * b.w};
}
template <typename T>
constexpr inline vec<T, 4> operator/(const vec<T, 4>& a, const vec<T, 4>& b) {
    return {a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w};
}
template <typename T, typename T1>
constexpr inline vec<T, 4> operator/(const vec<T, 4>& a, T1 b) {
    return {a.x / b, a.y / b, a.z / b, a.w / b};
}
template <typename T, typename T1>
constexpr inline vec<T, 4> operator/(T1 a, const vec<T, 4>& b) {
    return {a / b.x, a / b.y, a / b.z, a / b.w};
}

// Vector assignments
template <typename T, int N>
constexpr inline vec<T, N>& operator+=(vec<T, N>& a, const vec<T, N>& b) {
    return a = a + b;
}
template <typename T, int N, typename T1>
constexpr inline vec<T, N>& operator+=(vec<T, N>& a, T1 b) {
    return a = a + b;
}
template <typename T, int N>
constexpr inline vec<T, N>& operator-=(vec<T, N>& a, const vec<T, N>& b) {
    return a = a - b;
}
template <typename T, int N, typename T1>
constexpr inline vec<T, N>& operator-=(vec<T, N>& a, T1 b) {
    return a = a - b;
}
template <typename T, int N>
constexpr inline vec<T, N>& operator*=(vec<T, N>& a, const vec<T, N>& b) {
    return a = a * b;
}
template <typename T, int N, typename T1>
constexpr inline vec<T, N>& operator*=(vec<T, N>& a, T1 b) {
    return a = a * b;
}
template <typename T, int N>
constexpr inline vec<T, N>& operator/=(vec<T, N>& a, const vec<T, N>& b) {
    return a = a / b;
}
template <typename T, int N, typename T1>
constexpr inline vec<T, N>& operator/=(vec<T, N>& a, T1 b) {
    return a = a / b;
}

// Vector products and lengths.
template <typename T>
constexpr inline T dot(const vec<T, 1>& a, const vec<T, 1>& b) {
    return a.x * b.x;
}
template <typename T>
constexpr inline T dot(const vec<T, 2>& a, const vec<T, 2>& b) {
    return a.x * b.x + a.y * b.y;
}
template <typename T>
constexpr inline T cross(const vec<T, 2>& a, const vec<T, 2>& b) {
    return a.x * b.y - a.y * b.x;
}
template <typename T>
constexpr inline T dot(const vec<T, 3>& a, const vec<T, 3>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
template <typename T>
constexpr inline vec<T, 3> cross(const vec<T, 3>& a, const vec<T, 3>& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
template <typename T>
constexpr inline T dot(const vec<T, 4>& a, const vec<T, 4>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

template <typename T, int N>
inline T length(const vec<T, N>& a) {
    return sqrt(dot(a, a));
}
template <typename T, int N>
constexpr inline T length_squared(const vec<T, N>& a) {
    return dot(a, a);
}
template <typename T, int N>
inline vec<T, N> normalize(const vec<T, N>& a) {
    auto l = length(a);
    return (l) ? a / l : a;
}
template <typename T, int N>
inline T distance(const vec<T, N>& a, const vec<T, N>& b) {
    return length(a - b);
}
template <typename T, int N>
inline T distance_squared(const vec<T, N>& a, const vec<T, N>& b) {
    return length_squared(a - b);
}

// Vector angles and slerps.
template <typename T>
inline T angle(const vec<T, 3>& a, const vec<T, 3>& b) {
    return acos(clamp(dot(normalize(a), normalize(b)), -1.0f, 1.0f));
}
template <typename T, typename T1>
inline vec<T, 4> slerp(const vec<T, 4>& a, const vec<T, 4>& b, T1 u) {
    // https://en.wikipedia.org/wiki/Slerp
    auto an = normalize(a), bn = normalize(b);
    auto d = dot(an, bn);
    if (d < 0) {
        bn = -bn;
        d  = -d;
    }
    if (d > 0.9995f) return normalize(an + u * (bn - an));
    auto th = acos(clamp(d, -1.0f, 1.0f));
    if (!th) return an;
    return an * (sin(th * (1 - u)) / sin(th)) + bn * (sin(th * u) / sin(th));
}

// Orthogonal vectors.
template <typename T>
inline vec<T, 3> orthogonal(const vec<T, 3>& v) {
    // http://lolengine.net/blog/2013/09/21/picking-orthogonal-vector-combing-coconuts)
    return fabs(v.x) > fabs(v.z) ? vec<T, 3>{-v.y, v.x, 0} :
                                   vec<T, 3>{0, -v.z, v.y};
}
template <typename T>
inline vec<T, 3> orthonormalize(const vec<T, 3>& a, const vec<T, 3>& b) {
    return normalize(a - b * dot(a, b));
}

// Reflected and refracted vector.
template <typename T>
inline vec<T, 3> reflect(const vec<T, 3>& direction, const vec<T, 3>& normal) {
    return -direction + 2 * dot(normal, direction) * normal;
}
template <typename T, typename T1>
inline vec<T, 3> refract(
    const vec<T, 3>& direction, const vec<T, 3>& normal, T1 eta) {
    // auto k = 1.0 - eta * eta * (1.0 - dot(n, w) * dot(n, w));
    auto k = 1 -
             eta * eta *
                 max(0.0f, 1 - dot(normal, direction) * dot(normal, direction));
    if (k < 0) return {0, 0, 0};  // tir
    return -direction * eta + (eta * dot(normal, direction) - sqrt(k)) * normal;
}

// Max element and clamp.
template <typename T, typename T1, typename T2>
constexpr inline vec<T, 2> clamp(const vec<T, 2>& value, T1 min, T2 max) {
    return {clamp(value.x, min, max), clamp(value.y, min, max)};
}
template <typename T, typename T1, typename T2>
constexpr inline vec<T, 3> clamp(const vec<T, 3>& value, T1 min, T2 max) {
    return {clamp(value.x, min, max), clamp(value.y, min, max),
        clamp(value.z, min, max)};
}
template <typename T, typename T1, typename T2>
constexpr inline vec<T, 4> clamp(const vec<T, 4>& value, T1 min, T2 max) {
    return {clamp(value.x, min, max), clamp(value.y, min, max),
        clamp(value.z, min, max), clamp(value.w, min, max)};
}
template <typename T>
constexpr inline T max(const vec<T, 2>& a) {
    return max(a.x, a.y);
}
template <typename T>
constexpr inline T max(const vec<T, 3>& a) {
    return max(max(a.x, a.y), a.z);
}
template <typename T>
constexpr inline T max(const vec<T, 4>& a) {
    return max(max(max(a.x, a.y), a.z), a.w);
}
template <typename T>
constexpr inline T min(const vec<T, 2>& a) {
    return min(a.x, a.y);
}
template <typename T>
constexpr inline T min(const vec<T, 3>& a) {
    return min(min(a.x, a.y), a.z);
}
template <typename T>
constexpr inline T min(const vec<T, 4>& a) {
    return min(min(min(a.x, a.y), a.z), a.w);
}
template <typename T>
constexpr inline T mean(const vec<T, 2>& a) {
    return (a.x + a.y) / 2;
}
template <typename T>
constexpr inline T mean(const vec<T, 3>& a) {
    return (a.x + a.y + a.z) / 3;
}
template <typename T>
constexpr inline T mean(const vec<T, 4>& a) {
    return (a.x + a.y + a.z + a.w) / 4;
}

// Quaternion operatons represented as xi + yj + zk + w
const auto identity_quat4f = vec4f{0, 0, 0, 1};
template <typename T>
constexpr inline vec<T, 4> quat_mul(const vec<T, 4>& a, float b) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
}
template <typename T>
constexpr inline vec<T, 4> quat_mul(const vec<T, 4>& a, const vec<T, 4>& b) {
    return {a.x * b.w + a.w * b.x + a.y * b.w - a.z * b.y,
        a.y * b.w + a.w * b.y + a.z * b.x - a.x * b.z,
        a.z * b.w + a.w * b.z + a.x * b.y - a.y * b.x,
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z};
}
template <typename T>
constexpr inline vec<T, 4> quat_conjugate(const vec<T, 4>& a) {
    return {-a.x, -a.y, -a.z, a.w};
}
template <typename T>
constexpr inline vec<T, 4> quat_inverse(const vec<T, 4>& a) {
    return quat_conjugate(a) / dot(a, a);
}

}  // namespace ygl

namespace std {

// Hash functor for vector for use with unordered_map
template <typename T>
struct hash<ygl::vec<T, 2>> {
    size_t operator()(const ygl::vec<T, 2>& v) const {
        auto vh = hash<T>();
        auto h  = (size_t)0;
        for (auto i = 0; i < 2; i++)
            h ^= vh((&v.x)[i]) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};
// Hash functor for vector for use with unordered_map
template <typename T>
struct hash<ygl::vec<T, 3>> {
    size_t operator()(const ygl::vec<T, 3>& v) const {
        auto vh = hash<int>();
        auto h  = (size_t)0;
        for (auto i = 0; i < 3; i++)
            h ^= vh((&v.x)[i]) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};
// Hash functor for vector for use with unordered_map
template <typename T>
struct hash<ygl::vec<T, 4>> {
    size_t operator()(const ygl::vec<T, 4>& v) const {
        auto vh = hash<int>();
        auto h  = (size_t)0;
        for (auto i = 0; i < 4; i++)
            h ^= vh((&v.x)[i]) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

}  // namespace std

// -----------------------------------------------------------------------------
// MATRICES
// -----------------------------------------------------------------------------
namespace ygl {

// Small Fixed-size square matrices stored in column major format.
template <typename T, int N, int M>
struct mat;

// Small Fixed-size square matrices stored in column major format.
template <typename T, int N>
struct mat<T, N, 1> {
    vec<T, N> x = {};
};
template <typename T, int N>
struct mat<T, N, 2> {
    vec<T, N> x = {};
    vec<T, N> y = {};
};
template <typename T, int N>
struct mat<T, N, 3> {
    vec<T, 3> x = {};
    vec<T, 3> y = {};
    vec<T, 3> z = {};
};
template <typename T, int N>
struct mat<T, N, 4> {
    vec<T, 4> x = {};
    vec<T, 4> y = {};
    vec<T, 4> z = {};
    vec<T, 4> w = {};
};

// Type aliases.
using mat1f = mat<float, 1, 1>;
using mat2f = mat<float, 2, 2>;
using mat3f = mat<float, 3, 3>;
using mat4f = mat<float, 4, 4>;

// Identity matrices constants.
const auto identity_mat1f = mat1f{{1}};
const auto identity_mat2f = mat2f{{1, 0}, {0, 1}};
const auto identity_mat3f = mat3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
const auto identity_mat4f = mat4f{
    {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};

// Matrix comparisons.
template <typename T, int N>
constexpr inline bool operator==(const mat<T, N, 1>& a, const mat<T, N, 1>& b) {
    return a.x == b.x;
}
template <typename T, int N>
constexpr inline bool operator!=(const mat<T, N, 1>& a, const mat<T, N, 1>& b) {
    return !(a == b);
}
template <typename T, int N>
constexpr inline bool operator==(const mat<T, N, 2>& a, const mat<T, N, 2>& b) {
    return a.x == b.x && a.y == b.y;
}
template <typename T, int N>
constexpr inline bool operator!=(const mat<T, N, 2>& a, const mat<T, N, 2>& b) {
    return !(a == b);
}
template <typename T, int N>
constexpr inline bool operator==(const mat<T, N, 3>& a, const mat<T, N, 3>& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}
template <typename T, int N>
constexpr inline bool operator!=(const mat<T, N, 3>& a, const mat<T, N, 3>& b) {
    return !(a == b);
}
template <typename T, int N>
constexpr inline bool operator==(const mat<T, N, 4>& a, const mat<T, N, 4>& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}
template <typename T, int N>
constexpr inline bool operator!=(const mat<T, N, 4>& a, const mat<T, N, 4>& b) {
    return !(a == b);
}

// Matrix operations.
template <typename T, int N>
constexpr inline mat<T, N, 1> operator+(
    const mat<T, N, 1>& a, const mat<T, N, 1>& b) {
    return {a.x + b.x};
}
template <typename T, int N>
constexpr inline mat<T, N, 1> operator*(const mat<T, N, 1>& a, T b) {
    return {a.x * b};
}
template <typename T, int N>
constexpr inline vec<T, 1> operator*(const mat<T, N, 1>& a, const vec<T, 1>& b) {
    return a.x * b.x;
}
template <typename T, int N>
constexpr inline vec<T, 1> operator*(const vec<T, N>& a, const mat<T, N, 1>& b) {
    return {dot(a, b.x)};
}
template <typename T, int N, int M>
constexpr inline mat<T, N, 1> operator*(
    const mat<T, N, M>& a, const mat<T, M, 1>& b) {
    return {a * b.x};
}

// Matrix operations.
template <typename T, int N>
constexpr inline mat<T, N, 2> operator+(
    const mat<T, N, 2>& a, const mat<T, N, 2>& b) {
    return {a.x + b.x, a.y + b.y};
}
template <typename T, int N>
constexpr inline mat<T, N, 2> operator*(const mat<T, N, 2>& a, T b) {
    return {a.x * b, a.y * b};
}
template <typename T, int N>
constexpr inline vec<T, 2> operator*(const mat<T, N, 2>& a, const vec<T, 2>& b) {
    return a.x * b.x + a.y * b.y;
}
template <typename T, int N>
constexpr inline vec<T, 2> operator*(const vec<T, N>& a, const mat<T, N, 2>& b) {
    return {dot(a, b.x), dot(a, b.y)};
}
template <typename T, int N, int M>
constexpr inline mat<T, N, 2> operator*(
    const mat<T, N, M>& a, const mat<T, M, 2>& b) {
    return {a * b.x, a * b.y};
}

// Matrix operations.
template <typename T, int N>
constexpr inline mat<T, N, 3> operator+(
    const mat<T, N, 3>& a, const mat<T, N, 3>& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
template <typename T, int N>
constexpr inline mat<T, N, 3> operator*(const mat<T, N, 3>& a, T b) {
    return {a.x * b, a.y * b, a.z * b};
}
template <typename T, int N>
constexpr inline vec<T, 3> operator*(const mat<T, N, 3>& a, const vec<T, 3>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
template <typename T, int N>
constexpr inline vec<T, 3> operator*(const vec<T, N> a, const mat<T, N, 3>& b) {
    return {dot(a, b.x), dot(a, b.y), dot(a, b.z)};
}
template <typename T, int N, int M>
constexpr inline mat<T, N, 3> operator*(
    const mat<T, N, M>& a, const mat<T, M, 3>& b) {
    return {a * b.x, a * b.y, a * b.z};
}

// Matrix operations.
template <typename T, int N>
constexpr inline mat<T, N, 4> operator+(
    const mat<T, N, 4>& a, const mat<T, N, 4>& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}
template <typename T, int N>
constexpr inline mat<T, N, 4> operator*(const mat<T, N, 4>& a, T b) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
}
template <typename T, int N>
constexpr inline vec<T, 4> operator*(const mat<T, N, 4>& a, const vec<T, 4>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}
template <typename T, int N>
constexpr inline vec<T, 4> operator*(const vec<T, N>& a, const mat<T, N, 4>& b) {
    return {dot(a, b.x), dot(a, b.y), dot(a, b.z), dot(a, b.w)};
}
template <typename T, int N, int M>
constexpr inline mat<T, N, 4> operator*(
    const mat<T, N, M>& a, const mat<T, M, 4>& b) {
    return {a * b.x, a * b.y, a * b.z, a * b.w};
}

// Matrix assignments.
template <typename T, int N, int M>
constexpr inline mat<T, N, M>& operator+=(mat<T, N, M>& a, const mat<T, N, M>& b) {
    return a = a + b;
}
template <typename T, int N, int M>
constexpr inline mat<T, N, M>& operator*=(mat<T, N, M>& a, const mat<T, N, M>& b) {
    return a = a * b;
}
template <typename T, int N, int M, typename T1>
constexpr inline mat<T, N, M>& operator*=(mat<T, N, M>& a, T1 b) {
    return a = a * b;
}

// Matrix diagonals and transposes.
template <typename T>
constexpr inline vec<T, 1> diagonal(const mat<T, 1, 1>& a) {
    return {a.x.x};
}
template <typename T>
constexpr inline vec<T, 2> diagonal(const mat<T, 2, 2>& a) {
    return {a.x.x, a.y.y};
}
template <typename T>
constexpr inline vec<T, 3> diagonal(const mat<T, 3, 3>& a) {
    return {a.x.x, a.y.y, a.z.z};
}
template <typename T>
constexpr inline vec<T, 4> diagonal(const mat<T, 4, 4>& a) {
    return {a.x.x, a.y.y, a.z.z, a.w.w};
}
template <typename T>
constexpr inline mat<T, 1, 1> transpose(const mat<T, 1, 1>& a);
template <typename T>
constexpr inline mat<T, 2, 2> transpose(const mat<T, 2, 2>& a);
template <typename T>
constexpr inline mat<T, 3, 3> transpose(const mat<T, 3, 3>& a);
template <typename T>
constexpr inline mat<T, 4, 4> transpose(const mat<T, 4, 4>& a);

// Matrix adjugates, determinant and inverses.
template <typename T>
constexpr inline mat<T, 1, 1> adjugate(const mat<T, 1, 1>& a);
template <typename T>
constexpr inline mat<T, 2, 2> adjugate(const mat<T, 2, 2>& a);
template <typename T>
constexpr inline mat<T, 3, 3> adjugate(const mat<T, 3, 3>& a);
template <typename T>
constexpr inline mat<T, 4, 4> adjugate(const mat<T, 4, 4>& a);
template <typename T>
constexpr inline T determinant(const mat<T, 1, 1>& a);
template <typename T>
constexpr inline T determinant(const mat<T, 2, 2>& a);
template <typename T>
constexpr inline T determinant(const mat<T, 3, 3>& a);
template <typename T>
constexpr inline T determinant(const mat<T, 4, 4>& a);
template <typename T, int N>
constexpr inline mat<T, N, N> comatrix(const mat<T, N, N>& a) {
    return transpose(adjugate(a));
}
template <typename T, int N>
constexpr inline mat<T, N, N> inverse(const mat<T, N, N>& a) {
    return adjugate(a) * (1 / determinant(a));
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// RIGID BODY TRANSFORMS/FRAMES
// -----------------------------------------------------------------------------
namespace ygl {

// Rigid frames stored as a column-major affine transform matrix.
template <typename T, int N>
struct frame;

// Rigid frames stored as a column-major affine transform matrix.
template <typename T>
struct frame<T, 2> {
    vec<T, 2> x = {1, 0};
    vec<T, 2> y = {0, 1};
    vec<T, 2> o = {0, 0};
};
template <typename T>
struct frame<T, 3> {
    vec<T, 3> x = {1, 0, 0};
    vec<T, 3> y = {0, 1, 0};
    vec<T, 3> z = {0, 0, 1};
    vec<T, 3> o = {0, 0, 0};
};

// Type aliases.
using frame2f = frame<float, 2>;
using frame3f = frame<float, 3>;

// Indentity frames.
const auto identity_frame2f = frame2f{{1, 0}, {0, 1}, {0, 0}};
const auto identity_frame3f = frame3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};

// Frame construction from axis.
template <typename T>
constexpr inline frame<T, 3> make_frame_fromz(
    const vec<T, 3>& o, const vec<T, 3>& v) {
    auto z = normalize(v);
    auto x = normalize(orthogonal(z));
    auto y = normalize(cross(z, x));
    return {x, y, z, o};
}
template <typename T>
constexpr inline frame<T, 3> make_frame_fromzx(
    const vec<T, 3>& o, const vec<T, 3>& z_, const vec<T, 3>& x_) {
    auto z = normalize(z_);
    auto x = orthonormalize(x_, z);
    auto y = normalize(cross(z, x));
    return {x, y, z, o};
}

// Frame to matrix conversion.
template <typename T>
constexpr inline mat<T, 3, 3> frame_to_mat(const frame<T, 2>& a) {
    return {
        {a.x.x, a.x.y, 0},
        {a.y.x, a.y.y, 0},
        {a.o.x, a.o.y, 1},
    };
}
template <typename T>
constexpr inline frame<T, 2> mat_to_frame(const mat<T, 3, 3>& a) {
    return {
        {a.x.x, a.x.y},
        {a.y.x, a.y.y},
        {a.z.x, a.z.y},
    };
}
template <typename T>
constexpr inline mat<T, 4, 4> frame_to_mat(const frame<T, 3>& a) {
    return {
        {a.x.x, a.x.y, a.x.z, 0},
        {a.y.x, a.y.y, a.y.z, 0},
        {a.z.x, a.z.y, a.z.z, 0},
        {a.o.x, a.o.y, a.o.z, 1},
    };
}
template <typename T>
constexpr inline frame<T, 3> mat_to_frame(const mat<T, 4, 4>& a) {
    return {
        {a.x.x, a.x.y, a.x.z},
        {a.y.x, a.y.y, a.y.z},
        {a.z.x, a.z.y, a.z.z},
        {a.w.x, a.w.y, a.w.z},
    };
}

// Frame comparisons.
template <typename T>
constexpr inline bool operator==(const frame<T, 2>& a, const frame<T, 2>& b) {
    return a.x == b.x && a.y == b.y && a.o == b.o;
}
template <typename T>
constexpr inline bool operator!=(const frame<T, 2>& a, const frame<T, 2>& b) {
    return !(a == b);
}
template <typename T>
constexpr inline bool operator==(const frame<T, 3>& a, const frame<T, 3>& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.o == b.o;
}
template <typename T>
constexpr inline bool operator!=(const frame<T, 3>& a, const frame<T, 3>& b) {
    return !(a == b);
}

// Frame composition, equivalent to affine matrix product.
template <typename T>
constexpr inline frame<T, 2> operator*(
    const frame<T, 2>& a, const frame<T, 2>& b) {
    auto rot = mat<T, 2, 2>{a.x, a.y} * mat<T, 2, 2>{b.x, b.y};
    auto pos = mat<T, 2, 2>{a.x, a.y} * b.o + a.o;
    return {rot.x, rot.y, pos};
}
template <typename T>
constexpr inline frame<T, 3> operator*(
    const frame<T, 3>& a, const frame<T, 3>& b) {
    auto rot = mat<T, 3, 3>{a.x, a.y, a.z} * mat<T, 3, 3>{b.x, b.y, b.z};
    auto pos = mat<T, 3, 3>{a.x, a.y, a.z} * b.o + a.o;
    return {rot.x, rot.y, rot.z, pos};
}
// Frame inverse, equivalent to rigid affine inverse.
template <typename T>
constexpr inline frame<T, 2> inverse(const frame<T, 2>& a, bool is_rigid = true) {
    auto minv = (is_rigid) ? transpose(mat<T, 2, 2>{a.x, a.y}) :
                             inverse(mat<T, 2, 2>{a.x, a.y});
    return {minv.x, minv.y, -(minv * a.o)};
}
template <typename T>
constexpr inline frame<T, 3> inverse(const frame<T, 3>& a, bool is_rigid = true) {
    auto minv = (is_rigid) ? transpose(mat<T, 3, 3>{a.x, a.y, a.z}) :
                             inverse(mat<T, 3, 3>{a.x, a.y, a.z});
    return {minv.x, minv.y, minv.z, -(minv * a.o)};
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// AXIS ALIGNED BOUNDING BOXES
// -----------------------------------------------------------------------------
namespace ygl {

// DEscribes a range of values in N dimensions.
template <typename T, int N>
struct bbox;

// Range of values in 1D.
template <typename T>
struct bbox<T, 1> {
    T min = maxt<T>();
    T max = mint<T>();
};

// Axis aligned bounding box represented as a min/max vector pairs.
template <typename T>
struct bbox<T, 2> {
    vec<T, 2> min = {maxt<T>(), maxt<T>()};
    vec<T, 2> max = {mint<T>(), mint<T>()};
};

// Axis aligned bounding box represented as a min/max vector pairs.
template <typename T>
struct bbox<T, 3> {
    vec<T, 3> min = {maxt<T>(), maxt<T>(), maxt<T>()};
    vec<T, 3> max = {mint<T>(), mint<T>(), mint<T>()};
};

// Axis aligned bounding box represented as a min/max vector pairs.
template <typename T>
struct bbox<T, 4> {
    vec<T, 4> min = {maxt<T>(), maxt<T>(), maxt<T>(), maxt<T>()};
    vec<T, 4> max = {mint<T>(), mint<T>(), mint<T>(), mint<T>()};
};

// Type aliases
using bbox1f = bbox<float, 1>;
using bbox2f = bbox<float, 2>;
using bbox3f = bbox<float, 3>;
using bbox4f = bbox<float, 4>;

// Empty bbox constant.
const auto invalid_bbox1f = bbox1f();
const auto invalid_bbox2f = bbox2f();
const auto invalid_bbox3f = bbox3f();
const auto invalid_bbox4f = bbox4f();

// Bounding box comparisons.
template <typename T>
constexpr inline bool operator==(const bbox<T, 1>& a, const bbox<T, 1>& b) {
    return a.min == b.min && a.max == b.max;
}
template <typename T>
constexpr inline bool operator!=(const bbox<T, 1>& a, const bbox<T, 1>& b) {
    return a.min != b.min || a.max != b.max;
}
template <typename T>
constexpr inline bool operator==(const bbox<T, 2>& a, const bbox<T, 2>& b) {
    return a.min == b.min && a.max == b.max;
}
template <typename T>
constexpr inline bool operator!=(const bbox<T, 2>& a, const bbox<T, 2>& b) {
    return a.min != b.min || a.max != b.max;
}
template <typename T>
constexpr inline bool operator==(const bbox<T, 3>& a, const bbox<T, 3>& b) {
    return a.min == b.min && a.max == b.max;
}
template <typename T>
constexpr inline bool operator!=(const bbox<T, 3>& a, const bbox<T, 3>& b) {
    return a.min != b.min || a.max != b.max;
}
template <typename T>
constexpr inline bool operator==(const bbox<T, 4>& a, const bbox<T, 4>& b) {
    return a.min == b.min && a.max == b.max;
}
template <typename T>
constexpr inline bool operator!=(const bbox<T, 4>& a, const bbox<T, 4>& b) {
    return a.min != b.min || a.max != b.max;
}

// Bounding box expansions with points and other boxes.
template <typename T>
constexpr inline bbox<T, 1>& operator+=(bbox<T, 1>& a, T b) {
    a.min = min(a.min, b);
    a.max = max(a.max, b);
    return a;
}
template <typename T>
constexpr inline bbox<T, 1>& operator+=(bbox<T, 1>& a, const bbox<T, 1>& b) {
    a.min = min(a.min, b.min);
    a.max = max(a.max, b.max);
    return a;
}
// Bounding box expansions with points and other boxes.
template <typename T>
constexpr inline bbox<T, 2>& operator+=(bbox<T, 2>& a, const vec<T, 2>& b) {
    a.min = {min(a.min.x, b.x), min(a.min.y, b.y)};
    a.max = {max(a.max.x, b.x), max(a.max.y, b.y)};
    return a;
}
template <typename T>
constexpr inline bbox<T, 2>& operator+=(bbox<T, 2>& a, const bbox<T, 2>& b) {
    a.min = {min(a.min.x, b.min.x), min(a.min.y, b.min.y)};
    a.max = {max(a.max.x, b.max.x), max(a.max.y, b.max.y)};
    return a;
}
// Bounding box expansions with points and other boxes.
template <typename T>
constexpr inline bbox<T, 3>& operator+=(bbox<T, 3>& a, const vec<T, 3>& b) {
    a.min = {min(a.min.x, b.x), min(a.min.y, b.y), min(a.min.z, b.z)};
    a.max = {max(a.max.x, b.x), max(a.max.y, b.y), max(a.max.z, b.z)};
    return a;
}
template <typename T>
constexpr inline bbox<T, 3>& operator+=(bbox<T, 3>& a, const bbox<T, 3>& b) {
    a.min = {min(a.min.x, b.min.x), min(a.min.y, b.min.y), min(a.min.z, b.min.z)};
    a.max = {max(a.max.x, b.max.x), max(a.max.y, b.max.y), max(a.max.z, b.max.z)};
    return a;
}
// Bounding box expansions with points and other boxes.
template <typename T>
constexpr inline bbox<T, 4>& operator+=(bbox<T, 4>& a, const vec<T, 4>& b) {
    a.min = {min(a.min.x, b.x), min(a.min.y, b.y), min(a.min.z, b.z),
        min(a.min.w, b.w)};
    a.max = {max(a.max.x, b.x), max(a.max.y, b.y), max(a.max.z, b.z),
        max(a.max.w, b.w)};
    return a;
}
template <typename T>
constexpr inline bbox<T, 4>& operator+=(bbox<T, 4>& a, const bbox<T, 4>& b) {
    a.min = {min(a.min.x, b.min.x), min(a.min.y, b.min.y),
        min(a.min.z, b.min.z), min(a.min.w, b.min.w)};
    a.max = {max(a.max.x, b.max.x), max(a.max.y, b.max.y),
        max(a.max.z, b.max.z), max(a.max.w, b.max.w)};
    return a;
}

// Primitive bounds.
template <typename T>
constexpr inline bbox<T, 3> point_bounds(const vec<T, 3>& p, T r = 0) {
    auto bounds = bbox<T, 3>{};
    bounds += p - vec3f{r, r, r};
    bounds += p + vec3f{r, r, r};
    return bounds;
}
template <typename T>
constexpr inline bbox<T, 3> line_bounds(
    const vec<T, 3>& p0, const vec<T, 3>& p1, T r0 = 0, T r1 = 0) {
    auto bounds = bbox<T, 3>{};
    bounds += p0 - vec3f{r0, r0, r0};
    bounds += p0 + vec3f{r0, r0, r0};
    bounds += p1 - vec3f{r1, r1, r1};
    bounds += p1 + vec3f{r1, r1, r1};
    return bounds;
}
template <typename T>
constexpr inline bbox<T, 3> triangle_bounds(
    const vec<T, 3>& p0, const vec<T, 3>& p1, const vec<T, 3>& p2) {
    auto bounds = bbox<T, 3>{};
    bounds += p0;
    bounds += p1;
    bounds += p2;
    return bounds;
}
template <typename T>
constexpr inline bbox<T, 3> quad_bounds(const vec<T, 3>& p0,
    const vec<T, 3>& p1, const vec<T, 3>& p2, const vec<T, 3>& p3) {
    auto bounds = bbox<T, 3>{};
    bounds += p0;
    bounds += p1;
    bounds += p2;
    bounds += p3;
    return bounds;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// RAYS
// -----------------------------------------------------------------------------
namespace ygl {

// Rays with origin, direction and min/max t value.
template <typename T, int N>
struct ray;

// Rays with origin, direction and min/max t value.
template <typename T>
struct ray<T, 2> {
    vec<T, 2> o    = {0, 0};
    vec<T, 2> d    = {0, 1};
    T         tmin = 0;
    T         tmax = maxt<T>();
};
template <typename T>
struct ray<T, 3> {
    vec<T, 3> o    = {0, 0, 0};
    vec<T, 3> d    = {0, 0, 1};
    T         tmin = 0;
    T         tmax = maxt<T>();
};

// Type aliases.
using ray2f = ray<float, 2>;
using ray3f = ray<float, 3>;

// Construct a ray from direction or segments using a default epsilon.
template <typename T, int N>
constexpr inline ray<T, N> make_ray(
    const vec<T, N>& o, const vec<T, N>& d, T eps = 1e-4f) {
    return {o, d, eps, maxt<T>()};
}
template <typename T, int N>
constexpr inline ray<T, N> make_segment(
    const vec<T, N>& p1, const vec<T, N>& p2, T eps = 1e-4f) {
    return {p1, normalize(p2 - p1), eps, length(p2 - p1) - 2 * eps};
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// TRANSFORMS
// -----------------------------------------------------------------------------
namespace ygl {

// Transforms points, vectors and directions by matrices.
template <typename T>
constexpr inline vec<T, 2> transform_point(
    const mat<T, 3, 3>& a, const vec<T, 2>& b) {
    auto tvb = a * vec3f{b.x, b.y, 1};
    return vec<T, 2>{tvb.x, tvb.y} / tvb.z;
}
template <typename T>
constexpr inline vec<T, 3> transform_point(
    const mat<T, 4, 4>& a, const vec<T, 3>& b) {
    auto tvb = a * vec<T, 4>{b.x, b.y, b.z, 1};
    return vec<T, 3>{tvb.x, tvb.y, tvb.z} / tvb.w;
}
template <typename T>
constexpr inline vec<T, 2> transform_vector(
    const mat<T, 3, 3>& a, const vec<T, 2>& b) {
    auto tvb = a * vec<T, 3>{b.x, b.y, 0};
    return vec<T, 2>{tvb.x, tvb.y} / tvb.z;
}
template <typename T>
constexpr inline vec<T, 3> transform_vector(
    const mat<T, 3, 3>& a, const vec<T, 3>& b) {
    return a * b;
}
template <typename T>
constexpr inline vec<T, 3> transform_vector(
    const mat<T, 4, 4>& a, const vec<T, 3>& b) {
    auto tvb = a * vec4f{b.x, b.y, b.z, 0};
    return vec<T, 3>{tvb.x, tvb.y, tvb.z};
}
template <typename T, int N>
constexpr inline vec<T, N> transform_direction(
    const mat<T, 4, 4>& a, const vec<T, N>& b) {
    return normalize(transform_vector(a, b));
}

// Transforms points, vectors and directions by frames.
template <typename T>
constexpr inline vec<T, 2> transform_point(
    const frame<T, 2>& a, const vec<T, 2>& b) {
    return a.x * b.x + a.y * b.y + a.o;
}
template <typename T>
constexpr inline vec<T, 3> transform_point(
    const frame<T, 3>& a, const vec<T, 3>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.o;
}
template <typename T>
constexpr inline vec<T, 2> transform_vector(
    const frame<T, 2>& a, const vec<T, 2>& b) {
    return a.x * b.x + a.y * b.y;
}
template <typename T>
constexpr inline vec<T, 3> transform_vector(
    const frame<T, 3>& a, const vec<T, 3>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
template <typename T, int N>
constexpr inline vec<T, N> transform_direction(
    const frame<T, N>& a, const vec<T, N>& b) {
    return normalize(transform_vector(a, b));
}

// Transforms rays and bounding boxes by matrices.
template <typename T, int N>
constexpr inline ray<T, N> transform_ray(
    const frame<T, N>& a, const ray<T, N>& b) {
    return {transform_point(a, b.o), transform_vector(a, b.d), b.tmin, b.tmax};
}
template <typename T, int N>
constexpr inline ray<T, N> transform_ray(
    const mat<T, N + 1, N + 1>& a, const ray<T, N>& b) {
    return {transform_point(a, b.o), transform_vector(a, b.d), b.tmin, b.tmax};
}
template <typename T>
constexpr inline bbox<T, 3> transform_bbox(
    const frame<T, 3>& a, const bbox<T, 3>& b) {
    auto corners = {vec3f{b.min.x, b.min.y, b.min.z},
        vec3f{b.min.x, b.min.y, b.max.z}, vec3f{b.min.x, b.max.y, b.min.z},
        vec3f{b.min.x, b.max.y, b.max.z}, vec3f{b.max.x, b.min.y, b.min.z},
        vec3f{b.max.x, b.min.y, b.max.z}, vec3f{b.max.x, b.max.y, b.min.z},
        vec3f{b.max.x, b.max.y, b.max.z}};
    auto xformed = bbox<T, 3>();
    for (auto& corner : corners) xformed += transform_point(a, corner);
    return xformed;
}
template <typename T>
constexpr inline bbox<T, 3> transform_bbox(
    const mat<T, 4, 4>& a, const bbox<T, 3>& b) {
    auto corners = {vec3f{b.min.x, b.min.y, b.min.z},
        vec3f{b.min.x, b.min.y, b.max.z}, vec3f{b.min.x, b.max.y, b.min.z},
        vec3f{b.min.x, b.max.y, b.max.z}, vec3f{b.max.x, b.min.y, b.min.z},
        vec3f{b.max.x, b.min.y, b.max.z}, vec3f{b.max.x, b.max.y, b.min.z},
        vec3f{b.max.x, b.max.y, b.max.z}};
    auto xformed = bbox<T, 3>();
    for (auto& corner : corners) xformed += transform_point(a, corner);
    return xformed;
}

// Inverse transforms by frames, assuming they are rigid transforms.
template <typename T>
constexpr inline vec<T, 2> transform_point_inverse(
    const frame<T, 2>& a, const vec<T, 2>& b) {
    return {dot(b - a.o, a.x), dot(b - a.o, a.y)};
}
template <typename T>
constexpr inline vec3f transform_point_inverse(
    const frame<T, 3>& a, const vec<T, 3>& b) {
    return {dot(b - a.o, a.x), dot(b - a.o, a.y), dot(b - a.o, a.z)};
}
template <typename T>
constexpr inline vec<T, 2> transform_vector_inverse(
    const frame<T, 2>& a, const vec<T, 2>& b) {
    return {dot(b, a.x), dot(b, a.y)};
}
template <typename T>
constexpr inline vec3f transform_vector_inverse(
    const frame<T, 3>& a, const vec<T, 3>& b) {
    return {dot(b, a.x), dot(b, a.y), dot(b, a.z)};
}
template <typename T, int N>
constexpr inline vec3f transform_direction_inverse(
    const frame<T, N>& a, const vec<T, N>& b) {
    return normalize(transform_vector_inverse(a, b));
}
template <typename T, int N>
constexpr inline ray<T, N> transform_ray_inverse(
    const frame<T, N>& a, const ray<T, N>& b) {
    return {transform_point_inverse(a, b.o),
        transform_direction_inverse(a, b.d), b.tmin, b.tmax};
}
template <typename T>
constexpr inline bbox<T, 3> transform_bbox_inverse(
    const frame<T, 3>& a, const bbox<T, 3>& b) {
    return transform_bbox(inverse(a), b);
}

// Translation, scaling and rotations transforms.
template <typename T>
constexpr inline frame<T, 3> translation_frame(const vec<T, 3>& a) {
    return {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, a};
}
template <typename T>
constexpr inline frame<T, 3> scaling_frame(const vec<T, 3>& a) {
    return {{a.x, 0, 0}, {0, a.y, 0}, {0, 0, a.z}, {0, 0, 0}};
}
template <typename T>
constexpr inline frame<T, 3> rotation_frame(const vec<T, 3>& axis, T angle) {
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
template <typename T>
constexpr inline frame<T, 3> rotation_frame(const vec<T, 4>& quat) {
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
template <typename T>
constexpr inline frame<T, 3> rotation_frame(const mat<T, 3, 3>& rot) {
    return {rot.x, rot.y, rot.z, {0, 0, 0}};
}

// Lookat frame. Z-axis can be inverted with inv_xz.
template <typename T>
constexpr inline frame<T, 3> lookat_frame(const vec<T, 3>& eye,
    const vec<T, 3>& center, const vec<T, 3>& up, bool inv_xz = false) {
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
template <typename T>
constexpr inline mat<T, 4, 4> frustum_mat(T l, T r, T b, T t, T n, T f) {
    return {{2 * n / (r - l), 0, 0, 0}, {0, 2 * n / (t - b), 0, 0},
        {(r + l) / (r - l), (t + b) / (t - b), -(f + n) / (f - n), -1},
        {0, 0, -2 * f * n / (f - n), 0}};
}
template <typename T>
constexpr inline mat<T, 4, 4> orthographic_mat(T l, T r, T b, T t, T n, T f) {
    return {{2 / (r - l), 0, 0, 0}, {0, 2 / (t - b), 0, 0},
        {0, 0, -2 / (f - n), 0},
        {-(r + l) / (r - l), -(t + b) / (t - b), -(f + n) / (f - n), 1}};
}
template <typename T>
constexpr inline mat<T, 4, 4> orthographic2d_mat(
    T left, T right, T bottom, T top) {
    return orthographic_mat(left, right, bottom, top, -1, 1);
}
template <typename T>
constexpr inline mat<T, 4, 4> orthographic_mat(T xmag, T ymag, T near, T far) {
    return {{1 / xmag, 0, 0, 0}, {0, 1 / ymag, 0, 0},
        {0, 0, 2 / (near - far), 0}, {0, 0, (far + near) / (near - far), 1}};
}
template <typename T>
constexpr inline mat<T, 4, 4> perspective_mat(T fovy, T aspect, T near, T far) {
    auto tg = tan(fovy / 2);
    return {{1 / (aspect * tg), 0, 0, 0}, {0, 1 / tg, 0, 0},
        {0, 0, (far + near) / (near - far), -1},
        {0, 0, 2 * far * near / (near - far), 0}};
}
template <typename T>
constexpr inline mat<T, 4, 4> perspective_mat(T fovy, T aspect, T near) {
    auto tg = tan(fovy / 2);
    return {{1 / (aspect * tg), 0, 0, 0}, {0, 1 / tg, 0, 0}, {0, 0, -1, -1},
        {0, 0, 2 * near, 0}};
}

// Rotation conversions.
template <typename T>
constexpr inline tuple<vec<T, 3>, T> rotation_axisangle(const vec<T, 4>& quat) {
    return {normalize(vec3f{quat.x, quat.y, quat.z}), 2 * acos(quat.w)};
}
template <typename T>
constexpr inline vec<T, 4> rotation_quat(const vec<T, 3>& axis, T angle) {
    auto len = length(axis);
    if (!len) return {0, 0, 0, 1};
    return vec4f{sin(angle / 2) * axis.x / len, sin(angle / 2) * axis.y / len,
        sin(angle / 2) * axis.z / len, cos(angle / 2)};
}
template <typename T>
constexpr inline vec<T, 4> rotation_quat(const vec<T, 4>& axisangle) {
    return rotation_quat(
        vec3f{axisangle.x, axisangle.y, axisangle.z}, axisangle.w);
}

// Turntable and FPS Camera navigation.
inline void camera_turntable(vec3f& from, vec3f& to, vec3f& up,
    const vec2f& rotate, float dolly, const vec2f& pan);
inline void camera_turntable(frame3f& frame, float& focus, const vec2f& rotate,
    float dolly, const vec2f& pan);
inline void camera_fps(frame3f& frame, const vec3f& transl, const vec2f& rotate);

// Computes the image uv coordinates corresponding to the view parameters.
// Returns negative coordinates if out of the image.
inline vec2i get_image_coords(const vec2f& mouse_pos, const vec2f& center,
    float scale, const vec2i& txt_size) {
    auto xyf = (mouse_pos - center) / scale;
    return vec2i{(int)round(xyf.x + txt_size.x / 2.0f),
        (int)round(xyf.y + txt_size.y / 2.0f)};
}

// Center image and autofit.
inline void center_image(vec2f& center, float& scale, const vec2i& image_size,
    const vec2i& window_size, bool zoom_to_fit) {
    if (zoom_to_fit) {
        scale  = min(window_size.x / (float)image_size.x,
            window_size.y / (float)image_size.y);
        center = {(float)window_size.x / 2, (float)window_size.y / 2};
    } else {
        if (window_size.x >= image_size.x * scale) center.x = window_size.x / 2;
        if (window_size.y >= image_size.y * scale) center.y = window_size.y / 2;
    }
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// CONTAINER UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// check if a container contains a key
template <typename K, typename V>
inline bool contains(const unordered_map<K, V>& container, const K& value) {
    return container.find(value) != container.end();
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// CONCURRENCY UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// a simple concurrent queue that locks at every call
template <typename T>
struct concurrent_queue {
    concurrent_queue() {}
    concurrent_queue(const concurrent_queue& other) {
        if (!other._queue.empty()) log_error("cannot copy full queue");
        clear();
    }
    concurrent_queue& operator=(const concurrent_queue& other) {
        if (!other._queue.empty()) log_error("cannot copy full queue");
        clear();
    }

    bool empty() {
        lock_guard<mutex> lock(_mutex);
        return _queue.empty();
    }
    void clear() {
        lock_guard<mutex> lock(_mutex);
        _queue.clear();
    }
    void push(const T& value) {
        lock_guard<mutex> lock(_mutex);
        _queue.push_back(value);
    }
    bool try_pop(T& value) {
        lock_guard<mutex> lock(_mutex);
        if (_queue.empty()) return false;
        value = _queue.front();
        _queue.pop_front();
        return true;
    }

    mutex    _mutex;
    deque<T> _queue;
};

}  // namespace ygl

// -----------------------------------------------------------------------------
// RANDOM NUMBER GENERATION
// -----------------------------------------------------------------------------
namespace ygl {

// PCG random numbers from http://www.pcg-random.org/
struct rng_state {
    uint64_t state = 0x853c49e6748fea9bULL;
    uint64_t inc   = 0xda3e39cb94b95bdbULL;
};

// Next random number.
inline uint32_t advance_rng(rng_state& rng) {
    uint64_t oldstate   = rng.state;
    rng.state           = oldstate * 6364136223846793005ULL + rng.inc;
    uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
    uint32_t rot        = (uint32_t)(oldstate >> 59u);
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
}

// Init a random number generator with a state state from the sequence seq.
inline rng_state make_rng(uint64_t seed, uint64_t seq = 1) {
    auto rng  = rng_state();
    rng.state = 0U;
    rng.inc   = (seq << 1u) | 1u;
    advance_rng(rng);
    rng.state += seed;
    advance_rng(rng);
    return rng;
}

// Next random numbers: floats in [0,1), ints in [0,n).
inline int get_random_int(rng_state& rng, int n) {
    return advance_rng(rng) % n;
}
inline float get_random_float(rng_state& rng) {
    union {
        uint32_t u;
        float    f;
    } x;
    x.u = (advance_rng(rng) >> 9) | 0x3f800000u;
    return x.f - 1.0f;
    // alternate implementation
    // const static auto scale = (float)(1.0 / numeric_limits<uint32_t>::max());
    // return advance_rng(rng) * scale;
}
inline vec2f get_random_vec2f(rng_state& rng) {
    // force order of evaluation by using separate assignments.
    auto x = get_random_float(rng);
    auto y = get_random_float(rng);
    return {x, y};
}
inline vec3f get_random_vec3f(rng_state& rng) {
    // force order of evaluation by using separate assignments.
    auto x = get_random_float(rng);
    auto y = get_random_float(rng);
    auto z = get_random_float(rng);
    return {x, y, z};
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// MONETACARLO SAMPLING FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Sample an hemispherical direction with uniform distribution.
inline vec3f sample_hemisphere_direction(const vec2f& ruv) {
    auto z   = ruv.y;
    auto r   = sqrt(1 - z * z);
    auto phi = 2 * pif * ruv.x;
    return {r * cos(phi), r * sin(phi), z};
}
inline float sample_hemisphere_direction_pdf(const vec3f& direction) {
    return (direction.z <= 0) ? 0 : 1 / (2 * pif);
}

// Sample a spherical direction with uniform distribution.
inline vec3f sample_sphere_direction(const vec2f& ruv) {
    auto z   = 2 * ruv.y - 1;
    auto r   = sqrt(1 - z * z);
    auto phi = 2 * pif * ruv.x;
    return {r * cos(phi), r * sin(phi), z};
}
inline float sample_sphere_direction_pdf(const vec3f& w) {
    return 1 / (4 * pif);
}

// Sample an hemispherical direction with cosine distribution.
inline vec3f sample_hemisphere_direction_cosine(const vec2f& ruv) {
    auto z   = sqrt(ruv.y);
    auto r   = sqrt(1 - z * z);
    auto phi = 2 * pif * ruv.x;
    return {r * cos(phi), r * sin(phi), z};
}
inline float sample_hemisphere_direction_cosine_pdf(const vec3f& direction) {
    return (direction.z <= 0) ? 0 : direction.z / pif;
}

// Sample an hemispherical direction with cosine power distribution.
inline vec3f sample_hemisphere_direction_cospower(
    float exponent, const vec2f& ruv) {
    auto z   = pow(ruv.y, 1 / (exponent + 1));
    auto r   = sqrt(1 - z * z);
    auto phi = 2 * pif * ruv.x;
    return {r * cos(phi), r * sin(phi), z};
}
inline float sample_hemisphere_direction_cospower_pdf(
    float exponent, const vec3f& direction) {
    return (direction.z <= 0) ?
               0 :
               pow(direction.z, exponent) * (exponent + 1) / (2 * pif);
}

// Sample a point uniformly on a disk.
inline vec3f sample_disk_point(const vec2f& ruv) {
    auto r   = sqrt(ruv.y);
    auto phi = 2 * pif * ruv.x;
    return {cos(phi) * r, sin(phi) * r, 0};
}
inline float sample_disk_point_pdf() { return 1 / pif; }

// Sample a point uniformly on a cylinder, without caps.
inline vec3f sample_cylinder_point(const vec2f& ruv) {
    auto phi = 2 * pif * ruv.x;
    return {sin(phi), cos(phi), ruv.y * 2 - 1};
}
inline float sample_cylinder_point_pdf() { return 1 / pif; }

// Sample a point uniformly on a triangle returning the baricentric coordinates.
inline vec2f sample_triangle_coordinates(const vec2f& ruv) {
    return {1 - sqrt(ruv.x), ruv.y * sqrt(ruv.x)};
}

// Sample a point uniformly on a triangle.
inline vec3f sample_triangle_point(
    const vec3f& p0, const vec3f& p1, const vec3f& p2, const vec2f& ruv) {
    auto uv = sample_triangle_coordinates(ruv);
    return p0 * (1 - uv.x - uv.y) + p1 * uv.x + p2 * uv.y;
}
// Pdf for uniform triangle sampling, i.e. triangle area.
inline float sample_triangle_point_pdf(
    const vec3f& p0, const vec3f& p1, const vec3f& p2) {
    return 2 / length(cross(p1 - p0, p2 - p0));
}

// Sample an index with uniform distribution.
inline int sample_uniform_index(int size, float r) {
    return clamp((int)(r * size), 0, size - 1);
}
inline float sample_uniform_index_pdf(int size) { return 1.0f / size; }

// Sample an index with uniform distribution.
template <typename T>
inline T sample_uniform_element(const vector<T>& elements, float r) {
    if (elements.empty()) return {};
    auto size = (int)elements.size();
    return elements[clamp((int)(r * size), 0, size - 1)];
}
template <typename T>
inline float sample_uniform_element_pdf(const vector<T>& elements) {
    if (elements.empty()) return 0;
    return 1.0f / (int)elements.size();
}

// Sample a discrete distribution represented by its cdf.
inline int sample_discrete_distribution(const vector<float>& cdf, float r) {
    r        = clamp(r * cdf.back(), 0.0f, cdf.back() - 0.00001f);
    auto idx = (int)(std::upper_bound(cdf.data(), cdf.data() + cdf.size(), r) -
                     cdf.data());
    return clamp(idx, 0, (int)cdf.size() - 1);
}
// Pdf for uniform discrete distribution sampling.
inline float sample_discrete_distribution_pdf(const vector<float>& cdf, int idx) {
    if (idx == 0) return cdf.at(0);
    return cdf.at(idx) - cdf.at(idx - 1);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// PERLIN NOISE FUNCTION
// -----------------------------------------------------------------------------
namespace ygl {

// Compute the revised Perlin noise function. Wrap provides a wrapping noise
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
inline vec3f line_tangent(const vec3f& p0, const vec3f& p1) {
    return normalize(p1 - p0);
}
inline float line_length(const vec3f& p0, const vec3f& p1) {
    return length(p1 - p0);
}

// Triangle properties.
inline vec3f triangle_normal(const vec3f& p0, const vec3f& p1, const vec3f& p2) {
    return normalize(cross(p1 - p0, p2 - p0));
}
inline float triangle_area(const vec3f& p0, const vec3f& p1, const vec3f& p2) {
    return length(cross(p1 - p0, p2 - p0)) / 2;
}

// Quad propeties.
inline vec3f quad_normal(
    const vec3f& p0, const vec3f& p1, const vec3f& p2, const vec3f& p3) {
    return normalize(triangle_normal(p0, p1, p3) + triangle_normal(p2, p3, p1));
}
inline float quad_area(
    const vec3f& p0, const vec3f& p1, const vec3f& p2, const vec3f& p3) {
    return triangle_area(p0, p1, p3) + triangle_area(p2, p3, p1);
}

// Triangle tangent and bitangent from uv
inline tuple<vec3f, vec3f> triangle_tangents_fromuv(const vec3f& p0,
    const vec3f& p1, const vec3f& p2, const vec2f& uv0, const vec2f& uv1,
    const vec2f& uv2) {
    // Follows the definition in http://www.terathon.com/code/tangent.html and
    // https://gist.github.com/aras-p/2843984
    // normal points up from texture space
    auto p   = p1 - p0;
    auto q   = p2 - p0;
    auto s   = vec2f{uv1.x - uv0.x, uv2.x - uv0.x};
    auto t   = vec2f{uv1.y - uv0.y, uv2.y - uv0.y};
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

// Interpolates values over a line parameterized from a to b by u. Same as lerp.
template <typename T>
inline T interpolate_line(const T& p0, const T& p1, float u) {
    return p0 * (1 - u) + p1 * u;
}
// Interpolates values over a triangle parameterized by u and v along the
// (p1-p0) and (p2-p0) directions. Same as barycentric interpolation.
template <typename T>
inline T interpolate_triangle(
    const T& p0, const T& p1, const T& p2, const vec2f& uv) {
    return p0 * (1 - uv.x - uv.y) + p1 * uv.x + p2 * uv.y;
}
// Interpolates values over a quad parameterized by u and v along the
// (p1-p0) and (p2-p1) directions. Same as bilinear interpolation.
template <typename T>
inline T interpolate_quad(
    const T& p0, const T& p1, const T& p2, const T& p3, const vec2f& uv) {
#if YGL_QUADS_AS_TRIANGLES
    if (uv.x + uv.y <= 1) {
        return interpolate_triangle(p0, p1, p3, uv);
    } else {
        return interpolate_triangle(p2, p3, p1, 1 - uv);
    }
#else
    return p0 * (1 - uv.x) * (1 - uv.y) + p1 * uv.x * (1 - uv.y) +
           p2 * uv.x * uv.y + p3 * (1 - uv.x) * uv.y;
#endif
}

// Interpolates values along a cubic Bezier segment parametrized by u.
template <typename T>
inline T interpolate_bezier(
    const T& p0, const T& p1, const T& p2, const T& p3, float u) {
    return p0 * (1 - u) * (1 - u) * (1 - u) + p1 * 3 * u * (1 - u) * (1 - u) +
           p2 * 3 * u * u * (1 - u) + p3 * u * u * u;
}
// Computes the derivative of a cubic Bezier segment parametrized by u.
template <typename T>
inline T interpolate_bezier_derivative(
    const T& p0, const T& p1, const T& p2, const T& p3, float u) {
    return (p1 - p0) * 3 * (1 - u) * (1 - u) + (p2 - p1) * 6 * u * (1 - u) +
           (p3 - p2) * 3 * u * u;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// SHAPE UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Compute per-vertex normals/tangents for lines/triangles/quads.
vector<vec3f> compute_vertex_tangents(
    const vector<vec2i>& lines, const vector<vec3f>& positions);
vector<vec3f> compute_vertex_normals(
    const vector<vec3i>& triangles, const vector<vec3f>& positions);
vector<vec3f> compute_vertex_normals(
    const vector<vec4i>& quads, const vector<vec3f>& positions);

// Compute per-vertex tangent space for triangle meshes.
// Tangent space is defined by a four component vector.
// The first three components are the tangent with respect to the u texcoord.
// The fourth component is the sign of the tangent wrt the v texcoord.
// Tangent frame is useful in normal mapping.
vector<vec4f> compute_tangent_spaces(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec2f>& texturecoords);

// Apply skinning to vertex position and normals.
tuple<vector<vec3f>, vector<vec3f>> compute_skinning(
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec4f>& weights, const vector<vec4i>& joints,
    const vector<frame3f>& xforms);
// Apply skinning as specified in Khronos glTF.
tuple<vector<vec3f>, vector<vec3f>> compute_matrix_skinning(
    const vector<vec3f>& positions, const vector<vec3f>& normals,
    const vector<vec4f>& weights, const vector<vec4i>& joints,
    const vector<mat4f>& xforms);

// Dictionary to store edge information.
using edge_map = unordered_map<vec2i, vec2i>;

// Initialize an edge map with elements.
edge_map make_edge_map(const vector<vec3i>& triangles);
edge_map make_edge_map(const vector<vec4i>& quads);
// Insert an edge and return its index
int insert_edge(edge_map& emap, const vec2i& edge);
// Get the edge index / insertion count
int get_edge_index(const edge_map& emap, const vec2i& edge);
int get_edge_count(const edge_map& emap, const vec2i& edge);
// Get list of edges / boundary edges
vector<vec2i> get_edges(const edge_map& emap);
vector<vec2i> get_boundary(const edge_map& emap);

// Create an array of edges.
inline vector<vec2i> get_edges(const vector<vec3i>& triangles) {
    return get_edges(make_edge_map(triangles));
}
inline vector<vec2i> get_edges(const vector<vec4i>& quads) {
    return get_edges(make_edge_map(quads));
}

// Convert quads to triangles
vector<vec3i> convert_quads_to_triangles(const vector<vec4i>& quads);
// Convert quads to triangles with a diamond-like topology.
// Quads have to be consecutive one row after another.
vector<vec3i> convert_quads_to_triangles(
    const vector<vec4i>& quads, int row_length);
// Convert triangles to quads by creating degenerate quads
vector<vec4i> convert_triangles_to_quads(const vector<vec3i>& triangles);

// Convert beziers to lines using 3 lines for each bezier.
vector<vec2i> convert_bezier_to_lines(const vector<vec4i>& beziers);

// Convert face-varying data to single primitives. Returns the quads indices
// and face ids and filled vectors for pos, norm and texcoord. When used
// with ids, it also plits the faces per id.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>> convert_face_varying(
    const vector<vec4i>& quads_positions, const vector<vec4i>& quads_normals,
    const vector<vec4i>& quads_texturecoords, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texturecoords);

// Split primitives per id
vector<vector<vec2i>> ungroup_lines(
    const vector<vec2i>& lines, const vector<int>& ids);
vector<vector<vec3i>> ungroup_triangles(
    const vector<vec3i>& triangles, const vector<int>& ids);
vector<vector<vec4i>> ungroup_quads(
    const vector<vec4i>& quads, const vector<int>& ids);

// Subdivide lines by splitting each line in half.
template <typename T>
tuple<vector<vec2i>, vector<T>> subdivide_lines(
    const vector<vec2i>& lines, const vector<T>& vert);
// Subdivide triangle by splitting each triangle in four, creating new
// vertices for each edge.
template <typename T>
tuple<vector<vec3i>, vector<T>> subdivide_triangles(
    const vector<vec3i>& triangles, const vector<T>& vert);
// Subdivide quads by splitting each quads in four, creating new
// vertices for each edge and for each face.
template <typename T>
tuple<vector<vec4i>, vector<T>> subdivide_quads(
    const vector<vec4i>& quads, const vector<T>& vert);
// Subdivide beziers by splitting each segment in two.
template <typename T>
tuple<vector<vec4i>, vector<T>> subdivide_beziers(
    const vector<vec4i>& beziers, const vector<T>& vert);
// Subdivide quads using Carmull-Clark subdivision rules.
template <typename T>
tuple<vector<vec4i>, vector<T>> subdivide_catmullclark(const vector<vec4i>& quads,
    const vector<T>& vert, bool lock_boundary = false);

// Weld vertices within a threshold. For noe the implementation is O(n^2).
tuple<vector<vec3f>, vector<int>> weld_vertices(
    const vector<vec3f>& positions, float threshold);
tuple<vector<vec3i>, vector<vec3f>> weld_triangles(const vector<vec3i>& triangles,
    const vector<vec3f>& positions, float threshold);
tuple<vector<vec4i>, vector<vec3f>> weld_quads(const vector<vec4i>& quads,
    const vector<vec3f>& positions, float threshold);

// Pick a point in a point set uniformly.
inline int sample_points_element(int npoints, float re) {
    return sample_uniform_index(npoints, re);
}
inline vector<float> sample_points_element_cdf(int npoints) {
    auto cdf = vector<float>(npoints);
    for (auto i = 0; i < cdf.size(); i++) cdf[i] = 1 + (i ? cdf[i - 1] : 0);
    return cdf;
}
inline int sample_points_element(const vector<float>& cdf, float re) {
    return sample_discrete_distribution(cdf, re);
}

// Pick a point on lines uniformly.
inline vector<float> sample_lines_element_cdf(
    const vector<vec2i>& lines, const vector<vec3f>& positions) {
    auto cdf = vector<float>(lines.size());
    for (auto i = 0; i < cdf.size(); i++) {
        auto l = lines[i];
        auto w = line_length(positions[l.x], positions[l.y]);
        cdf[i] = w + (i ? cdf[i - 1] : 0);
    }
    return cdf;
}
inline tuple<int, float> sample_lines_element(
    const vector<float>& cdf, float re, float ru) {
    return {sample_discrete_distribution(cdf, re), ru};
}

// Pick a point on a triangle mesh uniformly.
inline vector<float> sample_triangles_element_cdf(
    const vector<vec3i>& triangles, const vector<vec3f>& positions) {
    auto cdf = vector<float>(triangles.size());
    for (auto i = 0; i < cdf.size(); i++) {
        auto t = triangles[i];
        auto w = triangle_area(positions[t.x], positions[t.y], positions[t.z]);
        cdf[i] = w + (i ? cdf[i - 1] : 0);
    }
    return cdf;
}
inline tuple<int, vec2f> sample_triangles_element(
    const vector<float>& cdf, float re, const vec2f& ruv) {
    return {sample_discrete_distribution(cdf, re),
        sample_triangle_coordinates(ruv)};
}

// Pick a point on a quad mesh uniformly.
inline vector<float> sample_quads_element_cdf(
    const vector<vec4i>& quads, const vector<vec3f>& positions) {
    auto cdf = vector<float>(quads.size());
    for (auto i = 0; i < cdf.size(); i++) {
        auto q = quads[i];
        auto w = quad_area(
            positions[q.x], positions[q.y], positions[q.z], positions[q.w]);
        cdf[i] = w + (i ? cdf[i - 1] : 0);
    }
    return cdf;
}
inline tuple<int, vec2f> sample_quads_element(
    const vector<float>& cdf, float re, const vec2f& ruv) {
    return {sample_discrete_distribution(cdf, re), ruv};
}
inline tuple<int, vec2f> sample_quads_element(const vector<vec4i>& quads,
    const vector<float>& cdf, float re, const vec2f& ruv) {
    auto element_id = sample_discrete_distribution(cdf, re);
    if (quads[element_id].z == quads[element_id].w) {
        return {element_id, sample_triangle_coordinates(ruv)};
    } else {
        return {element_id, ruv};
    }
}

// Samples a set of points over a triangle/quad mesh uniformly. Returns pos,
// norm and texcoord of the sampled points.
tuple<vector<vec3f>, vector<vec3f>, vector<vec2f>> sample_triangles_points(
    const vector<vec3i>& triangles, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texturecoords,
    int npoints, int seed = 7);
tuple<vector<vec3f>, vector<vec3f>, vector<vec2f>> sample_quads_points(
    const vector<vec4i>& quads, const vector<vec3f>& positions,
    const vector<vec3f>& normals, const vector<vec2f>& texturecoords,
    int npoints, int seed = 7);

}  // namespace ygl

// -----------------------------------------------------------------------------
// RAY INTERSECTION AND CLOSEST POINT FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Intersect a ray with a point (approximate).
// Based on http://geomalgorithms.com/a02-lines.html.
bool intersect_point(
    const ray3f& ray, const vec3f& p, float r, float& distance, vec2f& uv);

// Intersect a ray with a line (approximate).
// Based on http://geomalgorithms.com/a05-intersect-1.html and
// http://geomalgorithms.com/a07-distance.html#
//     dist3D_Segment_to_Segment
bool intersect_line(const ray3f& ray, const vec3f& p0, const vec3f& p1,
    float r0, float r1, float& distance, vec2f& uv);

// Intersect a ray with a triangle.
bool intersect_triangle(const ray3f& ray, const vec3f& p0, const vec3f& p1,
    const vec3f& p2, float& distance, vec2f& uv);

// Intersect a ray with a quad represented as two triangles (0,1,3) and
// (2,3,1), with the uv coordinates of the second triangle corrected by u =
// 1-u' and v = 1-v' to produce a quad parametrization where u and v go from 0
// to 1. This is equivalent to Intel's Embree.
bool intersect_quad(const ray3f& ray, const vec3f& p0, const vec3f& p1,
    const vec3f& p2, const vec3f& p3, float& distance, vec2f& uv);

// Intersect a ray with a axis-aligned bounding box.
bool intersect_bbox(const ray3f& ray, const bbox3f& bbox);

// Intersect a ray with a axis-aligned bounding box, implemented as
// "Robust BVH Ray Traversal" by T. Ize published at
// http://jcgt.org/published/0002/02/02/paper.pdf
bool intersect_bbox(const ray3f& ray, const vec3f& ray_dinv,
    const vec3i& ray_dsign, const bbox3f& bbox);

// Check if a point overlaps a position within a max distance.
bool overlap_point(const vec3f& pos, float dist_max, const vec3f& p0, float r0,
    float& distance, vec2f& uv);

// Find closest line point to a position.
float closestuv_line(const vec3f& pos, const vec3f& p0, const vec3f& p1);

// Check if a line overlaps a position within a max distance.
bool overlap_line(const vec3f& pos, float dist_max, const vec3f& p0,
    const vec3f& p1, float r0, float r1, float& distance, vec2f& uv);

// Find closest triangle point to a position.
vec2f closestuv_triangle(
    const vec3f& pos, const vec3f& p0, const vec3f& p1, const vec3f& p2);

// Check if a triangle overlaps a position within a max distance.
bool overlap_triangle(const vec3f& pos, float dist_max, const vec3f& p0,
    const vec3f& p1, const vec3f& p2, float r0, float r1, float r2,
    float& distance, vec2f& uv);

// Check if a quad overlaps a position within a max distance.
bool overlap_quad(const vec3f& pos, float dist_max, const vec3f& p0,
    const vec3f& p1, const vec3f& p2, const vec3f& p3, float r0, float r1,
    float r2, float r3, float& distance, vec2f& uv);

// Check if a bounding box overlaps a position within a max distance.
bool overlap_bbox(const vec3f& pos, float dist_max, const bbox3f& bbox);

// Check if two bounding boxes overlap.
bool overlap_bbox(const bbox3f& bbox1, const bbox3f& bbox2);

}  // namespace ygl

// -----------------------------------------------------------------------------
// BVH FOR RAY INTERSECTION AND CLOSEST ELEMENT
// -----------------------------------------------------------------------------
namespace ygl {

// Maximum number of primitives per BVH node.
const int bvh_max_prims = 4;

// BVH tree node containing its bounds, indices to the BVH arrays of either
// primitives or internal nodes, the node element type,
// and the split axis. Leaf and internal nodes are identical, except that
// indices refer to primitives for leaf nodes or other nodes for internal
// nodes. See bvh_scene for more details.
struct bvh_node {
    bbox3f   bbox;
    uint16_t num_primitives;
    bool     is_internal;
    uint8_t  split_axis;
    uint32_t primitive_ids[bvh_max_prims];
};

// BVH for shapes made of points, lines, triangles or quads. Only one primitive
// type can be used.
// To build, fill in the shape data, then call `build_shape_bvh()`.
// The BVH is stored as a node array with the tree structure is encoded using
// array indices. BVH nodes indices refer to either the node array,
// for internal nodes, or the primitive arrays, for leaf nodes.
struct bvh_shape {
    // data for shape BVH
    vector<vec3f> positions;
    vector<float> radius;
    vector<int>   points;
    vector<vec2i> lines;
    vector<vec3i> triangles;
    vector<vec4i> quads;

    // bvh internal nodes
    vector<bvh_node> nodes;

    // Embree opaque data
    void* embree_bvh = nullptr;
};

// Instance for a scene BVH.
struct bvh_instance {
    frame3f frame         = identity_frame3f;
    frame3f frame_inverse = identity_frame3f;
    int     shape_id      = -1;
    int     surface_id    = -1;
};

// BVH for scenes made of instances to shapes.
// To build, first build the shape BVHs, then fill in the instance data,
// then call `build_scene_bvh()`.
// The BVH is stored as a node array with the tree structure is encoded using
// array indices. BVH nodes indices refer to either the node array,
// for internal nodes, or the primitive arrays, for leaf nodes.
struct bvh_scene {
    // data for instance BVH
    vector<bvh_instance> instances;
    vector<bvh_shape>    shape_bvhs;
    vector<bvh_shape>    surface_bvhs;

    // bvh internal nodes
    vector<bvh_node> nodes;

    // Embree opaque data
    void* embree_bvh = nullptr;
};

// Build a BVH from the given set of primitives.
void build_shape_bvh(bvh_shape& bvh, bool high_quality = false);
void build_scene_bvh(bvh_scene& bvh, bool high_quality = false);
// Update the node bounds for a shape bvh.
void refit_shape_bvh(bvh_shape& bvh);
void refit_scene_bvh(bvh_scene& bvh);

// Build a BVH from the given set of primitives.
// Uses Embree if available and requested, otherwise the standard build.
void build_shape_bvh_embree(bvh_shape& bvh, bool high_quality = false);
void clear_shape_bvh_embree(bvh_shape& bvh);
void build_scene_bvh_embree(bvh_scene& bvh, bool high_quality = false);
void clear_scene_bvh_embree(bvh_scene& bvh);

// Intersect ray with a bvh returning either the first or any intersection
// depending on `find_any`. Returns the ray distance , the instance id,
// the shape element index and the element barycentric coordinates.
bool intersect_shape_bvh(const bvh_shape& bvh, const ray3f& ray, bool find_any,
    float& distance, int& element_id, vec2f& element_uv);
bool intersect_scene_bvh(const bvh_scene& bvh, const ray3f& ray, bool find_any,
    float& distance, int& instance_id, int& element_id, vec2f& element_uv);

// Find a shape element that overlaps a point within a given distance
// max distance, returning either the closest or any overlap depending on
// `find_any`. Returns the point distance, the instance id, the shape element
// index and the element barycentric coordinates.
bool overlap_shape_bvh(const bvh_shape& bvh, const vec3f& pos, float max_distance,
    bool find_any, float& distance, int& element_id, vec2f& element_uv);
bool overlap_scene_bvh(const bvh_scene& bvh, const vec3f& pos,
    float max_distance, bool find_any, float& distance, int& instance_id,
    int& element_id, vec2f& element_uv);

}  // namespace ygl

// -----------------------------------------------------------------------------
// SHAPE EXAMPLES
// -----------------------------------------------------------------------------
namespace ygl {

// Shape data returned by make_<shape> functions.
struct make_shape_data {
    // vertex data
    vector<vec3f> positions;
    vector<vec3f> normals;
    vector<vec2f> texturecoords;
    vector<float> radius;

    // elements data
    vector<int>   points;
    vector<vec2i> lines;
    vector<vec3i> triangles;
    vector<vec4i> quads;
    vector<vec4i> beziers;

    // face-varying elements data
    vector<vec4i> quads_positions;
    vector<vec4i> quads_normals;
    vector<vec4i> quads_texturecoords;
    vector<int>   quads_materials;
};

// Make examples shapes that are not watertight (besides quads).
// Return (triangles, quads, pos, norm, texcoord)
make_shape_data make_quad_shape(const vec2i& steps, const vec2f& size,
    const vec2f& uvsize, bool as_triangles);
make_shape_data make_quad_stack_shape(const vec3i& steps, const vec3f& size,
    const vec2f& uvsize, bool as_triangles);
make_shape_data make_floor_shape(const vec2i& steps, const vec2f& size,
    const vec2f& uvsize, bool as_triangles);
make_shape_data make_cube_shape(const vec3i& steps, const vec3f& size,
    const vec3f& uvsize, bool as_triangles);
make_shape_data make_cube_rounded_shape(const vec3i& steps, const vec3f& size,
    const vec3f& uvsize, float radius, bool as_triangles);
make_shape_data make_sphere_shape(
    const vec2i& steps, float size, const vec2f& uvsize, bool as_triangles);
make_shape_data make_sphere_cube_shape(
    int steps, float size, float uvsize, bool as_triangles);
make_shape_data make_sphere_flipcap_shape(const vec2i& steps, float size,
    const vec2f& uvsize, const vec2f& zflip, bool as_triangles);
make_shape_data make_disk_shape(
    const vec2i& steps, float size, const vec2f& uvsize, bool as_triangles);
make_shape_data make_disk_quad_shape(
    int steps, float size, float uvsize, bool as_triangles);
make_shape_data make_disk_bulged_shape(
    int steps, float size, float uvsize, float height, bool as_triangles);
make_shape_data make_cylinder_side_shape(const vec2i& steps, const vec2f& size,
    const vec2f& uvsize, bool as_triangles);
make_shape_data make_cylinder_shape(const vec3i& steps, const vec2f& size,
    const vec3f& uvsize, bool as_triangles);
make_shape_data make_cylinder_rounded_shape(const vec3i& steps,
    const vec2f& size, const vec3f& uvsize, float radius, bool as_triangles);
make_shape_data make_geodesic_sphere_shape(
    int tesselation, float size, bool as_triangles);

// Make examples shapes with are watertight (good for subdivs).
// Returns (triangles, quads, pos)
make_shape_data make_suzanne_shape(float size, bool as_triangles);
make_shape_data make_cube_shape(const vec3f& size, bool as_triangles);

// Make facevarying example shapes that are watertight (good for subdivs).
make_shape_data make_cube_facevarying_shape(
    const vec3i& steps, const vec3f& size, const vec3f& uvsize);
make_shape_data make_cube_posonly_shape(
    const vec3i& steps, const vec3f& size, const vec3f& uvsize);
make_shape_data make_cube_multiplematerials_shape(
    const vec3i& steps, const vec3f& size, const vec3f& uvsize);

// Generate lines set along a quad. Returns lines, pos, norm, texcoord, radius.
make_shape_data make_lines_shape(const vec2i& steps, const vec2f& size,
    const vec2f& uvsize, const vec2f& line_radius = {0.001f, 0.001f});

// Make point primitives. Returns points, pos, norm, texcoord, radius.
make_shape_data make_point_shape(float point_radius = 0.001f);
make_shape_data make_points_shape(
    int num, float uvsize, float point_radius = 0.001f);
make_shape_data make_random_points_shape(int num, const vec3f& size,
    float uvsize, float point_radius = 0.001f, uint64_t seed = 0);

// Make a bezier circle. Returns bezier, pos.
make_shape_data make_bezier_circle_shape(
    vector<vec4i>& beziers, vector<vec3f>& pos);

// Make a hair ball around a shape.  Returns lines, pos, norm, texcoord, radius.
// length: minimum and maximum length
// rad: minimum and maximum radius from base to tip
// noise: noise added to hair (strength/scale)
// clump: clump added to hair (number/strength)
// rotation: rotation added to hair (angle/strength)
make_shape_data make_hair_shape(const vec2i& steps,
    const vector<vec3i>& striangles, const vector<vec4i>& squads,
    const vector<vec3f>& spos, const vector<vec3f>& snorm,
    const vector<vec2f>& stexcoord, const vec2f& length = {0.1f, 0.1f},
    const vec2f& rad = {0.001f, 0.001f}, const vec2f& noise = zero2f,
    const vec2f& clump = zero2f, const vec2f& rotation = zero2f, int seed = 7);

// Helper to concatenated shape data for non-facevarying shapes.
make_shape_data merge_shape_data(const vector<make_shape_data>& shapes);

}  // namespace ygl

// -----------------------------------------------------------------------------
// COLOR CONVERSION UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Element-wise float to byte conversion.
inline vec4b float_to_byte(const vec4f& a) {
    return {(byte)clamp(int(a.x * 256), 0, 255),
        (byte)clamp(int(a.y * 256), 0, 255), (byte)clamp(int(a.z * 256), 0, 255),
        (byte)clamp(int(a.w * 256), 0, 255)};
}
inline vec4f byte_to_float(const vec4b& a) {
    return {a.x / 255.0f, a.y / 255.0f, a.z / 255.0f, a.w / 255.0f};
}

// Conversion between linear and gamma-encoded colors.
inline vec3f gamma_to_linear(const vec3f& srgb, float gamma = 2.2f) {
    return {pow(srgb.x, gamma), pow(srgb.y, gamma), pow(srgb.z, gamma)};
}
inline vec3f linear_to_gamma(const vec3f& lin, float gamma = 2.2f) {
    return {pow(lin.x, 1 / gamma), pow(lin.y, 1 / gamma), pow(lin.z, 1 / gamma)};
}
inline vec4f gamma_to_linear(const vec4f& srgb, float gamma = 2.2f) {
    return {pow(srgb.x, gamma), pow(srgb.y, gamma), pow(srgb.z, gamma), srgb.w};
}
inline vec4f linear_to_gamma(const vec4f& lin, float gamma = 2.2f) {
    return {pow(lin.x, 1 / gamma), pow(lin.y, 1 / gamma), pow(lin.z, 1 / gamma),
        lin.w};
}

// sRGB non-linear curve
inline float srgb_to_linear(float srgb) {
    if (srgb <= 0.04045) {
        return srgb / 12.92f;
    } else {
        return pow((srgb + 0.055f) / (1.0f + 0.055f), 2.4f);
    }
}
inline float linear_to_srgb(float lin) {
    if (lin <= 0.0031308f) {
        return 12.92f * lin;
    } else {
        return (1 + 0.055f) * pow(lin, 1 / 2.4f) - 0.055f;
    }
}

// Conversion between linear and srgb colors.
inline vec3f srgb_to_linear(const vec3f& srgb) {
    return {
        srgb_to_linear(srgb.x), srgb_to_linear(srgb.y), srgb_to_linear(srgb.z)};
}
inline vec3f linear_to_srgb(const vec3f& lin) {
    return {linear_to_srgb(lin.x), linear_to_srgb(lin.y), linear_to_srgb(lin.z)};
}
inline vec4f srgb_to_linear(const vec4f& srgb) {
    return {srgb_to_linear(srgb.x), srgb_to_linear(srgb.y),
        srgb_to_linear(srgb.z), srgb.w};
}
inline vec4f linear_to_srgb(const vec4f& lin) {
    return {linear_to_srgb(lin.x), linear_to_srgb(lin.y), linear_to_srgb(lin.z),
        lin.w};
}

// Approximate luminance estimate for sRGB primaries (better relative luminance)
inline float luminance(const vec3f& a) {
    return (0.2126f * a.x + 0.7152f * a.y + 0.0722 * a.z);
}
inline float luminance(const vec4f& a) {
    return (0.2126f * a.x + 0.7152f * a.y + 0.0722 * a.z);
}

// Fitted ACES tonemapping curve.
inline float tonemap_filmic(float hdr) {
    // https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
    // hdr *= 0.6; // brings it back to ACES range
    return (hdr * hdr * 2.51f + hdr * 0.03f) /
           (hdr * hdr * 2.43f + hdr * 0.59f + 0.14f);
}
// Apply ACES fitted curve.
inline vec4f tonemap_filmic(const vec4f& hdr) {
    return {tonemap_filmic(hdr.x), tonemap_filmic(hdr.y), tonemap_filmic(hdr.z),
        hdr.w};
}

// Tonemap a color value according to an exposure-gamma tone mapper, with
// an optional filmic curve.
inline vec4f tonemap_filmic(
    const vec4f& hdr, float exposure, bool filmic, bool srgb) {
    auto scale = pow(2.0f, exposure);
    auto ldr   = vec4f{hdr.x * scale, hdr.y * scale, hdr.z * scale, hdr.w};
    if (filmic) ldr = tonemap_filmic(ldr);
    if (srgb) ldr = linear_to_srgb(ldr);
    return ldr;
}

// Converts HSV to RGB.
vec3f hsv_to_rgb(const vec3f& hsv);
vec3f rgb_to_hsv(const vec3f& rgb);
// Convert between CIE XYZ and xyY
vec3f xyz_to_xyY(const vec3f& xyz);
vec3f xyY_to_xyz(const vec3f& xyY);
// Convert between CIE XYZ and RGB
vec3f xyz_to_rgb(const vec3f& xyz);
vec3f rgb_to_xyz(const vec3f& rgb);

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMAGE UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Image container.
template <typename T>
struct image {
    vec2i     size   = {0, 0};
    vector<T> pixels = {};
};

// Image onstructors
template <typename T>
inline image<T> make_image(const vec2i& size, const T& v = T{}) {
    return image<T>{size, vector<T>((size_t)(size.x * size.y), v)};
}
template <typename T>
inline image<T> make_image(const vec2i& size, const T* v) {
    return image<T>{size, vector<T>(v, v + size.x * size.y)};
}

// Element access.
template <typename T>
inline T& at(image<T>& img, const vec2i& ij) {
    return img.pixels[ij.y * img.size.x + ij.x];
}
template <typename T>
inline const T& at(const image<T>& img, const vec2i& ij) {
    return img.pixels[ij.y * img.size.x + ij.x];
}

// Size
template <typename T>
inline float get_image_aspect(const image<T>& img) {
    return (float)img.size.x / (float)img.size.y;
}

// Image region defined by its corner at x,y and with size width x height
struct image_region {
    vec2i offset = {0, 0};
    vec2i size   = {0, 0};
};

// Splits an image into an array of regions
vector<image_region> make_image_regions(
    const vec2i& image_size, int region_size = 32);

// Gets pixels in an image region
template <typename T>
inline image<T> get_image_region(const image<T>& img, const image_region& region);

// Gets an image size from a suggested size and an aspect ratio. The suggested
// size may have zeros in either components. In which case, we use the aspect
// ration to compute the other.
vec2i get_image_size(const vec2i& size, float aspect);

// Conversion from/to floats.
image<vec4f> byte_to_float(const image<vec4b>& bt);
image<vec4b> float_to_byte(const image<vec4f>& fl);

// Conversion between linear and gamma-encoded images.
image<vec4f> gamma_to_linear(const image<vec4f>& srgb, float gamma);
image<vec4f> linear_to_gamma(const image<vec4f>& lin, float gamma);

// Conversion between linear and sRGB images.
image<vec4f> srgb_to_linear(const image<vec4f>& srgb);
image<vec4f> linear_to_srgb(const image<vec4f>& lin);

// Apply exposure and filmic tone mapping
image<vec4f> tonemap_image(
    const image<vec4f>& hdr, float exposure, bool filmic, bool srgb);
void tonemap_image_region(image<vec4f>& ldr, const image_region& region,
    const image<vec4f>& hdr, float exposure, bool filmic, bool srgb);

// Resize an image.
image<vec4f> resize_image(const image<vec4f>& img, const vec2i& size);

}  // namespace ygl

// -----------------------------------------------------------------------------
// EXAMPLE IMAGES
// -----------------------------------------------------------------------------
namespace ygl {

// Make example images.
image<vec4f> make_grid_image(const vec2i& size, int tile = 8,
    const vec4f& c0 = {0.2f, 0.2f, 0.2f, 1},
    const vec4f& c1 = {0.8f, 0.8f, 0.8f, 1});
image<vec4f> make_checker_image(const vec2i& size, int tile = 8,
    const vec4f& c0 = {0.2f, 0.2f, 0.2f, 1},
    const vec4f& c1 = {0.8f, 0.8f, 0.8f, 1});
image<vec4f> make_bumpdimple_image(const vec2i& size, int tile = 8);
image<vec4f> make_ramp_image(
    const vec2i& size, const vec4f& c0, const vec4f& c1, float srgb = false);
image<vec4f> make_gammaramp_image(const vec2i& size);
image<vec4f> make_uvramp_image(const vec2i& size);
image<vec4f> make_uvgrid_image(
    const vec2i& size, int tile = 8, bool colored = true);

// Comvert a bump map to a normal map.
image<vec4f> bump_to_normal_map(const image<vec4f>& img, float scale = 1);

// Make a sunsky HDR model with sun at theta elevation in [0,pif/2], turbidity
// in [1.7,10] with or without sun.
image<vec4f> make_sunsky_image(const vec2i& size, float thetaSun,
    float turbidity = 3, bool has_sun = false,
    const vec3f& ground_albedo = {0.7f, 0.7f, 0.7f});
// Make an image of multiple lights.
image<vec4f> make_lights_image(const vec2i& size, const vec3f& le = {1, 1, 1},
    int nlights = 4, float langle = pif / 4, float lwidth = pif / 16,
    float lheight = pif / 16);

// Make a noise image. Wrap works only if both resx and resy are powers of two.
image<vec4f> make_noise_image(
    const vec2i& size, float scale = 1, bool wrap = true);
image<vec4f> make_fbm_image(const vec2i& size, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, int octaves = 6, bool wrap = true);
image<vec4f> make_ridge_image(const vec2i& size, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, float offset = 1.0f,
    int octaves = 6, bool wrap = true);
image<vec4f> make_turbulence_image(const vec2i& size, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, int octaves = 6, bool wrap = true);

}  // namespace ygl

// -----------------------------------------------------------------------------
// ANIMATION UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Find the first keyframe value that is greater than the argument.
inline int evaluate_keyframed_index(
    const vector<float>& times, const float& time) {
    for (auto i = 0; i < times.size(); i++)
        if (times[i] > time) return i;
    return (int)times.size();
}

// Evaluates a keyframed value using step interpolation.
template <typename T>
inline T evaluate_keyframed_step(
    const vector<float>& times, const vector<T>& vals, float time) {
    if (time <= times.front()) return vals.front();
    if (time >= times.back()) return vals.back();
    time     = clamp(time, times.front(), times.back() - 0.001f);
    auto idx = evaluate_keyframed_index(times, time);
    return vals.at(idx - 1);
}

// Evaluates a keyframed value using linear interpolation.
template <typename T>
inline vec4f evaluate_keyframed_slerp(
    const vector<float>& times, const vector<vec4f>& vals, float time) {
    if (time <= times.front()) return vals.front();
    if (time >= times.back()) return vals.back();
    time     = clamp(time, times.front(), times.back() - 0.001f);
    auto idx = evaluate_keyframed_index(times, time);
    auto t   = (time - times.at(idx - 1)) / (times.at(idx) - times.at(idx - 1));
    return slerp(vals.at(idx - 1), vals.at(idx), t);
}

// Evaluates a keyframed value using linear interpolation.
template <typename T>
inline T evaluate_keyframed_linear(
    const vector<float>& times, const vector<T>& vals, float time) {
    if (time <= times.front()) return vals.front();
    if (time >= times.back()) return vals.back();
    time     = clamp(time, times.front(), times.back() - 0.001f);
    auto idx = evaluate_keyframed_index(times, time);
    auto t   = (time - times.at(idx - 1)) / (times.at(idx) - times.at(idx - 1));
    return vals.at(idx - 1) * (1 - t) + vals.at(idx) * t;
}

// Evaluates a keyframed value using Bezier interpolation.
template <typename T>
inline T evaluate_keyframed_bezier(
    const vector<float>& times, const vector<T>& vals, float time) {
    if (time <= times.front()) return vals.front();
    if (time >= times.back()) return vals.back();
    time     = clamp(time, times.front(), times.back() - 0.001f);
    auto idx = evaluate_keyframed_index(times, time);
    auto t   = (time - times.at(idx - 1)) / (times.at(idx) - times.at(idx - 1));
    return interpolate_bezier(
        vals.at(idx - 3), vals.at(idx - 2), vals.at(idx - 1), vals.at(idx), t);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// VOLUME TYPE AND UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Volume container.
template <typename T>
struct volume {
    vec3i     size   = {0, 0, 0};
    vector<T> voxels = {};
};

// Volume onstructors
template <typename T>
inline volume<T> make_volume(const vec3i& size, const T& v = T{}) {
    return volume<T>{size, vector<T>((size_t)(size.x * size.y * size.z), v)};
}
template <typename T>
inline volume<T> make_volume(const vec3i& size, const T* v) {
    return volume<T>{size, vector<T>(v, v + size.x * size.y * size.z)};
}

// Element access
template <typename T>
T& at(volume<T>& vol, const vec3i& ijk) {
    return vol.voxels[ijk.z * vol.size.x * vol.size.y + ijk.y * vol.size.x + ijk.x];
}
template <typename T>
const T& at(const volume<T>& vol, const vec3i& ijk) {
    return vol.voxels[ijk.z * vol.size.x * vol.size.y + ijk.y * vol.size.x + ijk.x];
}

// make a simple example volume
volume<float> make_test_volume1f(
    const vec3i& size, float scale = 10, float exponent = 6);

}  // namespace ygl

// -----------------------------------------------------------------------------
// SCENE DATA
// -----------------------------------------------------------------------------
namespace ygl {

// Camera based on a simple lens model. The camera is placed using a frame.
// Camera projection is described in photorgaphics terms. In particular,
// we specify fil size (35mm by default), the focal lengthm the focus
// distance and the lens_aperture. All values are in meters.
struct yocto_camera {
    string  name           = "";
    frame3f frame          = identity_frame3f;
    bool    orthographic   = false;
    vec2f   film_size      = {0.036f, 0.024f};
    float   focal_length   = 0.050f;
    float   focus_distance = maxf;
    float   lens_aperture  = 0;
};

// Texture containing either an LDR or HDR image. Textures are rendered
// using linear interpolation (unless `no_interoilation` is set) and
// weith tiling (unless `clamp_to_edge` is set). HdR images are encoded
// in linear color space, while LDRs are encoded as sRGB. The latter
// conversion can be disabled with `ldr_as_linear` for example to render
// normal maps.
struct yocto_texture {
    string       name             = "";
    string       filename         = "";
    image<vec4f> hdr_image        = {};
    image<vec4b> ldr_image        = {};
    bool         clamp_to_edge    = false;
    bool         no_interpolation = false;
    float        height_scale     = 1;
    bool         ldr_as_linear    = false;
    bool         has_opacity      = false;
};

// Volumetric texture containing a float only volume data. See texture
// above for other propoerties.
struct yocto_voltexture {
    string        name             = "";
    string        filename         = "";
    volume<float> volume_data      = {};
    bool          clamp_to_edge    = false;
    bool          no_interpolation = false;
};

// Material for surfaces, lines and triangles.
// For surfaces, uses a microfacet model with thin sheet transmission.
// The model is based on OBJ, but contains glTF compatibility.
// For the documentation on the values, please see the OBJ format.
struct yocto_material {
    string name          = "";
    bool   base_metallic = false;  // base-metallic parametrization
    bool   gltf_textures = false;  // glTF packed textures
    bool   double_sided  = false;  // double sided rendering

    // base values
    vec3f emission     = {0, 0, 0};
    vec3f diffuse      = {0, 0, 0};
    vec3f specular     = {0, 0, 0};
    vec3f transmission = {0, 0, 0};
    float roughness    = 0.0001;
    float opacity      = 1;
    bool  fresnel      = true;
    bool  refract      = false;

    // textures
    int emission_texture     = -1;
    int diffuse_texture      = -1;
    int specular_texture     = -1;
    int transmission_texture = -1;
    int roughness_texture    = -1;
    int opacity_texture      = -1;
    int occlusion_texture    = -1;
    int bump_texture         = -1;
    int displacement_texture = -1;
    int normal_texture       = -1;

    // volume properties
    // albedo = scattering / (absorption + scattering)
    // density = absorption + scattering
    vec3f volume_emission = {0, 0, 0};
    vec3f volume_albedo   = {0, 0, 0};
    vec3f volume_density  = {0, 0, 0};
    float volume_phaseg   = 0;

    // volume textures
    int volume_density_texture = -1;
};

// Shape data represented as an indexed meshes of elements.
// May contain either points, lines, triangles and quads.
struct yocto_shape {
    // shape data
    string name     = "";
    string filename = "";
    int    material = -1;

    // subdision properties
    int  subdivision_level      = 0;
    bool catmull_clark          = false;
    bool compute_vertex_normals = false;

    // primitives
    vector<int>   points    = {};
    vector<vec2i> lines     = {};
    vector<vec3i> triangles = {};
    vector<vec4i> quads     = {};

    // vertex data
    vector<vec3f> positions     = {};
    vector<vec3f> normals       = {};
    vector<vec2f> texturecoords = {};
    vector<vec4f> colors        = {};
    vector<float> radius        = {};
    vector<vec4f> tangentspaces = {};
};

// Shape data represented as an indexed meshes of elements with face-varying
// data. Each face maintains different topologies for positions, normals and
// texture coordinates.
struct yocto_surface {
    // shape data
    string      name      = "";
    string      filename  = "";
    vector<int> materials = {};

    // subdision properties
    int  subdivision_level      = 0;
    bool catmull_clark          = false;
    bool compute_vertex_normals = false;

    // face-varying primitives
    vector<vec4i> quads_positions     = {};
    vector<vec4i> quads_normals       = {};
    vector<vec4i> quads_texturecoords = {};
    vector<int>   quads_materials     = {};

    // vertex data
    vector<vec3f> positions     = {};
    vector<vec3f> normals       = {};
    vector<vec2f> texturecoords = {};
};

// Instance of a visible object in the scene. For now, this can be either
// a shape or a surface.
struct yocto_instance {
    string  name    = "";
    frame3f frame   = identity_frame3f;
    int     shape   = -1;
    int     surface = -1;
};

// Environment map.
struct yocto_environment {
    string  name             = "";
    frame3f frame            = identity_frame3f;
    vec3f   emission         = {0, 0, 0};
    int     emission_texture = -1;
};

// Node in a transform hierarchy.
struct yocto_scene_node {
    string        name        = "";
    int           parent      = -1;
    frame3f       local       = identity_frame3f;
    vec3f         translation = {0, 0, 0};
    vec4f         rotation    = {0, 0, 0, 1};
    vec3f         scale       = {1, 1, 1};
    vector<float> weights     = {};
    int           camera      = -1;
    int           instance    = -1;
    int           environment = -1;

    // compute properties
    vector<int> children = {};
};

// Keyframe type.
enum struct yocto_interpolation_type { linear, step, bezier };

// Keyframe data.
struct yocto_animation {
    string                   name               = "";
    string                   filename           = "";
    string                   animation_group    = "";
    yocto_interpolation_type interpolation_type = yocto_interpolation_type::linear;
    vector<float>            keyframes_times    = {};
    vector<vec3f>            translation_keyframes   = {};
    vector<vec4f>            rotation_keyframes      = {};
    vector<vec3f>            scale_keyframes         = {};
    vector<vector<float>>    morph_weights_keyframes = {};
    vector<int>              node_targets            = {};
};

// Scene comprised an array of objects whose memory is owened by the scene.
// All members are optional,Scene objects (camera, instances, environments)
// have transforms defined internally. A scene can optionally contain a
// node hierarchy where each node might point to a camera, instance or
// environment. In that case, the element transforms are computed from
// the hierarchy. Animation is also optional, with keyframe data that
// updates node transformations only if defined.
struct yocto_scene {
    string                    name         = "";
    vector<yocto_camera>      cameras      = {};
    vector<yocto_shape>       shapes       = {};
    vector<yocto_surface>     surfaces     = {};
    vector<yocto_instance>    instances    = {};
    vector<yocto_material>    materials    = {};
    vector<yocto_texture>     textures     = {};
    vector<yocto_environment> environments = {};
    vector<yocto_voltexture>  voltextures  = {};
    vector<yocto_scene_node>  nodes        = {};
    vector<yocto_animation>   animations   = {};
};

}  // namespace ygl

// -----------------------------------------------------------------------------
// SCENE UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Print scene statistics.
void print_stats(const yocto_scene& scene);

// Merge scene into one another. Note that the objects are _moved_ from
// merge_from to merged_into, so merge_from will be empty after this function.
void merge_scene(yocto_scene& merge_into, yocto_scene& merge_from);

}  // namespace ygl

// -----------------------------------------------------------------------------
// EVALUATION OF SCENE PROPERTIES
// -----------------------------------------------------------------------------
namespace ygl {

// Update node transforms.
void update_transforms(
    yocto_scene& scene, float time = 0, const string& anim_group = "");

// Compute animation range.
vec2f compute_animation_range(
    const yocto_scene& scene, const string& anim_group = "");

// Computes shape/scene approximate bounds.
bbox3f compute_shape_bounds(const yocto_shape& shape);
bbox3f compute_scene_bounds(const yocto_scene& scene);

// Compute shape vertex normals
vector<vec3f> compute_shape_normals(const yocto_shape& shape);

// Updates/refits bvh.
bvh_shape make_shape_bvh(
    const yocto_shape& shape, bool high_quality, bool embree = false);
bvh_shape make_surface_bvh(
    const yocto_surface& surface, bool high_quality, bool embree = false);
bvh_scene make_scene_bvh(
    const yocto_scene& scene, bool high_quality, bool embree = false);
void refit_shape_bvh(const yocto_shape& shape, bvh_shape& bvh);
void refit_scene_bvh(const yocto_scene& scene, bvh_scene& bvh);

// Apply subdivision and displacement rules.
void tesselate_shapes_and_surfaces(yocto_scene& scene);

// Add missing names, normals, tangents and hierarchy.
void add_missing_names(yocto_scene& scene);
void add_missing_normals(yocto_scene& scene);
void add_missing_tangent_space(yocto_scene& scene);
void add_missing_materials(yocto_scene& scene);
void add_missing_cameras(yocto_scene& scene);

// Add a sky environment
void add_sky_environment(yocto_scene& scene, float sun_angle = pif / 4);

// Checks for validity of the scene.
vector<string> validate_scene(
    const yocto_scene& scene, bool skip_textures = false);
void log_validation_errors(const yocto_scene& scene, bool skip_textures = false);

// Queries on objects
bool is_shape_face_varying(const yocto_shape& shape);

// Scene intersection. Upron intersection we set the instance pointer,
// the shape element_id and element_uv and the inetrsection distance.
struct scene_intersection {
    int   instance_id = -1;
    int   element_id  = -1;
    vec2f element_uv  = zero2f;
    float distance    = maxf;
};

// Intersects a ray with the scene.
scene_intersection intersect_scene(const yocto_scene& scene,
    const bvh_scene& bvh, const ray3f& ray, bool find_any = false);
// Intersects a ray with a scene instance.
scene_intersection intersect_scene(const yocto_scene& scene, int instance_id,
    const bvh_scene& bvh, const ray3f& ray, bool find_any = false);

// Shape values interpolated using barycentric coordinates.
vec3f evaluate_shape_position(
    const yocto_shape& shape, int element_id, const vec2f& element_uv);
vec3f evaluate_shape_normal(
    const yocto_shape& shape, int element_id, const vec2f& element_uv);
vec2f evaluate_shape_texturecoord(
    const yocto_shape& shape, int element_id, const vec2f& element_uv);
vec4f evaluate_shape_color(
    const yocto_shape& shape, int element_id, const vec2f& element_uv);
float evaluate_shape_radius(
    const yocto_shape& shape, int element_id, const vec2f& element_uv);
vec4f evaluate_shape_tangentspace(
    const yocto_shape& shape, int element_id, const vec2f& element_uv);
vec3f evaluate_shape_tangentspace(const yocto_shape& shape, int element_id,
    const vec2f& element_uv, bool& left_handed);
// Shape element values.
vec3f evaluate_shape_element_normal(const yocto_shape& shape, int element_id);
vec4f evaluate_shape_element_tangentspace(
    const yocto_shape& shape, int element_id);

// Sample a shape element based on area/length.
vector<float>     compute_shape_elements_cdf(const yocto_shape& shape);
tuple<int, vec2f> sample_shape_element(const yocto_shape& shape,
    const vector<float>& elem_cdf, float re, const vec2f& ruv);
float             sample_shape_element_pdf(const yocto_shape& shape,
                const vector<float>& elem_cdf, int element_id, const vec2f& element_uv);

// Surface values interpolated using barycentric coordinates.
vec3f evaluate_surface_position(
    const yocto_surface& surface, int element_id, const vec2f& element_uv);
vec3f evaluate_surface_normal(
    const yocto_surface& surface, int element_id, const vec2f& element_uv);
vec2f evaluate_surface_texturecoord(
    const yocto_surface& surface, int element_id, const vec2f& element_uv);
// Surface element values.
vec3f evaluate_surface_element_normal(const yocto_surface& shape, int element_id);
// Per-element material.
int get_surface_element_material(const yocto_surface& surface, int element_id);

// Sample a surface element based on area.
vector<float>     compute_surface_elements_cdf(const yocto_surface& surface);
tuple<int, vec2f> sample_surface_element(const yocto_surface& surface,
    const vector<float>& elem_cdf, float re, const vec2f& ruv);
float             sample_surface_element_pdf(const yocto_surface& surface,
                const vector<float>& elem_cdf, int element_id, const vec2f& element_uv);

// Evaluate a texture.
vec2i evaluate_texture_size(const yocto_texture& texture);
vec4f lookup_texture(const yocto_texture& texture, const vec2i& ij);
vec4f evaluate_texture(const yocto_texture& texture, const vec2f& texcoord);
float lookup_voltexture(const yocto_voltexture& texture, const vec3i& ijk);
float evaluate_voltexture(const yocto_voltexture& texture, const vec3f& texcoord);

// Set and evaluate camera parameters. Setters take zeros as default values.
float get_camera_fovx(const yocto_camera& camera);
float get_camera_fovy(const yocto_camera& camera);
float get_camera_aspect(const yocto_camera& camera);
vec2i get_camera_image_size(const yocto_camera& camera, const vec2i& size);
void  set_camera_fovy(
     yocto_camera& camera, float fovy, float aspect, float width = 0.036f);
// Sets camera field of view to enclose all the bbox. Camera view direction
// fiom size and forcal lemgth can be overridden if we pass non zero values.
void set_camera_view(yocto_camera& camera, const bbox3f& bbox,
    const vec3f& view_direction = zero3f, const vec2f& film = zero2f,
    float focal = 0);

// Generates a ray from a camera image coordinate and lens coordinates.
ray3f evaluate_camera_ray(
    const yocto_camera& camera, const vec2f& image_uv, const vec2f& lens_uv);
// Generates a ray from a camera for pixel `image_ij`, the image size,
// the sub-pixel coordinates `pixel_uv` and the lens coordinates `lens_uv`
// and the image resolution `image_size`.
ray3f evaluate_camera_ray(const yocto_camera& camera, const vec2i& image_ij,
    const vec2i& image_size, const vec2f& pixel_uv, const vec2f& lens_uv);
// Generates a ray from a camera for pixel index `idx`, the image size,
// the sub-pixel coordinates `pixel_uv` and the lens coordinates `lens_uv`.
ray3f evaluate_camera_ray(const yocto_camera& camera, int idx,
    const vec2i& image_size, const vec2f& pixel_uv, const vec2f& lens_uv);

// Evaluates material parameters: emission, diffuse, specular, transmission,
// roughness and opacity.
vec3f evaluate_material_emission(const yocto_scene& scene,
    const yocto_material& material, const vec2f& texturecoord,
    const vec4f& shape_color = {1, 1, 1, 1});
vec3f evaluate_material_diffuse(const yocto_scene& scene,
    const yocto_material& material, const vec2f& texturecoord,
    const vec4f& shape_color = {1, 1, 1, 1});
vec3f evaluate_material_specular(const yocto_scene& scene,
    const yocto_material& material, const vec2f& texturecoord,
    const vec4f& shape_color = {1, 1, 1, 1});
vec3f evaluate_material_transmission(const yocto_scene& scene,
    const yocto_material& material, const vec2f& texturecoord,
    const vec4f& shape_color = {1, 1, 1, 1});
float evaluate_material_roughness(const yocto_scene& scene,
    const yocto_material& material, const vec2f& texturecoord,
    const vec4f& shape_color = {1, 1, 1, 1});
float evaluate_material_opacity(const yocto_scene& scene,
    const yocto_material& material, const vec2f& texturecoord,
    const vec4f& shape_color = {1, 1, 1, 1});
// Query material properties
bool is_material_emissive(const yocto_material& material);

// Material values packed into a convenience structure.
struct microfacet_brdf {
    vec3f diffuse      = zero3f;
    vec3f specular     = zero3f;
    vec3f transmission = zero3f;
    float roughness    = 1;
    bool  refract      = false;
};
microfacet_brdf evaluate_material_brdf(const yocto_scene& scene,
    const yocto_material& material, const vec2f& texturecoord,
    const vec4f& shape_color = {1, 1, 1, 1});
bool            is_bsdf_delta(const microfacet_brdf& f);

// Check volume properties.
bool is_material_volume_homogeneus(const yocto_material& vol);
bool is_material_volume_colored(const yocto_material& vol);

// Instance values interpolated using barycentric coordinates.
// Handles defaults if data is missing.
vec3f evaluate_instance_position(const yocto_scene& scene,
    const yocto_instance& instance, int element_id, const vec2f& element_uv);
vec3f evaluate_instance_normal(const yocto_scene& scene,
    const yocto_instance& instance, int element_id, const vec2f& element_uv);
vec2f evaluate_instance_texturecoord(const yocto_scene& scene,
    const yocto_instance& instance, int element_id, const vec2f& element_uv);
vec4f evaluate_instance_color(
    const yocto_instance& instance, int element_id, const vec2f& element_uv);
vec3f evaluate_instance_tangentspace(const yocto_scene& scene,
    const yocto_instance& instance, int element_id, const vec2f& element_uv,
    bool& left_handed);
// Instance element values.
vec3f evaluate_instance_element_normal(
    const yocto_scene& scene, const yocto_instance& instance, int element_id);
// Shading normals including material perturbations.
vec3f evaluate_instance_shading_normal(const yocto_scene& scene,
    const yocto_instance& instance, int element_id, const vec2f& element_uv,
    const vec3f& o);

// Material values
int   get_instance_material_id(const yocto_scene& scene,
      const yocto_instance& instance, int element_id, const vec2f& element_uv);
vec3f evaluate_instance_emission(const yocto_scene& scene,
    const yocto_instance& instance, int element_id, const vec2f& element_uv);
float evaluate_instance_opacity(const yocto_scene& scene,
    const yocto_instance& instance, int element_id, const vec2f& element_uv);
bool  is_instance_emissive(
     const yocto_scene& scene, const yocto_instance& instance);

// <aterial brdf
microfacet_brdf evaluate_instance_brdf(const yocto_scene& scene,
    const yocto_instance& instance, int element_id, const vec2f& element_uv);

// Environment texture coordinates from the incoming direction.
vec2f evaluate_environment_texturecoord(
    const yocto_environment& environment, const vec3f& direction);
// Evaluate the incoming direction from the element_uv.
vec3f evaluate_environment_direction(
    const yocto_environment& environment, const vec2f& environment_uv);
// Evaluate the environment emission.
vec3f evaluate_environment_emission(const yocto_scene& scene,
    const yocto_environment& environment, const vec3f& direction);
// Evaluate all environment emission.
vec3f evaluate_environment_emission(
    const yocto_scene& scene, const vec3f& direction);

// Sample an environment based on either texel values of uniform
vector<float> compute_environment_texels_cdf(
    const yocto_scene& scene, const yocto_environment& environment);
vec3f sample_environment_direction(const yocto_scene& scene,
    const yocto_environment& environment, const vector<float>& texels_cdf,
    float re, const vec2f& ruv);
float sample_environment_direction_pdf(const yocto_scene& scene,
    const yocto_environment& environment, const vector<float>& texels_cdf,
    const vec3f& direction);

}  // namespace ygl

// -----------------------------------------------------------------------------
// PATH TRACING
// -----------------------------------------------------------------------------
namespace ygl {

// Default trace seed
const auto trace_default_seed = 961748941ull;

// Trace lights used during rendering.
struct trace_lights {
    vector<int>           instances               = {};
    vector<int>           environments            = {};
    vector<vector<float>> shape_elements_cdf      = {};
    vector<vector<float>> surface_elements_cdf    = {};
    vector<vector<float>> environment_texture_cdf = {};
};

// Initialize lights.
trace_lights make_trace_lights(const yocto_scene& scene);
inline bool  empty(const trace_lights& lights) {
    return lights.instances.empty() && lights.environments.empty();
}

// Initialize state of the renderer.
image<rng_state> make_trace_rngs(
    const vec2i& image_size, uint64_t random_seed = trace_default_seed);

// Type of tracing algorithm to use
enum struct trace_sampler_type {
    path,               // path tracing
    direct,             // direct illumination
    environment,        // environment illumination only
    eyelight,           // eyelight rendering
    path_nomis,         // path tracer without mis
    path_naive,         // naive path tracing
    direct_nomis,       // direct illumition without mis
    debug_normal,       // debug - normal
    debug_albedo,       // debug - albedo
    debug_texcoord,     // debug - texcoord
    debug_frontfacing,  // debug - faceforward
    debug_diffuse,      // debug - diffuse
    debug_specular,     // debug - specular
    debug_roughness,    // debug - roughness
};

const auto trace_sampler_type_names = vector<string>{"path", "direct",
    "environment", "eyelight", "path_nomis", "path_naive", "direct_nomis",
    "debug_normal", "debug_albedo", "debug_texcoord", "debug_frontfacing",
    "debug_diffuse", "debug_specular", "debug_roughness"};

// Tracer function
using trace_sampler_func = function<vec4f(const yocto_scene& scene,
    const bvh_scene& bvh, const trace_lights& lights, const vec3f& position,
    const vec3f& direction, rng_state& rng, int max_bounces)>;
trace_sampler_func get_trace_sampler_func(trace_sampler_type type);

// Progressively compute an image by calling trace_samples multiple times.
void trace_image(image<vec4f>& rendered_image, const yocto_scene& scene,
    const yocto_camera& camera, const bvh_scene& bvh, const trace_lights& lights,
    const trace_sampler_func& trace_sampler, int num_samples, int max_bounces,
    float pixel_clamp = 100, bool no_parallel = false);

// Progressively compute an image by calling trace_samples multiple times.
// Start with an empty state and then successively call this function to
// render the next batch of samples.
void trace_samples(image<vec4f>& rendered_image, const yocto_scene& scene,
    const yocto_camera& camera, const bvh_scene& bvh,
    const trace_lights& lights, const trace_sampler_func& trace_sampler,
    int current_sample, int num_samples, int max_bounces,
    image<rng_state>& rngs, float pixel_clamp = 100, bool no_parallel = false);

// Starts an anyncrhounous renderer. The function will keep a reference to
// params.
void trace_async_start(image<vec4f>& rendered_image, const yocto_scene& scene,
    const yocto_camera& camera, const bvh_scene& bvh,
    const trace_lights& lights, const trace_sampler_func& trace_sampler,
    int num_samples, int max_bounces, image<rng_state>& rngs,
    vector<thread>& threads, bool& stop_flag, int& current_sample,
    concurrent_queue<image_region>& queue, float pixel_clamp = 100);
// Stop the asynchronous renderer.
void trace_async_stop(vector<thread>& threads, bool& stop_flag,
    concurrent_queue<image_region>& queue);

// Trace statistics for last run used for fine tuning implementation.
// For now returns number of paths and number of rays.
tuple<uint64_t, uint64_t> get_trace_stats();
void                      reset_trace_stats();

}  // namespace ygl

// -----------------------------------------------------------------------------
// PATH TRACING SUPPORT FUNCTION
// -----------------------------------------------------------------------------
namespace ygl {

// Phong exponent to roughness.
float specular_exponent_to_roughness(float n);

// Specular to fresnel eta.
void  specular_fresnel_from_ks(const vec3f& ks, vec3f& es, vec3f& esk);
float specular_to_eta(const vec3f& ks);
// Compute the fresnel term for dielectrics.
vec3f fresnel_dielectric(float cosw, const vec3f& eta_);
// Compute the fresnel term for metals.
vec3f fresnel_metal(float cosw, const vec3f& eta, const vec3f& etak);
// Schlick approximation of Fresnel term, optionally weighted by rs;
vec3f fresnel_schlick(const vec3f& ks, float cosw);
vec3f fresnel_schlick(const vec3f& ks, float cosw, float rs);
vec3f fresnel_schlick(const vec3f& ks, const vec3f& h, const vec3f& o);
vec3f fresnel_schlick(const vec3f& ks, const vec3f& h, const vec3f& o, float rs);

// Evaluates the GGX distribution and geometric term.
float evaluate_ggx(float rs, float ndh, float ndi, float ndo);
// Sample the GGX distribution.
vec3f sample_ggx(float rs, const vec2f& rn);
float sample_ggx_pdf(float rs, float ndh);

// Evaluates the GGX distribution and geometric term.
float evaluate_ggx_dist(float rs, const vec3f& n, const vec3f& h);
float evaluate_ggx_sm(float rs, const vec3f& n, const vec3f& o, const vec3f& i);

// Evaluate and sample volume phase function.
vec3f sample_phase_function(float vg, const vec2f& u);
float evaluate_phase_function(float cos_theta, float vg);

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR MATRICES
// -----------------------------------------------------------------------------
namespace ygl {

// Matrix diagonals and transposes.
template <typename T>
constexpr inline mat<T, 1, 1> transpose(const mat<T, 1, 1>& a) {
    return {{a.x}};
}
template <typename T>
constexpr inline mat<T, 2, 2> transpose(const mat<T, 2, 2>& a) {
    return {{a.x.x, a.y.x}, {a.x.y, a.y.y}};
}
template <typename T>
constexpr inline mat<T, 3, 3> transpose(const mat<T, 3, 3>& a) {
    return {
        {a.x.x, a.y.x, a.z.x},
        {a.x.y, a.y.y, a.z.y},
        {a.x.z, a.y.z, a.z.z},
    };
}
template <typename T>
constexpr inline mat<T, 4, 4> transpose(const mat<T, 4, 4>& a) {
    return {
        {a.x.x, a.y.x, a.z.x, a.w.x},
        {a.x.y, a.y.y, a.z.y, a.w.y},
        {a.x.z, a.y.z, a.z.z, a.w.z},
        {a.x.w, a.y.w, a.z.w, a.w.w},
    };
}

// Matrix adjugates, determinant and inverses.
template <typename T>
constexpr inline mat<T, 1, 1> adjugate(const mat<T, 1, 1>& a) {
    return {{a.x}};
}
template <typename T>
constexpr inline mat<T, 2, 2> adjugate(const mat<T, 2, 2>& a) {
    return {{a.y.y, -a.x.y}, {-a.y.x, a.x.x}};
}
template <typename T>
constexpr inline mat<T, 3, 3> adjugate(const mat<T, 3, 3>& a) {
    return {
        {
            a.y.y * a.z.z - a.z.y * a.y.z,
            a.z.y * a.x.z - a.x.y * a.z.z,
            a.x.y * a.y.z - a.y.y * a.x.z,
        },
        {
            a.y.z * a.z.x - a.z.z * a.y.x,
            a.z.z * a.x.x - a.x.z * a.z.x,
            a.x.z * a.y.x - a.y.z * a.x.x,
        },
        {
            a.y.x * a.z.y - a.z.x * a.y.y,
            a.z.x * a.x.y - a.x.x * a.z.y,
            a.x.x * a.y.y - a.y.x * a.x.y,
        },
    };
}
template <typename T>
constexpr inline mat<T, 4, 4> adjugate(const mat<T, 4, 4>& a) {
    return {
        {
            a.y.y * a.z.z * a.w.w + a.w.y * a.y.z * a.z.w +
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
                a.z.y * a.x.z * a.y.w - a.y.y * a.z.z * a.x.w,
        },
        {
            a.y.z * a.w.w * a.z.x + a.z.z * a.y.w * a.w.x +
                a.w.z * a.z.w * a.y.x - a.y.z * a.z.w * a.w.x -
                a.w.z * a.y.w * a.z.x - a.z.z * a.w.w * a.y.x,
            a.x.z * a.z.w * a.w.x + a.w.z * a.x.w * a.z.x +
                a.z.z * a.w.w * a.x.x - a.x.z * a.w.w * a.z.x -
                a.z.z * a.x.w * a.w.x - a.w.z * a.z.w * a.x.x,
            a.x.z * a.w.w * a.y.x + a.y.z * a.x.w * a.w.x +
                a.w.z * a.y.w * a.x.x - a.x.z * a.y.w * a.w.x -
                a.w.z * a.x.w * a.y.x - a.y.z * a.w.w * a.x.x,
            a.x.z * a.y.w * a.z.x + a.z.z * a.x.w * a.y.x +
                a.y.z * a.z.w * a.x.x - a.x.z * a.z.w * a.y.x -
                a.y.z * a.x.w * a.z.x - a.z.z * a.y.w * a.x.x,
        },
        {
            a.y.w * a.z.x * a.w.y + a.w.w * a.y.x * a.z.y +
                a.z.w * a.w.x * a.y.y - a.y.w * a.w.x * a.z.y -
                a.z.w * a.y.x * a.w.y - a.w.w * a.z.x * a.y.y,
            a.x.w * a.w.x * a.z.y + a.z.w * a.x.x * a.w.y +
                a.w.w * a.z.x * a.x.y - a.x.w * a.z.x * a.w.y -
                a.w.w * a.x.x * a.z.y - a.z.w * a.w.x * a.x.y,
            a.x.w * a.y.x * a.w.y + a.w.w * a.x.x * a.y.y +
                a.y.w * a.w.x * a.x.y - a.x.w * a.w.x * a.y.y -
                a.y.w * a.x.x * a.w.y - a.w.w * a.y.x * a.x.y,
            a.x.w * a.z.x * a.y.y + a.y.w * a.x.x * a.z.y +
                a.z.w * a.y.x * a.x.y - a.x.w * a.y.x * a.z.y -
                a.z.w * a.x.x * a.y.y - a.y.w * a.z.x * a.x.y,
        },
        {
            a.y.x * a.w.y * a.z.z + a.z.x * a.y.y * a.w.z +
                a.w.x * a.z.y * a.y.z - a.y.x * a.z.y * a.w.z -
                a.w.x * a.y.y * a.z.z - a.z.x * a.w.y * a.y.z,
            a.x.x * a.z.y * a.w.z + a.w.x * a.x.y * a.z.z +
                a.z.x * a.w.y * a.x.z - a.x.x * a.w.y * a.z.z -
                a.z.x * a.x.y * a.w.z - a.w.x * a.z.y * a.x.z,
            a.x.x * a.w.y * a.y.z + a.y.x * a.x.y * a.w.z +
                a.w.x * a.y.y * a.x.z - a.x.x * a.y.y * a.w.z -
                a.w.x * a.x.y * a.y.z - a.y.x * a.w.y * a.x.z,
            a.x.x * a.y.y * a.z.z + a.z.x * a.x.y * a.y.z +
                a.y.x * a.z.y * a.x.z - a.x.x * a.z.y * a.y.z -
                a.y.x * a.x.y * a.z.z - a.z.x * a.y.y * a.x.z,
        },
    };
}
template <typename T>
constexpr inline T determinant(const mat<T, 1, 1>& a) {
    return a.x;
}
template <typename T>
constexpr inline T determinant(const mat<T, 2, 2>& a) {
    return a.x.x * a.y.y - a.x.y * a.y.x;
}
template <typename T>
constexpr inline T determinant(const mat<T, 3, 3>& a) {
    return a.x.x * (a.y.y * a.z.z - a.z.y * a.y.z) +
           a.x.y * (a.y.z * a.z.x - a.z.z * a.y.x) +
           a.x.z * (a.y.x * a.z.y - a.z.x * a.y.y);
}
template <typename T>
constexpr inline T determinant(const mat<T, 4, 4>& a) {
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
// IMPLEMENTATION FOR UI UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Turntable for UI navigation.
inline void camera_turntable(vec3f& from, vec3f& to, vec3f& up,
    const vec2f& rotate, float dolly, const vec2f& pan) {
    // rotate if necessary
    if (rotate.x || rotate.y) {
        auto z     = normalize(to - from);
        auto lz    = length(to - from);
        auto phi   = atan2(z.z, z.x) + rotate.x;
        auto theta = acos(z.y) + rotate.y;
        theta      = clamp(theta, 0.001f, pif - 0.001f);
        auto nz    = vec3f{sin(theta) * cos(phi) * lz, cos(theta) * lz,
            sin(theta) * sin(phi) * lz};
        from       = to - nz;
    }

    // dolly if necessary
    if (dolly) {
        auto z  = normalize(to - from);
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

// Turntable for UI navigation.
inline void camera_turntable(frame3f& frame, float& focus, const vec2f& rotate,
    float dolly, const vec2f& pan) {
    // rotate if necessary
    if (rotate != zero2f) {
        auto phi   = atan2(frame.z.z, frame.z.x) + rotate.x;
        auto theta = acos(frame.z.y) + rotate.y;
        theta      = clamp(theta, 0.001f, pif - 0.001f);
        auto new_z = vec3f{
            sin(theta) * cos(phi), cos(theta), sin(theta) * sin(phi)};
        auto new_center = frame.o - frame.z * focus;
        auto new_o      = new_center + new_z * focus;
        frame           = lookat_frame(new_o, new_center, {0, 1, 0});
        focus           = length(new_o - new_center);
    }

    // pan if necessary
    if (dolly) {
        auto c  = frame.o - frame.z * focus;
        focus   = max(focus * (1 + dolly), 0.001f);
        frame.o = c + frame.z * focus;
    }

    // pan if necessary
    if (pan.x || pan.y) {
        frame.o += frame.x * pan.x + frame.y * pan.y;
    }
}

// FPS camera for UI navigation for a frame parametrization.
inline void camera_fps(frame3f& frame, vec3f transl, vec2f rotate) {
    // https://gamedev.stackexchange.com/questions/30644/how-to-keep-my-quaternion-using-fps-camera-from-tilting-and-messing-up
    auto y = vec3f{0, 1, 0};
    auto z = orthonormalize(frame.z, y);
    auto x = cross(y, z);

    auto rot = rotation_frame(vec3f{1, 0, 0}, rotate.y) *
               frame3f{frame.x, frame.y, frame.z, vec3f{0, 0, 0}} *
               rotation_frame(vec3f{0, 1, 0}, rotate.x);
    auto pos = frame.o + transl.x * x + transl.y * y + transl.z * z;

    frame = {rot.x, rot.y, rot.z, pos};
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF IMAGE UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Gets pixels in an image region
template <typename T>
inline image<T> get_image_region(const image<T>& img, const image_region& region) {
    auto clipped = make_image<T>(region.size);
    for (auto j = 0; j < region.size.y; j++) {
        for (auto i = 0; i < region.size.x; i++) {
            at(clipped, {i, j}) = at(
                img, {i + region.offset.x, j + region.offset.y});
        }
    }
    return clipped;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF STRING/TIME UTILITIES FOR CLI APPLICATIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Prints basic types
inline bool print_value(string& str, const string& value) {
    str += value;
    return true;
}
inline bool print_value(string& str, const char* value) {
    str += value;
    return true;
}
inline bool print_value(string& str, int value) {
    str += std::to_string(value);
    return true;
}
inline bool print_value(string& str, float value) {
    str += std::to_string(value);
    return true;
}
inline bool print_value(string& str, double value) {
    str += std::to_string(value);
    return true;
}
template <typename T>
inline bool print_value(string& str, const T* value) {
    char buffer[512];
    sprintf(buffer, "%p", value);
    str += buffer;
    return true;
}

template <typename T, size_t N>
inline bool print_value(string& str, const array<T, N>& value) {
    for (auto i = 0; i < N; i++) {
        if (i) str += " ";
        str += std::to_string(value[i]);
    }
    return true;
}

// Print compound types.
template <typename T, int N>
inline bool print_value(string& str, const vec<T, N>& v) {
    return print_value(str, (const array<T, N>&)v);
}
template <typename T, int N, int M>
inline bool print_value(string& str, const mat<T, N, M>& v) {
    return print_value(str, (const array<T, N * M>&)v);
}
template <typename T, int N>
inline bool print_value(string& str, const frame<T, N>& v) {
    return print_value(str, (const array<T, N*(N + 1)>&)v);
}
template <typename T, int N>
inline bool print_value(string& str, const bbox<T, N>& v) {
    return print_value(str, (const array<T, N * 2>&)v);
}
template <typename T, int N>
inline bool print_value(string& str, const ray<T, N>& v) {
    return print_value(str, (const array<T, N * 2 + 2>&)v);
}
inline bool print_value(string& str, const image_region& v) {
    return print_value(str, (const array<int, 4>&)v);
}

// Prints a string.
inline bool print_next(string& str, const string& fmt) {
    return print_value(str, fmt);
}
template <typename Arg, typename... Args>
inline bool print_next(
    string& str, const string& fmt, const Arg& arg, const Args&... args) {
    auto pos = fmt.find("{}");
    if (pos == string::npos) return print_value(str, fmt);
    if (!print_value(str, fmt.substr(0, pos))) return false;
    if (!print_value(str, arg)) return false;
    return print_next(str, fmt.substr(pos + 2), args...);
}

// Formats a string `fmt` with values taken from `args`. Uses `{}` as
// placeholder.
template <typename... Args>
inline string format(const string& fmt, const Args&... args) {
    auto str = string();
    print_next(str, fmt, args...);
    return str;
}

// Prints a string.
template <typename... Args>
inline bool print(FILE* fs, const string& fmt, const Args&... args) {
    auto str = format(fmt, args...);
    return fprintf(fs, "%s", str.c_str()) >= 0;
}

// Converts to string.
template <typename T>
inline string to_string(const T& value) {
    auto str = string();
    print_value(str, value);
    return str;
}

// Trivial wrapper used for simplicity
struct parse_string_view {
    const char* str = nullptr;
};

// Prints basic types to string
inline bool parse_value(parse_string_view& str, string& value) {
    auto n = 0;
    char buffer[4096];
    if (sscanf(str.str, "%4095s%n", buffer, &n) != 1) return false;
    value = buffer;
    str.str += n;
    return true;
}
inline bool parse_value(parse_string_view& str, bool& value) {
    auto n = 0;
    auto v = 0;
    if (sscanf(str.str, "%d%n", &v, &n) != 1) return false;
    str.str += n;
    value = (bool)v;
    return true;
}
inline bool parse_value(parse_string_view& str, int& value) {
    auto n = 0;
    if (sscanf(str.str, "%d%n", &value, &n) != 1) return false;
    str.str += n;
    return true;
}
inline bool parse_value(parse_string_view& str, float& value) {
    auto n = 0;
    if (sscanf(str.str, "%f%n", &value, &n) != 1) return false;
    str.str += n;
    return true;
}
inline bool parse_value(parse_string_view& str, double& value) {
    auto n = 0;
    if (sscanf(str.str, "%lf%n", &value, &n) != 1) return false;
    str.str += n;
    return true;
}

// Print compound types
template <typename T, size_t N>
inline bool parse_value(parse_string_view& str, array<T, N>& value) {
    for (auto i = 0; i < N; i++) {
        if (!parse_value(str, value[i])) return false;
    }
    return true;
}
// Data acess
template <typename T, int N>
inline bool parse_value(parse_string_view& str, vec<T, N>& v) {
    return parse_value(str, (array<T, N>&)v);
}
template <typename T, int N, int M>
inline bool parse_value(parse_string_view& str, mat<T, N, M>& v) {
    return parse_value(str, (array<T, N * M>&)v);
}
template <typename T, int N>
inline bool parse_value(parse_string_view& str, frame<T, N>& v) {
    return parse_value(str, (array<T, N*(N + 1)>&)v);
}
template <typename T, int N>
inline bool parse_value(parse_string_view& str, bbox<T, N>& v) {
    return parse_value(str, (array<T, N * 2>&)v);
}
template <typename T, int N>
inline bool parse_value(parse_string_view& str, ray<T, N>& v) {
    return parse_value(str, (array<T, N * 2 + 2>&)v);
}
inline bool parse_value(parse_string_view& str, image_region& v) {
    return parse_value(str, (array<int, 4>&)v);
}

// Prints a string.
inline bool parse_next(parse_string_view& str) { return true; }
template <typename Arg, typename... Args>
inline bool parse_next(parse_string_view& str, Arg& arg, Args&... args) {
    if (!parse_value(str, arg)) return false;
    return parse_next(str, args...);
}

// Returns trus if this is white space
inline bool is_whitespace(parse_string_view str) {
    while (*str.str && isspace((unsigned char)*str.str)) str.str++;
    return *str.str == 0;
}

// Parse a list of space separated values.
template <typename... Args>
inline bool parse(const string& str, Args&... args) {
    auto view = parse_string_view{str.c_str()};
    if (!parse_next(view, args...)) return false;
    return is_whitespace(view);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF LOGGING UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Logging configutation
inline bool& _log_console() {
    static auto _log_console = true;
    return _log_console;
}
inline FILE*& _log_filestream() {
    static auto _log_filestream = (FILE*)nullptr;
    return _log_filestream;
}

// Logs a message
inline void log_message(const char* lbl, const char* msg) {
    if (_log_console()) {
        printf("%s\n", msg);
        fflush(stdout);
    }
    if (_log_filestream()) {
        fprintf(_log_filestream(), "%s %s\n", lbl, msg);
        fflush(_log_filestream());
    }
}

// Log info/error/fatal/trace message
template <typename... Args>
inline void log_info(const string& fmt, const Args&... args) {
    log_message("INFO ", format(fmt, args...).c_str());
}
template <typename... Args>
inline void log_error(const string& fmt, const Args&... args) {
    log_message("ERROR", format(fmt, args...).c_str());
}
template <typename... Args>
inline void log_fatal(const string& fmt, const Args&... args) {
    log_message("FATAL", format(fmt, args...).c_str());
    exit(1);
}

// Log traces for timing and program debugging
struct log_scope {
    string  message    = "";
    int64_t start_time = -1;
    bool    scoped     = false;
    ~log_scope();
};
template <typename... Args>
inline void log_trace(const string& fmt, const Args&... args) {
    log_message("TRACE", format(fmt, args...).c_str());
}
template <typename... Args>
inline log_scope log_trace_begin(const string& fmt, const Args&... args) {
    auto message = format(fmt, args...);
    log_trace(message + " [started]");
    return {message, get_time(), false};
}
template <typename... Args>
inline void log_trace_end(log_scope& scope) {
    if (scope.start_time >= 0) {
        log_trace(scope.message + " [ended: " +
                  format_duration(get_time() - scope.start_time) + "]");
    } else {
        log_trace(scope.message + " [ended]");
    }
}
template <typename... Args>
inline log_scope log_trace_scoped(const string& fmt, const Args&... args) {
    auto message = format(fmt, args...);
    log_trace(message + " [started]");
    return {message, get_time(), true};
}
inline log_scope::~log_scope() {
    if (scoped) log_trace_end(*this);
}

// Configure the logging
inline void set_log_console(bool enabled) { _log_console() = enabled; }
inline void set_log_file(const string& filename, bool append) {
    if (_log_filestream()) {
        fclose(_log_filestream());
        _log_filestream() = nullptr;
    }
    if (filename.empty()) return;
    _log_filestream() = fopen(filename.c_str(), append ? "at" : "wt");
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF STRING FORMAT UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Format duration string from nanoseconds
inline string format_duration(int64_t duration) {
    auto elapsed = duration / 1000000;  // milliseconds
    auto hours   = (int)(elapsed / 3600000);
    elapsed %= 3600000;
    auto mins = (int)(elapsed / 60000);
    elapsed %= 60000;
    auto secs  = (int)(elapsed / 1000);
    auto msecs = (int)(elapsed % 1000);
    char buffer[256];
    sprintf(buffer, "%02d:%02d:%02d.%03d", hours, mins, secs, msecs);
    return buffer;
}
// Format a large integer number in human readable form
inline string format_num(uint64_t num) {
    auto rem = num % 1000;
    auto div = num / 1000;
    if (div > 0) return format_num(div) + "," + std::to_string(rem);
    return std::to_string(rem);
}

}  // namespace ygl

#endif
