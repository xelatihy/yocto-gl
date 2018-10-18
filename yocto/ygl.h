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
// Yocto/GL tries to follow a "data-driven programming model" that makes data
// explicit. Data is stored in simple structs and access with free functions
// or directly. All data is public, so we make no attempt at encapsulation.
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
// in graphics. In particular, we support 2-4 dimensional vectors `vec<T, 2>`,
// `vec<T, 3>`, `vec<T, 4>`.
//
// We support 2-4 dimensional generic matrices `mat<T, N, 2>`, `mat<T, N, 3>`,
// `mat<T, N, 4>`, with matrix-matrix and matrix-vector products, transposes and
// inverses. Matrices are stored in column-major ordered and are accessed and
// constructed by column.
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
// 4. compute smooth normals and tangents with `compute_normals()`
//   `compute_tangents()`
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
// 9.  sample a could of point over a surface with `sample_triangles_points()`
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
//     4. generate random integers in an interval with `rand1i()`
//     5. generate random floats and double in the [0,1) range with
//        `rand1f()`, `rand2f()`, `rand3f()`, `next_rand1d()`
// 2. Perlin noise: `perlin_noise()` to generate Perlin noise with optional
//    wrapping, with fractal variations `perlin_ridge_noise()`,
//    `perlin_fbm_noise()`, `perlin_turbulence_noise()`
// 3. Monte Carlo support: warp functions from [0,1)^k domains to domains
//    commonly used in path tracing. In particular, use `sample_hemisphere()`,
//    `sample_sphere()`, `sample_hemisphere_cosine()`,
//    `sample_hemisphere_cospower()`. `sample_disk()`. `sample_cylinder()`.
//    `sample_triangle()`, `sample_discrete()`. For each warp, you can compute
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
// 2. build the BVH with `build_bvh()`
// 3. perform ray-element intersection with `intersect_bvh()`
// 4. perform point overlap queries with `overlap_bvh()`
// 5. refit the BVH with `refit_bvh()` after updating internal data
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
// 5. make various image examples with the `make_XXX_image4f()` functions
// 6. create procedural sun-sky images with `make_sunsky_image4f()`
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
// 3. use compute_bbox()` to compute element bounds
// 4. can merge scene together with `merge_into()`
// 5. make scene elements with `make_XXX()` functions
// 8. for ray-intersection and closest point queries, a BVH can be created with
//    `update_bvh()` and refit with `refit_bvh()`
// 9. compute interpolated values over scene elements with `eval_XXX()`
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
// 1. prepare the ray-tracing acceleration structure with `build_bvh()`
// 2. prepare lights for rendering with `make_trace_lights()`
// 3. create the image buffer and random number generators `make_trace_state()`
// 4. render blocks of samples with `trace_samples()`
// 5. you can also start an asynchronous renderer with `trace_asynch_start()`
//
//
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

#include <algorithm>  // for std::upper_bound
#include <array>
#include <atomic>
#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>  // for std::hash
#include <limits>
#include <map>
#include <memory>
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
using std::function;
using std::make_unique;
using std::map;
using std::runtime_error;
using std::string;
using std::thread;
using std::tie;
using std::tuple;
using std::unique_ptr;
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
inline T clamp(const T& x, const T& min_, const T& max_) {
    return min(max(x, min_), max_);
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
inline T mint() {
    return std::numeric_limits<T>::lowest();
}
template <typename T>
inline T maxt() {
    return std::numeric_limits<T>::max();
}
template <typename T>
inline T epst() {
    return std::numeric_limits<T>::epsilon();
}
const auto maxf = maxt<float>();
const auto epsf = epst<float>();

template <class T>
const T    pi  = (T)3.14159265358979323846;
const auto pif = 3.14159265f;

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
using vec2f = vec<float, 2>;
using vec3f = vec<float, 3>;
using vec4f = vec<float, 4>;
using vec2i = vec<int, 2>;
using vec3i = vec<int, 3>;
using vec4i = vec<int, 4>;
using vec4b = vec<byte, 4>;

// Zero vector constants.
const auto zero2f = vec2f{0, 0};
const auto zero3f = vec3f{0, 0, 0};
const auto zero4f = vec4f{0, 0, 0, 0};
const auto zero2i = vec2i{0, 0};
const auto zero3i = vec3i{0, 0, 0};
const auto zero4i = vec4i{0, 0, 0, 0};
const auto zero4b = vec4b{0, 0, 0, 0};

// Access component by index.
template <typename T>
inline const T& element_at(const vec<T, 3>& v, int i) {
    return *(&v.x + i);
}
template <typename T>
inline T& element_at(vec<T, 3>& v, int i) {
    return *(&v.x + i);
}

// Access xyz component of a vec4 typically used for color operation.
template <typename T>
inline vec<T, 3>& xyz(const vec<T, 4>& a) {
    return (vec<T, 3>&)a;
}
template <typename T>
inline vec<T, 3>& xyz(vec<T, 4>& a) {
    return (vec<T, 3>&)a;
}

// Vector comparison operations.
template <typename T>
inline bool operator==(const vec<T, 2>& a, const vec<T, 2>& b) {
    return a.x == b.x && a.y == b.y;
}
template <typename T>
inline bool operator!=(const vec<T, 2>& a, const vec<T, 2>& b) {
    return a.x != b.x || a.y != b.y;
}
template <typename T, typename T1>
inline bool operator==(const vec<T, 2>& a, T1 b) {
    return a.x == b && a.y == b;
}
template <typename T, typename T1>
inline bool operator!=(const vec<T, 2>& a, T1 b) {
    return a.x != b || a.y != b;
}
template <typename T>
inline bool operator==(const vec<T, 3>& a, const vec<T, 3>& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}
template <typename T>
inline bool operator!=(const vec<T, 3>& a, const vec<T, 3>& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z;
}
template <typename T, typename T1>
inline bool operator==(const vec<T, 3>& a, T1 b) {
    return a.x == b && a.y == b && a.z == b;
}
template <typename T, typename T1>
inline bool operator!=(const vec<T, 3>& a, T1 b) {
    return a.x != b || a.y != b || a.z != b;
}
template <typename T>
inline bool operator==(const vec<T, 4>& a, const vec<T, 4>& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}
template <typename T>
inline bool operator!=(const vec<T, 4>& a, const vec<T, 4>& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w;
}
template <typename T, typename T1>
inline bool operator==(const vec<T, 4>& a, T1 b) {
    return a.x == b && a.y == b && a.z == b && a.w == b;
}
template <typename T, typename T1>
inline bool operator!=(const vec<T, 4>& a, T1 b) {
    return a.x != b || a.y != b || a.z != b || a.w != b;
}

// Vector operations.
template <typename T>
inline vec<T, 2> operator-(const vec<T, 2>& a) {
    return {-a.x, -a.y};
}
template <typename T>
inline vec<T, 2> operator+(const vec<T, 2>& a, const vec<T, 2>& b) {
    return {a.x + b.x, a.y + b.y};
}
template <typename T, typename T1>
inline vec<T, 2> operator+(const vec<T, 2>& a, T1 b) {
    return {a.x + b, a.y + b};
}
template <typename T, typename T1>
inline vec<T, 2> operator+(T1 a, const vec<T, 2>& b) {
    return {a + b.x, a + b.y};
}
template <typename T>
inline vec<T, 2> operator-(const vec<T, 2>& a, const vec<T, 2>& b) {
    return {a.x - b.x, a.y - b.y};
}
template <typename T, typename T1>
inline vec<T, 2> operator-(const vec<T, 2>& a, T1 b) {
    return {a.x - b, a.y - b};
}
template <typename T, typename T1>
inline vec<T, 2> operator-(T1 a, const vec<T, 2>& b) {
    return {a - b.x, a - b.y};
}
template <typename T>
inline vec<T, 2> operator*(const vec<T, 2>& a, const vec<T, 2>& b) {
    return {a.x * b.x, a.y * b.y};
}
template <typename T, typename T1>
inline vec<T, 2> operator*(const vec<T, 2>& a, T1 b) {
    return {a.x * b, a.y * b};
}
template <typename T, typename T1>
inline vec<T, 2> operator*(T1 a, const vec<T, 2>& b) {
    return {a * b.x, a * b.y};
}
template <typename T>
inline vec<T, 2> operator/(const vec<T, 2>& a, const vec<T, 2>& b) {
    return {a.x / b.x, a.y / b.y};
}
template <typename T, typename T1>
inline vec<T, 2> operator/(const vec<T, 2>& a, T1 b) {
    return {a.x / b, a.y / b};
}
template <typename T, typename T1>
inline vec<T, 2> operator/(T1 a, const vec<T, 2>& b) {
    return {a / b.x, a / b.y};
}

// Vector operations.
template <typename T>
inline vec<T, 3> operator+(const vec<T, 3>& a) {
    return a;
}
template <typename T>
inline vec<T, 3> operator-(const vec<T, 3>& a) {
    return {-a.x, -a.y, -a.z};
}
template <typename T>
inline vec<T, 3> operator+(const vec<T, 3>& a, const vec<T, 3>& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
template <typename T, typename T1>
inline vec<T, 3> operator+(const vec<T, 3>& a, T1 b) {
    return {a.x + b, a.y + b, a.z + b};
}
template <typename T, typename T1>
inline vec<T, 3> operator+(T1 a, const vec<T, 3>& b) {
    return {a + b.x, a + b.y, a + b.z};
}
template <typename T>
inline vec<T, 3> operator-(const vec<T, 3>& a, const vec<T, 3>& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
template <typename T, typename T1>
inline vec<T, 3> operator-(const vec<T, 3>& a, T1 b) {
    return {a.x - b, a.y - b, a.z - b};
}
template <typename T, typename T1>
inline vec<T, 3> operator-(T1 a, const vec<T, 3>& b) {
    return {a - b.x, a - b.y, a - b.z};
}
template <typename T>
inline vec<T, 3> operator*(const vec<T, 3>& a, const vec<T, 3>& b) {
    return {a.x * b.x, a.y * b.y, a.z * b.z};
}
template <typename T, typename T1>
inline vec<T, 3> operator*(const vec<T, 3>& a, T1 b) {
    return {a.x * b, a.y * b, a.z * b};
}
template <typename T, typename T1>
inline vec<T, 3> operator*(T1 a, const vec<T, 3>& b) {
    return {a * b.x, a * b.y, a * b.z};
}
template <typename T>
inline vec<T, 3> operator/(const vec<T, 3>& a, const vec<T, 3>& b) {
    return {a.x / b.x, a.y / b.y, a.z / b.z};
}
template <typename T, typename T1>
inline vec<T, 3> operator/(const vec<T, 3>& a, T1 b) {
    return {a.x / b, a.y / b, a.z / b};
}
template <typename T, typename T1>
inline vec<T, 3> operator/(T1 a, const vec<T, 3>& b) {
    return {a / b.x, a / b.y, a / b.z};
}

// Vector operations.
template <typename T>
inline vec<T, 4> operator-(const vec<T, 4>& a) {
    return {-a.x, -a.y, -a.z, -a.w};
}
template <typename T>
inline vec<T, 4> operator+(const vec<T, 4>& a, const vec<T, 4>& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}
template <typename T, typename T1>
inline vec<T, 4> operator+(const vec<T, 4>& a, T1 b) {
    return {a.x + b, a.y + b, a.z + b, a.w + b};
}
template <typename T, typename T1>
inline vec<T, 4> operator+(T1 a, const vec<T, 4>& b) {
    return {a + b.x, a + b.y, a + b.z, a + b.w};
}
template <typename T>
inline vec<T, 4> operator-(const vec<T, 4>& a, const vec<T, 4>& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w};
}
template <typename T, typename T1>
inline vec<T, 4> operator-(const vec<T, 4>& a, T1 b) {
    return {a.x - b, a.y - b, a.z - b, a.w - b};
}
template <typename T, typename T1>
inline vec<T, 4> operator-(T1 a, const vec<T, 4>& b) {
    return {a - b.x, a - b.y, a - b.z, a - b.w};
}
template <typename T>
inline vec<T, 4> operator*(const vec<T, 4>& a, const vec<T, 4>& b) {
    return {a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w};
}
template <typename T, typename T1>
inline vec<T, 4> operator*(const vec<T, 4>& a, T1 b) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
}
template <typename T, typename T1>
inline vec<T, 4> operator*(T1 a, const vec<T, 4>& b) {
    return {a * b.x, a * b.y, a * b.z, a * b.w};
}
template <typename T>
inline vec<T, 4> operator/(const vec<T, 4>& a, const vec<T, 4>& b) {
    return {a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w};
}
template <typename T, typename T1>
inline vec<T, 4> operator/(const vec<T, 4>& a, T1 b) {
    return {a.x / b, a.y / b, a.z / b, a.w / b};
}
template <typename T, typename T1>
inline vec<T, 4> operator/(T1 a, const vec<T, 4>& b) {
    return {a / b.x, a / b.y, a / b.z, a / b.w};
}

// Vector assignments
template <typename T, int N>
inline vec<T, N>& operator+=(vec<T, N>& a, const vec<T, N>& b) {
    return a = a + b;
}
template <typename T, int N, typename T1>
inline vec<T, N>& operator+=(vec<T, N>& a, T1 b) {
    return a = a + b;
}
template <typename T, int N>
inline vec<T, N>& operator-=(vec<T, N>& a, const vec<T, N>& b) {
    return a = a - b;
}
template <typename T, int N, typename T1>
inline vec<T, N>& operator-=(vec<T, N>& a, T1 b) {
    return a = a - b;
}
template <typename T, int N>
inline vec<T, N>& operator*=(vec<T, N>& a, const vec<T, N>& b) {
    return a = a * b;
}
template <typename T, int N, typename T1>
inline vec<T, N>& operator*=(vec<T, N>& a, T1 b) {
    return a = a * b;
}
template <typename T, int N>
inline vec<T, N>& operator/=(vec<T, N>& a, const vec<T, N>& b) {
    return a = a / b;
}
template <typename T, int N, typename T1>
inline vec<T, N>& operator/=(vec<T, N>& a, T1 b) {
    return a = a / b;
}

// Vector products and lengths.
template <typename T>
inline T dot(const vec<T, 2>& a, const vec<T, 2>& b) {
    return a.x * b.x + a.y * b.y;
}
template <typename T>
inline T cross(const vec<T, 2>& a, const vec<T, 2>& b) {
    return a.x * b.y - a.y * b.x;
}
template <typename T>
inline T dot(const vec<T, 3>& a, const vec<T, 3>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
template <typename T>
inline vec<T, 3> cross(const vec<T, 3>& a, const vec<T, 3>& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
template <typename T>
inline T dot(const vec<T, 4>& a, const vec<T, 4>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

template <typename T, int N>
inline T length(const vec<T, N>& a) {
    return sqrt(dot(a, a));
}
template <typename T, int N>
inline vec<T, N> normalize(const vec<T, N>& a) {
    auto l = length(a);
    return (l) ? a / l : a;
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
inline vec<T, 3> reflect(const vec<T, 3>& w, const vec<T, 3>& n) {
    return -w + 2 * dot(n, w) * n;
}
template <typename T, typename T1>
inline vec<T, 3> refract(const vec<T, 3>& w, const vec<T, 3>& n, T1 eta) {
    // auto k = 1.0 - eta * eta * (1.0 - dot(n, w) * dot(n, w));
    auto k = 1 - eta * eta * max(0.0f, 1 - dot(n, w) * dot(n, w));
    if (k < 0) return {0, 0, 0};  // tir
    return -w * eta + (eta * dot(n, w) - sqrt(k)) * n;
}

// Max element and clamp.
template <typename T, typename T1, typename T2>
inline vec<T, 2> clamp(const vec<T, 2>& x, T1 min, T2 max) {
    return {clamp(x.x, min, max), clamp(x.y, min, max)};
}
template <typename T, typename T1, typename T2>
inline vec<T, 3> clamp(const vec<T, 3>& x, T1 min, T2 max) {
    return {clamp(x.x, min, max), clamp(x.y, min, max), clamp(x.z, min, max)};
}
template <typename T, typename T1, typename T2>
inline vec<T, 4> clamp(const vec<T, 4>& x, T1 min, T2 max) {
    return {clamp(x.x, min, max), clamp(x.y, min, max), clamp(x.z, min, max),
        clamp(x.w, min, max)};
}
template <typename T>
inline T max(const vec<T, 2>& a) {
    return max(a.x, a.y);
}
template <typename T>
inline T max(const vec<T, 3>& a) {
    return max(max(a.x, a.y), a.z);
}
template <typename T>
inline T max(const vec<T, 4>& a) {
    return max(max(max(a.x, a.y), a.z), a.w);
}
template <typename T>
inline T min(const vec<T, 2>& a) {
    return min(a.x, a.y);
}
template <typename T>
inline T min(const vec<T, 3>& a) {
    return min(min(a.x, a.y), a.z);
}
template <typename T>
inline T min(const vec<T, 4>& a) {
    return min(min(min(a.x, a.y), a.z), a.w);
}
template <typename T>
inline T mean(const vec<T, 2>& a) {
    return (a.x + a.y) / 2;
}
template <typename T>
inline T mean(const vec<T, 3>& a) {
    return (a.x + a.y + a.z) / 3;
}
template <typename T>
inline T mean(const vec<T, 4>& a) {
    return (a.x + a.y + a.z + a.w) / 4;
}

// Quaternion operatons represented as xi + yj + zk + w
const auto identity_quat4f = vec4f{0, 0, 0, 1};
template <typename T>
inline vec<T, 4> quat_mul(const vec<T, 4>& a, float b) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
}
template <typename T>
inline vec<T, 4> quat_mul(const vec<T, 4>& a, const vec<T, 4>& b) {
    return {a.x * b.w + a.w * b.x + a.y * b.w - a.z * b.y,
        a.y * b.w + a.w * b.y + a.z * b.x - a.x * b.z,
        a.z * b.w + a.w * b.z + a.x * b.y - a.y * b.x,
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z};
}
template <typename T>
inline vec<T, 4> quat_conjugate(const vec<T, 4>& a) {
    return {-a.x, -a.y, -a.z, a.w};
}
template <typename T>
inline vec<T, 4> quat_inverse(const vec<T, 4>& a) {
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
struct mat<T, N, 2> {
    vec<T, N> x = {0, 0};
    vec<T, N> y = {0, 0};
};
template <typename T, int N>
struct mat<T, N, 3> {
    vec<T, 3> x = {0, 0, 0};
    vec<T, 3> y = {0, 0, 0};
    vec<T, 3> z = {0, 0, 0};
};
template <typename T, int N>
struct mat<T, N, 4> {
    vec<T, 4> x = {0, 0, 0, 0};
    vec<T, 4> y = {0, 0, 0, 0};
    vec<T, 4> z = {0, 0, 0, 0};
    vec<T, 4> w = {0, 0, 0, 0};
};

// Type aliases.
using mat2f = mat<float, 2, 2>;
using mat3f = mat<float, 3, 3>;
using mat4f = mat<float, 4, 4>;

// Identity matrices constants.
const auto identity_mat2f = mat2f{{1, 0}, {0, 1}};
const auto identity_mat3f = mat3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
const auto identity_mat4f = mat4f{
    {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};

// Matrix comparisons.
template <typename T, int N>
inline bool operator==(const mat<T, N, 2>& a, const mat<T, N, 2>& b) {
    return a.x == b.x && a.y == b.y;
}
template <typename T, int N>
inline bool operator!=(const mat<T, N, 2>& a, const mat<T, N, 2>& b) {
    return !(a == b);
}
template <typename T, int N>
inline bool operator==(const mat<T, N, 3>& a, const mat<T, N, 3>& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}
template <typename T, int N>
inline bool operator!=(const mat<T, N, 3>& a, const mat<T, N, 3>& b) {
    return !(a == b);
}
template <typename T, int N>
inline bool operator==(const mat<T, N, 4>& a, const mat<T, N, 4>& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}
template <typename T, int N>
inline bool operator!=(const mat<T, N, 4>& a, const mat<T, N, 4>& b) {
    return !(a == b);
}

// Matrix operations.
template <typename T, int N>
inline mat<T, N, 2> operator+(const mat<T, N, 2>& a, const mat<T, N, 2>& b) {
    return {a.x + b.x, a.y + b.y};
}
template <typename T, int N, typename T1>
inline mat<T, N, 2> operator*(const mat<T, N, 2>& a, T1 b) {
    return {a.x * b, a.y * b};
}
template <typename T, int N>
inline vec<T, 2> operator*(const mat<T, N, 2>& a, const vec2f& b) {
    return a.x * b.x + a.y * b.y;
}
template <typename T, int N>
inline vec<T, 2> operator*(const vec2f& a, const mat<T, N, 2>& b) {
    return {dot(a, b.x), dot(a, b.y)};
}
template <typename T, int N>
inline mat<T, N, 2> operator*(const mat<T, N, 2>& a, const mat<T, N, 2>& b) {
    return {a * b.x, a * b.y};
}

// Matrix operations.
template <typename T, int N>
inline mat<T, N, 3> operator+(const mat<T, N, 3>& a, const mat<T, N, 3>& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
template <typename T, int N, typename T1>
inline mat<T, N, 3> operator*(const mat<T, N, 3>& a, T1 b) {
    return {a.x * b, a.y * b, a.z * b};
}
template <typename T, int N>
inline vec<T, 3> operator*(const mat<T, N, 3>& a, const vec<T, 3>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
template <typename T, int N>
inline vec3f operator*(const vec3f& a, const mat<T, N, 3>& b) {
    return {dot(a, b.x), dot(a, b.y), dot(a, b.z)};
}
template <typename T, int N>
inline mat<T, N, 3> operator*(const mat<T, N, 3>& a, const mat<T, N, 3>& b) {
    return {a * b.x, a * b.y, a * b.z};
}

// Matrix operations.
template <typename T, int N>
inline mat<T, N, 4> operator+(const mat<T, N, 4>& a, const mat<T, N, 4>& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}
template <typename T, int N, typename T1>
inline mat<T, N, 4> operator*(const mat<T, N, 4>& a, T1 b) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
}
template <typename T, int N>
inline vec<T, 4> operator*(const mat<T, N, 4>& a, const vec<T, 4>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}
template <typename T, int N>
inline vec<T, 4> operator*(const vec<T, 4>& a, const mat<T, N, 4>& b) {
    return {dot(a, b.x), dot(a, b.y), dot(a, b.z), dot(a, b.w)};
}
template <typename T, int N>
inline mat<T, N, 4> operator*(const mat<T, N, 4>& a, const mat<T, N, 4>& b) {
    return {a * b.x, a * b.y, a * b.z, a * b.w};
}

// Matrix assignments.
template <typename T, int N, int M>
inline mat<T, N, M>& operator+=(mat<T, N, M>& a, const mat<T, N, M>& b) {
    return a = a + b;
}
template <typename T, int N, int M>
inline mat<T, N, M>& operator*=(mat<T, N, M>& a, const mat<T, N, M>& b) {
    return a = a * b;
}
template <typename T, int N, int M, typename T1>
inline mat<T, N, M>& operator*=(mat<T, N, M>& a, T1 b) {
    return a = a * b;
}

// Matrix diagonals and transposes.
template <typename T>
inline vec<T, 2> diagonal(const mat<T, 2, 2>& a) {
    return {a.x.x, a.y.y};
}
template <typename T>
inline vec<T, 3> diagonal(const mat<T, 3, 3>& a) {
    return {a.x.x, a.y.y, a.z.z};
}
template <typename T>
inline vec<T, 4> diagonal(const mat<T, 4, 4>& a) {
    return {a.x.x, a.y.y, a.z.z, a.w.w};
}
template <typename T>
inline mat<T, 2, 2> transpose(const mat<T, 2, 2>& a);
template <typename T>
inline mat<T, 3, 3> transpose(const mat<T, 3, 3>& a);
template <typename T>
inline mat<T, 4, 4> transpose(const mat<T, 4, 4>& a);

// Matrix adjugates, determinant and inverses.
template <typename T>
inline mat<T, 2, 2> adjugate(const mat<T, 2, 2>& a);
template <typename T>
inline mat<T, 3, 3> adjugate(const mat<T, 3, 3>& a);
template <typename T>
inline mat<T, 4, 4> adjugate(const mat<T, 4, 4>& a);
template <typename T>
inline T determinant(const mat<T, 2, 2>& a);
template <typename T>
inline T determinant(const mat<T, 3, 3>& a);
template <typename T>
inline T determinant(const mat<T, 4, 4>& a);
template <typename T, int N>
inline mat<T, N, N> inverse(const mat<T, N, N>& a) {
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
inline frame<T, 3> make_frame_fromz(const vec<T, 3>& o, const vec<T, 3>& v) {
    auto z = normalize(v);
    auto x = normalize(orthogonal(z));
    auto y = normalize(cross(z, x));
    return {x, y, z, o};
}
template <typename T>
inline frame<T, 3> make_frame_fromzx(
    const vec<T, 3>& o, const vec<T, 3>& z_, const vec<T, 3>& x_) {
    auto z = normalize(z_);
    auto x = orthonormalize(x_, z);
    auto y = normalize(cross(z, x));
    return {x, y, z, o};
}

// Frame to matrix conversion.
template <typename T>
inline mat<T, 3, 3> frame_to_mat(const frame<T, 2>& a) {
    return {
        {a.x.x, a.x.y, 0},
        {a.y.x, a.y.y, 0},
        {a.o.x, a.o.y, 1},
    };
}
template <typename T>
inline frame<T, 2> mat_to_frame(const mat<T, 3, 3>& a) {
    return {
        {a.x.x, a.x.y},
        {a.y.x, a.y.y},
        {a.z.x, a.z.y},
    };
}
template <typename T>
inline mat<T, 4, 4> frame_to_mat(const frame<T, 3>& a) {
    return {
        {a.x.x, a.x.y, a.x.z, 0},
        {a.y.x, a.y.y, a.y.z, 0},
        {a.z.x, a.z.y, a.z.z, 0},
        {a.o.x, a.o.y, a.o.z, 1},
    };
}
template <typename T>
inline frame<T, 3> mat_to_frame(const mat<T, 4, 4>& a) {
    return {
        {a.x.x, a.x.y, a.x.z},
        {a.y.x, a.y.y, a.y.z},
        {a.z.x, a.z.y, a.z.z},
        {a.w.x, a.w.y, a.w.z},
    };
}

// Frame comparisons.
template <typename T>
inline bool operator==(const frame<T, 2>& a, const frame<T, 2>& b) {
    return a.x == b.x && a.y == b.y && a.o == b.o;
}
template <typename T>
inline bool operator!=(const frame<T, 2>& a, const frame<T, 2>& b) {
    return !(a == b);
}
template <typename T>
inline bool operator==(const frame<T, 3>& a, const frame<T, 3>& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.o == b.o;
}
template <typename T>
inline bool operator!=(const frame<T, 3>& a, const frame<T, 3>& b) {
    return !(a == b);
}

// Frame composition, equivalent to affine matrix product.
template <typename T>
inline frame<T, 2> operator*(const frame<T, 2>& a, const frame<T, 2>& b) {
    auto rot = mat<T, 2, 2>{a.x, a.y} * mat<T, 2, 2>{b.x, b.y};
    auto pos = mat<T, 2, 2>{a.x, a.y} * b.o + a.o;
    return {rot.x, rot.y, pos};
}
template <typename T>
inline frame<T, 3> operator*(const frame<T, 3>& a, const frame<T, 3>& b) {
    auto rot = mat<T, 3, 3>{a.x, a.y, a.z} * mat<T, 3, 3>{b.x, b.y, b.z};
    auto pos = mat<T, 3, 3>{a.x, a.y, a.z} * b.o + a.o;
    return {rot.x, rot.y, rot.z, pos};
}
// Frame inverse, equivalent to rigid affine inverse.
template <typename T>
inline frame<T, 2> inverse(const frame<T, 2>& a, bool is_rigid = true) {
    auto minv = (is_rigid) ? transpose(mat<T, 2, 2>{a.x, a.y}) :
                             inverse(mat<T, 2, 2>{a.x, a.y});
    return {minv.x, minv.y, -(minv * a.o)};
}
template <typename T>
inline frame<T, 3> inverse(const frame<T, 3>& a, bool is_rigid = true) {
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
inline bool operator==(const bbox<T, 1>& a, const bbox<T, 1>& b) {
    return a.min == b.min && a.max == b.max;
}
template <typename T>
inline bool operator!=(const bbox<T, 1>& a, const bbox<T, 1>& b) {
    return a.min != b.min || a.max != b.max;
}
template <typename T>
inline bool operator==(const bbox<T, 2>& a, const bbox<T, 2>& b) {
    return a.min == b.min && a.max == b.max;
}
template <typename T>
inline bool operator!=(const bbox<T, 2>& a, const bbox<T, 2>& b) {
    return a.min != b.min || a.max != b.max;
}
template <typename T>
inline bool operator==(const bbox<T, 3>& a, const bbox<T, 3>& b) {
    return a.min == b.min && a.max == b.max;
}
template <typename T>
inline bool operator!=(const bbox<T, 3>& a, const bbox<T, 3>& b) {
    return a.min != b.min || a.max != b.max;
}
template <typename T>
inline bool operator==(const bbox<T, 4>& a, const bbox<T, 4>& b) {
    return a.min == b.min && a.max == b.max;
}
template <typename T>
inline bool operator!=(const bbox<T, 4>& a, const bbox<T, 4>& b) {
    return a.min != b.min || a.max != b.max;
}

// Bounding box expansions with points and other boxes.
template <typename T>
inline bbox<T, 1>& operator+=(bbox<T, 1>& a, T b) {
    a.min = min(a.min, b);
    a.max = max(a.max, b);
    return a;
}
template <typename T>
inline bbox<T, 1>& operator+=(bbox<T, 1>& a, const bbox<T, 1>& b) {
    a.min = min(a.min, b.min);
    a.max = max(a.max, b.max);
    return a;
}
// Bounding box expansions with points and other boxes.
template <typename T>
inline bbox<T, 2>& operator+=(bbox<T, 2>& a, const vec<T, 2>& b) {
    a.min = {min(a.min.x, b.x), min(a.min.y, b.y)};
    a.max = {max(a.max.x, b.x), max(a.max.y, b.y)};
    return a;
}
template <typename T>
inline bbox<T, 2>& operator+=(bbox<T, 2>& a, const bbox<T, 2>& b) {
    a.min = {min(a.min.x, b.min.x), min(a.min.y, b.min.y)};
    a.max = {max(a.max.x, b.max.x), max(a.max.y, b.max.y)};
    return a;
}
// Bounding box expansions with points and other boxes.
template <typename T>
inline bbox<T, 3>& operator+=(bbox<T, 3>& a, const vec<T, 3>& b) {
    a.min = {min(a.min.x, b.x), min(a.min.y, b.y), min(a.min.z, b.z)};
    a.max = {max(a.max.x, b.x), max(a.max.y, b.y), max(a.max.z, b.z)};
    return a;
}
template <typename T>
inline bbox<T, 3>& operator+=(bbox<T, 3>& a, const bbox<T, 3>& b) {
    a.min = {min(a.min.x, b.min.x), min(a.min.y, b.min.y), min(a.min.z, b.min.z)};
    a.max = {max(a.max.x, b.max.x), max(a.max.y, b.max.y), max(a.max.z, b.max.z)};
    return a;
}
// Bounding box expansions with points and other boxes.
template <typename T>
inline bbox<T, 4>& operator+=(bbox<T, 4>& a, const vec<T, 4>& b) {
    a.min = {min(a.min.x, b.x), min(a.min.y, b.y), min(a.min.z, b.z),
        min(a.min.w, b.w)};
    a.max = {max(a.max.x, b.x), max(a.max.y, b.y), max(a.max.z, b.z),
        max(a.max.w, b.w)};
    return a;
}
template <typename T>
inline bbox<T, 4>& operator+=(bbox<T, 4>& a, const bbox<T, 4>& b) {
    a.min = {min(a.min.x, b.min.x), min(a.min.y, b.min.y),
        min(a.min.z, b.min.z), min(a.min.w, b.min.w)};
    a.max = {max(a.max.x, b.max.x), max(a.max.y, b.max.y),
        max(a.max.z, b.max.z), max(a.max.w, b.max.w)};
    return a;
}

// Primitive bounds.
template <typename T>
inline bbox<T, 3> point_bbox(const vec<T, 3>& p, T r = 0) {
    auto bounds = bbox<T, 3>{};
    bounds += p - vec3f{r, r, r};
    bounds += p + vec3f{r, r, r};
    return bounds;
}
template <typename T>
inline bbox<T, 3> line_bbox(
    const vec<T, 3>& p0, const vec<T, 3>& p1, T r0 = 0, T r1 = 0) {
    auto bounds = bbox<T, 3>{};
    bounds += p0 - vec3f{r0, r0, r0};
    bounds += p0 + vec3f{r0, r0, r0};
    bounds += p1 - vec3f{r1, r1, r1};
    bounds += p1 + vec3f{r1, r1, r1};
    return bounds;
}
template <typename T>
inline bbox<T, 3> triangle_bbox(
    const vec<T, 3>& p0, const vec<T, 3>& p1, const vec<T, 3>& p2) {
    auto bounds = bbox<T, 3>{};
    bounds += p0;
    bounds += p1;
    bounds += p2;
    return bounds;
}
template <typename T>
inline bbox<T, 3> quad_bbox(const vec<T, 3>& p0, const vec<T, 3>& p1,
    const vec<T, 3>& p2, const vec<T, 3>& p3) {
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
inline ray<T, N> make_ray(const vec<T, N>& o, const vec<T, N>& d, T eps = 1e-4f) {
    return {o, d, eps, maxt<T>()};
}
template <typename T, int N>
inline ray<T, N> make_segment(
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
inline vec<T, 2> transform_point(const mat<T, 3, 3>& a, const vec<T, 2>& b) {
    auto tvb = a * vec3f{b.x, b.y, 1};
    return vec<T, 2>{tvb.x, tvb.y} / tvb.z;
}
template <typename T>
inline vec<T, 3> transform_point(const mat<T, 4, 4>& a, const vec<T, 3>& b) {
    auto tvb = a * vec<T, 4>{b.x, b.y, b.z, 1};
    return vec<T, 3>{tvb.x, tvb.y, tvb.z} / tvb.w;
}
template <typename T>
inline vec<T, 2> transform_vector(const mat<T, 3, 3>& a, const vec<T, 2>& b) {
    auto tvb = a * vec<T, 3>{b.x, b.y, 0};
    return vec<T, 2>{tvb.x, tvb.y} / tvb.z;
}
template <typename T>
inline vec<T, 3> transform_vector(const mat<T, 3, 3>& a, const vec<T, 3>& b) {
    return a * b;
}
template <typename T>
inline vec<T, 3> transform_vector(const mat<T, 4, 4>& a, const vec<T, 3>& b) {
    auto tvb = a * vec4f{b.x, b.y, b.z, 0};
    return vec<T, 3>{tvb.x, tvb.y, tvb.z};
}
template <typename T, int N>
inline vec<T, N> transform_direction(const mat<T, 4, 4>& a, const vec<T, N>& b) {
    return normalize(transform_vector(a, b));
}

// Transforms points, vectors and directions by frames.
template <typename T>
inline vec<T, 2> transform_point(const frame<T, 2>& a, const vec<T, 2>& b) {
    return a.x * b.x + a.y * b.y + a.o;
}
template <typename T>
inline vec<T, 3> transform_point(const frame<T, 3>& a, const vec<T, 3>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.o;
}
template <typename T>
inline vec<T, 2> transform_vector(const frame<T, 2>& a, const vec<T, 2>& b) {
    return a.x * b.x + a.y * b.y;
}
template <typename T>
inline vec<T, 3> transform_vector(const frame<T, 3>& a, const vec<T, 3>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
template <typename T, int N>
inline vec<T, N> transform_direction(const frame<T, N>& a, const vec<T, N>& b) {
    return normalize(transform_vector(a, b));
}

// Transforms rays and bounding boxes by matrices.
template <typename T, int N>
inline ray<T, N> transform_ray(const frame<T, N>& a, const ray<T, N>& b) {
    return {transform_point(a, b.o), transform_vector(a, b.d), b.tmin, b.tmax};
}
template <typename T, int N>
inline ray<T, N> transform_ray(const mat<T, N + 1, N + 1>& a, const ray<T, N>& b) {
    return {transform_point(a, b.o), transform_vector(a, b.d), b.tmin, b.tmax};
}
template <typename T>
inline bbox<T, 3> transform_bbox(const frame<T, 3>& a, const bbox<T, 3>& b) {
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
inline bbox<T, 3> transform_bbox(const mat<T, 4, 4>& a, const bbox<T, 3>& b) {
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
inline vec<T, 2> transform_point_inverse(
    const frame<T, 2>& a, const vec<T, 2>& b) {
    return {dot(b - a.o, a.x), dot(b - a.o, a.y)};
}
template <typename T>
inline vec3f transform_point_inverse(const frame<T, 3>& a, const vec<T, 3>& b) {
    return {dot(b - a.o, a.x), dot(b - a.o, a.y), dot(b - a.o, a.z)};
}
template <typename T>
inline vec<T, 2> transform_vector_inverse(
    const frame<T, 2>& a, const vec<T, 2>& b) {
    return {dot(b, a.x), dot(b, a.y)};
}
template <typename T>
inline vec3f transform_vector_inverse(const frame<T, 3>& a, const vec<T, 3>& b) {
    return {dot(b, a.x), dot(b, a.y), dot(b, a.z)};
}
template <typename T, int N>
inline vec3f transform_direction_inverse(
    const frame<T, N>& a, const vec<T, N>& b) {
    return normalize(transform_vector_inverse(a, b));
}
template <typename T, int N>
inline ray<T, N> transform_ray_inverse(const frame<T, N>& a, const ray<T, N>& b) {
    return {transform_point_inverse(a, b.o),
        transform_direction_inverse(a, b.d), b.tmin, b.tmax};
}
template <typename T>
inline bbox<T, 3> transform_bbox_inverse(
    const frame<T, 3>& a, const bbox<T, 3>& b) {
    return transform_bbox(inverse(a), b);
}

// Translation, scaling and rotations transforms.
template <typename T>
inline frame<T, 3> translation_frame(const vec<T, 3>& a) {
    return {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, a};
}
template <typename T>
inline frame<T, 3> scaling_frame(const vec<T, 3>& a) {
    return {{a.x, 0, 0}, {0, a.y, 0}, {0, 0, a.z}, {0, 0, 0}};
}
template <typename T>
inline frame<T, 3> rotation_frame(const vec<T, 3>& axis, T angle) {
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
inline frame<T, 3> rotation_frame(const vec<T, 4>& quat) {
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
inline frame<T, 3> rotation_frame(const mat<T, 3, 3>& rot) {
    return {rot.x, rot.y, rot.z, {0, 0, 0}};
}

// Lookat frame. Z-axis can be inverted with inv_xz.
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

// OpenGL frustum, ortho and perspecgive matrices.
template <typename T>
inline mat<T, 4, 4> frustum_mat(T l, T r, T b, T t, T n, T f) {
    return {{2 * n / (r - l), 0, 0, 0}, {0, 2 * n / (t - b), 0, 0},
        {(r + l) / (r - l), (t + b) / (t - b), -(f + n) / (f - n), -1},
        {0, 0, -2 * f * n / (f - n), 0}};
}
template <typename T>
inline mat<T, 4, 4> ortho_mat(T l, T r, T b, T t, T n, T f) {
    return {{2 / (r - l), 0, 0, 0}, {0, 2 / (t - b), 0, 0},
        {0, 0, -2 / (f - n), 0},
        {-(r + l) / (r - l), -(t + b) / (t - b), -(f + n) / (f - n), 1}};
}
template <typename T>
inline mat<T, 4, 4> ortho2d_mat(T left, T right, T bottom, T top) {
    return ortho_mat(left, right, bottom, top, -1, 1);
}
template <typename T>
inline mat<T, 4, 4> ortho_mat(T xmag, T ymag, T near, T far) {
    return {{1 / xmag, 0, 0, 0}, {0, 1 / ymag, 0, 0},
        {0, 0, 2 / (near - far), 0}, {0, 0, (far + near) / (near - far), 1}};
}
template <typename T>
inline mat<T, 4, 4> perspective_mat(T fovy, T aspect, T near, T far) {
    auto tg = tan(fovy / 2);
    return {{1 / (aspect * tg), 0, 0, 0}, {0, 1 / tg, 0, 0},
        {0, 0, (far + near) / (near - far), -1},
        {0, 0, 2 * far * near / (near - far), 0}};
}
template <typename T>
inline mat<T, 4, 4> perspective_mat(T fovy, T aspect, T near) {
    auto tg = tan(fovy / 2);
    return {{1 / (aspect * tg), 0, 0, 0}, {0, 1 / tg, 0, 0}, {0, 0, -1, -1},
        {0, 0, 2 * near, 0}};
}

// Rotation conversions.
template <typename T>
inline tuple<vec<T, 3>, T> rotation_axisangle(const vec<T, 4>& quat) {
    return {normalize(vec3f{quat.x, quat.y, quat.z}), 2 * acos(quat.w)};
}
template <typename T>
inline vec<T, 4> rotation_quat(const vec<T, 3>& axis, T angle) {
    auto len = length(axis);
    if (!len) return {0, 0, 0, 1};
    return vec4f{sin(angle / 2) * axis.x / len, sin(angle / 2) * axis.y / len,
        sin(angle / 2) * axis.z / len, cos(angle / 2)};
}
template <typename T>
inline vec<T, 4> rotation_quat(const vec<T, 4>& axisangle) {
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
inline void center_image4f(vec2f& center, float& scale, const vec2i& imsize,
    const vec2i& winsize, bool zoom_to_fit) {
    if (zoom_to_fit) {
        scale  = min(winsize.x / (float)imsize.x, winsize.y / (float)imsize.y);
        center = {(float)winsize.x / 2, (float)winsize.y / 2};
    } else {
        if (winsize.x >= imsize.x * scale) center.x = winsize.x / 2;
        if (winsize.y >= imsize.y * scale) center.y = winsize.y / 2;
    }
}

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
inline int   rand1i(rng_state& rng, int n) { return advance_rng(rng) % n; }
inline float rand1f(rng_state& rng) {
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
inline vec2f rand2f(rng_state& rng) { return {rand1f(rng), rand1f(rng)}; }
inline vec3f rand3f(rng_state& rng) {
    return {rand1f(rng), rand1f(rng), rand1f(rng)};
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// MONETACARLO SAMPLING FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Sample an hemispherical direction with uniform distribution.
inline vec3f sample_hemisphere(const vec2f& ruv) {
    auto z   = ruv.y;
    auto r   = sqrt(1 - z * z);
    auto phi = 2 * pif * ruv.x;
    return {r * cos(phi), r * sin(phi), z};
}
inline float sample_hemisphere_pdf(const vec3f& w) {
    return (w.z <= 0) ? 0 : 1 / (2 * pif);
}

// Sample a spherical direction with uniform distribution.
inline vec3f sample_sphere(const vec2f& ruv) {
    auto z   = 2 * ruv.y - 1;
    auto r   = sqrt(1 - z * z);
    auto phi = 2 * pif * ruv.x;
    return {r * cos(phi), r * sin(phi), z};
}
inline float sample_sphere_pdf(const vec3f& w) { return 1 / (4 * pif); }

// Sample spherical coordinates uniformly.
inline vec2f sample_spherical(const vec2f& ruv) {
    // BUG: FIXME this is not uniform at all!!!!
    return {ruv.x, ruv.y};
}
inline float sample_spherical_pdf(const vec2f& w) { return 1 / (4 * pif); }

// Sample an hemispherical direction with cosine distribution.
inline vec3f sample_hemisphere_cosine(const vec2f& ruv) {
    auto z   = sqrt(ruv.y);
    auto r   = sqrt(1 - z * z);
    auto phi = 2 * pif * ruv.x;
    return {r * cos(phi), r * sin(phi), z};
}
inline float sample_hemisphere_cosine_pdf(const vec3f& w) {
    return (w.z <= 0) ? 0 : w.z / pif;
}

// Sample an hemispherical direction with cosine power distribution.
inline vec3f sample_hemisphere_cospower(float n, const vec2f& ruv) {
    auto z   = pow(ruv.y, 1 / (n + 1));
    auto r   = sqrt(1 - z * z);
    auto phi = 2 * pif * ruv.x;
    return {r * cos(phi), r * sin(phi), z};
}
inline float sample_hemisphere_cospower_pdf(float n, const vec3f& w) {
    return (w.z <= 0) ? 0 : pow(w.z, n) * (n + 1) / (2 * pif);
}

// Sample a point uniformly on a disk.
inline vec3f sample_disk(const vec2f& ruv) {
    auto r   = sqrt(ruv.y);
    auto phi = 2 * pif * ruv.x;
    return {cos(phi) * r, sin(phi) * r, 0};
}
inline float sample_disk_pdf() { return 1 / pif; }

// Sample a point uniformly on a cylinder, without caps.
inline vec3f sample_cylinder(const vec2f& ruv) {
    auto phi = 2 * pif * ruv.x;
    return {sin(phi), cos(phi), ruv.y * 2 - 1};
}
inline float sample_cylinder_pdf() { return 1 / pif; }

// Sample a point uniformly on a triangle.
inline vec2f sample_triangle(const vec2f& ruv) {
    return {1 - sqrt(ruv.x), ruv.y * sqrt(ruv.x)};
}
inline vec3f sample_triangle(
    const vec3f& p0, const vec3f& p1, const vec3f& p2, const vec2f& ruv) {
    auto uv = sample_triangle(ruv);
    return p0 * (1 - uv.x - uv.y) + p1 * uv.x + p2 * uv.y;
}
// Pdf for uniform triangle sampling, i.e. triangle area.
inline float sample_triangle_pdf(
    const vec3f& p0, const vec3f& p1, const vec3f& p2) {
    return 2 / length(cross(p1 - p0, p2 - p0));
}

// Sample an index with uniform distribution.
inline int sample_index(int size, float r) {
    return clamp((int)(r * size), 0, size - 1);
}
inline float sample_index_pdf(int size) { return 1.0f / size; }

// Sample a discrete distribution represented by its cdf.
inline int sample_discrete(const vector<float>& cdf, float r) {
    r        = clamp(r * cdf.back(), 0.0f, cdf.back() - 0.00001f);
    auto idx = (int)(std::upper_bound(cdf.data(), cdf.data() + cdf.size(), r) -
                     cdf.data());
    return clamp(idx, 0, (int)cdf.size() - 1);
}
// Pdf for uniform discrete distribution sampling.
inline float sample_discrete_pdf(const vector<float>& cdf, int idx) {
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
vector<vec3f> compute_tangents(
    const vector<vec2i>& lines, const vector<vec3f>& pos);
vector<vec3f> compute_normals(
    const vector<vec3i>& triangles, const vector<vec3f>& pos);
vector<vec3f> compute_normals(
    const vector<vec4i>& quads, const vector<vec3f>& pos);

// Compute per-vertex tangent space for triangle meshes.
// Tangent space is defined by a four component vector.
// The first three components are the tangent with respect to the u texcoord.
// The fourth component is the sign of the tangent wrt the v texcoord.
// Tangent frame is useful in normal mapping.
vector<vec4f> compute_tangent_space(const vector<vec3i>& triangles,
    const vector<vec3f>& pos, const vector<vec3f>& norm,
    const vector<vec2f>& texcoord);

// Apply skinning to vertex position and normals.
tuple<vector<vec3f>, vector<vec3f>> compute_skinning(const vector<vec3f>& pos,
    const vector<vec3f>& norm, const vector<vec4f>& weights,
    const vector<vec4i>& joints, const vector<frame3f>& xforms);
// Apply skinning as specified in Khronos glTF.
tuple<vector<vec3f>, vector<vec3f>> compute_matrix_skinning(
    const vector<vec3f>& pos, const vector<vec3f>& norm,
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
// and face ids and filled vectors for pos, norm and texcoord.
void convert_face_varying(vector<vec4i>& qquads, vector<vec3f>& qpos,
    vector<vec3f>& qnorm, vector<vec2f>& qtexcoord, vector<vec4f>& qcolor,
    const vector<vec4i>& quads_pos, const vector<vec4i>& quads_norm,
    const vector<vec4i>& quads_texcoord, const vector<vec4i>& quads_color,
    const vector<vec3f>& pos, const vector<vec3f>& norm,
    const vector<vec2f>& texcoord, const vector<vec4f>& color);

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
    const vector<vec3f>& pos, float threshold);
tuple<vector<vec3i>, vector<vec3f>> weld_triangles(
    const vector<vec3i>& triangles, const vector<vec3f>& pos, float threshold);
tuple<vector<vec4i>, vector<vec3f>> weld_quads(
    const vector<vec4i>& quads, const vector<vec3f>& pos, float threshold);

// Pick a point in a point set uniformly.
inline int sample_points(int npoints, float re) {
    return sample_index(npoints, re);
}
inline vector<float> sample_points_cdf(int npoints) {
    auto cdf = vector<float>(npoints);
    for (auto i = 0; i < cdf.size(); i++) cdf[i] = 1 + (i ? cdf[i - 1] : 0);
    return cdf;
}
inline int sample_points(const vector<float>& cdf, float re) {
    return sample_discrete(cdf, re);
}

// Pick a point on lines uniformly.
inline vector<float> sample_lines_cdf(
    const vector<vec2i>& lines, const vector<vec3f>& pos) {
    auto cdf = vector<float>(lines.size());
    for (auto i = 0; i < cdf.size(); i++) {
        auto l = lines[i];
        auto w = line_length(pos[l.x], pos[l.y]);
        cdf[i] = w + (i ? cdf[i - 1] : 0);
    }
    return cdf;
}
inline tuple<int, float> sample_lines(
    const vector<float>& cdf, float re, float ru) {
    return {sample_discrete(cdf, re), ru};
}

// Pick a point on a triangle mesh uniformly.
inline vector<float> sample_triangles_cdf(
    const vector<vec3i>& triangles, const vector<vec3f>& pos) {
    auto cdf = vector<float>(triangles.size());
    for (auto i = 0; i < cdf.size(); i++) {
        auto t = triangles[i];
        auto w = triangle_area(pos[t.x], pos[t.y], pos[t.z]);
        cdf[i] = w + (i ? cdf[i - 1] : 0);
    }
    return cdf;
}
inline tuple<int, vec2f> sample_triangles(
    const vector<float>& cdf, float re, const vec2f& ruv) {
    return {sample_discrete(cdf, re), sample_triangle(ruv)};
}

// Pick a point on a quad mesh uniformly.
inline vector<float> sample_quads_cdf(
    const vector<vec4i>& quads, const vector<vec3f>& pos) {
    auto cdf = vector<float>(quads.size());
    for (auto i = 0; i < cdf.size(); i++) {
        auto q = quads[i];
        auto w = quad_area(pos[q.x], pos[q.y], pos[q.z], pos[q.w]);
        cdf[i] = w + (i ? cdf[i - 1] : 0);
    }
    return cdf;
}
inline tuple<int, vec2f> sample_quads(
    const vector<float>& cdf, float re, const vec2f& ruv) {
    return {sample_discrete(cdf, re), ruv};
}
inline tuple<int, vec2f> sample_quads(const vector<vec4i>& quads,
    const vector<float>& cdf, float re, const vec2f& ruv) {
    auto ei = sample_discrete(cdf, re);
    if (quads[ei].z == quads[ei].w) {
        return {ei, sample_triangle(ruv)};
    } else {
        return {ei, ruv};
    }
}

// Samples a set of points over a triangle mesh uniformly. Returns pos, norm
// and texcoord of the sampled points.
tuple<vector<vec3f>, vector<vec3f>, vector<vec2f>> sample_triangles_points(
    const vector<vec3i>& triangles, const vector<vec3f>& pos,
    const vector<vec3f>& norm, const vector<vec2f>& texcoord, int npoints,
    int seed = 7);

}  // namespace ygl

// -----------------------------------------------------------------------------
// RAY INTERSECTION AND CLOSEST POINT FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Intersect a ray with a point (approximate).
// Based on http://geomalgorithms.com/a02-lines.html.
bool intersect_point(
    const ray3f& ray, const vec3f& p, float r, float& dist, vec2f& uv);

// Intersect a ray with a line (approximate).
// Based on http://geomalgorithms.com/a05-intersect-1.html and
// http://geomalgorithms.com/a07-distance.html#
//     dist3D_Segment_to_Segment
bool intersect_line(const ray3f& ray, const vec3f& p0, const vec3f& p1,
    float r0, float r1, float& dist, vec2f& uv);

// Intersect a ray with a triangle.
bool intersect_triangle(const ray3f& ray, const vec3f& p0, const vec3f& p1,
    const vec3f& p2, float& dist, vec2f& uv);

// Intersect a ray with a quad represented as two triangles (0,1,3) and
// (2,3,1), with the uv coordinates of the second triangle corrected by u =
// 1-u' and v = 1-v' to produce a quad parametrization where u and v go from 0
// to 1. This is equivalent to Intel's Embree.
bool intersect_quad(const ray3f& ray, const vec3f& p0, const vec3f& p1,
    const vec3f& p2, const vec3f& p3, float& dist, vec2f& uv);

// Intersect a ray with a axis-aligned bounding box.
bool intersect_bbox(const ray3f& ray, const bbox3f& bbox);

// Intersect a ray with a axis-aligned bounding box, implemented as
// "Robust BVH Ray Traversal" by T. Ize published at
// http://jcgt.org/published/0002/02/02/paper.pdf
bool intersect_bbox(const ray3f& ray, const vec3f& ray_dinv,
    const vec3i& ray_dsign, const bbox3f& bbox);

// Check if a point overlaps a position within a max distance.
bool overlap_point(const vec3f& pos, float dist_max, const vec3f& p0, float r0,
    float& dist, vec2f& uv);

// Find closest line point to a position.
float closestuv_line(const vec3f& pos, const vec3f& p0, const vec3f& p1);

// Check if a line overlaps a position within a max distance.
bool overlap_line(const vec3f& pos, float dist_max, const vec3f& p0,
    const vec3f& p1, float r0, float r1, float& dist, vec2f& uv);

// Find closest triangle point to a position.
vec2f closestuv_triangle(
    const vec3f& pos, const vec3f& p0, const vec3f& p1, const vec3f& p2);

// Check if a triangle overlaps a position within a max distance.
bool overlap_triangle(const vec3f& pos, float dist_max, const vec3f& p0,
    const vec3f& p1, const vec3f& p2, float r0, float r1, float r2, float& dist,
    vec2f& uv);

// Check if a quad overlaps a position within a max distance.
bool overlap_quad(const vec3f& pos, float dist_max, const vec3f& p0,
    const vec3f& p1, const vec3f& p2, const vec3f& p3, float r0, float r1,
    float r2, float r3, float& dist, vec2f& uv);

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
// nodes. See bvh_tree for more details.
struct bvh_node {
    bbox3f   bbox;
    uint16_t num_primitives;
    bool     is_internal;
    uint8_t  split_axis;
    uint32_t primitive_ids[bvh_max_prims];
};

// Instance for a scene BVH.
struct bvh_instance {
    frame3f frame         = identity_frame3f;
    frame3f frame_inverse = identity_frame3f;
    int     shape_id      = -1;
};

// BVH tree, stored as a node array. The tree structure is encoded using array
// indices instead of pointers, both for speed but also to simplify code.
// BVH nodes indices refer to either the node array, for internal nodes,
// or the primitive arrays, for leaf nodes. BVH trees may contain only one type
// of geometric primitive, like points, lines, triangle or instances of other
// BVHs. To handle multiple primitive types and transformed primitives, build
// a two-level hierarchy with the outer BVH, the scene BVH, containing inner
// BVHs, shape BVHs, each of which of a uniform primitive type.
// To build a BVH, first fill in either the shape or instance data, then
// call `build_bvh()`.
struct bvh_tree {
    // data for shape BVH
    vector<vec3f> positions;
    vector<float> radius;
    vector<int>   points;
    vector<vec2i> lines;
    vector<vec3i> triangles;
    vector<vec4i> quads;

    // data for instance BVH
    vector<bvh_instance> instances;
    vector<bvh_tree*>    shape_bvhs;

    // optional application specific data to key from a pointer to internal ids
    unordered_map<void*, int> instance_ids;
    unordered_map<void*, int> shape_ids;

    // bvh internal nodes
    vector<bvh_node> nodes;

    // Embree opaque data
    void* embree_bvh = nullptr;

    // cleanup
    ~bvh_tree();
};

// Build a BVH from the given set of primitives.
void build_bvh(bvh_tree* bvh, bool high_quality = false);
// Update the node bounds for a shape bvh.
void refit_bvh(bvh_tree* bvh);

// Build a BVH from the given set of primitives.
// Uses Embree if available and requested, otherwise the standard build.
void build_bvh_embree(bvh_tree* bvh, bool high_quality = false);

// Intersect ray with a bvh returning either the first or any intersection
// depending on `find_any`. Returns the ray distance `dist`, the instance
// id `iid`, the shape id `sid`, the shape element index `eid` and the
// shape barycentric coordinates `uv`.
bool intersect_bvh(const bvh_tree* bvh, const ray3f& ray, bool find_any,
    float& dist, int& iid, int& eid, vec2f& uv);

// Find a shape element that overlaps a point within a given distance
// `max_dist`, returning either the closest or any overlap depending on
// `find_any`. Returns the point distance `dist`, the instance id `iid`, the
// shape id `sid`, the shape element index `eid` and the shape barycentric
// coordinates `uv`.
bool overlap_bvh(const bvh_tree* bvh, const vec3f& pos, float max_dist,
    bool find_any, float& dist, int& iid, int& eid, vec2f& uv);

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
};

// Shape data returned by make_fv<shape> functions.
struct make_fvshape_data {
    // Vertex data
    vector<vec3f> positions;
    vector<vec3f> normals;
    vector<vec2f> texturecoords;

    // Faces swith different topology for each data
    vector<vec4i> positions_quads;
    vector<vec4i> normals_quads;
    vector<vec4i> quads_texcoord;
};

// Make examples shapes that are not watertight (besides quads).
// Return (triangles, quads, pos, norm, texcoord)
make_shape_data make_quad(const vec2i& steps, const vec2f& size,
    const vec2f& uvsize, bool as_triangles);
make_shape_data make_quad_stack(const vec3i& steps, const vec3f& size,
    const vec2f& uvsize, bool as_triangles);
make_shape_data make_floor(const vec2i& steps, const vec2f& size,
    const vec2f& uvsize, bool as_triangles);
make_shape_data make_cube(const vec3i& steps, const vec3f& size,
    const vec3f& uvsize, bool as_triangles);
make_shape_data make_cube_rounded(const vec3i& steps, const vec3f& size,
    const vec3f& uvsize, float radius, bool as_triangles);
make_shape_data make_sphere(
    const vec2i& steps, float size, const vec2f& uvsize, bool as_triangles);
make_shape_data make_sphere_cube(
    int steps, float size, float uvsize, bool as_triangles);
make_shape_data make_sphere_flipcap(const vec2i& steps, float size,
    const vec2f& uvsize, const vec2f& zflip, bool as_triangles);
make_shape_data make_disk(
    const vec2i& steps, float size, const vec2f& uvsize, bool as_triangles);
make_shape_data make_disk_quad(
    int steps, float size, float uvsize, bool as_triangles);
make_shape_data make_disk_bulged(
    int steps, float size, float uvsize, float height, bool as_triangles);
make_shape_data make_cylinder_side(const vec2i& steps, const vec2f& size,
    const vec2f& uvsize, bool as_triangles);
make_shape_data make_cylinder(const vec3i& steps, const vec2f& size,
    const vec3f& uvsize, bool as_triangles);
make_shape_data make_cylinder_rounded(const vec3i& steps, const vec2f& size,
    const vec3f& uvsize, float radius, bool as_triangles);
make_shape_data make_geodesic_sphere(
    int tesselation, float size, bool as_triangles);

// Make examples shapes with are watertight (good for subdivs).
// Returns (triangles, quads, pos)
make_shape_data make_suzanne(float size, bool as_triangles);
make_shape_data make_cube(const vec3f& size, bool as_triangles);

// Make facevarying example shapes that are watertight (good for subdivs).
make_fvshape_data make_fvcube(
    const vec3i& steps, const vec3f& size, const vec3f& uvsize);

// Generate lines set along a quad. Returns lines, pos, norm, texcoord, radius.
make_shape_data make_lines(const vec2i& steps, const vec2f& size,
    const vec2f& uvsize, const vec2f& line_radius = {0.001f, 0.001f});

// Make point primitives. Returns points, pos, norm, texcoord, radius.
make_shape_data make_point(float point_radius = 0.001f);
make_shape_data make_points(int num, float uvsize, float point_radius = 0.001f);
make_shape_data make_random_points(int num, const vec3f& size, float uvsize,
    float point_radius = 0.001f, uint64_t seed = 0);

// Make a bezier circle. Returns bezier, pos.
make_shape_data make_bezier_circle(vector<vec4i>& beziers, vector<vec3f>& pos);

// Make a hair ball around a shape.  Returns lines, pos, norm, texcoord, radius.
// length: minimum and maximum length
// rad: minimum and maximum radius from base to tip
// noise: noise added to hair (strength/scale)
// clump: clump added to hair (number/strength)
// rotation: rotation added to hair (angle/strength)
make_shape_data make_hair(const vec2i& steps, const vector<vec3i>& striangles,
    const vector<vec4i>& squads, const vector<vec3f>& spos,
    const vector<vec3f>& snorm, const vector<vec2f>& stexcoord,
    const vec2f& length = {0.1f, 0.1f}, const vec2f& rad = {0.001f, 0.001f},
    const vec2f& noise = zero2f, const vec2f& clump = zero2f,
    const vec2f& rotation = zero2f, int seed = 7);

// Helper to concatenated shape data for non-facevarying shapes.
make_shape_data merge_shape_data(const vector<make_shape_data>& shapes);

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMAGE TYPE
// -----------------------------------------------------------------------------
namespace ygl {

// Image container.
template <typename T>
struct image {
    int       width  = 0;
    int       height = 0;
    vector<T> pixels = {};

    // constructors
    image() : width{0}, height{0}, pixels() {}
    image(int w, int h, const T& v = T{})
        : width{w}, height{h}, pixels(w * h, v) {}
    image(int w, int h, const T* v)
        : width{w}, height{h}, pixels(v, v + w * h) {}
};

// Element access.
template <typename T>
inline T& pixel_at(image<T>& img, int i, int j) {
    return img.pixels[j * img.width + i];
}
template <typename T>
inline const T& pixel_at(const image<T>& img, int i, int j) {
    return img.pixels[j * img.width + i];
}

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
image<vec4f> tonemap_filmic(
    const image<vec4f>& hdr, float exposure, bool filmic, bool srgb);

// Resize an image.
image<vec4f> resize_image(const image<vec4f>& img, int width, int height);

}  // namespace ygl

// -----------------------------------------------------------------------------
// EXAMPLE IMAGES
// -----------------------------------------------------------------------------
namespace ygl {

// Make example images.
image<vec4f> make_grid_image4f(int width, int height, int tile = 8,
    const vec4f& c0 = {0.2f, 0.2f, 0.2f, 1},
    const vec4f& c1 = {0.8f, 0.8f, 0.8f, 1});
image<vec4f> make_checker_image4f(int width, int height, int tile = 8,
    const vec4f& c0 = {0.2f, 0.2f, 0.2f, 1},
    const vec4f& c1 = {0.8f, 0.8f, 0.8f, 1});
image<vec4f> make_bumpdimple_image4f(int width, int height, int tile = 8);
image<vec4f> make_ramp_image4f(int width, int height, const vec4f& c0,
    const vec4f& c1, float srgb = false);
image<vec4f> make_gammaramp_image4f(int width, int height);
image<vec4f> make_uvramp_image4f(int width, int height);
image<vec4f> make_uvgrid_image4f(
    int width, int height, int tile = 8, bool colored = true);

// Comvert a bump map to a normal map.
image<vec4f> bump_to_normal_map(const image<vec4f>& img, float scale = 1);

// Make a sunsky HDR model with sun at theta elevation in [0,pif/2], turbidity
// in [1.7,10] with or without sun.
image<vec4f> make_sunsky_image4f(int width, int height, float thetaSun,
    float turbidity = 3, bool has_sun = false,
    const vec3f& ground_albedo = {0.7f, 0.7f, 0.7f});
// Make an image of multiple lights.
image<vec4f> make_lights_image4f(int width, int height,
    const vec3f& le = {1, 1, 1}, int nlights = 4, float langle = pif / 4,
    float lwidth = pif / 16, float lheight = pif / 16);

// Make a noise image. Wrap works only if both resx and resy are powers of two.
image<vec4f> make_noise_image4f(
    int width, int height, float scale = 1, bool wrap = true);
image<vec4f> make_fbm_image4f(int width, int height, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, int octaves = 6, bool wrap = true);
image<vec4f> make_ridge_image4f(int width, int height, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, float offset = 1.0f,
    int octaves = 6, bool wrap = true);
image<vec4f> make_turbulence_image4f(int width, int height, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, int octaves = 6, bool wrap = true);

}  // namespace ygl

// -----------------------------------------------------------------------------
// ANIMATION UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Find the first keyframe value that is greater than the argument.
inline int eval_keyframed_index(const vector<float>& times, const float& time) {
    for (auto i = 0; i < times.size(); i++)
        if (times[i] > time) return i;
    return (int)times.size();
}

// Evaluates a keyframed value using step interpolation.
template <typename T>
inline T eval_keyframed_step(
    const vector<float>& times, const vector<T>& vals, float time) {
    if (time <= times.front()) return vals.front();
    if (time >= times.back()) return vals.back();
    time     = clamp(time, times.front(), times.back() - 0.001f);
    auto idx = eval_keyframed_index(times, time);
    return vals.at(idx - 1);
}

// Evaluates a keyframed value using linear interpolation.
template <typename T>
inline vec4f eval_keyframed_slerp(
    const vector<float>& times, const vector<vec4f>& vals, float time) {
    if (time <= times.front()) return vals.front();
    if (time >= times.back()) return vals.back();
    time     = clamp(time, times.front(), times.back() - 0.001f);
    auto idx = eval_keyframed_index(times, time);
    auto t   = (time - times.at(idx - 1)) / (times.at(idx) - times.at(idx - 1));
    return slerp(vals.at(idx - 1), vals.at(idx), t);
}

// Evaluates a keyframed value using linear interpolation.
template <typename T>
inline T eval_keyframed_linear(
    const vector<float>& times, const vector<T>& vals, float time) {
    if (time <= times.front()) return vals.front();
    if (time >= times.back()) return vals.back();
    time     = clamp(time, times.front(), times.back() - 0.001f);
    auto idx = eval_keyframed_index(times, time);
    auto t   = (time - times.at(idx - 1)) / (times.at(idx) - times.at(idx - 1));
    return vals.at(idx - 1) * (1 - t) + vals.at(idx) * t;
}

// Evaluates a keyframed value using Bezier interpolation.
template <typename T>
inline T eval_keyframed_bezier(
    const vector<float>& times, const vector<T>& vals, float time) {
    if (time <= times.front()) return vals.front();
    if (time >= times.back()) return vals.back();
    time     = clamp(time, times.front(), times.back() - 0.001f);
    auto idx = eval_keyframed_index(times, time);
    auto t   = (time - times.at(idx - 1)) / (times.at(idx) - times.at(idx - 1));
    return interpolate_bezier(
        vals.at(idx - 3), vals.at(idx - 2), vals.at(idx - 1), vals.at(idx), t);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// VOLUME TYPE
// -----------------------------------------------------------------------------
namespace ygl {

// Volume container.
template <typename T>
struct volume {
    int       width  = 0;
    int       height = 0;
    int       depth  = 0;
    vector<T> voxels = {};

    // constructors
    volume() : width{0}, height{0}, depth{0}, voxels() {}
    volume(int w, int h, int d, const T& v = T{})
        : width{w}, height{h}, depth{d}, voxels(w * h * d, v) {}
    volume(int w, int h, int d, const T* v)
        : width{w}, height{h}, depth{d}, voxels(v, v + w * h * d) {}
};

// Element access
template <typename T>
T& voxel_at(volume<T>& vol, int i, int j, int k) {
    return vol.voxels[k * vol.width * vol.height + j * vol.width + i];
}
template <typename T>
const T& voxel_at(const volume<T>& vol, int i, int j, int k) {
    return vol.voxels[k * vol.width * vol.height + j * vol.width + i];
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// SCENE DATA
// -----------------------------------------------------------------------------
namespace ygl {

// make a simple example volume
volume<float> make_test_volume1f(
    int width, int height, int depth, float scale = 10, float exponent = 6);

}  // namespace ygl

// -----------------------------------------------------------------------------
// SCENE DATA
// -----------------------------------------------------------------------------
namespace ygl {

// forward declaration
struct bvh_tree;

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
    yocto_texture* emission_texture     = nullptr;
    yocto_texture* diffuse_texture      = nullptr;
    yocto_texture* specular_texture     = nullptr;
    yocto_texture* transmission_texture = nullptr;
    yocto_texture* roughness_texture    = nullptr;
    yocto_texture* opacity_texture      = nullptr;
    yocto_texture* occlusion_texture    = nullptr;
    yocto_texture* bump_texture         = nullptr;
    yocto_texture* displacement_texture = nullptr;
    yocto_texture* normal_texture       = nullptr;

    // volume properties
    // albedo = scattering / (absorption + scattering)
    // density = absorption + scattering
    vec3f volume_emission = {0, 0, 0};
    vec3f volume_albedo   = {0, 0, 0};
    vec3f volume_density  = {0, 0, 0};
    float volume_phaseg   = 0;

    // volume textures
    yocto_voltexture* volume_density_texture = nullptr;
};

// Shape data represented as an indexed meshes of elements.
// May contain either points, lines, triangles and quads.
struct yocto_shape {
    string name     = "";
    string filename = "";

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

// Subdivision surface.
struct yocto_surface {
    string name              = "";
    string filename          = "";
    int    subdivision_level = 0;
    bool   catmull_clark     = true;
    bool   compute_normals   = true;

    // primitives for each vertex propoerty
    vector<vec4i> positions_quads     = {};
    vector<vec4i> texturecoords_quads = {};
    vector<vec4i> colors_quads        = {};

    // creases
    vector<vec3i> positions_creases     = {};
    vector<vec3i> texturecoords_creases = {};

    // vertex data
    vector<vec3f> positions     = {};
    vector<vec2f> texturecoords = {};
    vector<vec4f> colors        = {};
};

// Shape instance.
struct yocto_instance {
    string          name     = "";
    frame3f         frame    = identity_frame3f;
    yocto_shape*    shape    = nullptr;
    yocto_material* material = nullptr;
    yocto_surface*  surface  = nullptr;
};

// Environment map.
struct yocto_environment {
    string         name             = "";
    frame3f        frame            = identity_frame3f;
    vec3f          emission         = {0, 0, 0};
    yocto_texture* emission_texture = nullptr;
};

// Node in a transform hierarchy.
struct yocto_scene_node {
    string             name        = "";
    yocto_scene_node*  parent      = nullptr;
    frame3f            local       = identity_frame3f;
    vec3f              translation = {0, 0, 0};
    vec4f              rotation    = {0, 0, 0, 1};
    vec3f              scale       = {1, 1, 1};
    vector<float>      weights     = {};
    yocto_camera*      camera      = nullptr;
    yocto_instance*    instance    = nullptr;
    yocto_environment* environment = nullptr;

    // compute properties
    vector<yocto_scene_node*> children = {};
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
    vector<yocto_scene_node*> node_targets           = {};
};

// Scene comprised an array of objects whose memory is owened by the scene.
// All members are optional,Scene objects (camera, instances, environments)
// have transforms defined internally. A scene can optionally contain a
// node hierarchy where each node might point to a camera, instance or
// environment. In that case, the element transforms are computed from
// the hierarchy. Animation is also optional, with keyframe data that
// updates node transformations only if defined.
struct yocto_scene {
    string                     name         = "";
    vector<yocto_camera*>      cameras      = {};
    vector<yocto_shape*>       shapes       = {};
    vector<yocto_surface*>     surfaces     = {};
    vector<yocto_instance*>    instances    = {};
    vector<yocto_material*>    materials    = {};
    vector<yocto_texture*>     textures     = {};
    vector<yocto_environment*> environments = {};
    vector<yocto_voltexture*>  voltextures  = {};
    vector<yocto_scene_node*>  nodes        = {};
    vector<yocto_animation*>   animations   = {};

    // cleanup
    ~yocto_scene();
};

}  // namespace ygl

// -----------------------------------------------------------------------------
// SCENE UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Print scene statistics.
void print_stats(const yocto_scene* scene);

// Merge scene into one another. Note that the objects are _moved_ from
// merge_from to merged_into, so merge_from will be empty after this function.
void merge_into(const yocto_scene* merge_into, yocto_scene* merge_from);

}  // namespace ygl

// -----------------------------------------------------------------------------
// UPDATES TO COMPUTED PROPERTIES
// -----------------------------------------------------------------------------
namespace ygl {

// Update node transforms.
void update_transforms(
    yocto_scene* scene, float time = 0, const string& anim_group = "");
// Compute animation range.
vec2f compute_animation_range(
    const yocto_scene* scene, const string& anim_group = "");

// Computes shape/scene approximate bounds.
bbox3f compute_bbox(const yocto_shape* shape);
bbox3f compute_bbox(const yocto_scene* scene);

// Generate a distribution for sampling a shape uniformly based on area/length.
vector<float> compute_shape_cdf(const yocto_shape* shape);
// Generate a distribution for sampling an environment texture uniformly
// based on angle and texture intensity.
vector<float> compute_environment_cdf(const yocto_environment* environment);

// Updates/refits bvh.
bvh_tree* build_bvh(
    const yocto_shape* shape, bool high_quality, bool embree = false);
bvh_tree* build_bvh(
    const yocto_scene* scene, bool high_quality, bool embree = false);
void refit_bvh(const yocto_shape* shape, bvh_tree* bvh);
void refit_bvh(const yocto_scene* scene, bvh_tree* bvh);

// Updates tesselation.
void tesselate_subdiv(const yocto_surface* surface, yocto_shape* shape);
void tesselate_subdivs(yocto_scene* scene);

// Add missing names, normals, tangents and hierarchy.
void add_missing_names(yocto_scene* scene);
void add_missing_normals(yocto_scene* scene);
void add_missing_tangent_space(yocto_scene* scene);
void add_missing_materials(yocto_scene* scene);
void add_missing_cameras(yocto_scene* scene);
// Checks for validity of the scene.
vector<string> validate(const yocto_scene* scene, bool skip_textures = false);

// make camera
yocto_camera* make_bbox_camera(const string& name, const bbox3f& bbox,
    const vec2f& film = {0.036f, 0.024f}, float focal = 0.050f);
// make default material
inline yocto_material* make_default_material(const string& name) {
    auto mat     = new yocto_material();
    mat->name    = name;
    mat->diffuse = {0.2f, 0.2f, 0.2f};
    return mat;
}

// Add a sky environment
inline yocto_environment* make_sky_environment(
    const string& name, float sun_angle = pif / 4) {
    auto texture                  = new yocto_texture();
    texture->name                 = name;
    texture->filename             = "textures/" + name + ".hdr";
    texture->hdr_image            = make_sunsky_image4f(1024, 512, sun_angle);
    auto environment              = new yocto_environment();
    environment->name             = name;
    environment->emission         = {1, 1, 1};
    environment->emission_texture = texture;
    return environment;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// EVALUATION OF SCENE PROPERTIES
// -----------------------------------------------------------------------------
namespace ygl {

// Scene intersection. Upron intersection we set the instance pointer,
// the shape element_id and element_uv and the inetrsection distance.
struct scene_intersection {
    yocto_instance* instance   = nullptr;
    int             element_id = 0;
    vec2f           element_uv = zero2f;
    float           distance   = maxf;
};

// Intersects a ray with an instance. The bvh refers is the shape bvh.
scene_intersection intersect_ray(const yocto_instance* instance,
    const bvh_tree* sbvh, const ray3f& ray, bool find_any = false);
// Intersects a ray with the scene.
scene_intersection intersect_ray(const yocto_scene* scene, const bvh_tree* bvh,
    const ray3f& ray, bool find_any = false);

// Shape values interpolated using barycentric coordinates.
vec3f eval_position(const yocto_shape* shape, int ei, const vec2f& uv);
vec3f eval_normal(const yocto_shape* shape, int ei, const vec2f& uv);
vec2f eval_texturecoord(const yocto_shape* shape, int ei, const vec2f& uv);
vec4f eval_color(const yocto_shape* shape, int ei, const vec2f& uv);
float eval_radius(const yocto_shape* shape, int ei, const vec2f& uv);
vec4f eval_tangentspace(const yocto_shape* shape, int ei, const vec2f& uv);
vec3f eval_tangentspace(
    const yocto_shape* shape, int ei, const vec2f& uv, bool& left_handed);
// Shape element values.
vec3f eval_element_normal(const yocto_shape* shape, int ei);
vec4f eval_element_tangentspace(const yocto_shape* shape, int ei);

// Instance values interpolated using barycentric coordinates.
// Handles defaults if data is missing.
vec3f eval_position(const yocto_instance* instance, int ei, const vec2f& uv);
vec3f eval_normal(const yocto_instance* instance, int ei, const vec2f& uv);
vec2f eval_texturecoord(const yocto_instance* instance, int ei, const vec2f& uv);
vec4f eval_color(const yocto_instance* instance, int ei, const vec2f& uv);
float eval_radius(const yocto_instance* instance, int ei, const vec2f& uv);
vec3f eval_tangentspace(
    const yocto_instance* instance, int ei, const vec2f& uv, bool& left_handed);
// Instance element values.
vec3f eval_element_normal(const yocto_instance* instance, int ei);
// Shading normals including material perturbations.
vec3f eval_shading_normal(
    const yocto_instance* instance, int ei, const vec2f& uv, const vec3f& o);

// Environment texture coordinates from the incoming direction.
vec2f eval_texcoord(const yocto_environment* environment, const vec3f& i);
// Evaluate the incoming direction from the uv.
vec3f eval_direction(const yocto_environment* environment, const vec2f& uv);
// Evaluate the environment emission.
vec3f eval_emission(const yocto_environment* environment, const vec3f& i);
// Evaluate all environment emission.
vec3f eval_emission(const yocto_scene* scene, const vec3f& i);

// Evaluate a texture.
vec2i eval_texture_size(const yocto_texture* texture);
vec4f lookup_texture(const yocto_texture* texture, int i, int j);
vec4f eval_texture(const yocto_texture* texture, const vec2f& texcoord);
float lookup_voltexture(const yocto_voltexture* texture, int i, int j, int k);
float eval_voltexture(const yocto_voltexture* texture, const vec3f& texcoord);

// Set and evaluate camera parameters. Setters take zeros as default values.
float eval_camera_fovy(const yocto_camera* camera);
float eval_camera_aspect(const yocto_camera* camera);
void  set_camera_fovy(
     yocto_camera* camera, float fovy, float aspect, float width = 0.036f);
vec2i eval_image_size(const yocto_camera* camera, int yresolution);

// Generates a ray from a camera image coordinate `uv` and lens coordinates
// `luv`.
ray3f eval_camera_ray(
    const yocto_camera* camera, const vec2f& uv, const vec2f& luv);
// Generates a ray from a camera for pixel coordinates `ij`, the image size
// `imsize`, the sub-pixel coordinates `puv` and the lens coordinates `luv` and
// the image resolution `res`.
ray3f eval_camera_ray(const yocto_camera* camera, const vec2i& ij,
    const vec2i& imsize, const vec2f& puv, const vec2f& luv);
// Generates a ray from a camera for pixel index `idx`, the image size
// `imsize`, the sub-pixel coordinates `puv` and the lens coordinates `luv` and
// the image resolution `res`.
ray3f eval_camera_ray(const yocto_camera* camera, int idx, const vec2i& imsize,
    const vec2f& puv, const vec2f& luv);

// Evaluates material parameters: emission, diffuse, specular, transmission,
// roughness and opacity.
vec3f eval_emission(const yocto_instance* instance, int ei, const vec2f& uv);
vec3f eval_diffuse(const yocto_instance* instance, int ei, const vec2f& uv);
vec3f eval_specular(const yocto_instance* instance, int ei, const vec2f& uv);
vec3f eval_transmission(const yocto_instance* instance, int ei, const vec2f& uv);
float eval_roughness(const yocto_instance* instance, int ei, const vec2f& uv);
float eval_opacity(const yocto_instance* instance, int ei, const vec2f& uv);

// Material values packed into a convenience structure.
struct bsdf {
    vec3f kd      = zero3f;  // diffuse
    vec3f ks      = zero3f;  // specular
    vec3f kt      = zero3f;  // transmission
    float rs      = 1;       // roughness
    bool  refract = false;   // whether to use refraction in transmission
};
bsdf eval_bsdf(const yocto_instance* instance, int ei, const vec2f& uv);
bool is_delta_bsdf(const bsdf& f);

// Check volume properties.
bool is_volume_homogeneus(const yocto_material* vol);
bool has_volume_color(const yocto_material* vol);

}  // namespace ygl

// -----------------------------------------------------------------------------
// PATH TRACING
// -----------------------------------------------------------------------------
namespace ygl {

// Default trace seed
const auto trace_default_seed = 961748941;

// Type of tracing algorithm to use
enum struct trace_type {
    path,               // path tracing
    volpath,            // volumetric path tracing
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

const auto trace_type_names = vector<string>{"path", "volpath", "direct",
    "environment", "eyelight", "path_nomis", "path_naive", "direct_nomis",
    "debug_normal", "debug_albedo", "debug_texcoord", "debug_frontfacing",
    "debug_diffuse", "debug_specular", "debug_roughness"};

// Trace options
struct trace_params {
    int        camera_id           = 0;
    int        vertical_resolution = 256;
    trace_type sample_tracer       = trace_type::path;
    int        num_samples         = 256;
    int        max_bounces         = 8;
    float      pixel_clamp         = 100;
    int        samples_per_batch   = 16;
    bool       no_parallel         = false;
    int        preview_ratio       = 8;
    float      display_exposure    = 0;
    bool       display_filmic      = false;
    bool       display_srgb        = true;
    int        random_seed         = trace_default_seed;
};

// Trace lights used during rendering.
struct trace_lights {
    vector<yocto_instance*>                          instances;
    vector<yocto_environment*>                       environments;
    unordered_map<yocto_shape*, vector<float>>       shapes_cdfs;
    unordered_map<yocto_environment*, vector<float>> environment_cdfs;
};

// Trace data used during rendering. Initialize with `make_trace_state()`
struct trace_state {
    image<vec4f> rendered_image = {};
    image<vec4f> display_image  = {};

    // internal data used during rendering
    image<vec4f>     accumulation_buffer      = {};
    image<int>       samples_per_pixel        = {};
    image<rng_state> random_number_generators = {};
    int              current_sample           = 0;
    vector<thread>   async_threads;
    bool             async_stop_flag = false;
};

// Initialize lights.
trace_lights* make_trace_lights(
    const yocto_scene* scene, const trace_params& params);

// Initialize state of the renderer.
trace_state* make_trace_state(
    const yocto_scene* scene, const trace_params& params);

// Progressively compute an image by calling trace_samples multiple times.
image<vec4f> trace_image4f(const yocto_scene* scene, const bvh_tree* bvh,
    const trace_lights* lights, const trace_params& params);

// Progressively compute an image by calling trace_samples multiple times.
// Start with an empty state and then successively call this function to
// render the next batch of samples.
bool trace_samples(trace_state* state, const yocto_scene* scene,
    const bvh_tree* bvh, const trace_lights* lights, const trace_params& params);

// Starts an anyncrhounous renderer. The function will keep a reference to
// params.
void trace_async_start(trace_state* state, const yocto_scene* scene,
    const bvh_tree* bvh, const trace_lights* lights, const trace_params& params);
// Stop the asynchronous renderer.
void trace_async_stop(trace_state* state);

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
float eval_ggx(float rs, float ndh, float ndi, float ndo);
// Sample the GGX distribution.
vec3f sample_ggx(float rs, const vec2f& rn);
float sample_ggx_pdf(float rs, float ndh);

// Evaluates the GGX distribution and geometric term.
float eval_ggx_dist(float rs, const vec3f& n, const vec3f& h);
float eval_ggx_sm(float rs, const vec3f& n, const vec3f& o, const vec3f& i);

// Evaluate and sample volume phase function.
vec3f sample_phase_function(float vg, const vec2f& u);
float eval_phase_function(float cos_theta, float vg);

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR MATRICES
// -----------------------------------------------------------------------------
namespace ygl {

// Matrix diagonals and transposes.
template <typename T>
inline mat<T, 2, 2> transpose(const mat<T, 2, 2>& a) {
    return {{a.x.x, a.y.x}, {a.x.y, a.y.y}};
}
template <typename T>
inline mat<T, 3, 3> transpose(const mat<T, 3, 3>& a) {
    return {
        {a.x.x, a.y.x, a.z.x},
        {a.x.y, a.y.y, a.z.y},
        {a.x.z, a.y.z, a.z.z},
    };
}
template <typename T>
inline mat<T, 4, 4> transpose(const mat<T, 4, 4>& a) {
    return {
        {a.x.x, a.y.x, a.z.x, a.w.x},
        {a.x.y, a.y.y, a.z.y, a.w.y},
        {a.x.z, a.y.z, a.z.z, a.w.z},
        {a.x.w, a.y.w, a.z.w, a.w.w},
    };
}

// Matrix adjugates, determinant and inverses.
template <typename T>
inline mat<T, 2, 2> adjugate(const mat<T, 2, 2>& a) {
    return {{a.y.y, -a.x.y}, {-a.y.x, a.x.x}};
}
template <typename T>
inline mat<T, 3, 3> adjugate(const mat<T, 3, 3>& a) {
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
inline mat<T, 4, 4> adjugate(const mat<T, 4, 4>& a) {
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
inline T determinant(const mat<T, 2, 2>& a) {
    return a.x.x * a.y.y - a.x.y * a.y.x;
}
template <typename T>
inline T determinant(const mat<T, 3, 3>& a) {
    return a.x.x * (a.y.y * a.z.z - a.z.y * a.y.z) +
           a.x.y * (a.y.z * a.z.x - a.z.z * a.y.x) +
           a.x.z * (a.y.x * a.z.y - a.z.x * a.y.y);
}
template <typename T>
inline T determinant(const mat<T, 4, 4>& a) {
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
    if (pan.x || pan.y) { frame.o += frame.x * pan.x + frame.y * pan.y; }
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

#endif
