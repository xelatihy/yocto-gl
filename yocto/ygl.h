//
// # Yocto/GL: Tiny C++ Library for Physically-Based Rendering
//
//
// Yocto/GL is a collection of utilities for creating a simple path tracer.
// The list of utilities is described below.
//
//
// ## Small Vectors and Matrices, Frames, Bounding Boxes and Transforms
//
// We provide common operations for small vectors and matrices typically used
// in graphics. In particular, we support 2-4 dimensional vectors `vec<T, 2>`,
// `vec<T, 3>`, `vec<T, 4>`.
//
// We support 2-4 dimensional generic matrices `mat<T, 2>`, `mat<T, 3>`,
// `mat<T, 4>`, with matrix-matrix and matrix-vector products, transposes and
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
// expansion operations for points and otehr bboxes. We provide operations to
// compute bounds for points, lines, triangles and quads.
//
// For both matrices and frames we support transform operations for points,
// vectors and directions (`trasform_point()`, `trasform_vector()`,
// `trasform_direction()`). For frames we also the support inverse operations
// (`transform_xxx_inverse()`). Transform matrices and frames can be
// constructed from basic translation, rotation and scaling, e.g. with
// `translation_mat<T, 4>()` or `translation_frame<T, 3>()` respectively, etc.
// For rotation we support axis-angle and quaternions, with slerp.
//
//
// ## Geometry functions
//
// The library supports basic geomtry functions such as computing
// line/triangle/quad normals and areas, picking points on triangles
// and the like. In these functions triangles are parametrized with us written
// w.r.t the (v1-v0) and (v2-v0) axis respetively. Quads are internally handled
// as pairs of two triangles v0,v1,v3 and v2,v3,v1, with the u/v coordinates
// of the second triangle corrected as 1-u and 1-v to produce a quad
// parametrization where u and v go from 0 to 1. Degenerate quads with v2==v3
// represent triangles correctly, an this convention is used throught the
// library. This is equivalent to Intel's Embree.
//
//
// ## Shape functions
//
// We provide a small number of utlities for shape manipulation for index
// triangle and quad meshes, indexed line and point sets and indexed beziers.
// The utlities collected here are written to support a global illujmination
// rendering and not for generic geometry processing. We support operation for
// shape smoothing, shape subdivision (including Catmull-Clark subdivs), and
// example shape creation.
//
// 1. compute line tangents, and triangle and quad areas and normals with
//    `compute_line_tangents()`, `compute_triangle_normals()`,
//    `compute_quad_normals()`
// 2. interpolate values over primitives with `eval_line()`,
//    `eval_triangle()` and `eval_quad()`
// 3. evaluate Bezier curves and derivatives with `eval_bezier()` and
//    `eval_bezier_derivative()`
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
// 9.  samnple a could of point over a surface with `sample_triangles_points()`
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
// ionstancing for large scenes and easy BVH refitting for interactive
// applictions.
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
// color utilitis, example image creation, tonempping, image resizing, and
// sunsky procedural images. Yocto/Image is written to support the need of a
// minimal, but fully-featured, global illumination renderer, rather than the
// need of generic image editing.
//
// 0. load and save image with Yocto/GLIO
// 1. create images with `image<T>` data structure
// 2. resize images with `resize_image()`
// 3. tonemap images with `tonemap_image()`
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
// 3. use `update_bbox()` to compute element bounds
// 4. can merge scene together with `merge_into()`
// 5. make scene elements with `make_XXX()` functions
// 6. make procedural elements and scenes with `make_proc_XXX()` functions
// 7. for ray-intersection and closest point queries, a BVH can be created with
//    `update_bvh()` and refit with `refit_bvh()`
// 8. compute interpolated values over scene elements with `eval_XXX()`
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
// env is used.
//
// 1. prepare the scene for tracing
//    - build the ray-tracing acceleration structure with `update_bvh()`
//     - prepare lights for rendering with `update_lights()`
// 2. create the inmage buffer and random number generators `make_trace_rngs()`
// 3. render blocks of samples with `trace_samples()`
// 4. you can also start an asynchronous renderer with `trace_asynch_start()`
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

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <algorithm>  // for std::upper_bound
#include <cctype>
#include <cfloat>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>  // for std::hash
#include <iostream>
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

using std::acos;
using std::asin;
using std::atan;
using std::atan2;
using std::cos;
using std::exp;
using std::fabs;
using std::floor;
using std::isfinite;
using std::log;
using std::pow;
using std::round;
using std::sin;
using std::sqrt;
using std::tan;

using byte = unsigned char;
using uint = unsigned int;

// const auto pi_d = 3.14159265358979323846;
const auto pi = 3.14159265f;
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
template <typename T, int N>
struct vec;

// Small size vectors.
template <typename T>
struct vec<T, 2> {
    // constructurs
    vec() : x{0}, y{0} {}
    vec(T x_, T y_) : x{x_}, y{y_} {}

    // elements
    T x = 0;
    T y = 0;
};
template <typename T>
struct vec<T, 3> {
    // constructurs
    vec() : x{0}, y{0}, z{0} {}
    vec(T x_, T y_, T z_) : x{x_}, y{y_}, z{z_} {}

    // elements
    T x = 0;
    T y = 0;
    T z = 0;
};
template <typename T>
struct vec<T, 4> {
    // constructurs
    vec() : x{0}, y{0}, z{0}, w{0} {}
    vec(T x_, T y_, T z_, T w_) : x{x_}, y{y_}, z{z_}, w{w_} {}

    // elements
    T x = 0;
    T y = 0;
    T z = 0;
    T w = 0;
};

// Type aliases
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

// Access xyz component of a vec4 typically used for color operation.
template <typename T>
inline const vec<T, 3>& xyz(const vec<T, 4>& a) {
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
template <typename T>
inline bool operator==(const vec<T, 3>& a, const vec<T, 3>& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}
template <typename T>
inline bool operator!=(const vec<T, 3>& a, const vec<T, 3>& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z;
}
template <typename T>
inline bool operator==(const vec<T, 4>& a, const vec<T, 4>& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}
template <typename T>
inline bool operator!=(const vec<T, 4>& a, const vec<T, 4>& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w;
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
template <typename T>
inline vec<T, 2> operator-(const vec<T, 2>& a, const vec<T, 2>& b) {
    return {a.x - b.x, a.y - b.y};
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
template <typename T>
inline vec<T, 3> operator-(const vec<T, 3>& a, const vec<T, 3>& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
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
template <typename T>
inline vec<T, 4> operator-(const vec<T, 4>& a, const vec<T, 4>& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w};
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
template <typename T, int N>
inline vec<T, N>& operator-=(vec<T, N>& a, const vec<T, N>& b) {
    return a = a - b;
}
template <typename T, int N, typename T1>
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
inline float dot(const vec<T, 2>& a, const vec<T, 2>& b) {
    return a.x * b.x + a.y * b.y;
}
template <typename T>
inline float dot(const vec<T, 3>& a, const vec<T, 3>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
template <typename T>
inline float dot(const vec<T, 4>& a, const vec<T, 4>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}
template <typename T>
inline float cross(const vec<T, 2>& a, const vec<T, 2>& b) {
    return a.x * b.y - a.y * b.x;
}
template <typename T>
inline vec<T, 3> cross(const vec<T, 3>& a, const vec<T, 3>& b) {
    return {
        a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
template <typename T, int N>
inline T length(const vec<T, N>& a) {
    return sqrt(dot(a, a));
}
template <typename T, int N>
inline T length_sqr(const vec<T, N>& a) {
    return dot(a, a);
}
template <typename T, int N>
inline vec<T, N> normalize(const vec<T, N>& a) {
    return length(a) ? a / length(a) : a;
}
template <typename T, int N>
inline T adot(const vec<T, N>& a, const vec<T, N>& b) {
    return fabs(dot(a, b));
}

// Vecror angles and slerps.
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
        d = -d;
    }
    if (d > (T)0.9995) return normalize(an + u * (bn - an));
    auto th = acos(clamp(d, (T)-1, (T)1));
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
    auto k = 1 - eta * eta * max((T)0, (T)1 - dot(n, w) * dot(n, w));
    if (k < 0) return vec<T, 3>();  // tir
    return -w * eta + (eta * dot(n, w) - sqrt(k)) * n;
}

// Max element and clamp.
template <typename T, typename T1>
inline vec<T, 2> clamp(const vec<T, 2>& x, T1 min, T1 max) {
    return {clamp(x.x, min, max), clamp(x.y, min, max)};
}
template <typename T, typename T1>
inline vec<T, 3> clamp(const vec<T, 3>& x, T1 min, T1 max) {
    return {clamp(x.x, min, max), clamp(x.y, min, max), clamp(x.z, min, max)};
}
template <typename T, typename T1>
inline vec<T, 4> clamp(const vec<T, 4>& x, T1 min, T1 max) {
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

// Quaternion operatons represented as xi + yj + zk + w
const auto identity_quat4f = vec<float, 4>{0, 0, 0, 1};
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
    return quat_conjugate(a) / length_sqr(a);
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const vec<T, 2>& a) {
    return os << a.x << " " << a.y;
}
template <typename T>
std::ostream& operator<<(std::ostream& os, const vec<T, 3>& a) {
    return os << a.x << " " << a.y << " " << a.z;
}
template <typename T>
std::ostream& operator<<(std::ostream& os, const vec<T, 4>& a) {
    return os << a.x << " " << a.y << " " << a.z << " " << a.w;
}
template <typename T>
std::istream& operator>>(std::istream& is, vec<T, 2>& a) {
    return is >> a.x >> a.y;
}
template <typename T>
std::istream& operator>>(std::istream& is, vec<T, 3>& a) {
    return is >> a.x >> a.y >> a.z;
}
template <typename T>
std::istream& operator>>(std::istream& is, vec<T, 4>& a) {
    return is >> a.x >> a.y >> a.z >> a.w;
}

}  // namespace ygl

namespace std {

// Hash functor for vector for use with unordered_map
template <typename T, int N>
struct hash<ygl::vec<T, N>> {
    size_t operator()(const ygl::vec<T, N>& v) const {
        auto vh = hash<T>();
        auto h = (size_t)0;
        for (auto i = 0; i < N; i++)
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
template <typename T, int N>
struct mat;

// Small Fixed-size square matrices stored in column major format.
template <typename T>
struct mat<T, 2> {
    // constructors
    mat() : x{1, 0}, y{0, 1} {}
    mat(const vec<T, 2>& x_, const vec<T, 2>& y_) : x{x_}, y{y_} {}

    // elements
    vec<T, 2> x = {1, 0};
    vec<T, 2> y = {0, 1};
};
template <typename T>
struct mat<T, 3> {
    // constructors
    mat() : x{1, 0, 0}, y{0, 1, 0}, z{0, 0, 1} {}
    mat(const vec<T, 3>& x_, const vec<T, 3>& y_, const vec<T, 3>& z_)
        : x{x_}, y{y_}, z{z_} {}

    // elements
    vec<T, 3> x = {1, 0, 0};
    vec<T, 3> y = {0, 1, 0};
    vec<T, 3> z = {0, 0, 1};
};
template <typename T>
struct mat<T, 4> {
    // constructors
    mat() : x{1, 0, 0, 0}, y{0, 1, 0, 0}, z{0, 0, 1, 0}, w{0, 0, 0, 1} {}
    mat(const vec<T, 4>& x_, const vec<T, 4>& y_, const vec<T, 4>& z_,
        const vec<T, 4>& w_)
        : x{x_}, y{y_}, z{z_}, w{w_} {}

    // elements
    vec<T, 4> x = {1, 0, 0, 0};
    vec<T, 4> y = {0, 1, 0, 0};
    vec<T, 4> z = {0, 0, 1, 0};
    vec<T, 4> w = {0, 0, 0, 1};
};

// Type aliases.
using mat2f = mat<float, 2>;
using mat3f = mat<float, 3>;
using mat4f = mat<float, 4>;

// Identity matrices constants.
const auto identity_mat2f = mat2f();
const auto identity_mat3f = mat3f();
const auto identity_mat4f = mat4f();

// Matrix comparisons.
template <typename T>
inline bool operator==(const mat<T, 2>& a, const mat<T, 2>& b) {
    return a.x == b.x && a.y == b.y;
}
template <typename T>
inline bool operator!=(const mat<T, 2>& a, const mat<T, 2>& b) {
    return !(a == b);
}
template <typename T>
inline bool operator==(const mat<T, 3>& a, const mat<T, 3>& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}
template <typename T>
inline bool operator!=(const mat<T, 3>& a, const mat<T, 3>& b) {
    return !(a == b);
}
template <typename T>
inline bool operator==(const mat<T, 4>& a, const mat<T, 4>& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}
template <typename T>
inline bool operator!=(const mat<T, 4>& a, const mat<T, 4>& b) {
    return !(a == b);
}

// Matrix operations.
template <typename T>
inline mat<T, 2> operator+(const mat<T, 2>& a, const mat<T, 2>& b) {
    return {a.x + b.x, a.y + b.y};
}
template <typename T, typename T1>
inline mat<T, 2> operator*(const mat<T, 2>& a, T1 b) {
    return {a.x * b, a.y * b};
}
template <typename T, typename T1>
inline mat<T, 2> operator/(const mat<T, 2>& a, T1 b) {
    return {a.x / b, a.y / b};
}
template <typename T>
inline vec<T, 2> operator*(const mat<T, 2>& a, const vec<T, 2>& b) {
    return a.x * b.x + a.y * b.y;
}
template <typename T>
inline vec<T, 2> operator*(const vec<T, 2>& a, const mat<T, 2>& b) {
    return {dot(a, b.x), dot(a, b.y)};
}
template <typename T>
inline mat<T, 2> operator*(const mat<T, 2>& a, const mat<T, 2>& b) {
    return {a * b.x, a * b.y};
}

// Matrix operations.
template <typename T>
inline mat<T, 3> operator+(const mat<T, 3>& a, const mat<T, 3>& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
template <typename T, typename T1>
inline mat<T, 3> operator*(const mat<T, 3>& a, T1 b) {
    return {a.x * b, a.y * b, a.z * b};
}
template <typename T, typename T1>
inline mat<T, 3> operator/(const mat<T, 3>& a, T1 b) {
    return {a.x / b, a.y / b, a.z / b};
}
template <typename T>
inline vec<T, 3> operator*(const mat<T, 3>& a, const vec<T, 3>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
template <typename T>
inline vec<T, 3> operator*(const vec<T, 3>& a, const mat<T, 3>& b) {
    return {dot(a, b.x), dot(a, b.y), dot(a, b.z)};
}
template <typename T>
inline mat<T, 3> operator*(const mat<T, 3>& a, const mat<T, 3>& b) {
    return {a * b.x, a * b.y, a * b.z};
}

// Matrix operations.
template <typename T>
inline mat<T, 4> operator+(const mat<T, 4>& a, const mat<T, 4>& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}
template <typename T, typename T1>
inline mat<T, 4> operator*(const mat<T, 4>& a, T1 b) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
}
template <typename T>
inline vec<T, 4> operator*(const mat<T, 4>& a, const vec<T, 4>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}
template <typename T>
inline vec<T, 4> operator*(const vec<T, 4>& a, const mat<T, 4>& b) {
    return {dot(a, b.x), dot(a, b.y), dot(a, b.z), dot(a, b.w)};
}
template <typename T>
inline mat<T, 4> operator*(const mat<T, 4>& a, const mat<T, 4>& b) {
    return {a * b.x, a * b.y, a * b.z, a * b.w};
}

// Matrix assignments.
template <typename T, int N>
inline mat<T, N>& operator+=(mat<T, N>& a, const mat<T, N>& b) {
    return a = a + b;
}
template <typename T, int N>
inline mat<T, N>& operator*=(mat<T, N>& a, const mat<T, N>& b) {
    return a = a * b;
}
template <typename T, int N, typename T1>
inline mat<T, N>& operator*=(mat<T, N>& a, T1 b) {
    return a = a * b;
}

// Matrix diagonals and transposes.
template <typename T>
inline vec<T, 2> diagonal(const mat<T, 2>& a) {
    return {a.x.x, a.y.y};
}
template <typename T>
inline vec<T, 3> diagonal(const mat<T, 3>& a) {
    return {a.x.x, a.y.y, a.z.z};
}
template <typename T>
inline vec<T, 4> diagonal(const mat<T, 4>& a) {
    return {a.x.x, a.y.y, a.z.z, a.w.w};
}
template <typename T>
inline mat<T, 2> transpose(const mat<T, 2>& a) {
    return {{a.x.x, a.y.x}, {a.x.y, a.y.y}};
}
template <typename T>
inline mat<T, 3> transpose(const mat<T, 3>& a) {
    return {
        {a.x.x, a.y.x, a.z.x}, {a.x.y, a.y.y, a.z.y}, {a.x.z, a.y.z, a.z.z}};
}
template <typename T>
inline mat<T, 4> transpose(const mat<T, 4>& a) {
    return {{a.x.x, a.y.x, a.z.x, a.w.x}, {a.x.y, a.y.y, a.z.y, a.w.y},
        {a.x.z, a.y.z, a.z.z, a.w.z}, {a.x.w, a.y.w, a.z.w, a.w.w}};
}

// Matrix adjugates, determinant and inverses.
template <typename T>
inline mat<T, 2> adjugate(const mat<T, 2>& a);
template <typename T>
inline mat<T, 3> adjugate(const mat<T, 3>& a);
template <typename T>
inline mat<T, 4> adjugate(const mat<T, 4>& a);
template <typename T>
inline T determinant(const mat<T, 2>& a);
template <typename T>
inline T determinant(const mat<T, 3>& a);
template <typename T>
inline T determinant(const mat<T, 4>& a);
template <typename T>
inline mat<T, 2> inverse(const mat<T, 2>& a) {
    return adjugate(a) * (1 / determinant(a));
}
template <typename T>
inline mat<T, 3> inverse(const mat<T, 3>& a) {
    return adjugate(a) * (1 / determinant(a));
}
template <typename T>
inline mat<T, 4> inverse(const mat<T, 4>& a) {
    return adjugate(a) * (1 / determinant(a));
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const mat<T, 2>& a) {
    return os << a.x << " " << a.y;
}
template <typename T>
std::ostream& operator<<(std::ostream& os, const mat<T, 3>& a) {
    return os << a.x << " " << a.y << " " << a.z;
}
template <typename T>
std::ostream& operator<<(std::ostream& os, const mat<T, 4>& a) {
    return os << a.x << " " << a.y << " " << a.z << " " << a.w;
}
template <typename T>
std::istream& operator>>(std::istream& is, mat<T, 2>& a) {
    return is >> a.x >> a.y;
}
template <typename T>
std::istream& operator>>(std::istream& is, mat<T, 3>& a) {
    return is >> a.x >> a.y >> a.z;
}
template <typename T>
std::istream& operator>>(std::istream& is, mat<T, 4>& a) {
    return is >> a.x >> a.y >> a.z >> a.w;
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
    // constructors
    frame() : x{1, 0}, y{0, 1}, o{0, 0} {}
    frame(const vec<T, 2>& x_, const vec<T, 2>& y_, const vec<T, 2>& o_)
        : x{x_}, y{y_}, o{o_} {}

    // elements
    vec<T, 2> x = {1, 0};
    vec<T, 2> y = {0, 1};
    vec<T, 2> o = {0, 0};
};
template <typename T>
struct frame<T, 3> {
    // constructors
    frame() : x{1, 0, 0}, y{0, 1, 0}, z{0, 0, 1}, o{0, 0, 0} {}
    frame(const vec<T, 3>& x_, const vec<T, 3>& y_, const vec<T, 3>& z_,
        const vec<T, 3>& o_)
        : x{x_}, y{y_}, z{z_}, o{o_} {}

    // elements
    vec<T, 3> x = {1, 0, 0};
    vec<T, 3> y = {0, 1, 0};
    vec<T, 3> z = {0, 0, 1};
    vec<T, 3> o = {0, 0, 0};
};

// Type aliases
using frame2f = frame<float, 2>;
using frame3f = frame<float, 3>;

// Indentity frames.
const auto identity_frame2f = frame2f{{1, 0}, {0, 1}, {0, 0}};
const auto identity_frame3f =
    frame3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};

// Frame construction from axis.
template <typename T>
inline frame<T, 3> make_frame_fromz(const vec<T, 3>& o, const vec<T, 3>& z_) {
    auto z = normalize(z_);
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
inline mat<T, 4> frame_to_mat(const frame<T, 3>& a) {
    return {{a.x.x, a.x.y, a.x.z, 0}, {a.y.x, a.y.y, a.y.z, 0},
        {a.z.x, a.z.y, a.z.z, 0}, {a.o.x, a.o.y, a.o.z, 1}};
}
template <typename T>
inline frame<T, 3> mat_to_frame(const mat<T, 4>& a) {
    return {{a.x.x, a.x.y, a.x.z}, {a.y.x, a.y.y, a.y.z}, {a.z.x, a.z.y, a.z.z},
        {a.w.x, a.w.y, a.w.z}};
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
    auto rot = mat<T, 2>{a.x, a.y} * mat<T, 2>{b.x, b.y};
    auto pos = mat<T, 2>{a.x, a.y} * b.o + a.o;
    return {rot.x, rot.y, pos};
}
template <typename T>
inline frame<T, 3> operator*(const frame<T, 3>& a, const frame<T, 3>& b) {
    auto rot = mat<T, 3>{a.x, a.y, a.z} * mat<T, 3>{b.x, b.y, b.z};
    auto pos = mat<T, 3>{a.x, a.y, a.z} * b.o + a.o;
    return {rot.x, rot.y, rot.z, pos};
}
// Frame inverse, equivalent to rigid affine inverse.
template <typename T>
inline frame<T, 2> inverse(const frame<T, 2>& a, bool is_rigid = true) {
    auto minv = (is_rigid) ? transpose(mat<T, 2>{a.x, a.y}) :
                             inverse(mat<T, 2>{a.x, a.y});
    return {minv.x, minv.y, -(minv * a.o)};
}
template <typename T>
inline frame<T, 3> inverse(const frame<T, 3>& a, bool is_rigid = true) {
    auto minv = (is_rigid) ? transpose(mat<T, 3>{a.x, a.y, a.z}) :
                             inverse(mat<T, 3>{a.x, a.y, a.z});
    return {minv.x, minv.y, minv.z, -(minv * a.o)};
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const frame<T, 2>& a) {
    return os << a.x << " " << a.y << " " << a.o;
}
template <typename T>
std::ostream& operator<<(std::ostream& os, const frame<T, 3>& a) {
    return os << a.x << " " << a.y << " " << a.z << " " << a.o;
}
template <typename T>
std::istream& operator>>(std::istream& is, frame<T, 2>& a) {
    return is >> a.x >> a.y >> a.o;
}
template <typename T>
std::istream& operator>>(std::istream& is, frame<T, 3>& a) {
    return is >> a.x >> a.y >> a.z >> a.o;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// AXIS ALIGNED BOUNDING BOXES
// -----------------------------------------------------------------------------
namespace ygl {

// Axis aligned bounding box represented as a min/max vector pairs.
template <typename T, int N>
struct bbox;

// Range of values in 1D.
template <typename T>
struct bbox<T, 1> {
    // constructors
    bbox()
        : min{std::numeric_limits<T>::max()}
        , max{std::numeric_limits<T>::lowest()} {}
    bbox(T min_, T max_) : min{min_}, max{max_} {}

    // elements
    T min = std::numeric_limits<T>::max();
    T max = std::numeric_limits<T>::lowest();
};

// Axis aligned bounding box represented as a min/max vector pairs.
template <typename T>
struct bbox<T, 2> {
    // constructors
    bbox()
        : min{std::numeric_limits<T>::max(), std::numeric_limits<T>::max()}
        , max{std::numeric_limits<T>::lowest(),
              std::numeric_limits<T>::lowest()} {}
    bbox(const vec<T, 2>& min_, const vec<T, 2>& max_) : min{min_}, max{max_} {}

    // elements
    vec<T, 2> min = {
        std::numeric_limits<T>::max(), std::numeric_limits<T>::max()};
    vec<T, 2> max = {
        std::numeric_limits<T>::lowest(), std::numeric_limits<T>::lowest()};
};

// Axis aligned bounding box represented as a min/max vector pairs.
template <typename T>
struct bbox<T, 3> {
    // constructors
    bbox()
        : min{std::numeric_limits<T>::max(), std::numeric_limits<T>::max(),
              std::numeric_limits<T>::max()}
        , max{std::numeric_limits<T>::lowest(),
              std::numeric_limits<T>::lowest(),
              std::numeric_limits<T>::lowest()} {}
    bbox(const vec<T, 3>& min_, const vec<T, 3>& max_) : min{min_}, max{max_} {}

    // elements
    vec<T, 3> min = {std::numeric_limits<T>::max(),
        std::numeric_limits<T>::max(), std::numeric_limits<T>::max()};
    vec<T, 3> max = {std::numeric_limits<T>::lowest(),
        std::numeric_limits<T>::lowest(), std::numeric_limits<T>::lowest()};
};

// Axis aligned bounding box represented as a min/max vector pairs.
template <typename T>
struct bbox<T, 4> {
    // constructors
    bbox()
        : min{std::numeric_limits<T>::max(), std::numeric_limits<T>::max(),
              std::numeric_limits<T>::max(), std::numeric_limits<T>::max()}
        , max{std::numeric_limits<T>::lowest(),
              std::numeric_limits<T>::lowest(),
              std::numeric_limits<T>::lowest(),
              std::numeric_limits<T>::lowest()} {}
    bbox(const vec<T, 4>& min_, const vec<T, 4>& max_) : min{min_}, max{max_} {}

    vec<T, 4> min = {std::numeric_limits<T>::max(),
        std::numeric_limits<T>::max(), std::numeric_limits<T>::max(),
        std::numeric_limits<T>::max()};
    vec<T, 4> max = {std::numeric_limits<T>::lowest(),
        std::numeric_limits<T>::lowest(), std::numeric_limits<T>::lowest(),
        std::numeric_limits<T>::lowest()};
};

// Type alias
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
    a.min = {
        min(a.min.x, b.min.x), min(a.min.y, b.min.y), min(a.min.z, b.min.z)};
    a.max = {
        max(a.max.x, b.max.x), max(a.max.y, b.max.y), max(a.max.z, b.max.z)};
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
template <typename T, typename T1>
inline bbox<T, 3> point_bbox(const vec<T, 3>& p, T1 r = 0) {
    auto bbox = ygl::bbox<T, 3>{};
    bbox += p - vec<T, 3>{r, r, r};
    bbox += p + vec<T, 3>{r, r, r};
    return bbox;
}
template <typename T, typename T1>
inline bbox<T, 3> line_bbox(
    const vec<T, 3>& v0, const vec<T, 3>& v1, T1 r0 = 0, T1 r1 = 0) {
    auto bbox = ygl::bbox<T, 3>{};
    bbox += v0 - vec<T, 3>{r0, r0, r0};
    bbox += v0 + vec<T, 3>{r0, r0, r0};
    bbox += v1 - vec<T, 3>{r1, r1, r1};
    bbox += v1 + vec<T, 3>{r1, r1, r1};
    return bbox;
}
template <typename T>
inline bbox<T, 3> triangle_bbox(
    const vec<T, 3>& v0, const vec<T, 3>& v1, const vec<T, 3>& v2) {
    auto bbox = ygl::bbox<T, 3>{};
    bbox += v0;
    bbox += v1;
    bbox += v2;
    return bbox;
}
template <typename T>
inline bbox<T, 3> quad_bbox(const vec<T, 3>& v0, const vec<T, 3>& v1,
    const vec<T, 3>& v2, const vec<T, 3>& v3) {
    auto bbox = ygl::bbox<T, 3>{};
    bbox += v0;
    bbox += v1;
    bbox += v2;
    bbox += v3;
    return bbox;
}
template <typename T>
inline bbox<T, 3> tetrahedron_bbox(const vec<T, 3>& v0, const vec<T, 3>& v1,
    const vec<T, 3>& v2, const vec<T, 3>& v3) {
    auto bbox = ygl::bbox<T, 3>{};
    bbox += v0;
    bbox += v1;
    bbox += v2;
    bbox += v3;
    return bbox;
}

template <typename T, int N>
std::ostream& operator<<(std::ostream& os, const bbox<T, N>& a) {
    return os << a.min << " " << a.max;
}
template <typename T, int N>
std::istream& operator>>(std::istream& is, bbox<T, N>& a) {
    return is >> a.min >> a.max;
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
    // constructors
    ray() : o{0, 0}, d{0, 1}, tmin{0}, tmax{std::numeric_limits<T>::max()} {}
    ray(const vec<T, 2>& o_, const vec<T, 2>& d_, T tmin_, T tmax_)
        : o{o_}, d{d_}, tmin{tmin_}, tmax{tmax_} {}

    // elements
    vec<T, 2> o = {0, 0};
    vec<T, 2> d = {0, 1};
    T tmin = 0;
    T tmax = std::numeric_limits<T>::max();
};

// Rays with origin, direction and min/max t value.
template <typename T>
struct ray<T, 3> {
    // constructors
    ray()
        : o{0, 0, 0}
        , d{0, 0, 1}
        , tmin{0}
        , tmax{std::numeric_limits<T>::max()} {}
    ray(const vec<T, 3>& o_, const vec<T, 3>& d_, T tmin_, T tmax_)
        : o{o_}, d{d_}, tmin{tmin_}, tmax{tmax_} {}

    // elements
    vec<T, 3> o = {0, 0, 0};
    vec<T, 3> d = {0, 0, 1};
    T tmin = 0;
    T tmax = std::numeric_limits<T>::max();
};

// Type aliases
using ray2f = ray<float, 2>;
using ray3f = ray<float, 3>;

// Construct a ray from dirction or segments using a default epsilon.
template <typename T, int N>
inline ray<T, N> make_ray(
    const vec<T, N>& o, const vec<T, N>& d, T eps = (T)1e-4) {
    return ray<T, N>{o, d, eps, std::numeric_limits<T>::max()};
}
template <typename T, int N>
inline ray<T, N> make_segment(
    const vec<T, N>& p1, const vec<T, N>& p2, T eps = (T)1e-4) {
    return ray<T, N>{p1, normalize(p2 - p1), eps, length(p2 - p1) - 2 * eps};
}

template <typename T, int N>
std::ostream& operator<<(std::ostream& os, const ray<T, N>& a) {
    return os << a.o << " " << a.d << a.tmin << a.tmax;
}
template <typename T, int N>
std::istream& operator>>(std::istream& is, ray<T, N>& a) {
    return is >> a.o >> a.d >> a.tmin >> a.tmax;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// TRANSFORMS
// -----------------------------------------------------------------------------
namespace ygl {

// Transforms points, vectors and directions by matrices.
template <typename T>
inline vec<T, 2> transform_point(const mat<T, 3>& a, const vec<T, 2>& b) {
    auto tvb = a * vec<T, 3>{b.x, b.y, 1};
    return vec<T, 2>{tvb.x, tvb.y} / tvb.z;
}
template <typename T>
inline vec<T, 3> transform_point(const mat<T, 4>& a, const vec<T, 3>& b) {
    auto tvb = a * vec<T, 4>{b.x, b.y, b.z, 1};
    return vec<T, 3>{tvb.x, tvb.y, tvb.z} / tvb.w;
}
template <typename T>
inline vec<T, 2> transform_vector(const mat<T, 3>& a, const vec<T, 2>& b) {
    auto tvb = a * vec<T, 3>{b.x, b.y, 0};
    return vec<T, 2>{tvb.x, tvb.y} / tvb.z;
}
template <typename T>
inline vec<T, 3> transform_vector(const mat<T, 3>& a, const vec<T, 3>& b) {
    return a * b;
}
template <typename T>
inline vec<T, 3> transform_vector(const mat<T, 4>& a, const vec<T, 3>& b) {
    auto tvb = a * vec<T, 4>{b.x, b.y, b.z, 0};
    return vec<T, 3>{tvb.x, tvb.y, tvb.z};
}
template <typename T>
inline vec<T, 3> transform_direction(const mat<T, 4>& a, const vec<T, 3>& b) {
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
template <typename T>
inline vec<T, 3> transform_direction(const frame<T, 3>& a, const vec<T, 3>& b) {
    return normalize(transform_vector(a, b));
}

// Transforms rays and bounding boxes by matrices.
template <typename T>
inline ray<T, 3> transform_ray(const frame<T, 3>& a, const ray<T, 3>& b) {
    return {transform_point(a, b.o), transform_vector(a, b.d), b.tmin, b.tmax};
}
template <typename T>
inline ray<T, 3> transform_ray(const mat<T, 4>& a, const ray<T, 3>& b) {
    return {transform_point(a, b.o), transform_vector(a, b.d), b.tmin, b.tmax};
}
template <typename T>
inline bbox<T, 3> transform_bbox(const frame<T, 3>& a, const bbox<T, 3>& b) {
    auto corners = {vec<T, 3>{b.min.x, b.min.y, b.min.z},
        vec<T, 3>{b.min.x, b.min.y, b.max.z},
        vec<T, 3>{b.min.x, b.max.y, b.min.z},
        vec<T, 3>{b.min.x, b.max.y, b.max.z},
        vec<T, 3>{b.max.x, b.min.y, b.min.z},
        vec<T, 3>{b.max.x, b.min.y, b.max.z},
        vec<T, 3>{b.max.x, b.max.y, b.min.z},
        vec<T, 3>{b.max.x, b.max.y, b.max.z}};
    auto xformed = bbox<T, 3>();
    for (auto& corner : corners) xformed += transform_point(a, corner);
    return xformed;
}
template <typename T>
inline bbox<T, 3> transform_bbox(const mat<T, 4>& a, const bbox<T, 3>& b) {
    auto corners = {vec<T, 3>{b.min.x, b.min.y, b.min.z},
        vec<T, 3>{b.min.x, b.min.y, b.max.z},
        vec<T, 3>{b.min.x, b.max.y, b.min.z},
        vec<T, 3>{b.min.x, b.max.y, b.max.z},
        vec<T, 3>{b.max.x, b.min.y, b.min.z},
        vec<T, 3>{b.max.x, b.min.y, b.max.z},
        vec<T, 3>{b.max.x, b.max.y, b.min.z},
        vec<T, 3>{b.max.x, b.max.y, b.max.z}};
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
inline vec<T, 3> transform_point_inverse(
    const frame<T, 3>& a, const vec<T, 3>& b) {
    return {dot(b - a.o, a.x), dot(b - a.o, a.y), dot(b - a.o, a.z)};
}
template <typename T>
inline vec<T, 2> transform_vector_inverse(
    const frame<T, 2>& a, const vec<T, 2>& b) {
    return {dot(b, a.x), dot(b, a.y)};
}
template <typename T>
inline vec<T, 3> transform_vector_inverse(
    const frame<T, 3>& a, const vec<T, 3>& b) {
    return {dot(b, a.x), dot(b, a.y), dot(b, a.z)};
}
template <typename T>
inline vec<T, 3> transform_direction_inverse(
    const frame<T, 3>& a, const vec<T, 3>& b) {
    return normalize(transform_vector_inverse(a, b));
}
template <typename T>
inline ray<T, 3> transform_ray_inverse(
    const frame<T, 3>& a, const ray<T, 3>& b) {
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
inline frame<T, 3> rotation_frame(const mat<T, 3>& rot) {
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
inline mat<T, 4> frustum_mat(T l, T r, T b, T t, T n, T f) {
    return {{2 * n / (r - l), 0, 0, 0}, {0, 2 * n / (t - b), 0, 0},
        {(r + l) / (r - l), (t + b) / (t - b), -(f + n) / (f - n), -1},
        {0, 0, -2 * f * n / (f - n), 0}};
}
template <typename T>
inline mat<T, 4> ortho_mat(T l, T r, T b, T t, T n, T f) {
    return {{2 / (r - l), 0, 0, 0}, {0, 2 / (t - b), 0, 0},
        {0, 0, -2 / (f - n), 0},
        {-(r + l) / (r - l), -(t + b) / (t - b), -(f + n) / (f - n), 1}};
}
template <typename T>
inline mat<T, 4> ortho2d_mat(T left, T right, T bottom, T top) {
    return ortho_mat(left, right, bottom, top, -1, 1);
}
template <typename T>
inline mat<T, 4> ortho_mat(T xmag, T ymag, T near, T far) {
    return {{1 / xmag, 0, 0, 0}, {0, 1 / ymag, 0, 0},
        {0, 0, 2 / (near - far), 0}, {0, 0, (far + near) / (near - far), 1}};
}
template <typename T>
inline mat<T, 4> perspective_mat(T fovy, T aspect, T near, T far) {
    auto tg = tan(fovy / 2);
    return {{1 / (aspect * tg), 0, 0, 0}, {0, 1 / tg, 0, 0},
        {0, 0, (far + near) / (near - far), -1},
        {0, 0, 2 * far * near / (near - far), 0}};
}
template <typename T>
inline mat<T, 4> perspective_mat(T fovy, T aspect, T near) {
    auto tg = tan(fovy / 2);
    return {{1 / (aspect * tg), 0, 0, 0}, {0, 1 / tg, 0, 0}, {0, 0, -1, -1},
        {0, 0, 2 * near, 0}};
}

// Rotation conversions.
template <typename T>
inline std::pair<vec<T, 3>, T> rotation_axisangle(const vec<T, 4>& quat) {
    return {normalize(vec<T, 3>{quat.x, quat.y, quat.z}), 2 * acos(quat.w)};
}
template <typename T, typename T1>
inline vec<T, 4> rotation_quat(const vec<T, 3>& axis, T1 angle) {
    auto len = length(axis);
    if (!len) return {0, 0, 0, 1};
    return vec<T, 4>{sin(angle / 2) * axis.x / len,
        sin(angle / 2) * axis.y / len, sin(angle / 2) * axis.z / len,
        cos(angle / 2)};
}
template <typename T>
inline vec<T, 4> rotation_quat(const vec<T, 4>& axisangle) {
    return rotation_quat(
        vec<T, 3>{axisangle.x, axisangle.y, axisangle.z}, axisangle.w);
}

// Turntable and FPS Camera navigation.
template <typename T, typename T1>
inline void camera_turntable(vec<T, 3>& from, vec<T, 3>& to, vec<T, 3>& up,
    const vec<T, 2>& rotate, T1 dolly, const vec<T, 2>& pan);
template <typename T, typename T1>
inline void camera_turntable(frame<T, 3>& frame, float& focus,
    const vec<T, 2>& rotate, T1 dolly, const vec<T, 2>& pan);
template <typename T>
inline void camera_fps(
    frame<T, 3>& frame, const vec<T, 3>& transl, const vec<T, 2>& rotate);

// Computes the image uv coordinates corresponding to the view parameters.
// Returns negative coordinates if out of the image.
template <typename T>
inline vec<int, 2> get_image_coords(const vec<T, 2>& mouse_pos,
    const frame<T, 2>& frame, const vec<int, 2>& txt_size) {
    // assume an affine without rotation
    auto xyf = (mouse_pos - frame.o) / vec2f{frame.x.x, frame.y.y};
    return vec2i{(int)round(xyf.x + txt_size.x / 2.0f),
        (int)round(xyf.y + txt_size.y / 2.0f)};
}

// Center image and autofit.
template <typename T>
inline void center_image(frame<T, 2>& frame, const vec<int, 2>& imsize,
    const vec<int, 2>& winsize, bool zoom_to_fit) {
    if (zoom_to_fit) {
        frame.x.x = frame.y.y =
            ygl::min(winsize.x / (float)imsize.x, winsize.y / (float)imsize.y);
        frame.o = {(float)winsize.x / 2, (float)winsize.y / 2};
    } else {
        if (winsize.x >= imsize.x * frame.x.x) frame.o.x = winsize.x / 2;
        if (winsize.y >= imsize.y * frame.y.y) frame.o.y = winsize.y / 2;
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
inline rng_state make_rng(uint64_t seed, uint64_t seq = 1) {
    auto rng = rng_state();
    rng.state = 0U;
    rng.inc = (seq << 1u) | 1u;
    advance_rng(rng);
    rng.state += seed;
    advance_rng(rng);
    return rng;
}

// Next random numbers: floats in [0,1), ints in [0,n).
inline int rand1i(rng_state& rng, int n) { return advance_rng(rng) % n; }
inline float rand1f(rng_state& rng) {
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
template <typename T>
inline vec<T, 3> sample_hemisphere(const vec<T, 2>& ruv) {
    auto z = ruv.y;
    auto r = sqrt(1 - z * z);
    auto phi = 2 * pi * ruv.x;
    return {r * cos(phi), r * sin(phi), z};
}
template <typename T>
inline T sample_hemisphere_pdf(const vec<T, 3>& w) {
    return (w.z <= 0) ? 0 : 1 / (2 * pi);
}

// Sample a spherical direction with uniform distribution.
template <typename T>
inline vec<T, 3> sample_sphere(const vec<T, 2>& ruv) {
    auto z = 2 * ruv.y - 1;
    auto r = sqrt(1 - z * z);
    auto phi = 2 * pi * ruv.x;
    return {r * cos(phi), r * sin(phi), z};
}
template <typename T>
inline T sample_sphere_pdf(const vec<T, 3>& w) {
    return 1 / (4 * pi);
}

// Sample spherical coordinates uniformly.
template <typename T>
inline vec<T, 2> sample_spherical(const vec<T, 2>& ruv) {
    // BUG: FIXME this is not uniform at all!!!!
    return {ruv.x, ruv.y};
}
template <typename T>
inline T sample_spherical_pdf(const vec<T, 2>& w) {
    return 1 / (4 * pi);
}

// Sample an hemispherical direction with cosine distribution.
template <typename T>
inline vec<T, 3> sample_hemisphere_cosine(const vec<T, 2>& ruv) {
    auto z = sqrt(ruv.y);
    auto r = sqrt(1 - z * z);
    auto phi = 2 * pi * ruv.x;
    return {r * cos(phi), r * sin(phi), z};
}
template <typename T>
inline T sample_hemisphere_cosine_pdf(const vec<T, 3>& w) {
    return (w.z <= 0) ? 0 : w.z / pi;
}

// Sample an hemispherical direction with cosine power distribution.
template <typename T>
inline vec<T, 3> sample_hemisphere_cospower(T n, const vec<T, 2>& ruv) {
    auto z = pow(ruv.y, 1 / (n + 1));
    auto r = sqrt(1 - z * z);
    auto phi = 2 * pi * ruv.x;
    return {r * cos(phi), r * sin(phi), z};
}
template <typename T>
inline float sample_hemisphere_cospower_pdf(T n, const vec<T, 3>& w) {
    return (w.z <= 0) ? 0 : pow(w.z, n) * (n + 1) / (2 * pi);
}

// Sample a point uniformly on a disk.
template <typename T>
inline vec<T, 3> sample_disk(const vec<T, 2>& ruv) {
    auto r = sqrt(ruv.y);
    auto phi = 2 * pi * ruv.x;
    return {cos(phi) * r, sin(phi) * r, 0};
}
template <typename T>
inline T sample_disk_pdf() {
    return 1 / pi;
}

// Sample a point uniformly on a cylinder, without caps.
template <typename T>
inline vec<T, 3> sample_cylinder(const vec<T, 2>& ruv) {
    auto phi = 2 * pi * ruv.x;
    return {sin(phi), cos(phi), ruv.y * 2 - 1};
}
template <typename T>
inline T sample_cylinder_pdf() {
    return 1 / pi;
}

// Sample a point uniformly on a triangle.
template <typename T>
inline vec<T, 2> sample_triangle(const vec<T, 2>& ruv) {
    return {1 - sqrt(ruv.x), ruv.y * sqrt(ruv.x)};
}
template <typename T>
inline vec<T, 3> sample_triangle(const vec<T, 3>& v0, const vec<T, 3>& v1,
    const vec<T, 3>& v2, const vec<T, 2>& ruv) {
    auto uv = sample_triangle(ruv);
    return v0 * (1 - uv.x - uv.y) + v1 * uv.x + v2 * uv.y;
}
// Pdf for uniform triangle sampling, i.e. triangle area.
template <typename T>
inline T sample_triangle_pdf(
    const vec<T, 3>& v0, const vec<T, 3>& v1, const vec<T, 3>& v2) {
    return 2 / length(cross(v1 - v0, v2 - v0));
}

// Sample an index with uniform distribution.
template <typename T>
inline int sample_index(int size, T r) {
    return clamp((int)(r * size), 0, size - 1);
}
template <typename T>
inline T sample_index_pdf(int size) {
    return (T)1 / (T)size;
}

// Sample a discrete distribution represented by its cdf.
template <typename T>
inline int sample_discrete(const std::vector<T>& cdf, T r) {
    r = clamp(r * cdf.back(), (T)0.0, cdf.back() - (T)0.00001);
    auto idx = (int)(std::upper_bound(cdf.data(), cdf.data() + cdf.size(), r) -
                     cdf.data());
    return clamp(idx, 0, (int)cdf.size() - 1);
}
// Pdf for uniform discrete distribution sampling.
template <typename T>
inline T sample_discrete_pdf(const std::vector<T>& cdf, int idx) {
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
inline float perlin_noise(vec3f p, vec3i wrap = zero3i);
inline float perlin_ridge_noise(vec3f p, float lacunarity = 2.0f,
    float gain = 0.5f, float offset = 1.0f, int octaves = 6,
    vec3i wrap = zero3i);
inline float perlin_fbm_noise(vec3f p, float lacunarity = 2.0f,
    float gain = 0.5f, int octaves = 6, vec3i wrap = zero3i);
inline float perlin_turbulence_noise(vec3f p, float lacunarity = 2.0f,
    float gain = 0.5f, int octaves = 6, vec3i wrap = zero3i);

}  // namespace ygl

// -----------------------------------------------------------------------------
// GEOMETRY UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Line properties.
template <typename T>
inline vec<T, 3> line_tangent(const vec<T, 3>& v0, const vec<T, 3>& v1) {
    return normalize(v1 - v0);
}
template <typename T>
inline T line_length(const vec<T, 3>& v0, const vec<T, 3>& v1) {
    return length(v1 - v0);
}

// Triangle properties.
template <typename T>
inline vec<T, 3> triangle_normal(
    const vec<T, 3>& v0, const vec<T, 3>& v1, const vec<T, 3>& v2) {
    return normalize(cross(v1 - v0, v2 - v0));
}
template <typename T>
inline T triangle_area(
    const vec<T, 3>& v0, const vec<T, 3>& v1, const vec<T, 3>& v2) {
    return length(cross(v1 - v0, v2 - v0)) / 2;
}

// Quad propeties.
template <typename T>
inline vec<T, 3> quad_normal(const vec<T, 3>& v0, const vec<T, 3>& v1,
    const vec<T, 3>& v2, const vec<T, 3>& v3) {
    return normalize(triangle_normal(v0, v1, v3) + triangle_normal(v2, v3, v1));
}
template <typename T>
inline T quad_area(const vec<T, 3>& v0, const vec<T, 3>& v1,
    const vec<T, 3>& v2, const vec<T, 3>& v3) {
    return triangle_area(v0, v1, v3) + triangle_area(v2, v3, v1);
}

// Triangle tangent and bitangent from uv
template <typename T>
inline std::pair<vec<T, 3>, vec<T, 3>> triangle_tangents_fromuv(
    const vec<T, 3>& v0, const vec<T, 3>& v1, const vec<T, 3>& v2,
    const vec<T, 2>& uv0, const vec<T, 2>& uv1, const vec<T, 2>& uv2) {
    // Follows the definition in http://www.terathon.com/code/tangent.html and
    // https://gist.github.com/aras-p/2843984
    // normal points up from texture space
    auto p = v1 - v0;
    auto q = v2 - v0;
    auto s = vec<T, 2>{uv1.x - uv0.x, uv2.x - uv0.x};
    auto t = vec<T, 2>{uv1.y - uv0.y, uv2.y - uv0.y};
    auto div = s.x * t.y - s.y * t.x;

    if (div != 0) {
        auto tu = vec<T, 3>{t.y * p.x - t.x * q.x, t.y * p.y - t.x * q.y,
                      t.y * p.z - t.x * q.z} /
                  div;
        auto tv = vec<T, 3>{s.x * q.x - s.y * p.x, s.x * q.y - s.y * p.y,
                      s.x * q.z - s.y * p.z} /
                  div;
        return {tu, tv};
    } else {
        return {{1, 0, 0}, {0, 1, 0}};
    }
}

// Copies of point value. Here only for completeness.
template <typename TT>
inline TT interpolate_point(const std::vector<TT>& vals, int p) {
    if (vals.empty()) return TT();
    return vals[p];
}

// Interpolates values over a line parametrized from a to b by u. Same as lerp.
template <typename TT, typename T>
inline TT interpolate_line(const TT& v0, const TT& v1, T u) {
    return v0 * (1 - u) + v1 * u;
}
template <typename TT, typename T>
inline TT interpolate_line(
    const std::vector<TT>& vals, const vec<int, 2>& l, T u) {
    if (vals.empty()) return TT();
    return interpolate_line(vals[l.x], vals[l.y], u);
}

// Interpolates values over a triangle parametrized by u and v along the
// (v1-v0) and (v2-v0) directions. Same as barycentric interpolation.
template <typename TT, typename T>
inline TT interpolate_triangle(
    const TT& v0, const TT& v1, const TT& v2, const vec<T, 2>& uv) {
    return v0 * (1 - uv.x - uv.y) + v1 * uv.x + v2 * uv.y;
}
template <typename TT, typename T>
inline TT interpolate_triangle(
    const std::vector<TT>& vals, const vec<int, 3>& t, const vec<T, 2>& uv) {
    if (vals.empty()) return TT();
    return interpolate_triangle(vals[t.x], vals[t.y], vals[t.z], uv);
}
// Interpolates values over a quad parametrized by u and v along the
// (v1-v0) and (v2-v1) directions. Same as bilear interpolation.
template <typename TT, typename T>
inline TT interpolate_quad(const TT& v0, const TT& v1, const TT& v2,
    const TT& v3, const vec<T, 2>& uv) {
    return v0 * (1 - uv.x) * (1 - uv.y) + v1 * uv.x * (1 - uv.y) +
           v2 * uv.x * uv.y + v3 * (1 - uv.x) * uv.y;
}
template <typename TT, typename T>
inline TT interpolate_quad(
    const std::vector<TT>& vals, const vec<int, 4>& q, const vec<T, 2>& uv) {
    if (vals.empty()) return T();
    if (q.z == q.w)
        return interpolate_triangle(vals[q.x], vals[q.y], vals[q.z]);
    return interpolate_quad(vals[q.x], vals[q.y], vals[q.z], vals[q.w], uv);
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
inline T interpolate_bezier(const std::vector<T>& vals, vec4i b, float u) {
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
    const std::vector<T>& vals, vec4i b, float u) {
    if (vals.empty()) return T();
    return interpolate_bezier_derivative(
        vals[b.x], vals[b.y], vals[b.z], vals[b.w], u);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// SHAPE UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Compute per-vertex normals/tangents for lines/triangles/quads.
void compute_tangents(const std::vector<vec2i>& lines,
    const std::vector<vec3f>& pos, std::vector<vec3f>& tang);
void compute_normals(const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, std::vector<vec3f>& norm);
void compute_normals(const std::vector<vec4i>& quads,
    const std::vector<vec3f>& pos, std::vector<vec3f>& norm);

// Compute per-vertex tangent space for triangle meshes.
// Tangent space is defined by a four component vector.
// The first three components are the tangent with respect to the u texcoord.
// The fourth component is the sign of the tangent wrt the v texcoord.
// Tangent frame is useful in normal mapping.
void compute_tangent_space(const std::vector<vec3i>& triangles,
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
std::vector<vec2i> get_edges(const std::vector<vec3i>& triangles);
std::vector<vec2i> get_edges(const std::vector<vec4i>& quads);

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
void convert_face_varying(std::vector<vec4i>& qquads, std::vector<vec3f>& qpos,
    std::vector<vec3f>& qnorm, std::vector<vec2f>& qtexcoord,
    std::vector<vec4f>& qcolor, const std::vector<vec4i>& quads_pos,
    const std::vector<vec4i>& quads_norm,
    const std::vector<vec4i>& quads_texcoord,
    const std::vector<vec4i>& quads_color, const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec2f>& texcoord,
    const std::vector<vec4f>& color);

// Subdivide lines by splitting each line in half.
template <typename T>
void subdivide_lines(
    std::vector<vec2i>& lines, std::vector<T>& vert, int level);
// Subdivide triangle by splitting each triangle in four, creating new
// vertices for each edge.
template <typename T>
void subdivide_triangles(
    std::vector<vec3i>& triangles, std::vector<T>& vert, int level);
// Subdivide quads by splitting each quads in four, creating new
// vertices for each edge and for each face.
template <typename T>
void subdivide_quads(
    std::vector<vec4i>& quads, std::vector<T>& vert, int level);
// Subdivide beziers by splitting each segment in two.
template <typename T>
void subdivide_beziers(
    std::vector<vec4i>& beziers, std::vector<T>& vert, int level);
// Subdivide quads using Carmull-Clark subdivision rules.
template <typename T>
void subdivide_catmullclark(std::vector<vec4i>& quads, std::vector<T>& vert,
    int level, bool lock_boundary = false);

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

// Weld vertices within a threshold. For noe the implementation is O(n^2).
std::vector<int> weld_vertices(std::vector<vec3f>& pos, float threshold);
void weld_triangles(
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos, float threshold);
void weld_quads(
    std::vector<vec4i>& quads, std::vector<vec3f>& pos, float threshold);

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
    const std::vector<vec2f>& texcoord, int npoints, int seed = 7);

}  // namespace ygl

// -----------------------------------------------------------------------------
// RAY INTERSECTION AND CLOSEST POINT FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Intersect a ray with a point (approximate).
// Based on http://geomalgorithms.com/a02-lines.html.
bool intersect_point(
    const ray3f& ray, vec3f p, float r, float& dist, vec2f& uv);

// Intersect a ray with a line (approximate).
// Based on http://geomalgorithms.com/a05-intersect-1.html and
// http://geomalgorithms.com/a07-distance.html#
//     dist3D_Segment_to_Segment
bool intersect_line(const ray3f& ray, const vec3f& v0, const vec3f& v1,
    float r0, float r1, float& dist, vec2f& uv);

// Intersect a ray with a triangle.
bool intersect_triangle(const ray3f& ray, const vec3f& v0, const vec3f& v1,
    const vec3f& v2, float& dist, vec2f& uv);

// Intersect a ray with a quad represented as two triangles (0,1,3) and
// (2,3,1), with the uv coordinates of the second triangle corrected by u =
// 1-u' and v = 1-v' to produce a quad parametrization where u and v go from 0
// to 1. This is equivalent to Intel's Embree.
bool intersect_quad(const ray3f& ray, const vec3f& v0, const vec3f& v1,
    const vec3f& v2, const vec3f& v3, float& dist, vec2f& uv);

// Intersect a ray with a axis-aligned bounding box.
bool intersect_bbox(const ray3f& ray, const bbox3f& bbox);

// Intersect a ray with a axis-aligned bounding box, implemented as
// "Robust BVH Ray Traversal" by T. Ize published at
// http://jcgt.org/published/0002/02/02/paper.pdf
bool intersect_bbox(
    const ray3f& ray, vec3f ray_dinv, vec3i ray_dsign, const bbox3f& bbox);

// Check if a point overlaps a position within a max distance.
bool overlap_point(const vec3f& pos, float dist_max, const vec3f& v0, float r0,
    float& dist, vec2f& uv);

// Find closest line point to a position.
float closestuv_line(const vec3f& pos, const vec3f& v0, const vec3f& v1);

// Check if a line overlaps a position within a max distance.
bool overlap_line(const vec3f& pos, float dist_max, const vec3f& v0,
    const vec3f& v1, float r0, float r1, float& dist, vec2f& uv);

// Find closest triangle point to a position.
vec2f closestuv_triangle(
    const vec3f& pos, const vec3f& v0, const vec3f& v1, const vec3f& v2);

// Check if a triangle overlaps a position within a max distance.
bool overlap_triangle(const vec3f& pos, float dist_max, const vec3f& v0,
    const vec3f& v1, const vec3f& v2, float r0, float r1, float r2, float& dist,
    vec2f& uv);

// Check if a quad overlaps a position within a max distance.
bool overlap_quad(const vec3f& pos, float dist_max, const vec3f& v0,
    const vec3f& v1, const vec3f& v2, const vec3f& v3, float r0, float r1,
    float r2, float r3, float& dist, vec2f& uv);

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
    internal,
    point,
    line,
    triangle,
    vertex,
    instance
};

// Maximum number of primitives per BVH node.
const int bvh_max_prims = 4;

// BVH tree node containing its bounds, indices to the BVH arrays of either
// primitives or internal nodes, the node element type,
// and the split axis. Leaf and internal nodes are identical, except that
// indices refer to primitives for leaf nodes or other nodes for internal
// nodes. See bvh_tree for more details.
struct bvh_node {
    bbox3f bbox;                    // bouds
    uint32_t prims[bvh_max_prims];  // primitives
    uint16_t count;                 // number of prims
    bvh_node_type type;             // node type
    uint8_t split_axis;             // split axis
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
    std::vector<vec3f> pos;        // Positions for shape BVHs.
    std::vector<float> radius;     // Radius for shape BVHs.
    std::vector<int> points;       // Points for shape BVHs.
    std::vector<vec2i> lines;      // Lines for shape BVHs.
    std::vector<vec3i> triangles;  // Triangles for shape BVHs.

    // data for instance BVH
    std::vector<frame3f> ist_frames;                  // instance frames
    std::vector<frame3f> ist_inv_frames;              // instance inverse frames
    std::vector<std::shared_ptr<bvh_tree>> ist_bvhs;  // instance shape bvhs

    // bvh nodes
    std::vector<bvh_node> nodes;  // Internal nodes.
};

// Build a BVH from the given set of primitives.
void build_bvh(const std::shared_ptr<bvh_tree>& bvh, bool equalsize);
// Update the node bounds for a shape bvh.
void refit_bvh(const std::shared_ptr<bvh_tree>& bvh);

// Intersect ray with a bvh returning either the first or any intersection
// depending on `find_any`. Returns the ray distance `dist`, the instance
// id `iid`, the shape id `sid`, the shape element index `eid` and the
// shape barycentric coordinates `uv`.
bool intersect_bvh(const std::shared_ptr<bvh_tree>& bvh, const ray3f& ray,
    bool find_any, float& dist, int& iid, int& eid, vec2f& uv);

// Find a shape element that overlaps a point within a given distance
// `max_dist`, returning either the closest or any overlap depending on
// `find_any`. Returns the point distance `dist`, the instance id `iid`, the
// shape id `sid`, the shape element index `eid` and the shape barycentric
// coordinates `uv`.
bool overlap_bvh(const std::shared_ptr<bvh_tree>& bvh, const vec3f& pos,
    float max_dist, bool find_any, float& dist, int& iid, int& eid, vec2f& uv);

}  // namespace ygl

// -----------------------------------------------------------------------------
// SHAPE EXAMPLES
// -----------------------------------------------------------------------------
namespace ygl {

// Make examples triangle shapes with shared vertices (not watertight).
void make_quad(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec2i& steps,
    const vec2f& size, const vec2f& uvsize);
void make_quad_stack(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec3i& steps,
    const vec3f& size, const vec2f& uvsize);
void make_cube(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec3i& steps,
    const vec3f& size, const vec3f& uvsize);
void make_cube_rounded(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec3i& steps,
    const vec3f& size, const vec3f& uvsize, float radius);
void make_sphere(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec2i& steps,
    float size, const vec2f& uvsize);
void make_sphere_cube(std::vector<vec4i>& quadss, std::vector<vec3f>& pos,
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
void make_disk_bulged(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, int steps,
    float size, float uvsize, float height);
void make_cylinder_side(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec2i& steps,
    const vec2f& size, const vec2f& uvsize);
void make_cylinder(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec3i& steps,
    const vec2f& size, const vec3f& uvsize);
void make_cylinder_rounded(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec3i& steps,
    const vec2f& size, const vec3f& uvsize, float radius);
void make_geodesic_sphere(std::vector<vec3i>& triangles,
    std::vector<vec3f>& pos, int tesselation, float size);

// Make example watertight quad meshes for subdivision surfaces.
void make_fvcube(std::vector<vec4i>& quads_pos, std::vector<vec3f>& pos,
    std::vector<vec4i>& quads_norm, std::vector<vec3f>& norm,
    std::vector<vec4i>& quads_texcoord, std::vector<vec2f>& texcoord,
    vec3i steps, const vec3f& size, const vec3f& uvsize);
void make_suzanne(std::vector<vec4i>& quads, std::vector<vec3f>& pos);
void make_cube(std::vector<vec4i>& quads, std::vector<vec3f>& pos);

// Generate lines set along a quad.
void make_lines(std::vector<vec2i>& lines, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<float>& radius, const vec2i& steps, const vec2f& size,
    const vec2f& uvsize, vec2f line_radius = {0.001f, 0.001f});

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

// Make a hair ball around a shape.
// length: minimum and maximum length
// rad: minimum and maximum radius from base to tip
// noise: noise added to hair (strength/scale)
// clump: clump added to hair (number/strength)
// rotation: rotation added to hair (angle/strength)
void make_hair(std::vector<vec2i>& lines, std::vector<vec3f>& pos,
    std::vector<vec3f>& tang, std::vector<vec2f>& texcoord,
    std::vector<float>& radius, const vec2i& steps,
    const std::vector<vec3i>& striangles, const std::vector<vec3f>& spos,
    const std::vector<vec3f>& snorm, const std::vector<vec2f>& stexcoord,
    const vec2f& length = {0.1f, 0.1f}, const vec2f& rad = {0.001f, 0.001f},
    const vec2f& noise = zero2f, const vec2f& clump = zero2f,
    const vec2f& rotation = zero2f, int seed = 7);

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMAGE TYPE
// -----------------------------------------------------------------------------
namespace ygl {

// Image container.
template <typename T>
struct image {
    // constructors
    image() : w(0), h(0), pxl() {}
    image(int width, int height, const T& c = T())
        : w(width), h(height), pxl(width * height, c) {}
    image(int width, int height, const std::vector<T>& p)
        : w(width), h(height), pxl(p) {}

    // size checks
    int width() const { return w; };
    int height() const { return h; };
    int size() const { return w * h; };
    bool empty() const { return !w || !h; };

    // pixel access
    T& operator[](int idx) { return pxl[idx]; }
    const T& operator[](int idx) const { return pxl[idx]; }
    T& operator[](vec2i ij) { return pxl[ij.y * w + ij.x]; }
    const T& operator[](vec2i ij) const { return pxl[ij.y * w + ij.x]; }
    T& at(int idx) { return pxl.at(idx); }
    const T& at(int idx) const { return pxl.at(idx); }
    T& at(vec2i ij) { return pxl.at(ij.y * w + ij.x); }
    const T& at(vec2i ij) const { return pxl.at(ij.y * w + ij.x); }

    // data access
    T* data() { return pxl.data(); }
    const T* data() const { return pxl.data(); }
    std::vector<T>& pixels() { return pxl; }
    const std::vector<T>& pixels() const { return pxl; }
    T* begin() { return pxl.data(); }
    T* end() { return pxl.data() + pxl.size(); }
    const T* begin() const { return pxl.data(); }
    const T* end() const { return pxl.data() + pxl.size(); }

   private:
    int w = 0;
    int h = 0;
    std::vector<T> pxl;
};

// Type aliases
using image4f = image<vec4f>;
using image4b = image<vec4b>;

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMAGE UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Conversion from/to floats.
image4f byte_to_float(const image4b& bt);
image4b float_to_byte(const image4f& fl);

// Conversion between linear and gamma-encoded images.
image4f gamma_to_linear(const image4f& srgb, float gamma = 2.2f);
image4f linear_to_gamma(const image4f& lin, float gamma = 2.2f);

// Apply exposure and filmic tone mapping
image4f expose_image(const image4f& hdr, float exposure);
image4f filmic_tonemap_image(const image4f& hdr);

// Resize an image.
image4f resize_image(const image4f& img, int res_width, int res_height);

}  // namespace ygl

// -----------------------------------------------------------------------------
// COLOR CONVERSION UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Element-wise float to byte conversion.
inline vec4b float_to_byte(const vec4f& a) {
    return {(byte)clamp(int(a.x * 256), 0, 255),
        (byte)clamp(int(a.y * 256), 0, 255),
        (byte)clamp(int(a.z * 256), 0, 255),
        (byte)clamp(int(a.w * 256), 0, 255)};
}
inline vec4f byte_to_float(const vec4b& a) {
    return {a.x / 255.0f, a.y / 255.0f, a.z / 255.0f, a.w / 255.0f};
}

// Conversion between linear and gamma-encoded images.
inline vec3f gamma_to_linear(const vec3f& srgb, float gamma = 2.2f) {
    return {pow(srgb.x, gamma), pow(srgb.y, gamma), pow(srgb.z, gamma)};
}
inline vec3f linear_to_gamma(const vec3f& lin, float gamma = 2.2f) {
    return {
        pow(lin.x, 1 / gamma), pow(lin.y, 1 / gamma), pow(lin.z, 1 / gamma)};
}
inline vec4f gamma_to_linear(const vec4f& srgb, float gamma = 2.2f) {
    return {pow(srgb.x, gamma), pow(srgb.y, gamma), pow(srgb.z, gamma), srgb.w};
}
inline vec4f linear_to_gamma(const vec4f& lin, float gamma = 2.2f) {
    return {pow(lin.x, 1 / gamma), pow(lin.y, 1 / gamma), pow(lin.z, 1 / gamma),
        lin.w};
}

// Approximate luminance estimate
inline float luminance(const vec3f& a) { return (a.x + a.y + a.z) / 3; }

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
// EXAMPLE IMAGES
// -----------------------------------------------------------------------------
namespace ygl {

// Make example images.
image4f make_grid_image(int width, int height, int tile = 8,
    const vec4f& c0 = {0.5f, 0.5f, 0.5f, 1},
    const vec4f& c1 = {0.8f, 0.8f, 0.8f, 1});
image4f make_checker_image(int width, int height, int tile = 8,
    const vec4f& c0 = {0.5f, 0.5f, 0.5f, 1},
    const vec4f& c1 = {0.8f, 0.8f, 0.8f, 1});
image4f make_bumpdimple_image(int width, int height, int tile = 8);
image4f make_ramp_image(int width, int height, const vec4f& c0, const vec4f& c1,
    float srgb = false);
image4f make_gammaramp_image(int width, int height);
image4f make_uv_image(int width, int height);
image4f make_uvgrid_image(
    int width, int height, int tile = 8, bool colored = true);

// Comvert a bump map to a normal map.
image4f bump_to_normal_map(const image4f& img, float scale = 1);

// Make a sunsky HDR model with sun at theta elevation in [0,pi/2], turbidity
// in [1.7,10] with or without sun.
image4f make_sunsky_image(int width, int height, float thetaSun,
    float turbidity = 3, bool has_sun = false,
    vec3f ground_albedo = {0.7f, 0.7f, 0.7f});
// Make an image of multiple lights.
image4f make_lights_image(int width, int height, const vec3f& le = {1, 1, 1},
    int nlights = 4, float langle = pi / 4, float lwidth = pi / 16,
    float lheight = pi / 16);

// Make a noise image. Wrap works only if both resx and resy are powers of two.
image4f make_noise_image(
    int width, int height, float scale = 1, bool wrap = true);
image4f make_fbm_image(int width, int height, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, int octaves = 6, bool wrap = true);
image4f make_ridge_image(int width, int height, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, float offset = 1.0f,
    int octaves = 6, bool wrap = true);
image4f make_turbulence_image(int width, int height, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, int octaves = 6, bool wrap = true);

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
template <typename T, typename T1>
inline T eval_keyframed_step(
    const std::vector<T1>& times, const std::vector<T>& vals, T1 time) {
    if (time <= times.front()) return vals.front();
    if (time >= times.back()) return vals.back();
    time = clamp(time, times.front(), times.back() - (T1)0.001);
    auto idx = eval_keyframed_index(times, time);
    return vals.at(idx - 1);
}

// Evalautes a keyframed value using linear interpolation.
template <typename T, typename T1>
inline vec<T, 4> eval_keyframed_slerp(
    const std::vector<T1>& times, const std::vector<vec<T, 4>>& vals, T1 time) {
    if (time <= times.front()) return vals.front();
    if (time >= times.back()) return vals.back();
    time = clamp(time, times.front(), times.back() - (T1)0.001);
    auto idx = eval_keyframed_index(times, time);
    auto t = (time - times.at(idx - 1)) / (times.at(idx) - times.at(idx - 1));
    return slerp(vals.at(idx - 1), vals.at(idx), t);
}

// Evalautes a keyframed value using linear interpolation.
template <typename T, typename T1>
inline T eval_keyframed_linear(
    const std::vector<T1>& times, const std::vector<T>& vals, T1 time) {
    if (time <= times.front()) return vals.front();
    if (time >= times.back()) return vals.back();
    time = clamp(time, times.front(), times.back() - (T1)0.001);
    auto idx = eval_keyframed_index(times, time);
    auto t = (time - times.at(idx - 1)) / (times.at(idx) - times.at(idx - 1));
    return vals.at(idx - 1) * (1 - t) + vals.at(idx) * t;
}

// Evalautes a keyframed value using Bezier interpolation.
template <typename T, typename T1>
inline T eval_keyframed_bezier(
    const std::vector<T1>& times, const std::vector<T>& vals, T1 time) {
    if (time <= times.front()) return vals.front();
    if (time >= times.back()) return vals.back();
    time = clamp(time, times.front(), times.back() - (T1)0.001);
    auto idx = eval_keyframed_index(times, time);
    auto t = (time - times.at(idx - 1)) / (times.at(idx) - times.at(idx - 1));
    return interpolate_bezier(
        vals.at(idx - 3), vals.at(idx - 2), vals.at(idx - 1), vals.at(idx), t);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// SCENE DATA
// -----------------------------------------------------------------------------
namespace ygl {

// forward declaration
struct bvh_tree;

// Camera.
struct camera {
    std::string name = "";             // name
    frame3f frame = identity_frame3f;  // transform frame
    bool ortho = false;                // orthographic
    float width = 0.036f;              // film width (default to 35mm)
    float height = 0.024f;             // film height (default to 35mm)
    float focal = 0.050f;              // focal length (defaul to 50 mm)
    float focus = flt_max;             // focal distance (default to inf focus)
    float aperture = 0;                // lens aperture
    float near = 0.01f;                // near plane distance
    float far = 10000;                 // far plane distance
};

// Texture containing either an LDR or HDR image.
struct texture {
    std::string name = "";  // name
    std::string path = "";  // file path
    image4f img = {};       // image
    bool clamp = false;     // clamp textures coordinates
    float scale = 1;        // scale for occ, normal, bumps
    uint gl_txt = 0;        // unmanaged data for OpenGL viewer
};

// Material for surfaces, lines and triangles.
// For surfaces, uses a microfacet model with thin sheet transmission.
// The model is based on glTF for compatibility and adapted to OBJ.
// For lines, uses Kajija-Kay model. For points, a hacked up shading.
struct material {
    std::string name = "";       // name
    bool base_metallic = false;  // base-metallic parametrization
    bool gltf_textures = false;  // glTF packed textures
    bool double_sided = false;   // double sided rendering

    // base values
    vec3f ke = {0, 0, 0};  // emission color
    vec3f kd = {0, 0, 0};  // diffuse/base color
    vec3f ks = {0, 0, 0};  // specular color / metallic factor
    vec3f kt = {0, 0, 0};  // transmission color
    float rs = 0.0001;     // roughness mapped as glTF
    float op = 1;          // opacity
    bool fresnel = true;   // whether to use fresnel in reflections/transmission
    bool refract = false;  // whether to use use refraction in tranmission

    // textures
    std::shared_ptr<texture> ke_txt;    // emission texture
    std::shared_ptr<texture> kd_txt;    // diffuse texture
    std::shared_ptr<texture> ks_txt;    // specular texture
    std::shared_ptr<texture> kt_txt;    // transmission texture
    std::shared_ptr<texture> rs_txt;    // roughness texture
    std::shared_ptr<texture> op_txt;    // opacity texture
    std::shared_ptr<texture> occ_txt;   // occlusion texture
    std::shared_ptr<texture> bump_txt;  // bump map texture (heighfield)
    std::shared_ptr<texture> disp_txt;  // displacement map texture (heighfield)
    std::shared_ptr<texture> norm_txt;  // normal texture
};

// Shape data represented as an indexed meshes of elements.
// May contain either tringles, lines or a set of vertices.
struct shape {
    std::string name = "";  // name
    std::string path = "";  // path for glTF buffers

    // primitives
    std::vector<int> points;       // points
    std::vector<vec2i> lines;      // lines
    std::vector<vec3i> triangles;  // triangles

    // vertex data
    std::vector<vec3f> pos;       // positions
    std::vector<vec3f> norm;      // normals/tangents
    std::vector<vec2f> texcoord;  // texcoord coordinates
    std::vector<vec4f> color;     // colors
    std::vector<float> radius;    // radia for lines/points
    std::vector<vec4f> tangsp;    // tangent space for triangles

    // computed properties
    bbox3f bbox = invalid_bbox3f;             // boudning box
    std::vector<float> elem_cdf = {};         // element cdf for sampling
    std::shared_ptr<bvh_tree> bvh = nullptr;  // bvh for ray intersection
    uint gl_pos = 0, gl_norm = 0, gl_texcoord = 0, gl_color = 0, gl_tangsp = 0,
         gl_points = 0, gl_lines = 0,
         gl_triangles = 0;  // unmanaged data for OpenGL viewer
};

// Subdivision surface.
struct subdiv {
    std::string name = "";        // name
    std::string path = "";        // path for glTF buffers
    int level = 0;                // subdivision level
    bool catmull_clark = true;    // catmull clark subdiv
    bool compute_normals = true;  // faceted subdivision

    // primitives
    std::vector<vec4i> quads_pos;       // quads for position
    std::vector<vec4i> quads_texcoord;  // quads for texture coordinates
    std::vector<vec4i> quads_color;     // quads for color

    // creases
    std::vector<vec3i> crease_pos;       // crease for position
    std::vector<vec3i> crease_texcoord;  // crease for texture coordinates

    // vertex data
    std::vector<vec3f> pos;       // positions
    std::vector<vec2f> texcoord;  // texcoord coordinates
    std::vector<vec4f> color;     // colors
};

// Shape instance.
struct instance {
    std::string name;                         // name
    frame3f frame = identity_frame3f;         // transform frame
    std::shared_ptr<shape> shp = nullptr;     // shape
    std::shared_ptr<material> mat = nullptr;  // material
    std::shared_ptr<subdiv> sbd = nullptr;    // subdivision shape

    // compute properties
    bbox3f bbox = invalid_bbox3f;  // boudning box
};

// Envinonment map.
struct environment {
    std::string name = "";                 // name
    frame3f frame = identity_frame3f;      // transform frame
    vec3f ke = {0, 0, 0};                  // emission color
    std::shared_ptr<texture> ke_txt = {};  // emission texture

    // computed properties
    std::vector<float> elem_cdf;  // element cdf for sampling
};

// Node in a transform hierarchy.
struct node {
    std::string name = "";                       // name
    std::shared_ptr<node> parent = nullptr;      // parent
    frame3f frame = identity_frame3f;            // transform frame
    vec3f translation = {0, 0, 0};               // translation
    vec4f rotation = {0, 0, 0, 1};               // rotation
    vec3f scale = {1, 1, 1};                     // scale
    std::vector<float> weights = {};             // morph weights
    std::shared_ptr<camera> cam = nullptr;       // camera
    std::shared_ptr<instance> ist = nullptr;     // instance
    std::shared_ptr<environment> env = nullptr;  // environment

    // compute properties
    std::vector<std::weak_ptr<node>> children = {};  // child nodes
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
    std::vector<std::shared_ptr<node>> targets;    // target nodes
};

// Scene comprised an array of objects whose memory is owened by the scene.
// All members are optional,Scene objects (camera, instances, environments)
// have transforms defined internally. A scene can optionally contain a
// node hierarchy where each node might point to a camera, instance or
// environment. In that case, the element transforms are computed from
// the hierarchy. Animation is also optional, with keyframe data that
// updates node transformations only if defined.
struct scene {
    std::string name;                                       // name
    std::vector<std::shared_ptr<camera>> cameras = {};      // cameras
    std::vector<std::shared_ptr<shape>> shapes = {};        // shapes
    std::vector<std::shared_ptr<subdiv>> subdivs = {};      // subdivs
    std::vector<std::shared_ptr<instance>> instances = {};  // instances
    std::vector<std::shared_ptr<material>> materials = {};  // materials
    std::vector<std::shared_ptr<texture>> textures = {};    // textures
    std::vector<std::shared_ptr<environment>> environments =
        {};  // environments

    std::vector<std::shared_ptr<node>> nodes = {};  // node hierarchy [optional]
    std::vector<std::shared_ptr<animation>> animations =
        {};  // animations [optional]

    // compute properties
    std::vector<std::shared_ptr<instance>> lights;
    bbox3f bbox = invalid_bbox3f;  // boudning box
    std::shared_ptr<bvh_tree> bvh = nullptr;
};

}  // namespace ygl

// -----------------------------------------------------------------------------
// SCENE UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Print scene statistics.
void print_stats(const std::shared_ptr<scene>& scn);

// Merge scene into one another. Note that the objects are _moved_ from
// merge_from to merged_into, so merge_from will be empty after this function.
void merge_into(const std::shared_ptr<scene>& merge_into,
    std::shared_ptr<scene> merge_from);

}  // namespace ygl

// -----------------------------------------------------------------------------
// UPDATES TO COMPUTED PROPERTIES
// -----------------------------------------------------------------------------
namespace ygl {

// Update the normals of a shape.
void update_normals(const std::shared_ptr<shape>& shp);

// Update node transforms.
void update_transforms(const std::shared_ptr<scene>& scn, float time = 0,
    const std::string& anim_group = "");
// Compute animation range.
vec2f compute_animation_range(
    const std::shared_ptr<scene>& scn, const std::string& anim_group = "");

// Computes shape/scene approximate bounds.
void udpate_bbox(const std::shared_ptr<shape>& shp);
void update_bbox(const std::shared_ptr<scene>& scn, bool do_shapes = true);

// Update lights.
void update_lights(const std::shared_ptr<scene>& scn, bool do_shapes = true,
    bool do_environments = false);
// Generate a distribution for sampling a shape uniformly based on area/length.
void update_shape_cdf(const std::shared_ptr<shape>& shp);
// Generate a distribution for sampling an environment texture uniformly
// based on angle and texture intensity.
void update_environment_cdf(std::shared_ptr<environment> env);

// Updates/refits bvh.
void update_bvh(const std::shared_ptr<shape>& shp, bool equalsize = false);
void update_bvh(const std::shared_ptr<scene>& scn, bool do_shapes = true,
    bool equalsize = false);
void refit_bvh(const std::shared_ptr<shape>& shp);
void refit_bvh(const std::shared_ptr<scene>& scn, bool do_shapes = true);

// Updates tesselation.
void update_tesselation(
    const std::shared_ptr<subdiv>& sbd, std::shared_ptr<shape> shp);
void update_tesselation(const std::shared_ptr<scene>& scn);

// Add missing names, normals, tangents and hierarchy.
void add_missing_names(const std::shared_ptr<scene>& scn);
void add_missing_normals(const std::shared_ptr<scene>& scn);
void add_missing_tangent_space(const std::shared_ptr<scene>& scn);
void add_missing_materials(const std::shared_ptr<scene>& scn);
// Checks for validity of the scene.
std::vector<std::string> validate(
    const std::shared_ptr<scene>& scn, bool skip_textures = false);

// make camera
std::shared_ptr<camera> make_bbox_camera(const std::string& name,
    const bbox3f& bbox, float width = 0.036f, float height = 0.024f,
    float focal = 0.050f);
// make default material
inline std::shared_ptr<material> make_default_material(
    const std::string& name) {
    auto mat = std::make_shared<material>();
    mat->name = name;
    mat->kd = {0.2f, 0.2f, 0.2f};
    return mat;
}

// Add a sky environment
inline std::shared_ptr<environment> make_sky_environment(
    const std::string& name, float sun_angle = pi / 4) {
    auto txt = std::make_shared<texture>();
    txt->name = name;
    txt->path = "textures/" + name + ".hdr";
    txt->img = make_sunsky_image(1024, 512, sun_angle);
    auto env = std::make_shared<environment>();
    env->name = name;
    env->ke = {1, 1, 1};
    env->ke_txt = txt;
    return env;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// INTERSECTION, EVAL AND SAMPLING FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Scene intersection.
struct scene_intersection {
    std::shared_ptr<instance> ist =
        nullptr;        // instance or null for no intersection
    int ei = 0;         // shape element index
    vec2f uv = zero2f;  // shape element coordinates
    float dist = 0;     // ray/point distance
};

// Intersects a ray with the scene.
scene_intersection intersect_ray(
    const std::shared_ptr<scene>& scn, const ray3f& ray, bool find_any = false);

// Shape values interpolated using barycentric coordinates.
vec3f eval_pos(const std::shared_ptr<shape>& shp, int ei, const vec2f& uv);
vec3f eval_norm(const std::shared_ptr<shape>& shp, int ei, const vec2f& uv);
vec2f eval_texcoord(const std::shared_ptr<shape>& shp, int ei, const vec2f& uv);
vec4f eval_color(const std::shared_ptr<shape>& shp, int ei, const vec2f& uv);
float eval_radius(const std::shared_ptr<shape>& shp, int ei, const vec2f& uv);
vec4f eval_tangsp(const std::shared_ptr<shape>& shp, int ei, const vec2f& uv);
vec3f eval_tangsp(const std::shared_ptr<shape>& shp, int ei, const vec2f& uv,
    bool& left_handed);
// Shape element values.
vec3f eval_elem_norm(const std::shared_ptr<shape>& shp, int ei);
vec4f eval_elem_tangsp(const std::shared_ptr<shape>& shp, int ei);

// Instance values interpolated using barycentric coordinates.
// Handles defaults if data is missing.
vec3f eval_pos(const std::shared_ptr<instance>& ist, int ei, const vec2f& uv);
vec3f eval_norm(const std::shared_ptr<instance>& ist, int ei, const vec2f& uv);
vec2f eval_texcoord(
    const std::shared_ptr<instance>& ist, int ei, const vec2f& uv);
vec4f eval_color(const std::shared_ptr<instance>& ist, int ei, const vec2f& uv);
float eval_radius(
    const std::shared_ptr<instance>& ist, int ei, const vec2f& uv);
vec3f eval_tangsp(const std::shared_ptr<instance>& ist, int ei, const vec2f& uv,
    bool& left_handed);
// Instance element values.
vec3f eval_elem_norm(const std::shared_ptr<instance>& ist, int ei);
// Shading normals including material perturbations.
vec3f eval_shading_norm(
    const std::shared_ptr<instance>& ist, int ei, const vec2f& uv, vec3f o);

// Environment texture coordinates from the incoming direction.
vec2f eval_texcoord(const std::shared_ptr<environment>& env, vec3f i);
// Evaluate the incoming direction from the uv.
vec3f eval_direction(const std::shared_ptr<environment>& env, const vec2f& uv);
// Evaluate the environment emission.
vec3f eval_environment(const std::shared_ptr<environment>& env, vec3f i);

// Evaluate a texture.
vec4f eval_texture(const std::shared_ptr<texture>& txt, const vec2f& texcoord);

// Set and evaluate camera parameters. Setters take zeros as default values.
float eval_camera_fovy(const std::shared_ptr<camera>& cam);
float eval_camera_aspect(const std::shared_ptr<camera>& cam);
void set_camera_fovy(const std::shared_ptr<camera>& cam, float fovy,
    float aspect, float width = 0.036f);

// Generates a ray from a camera image coordinate `uv` and lens coordinates
// `luv`.
ray3f eval_camera_ray(
    const std::shared_ptr<camera>& cam, const vec2f& uv, const vec2f& luv);
// Generates a ray from a camera for pixel coordinates `ij`, the resolution
// `res`, the sub-pixel coordinates `puv` and the lens coordinates `luv` and
// the image resolution `res`.
ray3f eval_camera_ray(const std::shared_ptr<camera>& cam, const vec2i& ij,
    const vec2i& imsize, const vec2f& puv, const vec2f& luv);

// Evaluates material parameters: emission, diffuse, specular, transmission,
// roughness and opacity.
vec3f eval_emission(
    const std::shared_ptr<instance>& ist, int ei, const vec2f& uv);
vec3f eval_diffuse(
    const std::shared_ptr<instance>& ist, int ei, const vec2f& uv);
vec3f eval_specular(
    const std::shared_ptr<instance>& ist, int ei, const vec2f& uv);
vec3f eval_transmission(
    const std::shared_ptr<instance>& ist, int ei, const vec2f& uv);
float eval_roughness(
    const std::shared_ptr<instance>& ist, int ei, const vec2f& uv);
float eval_opacity(
    const std::shared_ptr<instance>& ist, int ei, const vec2f& uv);

// Material values packed into a convenience structure.
struct bsdf {
    vec3f kd = zero3f;     // diffuse
    vec3f ks = zero3f;     // specular
    vec3f kt = zero3f;     // transmission
    float rs = 1;          // roughness
    bool refract = false;  // whether to use refraction in transmission
};
bsdf eval_bsdf(const std::shared_ptr<instance>& ist, int ei, const vec2f& uv);
bool is_delta_bsdf(const bsdf& f);

// Sample a shape based on a distribution.
std::pair<int, vec2f> sample_shape(
    const std::shared_ptr<shape>& shp, float re, const vec2f& ruv);

// Sample an environment uniformly.
vec2f sample_environment(
    const std::shared_ptr<environment>& env, const vec2f& ruv);

}  // namespace ygl

// -----------------------------------------------------------------------------
// PATH TRACING
// -----------------------------------------------------------------------------
namespace ygl {

// Trace evaluation function.
using trace_func = std::function<vec3f(const std::shared_ptr<scene>& scn,
    const ray3f& ray, rng_state& rng, int nbounces, bool* hit)>;

// Init a sequence of random number generators.
inline image<rng_state> make_trace_rngs(int w, int h, uint64_t seed) {
    auto rngs = image<rng_state>{w, h};
    int rseed = 1301081;  // large prime
    for (auto i = 0; i < w * h; i++) {
        rngs[i] = make_rng(seed, rseed + 1);
        rseed = (rseed * 1103515245 + 12345) & ((1U << 31) - 1);  // bsd rand
    }
    return rngs;
}

// Trace the next nsamples samples. Assumes that the
// image contains cur_samples already. Returns true when done.
void trace_samples(const std::shared_ptr<scene>& scn,
    const std::shared_ptr<camera>& cam, image4f& img, image<rng_state>& rngs,
    int cur_samples, int nsamples, trace_func tracer, int nbounces,
    float pixel_clamp = 100);
// Like before but with multiplthreading.
void trace_samples_mt(const std::shared_ptr<scene>& scn,
    const std::shared_ptr<camera>& cam, image4f& img, image<rng_state>& rngs,
    int cur_samples, int nsamples, trace_func tracer, int nbounces,
    float pixel_clamp = 100);

// Starts an anyncrhounous renderer.
void trace_async_start(const std::shared_ptr<scene>& scn,
    const std::shared_ptr<camera>& cam, image4f& img, image<rng_state>& rngs,
    int nsamples, trace_func tracer, int nbounces,
    std::vector<std::thread>& threads, bool& stop_flag, int& cur_sample,
    float pixel_clamp = 100);
// Stop the asynchronous renderer.
void trace_async_stop(std::vector<std::thread>& threads, bool& stop_flag);

// Trace function - path tracer.
vec3f trace_path(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit = nullptr);
// Trace function - path tracer without mis.
vec3f trace_path_nomis(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit = nullptr);
// Trace function - naive path tracer.
vec3f trace_path_naive(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit = nullptr);
// Trace function - direct illumination.
vec3f trace_direct(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit = nullptr);
// Trace function - direct illumination without mis.
vec3f trace_direct_nomis(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit = nullptr);
// Trace function - pure environment illumination with no shadows.
vec3f trace_environment(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit = nullptr);
// Trace function - eyelight rendering.
vec3f trace_eyelight(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit = nullptr);
// Trace function - normal debug visualization.
vec3f trace_debug_normal(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit = nullptr);
// Trace function - faceforward debug visualization.
vec3f trace_debug_frontfacing(const std::shared_ptr<scene>& scn,
    const ray3f& ray, rng_state& rng, int nbounces, bool* hit = nullptr);
// Trace function - albedo debug visualization.
vec3f trace_debug_albedo(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit = nullptr);
// Trace function - diffuse debug visualization.
vec3f trace_debug_diffuse(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit = nullptr);
// Trace function - specular debug visualization.
vec3f trace_debug_specular(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit = nullptr);
// Trace function - roughness debug visualization.
vec3f trace_debug_roughness(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit = nullptr);
// Trace function - texcoord debug visualization.
vec3f trace_debug_texcoord(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit = nullptr);

}  // namespace ygl

// -----------------------------------------------------------------------------
// PATH TRACING SUPPORT FUNCTION
// -----------------------------------------------------------------------------
namespace ygl {

// Phong exponent to roughness.
float specular_exponent_to_roughness(float n);

// Specular to fresnel eta.
void specular_fresnel_from_ks(const vec3f& ks, vec3f& es, vec3f& esk);
float specular_to_eta(const vec3f& ks);
// Compute the fresnel term for dielectrics.
vec3f fresnel_dielectric(float cosw, const vec3f& eta_);
// Compute the fresnel term for metals.
vec3f fresnel_metal(float cosw, const vec3f& eta, const vec3f& etak);
// Schlick approximation of Fresnel term, optionally weighted by rs;
vec3f fresnel_schlick(const vec3f& ks, float cosw);
vec3f fresnel_schlick(const vec3f& ks, float cosw, float rs);
vec3f fresnel_schlick(const vec3f& ks, const vec3f& h, const vec3f& o);
vec3f fresnel_schlick(
    const vec3f& ks, const vec3f& h, const vec3f& o, float rs);

// Evaluates the GGX distribution and geometric term.
float eval_ggx(float rs, float ndh, float ndi, float ndo);
// Sample the GGX distribution.
vec3f sample_ggx(float rs, const vec2f& rn);
float sample_ggx_pdf(float rs, float ndh);

// Evaluates the GGX distribution and geometric term.
float eval_ggx_dist(float rs, const vec3f& n, const vec3f& h);
float eval_ggx_sm(float rs, const vec3f& n, const vec3f& o, const vec3f& i);

}  // namespace ygl

// -----------------------------------------------------------------------------
// TIMER UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// get the timer time in nanoseconds
inline int64_t get_time() {
    return std::chrono::high_resolution_clock::now().time_since_epoch().count();
}
// get duration string from nanoseconds
inline std::string format_duration(int64_t duration) {
    auto elapsed = duration / 1000000;  // milliseconds
    auto hours = (int)(elapsed / 3600000);
    elapsed %= 3600000;
    auto mins = (int)(elapsed / 60000);
    elapsed %= 60000;
    auto secs = (int)(elapsed / 1000);
    auto msecs = (int)(elapsed % 1000);
    char buf[256];
    sprintf(buf, "%02d:%02d:%02d.%03d", hours, mins, secs, msecs);
    return buf;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR MATRICES
// -----------------------------------------------------------------------------
namespace ygl {

// Matrix adjugates, determinant and inverses.
template <typename T>
inline mat<T, 2> adjugate(const mat<T, 2>& a) {
    return {{a.y.y, -a.x.y}, {-a.y.x, a.x.x}};
}
template <typename T>
inline mat<T, 3> adjugate(const mat<T, 3>& a) {
    return {{a.y.y * a.z.z - a.z.y * a.y.z, a.z.y * a.x.z - a.x.y * a.z.z,
                a.x.y * a.y.z - a.y.y * a.x.z},
        {a.y.z * a.z.x - a.z.z * a.y.x, a.z.z * a.x.x - a.x.z * a.z.x,
            a.x.z * a.y.x - a.y.z * a.x.x},
        {a.y.x * a.z.y - a.z.x * a.y.y, a.z.x * a.x.y - a.x.x * a.z.y,
            a.x.x * a.y.y - a.y.x * a.x.y}};
}
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
template <typename T>
inline T determinant(const mat<T, 2>& a) {
    return a.x.x * a.y.y - a.x.y * a.y.x;
}
template <typename T>
inline T determinant(const mat<T, 3>& a) {
    return a.x.x * (a.y.y * a.z.z - a.z.y * a.y.z) +
           a.x.y * (a.y.z * a.z.x - a.z.z * a.y.x) +
           a.x.z * (a.y.x * a.z.y - a.z.x * a.y.y);
}
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

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR UI UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Turntable for UI navigation.
template <typename T, typename T1>
inline void camera_turntable(vec<T, 3>& from, vec<T, 3>& to, vec<T, 3>& up,
    const vec<T, 2>& rotate, T1 dolly, const vec<T, 2>& pan) {
    // rotate if necessary
    if (rotate.x || rotate.y) {
        auto z = normalize(to - from);
        auto lz = length(to - from);
        auto phi = atan2(z.z, z.x) + rotate.x;
        auto theta = acos(z.y) + rotate.y;
        theta = clamp(theta, (T)0.001, pi - (T)0.001);
        auto nz = vec<T, 3>{sin(theta) * cos(phi) * lz, cos(theta) * lz,
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
        auto t = vec<T, 3>{pan.x * x.x + pan.y * y.x, pan.x * x.y + pan.y * y.y,
            pan.x * x.z + pan.y * y.z};
        from += t;
        to += t;
    }
}

// Turntable for UI navigation.
template <typename T, typename T1>
inline void camera_turntable(frame<T, 3>& frame, float& focus,
    const vec<T, 2>& rotate, T1 dolly, const vec<T, 2>& pan) {
    // rotate if necessary
    if (rotate != zero2f) {
        auto phi = atan2(frame.z.z, frame.z.x) + rotate.x;
        auto theta = acos(frame.z.y) + rotate.y;
        theta = clamp(theta, (T)0.001, pi - (T)0.001);
        auto new_z =
            vec<T, 3>{sin(theta) * cos(phi), cos(theta), sin(theta) * sin(phi)};
        auto new_center = frame.o - frame.z * focus;
        auto new_o = new_center + new_z * focus;
        frame = lookat_frame(new_o, new_center, {0, 1, 0});
        focus = length(new_o - new_center);
    }

    // pan if necessary
    if (dolly) {
        auto c = frame.o - frame.z * focus;
        focus = max(focus * (1 + dolly), 0.001f);
        frame.o = c + frame.z * focus;
    }

    // pan if necessary
    if (pan.x || pan.y) { frame.o += frame.x * pan.x + frame.y * pan.y; }
}

// FPS camera for UI navigation for a frame parametrization.
template <typename T>
inline void camera_fps(
    frame<T, 3>& frame, const vec<T, 3>& transl, const vec<T, 2>& rotate) {
    // https://gamedev.stackexchange.com/questions/30644/how-to-keep-my-quaternion-using-fps-camera-from-tilting-and-messing-up
    auto y = vec<T, 3>{0, 1, 0};
    auto z = orthonormalize(frame.z, y);
    auto x = cross(y, z);

    auto rot = rotation_frame(vec<T, 3>{1, 0, 0}, rotate.y) *
               ygl::frame<T, 3>{frame.x, frame.y, frame.z, vec<T, 3>{0, 0, 0}} *
               rotation_frame(vec<T, 3>{0, 1, 0}, rotate.x);
    auto pos = frame.o + transl.x * x + transl.y * y + transl.z * z;

    frame = {rot.x, rot.y, rot.z, pos};
}

}  // namespace ygl

#endif
