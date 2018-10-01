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
// in graphics. In particular, we support 2-4 dimensional vectors `vec2<T>`,
// `vec3f<T>`, `vec4<T>`.
//
// We support 2-4 dimensional generic matrices `mat2<T>`, `mat3<T>`,
// `mat4<T>`, with matrix-matrix and matrix-vector products, transposes and
// inverses. Matrices are stored in column-major ordered and are accessed and
// constructed by column.
//
// To represent transformations, most of the library facilities prefer the use
// coordinate frames, aka rigid transforms, represented as `frame2<T>`,
// `frame3<T>`. The structure store three coordinate axes and the origin.
// This is equivalent to a rigid transform written as a column-major affine
// matrix. Transform operations are better behaved with this representation.
//
// We represent coordinate bounds with axis-aligned bounding boxes with
// `bbox1<T>`, `bbox2<T>`, `bbox3<T>`, `bbox4<T>`, with support for
// expansion operations for points and other bounding boxes. We provide
// operations to compute bounds for points, lines, triangles and quads.
//
// For both matrices and frames we support transform operations for points,
// vectors and directions (`transform_point()`, `transform_vector()`,
// `transform_direction()`). For frames we also the support inverse operations
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
// minimal, but fully-featured, global illumination renderer, rather than the
// need of generic image editing.
//
// 0. load and save image with Yocto/GLIO
// 1. create images with `image<T>` data structure
// 2. resize images with `resize_image4f()`
// 3. tonemap images with `tonemap_image4f()`
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
//    - build the ray-tracing acceleration structure with `build_bvh()`
//     - prepare lights for rendering with `init_lights()`
// 2. create the image buffer and random number generators `make_trace_rngs()`
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

#ifndef YGL_EMBREE
#define YGL_EMBREE 1
#endif

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <algorithm>  // for std::upper_bound
#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <functional>  // for std::hash
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
using std::isfinite;
using std::log;
using std::pow;
using std::round;
using std::sin;
using std::sqrt;
using std::tan;

using byte = unsigned char;
using uint = unsigned int;

template<typename T>
inline T min(const T& x, const T& y) { return (x < y) ? x : y; }
template<typename T>
inline T max(const T& x, const T& y) { return (x > y) ? x : y; }
template<typename T>
inline T clamp(const T& x, const T& min_, const T& max_) { return min(max(x, min_), max_); }
template<typename T, typename T1>
inline T lerp(const T& a, const T& b, T1 u) { return a * (1 - u) + b * u; }

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
const T pi = (T)3.14159265358979323846;
const auto pif = 3.14159265f;

}  // namespace ygl

// -----------------------------------------------------------------------------
// VECTORS
// -----------------------------------------------------------------------------
namespace ygl {

// Small size vectors.
template <typename T>
struct vec2 {
    T x = 0;
    T y = 0;
};
template <typename T>
struct vec3 {
    T x = 0;
    T y = 0;
    T z = 0;
};
template <typename T>
struct vec4 {
    T x = 0;
    T y = 0;
    T z = 0;
    T w = 0;
};

// Type aliases.
using vec2f = vec2<float>;
using vec3f = vec3<float>;
using vec4f = vec4<float>;
using vec2i = vec2<int>;
using vec3i = vec3<int>;
using vec4i = vec4<int>;
using vec4b = vec4<byte>;

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
inline const T& at(const vec3<T>& v, int i) {
    return *(&v.x + i);
}
template <typename T>
inline T& at(vec3<T>& v, int i) {
    return *(&v.x + i);
}

// Access xyz component of a vec4 typically used for color operation.
template <typename T>
inline vec3<T>& xyz(const vec4<T>& a) {
    return (vec3<T>&)a;
}
template <typename T>
inline vec3<T>& xyz(vec4<T>& a) {
    return (vec3<T>&)a;
}

// Vector comparison operations.
template <typename T>
inline bool operator==(const vec2<T>& a, const vec2<T>& b) {
    return a.x == b.x && a.y == b.y;
}
template <typename T>
inline bool operator!=(const vec2<T>& a, const vec2<T>& b) {
    return a.x != b.x || a.y != b.y;
}
template <typename T, typename T1>
inline bool operator==(const vec2<T>& a, T1 b) {
    return a.x == b && a.y == b;
}
template <typename T, typename T1>
inline bool operator!=(const vec2<T>& a, T1 b) {
    return a.x != b || a.y != b;
}
template <typename T>
inline bool operator==(const vec3<T>& a, const vec3<T>& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}
template <typename T>
inline bool operator!=(const vec3<T>& a, const vec3<T>& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z;
}
template <typename T, typename T1>
inline bool operator==(const vec3<T>& a, T1 b) {
    return a.x == b && a.y == b && a.z == b;
}
template <typename T, typename T1>
inline bool operator!=(const vec3<T>& a, T1 b) {
    return a.x != b || a.y != b || a.z != b;
}
template <typename T>
inline bool operator==(const vec4<T>& a, const vec4<T>& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}
template <typename T>
inline bool operator!=(const vec4<T>& a, const vec4<T>& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w;
}
template <typename T, typename T1>
inline bool operator==(const vec4<T>& a, T1 b) {
    return a.x == b && a.y == b && a.z == b && a.w == b;
}
template <typename T, typename T1>
inline bool operator!=(const vec4<T>& a, T1 b) {
    return a.x != b || a.y != b || a.z != b || a.w != b;
}

// Vector operations.
template <typename T>
inline vec2<T> operator-(const vec2<T>& a) {
    return {-a.x, -a.y};
}
template <typename T>
inline vec2<T> operator+(const vec2<T>& a, const vec2<T>& b) {
    return {a.x + b.x, a.y + b.y};
}
template <typename T, typename T1>
inline vec2<T> operator+(const vec2<T>& a, T1 b) {
    return {a.x + b, a.y + b};
}
template <typename T, typename T1>
inline vec2<T> operator+(T1 a, const vec2<T>& b) {
    return {a + b.x, a + b.y};
}
template <typename T>
inline vec2<T> operator-(const vec2<T>& a, const vec2<T>& b) {
    return {a.x - b.x, a.y - b.y};
}
template <typename T, typename T1>
inline vec2<T> operator-(const vec2<T>& a, T1 b) {
    return {a.x - b, a.y - b};
}
template <typename T, typename T1>
inline vec2<T> operator-(T1 a, const vec2<T>& b) {
    return {a - b.x, a - b.y};
}
template <typename T>
inline vec2<T> operator*(const vec2<T>& a, const vec2<T>& b) {
    return {a.x * b.x, a.y * b.y};
}
template <typename T, typename T1>
inline vec2<T> operator*(const vec2<T>& a, T1 b) {
    return {a.x * b, a.y * b};
}
template <typename T, typename T1>
inline vec2<T> operator*(T1 a, const vec2<T>& b) {
    return {a * b.x, a * b.y};
}
template <typename T>
inline vec2<T> operator/(const vec2<T>& a, const vec2<T>& b) {
    return {a.x / b.x, a.y / b.y};
}
template <typename T, typename T1>
inline vec2<T> operator/(const vec2<T>& a, T1 b) {
    return {a.x / b, a.y / b};
}
template <typename T, typename T1>
inline vec2<T> operator/(T1 a, const vec2<T>& b) {
    return {a / b.x, a / b.y};
}

// Vector operations.
template <typename T>
inline vec3<T> operator+(const vec3<T>& a) {
    return a;
}
template <typename T>
inline vec3<T> operator-(const vec3<T>& a) {
    return {-a.x, -a.y, -a.z};
}
template <typename T>
inline vec3<T> operator+(const vec3<T>& a, const vec3<T>& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
template <typename T, typename T1>
inline vec3<T> operator+(const vec3<T>& a, T1 b) {
    return {a.x + b, a.y + b, a.z + b};
}
template <typename T, typename T1>
inline vec3<T> operator+(T1 a, const vec3<T>& b) {
    return {a + b.x, a + b.y, a + b.z};
}
template <typename T>
inline vec3<T> operator-(const vec3<T>& a, const vec3<T>& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
template <typename T, typename T1>
inline vec3<T> operator-(const vec3<T>& a, T1 b) {
    return {a.x - b, a.y - b, a.z - b};
}
template <typename T, typename T1>
inline vec3<T> operator-(T1 a, const vec3<T>& b) {
    return {a - b.x, a - b.y, a - b.z};
}
template <typename T>
inline vec3<T> operator*(const vec3<T>& a, const vec3<T>& b) {
    return {a.x * b.x, a.y * b.y, a.z * b.z};
}
template <typename T, typename T1>
inline vec3<T> operator*(const vec3<T>& a, T1 b) {
    return {a.x * b, a.y * b, a.z * b};
}
template <typename T, typename T1>
inline vec3<T> operator*(T1 a, const vec3<T>& b) {
    return {a * b.x, a * b.y, a * b.z};
}
template <typename T>
inline vec3<T> operator/(const vec3<T>& a, const vec3<T>& b) {
    return {a.x / b.x, a.y / b.y, a.z / b.z};
}
template <typename T, typename T1>
inline vec3<T> operator/(const vec3<T>& a, T1 b) {
    return {a.x / b, a.y / b, a.z / b};
}
template <typename T, typename T1>
inline vec3<T> operator/(T1 a, const vec3<T>& b) {
    return {a / b.x, a / b.y, a / b.z};
}

// Vector operations.
template <typename T>
inline vec4<T> operator-(const vec4<T>& a) {
    return {-a.x, -a.y, -a.z, -a.w};
}
template <typename T>
inline vec4<T> operator+(const vec4<T>& a, const vec4<T>& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}
template <typename T, typename T1>
inline vec4<T> operator+(const vec4<T>& a, T1 b) {
    return {a.x + b, a.y + b, a.z + b, a.w + b};
}
template <typename T, typename T1>
inline vec4<T> operator+(T1 a, const vec4<T>& b) {
    return {a + b.x, a + b.y, a + b.z, a + b.w};
}
template <typename T>
inline vec4<T> operator-(const vec4<T>& a, const vec4<T>& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w};
}
template <typename T, typename T1>
inline vec4<T> operator-(const vec4<T>& a, T1 b) {
    return {a.x - b, a.y - b, a.z - b, a.w - b};
}
template <typename T, typename T1>
inline vec4<T> operator-(T1 a, const vec4<T>& b) {
    return {a - b.x, a - b.y, a - b.z, a - b.w};
}
template <typename T>
inline vec4<T> operator*(const vec4<T>& a, const vec4<T>& b) {
    return {a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w};
}
template <typename T, typename T1>
inline vec4<T> operator*(const vec4<T>& a, T1 b) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
}
template <typename T, typename T1>
inline vec4<T> operator*(T1 a, const vec4<T>& b) {
    return {a * b.x, a * b.y, a * b.z, a * b.w};
}
template <typename T>
inline vec4<T> operator/(const vec4<T>& a, const vec4<T>& b) {
    return {a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w};
}
template <typename T, typename T1>
inline vec4<T> operator/(const vec4<T>& a, T1 b) {
    return {a.x / b, a.y / b, a.z / b, a.w / b};
}
template <typename T, typename T1>
inline vec4<T> operator/(T1 a, const vec4<T>& b) {
    return {a / b.x, a / b.y, a / b.z, a / b.w};
}

// Vector assignments
template <typename T>
inline vec2<T>& operator+=(vec2<T>& a, const vec2<T>& b) {
    return a = a + b;
}
template <typename T, typename T1>
inline vec2<T>& operator+=(vec2<T>& a, T1 b) {
    return a = a + b;
}
template <typename T>
inline vec2<T>& operator-=(vec2<T>& a, const vec2<T>& b) {
    return a = a - b;
}
template <typename T, typename T1>
inline vec2<T>& operator-=(vec2<T>& a, T1 b) {
    return a = a - b;
}
template <typename T>
inline vec2<T>& operator*=(vec2<T>& a, const vec2<T>& b) {
    return a = a * b;
}
template <typename T, typename T1>
inline vec2<T>& operator*=(vec2<T>& a, T1 b) {
    return a = a * b;
}
template <typename T>
inline vec2<T>& operator/=(vec2<T>& a, const vec2<T>& b) {
    return a = a / b;
}
template <typename T, typename T1>
inline vec2<T>& operator/=(vec2<T>& a, T1 b) {
    return a = a / b;
}

// Vector assignments
template <typename T>
inline vec3<T>& operator+=(vec3<T>& a, const vec3<T>& b) {
    return a = a + b;
}
template <typename T>
inline vec3<T>& operator-=(vec3<T>& a, const vec3<T>& b) {
    return a = a - b;
}
template <typename T>
inline vec3<T>& operator*=(vec3<T>& a, const vec3<T>& b) {
    return a = a * b;
}
template <typename T, typename T1>
inline vec3<T>& operator*=(vec3<T>& a, T1 b) {
    return a = a * b;
}
template <typename T>
inline vec3<T>& operator/=(vec3<T>& a, const vec3<T>& b) {
    return a = a / b;
}
template <typename T, typename T1>
inline vec3<T>& operator/=(vec3<T>& a, T1 b) {
    return a = a / b;
}

// Vector assignments
template <typename T>
inline vec4<T>& operator+=(vec4<T>& a, const vec4<T>& b) {
    return a = a + b;
}
template <typename T>
inline vec4<T>& operator-=(vec4<T>& a, const vec4<T>& b) {
    return a = a - b;
}
template <typename T>
inline vec4<T>& operator*=(vec4<T>& a, const vec4<T>& b) {
    return a = a * b;
}
template <typename T, typename T1>
inline vec4<T>& operator*=(vec4<T>& a, T1 b) {
    return a = a * b;
}
template <typename T>
inline vec4<T>& operator/=(vec4<T>& a, const vec4<T>& b) {
    return a = a / b;
}
template <typename T, typename T1>
inline vec4<T>& operator/=(vec4<T>& a, T1 b) {
    return a = a / b;
}

// Vector products and lengths.
template <typename T>
inline T dot(const vec2<T>& a, const vec2<T>& b) {
    return a.x * b.x + a.y * b.y;
}
template <typename T>
inline T length(const vec2<T>& a) {
    return sqrt(a.x * a.x + a.y * a.y);
}
template <typename T>
inline vec2<T> normalize(const vec2<T>& a) {
    auto l = length(a);
    return (l) ? a / l : a;
}
template <typename T>
inline T cross(const vec2<T>& a, const vec2<T>& b) {
    return a.x * b.y - a.y * b.x;
}
template <typename T>
inline T dot(const vec3<T>& a, const vec3<T>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
template <typename T>
inline T length(const vec3<T>& a) {
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}
template <typename T>
inline vec3<T> normalize(const vec3<T>& a) {
    auto l = length(a);
    return (l) ? a / l : a;
}
template <typename T>
inline vec3<T> cross(const vec3<T>& a, const vec3<T>& b) {
    return {
        a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
template <typename T>
inline T dot(const vec4<T>& a, const vec4<T>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}
template <typename T>
inline T length(const vec4<T>& a) {
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z + a.w * a.w);
}
template <typename T>
inline vec4<T> normalize(const vec4<T>& a) {
    auto l = length(a);
    return (l) ? a / l : a;
}

// Vector angles and slerps.
template <typename T>
inline T angle(const vec3<T>& a, const vec3<T>& b) {
    return acos(clamp(dot(normalize(a), normalize(b)), -1.0f, 1.0f));
}
template <typename T, typename T1>
inline vec4<T> slerp(const vec4<T>& a, const vec4<T>& b, T1 u) {
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
template <typename T>
inline vec3<T> orthogonal(const vec3<T>& v) {
    // http://lolengine.net/blog/2013/09/21/picking-orthogonal-vector-combing-coconuts)
    return fabs(v.x) > fabs(v.z) ? vec3<T>{-v.y, v.x, 0} :
                                   vec3<T>{0, -v.z, v.y};
}
template <typename T>
inline vec3<T> orthonormalize(const vec3<T>& a, const vec3<T>& b) {
    return normalize(a - b * dot(a, b));
}

// Reflected and refracted vector.
template <typename T>
inline vec3<T> reflect(const vec3<T>& w, const vec3<T>& n) {
    return -w + 2 * dot(n, w) * n;
}
template <typename T, typename T1>
inline vec3<T> refract(const vec3<T>& w, const vec3<T>& n, T1 eta) {
    // auto k = 1.0 - eta * eta * (1.0 - dot(n, w) * dot(n, w));
    auto k = 1 - eta * eta * max(0.0f, 1 - dot(n, w) * dot(n, w));
    if (k < 0) return {0, 0, 0};  // tir
    return -w * eta + (eta * dot(n, w) - sqrt(k)) * n;
}

// Max element and clamp.
template <typename T, typename T1, typename T2>
inline vec2<T> clamp(const vec2<T>& x, T1 min, T2 max) {
    return {clamp(x.x, min, max), clamp(x.y, min, max)};
}
template <typename T, typename T1, typename T2>
inline vec3<T> clamp(const vec3<T>& x, T1 min, T2 max) {
    return {clamp(x.x, min, max), clamp(x.y, min, max), clamp(x.z, min, max)};
}
template <typename T, typename T1, typename T2>
inline vec4<T> clamp(const vec4<T>& x, T1 min, T2 max) {
    return {clamp(x.x, min, max), clamp(x.y, min, max), clamp(x.z, min, max),
        clamp(x.w, min, max)};
}
template <typename T>
inline T max(const vec2<T>& a) {
    return max(a.x, a.y);
}
template <typename T>
inline T max(const vec3<T>& a) {
    return max(max(a.x, a.y), a.z);
}
template <typename T>
inline T max(const vec4<T>& a) {
    return max(max(max(a.x, a.y), a.z), a.w);
}
template <typename T>
inline T min(const vec2<T>& a) {
    return min(a.x, a.y);
}
template <typename T>
inline T min(const vec3<T>& a) {
    return min(min(a.x, a.y), a.z);
}
template <typename T>
inline T min(const vec4<T>& a) {
    return min(min(min(a.x, a.y), a.z), a.w);
}

// Quaternion operatons represented as xi + yj + zk + w
const auto identity_quat4f = vec4f{0, 0, 0, 1};
template <typename T>
inline vec4<T> quat_mul(const vec4<T>& a, float b) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
}
template <typename T>
inline vec4<T> quat_mul(const vec4<T>& a, const vec4<T>& b) {
    return {a.x * b.w + a.w * b.x + a.y * b.w - a.z * b.y,
        a.y * b.w + a.w * b.y + a.z * b.x - a.x * b.z,
        a.z * b.w + a.w * b.z + a.x * b.y - a.y * b.x,
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z};
}
template <typename T>
inline vec4<T> quat_conjugate(const vec4<T>& a) {
    return {-a.x, -a.y, -a.z, a.w};
}
template <typename T>
inline vec4<T> quat_inverse(const vec4<T>& a) {
    return quat_conjugate(a) / dot(a, a);
}

}  // namespace ygl

namespace std {

// Hash functor for vector for use with unordered_map
template <typename T>
struct hash<ygl::vec2<T>> {
    size_t operator()(const ygl::vec2<T>& v) const {
        auto vh = hash<T>();
        auto h = (size_t)0;
        for (auto i = 0; i < 2; i++)
            h ^= vh((&v.x)[i]) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};
// Hash functor for vector for use with unordered_map
template <typename T>
struct hash<ygl::vec3<T>> {
    size_t operator()(const ygl::vec3<T>& v) const {
        auto vh = hash<int>();
        auto h = (size_t)0;
        for (auto i = 0; i < 3; i++)
            h ^= vh((&v.x)[i]) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};
// Hash functor for vector for use with unordered_map
template <typename T>
struct hash<ygl::vec4<T>> {
    size_t operator()(const ygl::vec4<T>& v) const {
        auto vh = hash<int>();
        auto h = (size_t)0;
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
template <typename T>
struct mat2 {
    vec2<T> x = {1, 0};
    vec2<T> y = {0, 1};
};
template <typename T>
struct mat3 {
    vec3<T> x = {1, 0, 0};
    vec3<T> y = {0, 1, 0};
    vec3<T> z = {0, 0, 1};
};
template <typename T>
struct mat4 {
    vec4<T> x = {1, 0, 0, 0};
    vec4<T> y = {0, 1, 0, 0};
    vec4<T> z = {0, 0, 1, 0};
    vec4<T> w = {0, 0, 0, 1};
};

// Type aliases.
using mat2f = mat2<float>;
using mat3f = mat3<float>;
using mat4f = mat4<float>;

// Identity matrices constants.
const auto identity_mat2f = mat2f();
const auto identity_mat3f = mat3f();
const auto identity_mat4f = mat4f();

// Matrix comparisons.
template <typename T>
inline bool operator==(const mat2<T>& a, const mat2<T>& b) {
    return a.x == b.x && a.y == b.y;
}
template <typename T>
inline bool operator!=(const mat2<T>& a, const mat2<T>& b) {
    return !(a == b);
}
template <typename T>
inline bool operator==(const mat3<T>& a, const mat3<T>& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}
template <typename T>
inline bool operator!=(const mat3<T>& a, const mat3<T>& b) {
    return !(a == b);
}
template <typename T>
inline bool operator==(const mat4<T>& a, const mat4<T>& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}
template <typename T>
inline bool operator!=(const mat4<T>& a, const mat4<T>& b) {
    return !(a == b);
}

// Matrix operations.
template <typename T>
inline mat2<T> operator+(const mat2<T>& a, const mat2<T>& b) {
    return {a.x + b.x, a.y + b.y};
}
template <typename T, typename T1>
inline mat2<T> operator*(const mat2<T>& a, T1 b) {
    return {a.x * b, a.y * b};
}
template <typename T>
inline vec2<T> operator*(const mat2<T>& a, const vec2f& b) {
    return a.x * b.x + a.y * b.y;
}
template <typename T>
inline vec2<T> operator*(const vec2f& a, const mat2<T>& b) {
    return {dot(a, b.x), dot(a, b.y)};
}
template <typename T>
inline mat2<T> operator*(const mat2<T>& a, const mat2<T>& b) {
    return {a * b.x, a * b.y};
}

// Matrix operations.
template <typename T>
inline mat3<T> operator+(const mat3<T>& a, const mat3<T>& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
template <typename T, typename T1>
inline mat3<T> operator*(const mat3<T>& a, T1 b) {
    return {a.x * b, a.y * b, a.z * b};
}
template <typename T>
inline vec3<T> operator*(const mat3<T>& a, const vec3<T>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
template <typename T>
inline vec3f operator*(const vec3f& a, const mat3<T>& b) {
    return {dot(a, b.x), dot(a, b.y), dot(a, b.z)};
}
template <typename T>
inline mat3<T> operator*(const mat3<T>& a, const mat3<T>& b) {
    return {a * b.x, a * b.y, a * b.z};
}

// Matrix operations.
template <typename T>
inline mat4<T> operator+(const mat4<T>& a, const mat4<T>& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}
template <typename T, typename T1>
inline mat4<T> operator*(const mat4<T>& a, T1 b) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
}
template <typename T>
inline vec4<T> operator*(const mat4<T>& a, const vec4<T>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}
template <typename T>
inline vec4<T> operator*(const vec4<T>& a, const mat4<T>& b) {
    return {dot(a, b.x), dot(a, b.y), dot(a, b.z), dot(a, b.w)};
}
template <typename T>
inline mat4<T> operator*(const mat4<T>& a, const mat4<T>& b) {
    return {a * b.x, a * b.y, a * b.z, a * b.w};
}

// Matrix assignments.
template <typename T>
inline mat2<T>& operator+=(mat2<T>& a, const mat2<T>& b) {
    return a = a + b;
}
template <typename T>
inline mat2<T>& operator*=(mat2<T>& a, const mat2<T>& b) {
    return a = a * b;
}
template <typename T, typename T1>
inline mat2<T>& operator*=(mat2<T>& a, T1 b) {
    return a = a * b;
}

// Matrix assignments.
template <typename T>
inline mat3<T>& operator+=(mat3<T>& a, const mat3<T>& b) {
    return a = a + b;
}
template <typename T>
inline mat3<T>& operator*=(mat3<T>& a, const mat3<T>& b) {
    return a = a * b;
}
template <typename T, typename T1>
inline mat3<T>& operator*=(mat3<T>& a, T1 b) {
    return a = a * b;
}

// Matrix assignments.
template <typename T>
inline mat4<T>& operator+=(mat4<T>& a, const mat4<T>& b) {
    return a = a + b;
}
template <typename T>
inline mat4<T>& operator*=(mat4<T>& a, const mat4<T>& b) {
    return a = a * b;
}
template <typename T, typename T1>
inline mat4<T>& operator*=(mat4<T>& a, float b) {
    return a = a * b;
}

// Matrix diagonals and transposes.
template <typename T>
inline vec2<T> diagonal(const mat2<T>& a) {
    return {a.x.x, a.y.y};
}
template <typename T>
inline vec3<T> diagonal(const mat3<T>& a) {
    return {a.x.x, a.y.y, a.z.z};
}
template <typename T>
inline vec4<T> diagonal(const mat4<T>& a) {
    return {a.x.x, a.y.y, a.z.z, a.w.w};
}
template <typename T>
inline mat2<T> transpose(const mat2<T>& a);
template <typename T>
inline mat3<T> transpose(const mat3<T>& a);
template <typename T>
inline mat4<T> transpose(const mat4<T>& a);

// Matrix adjugates, determinant and inverses.
template <typename T>
inline mat2<T> adjugate(const mat2<T>& a);
template <typename T>
inline mat3<T> adjugate(const mat3<T>& a);
template <typename T>
inline mat4<T> adjugate(const mat4<T>& a);
template <typename T>
inline T determinant(const mat2<T>& a);
template <typename T>
inline T determinant(const mat3<T>& a);
template <typename T>
inline T determinant(const mat4<T>& a);
template <typename T>
inline mat2<T> inverse(const mat2<T>& a) {
    return adjugate(a) * (1 / determinant(a));
}
template <typename T>
inline mat3<T> inverse(const mat3<T>& a) {
    return adjugate(a) * (1 / determinant(a));
}
template <typename T>
inline mat4<T> inverse(const mat4<T>& a) {
    return adjugate(a) * (1 / determinant(a));
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// RIGID BODY TRANSFORMS/FRAMES
// -----------------------------------------------------------------------------
namespace ygl {

// Rigid frames stored as a column-major affine transform matrix.
template <typename T>
struct frame2 {
    vec2<T> x = {1, 0};
    vec2<T> y = {0, 1};
    vec2<T> o = {0, 0};
};
template <typename T>
struct frame3 {
    vec3<T> x = {1, 0, 0};
    vec3<T> y = {0, 1, 0};
    vec3<T> z = {0, 0, 1};
    vec3<T> o = {0, 0, 0};
};

// Type aliases.
using frame2f = frame2<float>;
using frame3f = frame3<float>;

// Indentity frames.
const auto identity_frame2f = frame2f{{1, 0}, {0, 1}, {0, 0}};
const auto identity_frame3f =
    frame3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};

// Frame construction from axis.
template <typename T>
inline frame3<T> make_frame_fromz(const vec3<T>& o, const vec3<T>& v) {
    auto z = normalize(v);
    auto x = normalize(orthogonal(z));
    auto y = normalize(cross(z, x));
    return {x, y, z, o};
}
template <typename T>
inline frame3<T> make_frame_fromzx(
    const vec3<T>& o, const vec3<T>& z_, const vec3<T>& x_) {
    auto z = normalize(z_);
    auto x = orthonormalize(x_, z);
    auto y = normalize(cross(z, x));
    return {x, y, z, o};
}

// Frame to matrix conversion.
template <typename T>
inline mat4<T> frame_to_mat(const frame3<T>& a) {
    return {
        {a.x.x, a.x.y, a.x.z, 0},
        {a.y.x, a.y.y, a.y.z, 0},
        {a.z.x, a.z.y, a.z.z, 0},
        {a.o.x, a.o.y, a.o.z, 1},
    };
}
template <typename T>
inline frame3<T> mat_to_frame(const mat4<T>& a) {
    return {
        {a.x.x, a.x.y, a.x.z},
        {a.y.x, a.y.y, a.y.z},
        {a.z.x, a.z.y, a.z.z},
        {a.w.x, a.w.y, a.w.z},
    };
}

// Frame comparisons.
template <typename T>
inline bool operator==(const frame2<T>& a, const frame2<T>& b) {
    return a.x == b.x && a.y == b.y && a.o == b.o;
}
template <typename T>
inline bool operator!=(const frame2<T>& a, const frame2<T>& b) {
    return !(a == b);
}
template <typename T>
inline bool operator==(const frame3<T>& a, const frame3<T>& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.o == b.o;
}
template <typename T>
inline bool operator!=(const frame3<T>& a, const frame3<T>& b) {
    return !(a == b);
}

// Frame composition, equivalent to affine matrix product.
template <typename T>
inline frame2<T> operator*(const frame2<T>& a, const frame2<T>& b) {
    auto rot = mat2<T>{a.x, a.y} * mat2<T>{b.x, b.y};
    auto pos = mat2<T>{a.x, a.y} * b.o + a.o;
    return {rot.x, rot.y, pos};
}
template <typename T>
inline frame3<T> operator*(const frame3<T>& a, const frame3<T>& b) {
    auto rot = mat3<T>{a.x, a.y, a.z} * mat3<T>{b.x, b.y, b.z};
    auto pos = mat3<T>{a.x, a.y, a.z} * b.o + a.o;
    return {rot.x, rot.y, rot.z, pos};
}
// Frame inverse, equivalent to rigid affine inverse.
template <typename T>
inline frame2<T> inverse(const frame2<T>& a, bool is_rigid = true) {
    auto minv =
        (is_rigid) ? transpose(mat2<T>{a.x, a.y}) : inverse(mat2<T>{a.x, a.y});
    return {minv.x, minv.y, -(minv * a.o)};
}
template <typename T>
inline frame3<T> inverse(const frame3<T>& a, bool is_rigid = true) {
    auto minv = (is_rigid) ? transpose(mat3<T>{a.x, a.y, a.z}) :
                             inverse(mat3<T>{a.x, a.y, a.z});
    return {minv.x, minv.y, minv.z, -(minv * a.o)};
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// AXIS ALIGNED BOUNDING BOXES
// -----------------------------------------------------------------------------
namespace ygl {

// Range of values in 1D.
template <typename T>
struct bbox1 {
    T min = maxt<T>();
    T max = mint<T>();
};

// Axis aligned bounding box represented as a min/max vector pairs.
template <typename T>
struct bbox2 {
    vec2<T> min = {maxt<T>(), maxt<T>()};
    vec2<T> max = {mint<T>(), mint<T>()};
};

// Axis aligned bounding box represented as a min/max vector pairs.
template <typename T>
struct bbox3 {
    vec3<T> min = {maxt<T>(), maxt<T>(), maxt<T>()};
    vec3<T> max = {mint<T>(), mint<T>(), mint<T>()};
};

// Axis aligned bounding box represented as a min/max vector pairs.
template <typename T>
struct bbox4 {
    vec4<T> min = {maxt<T>(), maxt<T>(), maxt<T>(), maxt<T>()};
    vec4<T> max = {mint<T>(), mint<T>(), mint<T>(), mint<T>()};
};

// Type aliases
using bbox1f = bbox1<float>;
using bbox2f = bbox2<float>;
using bbox3f = bbox3<float>;
using bbox4f = bbox4<float>;

// Empty bbox constant.
const auto invalid_bbox1f = bbox1f();
const auto invalid_bbox2f = bbox2f();
const auto invalid_bbox3f = bbox3f();
const auto invalid_bbox4f = bbox4f();

// Bounding box comparisons.
template <typename T>
inline bool operator==(const bbox1<T>& a, const bbox1<T>& b) {
    return a.min == b.min && a.max == b.max;
}
template <typename T>
inline bool operator!=(const bbox1<T>& a, const bbox1<T>& b) {
    return a.min != b.min || a.max != b.max;
}
template <typename T>
inline bool operator==(const bbox2<T>& a, const bbox2<T>& b) {
    return a.min == b.min && a.max == b.max;
}
template <typename T>
inline bool operator!=(const bbox2<T>& a, const bbox2<T>& b) {
    return a.min != b.min || a.max != b.max;
}
template <typename T>
inline bool operator==(const bbox3<T>& a, const bbox3<T>& b) {
    return a.min == b.min && a.max == b.max;
}
template <typename T>
inline bool operator!=(const bbox3<T>& a, const bbox3<T>& b) {
    return a.min != b.min || a.max != b.max;
}
template <typename T>
inline bool operator==(const bbox4<T>& a, const bbox4<T>& b) {
    return a.min == b.min && a.max == b.max;
}
template <typename T>
inline bool operator!=(const bbox4<T>& a, const bbox4<T>& b) {
    return a.min != b.min || a.max != b.max;
}

// Bounding box expansions with points and other boxes.
template <typename T>
inline bbox1<T>& operator+=(bbox1<T>& a, T b) {
    a.min = min(a.min, b);
    a.max = max(a.max, b);
    return a;
}
template <typename T>
inline bbox1<T>& operator+=(bbox1<T>& a, const bbox1<T>& b) {
    a.min = min(a.min, b.min);
    a.max = max(a.max, b.max);
    return a;
}
// Bounding box expansions with points and other boxes.
template <typename T>
inline bbox2<T>& operator+=(bbox2<T>& a, const vec2<T>& b) {
    a.min = {min(a.min.x, b.x), min(a.min.y, b.y)};
    a.max = {max(a.max.x, b.x), max(a.max.y, b.y)};
    return a;
}
template <typename T>
inline bbox2<T>& operator+=(bbox2<T>& a, const bbox2<T>& b) {
    a.min = {min(a.min.x, b.min.x), min(a.min.y, b.min.y)};
    a.max = {max(a.max.x, b.max.x), max(a.max.y, b.max.y)};
    return a;
}
// Bounding box expansions with points and other boxes.
template <typename T>
inline bbox3<T>& operator+=(bbox3<T>& a, const vec3<T>& b) {
    a.min = {min(a.min.x, b.x), min(a.min.y, b.y), min(a.min.z, b.z)};
    a.max = {max(a.max.x, b.x), max(a.max.y, b.y), max(a.max.z, b.z)};
    return a;
}
template <typename T>
inline bbox3<T>& operator+=(bbox3<T>& a, const bbox3<T>& b) {
    a.min = {
        min(a.min.x, b.min.x), min(a.min.y, b.min.y), min(a.min.z, b.min.z)};
    a.max = {
        max(a.max.x, b.max.x), max(a.max.y, b.max.y), max(a.max.z, b.max.z)};
    return a;
}
// Bounding box expansions with points and other boxes.
template <typename T>
inline bbox4<T>& operator+=(bbox4<T>& a, const vec4<T>& b) {
    a.min = {min(a.min.x, b.x), min(a.min.y, b.y), min(a.min.z, b.z),
        min(a.min.w, b.w)};
    a.max = {max(a.max.x, b.x), max(a.max.y, b.y), max(a.max.z, b.z),
        max(a.max.w, b.w)};
    return a;
}
template <typename T>
inline bbox4<T>& operator+=(bbox4<T>& a, const bbox4<T>& b) {
    a.min = {min(a.min.x, b.min.x), min(a.min.y, b.min.y),
        min(a.min.z, b.min.z), min(a.min.w, b.min.w)};
    a.max = {max(a.max.x, b.max.x), max(a.max.y, b.max.y),
        max(a.max.z, b.max.z), max(a.max.w, b.max.w)};
    return a;
}

// Primitive bounds.
template <typename T>
inline bbox3<T> point_bbox(const vec3<T>& p, T r = 0) {
    auto bbox = bbox3<T>{};
    bbox += p - vec3f{r, r, r};
    bbox += p + vec3f{r, r, r};
    return bbox;
}
template <typename T>
inline bbox3<T> line_bbox(
    const vec3<T>& p0, const vec3<T>& p1, T r0 = 0, T r1 = 0) {
    auto bbox = bbox3<T>{};
    bbox += p0 - vec3f{r0, r0, r0};
    bbox += p0 + vec3f{r0, r0, r0};
    bbox += p1 - vec3f{r1, r1, r1};
    bbox += p1 + vec3f{r1, r1, r1};
    return bbox;
}
template <typename T>
inline bbox3<T> triangle_bbox(
    const vec3<T>& p0, const vec3<T>& p1, const vec3<T>& p2) {
    auto bbox = bbox3<T>{};
    bbox += p0;
    bbox += p1;
    bbox += p2;
    return bbox;
}
template <typename T>
inline bbox3<T> quad_bbox(const vec3<T>& p0, const vec3<T>& p1,
    const vec3<T>& p2, const vec3<T>& p3) {
    auto bbox = bbox3<T>{};
    bbox += p0;
    bbox += p1;
    bbox += p2;
    bbox += p3;
    return bbox;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// RAYS
// -----------------------------------------------------------------------------
namespace ygl {

// Rays with origin, direction and min/max t value.
template <typename T>
struct ray2 {
    vec2<T> o = {0, 0};
    vec2<T> d = {0, 1};
    T tmin = 0;
    T tmax = maxt<T>();
};
template <typename T>
struct ray3 {
    vec3<T> o = {0, 0, 0};
    vec3<T> d = {0, 0, 1};
    T tmin = 0;
    T tmax = maxt<T>();
};

// Type aliases.
using ray2f = ray2<float>;
using ray3f = ray3<float>;

// Construct a ray from direction or segments using a default epsilon.
template <typename T>
inline ray3<T> make_ray(const vec3<T>& o, const vec3<T>& d, T eps = 1e-4f) {
    return {o, d, eps, maxt<T>()};
}
template <typename T>
inline ray3<T> make_segment(
    const vec3<T>& p1, const vec3<T>& p2, T eps = 1e-4f) {
    return {p1, normalize(p2 - p1), eps, length(p2 - p1) - 2 * eps};
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// TRANSFORMS
// -----------------------------------------------------------------------------
namespace ygl {

// Transforms points, vectors and directions by matrices.
template <typename T>
inline vec2<T> transform_point(const mat3<T>& a, const vec2<T>& b) {
    auto tvb = a * vec3f{b.x, b.y, 1};
    return vec2<T>{tvb.x, tvb.y} / tvb.z;
}
template <typename T>
inline vec3<T> transform_point(const mat4<T>& a, const vec3<T>& b) {
    auto tvb = a * vec4f{b.x, b.y, b.z, 1};
    return vec3<T>{tvb.x, tvb.y, tvb.z} / tvb.w;
}
template <typename T>
inline vec2<T> transform_vector(const mat3<T>& a, const vec2<T>& b) {
    auto tvb = a * vec3<T>{b.x, b.y, 0};
    return vec2<T>{tvb.x, tvb.y} / tvb.z;
}
template <typename T>
inline vec3<T> transform_vector(const mat3<T>& a, const vec3<T>& b) {
    return a * b;
}
template <typename T>
inline vec3<T> transform_vector(const mat4<T>& a, const vec3<T>& b) {
    auto tvb = a * vec4f{b.x, b.y, b.z, 0};
    return vec3<T>{tvb.x, tvb.y, tvb.z};
}
template <typename T>
inline vec3<T> transform_direction(const mat4<T>& a, const vec3<T>& b) {
    return normalize(transform_vector(a, b));
}

// Transforms points, vectors and directions by frames.
template <typename T>
inline vec2<T> transform_point(const frame2<T>& a, const vec2<T>& b) {
    return a.x * b.x + a.y * b.y + a.o;
}
template <typename T>
inline vec3<T> transform_point(const frame3<T>& a, const vec3<T>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.o;
}
template <typename T>
inline vec2<T> transform_vector(const frame2<T>& a, const vec2<T>& b) {
    return a.x * b.x + a.y * b.y;
}
template <typename T>
inline vec3<T> transform_vector(const frame3<T>& a, const vec3<T>& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
template <typename T>
inline vec3<T> transform_direction(const frame3<T>& a, const vec3<T>& b) {
    return normalize(transform_vector(a, b));
}

// Transforms rays and bounding boxes by matrices.
template <typename T>
inline ray3<T> transform_ray(const frame3<T>& a, const ray3<T>& b) {
    return {transform_point(a, b.o), transform_vector(a, b.d), b.tmin, b.tmax};
}
template <typename T>
inline ray3<T> transform_ray(const mat4<T>& a, const ray3<T>& b) {
    return {transform_point(a, b.o), transform_vector(a, b.d), b.tmin, b.tmax};
}
template <typename T>
inline bbox3<T> transform_bbox(const frame3<T>& a, const bbox3<T>& b) {
    auto corners = {vec3f{b.min.x, b.min.y, b.min.z},
        vec3f{b.min.x, b.min.y, b.max.z}, vec3f{b.min.x, b.max.y, b.min.z},
        vec3f{b.min.x, b.max.y, b.max.z}, vec3f{b.max.x, b.min.y, b.min.z},
        vec3f{b.max.x, b.min.y, b.max.z}, vec3f{b.max.x, b.max.y, b.min.z},
        vec3f{b.max.x, b.max.y, b.max.z}};
    auto xformed = bbox3<T>();
    for (auto& corner : corners) xformed += transform_point(a, corner);
    return xformed;
}
template <typename T>
inline bbox3<T> transform_bbox(const mat4<T>& a, const bbox3<T>& b) {
    auto corners = {vec3f{b.min.x, b.min.y, b.min.z},
        vec3f{b.min.x, b.min.y, b.max.z}, vec3f{b.min.x, b.max.y, b.min.z},
        vec3f{b.min.x, b.max.y, b.max.z}, vec3f{b.max.x, b.min.y, b.min.z},
        vec3f{b.max.x, b.min.y, b.max.z}, vec3f{b.max.x, b.max.y, b.min.z},
        vec3f{b.max.x, b.max.y, b.max.z}};
    auto xformed = bbox3<T>();
    for (auto& corner : corners) xformed += transform_point(a, corner);
    return xformed;
}

// Inverse transforms by frames, assuming they are rigid transforms.
template <typename T>
inline vec2<T> transform_point_inverse(const frame2<T>& a, const vec2<T>& b) {
    return {dot(b - a.o, a.x), dot(b - a.o, a.y)};
}
template <typename T>
inline vec3f transform_point_inverse(const frame3<T>& a, const vec3f& b) {
    return {dot(b - a.o, a.x), dot(b - a.o, a.y), dot(b - a.o, a.z)};
}
template <typename T>
inline vec2<T> transform_vector_inverse(const frame2<T>& a, const vec2<T>& b) {
    return {dot(b, a.x), dot(b, a.y)};
}
template <typename T>
inline vec3f transform_vector_inverse(const frame3<T>& a, const vec3f& b) {
    return {dot(b, a.x), dot(b, a.y), dot(b, a.z)};
}
template <typename T>
inline vec3f transform_direction_inverse(const frame3<T>& a, const vec3f& b) {
    return normalize(transform_vector_inverse(a, b));
}
template <typename T>
inline ray3<T> transform_ray_inverse(const frame3<T>& a, const ray3<T>& b) {
    return {transform_point_inverse(a, b.o),
        transform_direction_inverse(a, b.d), b.tmin, b.tmax};
}
template <typename T>
inline bbox3<T> transform_bbox_inverse(const frame3<T>& a, const bbox3<T>& b) {
    return transform_bbox(inverse(a), b);
}

// Translation, scaling and rotations transforms.
template <typename T>
inline frame3<T> translation_frame(const vec3<T>& a) {
    return {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, a};
}
template <typename T>
inline frame3<T> scaling_frame(const vec3<T>& a) {
    return {{a.x, 0, 0}, {0, a.y, 0}, {0, 0, a.z}, {0, 0, 0}};
}
template <typename T>
inline frame3<T> rotation_frame(const vec3<T>& axis, T angle) {
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
inline frame3<T> rotation_frame(const vec4<T>& quat) {
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
inline frame3<T> rotation_frame(const mat3<T>& rot) {
    return {rot.x, rot.y, rot.z, {0, 0, 0}};
}

// Lookat frame. Z-axis can be inverted with inv_xz.
template <typename T>
inline frame3<T> lookat_frame(const vec3<T>& eye, const vec3<T>& center,
    const vec3<T>& up, bool inv_xz = false) {
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
inline mat4<T> frustum_mat(T l, T r, T b, T t, T n, T f) {
    return {{2 * n / (r - l), 0, 0, 0}, {0, 2 * n / (t - b), 0, 0},
        {(r + l) / (r - l), (t + b) / (t - b), -(f + n) / (f - n), -1},
        {0, 0, -2 * f * n / (f - n), 0}};
}
template <typename T>
inline mat4<T> ortho_mat(T l, T r, T b, T t, T n, T f) {
    return {{2 / (r - l), 0, 0, 0}, {0, 2 / (t - b), 0, 0},
        {0, 0, -2 / (f - n), 0},
        {-(r + l) / (r - l), -(t + b) / (t - b), -(f + n) / (f - n), 1}};
}
template <typename T>
inline mat4<T> ortho2d_mat(T left, T right, T bottom, T top) {
    return ortho_mat(left, right, bottom, top, -1, 1);
}
template <typename T>
inline mat4<T> ortho_mat(T xmag, T ymag, T near, T far) {
    return {{1 / xmag, 0, 0, 0}, {0, 1 / ymag, 0, 0},
        {0, 0, 2 / (near - far), 0}, {0, 0, (far + near) / (near - far), 1}};
}
template <typename T>
inline mat4<T> perspective_mat(T fovy, T aspect, T near, T far) {
    auto tg = tan(fovy / 2);
    return {{1 / (aspect * tg), 0, 0, 0}, {0, 1 / tg, 0, 0},
        {0, 0, (far + near) / (near - far), -1},
        {0, 0, 2 * far * near / (near - far), 0}};
}
template <typename T>
inline mat4<T> perspective_mat(T fovy, T aspect, T near) {
    auto tg = tan(fovy / 2);
    return {{1 / (aspect * tg), 0, 0, 0}, {0, 1 / tg, 0, 0}, {0, 0, -1, -1},
        {0, 0, 2 * near, 0}};
}

// Rotation conversions.
template <typename T>
inline std::pair<vec3<T>, T> rotation_axisangle(const vec4<T>& quat) {
    return {normalize(vec3f{quat.x, quat.y, quat.z}), 2 * acos(quat.w)};
}
template <typename T>
inline vec4<T> rotation_quat(const vec3<T>& axis, T angle) {
    auto len = length(axis);
    if (!len) return {0, 0, 0, 1};
    return vec4f{sin(angle / 2) * axis.x / len, sin(angle / 2) * axis.y / len,
        sin(angle / 2) * axis.z / len, cos(angle / 2)};
}
template <typename T>
inline vec4<T> rotation_quat(const vec4<T>& axisangle) {
    return rotation_quat(
        vec3f{axisangle.x, axisangle.y, axisangle.z}, axisangle.w);
}

// Turntable and FPS Camera navigation.
inline void camera_turntable(vec3f& from, vec3f& to, vec3f& up,
    const vec2f& rotate, float dolly, const vec2f& pan);
inline void camera_turntable(frame3f& frame, float& focus, const vec2f& rotate,
    float dolly, const vec2f& pan);
inline void camera_fps(
    frame3f& frame, const vec3f& transl, const vec2f& rotate);

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
        scale = min(winsize.x / (float)imsize.x, winsize.y / (float)imsize.y);
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
inline vec3f sample_hemisphere(const vec2f& ruv) {
    auto z = ruv.y;
    auto r = sqrt(1 - z * z);
    auto phi = 2 * pif * ruv.x;
    return {r * cos(phi), r * sin(phi), z};
}
inline float sample_hemisphere_pdf(const vec3f& w) {
    return (w.z <= 0) ? 0 : 1 / (2 * pif);
}

// Sample a spherical direction with uniform distribution.
inline vec3f sample_sphere(const vec2f& ruv) {
    auto z = 2 * ruv.y - 1;
    auto r = sqrt(1 - z * z);
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
    auto z = sqrt(ruv.y);
    auto r = sqrt(1 - z * z);
    auto phi = 2 * pif * ruv.x;
    return {r * cos(phi), r * sin(phi), z};
}
inline float sample_hemisphere_cosine_pdf(const vec3f& w) {
    return (w.z <= 0) ? 0 : w.z / pif;
}

// Sample an hemispherical direction with cosine power distribution.
inline vec3f sample_hemisphere_cospower(float n, const vec2f& ruv) {
    auto z = pow(ruv.y, 1 / (n + 1));
    auto r = sqrt(1 - z * z);
    auto phi = 2 * pif * ruv.x;
    return {r * cos(phi), r * sin(phi), z};
}
inline float sample_hemisphere_cospower_pdf(float n, const vec3f& w) {
    return (w.z <= 0) ? 0 : pow(w.z, n) * (n + 1) / (2 * pif);
}

// Sample a point uniformly on a disk.
inline vec3f sample_disk(const vec2f& ruv) {
    auto r = sqrt(ruv.y);
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
inline int sample_discrete(const std::vector<float>& cdf, float r) {
    r = clamp(r * cdf.back(), 0.0f, cdf.back() - 0.00001f);
    auto idx = (int)(std::upper_bound(cdf.data(), cdf.data() + cdf.size(), r) -
                     cdf.data());
    return clamp(idx, 0, (int)cdf.size() - 1);
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
inline vec3f triangle_normal(
    const vec3f& p0, const vec3f& p1, const vec3f& p2) {
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
inline std::pair<vec3f, vec3f> triangle_tangents_fromuv(const vec3f& p0,
    const vec3f& p1, const vec3f& p2, const vec2f& uv0, const vec2f& uv1,
    const vec2f& uv2) {
    // Follows the definition in http://www.terathon.com/code/tangent.html and
    // https://gist.github.com/aras-p/2843984
    // normal points up from texture space
    auto p = p1 - p0;
    auto q = p2 - p0;
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
    return p0 * (1 - uv.x) * (1 - uv.y) + p1 * uv.x * (1 - uv.y) +
           p2 * uv.x * uv.y + p3 * (1 - uv.x) * uv.y;
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
std::vector<vec3f> compute_tangents(
    const std::vector<vec2i>& lines, const std::vector<vec3f>& pos);
std::vector<vec3f> compute_normals(
    const std::vector<vec3i>& triangles, const std::vector<vec3f>& pos);
std::vector<vec3f> compute_normals(
    const std::vector<vec4i>& quads, const std::vector<vec3f>& pos);

// Compute per-vertex tangent space for triangle meshes.
// Tangent space is defined by a four component vector.
// The first three components are the tangent with respect to the u texcoord.
// The fourth component is the sign of the tangent wrt the v texcoord.
// Tangent frame is useful in normal mapping.
std::vector<vec4f> compute_tangent_space(const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, const std::vector<vec3f>& norm,
    const std::vector<vec2f>& texcoord);

// Apply skinning to vertex position and normals.
std::pair<std::vector<vec3f>, std::vector<vec3f>> compute_skinning(
    const std::vector<vec3f>& pos, const std::vector<vec3f>& norm,
    const std::vector<vec4f>& weights, const std::vector<vec4i>& joints,
    const std::vector<frame3f>& xforms);
// Apply skinning as specified in Khronos glTF.
std::pair<std::vector<vec3f>, std::vector<vec3f>> compute_matrix_skinning(
    const std::vector<vec3f>& pos, const std::vector<vec3f>& norm,
    const std::vector<vec4f>& weights, const std::vector<vec4i>& joints,
    const std::vector<mat4f>& xforms);

// Dictionary to store edge information.
using edge_map = std::unordered_map<vec2i, vec2i>;

// Initialize an edge map with elements.
edge_map make_edge_map(const std::vector<vec3i>& triangles);
edge_map make_edge_map(const std::vector<vec4i>& quads);
// Insert an edge and return its index
int insert_edge(edge_map& emap, const vec2i& edge);
// Get the edge index / insertion count
int get_edge_index(const edge_map& emap, const vec2i& edge);
int get_edge_count(const edge_map& emap, const vec2i& edge);
// Get list of edges / boundary edges
std::vector<vec2i> get_edges(const edge_map& emap);
std::vector<vec2i> get_boundary(const edge_map& emap);

// Create an array of edges.
inline std::vector<vec2i> get_edges(const std::vector<vec3i>& triangles) {
    return get_edges(make_edge_map(triangles));
}
inline std::vector<vec2i> get_edges(const std::vector<vec4i>& quads) {
    return get_edges(make_edge_map(quads));
}

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
std::pair<std::vector<vec2i>, std::vector<T>> subdivide_lines(
    const std::vector<vec2i>& lines, const std::vector<T>& vert);
// Subdivide triangle by splitting each triangle in four, creating new
// vertices for each edge.
template <typename T>
std::pair<std::vector<vec3i>, std::vector<T>> subdivide_triangles(
    const std::vector<vec3i>& triangles, const std::vector<T>& vert);
// Subdivide quads by splitting each quads in four, creating new
// vertices for each edge and for each face.
template <typename T>
std::pair<std::vector<vec4i>, std::vector<T>> subdivide_quads(
    const std::vector<vec4i>& quads, const std::vector<T>& vert);
// Subdivide beziers by splitting each segment in two.
template <typename T>
std::pair<std::vector<vec4i>, std::vector<T>> subdivide_beziers(
    const std::vector<vec4i>& beziers, const std::vector<T>& vert);
// Subdivide quads using Carmull-Clark subdivision rules.
template <typename T>
std::pair<std::vector<vec4i>, std::vector<T>> subdivide_catmullclark(
    const std::vector<vec4i>& quads, const std::vector<T>& vert,
    bool lock_boundary = false);

// Weld vertices within a threshold. For noe the implementation is O(n^2).
std::pair<std::vector<vec3f>, std::vector<int>> weld_vertices(
    const std::vector<vec3f>& pos, float threshold);
std::pair<std::vector<vec3i>, std::vector<vec3f>> weld_triangles(
    const std::vector<vec3i>& triangles, const std::vector<vec3f>& pos,
    float threshold);
std::pair<std::vector<vec4i>, std::vector<vec3f>> weld_quads(
    const std::vector<vec4i>& quads, const std::vector<vec3f>& pos,
    float threshold);

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
// and texcoord of the sampled points.
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
    bbox3f bbox;                    // bounds
    uint32_t prims[bvh_max_prims];  // primitives
    uint16_t count;                 // number of prims
    bool internal;                  // whether it is an internal node or a leaf
    uint8_t split_axis;             // split axis
};

// Instance for a scene BVH.
struct bvh_instance {
    frame3f frame = identity_frame3f;      // frame
    frame3f frame_inv = identity_frame3f;  // frame inverse
    int sid = -1;                          // shape index
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
    std::vector<vec4i> quads;      // Quads for shape BVHs.

    // data for instance BVH
    std::vector<bvh_instance> instances;  // instances
    std::vector<bvh_tree*> shape_bvhs;    // shape bvhs

    // optional application specific data to key from a pointer to internal ids
    std::unordered_map<void*, int> instance_ids;
    std::unordered_map<void*, int> shape_ids;

    // bvh nodes
    std::vector<bvh_node> nodes;  // Internal nodes.

    // Embree opaque data
    void* embree_bvh = nullptr;

    // cleanup
    ~bvh_tree();
};

// Build a BVH from the given set of primitives.
void build_bvh(bvh_tree* bvh, bool high_quality);
// Update the node bounds for a shape bvh.
void refit_bvh(bvh_tree* bvh);

// Build a BVH from the given set of primitives.
// Uses Embree if available and requested, otherwise the standard build.
void build_bvh(bvh_tree* bvh, bool high_quality, bool embree);

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
    std::vector<vec3f> pos;       // positions
    std::vector<vec3f> norm;      // normals/tangents
    std::vector<vec2f> texcoord;  // texture coordinates
    std::vector<float> radius;    // radius for lines and points

    std::vector<int> points;       // points
    std::vector<vec2i> lines;      // lines
    std::vector<vec3i> triangles;  // triangles
    std::vector<vec4i> quads;      // quads
    std::vector<vec4i> beziers;    // beziers

    std::vector<vec4i> quads_pos;       // facevarying quads for pos
    std::vector<vec4i> quads_norm;      // facevarying quads for norm
    std::vector<vec4i> quads_texcoord;  // facevarying quads for texcoord
};

// Make examples shapes that are not watertight (besides quads).
// Return (triangles, quads, pos, norm, texcoord)
make_shape_data* make_quad(const vec2i& steps, const vec2f& size,
    const vec2f& uvsize, bool as_triangles);
make_shape_data* make_quad_stack(const vec3i& steps, const vec3f& size,
    const vec2f& uvsize, bool as_triangles);
make_shape_data* make_floor(const vec2i& steps, const vec2f& size,
    const vec2f& uvsize, bool as_triangles);
make_shape_data* make_cube(const vec3i& steps, const vec3f& size,
    const vec3f& uvsize, bool as_triangles);
make_shape_data* make_cube_rounded(const vec3i& steps, const vec3f& size,
    const vec3f& uvsize, float radius, bool as_triangles);
make_shape_data* make_sphere(
    const vec2i& steps, float size, const vec2f& uvsize, bool as_triangles);
make_shape_data* make_sphere_cube(
    int steps, float size, float uvsize, bool as_triangles);
make_shape_data* make_sphere_flipcap(const vec2i& steps, float size,
    const vec2f& uvsize, const vec2f& zflip, bool as_triangles);
make_shape_data* make_disk(
    const vec2i& steps, float size, const vec2f& uvsize, bool as_triangles);
make_shape_data* make_disk_quad(
    int steps, float size, float uvsize, bool as_triangles);
make_shape_data* make_disk_bulged(
    int steps, float size, float uvsize, float height, bool as_triangles);
make_shape_data* make_cylinder_side(const vec2i& steps, const vec2f& size,
    const vec2f& uvsize, bool as_triangles);
make_shape_data* make_cylinder(const vec3i& steps, const vec2f& size,
    const vec3f& uvsize, bool as_triangles);
make_shape_data* make_cylinder_rounded(const vec3i& steps, const vec2f& size,
    const vec3f& uvsize, float radius, bool as_triangles);
make_shape_data* make_geodesic_sphere(
    int tesselation, float size, bool as_triangles);

// Make examples shapes with are watertight (good for subdivs).
// Returns (triangles, quads, pos)
make_shape_data* make_suzanne(float size, bool as_triangles);
make_shape_data* make_cube(const vec3f& size, bool as_triangles);

// Make facevarying example shapes that are watertight (good for subdivs).
make_shape_data* make_fvcube(
    const vec3i& steps, const vec3f& size, const vec3f& uvsize);

// Generate lines set along a quad. Returns lines, pos, norm, texcoord, radius.
make_shape_data* make_lines(const vec2i& steps, const vec2f& size,
    const vec2f& uvsize, const vec2f& line_radius = {0.001f, 0.001f});

// Make point primitives. Returns points, pos, norm, texcoord, radius.
make_shape_data* make_point(float point_radius = 0.001f);
make_shape_data* make_points(
    int num, float uvsize, float point_radius = 0.001f);
make_shape_data* make_random_points(int num, const vec3f& size, float uvsize,
    float point_radius = 0.001f, uint64_t seed = 0);

// Make a bezier circle. Returns bezier, pos.
make_shape_data* make_bezier_circle(
    std::vector<vec4i>& beziers, std::vector<vec3f>& pos);

// Make a hair ball around a shape.  Returns lines, pos, norm, texcoord, radius.
// length: minimum and maximum length
// rad: minimum and maximum radius from base to tip
// noise: noise added to hair (strength/scale)
// clump: clump added to hair (number/strength)
// rotation: rotation added to hair (angle/strength)
make_shape_data* make_hair(const vec2i& steps,
    const std::vector<vec3i>& striangles, const std::vector<vec3f>& spos,
    const std::vector<vec3f>& snorm, const std::vector<vec2f>& stexcoord,
    const vec2f& length = {0.1f, 0.1f}, const vec2f& rad = {0.001f, 0.001f},
    const vec2f& noise = zero2f, const vec2f& clump = zero2f,
    const vec2f& rotation = zero2f, int seed = 7);

// Helper to concatenated shape data for non-facevarying shapes.
make_shape_data* merge_shape_data(const std::vector<make_shape_data*>& shapes);

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMAGE TYPE
// -----------------------------------------------------------------------------
namespace ygl {

// Image container.
template <typename T>
struct image {
    // constructors
    image() : _size{0, 0}, _data() {}
    image(const vec2i& wh, const T& v = T{})
        : _size{wh}, _data(wh.x * wh.y, v) {}
    image(const vec2i& wh, const T* v) : _size{wh}, _data(v, v + wh.x * wh.y) {}

    // size
    vec2i size() const { return _size; }
    bool empty() const { return _data.empty(); }

    // pixel access
    T& operator[](const vec2i& ij) { return _data[ij.y * _size.x + ij.x]; }
    const T& operator[](const vec2i& ij) const {
        return _data[ij.y * _size.x + ij.x];
    }
    // T& at(int i, int j) { return _data.at(j * _size.x + i); }
    // const T& at(int i, int j) const { return _data.at(j * _size.x + i); }

    // data acess
    T* data() { return _data.data(); }
    const T* data() const { return _data.data(); }
    T* begin() { return _data.data(); }
    T* end() { return _data.data() + _data.size(); }
    const T* begin() const { return _data.data(); }
    const T* end() const { return _data.data() + _data.size(); }
    const std::vector<T>& dataref() const { return _data; }

    // private data
   private:
    vec2i _size = {0, 0};
    std::vector<T> _data = {};
};

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
inline float luminance(const vec4f& a) { return (a.x + a.y + a.z) / 3; }

// Fitted ACES tonemapping curve.
inline float tonemap_filmic_curve(float hdr) {
    // https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
    // hdr *= 0.6; // brings it back to ACES range
    return (hdr * hdr * 2.51f + hdr * 0.03f) /
           (hdr * hdr * 2.43f + hdr * 0.59f + 0.14f);
}
// Apply ACES fitted curve.
inline vec4f tonemap_filmic_curve(const vec4f& hdr) {
    return {tonemap_filmic_curve(hdr.x), tonemap_filmic_curve(hdr.y),
        tonemap_filmic_curve(hdr.z), hdr.w};
}

// Tonemap a color value according to an exposure-gamma tone mapper, with
// an optional filmic curve.
inline vec4f tonemap_exposuregamma(
    const vec4f& hdr, float exposure, float gamma, bool filmic) {
    auto scale = pow(2.0f, exposure);
    auto ldr = vec4f{hdr.x * scale, hdr.y * scale, hdr.z * scale, hdr.w};
    if (filmic) ldr = tonemap_filmic_curve(ldr);
    return linear_to_gamma(ldr, gamma);
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
image<vec4f> gamma_to_linear(const image<vec4f>& srgb, float gamma = 2.2f);
image<vec4f> linear_to_gamma(const image<vec4f>& lin, float gamma = 2.2f);

// Apply exposure and filmic tone mapping
image<vec4f> tonemap_exposuregamma(
    const image<vec4f>& hdr, float exposure, float gamma, bool filmic);

// Resize an image.
image<vec4f> resize_image(const image<vec4f>& img, const vec2i& size);

}  // namespace ygl

// -----------------------------------------------------------------------------
// EXAMPLE IMAGES
// -----------------------------------------------------------------------------
namespace ygl {

// Make example images.
image<vec4f> make_grid_image4f(const vec2i& size, int tile = 8,
    const vec4f& c0 = {0.5f, 0.5f, 0.5f, 1},
    const vec4f& c1 = {0.8f, 0.8f, 0.8f, 1});
image<vec4f> make_checker_image4f(const vec2i& size, int tile = 8,
    const vec4f& c0 = {0.5f, 0.5f, 0.5f, 1},
    const vec4f& c1 = {0.8f, 0.8f, 0.8f, 1});
image<vec4f> make_bumpdimple_image4f(const vec2i& size, int tile = 8);
image<vec4f> make_ramp_image4f(
    const vec2i& size, const vec4f& c0, const vec4f& c1, float srgb = false);
image<vec4f> make_gammaramp_image4f(const vec2i& size);
image<vec4f> make_uvramp_image4f(const vec2i& size);
image<vec4f> make_uvgrid_image4f(
    const vec2i& size, int tile = 8, bool colored = true);

// Comvert a bump map to a normal map.
image<vec4f> bump_to_normal_map(const image<vec4f>& img, float scale = 1);

// Make a sunsky HDR model with sun at theta elevation in [0,pif/2], turbidity
// in [1.7,10] with or without sun.
image<vec4f> make_sunsky_image4f(const vec2i& size, float thetaSun,
    float turbidity = 3, bool has_sun = false,
    const vec3f& ground_albedo = {0.7f, 0.7f, 0.7f});
// Make an image of multiple lights.
image<vec4f> make_lights_image4f(const vec2i& size, const vec3f& le = {1, 1, 1},
    int nlights = 4, float langle = pif / 4, float lwidth = pif / 16,
    float lheight = pif / 16);

// Make a noise image. Wrap works only if both resx and resy are powers of two.
image<vec4f> make_noise_image4f(
    const vec2i& size, float scale = 1, bool wrap = true);
image<vec4f> make_fbm_image4f(const vec2i& size, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, int octaves = 6, bool wrap = true);
image<vec4f> make_ridge_image4f(const vec2i& size, float scale = 1,
    float lacunarity = 2, float gain = 0.5f, float offset = 1.0f,
    int octaves = 6, bool wrap = true);
image<vec4f> make_turbulence_image4f(const vec2i& size, float scale = 1,
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

// Evaluates a keyframed value using step interpolation.
template <typename T>
inline T eval_keyframed_step(
    const std::vector<float>& times, const std::vector<T>& vals, float time) {
    if (time <= times.front()) return vals.front();
    if (time >= times.back()) return vals.back();
    time = clamp(time, times.front(), times.back() - 0.001f);
    auto idx = eval_keyframed_index(times, time);
    return vals.at(idx - 1);
}

// Evaluates a keyframed value using linear interpolation.
template <typename T>
inline vec4f eval_keyframed_slerp(const std::vector<float>& times,
    const std::vector<vec4f>& vals, float time) {
    if (time <= times.front()) return vals.front();
    if (time >= times.back()) return vals.back();
    time = clamp(time, times.front(), times.back() - 0.001f);
    auto idx = eval_keyframed_index(times, time);
    auto t = (time - times.at(idx - 1)) / (times.at(idx) - times.at(idx - 1));
    return slerp(vals.at(idx - 1), vals.at(idx), t);
}

// Evaluates a keyframed value using linear interpolation.
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

// Evaluates a keyframed value using Bezier interpolation.
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
// VOLUME TYPE
// -----------------------------------------------------------------------------
namespace ygl {

// Volume container.
template <typename T>
struct volume {
    // constructors
    volume() : _size{0, 0}, _data() {}
    volume(const vec3i& size, const T& v = T{})
        : _size{size}, _data(size.x * size.y, v) {}
    volume(const vec3i& size, const T* v)
        : _size{size}, _data(v, v + size.x * size.y * size.z) {}

    // size
    vec3i size() const { return _size; }
    bool empty() const { return _data.empty(); }

    // pixel access
    T& operator[](const vec3i& ijk) {
        return _data[ijk.z * _size.x * _size.y + ijk.y * _size.x + ijk.x];
    }
    const T& operator[](const vec3i& ijk) const {
        return _data[ijk.z * _size.x * _size.y + ijk.y * _size.x + ijk.x];
    }
    // T& at(int i, int j) { return _data.at(ij.z * _size.x * _size.y + j *
    // _size.x + i); } const T& at(int i, int j) const { return _data.at(ij.z *
    // _size.x * _size.y + j * _size.x + i); }

    // data acess
    T* data() { return _data.data(); }
    const T* data() const { return _data.data(); }
    T* begin() { return _data.data(); }
    T* end() { return _data.data() + _data.size(); }
    const T* begin() const { return _data.data(); }
    const T* end() const { return _data.data() + _data.size(); }
    const std::vector<T>& dataref() const { return _data; }

    // private data
   private:
    vec3i _size = {0, 0};
    std::vector<T> _data = {};
};

}  // namespace ygl

// -----------------------------------------------------------------------------
// SCENE DATA
// -----------------------------------------------------------------------------
namespace ygl {

// make a simple example volume
volume<float> make_test_volume1f(
    const vec3i& size, float scale = 10, float exponent = 6);

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
    vec2f film = {0.036f, 0.024f};     // film size (default: 35mm)
    float focal = 0.050f;              // focal length (defaut: 50 mm)
    float focus = maxf;                // focal distance (default: infinite)
    float aperture = 0;                // lens aperture
};

// Texture containing either an LDR or HDR image.
struct texture {
    std::string name = "";     // name
    std::string path = "";     // file path
    image<vec4f> img = {};     // image
    bool clamp = false;        // clamp textures coordinates
    bool linear = true;        // use bilinear interpolation
    float scale = 1;           // scale for occ, normal, bumps
    float gamma = 2.2f;        // gamma correction for ldr textures in IO
    bool has_opacity = false;  // check whether alpha != 0
};

// Volumetric texture containing either an HDR image.
struct voltexture {
    std::string name = "";   // name
    std::string path = "";   // file path
    volume<float> vol = {};  // volume
    bool clamp = false;      // clamp textures coordinates
    bool linear = true;      // use trilinear interpolation
};

// Material for surfaces, lines and triangles.
// For surfaces, uses a microfacet model with thin sheet transmission.
// The model is based on glTF for compatibility and adapted to OBJ.
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
    texture* ke_txt = nullptr;    // emission texture
    texture* kd_txt = nullptr;    // diffuse texture
    texture* ks_txt = nullptr;    // specular texture
    texture* kt_txt = nullptr;    // transmission texture
    texture* rs_txt = nullptr;    // roughness texture
    texture* op_txt = nullptr;    // opacity texture
    texture* occ_txt = nullptr;   // occlusion texture
    texture* bump_txt = nullptr;  // bump map texture (heighfield)
    texture* disp_txt = nullptr;  // displacement map texture (heighfield)
    texture* norm_txt = nullptr;  // normal texture

    // volume properties
    vec3f ve = {0, 0, 0};  // volume emission
    vec3f va = {0, 0, 0};  // albedo: scattering / (absorption + scattering)
    vec3f vd = {0, 0, 0};  // density: absorption + scattering
    float vg = 0;          // phase function shape

    // volume textures
    voltexture* vd_txt = nullptr;  // density
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
    std::string name = "";             // name
    frame3f frame = identity_frame3f;  // transform frame
    shape* shp = nullptr;              // shape
    material* mat = nullptr;           // material
    subdiv* sbd = nullptr;             // subdivision shape
};

// Environment map.
struct environment {
    std::string name = "";             // name
    frame3f frame = identity_frame3f;  // transform frame
    vec3f ke = {0, 0, 0};              // emission color
    texture* ke_txt = nullptr;         // emission texture
};

// Node in a transform hierarchy.
struct node {
    std::string name = "";             // name
    node* parent = nullptr;            // parent
    frame3f local = identity_frame3f;  // transform frame
    vec3f translation = {0, 0, 0};     // translation
    vec4f rotation = {0, 0, 0, 1};     // rotation
    vec3f scale = {1, 1, 1};           // scale
    std::vector<float> weights = {};   // morph weights
    camera* cam = nullptr;             // camera
    instance* ist = nullptr;           // instance
    environment* env = nullptr;        // environment

    // compute properties
    std::vector<node*> children = {};  // child nodes
};

// Keyframe type.
enum struct animation_type { linear, step, bezier };

// Keyframe data.
struct animation {
    std::string name = "";                         // name
    std::string path = "";                         // path for glTF buffer
    std::string group = "";                        // group
    animation_type type = animation_type::linear;  // type
    std::vector<float> times;                      // keyframe times
    std::vector<vec3f> translation;                // translation keyframes
    std::vector<vec4f> rotation;                   // rotation keyframes
    std::vector<vec3f> scale;                      // scale keyframes
    std::vector<std::vector<float>> weights;       // morph weight keyframes
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
    std::string name;                             // name
    std::vector<camera*> cameras = {};            // cameras
    std::vector<shape*> shapes = {};              // shapes
    std::vector<subdiv*> subdivs = {};            // subdivs
    std::vector<instance*> instances = {};        // instances
    std::vector<material*> materials = {};        // materials
    std::vector<texture*> textures = {};          // textures
    std::vector<environment*> environments = {};  // environments
    std::vector<voltexture*> voltextures = {};    // volume textures

    std::vector<node*> nodes = {};            // node hierarchy [optional]
    std::vector<animation*> animations = {};  // animations [optional]

    // cleanup
    ~scene();
};

}  // namespace ygl

// -----------------------------------------------------------------------------
// SCENE UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Print scene statistics.
void print_stats(const scene* scn);

// Merge scene into one another. Note that the objects are _moved_ from
// merge_from to merged_into, so merge_from will be empty after this function.
void merge_into(const scene* merge_into, scene* merge_from);

}  // namespace ygl

// -----------------------------------------------------------------------------
// UPDATES TO COMPUTED PROPERTIES
// -----------------------------------------------------------------------------
namespace ygl {

// Update node transforms.
void update_transforms(
    scene* scn, float time = 0, const std::string& anim_group = "");
// Compute animation range.
vec2f compute_animation_range(
    const scene* scn, const std::string& anim_group = "");

// Computes shape/scene approximate bounds.
bbox3f compute_bbox(const shape* shp);
bbox3f compute_bbox(const scene* scn);

// Generate a distribution for sampling a shape uniformly based on area/length.
std::vector<float> compute_shape_cdf(const shape* shp);
// Generate a distribution for sampling an environment texture uniformly
// based on angle and texture intensity.
std::vector<float> compute_environment_cdf(const environment* env);

// Updates/refits bvh.
bvh_tree* build_bvh(const shape* shp, bool high_quality, bool embree = false);
bvh_tree* build_bvh(const scene* scn, bool high_quality, bool embree = false);
void refit_bvh(const shape* shp, bvh_tree* bvh);
void refit_bvh(const scene* scn, bvh_tree* bvh);

// Updates tesselation.
void tesselate_subdiv(const subdiv* sbd, shape* shp);
void tesselate_subdivs(scene* scn);

// Add missing names, normals, tangents and hierarchy.
void add_missing_names(scene* scn);
void add_missing_normals(scene* scn);
void add_missing_tangent_space(scene* scn);
void add_missing_materials(scene* scn);
// Checks for validity of the scene.
std::vector<std::string> validate(const scene* scn, bool skip_textures = false);

// make camera
camera* make_bbox_camera(const std::string& name, const bbox3f& bbox,
    const vec2f& film = {0.036f, 0.024f}, float focal = 0.050f);
// make default material
inline material* make_default_material(const std::string& name) {
    auto mat = new material();
    mat->name = name;
    mat->kd = {0.2f, 0.2f, 0.2f};
    return mat;
}

// Add a sky environment
inline environment* make_sky_environment(
    const std::string& name, float sun_angle = pif / 4) {
    auto txt = new texture();
    txt->name = name;
    txt->path = "textures/" + name + ".hdr";
    txt->img = make_sunsky_image4f({1024, 512}, sun_angle);
    auto env = new environment();
    env->name = name;
    env->ke = {1, 1, 1};
    env->ke_txt = txt;
    return env;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// EVALUATION OF SCENE PROPERTIES
// -----------------------------------------------------------------------------
namespace ygl {

// Scene intersection.
struct scene_intersection {
    instance* ist = nullptr;  // instance or null for no intersection
    int ei = 0;               // shape element index
    vec2f uv = zero2f;        // shape element coordinates
    float dist = maxf;        // ray/point distance
};

// Intersects a ray with an instance. The bvh refers is the shape bvh.
scene_intersection intersect_ray(const instance* ist, const bvh_tree* sbvh,
    const ray3f& ray, bool find_any = false);
// Intersects a ray with the scene.
scene_intersection intersect_ray(const scene* scn, const bvh_tree* bvh,
    const ray3f& ray, bool find_any = false);

// Shape values interpolated using barycentric coordinates.
vec3f eval_pos(const shape* shp, int ei, const vec2f& uv);
vec3f eval_norm(const shape* shp, int ei, const vec2f& uv);
vec2f eval_texcoord(const shape* shp, int ei, const vec2f& uv);
vec4f eval_color(const shape* shp, int ei, const vec2f& uv);
float eval_radius(const shape* shp, int ei, const vec2f& uv);
vec4f eval_tangsp(const shape* shp, int ei, const vec2f& uv);
vec3f eval_tangsp(const shape* shp, int ei, const vec2f& uv, bool& left_handed);
// Shape element values.
vec3f eval_elem_norm(const shape* shp, int ei);
vec4f eval_elem_tangsp(const shape* shp, int ei);

// Instance values interpolated using barycentric coordinates.
// Handles defaults if data is missing.
vec3f eval_pos(const instance* ist, int ei, const vec2f& uv);
vec3f eval_norm(const instance* ist, int ei, const vec2f& uv);
vec2f eval_texcoord(const instance* ist, int ei, const vec2f& uv);
vec4f eval_color(const instance* ist, int ei, const vec2f& uv);
float eval_radius(const instance* ist, int ei, const vec2f& uv);
vec3f eval_tangsp(
    const instance* ist, int ei, const vec2f& uv, bool& left_handed);
// Instance element values.
vec3f eval_elem_norm(const instance* ist, int ei);
// Shading normals including material perturbations.
vec3f eval_shading_norm(
    const instance* ist, int ei, const vec2f& uv, const vec3f& o);

// Environment texture coordinates from the incoming direction.
vec2f eval_texcoord(const environment* env, const vec3f& i);
// Evaluate the incoming direction from the uv.
vec3f eval_direction(const environment* env, const vec2f& uv);
// Evaluate the environment emission.
vec3f eval_environment(const environment* env, const vec3f& i);
// Evaluate all environment emission.
vec3f eval_environment(const scene* scn, const vec3f& i);

// Evaluate a texture.
vec4f eval_texture(const texture* txt, const vec2f& texcoord);
float eval_voltexture(const voltexture* txt, const vec3f& texcoord);

// Set and evaluate camera parameters. Setters take zeros as default values.
float eval_camera_fovy(const camera* cam);
float eval_camera_aspect(const camera* cam);
void set_camera_fovy(
    camera* cam, float fovy, float aspect, float width = 0.036f);
vec2i eval_image_size(const camera* cam, int yresolution);

// Generates a ray from a camera image coordinate `uv` and lens coordinates
// `luv`.
ray3f eval_camera_ray(const camera* cam, const vec2f& uv, const vec2f& luv);
// Generates a ray from a camera for pixel coordinates `ij`, the image size
// `imsize`, the sub-pixel coordinates `puv` and the lens coordinates `luv` and
// the image resolution `res`.
ray3f eval_camera_ray(const camera* cam, const vec2i& ij, const vec2i& imsize,
    const vec2f& puv, const vec2f& luv);
// Generates a ray from a camera for pixel index `idx`, the image size
// `imsize`, the sub-pixel coordinates `puv` and the lens coordinates `luv` and
// the image resolution `res`.
ray3f eval_camera_ray(const camera* cam, int idx, const vec2i& imsize,
    const vec2f& puv, const vec2f& luv);

// Evaluates material parameters: emission, diffuse, specular, transmission,
// roughness and opacity.
vec3f eval_emission(const instance* ist, int ei, const vec2f& uv);
vec3f eval_diffuse(const instance* ist, int ei, const vec2f& uv);
vec3f eval_specular(const instance* ist, int ei, const vec2f& uv);
vec3f eval_transmission(const instance* ist, int ei, const vec2f& uv);
float eval_roughness(const instance* ist, int ei, const vec2f& uv);
float eval_opacity(const instance* ist, int ei, const vec2f& uv);

// Material values packed into a convenience structure.
struct bsdf {
    vec3f kd = zero3f;     // diffuse
    vec3f ks = zero3f;     // specular
    vec3f kt = zero3f;     // transmission
    float rs = 1;          // roughness
    bool refract = false;  // whether to use refraction in transmission
};
bsdf eval_bsdf(const instance* ist, int ei, const vec2f& uv);
bool is_delta_bsdf(const bsdf& f);

// Check volume properties.
bool is_volume_homogeneus(const material* vol);
bool has_volume_color(const material* vol);

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

const auto trace_type_names = std::vector<std::string>{"path", "volpath",
    "direct", "environment", "eyelight", "path_nomis", "path_naive",
    "direct_nomis", "debug_normal", "debug_albedo", "debug_texcoord",
    "debug_frontfacing", "debug_diffuse", "debug_specular", "debug_roughness"};

// Trace options
struct trace_params {
    int camid = 0;                         // camera index
    int yresolution = 256;                 // vertical resolution
    trace_type tracer = trace_type::path;  // tracer type
    int nsamples = 256;                    // number of samples
    int nbounces = 8;                      // max number of bounces
    float pixel_clamp = 100;               // pixel clamping
    int nbatch = 16;                       // number of samples per batch
    bool noparallel = false;               // serial or parallel execution
    int preview_ratio = 8;                 // preview ratio for asycn rendering
    float exposure = 0;                    // tone mapping exposure
    float gamma = 2.2f;                    // tone mapping gamma
    bool filmic = false;                   // tone mapping filmic
    int seed = trace_default_seed;         // trace seed
};

// Trace lights used during rendering.
struct trace_lights {
    std::vector<instance*> lights;           // instance lights
    std::vector<environment*> environments;  // environments lights
    std::unordered_map<shape*, std::vector<float>> shape_cdf;      // shape cdfs
    std::unordered_map<environment*, std::vector<float>> env_cdf;  // env cdfs
};

// Trace data used during rendering. Initialize with `make_trace_state()`
struct trace_state {
    image<vec4f> img = {};      // image being rendered
    image<vec4f> display = {};  // image tone mapped for display

    // internal data used during rendering
    image<vec4f> acc = {};             // accumulation buffer
    image<int> samples = {};           // samples per pixel
    image<rng_state> rng = {};         // random number generators
    int sample = 0;                    // current sample being rendered
    std::vector<std::thread> threads;  // threads used during rendering
    bool stop = false;                 // stop flag for threads
};

// Initialize lights.
trace_lights* make_trace_lights(const scene* scn, const trace_params& params);

// Initialize state of the renderer.
trace_state* make_trace_state(const scene* scn, const trace_params& params);

// Progressively compute an image by calling trace_samples multiple times.
image<vec4f> trace_image4f(const scene* scn, const bvh_tree* bvh,
    const trace_lights* lights, const trace_params& params);

// Progressively compute an image by calling trace_samples multiple times.
// Start with an empty state and then successively call this function to
// render the next batch of samples.
bool trace_samples(trace_state* state, const scene* scn, const bvh_tree* bvh,
    const trace_lights* lights, const trace_params& params);

// Starts an anyncrhounous renderer. The function will keep a reference to
// params.
void trace_async_start(trace_state* state, const scene* scn,
    const bvh_tree* bvh, const trace_lights* lights,
    const trace_params& params);
// Stop the asynchronous renderer.
void trace_async_stop(trace_state* state);

// Trace statistics for last run used for fine tuning implementation.
// For now returns number of paths and number of rays.
std::pair<uint64_t, uint64_t> get_trace_stats();
void reset_trace_stats();

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
inline mat2<T> transpose(const mat2<T>& a) {
    return {{a.x.x, a.y.x}, {a.x.y, a.y.y}};
}
template <typename T>
inline mat3<T> transpose(const mat3<T>& a) {
    return {
        {a.x.x, a.y.x, a.z.x},
        {a.x.y, a.y.y, a.z.y},
        {a.x.z, a.y.z, a.z.z},
    };
}
template <typename T>
inline mat4<T> transpose(const mat4<T>& a) {
    return {
        {a.x.x, a.y.x, a.z.x, a.w.x},
        {a.x.y, a.y.y, a.z.y, a.w.y},
        {a.x.z, a.y.z, a.z.z, a.w.z},
        {a.x.w, a.y.w, a.z.w, a.w.w},
    };
}

// Matrix adjugates, determinant and inverses.
template <typename T>
inline mat2<T> adjugate(const mat2<T>& a) {
    return {{a.y.y, -a.x.y}, {-a.y.x, a.x.x}};
}
template <typename T>
inline mat3<T> adjugate(const mat3<T>& a) {
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
inline mat4<T> adjugate(const mat4<T>& a) {
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
inline T determinant(const mat2<T>& a) {
    return a.x.x * a.y.y - a.x.y * a.y.x;
}
template <typename T>
inline T determinant(const mat3<T>& a) {
    return a.x.x * (a.y.y * a.z.z - a.z.y * a.y.z) +
           a.x.y * (a.y.z * a.z.x - a.z.z * a.y.x) +
           a.x.z * (a.y.x * a.z.y - a.z.x * a.y.y);
}
template <typename T>
inline T determinant(const mat4<T>& a) {
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

// Turntable for UI navigation.
inline void camera_turntable(frame3f& frame, float& focus, const vec2f& rotate,
    float dolly, const vec2f& pan) {
    // rotate if necessary
    if (rotate != zero2f) {
        auto phi = atan2(frame.z.z, frame.z.x) + rotate.x;
        auto theta = acos(frame.z.y) + rotate.y;
        theta = clamp(theta, 0.001f, pif - 0.001f);
        auto new_z =
            vec3f{sin(theta) * cos(phi), cos(theta), sin(theta) * sin(phi)};
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
