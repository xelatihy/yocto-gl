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
// in graphics. In particular, we support 2-4 dimensional vectors `vec2f`,
// `vec3f`, `vec4f`.
//
// We support 2-4 dimensional generic matrices `mat2f`, `mat3f`,
// `mat4f`, with matrix-matrix and matrix-vector products, transposes and
// inverses. Matrices are stored in column-major ordered and are accessed and
// constructed by column.
//
// To represent transformations, most of the library facilities prefer the use
// coordinate frames, aka rigid transforms, represented as `frame2f`,
// `frame3f`. The structure store three coordinate axes and the origin.
// This is equivalent to a rigid transform written as a column-major affine
// matrix. Transform operations are better behaved with this representation.
//
// We represent coordinate bounds with axis-aligned bounding boxes with
// `bbox1f`, `bbox2f`, `bbox3f`, `bbox4f`, with support for
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
// w.r.t the (v1-v0) and (v2-v0) axis respectively. Quads are internally handled
// as pairs of two triangles v0,v1,v3 and v2,v3,v1, with the u/v coordinates
// of the second triangle corrected as 1-u and 1-v to produce a quad
// parametrization where u and v go from 0 to 1. Degenerate quads with v2==v3
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
struct vec2f {
    float x = 0;
    float y = 0;
};
struct vec3f {
    float x = 0;
    float y = 0;
    float z = 0;
};
struct vec4f {
    float x = 0;
    float y = 0;
    float z = 0;
    float w = 0;
};

// Zero vector constants.
const auto zero2f = vec2f{0, 0};
const auto zero3f = vec3f{0, 0, 0};
const auto zero4f = vec4f{0, 0, 0, 0};

// Access xyz component of a vec4 typically used for color operation.
inline vec3f& xyz(const vec4f& a) { return (vec3f&)a; }
inline vec3f& xyz(vec4f& a) { return (vec3f&)a; }

// Vector comparison operations.
inline bool operator==(const vec2f& a, const vec2f& b) {
    return a.x == b.x && a.y == b.y;
}
inline bool operator!=(const vec2f& a, const vec2f& b) {
    return a.x != b.x || a.y != b.y;
}
inline bool operator==(const vec3f& a, const vec3f& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}
inline bool operator!=(const vec3f& a, const vec3f& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z;
}
inline bool operator==(const vec4f& a, const vec4f& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}
inline bool operator!=(const vec4f& a, const vec4f& b) {
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

// Vector assignments
inline vec2f& operator+=(vec2f& a, const vec2f& b) { return a = a + b; }
inline vec2f& operator-=(vec2f& a, const vec2f& b) { return a = a - b; }
inline vec2f& operator*=(vec2f& a, const vec2f& b) { return a = a * b; }
inline vec2f& operator*=(vec2f& a, float b) { return a = a * b; }
inline vec2f& operator/=(vec2f& a, const vec2f& b) { return a = a / b; }
inline vec2f& operator/=(vec2f& a, float b) { return a = a / b; }

// Vector assignments
inline vec3f& operator+=(vec3f& a, const vec3f& b) { return a = a + b; }
inline vec3f& operator-=(vec3f& a, const vec3f& b) { return a = a - b; }
inline vec3f& operator*=(vec3f& a, const vec3f& b) { return a = a * b; }
inline vec3f& operator*=(vec3f& a, float b) { return a = a * b; }
inline vec3f& operator/=(vec3f& a, const vec3f& b) { return a = a / b; }
inline vec3f& operator/=(vec3f& a, float b) { return a = a / b; }

// Vector assignments
inline vec4f& operator+=(vec4f& a, const vec4f& b) { return a = a + b; }
inline vec4f& operator-=(vec4f& a, const vec4f& b) { return a = a - b; }
inline vec4f& operator*=(vec4f& a, const vec4f& b) { return a = a * b; }
inline vec4f& operator*=(vec4f& a, float b) { return a = a * b; }
inline vec4f& operator/=(vec4f& a, const vec4f& b) { return a = a / b; }
inline vec4f& operator/=(vec4f& a, float b) { return a = a / b; }

// Vector products and lengths.
// Vector products and lengths.
inline float dot(const vec2f& a, const vec2f& b) {
    return a.x * b.x + a.y * b.y;
}
inline float length(const vec2f& a) { return sqrt(a.x * a.x + a.y * a.y); }
inline vec2f normalize(const vec2f& a) {
    auto l = length(a);
    return (l) ? a / l : a;
}
inline float cross(const vec2f& a, const vec2f& b) {
    return a.x * b.y - a.y * b.x;
}
inline float dot(const vec3f& a, const vec3f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline float length(const vec3f& a) {
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}
inline vec3f normalize(const vec3f& a) {
    auto l = length(a);
    return (l) ? a / l : a;
}
inline vec3f cross(const vec3f& a, const vec3f& b) {
    return {
        a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline float dot(const vec4f& a, const vec4f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}
inline float length(const vec4f& a) {
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z + a.w * a.w);
}
inline vec4f normalize(const vec4f& a) {
    auto l = length(a);
    return (l) ? a / l : a;
}

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
    if (k < 0) return {0, 0, 0};  // tir
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
    return quat_conjugate(a) / dot(a, a);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// INT VECTORS
// -----------------------------------------------------------------------------
namespace ygl {

// Small size vectors.
struct vec2i {
    int x = 0;
    int y = 0;
};
struct vec3i {
    int x = 0;
    int y = 0;
    int z = 0;
};
struct vec4i {
    int x = 0;
    int y = 0;
    int z = 0;
    int w = 0;
};
struct vec4b {
    byte x = 0;
    byte y = 0;
    byte z = 0;
    byte w = 0;
};

// Zero vector constants.
const auto zero2i = vec2i{0, 0};
const auto zero3i = vec3i{0, 0, 0};
const auto zero4i = vec4i{0, 0, 0, 0};
const auto zero4b = vec4b{0, 0, 0, 0};

// Vector comparison operations.
inline bool operator==(const vec2i& a, const vec2i& b) {
    return a.x == b.x && a.y == b.y;
}
inline bool operator!=(const vec2i& a, const vec2i& b) {
    return a.x != b.x || a.y != b.y;
}
inline bool operator==(const vec3i& a, const vec3i& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}
inline bool operator!=(const vec3i& a, const vec3i& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z;
}
inline bool operator==(const vec4i& a, const vec4i& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}
inline bool operator!=(const vec4i& a, const vec4i& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w;
}

}  // namespace ygl

namespace std {

// Hash functor for vector for use with unordered_map
template <>
struct hash<ygl::vec2i> {
    size_t operator()(const ygl::vec2i& v) const {
        auto vh = hash<int>();
        auto h = (size_t)0;
        for (auto i = 0; i < 2; i++)
            h ^= vh((&v.x)[i]) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};
// Hash functor for vector for use with unordered_map
template <>
struct hash<ygl::vec3i> {
    size_t operator()(const ygl::vec3i& v) const {
        auto vh = hash<int>();
        auto h = (size_t)0;
        for (auto i = 0; i < 3; i++)
            h ^= vh((&v.x)[i]) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};
// Hash functor for vector for use with unordered_map
template <>
struct hash<ygl::vec4i> {
    size_t operator()(const ygl::vec4i& v) const {
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
struct mat2f {
    vec2f x = {1, 0};
    vec2f y = {0, 1};
};
struct mat3f {
    vec3f x = {1, 0, 0};
    vec3f y = {0, 1, 0};
    vec3f z = {0, 0, 1};
};
struct mat4f {
    vec4f x = {1, 0, 0, 0};
    vec4f y = {0, 1, 0, 0};
    vec4f z = {0, 0, 1, 0};
    vec4f w = {0, 0, 0, 1};
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

// Matrix assignments.
inline mat3f& operator+=(mat3f& a, const mat3f& b) { return a = a + b; }
inline mat3f& operator*=(mat3f& a, const mat3f& b) { return a = a * b; }
inline mat3f& operator*=(mat3f& a, float b) { return a = a * b; }

// Matrix assignments.
inline mat4f& operator+=(mat4f& a, const mat4f& b) { return a = a + b; }
inline mat4f& operator*=(mat4f& a, const mat4f& b) { return a = a * b; }
inline mat4f& operator*=(mat4f& a, float b) { return a = a * b; }

// Matrix diagonals and transposes.
inline vec2f diagonal(const mat2f& a) { return {a.x.x, a.y.y}; }
inline vec3f diagonal(const mat3f& a) { return {a.x.x, a.y.y, a.z.z}; }
inline vec4f diagonal(const mat4f& a) { return {a.x.x, a.y.y, a.z.z, a.w.w}; }
inline mat2f transpose(const mat2f& a);
inline mat3f transpose(const mat3f& a);
inline mat4f transpose(const mat4f& a);

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
    vec2f x = {1, 0};
    vec2f y = {0, 1};
    vec2f o = {0, 0};
};
struct frame3f {
    vec3f x = {1, 0, 0};
    vec3f y = {0, 1, 0};
    vec3f z = {0, 0, 1};
    vec3f o = {0, 0, 0};
};

// Indentity frames.
const auto identity_frame2f = frame2f{{1, 0}, {0, 1}, {0, 0}};
const auto identity_frame3f =
    frame3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};

// Frame construction from axis.
inline frame3f make_frame_fromz(const vec3f& o, const vec3f& v) {
    auto z = normalize(v);
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
    return {
        {a.x.x, a.x.y, a.x.z, 0},
        {a.y.x, a.y.y, a.y.z, 0},
        {a.z.x, a.z.y, a.z.z, 0},
        {a.o.x, a.o.y, a.o.z, 1},
    };
}
inline frame3f mat_to_frame(const mat4f& a) {
    return {
        {a.x.x, a.x.y, a.x.z},
        {a.y.x, a.y.y, a.y.z},
        {a.z.x, a.z.y, a.z.z},
        {a.w.x, a.w.y, a.w.z},
    };
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
inline frame2f operator*(const frame2f& a, const frame2f& b) {
    auto rot = mat2f{a.x, a.y} * mat2f{b.x, b.y};
    auto pos = mat2f{a.x, a.y} * b.o + a.o;
    return {rot.x, rot.y, pos};
}
inline frame3f operator*(const frame3f& a, const frame3f& b) {
    auto rot = mat3f{a.x, a.y, a.z} * mat3f{b.x, b.y, b.z};
    auto pos = mat3f{a.x, a.y, a.z} * b.o + a.o;
    return {rot.x, rot.y, rot.z, pos};
}
// Frame inverse, equivalent to rigid affine inverse.
inline frame2f inverse(const frame2f& a, bool is_rigid = true) {
    auto minv =
        (is_rigid) ? transpose(mat2f{a.x, a.y}) : inverse(mat2f{a.x, a.y});
    return {minv.x, minv.y, -(minv * a.o)};
}
inline frame3f inverse(const frame3f& a, bool is_rigid = true) {
    auto minv = (is_rigid) ? transpose(mat3f{a.x, a.y, a.z}) :
                             inverse(mat3f{a.x, a.y, a.z});
    return {minv.x, minv.y, minv.z, -(minv * a.o)};
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// AXIS ALIGNED BOUNDING BOXES
// -----------------------------------------------------------------------------
namespace ygl {

// Range of values in 1D.
struct bbox1f {
    float min = flt_max;
    float max = flt_min;
};

// Axis aligned bounding box represented as a min/max vector pairs.
struct bbox2f {
    vec2f min = {flt_max, flt_max};
    vec2f max = {flt_min, flt_min};
};

// Axis aligned bounding box represented as a min/max vector pairs.
struct bbox3f {
    vec3f min = {flt_max, flt_max, flt_max};
    vec3f max = {flt_min, flt_min, flt_min};
};

// Axis aligned bounding box represented as a min/max vector pairs.
struct bbox4f {
    vec4f min = {flt_max, flt_max, flt_max, flt_max};
    vec4f max = {flt_min, flt_min, flt_min, flt_min};
};

// Empty bbox constant.
const auto invalid_bbox1f = bbox1f();
const auto invalid_bbox2f = bbox2f();
const auto invalid_bbox3f = bbox3f();
const auto invalid_bbox4f = bbox4f();

// Bounding box comparisons.
inline bool operator==(const bbox1f& a, const bbox1f& b) {
    return a.min == b.min && a.max == b.max;
}
inline bool operator!=(const bbox1f& a, const bbox1f& b) {
    return a.min != b.min || a.max != b.max;
}
inline bool operator==(const bbox2f& a, const bbox2f& b) {
    return a.min == b.min && a.max == b.max;
}
inline bool operator!=(const bbox2f& a, const bbox2f& b) {
    return a.min != b.min || a.max != b.max;
}
inline bool operator==(const bbox3f& a, const bbox3f& b) {
    return a.min == b.min && a.max == b.max;
}
inline bool operator!=(const bbox3f& a, const bbox3f& b) {
    return a.min != b.min || a.max != b.max;
}
inline bool operator==(const bbox4f& a, const bbox4f& b) {
    return a.min == b.min && a.max == b.max;
}
inline bool operator!=(const bbox4f& a, const bbox4f& b) {
    return a.min != b.min || a.max != b.max;
}

// Bounding box expansions with points and other boxes.
inline bbox1f& operator+=(bbox1f& a, float b) {
    a.min = min(a.min, b);
    a.max = max(a.max, b);
    return a;
}
inline bbox1f& operator+=(bbox1f& a, const bbox1f& b) {
    a.min = min(a.min, b.min);
    a.max = max(a.max, b.max);
    return a;
}
// Bounding box expansions with points and other boxes.
inline bbox2f& operator+=(bbox2f& a, const vec2f& b) {
    a.min = {min(a.min.x, b.x), min(a.min.y, b.y)};
    a.max = {max(a.max.x, b.x), max(a.max.y, b.y)};
    return a;
}
inline bbox2f& operator+=(bbox2f& a, const bbox2f& b) {
    a.min = {min(a.min.x, b.min.x), min(a.min.y, b.min.y)};
    a.max = {max(a.max.x, b.max.x), max(a.max.y, b.max.y)};
    return a;
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
// Bounding box expansions with points and other boxes.
inline bbox4f& operator+=(bbox4f& a, const vec4f& b) {
    a.min = {min(a.min.x, b.x), min(a.min.y, b.y), min(a.min.z, b.z),
        min(a.min.w, b.w)};
    a.max = {max(a.max.x, b.x), max(a.max.y, b.y), max(a.max.z, b.z),
        max(a.max.w, b.w)};
    return a;
}
inline bbox4f& operator+=(bbox4f& a, const bbox4f& b) {
    a.min = {min(a.min.x, b.min.x), min(a.min.y, b.min.y),
        min(a.min.z, b.min.z), min(a.min.w, b.min.w)};
    a.max = {max(a.max.x, b.max.x), max(a.max.y, b.max.y),
        max(a.max.z, b.max.z), max(a.max.w, b.max.w)};
    return a;
}

// Primitive bounds.
inline bbox3f point_bbox(const vec3f& p, float r = 0) {
    auto bbox = ygl::bbox3f{};
    bbox += p - vec3f{r, r, r};
    bbox += p + vec3f{r, r, r};
    return bbox;
}
inline bbox3f line_bbox(
    const vec3f& v0, const vec3f& v1, float r0 = 0, float r1 = 0) {
    auto bbox = ygl::bbox3f{};
    bbox += v0 - vec3f{r0, r0, r0};
    bbox += v0 + vec3f{r0, r0, r0};
    bbox += v1 - vec3f{r1, r1, r1};
    bbox += v1 + vec3f{r1, r1, r1};
    return bbox;
}
inline bbox3f triangle_bbox(const vec3f& v0, const vec3f& v1, const vec3f& v2) {
    auto bbox = ygl::bbox3f{};
    bbox += v0;
    bbox += v1;
    bbox += v2;
    return bbox;
}
inline bbox3f quad_bbox(
    const vec3f& v0, const vec3f& v1, const vec3f& v2, const vec3f& v3) {
    auto bbox = ygl::bbox3f{};
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
struct ray2f {
    vec2f o = {0, 0};
    vec2f d = {0, 1};
    float tmin = 0;
    float tmax = flt_max;
};

// Rays with origin, direction and min/max t value.
struct ray3f {
    vec3f o = {0, 0, 0};
    vec3f d = {0, 0, 1};
    float tmin = 0;
    float tmax = flt_max;
};

// Construct a ray from dirction or segments using a default epsilon.
inline ray3f make_ray(const vec3f& o, const vec3f& d, float eps = 1e-4f) {
    return {o, d, eps, flt_max};
}
inline ray3f make_segment(const vec3f& p1, const vec3f& p2, float eps = 1e-4f) {
    return {p1, normalize(p2 - p1), eps, length(p2 - p1) - 2 * eps};
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
inline vec4f rotation_quat(const vec4f& axisangle) {
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
inline void center_image(vec2f& center, float& scale, const vec2i& imsize,
    vec2i winsize, bool zoom_to_fit) {
    if (zoom_to_fit) {
        scale =
            ygl::min(winsize.x / (float)imsize.x, winsize.y / (float)imsize.y);
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

// Interpolates values over a line parametrized from a to b by u. Same as lerp.
template <typename T>
inline T interpolate_line(const T& v0, const T& v1, float u) {
    return v0 * (1 - u) + v1 * u;
}
// Interpolates values over a triangle parametrized by u and v along the
// (v1-v0) and (v2-v0) directions. Same as barycentric interpolation.
template <typename T>
inline T interpolate_triangle(
    const T& v0, const T& v1, const T& v2, const vec2f& uv) {
    return v0 * (1 - uv.x - uv.y) + v1 * uv.x + v2 * uv.y;
}
// Interpolates values over a quad parametrized by u and v along the
// (v1-v0) and (v2-v1) directions. Same as bilear interpolation.
template <typename T>
inline T interpolate_quad(
    const T& v0, const T& v1, const T& v2, const T& v3, const vec2f& uv) {
    return v0 * (1 - uv.x) * (1 - uv.y) + v1 * uv.x * (1 - uv.y) +
           v2 * uv.x * uv.y + v3 * (1 - uv.x) * uv.y;
}

// Interpolates values along a cubic Bezier segment parametrized by u.
template <typename T>
inline T interpolate_bezier(
    const T& v0, const T& v1, const T& v2, const T& v3, float u) {
    return v0 * (1 - u) * (1 - u) * (1 - u) + v1 * 3 * u * (1 - u) * (1 - u) +
           v2 * 3 * u * u * (1 - u) + v3 * u * u * u;
}
// Computes the derivative of a cubic Bezier segment parametrized by u.
template <typename T>
inline T interpolate_bezier_derivative(
    const T& v0, const T& v1, const T& v2, const T& v3, float u) {
    return (v1 - v0) * 3 * (1 - u) * (1 - u) + (v2 - v1) * 6 * u * (1 - u) +
           (v3 - v2) * 3 * u * u;
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
    const ray3f& ray, const vec3f& p, float r, float& dist, vec2f& uv);

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
bool intersect_bbox(const ray3f& ray, const vec3f& ray_dinv,
    const vec3i& ray_dsign, const bbox3f& bbox);

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
    quad,
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
    std::vector<vec4i> quads;      // Quads for shape BVHs.

    // data for instance BVH
    std::vector<frame3f> ist_frames;      // instance frames
    std::vector<frame3f> ist_inv_frames;  // instance inverse frames
    std::vector<bvh_tree*> ist_bvhs;      // instance shape bvhs

    // bvh nodes
    std::vector<bvh_node> nodes;  // Internal nodes.
};

// Build a BVH from the given set of primitives.
void build_bvh(bvh_tree* bvh, bool sah);
// Update the node bounds for a shape bvh.
void refit_bvh(bvh_tree* bvh);

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
    int width = 0;
    int height = 0;
    std::vector<T> pxl = {};

    // pixel access
    T& at(int i, int j) { return pxl.at(j * width + i); }
    const T& at(int i, int j) const { return pxl.at(j * width + i); }
};

// Type aliases
using image4f = image<vec4f>;
using image4b = image<vec4b>;

// Image creation.
template <typename T>
inline image<T> make_image(int width, int height, T c = T{}) {
    auto img = image<T>{};
    img.width = width;
    img.height = height;
    img.pxl.resize(width * height, c);
    return img;
}
inline image4f make_image4f(
    int width, int height, const vec4f& c = {0, 0, 0, 0}) {
    return make_image<vec4f>(width, height, c);
}
inline image4b make_image4b(
    int width, int height, const vec4b& c = {0, 0, 0, 0}) {
    return make_image<vec4b>(width, height, c);
}

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
image4f tonemap_image(
    const image4f& hdr, float exposure, float gamma, bool filmic);

// Resize an image.
image4f resize_image(const image4f& img, int width, int height);

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

// Fitted ACES tonemapping
inline vec3f tonemap_filmic(const vec3f& hdr) {
    // https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
    // hdr *= 0.6; // brings it back to ACES range
    return (hdr * hdr * 2.51f + hdr * 0.03f) /
           (hdr * hdr * 2.43f + hdr * 0.59f + vec3f{0.14f, 0.14f, 0.14f});
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
image4f make_uvramp_image(int width, int height);
image4f make_uvgrid_image(
    int width, int height, int tile = 8, bool colored = true);

// Comvert a bump map to a normal map.
image4f bump_to_normal_map(const image4f& img, float scale = 1);

// Make a sunsky HDR model with sun at theta elevation in [0,pi/2], turbidity
// in [1.7,10] with or without sun.
image4f make_sunsky_image(int width, int height, float thetaSun,
    float turbidity = 3, bool has_sun = false,
    const vec3f& ground_albedo = {0.7f, 0.7f, 0.7f});
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
// VOLUME TYPE
// -----------------------------------------------------------------------------
namespace ygl {

// Volume container.
template <typename T>
struct volume {
    int width = 0;
    int height = 0;
    int depth = 0;
    std::vector<T> pxl = {};

    // pixel access
    T& at(int i, int j, int k) {
        return pxl.at(k * height * width + j * width + i);
    }
    const T& at(int i, int j, int k) const {
        return pxl.at(k * height * width + j * width + i);
    }
};

// Type aliases
using volume4f = volume<vec4f>;
using volume1f = volume<float>;

// Image creation.
template <typename T>
inline volume<T> make_volume(int width, int height, int depth, T c = T{}) {
    auto img = volume<T>{};
    img.width = width;
    img.height = height;
    img.depth = depth;
    img.pxl.resize(width * height * depth, c);
    return img;
}
inline volume4f make_volume4f(
    int width, int height, int depth, const vec4f& c = {0, 0, 0, 0}) {
    return make_volume<vec4f>(width, height, depth, c);
}
inline volume1f make_volume1f(int width, int height, int depth, float c = 0) {
    return make_volume<float>(width, height, depth, c);
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
    float width = 0.036f;              // film width (default: 35mm)
    float height = 0.024f;             // film height (default: 35mm)
    float focal = 0.050f;              // focal length (defaut: 50 mm)
    float focus = flt_max;             // focal distance (default: infinite)
    float aperture = 0;                // lens aperture
    float near = 0.01f;                // near plane distance
    float far = 10000;                 // far plane distance
};

// Texture containing either an LDR or HDR image.
struct texture {
    std::string name = "";     // name
    std::string path = "";     // file path
    image4f img = {};          // image
    volume1f vol = {};         // volume
    bool clamp = false;        // clamp textures coordinates
    float scale = 1;           // scale for occ, normal, bumps
    float gamma = 2.2f;        // gamma correction for ldr textures in IO
    bool has_opacity = false;  // check whether alpha != 0
    uint gl_txt = 0;           // unmanaged data for OpenGL viewer
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
    vec3f ve = zero3f;  // volume emission
    vec3f va = zero3f;  // albedo: scattering / (absorption + scattering)
    vec3f vd = zero3f;  // density: absorption + scattering
    float vg = 0;       // phase function shape

    // volume textures
    texture* vd_txt = nullptr;  // density
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
    std::vector<float> elem_cdf = {};  // element cdf for sampling
    bvh_tree* bvh = nullptr;           // bvh for ray intersection
    uint gl_pos = 0, gl_norm = 0, gl_texcoord = 0, gl_color = 0, gl_tangsp = 0,
         gl_points = 0, gl_lines = 0,
         gl_triangles = 0;       // unmanaged data for OpenGL viewer
    void* embree_bvh = nullptr;  // unmanaged data for Embree raytracer

    // cleanup
    ~shape();
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
    std::string name;                  // name
    frame3f frame = identity_frame3f;  // transform frame
    shape* shp = nullptr;              // shape
    material* mat = nullptr;           // material
    subdiv* sbd = nullptr;             // subdivision shape
};

// Envinonment map.
struct environment {
    std::string name = "";             // name
    frame3f frame = identity_frame3f;  // transform frame
    vec3f ke = {0, 0, 0};              // emission color
    texture* ke_txt = nullptr;         // emission texture

    // computed properties
    std::vector<float> elem_cdf;  // element cdf for sampling
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
    std::string name;                             // name
    std::vector<camera*> cameras = {};            // cameras
    std::vector<shape*> shapes = {};              // shapes
    std::vector<subdiv*> subdivs = {};            // subdivs
    std::vector<instance*> instances = {};        // instances
    std::vector<material*> materials = {};        // materials
    std::vector<texture*> textures = {};          // textures
    std::vector<environment*> environments = {};  // environments

    std::vector<node*> nodes = {};            // node hierarchy [optional]
    std::vector<animation*> animations = {};  // animations [optional]

    // compute properties
    std::vector<instance*> lights;
    bvh_tree* bvh = nullptr;
    void* embree_bvh = nullptr;     // unmanaged data for Embree raytracer
    void* embree_device = nullptr;  // unmanaged data for Embree raytracer

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

// Update lights.
void init_lights(
    scene* scn, bool do_shapes = true, bool do_environments = false);
// Generate a distribution for sampling a shape uniformly based on area/length.
void update_shape_cdf(shape* shp);
// Generate a distribution for sampling an environment texture uniformly
// based on angle and texture intensity.
void update_environment_cdf(environment* env);

// Updates/refits bvh.
void build_bvh(shape* shp, bool sah = true);
void build_bvh(scene* scn, bool sah = true);
void refit_bvh(shape* shp);
void refit_bvh(scene* scn);

#if YGL_EMBREE
// Build/update Embree BVH
void build_bvh_embree(scene* scn);
#endif

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
    float width = 0.036f, float height = 0.024f, float focal = 0.050f);
// make default material
inline material* make_default_material(const std::string& name) {
    auto mat = new material();
    mat->name = name;
    mat->kd = {0.2f, 0.2f, 0.2f};
    return mat;
}

// Add a sky environment
inline environment* make_sky_environment(
    const std::string& name, float sun_angle = pi / 4) {
    auto txt = new texture();
    txt->name = name;
    txt->path = "textures/" + name + ".hdr";
    txt->img = make_sunsky_image(1024, 512, sun_angle);
    auto env = new environment();
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
    instance* ist = nullptr;  // instance or null for no intersection
    int ei = 0;               // shape element index
    vec2f uv = zero2f;        // shape element coordinates
    float dist = 0;           // ray/point distance
};

// Intersects a ray with the scene.
scene_intersection intersect_ray(
    const scene* scn, const ray3f& ray, bool find_any = false);

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

// Evaluate a texture.
vec4f eval_texture(const texture* txt, const vec2f& texcoord);

// Set and evaluate camera parameters. Setters take zeros as default values.
float eval_camera_fovy(const camera* cam);
float eval_camera_aspect(const camera* cam);
void set_camera_fovy(
    camera* cam, float fovy, float aspect, float width = 0.036f);
int image_width(const camera* cam, int yresolution);
int image_height(const camera* cam, int yresolution);

// Generates a ray from a camera image coordinate `uv` and lens coordinates
// `luv`.
ray3f eval_camera_ray(const camera* cam, const vec2f& uv, const vec2f& luv);
// Generates a ray from a camera for pixel coordinates `ij`, the resolution
// `res`, the sub-pixel coordinates `puv` and the lens coordinates `luv` and
// the image resolution `res`.
ray3f eval_camera_ray(const camera* cam, int i, int j, int w, int h,
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

// Sample a shape based on a distribution.
std::pair<int, vec2f> sample_shape(
    const shape* shp, float re, const vec2f& ruv);

// Sample an environment uniformly.
vec2f sample_environment(const environment* env, const vec2f& ruv);

}  // namespace ygl

// -----------------------------------------------------------------------------
// PATH TRACING
// -----------------------------------------------------------------------------
namespace ygl {

// Default trace seed
const auto trace_default_seed = 961748941;

// Trace evaluation function.
using trace_func = std::function<vec3f(const scene* scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit)>;

// Progressively compute an image by calling trace_samples multiple times.
image4f trace_image(const scene* scn, const camera* cam, int yresolution,
    int nsamples, trace_func tracer, int nbounces = 8, float pixel_clamp = 100,
    bool noparallel = false, int seed = trace_default_seed);

// Initialize trace rngs
std::vector<rng_state> make_trace_rngs(
    int width, int height, uint64_t seed = trace_default_seed);

// Progressively compute an image by calling trace_samples multiple times.
// Start with an empty state and then successively call this function to
// render the next batch of samples.
void trace_samples(const scene* scn, const camera* cam, int nsamples,
    trace_func tracer, image4f& img, std::vector<rng_state>& rngs, int sample,
    int nbounces = 8, float pixel_clamp = 100, bool noparallel = false,
    int seed = trace_default_seed);

// Starts an anyncrhounous renderer.
void trace_async_start(const scene* scn, const camera* cam, int nsamples,
    trace_func tracer, image4f& img, image4f& display,
    std::vector<rng_state>& rngs, std::vector<std::thread>& threads, bool& stop,
    int& sample, float& exposure, float& gamma, bool& filmic,
    int preview_ratio = 8, int nbounces = 8, float pixel_clamp = 100,
    int seed = trace_default_seed);
// Stop the asynchronous renderer.
void trace_async_stop(std::vector<std::thread>& threads, bool& stop);

// Trace function - path tracer.
vec3f trace_path(const scene* scn, const ray3f& ray, rng_state& rng,
    int nbounces, bool* hit = nullptr);
// Trace function - path tracer without mis.
vec3f trace_path_nomis(const scene* scn, const ray3f& ray, rng_state& rng,
    int nbounces, bool* hit = nullptr);
// Trace function - naive path tracer.
vec3f trace_path_naive(const scene* scn, const ray3f& ray, rng_state& rng,
    int nbounces, bool* hit = nullptr);
// Trace function - direct illumination.
vec3f trace_direct(const scene* scn, const ray3f& ray, rng_state& rng,
    int nbounces, bool* hit = nullptr);
// Trace function - direct illumination without mis.
vec3f trace_direct_nomis(const scene* scn, const ray3f& ray, rng_state& rng,
    int nbounces, bool* hit = nullptr);
// Trace function - pure environment illumination with no shadows.
vec3f trace_environment(const scene* scn, const ray3f& ray, rng_state& rng,
    int nbounces, bool* hit = nullptr);
// Trace function - eyelight rendering.
vec3f trace_eyelight(const scene* scn, const ray3f& ray, rng_state& rng,
    int nbounces, bool* hit = nullptr);
// Trace function - normal debug visualization.
vec3f trace_debug_normal(const scene* scn, const ray3f& ray, rng_state& rng,
    int nbounces, bool* hit = nullptr);
// Trace function - faceforward debug visualization.
vec3f trace_debug_frontfacing(const scene* scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit = nullptr);
// Trace function - albedo debug visualization.
vec3f trace_debug_albedo(const scene* scn, const ray3f& ray, rng_state& rng,
    int nbounces, bool* hit = nullptr);
// Trace function - diffuse debug visualization.
vec3f trace_debug_diffuse(const scene* scn, const ray3f& ray, rng_state& rng,
    int nbounces, bool* hit = nullptr);
// Trace function - specular debug visualization.
vec3f trace_debug_specular(const scene* scn, const ray3f& ray, rng_state& rng,
    int nbounces, bool* hit = nullptr);
// Trace function - roughness debug visualization.
vec3f trace_debug_roughness(const scene* scn, const ray3f& ray, rng_state& rng,
    int nbounces, bool* hit = nullptr);
// Trace function - texcoord debug visualization.
vec3f trace_debug_texcoord(const scene* scn, const ray3f& ray, rng_state& rng,
    int nbounces, bool* hit = nullptr);

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

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR MATRICES
// -----------------------------------------------------------------------------
namespace ygl {

// Matrix diagonals and transposes.
inline mat2f transpose(const mat2f& a) {
    return {{a.x.x, a.y.x}, {a.x.y, a.y.y}};
}
inline mat3f transpose(const mat3f& a) {
    return {
        {a.x.x, a.y.x, a.z.x},
        {a.x.y, a.y.y, a.z.y},
        {a.x.z, a.y.z, a.z.z},
    };
}
inline mat4f transpose(const mat4f& a) {
    return {
        {a.x.x, a.y.x, a.z.x, a.w.x},
        {a.x.y, a.y.y, a.z.y, a.w.y},
        {a.x.z, a.y.z, a.z.z, a.w.z},
        {a.x.w, a.y.w, a.z.w, a.w.w},
    };
}

// Matrix adjugates, determinant and inverses.
inline mat2f adjugate(const mat2f& a) {
    return {{a.y.y, -a.x.y}, {-a.y.x, a.x.x}};
}
inline mat3f adjugate(const mat3f& a) {
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
inline mat4f adjugate(const mat4f& a) {
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
        theta = clamp(theta, 0.001f, pi - 0.001f);
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
        theta = clamp(theta, 0.001f, pi - 0.001f);
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
               ygl::frame3f{frame.x, frame.y, frame.z, vec3f{0, 0, 0}} *
               rotation_frame(vec3f{0, 1, 0}, rotate.x);
    auto pos = frame.o + transl.x * x + transl.y * y + transl.z * z;

    frame = {rot.x, rot.y, rot.z, pos};
}

}  // namespace ygl

#endif
