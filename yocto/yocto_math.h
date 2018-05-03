//
// # Yocto/Math: Tiny C++ Library of Math function for 3D Graphics
//
// Yocto/Math is a collection of simple math types and utilities for 3D
// graphics.
//
// ## Small Vectors and Matrices, Frames, Bounding Boxes and Transforms
//
// We provide common operations for small vectors and matrices typically used
// in graphics. In particular, we support 2-4 dimensional vectors for float
// (`vec2f`, `vec3f`, `vec4f`), int (`vec2i`, `vec3i`, `vec4i`) and bytes
// (`vec4b`). Vector operations are support for float types only.
//
// We support 2-4 dimensional generic matrices `mat2f`, `mat3f`,
// `mat4f`, with matrix-matrix and matrix-vector products, transposes and
// inverses. Matrices are stored in column-major ordered and are accessed and
// constructed by column.
//
// To represent transformations, most of the library facilities prefer the use
// coordinate frames, aka rigid transforms, represented as `frame3f`.
// The structure store three coordinate axis and the frame origin. This is
// equivalent to a rigid transform written as a column-major affine
// matrix. Transform operations are better behaved with this representation.
//
// We represent coordinate bounds with axis-aligned bounding boxes in 3
// dimensions with `bbox3f` with support for expansion operations for
// points and otehr bboxes. We provide operations to compute bounds for points,
// lines, triangles and quads.
//
// For both matrices and frames we support transform operations for points,
// vectors and directions (`trasform_point()`, `trasform_vector()`,
// `trasform_direction()`). For frames we also the support inverse operations
// (`transform_xxx_inverse()`). Transform matrices and frames can be
// constructed from basic translation, rotation and scaling, e.g. with
// `translation_mat4f()` or `translation_frame3f()` respectively, etc. For
// rotation we support axis-angle and quaternions, with slerp.
//
//
// ## Random Number Generation, Noise, and Monte Carlo support
//
// This library supports many facilities helpful in writing sampling
// functions targeting path tracing and shape generations.
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

#ifndef _YGL_MATH_H_
#define _YGL_MATH_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <algorithm>   // for std::upper_bound 
#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <functional>  // for std::hash
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
inline frame2f inverse(const frame2f& a) {
    auto minv = transpose(mat2f{a.x, a.y});
    return {minv.x, minv.y, -(minv * a.o)};
}
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
inline void camera_turntable(vec3f& from, vec3f& to, vec3f& up,
    const vec2f& rotate, float dolly, const vec2f& pan);
inline void camera_turntable(frame3f& frame, float& focus, const vec2f& rotate,
    float dolly, const vec2f& pan);
inline void camera_fps(
    frame3f& frame, const vec3f& transl, const vec2f& rotate);

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

// Init a sequence of random number generators.
inline std::vector<rng_state> make_rng_seq(int num, uint64_t seed) {
    auto rngs = std::vector<rng_state>(num);
    for (auto i = 0; i < num; i++) rngs[i] = make_rng(seed, 2 * i + 1);
    return rngs;
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
    auto idx = (int)(std::upper_bound(cdf.begin(), cdf.end(), r) - cdf.begin());
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
inline float perlin_noise(const vec3f& p, const vec3i& wrap = zero3i);
inline float perlin_ridge_noise(const vec3f& p, float lacunarity = 2.0f,
    float gain = 0.5f, float offset = 1.0f, int octaves = 6,
    const vec3i& wrap = zero3i);
inline float perlin_fbm_noise(const vec3f& p, float lacunarity = 2.0f,
    float gain = 0.5f, int octaves = 6, const vec3i& wrap = zero3i);
inline float perlin_turbulence_noise(const vec3f& p, float lacunarity = 2.0f,
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
inline void camera_fps(
    frame3f& frame, const vec3f& transl, const vec2f& rotate) {
    // https://gamedev.stackexchange.com/questions/30644/how-to-keep-my-quaternion-using-fps-camera-from-tilting-and-messing-up
    auto y = vec3f{0, 1, 0};
    auto z = orthonormalize(frame.z, y);
    auto x = cross(y, z);

    auto rot = rotation_frame(vec3f{1, 0, 0}, rotate.y) *
               frame3f{frame.x, frame.y, frame.z, zero3f} *
               rotation_frame(vec3f{0, 1, 0}, rotate.x);
    auto pos = frame.o + transl.x * x + transl.y * y + transl.z * z;

    frame = {rot.x, rot.y, rot.z, pos};
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PERLIN NOISE
// -----------------------------------------------------------------------------
namespace ygl {

// clang-format off
inline  float stb__perlin_lerp(float a, float b, float t)
{
   return a + (b-a) * t;
}

inline  int stb__perlin_fastfloor(float a)
{
	int ai = (int) a;
	return (a < ai) ? ai-1 : ai;
}

// different grad function from Perlin's, but easy to modify to match reference
inline  float stb__perlin_grad(int hash, float x, float y, float z)
{
   static float basis[12][4] =
   {
      {  1, 1, 0 },
      { -1, 1, 0 },
      {  1,-1, 0 },
      { -1,-1, 0 },
      {  1, 0, 1 },
      { -1, 0, 1 },
      {  1, 0,-1 },
      { -1, 0,-1 },
      {  0, 1, 1 },
      {  0,-1, 1 },
      {  0, 1,-1 },
      {  0,-1,-1 },
   };

   // perlin's gradient has 12 cases so some get used 1/16th of the time
   // and some 2/16ths. We reduce bias by changing those fractions
   // to 5/64ths and 6/64ths, and the same 4 cases get the extra weight.
   static unsigned char indices[64] =
   {
      0,1,2,3,4,5,6,7,8,9,10,11,
      0,9,1,11,
      0,1,2,3,4,5,6,7,8,9,10,11,
      0,1,2,3,4,5,6,7,8,9,10,11,
      0,1,2,3,4,5,6,7,8,9,10,11,
      0,1,2,3,4,5,6,7,8,9,10,11,
   };

   // if you use reference permutation table, change 63 below to 15 to match reference
   // (this is why the ordering of the table above is funky)
   float *grad = basis[indices[hash & 63]];
   return grad[0]*x + grad[1]*y + grad[2]*z;
}

inline float stb_perlin_noise3(float x, float y, float z, int x_wrap, int y_wrap, int z_wrap)
{
    // not same permutation table as Perlin's reference to avoid copyright issues;
    // Perlin's table can be found at http://mrl.nyu.edu/~perlin/noise/
    // @OPTIMIZE: should this be unsigned char instead of int for cache?
    static unsigned char stb__perlin_randtab[512] =
    {
    23, 125, 161, 52, 103, 117, 70, 37, 247, 101, 203, 169, 124, 126, 44, 123,
    152, 238, 145, 45, 171, 114, 253, 10, 192, 136, 4, 157, 249, 30, 35, 72,
    175, 63, 77, 90, 181, 16, 96, 111, 133, 104, 75, 162, 93, 56, 66, 240,
    8, 50, 84, 229, 49, 210, 173, 239, 141, 1, 87, 18, 2, 198, 143, 57,
    225, 160, 58, 217, 168, 206, 245, 204, 199, 6, 73, 60, 20, 230, 211, 233,
    94, 200, 88, 9, 74, 155, 33, 15, 219, 130, 226, 202, 83, 236, 42, 172,
    165, 218, 55, 222, 46, 107, 98, 154, 109, 67, 196, 178, 127, 158, 13, 243,
    65, 79, 166, 248, 25, 224, 115, 80, 68, 51, 184, 128, 232, 208, 151, 122,
    26, 212, 105, 43, 179, 213, 235, 148, 146, 89, 14, 195, 28, 78, 112, 76,
    250, 47, 24, 251, 140, 108, 186, 190, 228, 170, 183, 139, 39, 188, 244, 246,
    132, 48, 119, 144, 180, 138, 134, 193, 82, 182, 120, 121, 86, 220, 209, 3,
    91, 241, 149, 85, 205, 150, 113, 216, 31, 100, 41, 164, 177, 214, 153, 231,
    38, 71, 185, 174, 97, 201, 29, 95, 7, 92, 54, 254, 191, 118, 34, 221,
    131, 11, 163, 99, 234, 81, 227, 147, 156, 176, 17, 142, 69, 12, 110, 62,
    27, 255, 0, 194, 59, 116, 242, 252, 19, 21, 187, 53, 207, 129, 64, 135,
    61, 40, 167, 237, 102, 223, 106, 159, 197, 189, 215, 137, 36, 32, 22, 5,

    // and a second copy so we don't need an extra mask or static initializer
    23, 125, 161, 52, 103, 117, 70, 37, 247, 101, 203, 169, 124, 126, 44, 123,
    152, 238, 145, 45, 171, 114, 253, 10, 192, 136, 4, 157, 249, 30, 35, 72,
    175, 63, 77, 90, 181, 16, 96, 111, 133, 104, 75, 162, 93, 56, 66, 240,
    8, 50, 84, 229, 49, 210, 173, 239, 141, 1, 87, 18, 2, 198, 143, 57,
    225, 160, 58, 217, 168, 206, 245, 204, 199, 6, 73, 60, 20, 230, 211, 233,
    94, 200, 88, 9, 74, 155, 33, 15, 219, 130, 226, 202, 83, 236, 42, 172,
    165, 218, 55, 222, 46, 107, 98, 154, 109, 67, 196, 178, 127, 158, 13, 243,
    65, 79, 166, 248, 25, 224, 115, 80, 68, 51, 184, 128, 232, 208, 151, 122,
    26, 212, 105, 43, 179, 213, 235, 148, 146, 89, 14, 195, 28, 78, 112, 76,
    250, 47, 24, 251, 140, 108, 186, 190, 228, 170, 183, 139, 39, 188, 244, 246,
    132, 48, 119, 144, 180, 138, 134, 193, 82, 182, 120, 121, 86, 220, 209, 3,
    91, 241, 149, 85, 205, 150, 113, 216, 31, 100, 41, 164, 177, 214, 153, 231,
    38, 71, 185, 174, 97, 201, 29, 95, 7, 92, 54, 254, 191, 118, 34, 221,
    131, 11, 163, 99, 234, 81, 227, 147, 156, 176, 17, 142, 69, 12, 110, 62,
    27, 255, 0, 194, 59, 116, 242, 252, 19, 21, 187, 53, 207, 129, 64, 135,
    61, 40, 167, 237, 102, 223, 106, 159, 197, 189, 215, 137, 36, 32, 22, 5,
    };

   float u,v,w;
   float n000,n001,n010,n011,n100,n101,n110,n111;
   float n00,n01,n10,n11;
   float n0,n1;

   unsigned int x_mask = (x_wrap-1) & 255;
   unsigned int y_mask = (y_wrap-1) & 255;
   unsigned int z_mask = (z_wrap-1) & 255;
   int px = stb__perlin_fastfloor(x);
   int py = stb__perlin_fastfloor(y);
   int pz = stb__perlin_fastfloor(z);
   int x0 = px & x_mask, x1 = (px+1) & x_mask;
   int y0 = py & y_mask, y1 = (py+1) & y_mask;
   int z0 = pz & z_mask, z1 = (pz+1) & z_mask;
   int r0,r1, r00,r01,r10,r11;

   #define stb__perlin_ease(a)   (((a*6-15)*a + 10) * a * a * a)

   x -= px; u = stb__perlin_ease(x);
   y -= py; v = stb__perlin_ease(y);
   z -= pz; w = stb__perlin_ease(z);

   r0 = stb__perlin_randtab[x0];
   r1 = stb__perlin_randtab[x1];

   r00 = stb__perlin_randtab[r0+y0];
   r01 = stb__perlin_randtab[r0+y1];
   r10 = stb__perlin_randtab[r1+y0];
   r11 = stb__perlin_randtab[r1+y1];

   n000 = stb__perlin_grad(stb__perlin_randtab[r00+z0], x  , y  , z   );
   n001 = stb__perlin_grad(stb__perlin_randtab[r00+z1], x  , y  , z-1 );
   n010 = stb__perlin_grad(stb__perlin_randtab[r01+z0], x  , y-1, z   );
   n011 = stb__perlin_grad(stb__perlin_randtab[r01+z1], x  , y-1, z-1 );
   n100 = stb__perlin_grad(stb__perlin_randtab[r10+z0], x-1, y  , z   );
   n101 = stb__perlin_grad(stb__perlin_randtab[r10+z1], x-1, y  , z-1 );
   n110 = stb__perlin_grad(stb__perlin_randtab[r11+z0], x-1, y-1, z   );
   n111 = stb__perlin_grad(stb__perlin_randtab[r11+z1], x-1, y-1, z-1 );

   n00 = stb__perlin_lerp(n000,n001,w);
   n01 = stb__perlin_lerp(n010,n011,w);
   n10 = stb__perlin_lerp(n100,n101,w);
   n11 = stb__perlin_lerp(n110,n111,w);

   n0 = stb__perlin_lerp(n00,n01,v);
   n1 = stb__perlin_lerp(n10,n11,v);

   return stb__perlin_lerp(n0,n1,u);
}

inline float stb_perlin_ridge_noise3(float x, float y, float z,float lacunarity, float gain, float offset, int octaves,int x_wrap, int y_wrap, int z_wrap)
{
   int i;
   float frequency = 1.0f;
   float prev = 1.0f;
   float amplitude = 0.5f;
   float sum = 0.0f;

   for (i = 0; i < octaves; i++) {
      float r = (float)(stb_perlin_noise3(x*frequency,y*frequency,z*frequency,x_wrap,y_wrap,z_wrap));
      r = r<0 ? -r : r; // fabs()
      r = offset - r;
      r = r*r;
      sum += r*amplitude*prev;
      prev = r;
      frequency *= lacunarity;
      amplitude *= gain;
   }
   return sum;
}

inline float stb_perlin_fbm_noise3(float x, float y, float z,float lacunarity, float gain, int octaves,int x_wrap, int y_wrap, int z_wrap)
{
   int i;
   float frequency = 1.0f;
   float amplitude = 1.0f;
   float sum = 0.0f;

   for (i = 0; i < octaves; i++) {
      sum += stb_perlin_noise3(x*frequency,y*frequency,z*frequency,x_wrap,y_wrap,z_wrap)*amplitude;
      frequency *= lacunarity;
      amplitude *= gain;
   }
   return sum;
}

inline float stb_perlin_turbulence_noise3(float x, float y, float z, float lacunarity, float gain, int octaves,int x_wrap, int y_wrap, int z_wrap)
{
   int i;
   float frequency = 1.0f;
   float amplitude = 1.0f;
   float sum = 0.0f;

   for (i = 0; i < octaves; i++) {
      float r = stb_perlin_noise3(x*frequency,y*frequency,z*frequency,x_wrap,y_wrap,z_wrap)*amplitude;
      r = r<0 ? -r : r; // fabs()
      sum += r;
      frequency *= lacunarity;
      amplitude *= gain;
   }
   return sum;
}
// clang-format on

// adapeted  stb_perlin.h
inline float perlin_noise(const vec3f& p, const vec3i& wrap) {
    return stb_perlin_noise3(p.x, p.y, p.z, wrap.x, wrap.y, wrap.z);
}

// adapeted  stb_perlin.h
inline float perlin_ridge_noise(const vec3f& p, float lacunarity, float gain,
    float offset, int octaves, const vec3i& wrap) {
    return stb_perlin_ridge_noise3(p.x, p.y, p.z, lacunarity, gain, offset,
        octaves, wrap.x, wrap.y, wrap.z);
}

// adapeted  stb_perlin.h
inline float perlin_fbm_noise(const vec3f& p, float lacunarity, float gain,
    int octaves, const vec3i& wrap) {
    return stb_perlin_fbm_noise3(
        p.x, p.y, p.z, lacunarity, gain, octaves, wrap.x, wrap.y, wrap.z);
}

// adapeted  stb_perlin.h
inline float perlin_turbulence_noise(const vec3f& p, float lacunarity,
    float gain, int octaves, const vec3i& wrap) {
    return stb_perlin_turbulence_noise3(
        p.x, p.y, p.z, lacunarity, gain, octaves, wrap.x, wrap.y, wrap.z);
}

}  // namespace ygl

#endif
