//
// # Yocto/Math: Tiny C++ Library of Math function for 3D Graphics
//
// Yocto/Math is a collection of simple math types and utilities for 3D
// graphics.
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

#include <algorithm>  // for std::upper_bound
#include <cctype>
#include <cfloat>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <functional>  // for std::hash
#include <iostream>
#include <limits>
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
template <typename T, int N>
struct mat;

// Small Fixed-size square matrices stored in column major format.
template <typename T>
struct mat<T, 2> {
    vec<T, 2> x = {1, 0};
    vec<T, 2> y = {0, 1};
};
template <typename T>
struct mat<T, 3> {
    vec<T, 3> x = {1, 0, 0};
    vec<T, 3> y = {0, 1, 0};
    vec<T, 3> z = {0, 0, 1};
};
template <typename T>
struct mat<T, 4> {
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
    T min = std::numeric_limits<T>::max();
    T max = std::numeric_limits<T>::lowest();
};

// Axis aligned bounding box represented as a min/max vector pairs.
template <typename T>
struct bbox<T, 2> {
    vec<T, 2> min = {
        std::numeric_limits<T>::max(), std::numeric_limits<T>::max()};
    vec<T, 2> max = {
        std::numeric_limits<T>::lowest(), std::numeric_limits<T>::lowest()};
};

// Axis aligned bounding box represented as a min/max vector pairs.
template <typename T>
struct bbox<T, 3> {
    vec<T, 3> min = {std::numeric_limits<T>::max(),
        std::numeric_limits<T>::max(), std::numeric_limits<T>::max()};
    vec<T, 3> max = {std::numeric_limits<T>::lowest(),
        std::numeric_limits<T>::lowest(), std::numeric_limits<T>::lowest()};
};

// Axis aligned bounding box represented as a min/max vector pairs.
template <typename T>
struct bbox<T, 4> {
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
inline bbox<T, 3> line_bbox(const vec<T, 3>& v0, const vec<T, 3>& v1, T1 r0 = 0, T1 r1 = 0) {
    auto bbox = ygl::bbox<T, 3>{};
    bbox += v0 - vec<T, 3>{r0, r0, r0};
    bbox += v0 + vec<T, 3>{r0, r0, r0};
    bbox += v1 - vec<T, 3>{r1, r1, r1};
    bbox += v1 + vec<T, 3>{r1, r1, r1};
    return bbox;
}
template <typename T>
inline bbox<T, 3> triangle_bbox(const vec<T, 3>& v0, const vec<T, 3>& v1, const vec<T, 3>& v2) {
    auto bbox = ygl::bbox<T, 3>{};
    bbox += v0;
    bbox += v1;
    bbox += v2;
    return bbox;
}
template <typename T>
inline bbox<T, 3> quad_bbox(
    const vec<T, 3>& v0, const vec<T, 3>& v1, const vec<T, 3>& v2, const vec<T, 3>& v3) {
    auto bbox = ygl::bbox<T, 3>{};
    bbox += v0;
    bbox += v1;
    bbox += v2;
    bbox += v3;
    return bbox;
}
template <typename T>
inline bbox<T, 3> tetrahedron_bbox(
    const vec<T, 3>& v0, const vec<T, 3>& v1, const vec<T, 3>& v2, const vec<T, 3>& v3) {
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
    vec<T, 2> o = {0, 0};
    vec<T, 2> d = {0, 1};
    T tmin = 0;
    T tmax = flt_max;
};

// Rays with origin, direction and min/max t value.
template <typename T>
struct ray<T, 3> {
    vec<T, 3> o = {0, 0, 0};
    vec<T, 3> d = {0, 0, 1};
    T tmin = 0;
    T tmax = flt_max;
};

// Type aliases
using ray2f = ray<float, 2>;
using ray3f = ray<float, 3>;

// Construct a ray from dirction or segments using a default epsilon.
template <typename T, int N>
inline ray<T, N> make_ray(const vec<T, N>& o, const vec<T, N>& d, T eps = (T)1e-4) {
    return ray<T, N>{o, d, eps, std::numeric_limits<T>::max()};
}
template <typename T, int N>
inline ray<T, N> make_segment(const vec<T, N>& p1, const vec<T, N>& p2, T eps = (T)1e-4) {
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
inline frame<T, 3> lookat_frame(
    const vec<T, 3>& eye, const vec<T, 3>& center, const vec<T, 3>& up, bool inv_xz = false) {
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
inline void camera_turntable(frame<T, 3>& frame, float& focus, const vec<T, 2>& rotate,
    T1 dolly, const vec<T, 2>& pan);
template <typename T>
inline void camera_fps(frame<T, 3>& frame, const vec<T, 3>& transl, const vec<T, 2>& rotate);

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
inline vec<T, 3> sample_triangle(
    const vec<T, 3>& v0, const vec<T, 3>& v1, const vec<T, 3>& v2, const vec<T, 2>& ruv) {
    auto uv = sample_triangle(ruv);
    return v0 * (1 - uv.x - uv.y) + v1 * uv.x + v2 * uv.y;
}
// Pdf for uniform triangle sampling, i.e. triangle area.
template <typename T>
inline T sample_triangle_pdf(const vec<T, 3>& v0, const vec<T, 3>& v1, const vec<T, 3>& v2) {
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
inline vec<T, 3> triangle_normal(const vec<T, 3>& v0, const vec<T, 3>& v1, const vec<T, 3>& v2) {
    return normalize(cross(v1 - v0, v2 - v0));
}
template <typename T>
inline T triangle_area(const vec<T, 3>& v0, const vec<T, 3>& v1, const vec<T, 3>& v2) {
    return length(cross(v1 - v0, v2 - v0)) / 2;
}

// Quad propeties.
template <typename T>
inline vec<T, 3> quad_normal(
    const vec<T, 3>& v0, const vec<T, 3>& v1, const vec<T, 3>& v2, const vec<T, 3>& v3) {
    return normalize(triangle_normal(v0, v1, v3) + triangle_normal(v2, v3, v1));
}
template <typename T>
inline T quad_area(const vec<T, 3>& v0, const vec<T, 3>& v1, const vec<T, 3>& v2, const vec<T, 3>& v3) {
    return triangle_area(v0, v1, v3) + triangle_area(v2, v3, v1);
}

// Triangle tangent and bitangent from uv
template <typename T>
inline std::pair<vec<T, 3>, vec<T, 3>> triangle_tangents_fromuv(const vec<T, 3>& v0,
    const vec<T, 3>& v1, const vec<T, 3>& v2, const vec<T, 2>& uv0, const vec<T, 2>& uv1, const vec<T, 2>& uv2) {
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
inline TT interpolate_line(const std::vector<TT>& vals, const vec<int, 2>& l, T u) {
    if (vals.empty()) return TT();
    return vals[l.x] * (1 - u) + vals[l.y] * u;
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
    return vals[t.x] * (1 - uv.x - uv.y) + vals[t.y] * uv.x + vals[t.z] * uv.y;
}
// Interpolates values over a quad parametrized by u and v along the
// (v1-v0) and (v2-v1) directions. Same as bilear interpolation.
template <typename T>
inline T interpolate_quad(
    const T& v0, const T& v1, const T& v2, const T& v3, const vec<T, 2>& uv) {
    return v0 * (1 - uv.x) * (1 - uv.y) + v1 * uv.x * (1 - uv.y) +
           v2 * uv.x * uv.y + v3 * (1 - uv.x) * uv.y;
}
template <typename T>
inline T interpolate_quad(const std::vector<T>& vals, vec4i t, const vec<T, 2>& uv) {
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
inline void camera_turntable(frame<T, 3>& frame, float& focus, const vec<T, 2>& rotate,
    T1 dolly, const vec<T, 2>& pan) {
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
inline void camera_fps(frame<T, 3>& frame, const vec<T, 3>& transl, const vec<T, 2>& rotate) {
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
inline float perlin_noise(vec3f p, vec3i wrap) {
    return stb_perlin_noise3(p.x, p.y, p.z, wrap.x, wrap.y, wrap.z);
}

// adapeted  stb_perlin.h
inline float perlin_ridge_noise(vec3f p, float lacunarity, float gain,
    float offset, int octaves, vec3i wrap) {
    return stb_perlin_ridge_noise3(p.x, p.y, p.z, lacunarity, gain, offset,
        octaves, wrap.x, wrap.y, wrap.z);
}

// adapeted  stb_perlin.h
inline float perlin_fbm_noise(
    vec3f p, float lacunarity, float gain, int octaves, vec3i wrap) {
    return stb_perlin_fbm_noise3(
        p.x, p.y, p.z, lacunarity, gain, octaves, wrap.x, wrap.y, wrap.z);
}

// adapeted  stb_perlin.h
inline float perlin_turbulence_noise(
    vec3f p, float lacunarity, float gain, int octaves, vec3i wrap) {
    return stb_perlin_turbulence_noise3(
        p.x, p.y, p.z, lacunarity, gain, octaves, wrap.x, wrap.y, wrap.z);
}

}  // namespace ygl

#endif
