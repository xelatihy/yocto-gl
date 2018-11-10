//
// # Yocto/Math: Tiny library for math support in graphics applications.
//
// Yocto/Math provides the basic math primitives used in grahics, including
// small-sized vectors and matrixes, frames, bounding boxes and transforms.
// We also include support for random number generation, and primitives for 
// Monte Carlo sampling.
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
// LICENCE OF INCLUDED SOFTWARE FOR PERLIN NOISE
// https://github.com/nothings/stb/blob/master/stb_perlin.h
//
// ------------------------------------------------------------------------------
// ALTERNATIVE B - Public Domain (www.unlicense.org)
// This is free and unencumbered software released into the public domain.
// Anyone is free to copy, modify, publish, use, compile, sell, or distribute
// this software, either in source code form or as a compiled binary, for any
// purpose, commercial or non-commercial, and by any means. In jurisdictions
// that recognize copyright laws, the author or authors of this software
// dedicate any and all copyright interest in the software to the public domain.
// We make this dedication for the benefit of the public at large and to the
// detriment of our heirs and successors. We intend this dedication to be an
// overt act of relinquishment in perpetuity of all present and future rights to
// this software under copyright law.
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
// ------------------------------------------------------------------------------
//
//

#ifndef _YOCTO_MATH_H_
#define _YOCTO_MATH_H_

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
#include <future>
#include <limits>
#include <mutex>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

// -----------------------------------------------------------------------------
// MATH CONSTANTS AND FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

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
using std::get;
using std::ignore;
using std::lock_guard;
using std::make_unique;
using std::mutex;
using std::packaged_task;
using std::pair;
using std::runtime_error;
using std::shared_future;
using std::string;
using std::thread;
using std::tie;
using std::tuple;
using std::unordered_map;
using std::vector;
using namespace std::string_literals;
using namespace std::chrono_literals;

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

}  // namespace yocto

// -----------------------------------------------------------------------------
// VECTORS
// -----------------------------------------------------------------------------
namespace yocto {

// Small size vectors.
template <typename T, int N>
struct vec;

// Small size vectors.
template <typename T>
struct vec<T, 1> {
    T x = 0;

    constexpr T&       operator[](int idx) { return *(&x + idx); }
    constexpr const T& operator[](int idx) const { return *(&x + idx); }
};
template <typename T>
struct vec<T, 2> {
    T x = 0;
    T y = 0;

    constexpr T&       operator[](int idx) { return *(&x + idx); }
    constexpr const T& operator[](int idx) const { return *(&x + idx); }
};
template <typename T>
struct vec<T, 3> {
    T x = 0;
    T y = 0;
    T z = 0;

    constexpr T&       operator[](int idx) { return *(&x + idx); }
    constexpr const T& operator[](int idx) const { return *(&x + idx); }
};
template <typename T>
struct vec<T, 4> {
    T x = 0;
    T y = 0;
    T z = 0;
    T w = 0;

    constexpr T&       operator[](int idx) { return *(&x + idx); }
    constexpr const T& operator[](int idx) const { return *(&x + idx); }
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
    return &v[0];
}
template <typename T, int N>
constexpr const T* begin(const vec<T, N>& v) {
    return &v[0];
}
template <typename T, int N>
constexpr T* end(vec<T, N>& v) {
    return &v[0] + N;
}
template <typename T, int N>
constexpr const T* end(const vec<T, N>& v) {
    return &v[0] + N;
}
template <typename T, int N>
constexpr T* data(vec<T, N>& v) {
    return &v[0];
}
template <typename T, int N>
constexpr const T* data(const vec<T, N>& v) {
    return &v[0];
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

}  // namespace yocto

namespace std {

// Hash functor for vector for use with unordered_map
template <typename T, int N>
struct hash<yocto::vec<T, N>> {
    size_t operator()(const yocto::vec<T, N>& v) const {
        auto vh = hash<T>();
        auto h  = (size_t)0;
        for (auto i = 0; i < N; i++)
            h ^= vh(v[i]) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

}  // namespace std

// -----------------------------------------------------------------------------
// MATRICES
// -----------------------------------------------------------------------------
namespace yocto {

// Small Fixed-size square matrices stored in column major format.
template <typename T, int N, int M>
struct mat;

// Small Fixed-size square matrices stored in column major format.
template <typename T, int N>
struct mat<T, N, 1> {
    vec<T, N> x = {};

    constexpr vec<T, N>&       operator[](int idx) { return *(&x + idx); }
    constexpr const vec<T, N>& operator[](int idx) const { return *(&x + idx); }
};
template <typename T, int N>
struct mat<T, N, 2> {
    vec<T, N> x = {};
    vec<T, N> y = {};

    constexpr vec<T, N>&       operator[](int idx) { return *(&x + idx); }
    constexpr const vec<T, N>& operator[](int idx) const { return *(&x + idx); }
};
template <typename T, int N>
struct mat<T, N, 3> {
    vec<T, 3> x = {};
    vec<T, 3> y = {};
    vec<T, 3> z = {};

    constexpr vec<T, N>&       operator[](int idx) { return *(&x + idx); }
    constexpr const vec<T, N>& operator[](int idx) const { return *(&x + idx); }
};
template <typename T, int N>
struct mat<T, N, 4> {
    vec<T, 4> x = {};
    vec<T, 4> y = {};
    vec<T, 4> z = {};
    vec<T, 4> w = {};

    constexpr vec<T, N>&       operator[](int idx) { return *(&x + idx); }
    constexpr const vec<T, N>& operator[](int idx) const { return *(&x + idx); }
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

}  // namespace yocto

// -----------------------------------------------------------------------------
// RIGID BODY TRANSFORMS/FRAMES
// -----------------------------------------------------------------------------
namespace yocto {

// Rigid frames stored as a column-major affine transform matrix.
template <typename T, int N>
struct frame;

// Rigid frames stored as a column-major affine transform matrix.
template <typename T>
struct frame<T, 2> {
    vec<T, 2> x = {1, 0};
    vec<T, 2> y = {0, 1};
    vec<T, 2> o = {0, 0};

    constexpr vec<T, 2>&       operator[](int idx) { return *(&x + idx); }
    constexpr const vec<T, 2>& operator[](int idx) const { return *(&x + idx); }
};
template <typename T>
struct frame<T, 3> {
    vec<T, 3> x = {1, 0, 0};
    vec<T, 3> y = {0, 1, 0};
    vec<T, 3> z = {0, 0, 1};
    vec<T, 3> o = {0, 0, 0};

    constexpr vec<T, 3>&       operator[](int idx) { return *(&x + idx); }
    constexpr const vec<T, 3>& operator[](int idx) const { return *(&x + idx); }
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

}  // namespace yocto

// -----------------------------------------------------------------------------
// AXIS ALIGNED BOUNDING BOXES
// -----------------------------------------------------------------------------
namespace yocto {

// DEscribes a range of values in N dimensions.
template <typename T, int N>
struct bbox;

// Range of values in 1D.
template <typename T>
struct bbox<T, 1> {
    T min = maxt<T>();
    T max = mint<T>();

    constexpr T&       operator[](int idx) { return *(&min + idx); }
    constexpr const T& operator[](int idx) const { return *(&min + idx); }
};

// Axis aligned bounding box represented as a min/max vector pairs.
template <typename T>
struct bbox<T, 2> {
    vec<T, 2> min = {maxt<T>(), maxt<T>()};
    vec<T, 2> max = {mint<T>(), mint<T>()};

    constexpr T&       operator[](int idx) { return *(&min + idx); }
    constexpr const T& operator[](int idx) const { return *(&min + idx); }
};

// Axis aligned bounding box represented as a min/max vector pairs.
template <typename T>
struct bbox<T, 3> {
    vec<T, 3> min = {maxt<T>(), maxt<T>(), maxt<T>()};
    vec<T, 3> max = {mint<T>(), mint<T>(), mint<T>()};

    constexpr T&       operator[](int idx) { return *(&min + idx); }
    constexpr const T& operator[](int idx) const { return *(&min + idx); }
};

// Axis aligned bounding box represented as a min/max vector pairs.
template <typename T>
struct bbox<T, 4> {
    vec<T, 4> min = {maxt<T>(), maxt<T>(), maxt<T>(), maxt<T>()};
    vec<T, 4> max = {mint<T>(), mint<T>(), mint<T>(), mint<T>()};

    constexpr T&       operator[](int idx) { return *(&min + idx); }
    constexpr const T& operator[](int idx) const { return *(&min + idx); }
};

// Type aliases
using bbox1f = bbox<float, 1>;
using bbox2f = bbox<float, 2>;
using bbox3f = bbox<float, 3>;
using bbox4f = bbox<float, 4>;
using bbox1i = bbox<int, 1>;
using bbox2i = bbox<int, 2>;
using bbox3i = bbox<int, 3>;
using bbox4i = bbox<int, 4>;

// Empty bbox constant.
const auto invalid_bbox1f = bbox1f();
const auto invalid_bbox2f = bbox2f();
const auto invalid_bbox3f = bbox3f();
const auto invalid_bbox4f = bbox4f();
const auto invalid_bbox1i = bbox1i();
const auto invalid_bbox2i = bbox2i();
const auto invalid_bbox3i = bbox3i();
const auto invalid_bbox4i = bbox4i();

// Bounding box size and center
template <typename T, int N>
vec<T, N> bbox_size(const bbox<T, N>& a) {
    return a.max - a.min;
}
template <typename T, int N>
vec<T, N> bbox_center(const bbox<T, N>& a) {
    return (a.min + a.max) / 2;
}

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

}  // namespace yocto

// -----------------------------------------------------------------------------
// RAYS
// -----------------------------------------------------------------------------
namespace yocto {

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

// Default ray epsilon
const auto default_ray_epsf = 1e-4f;

// Construct a ray from direction or segments using a default epsilon.
template <typename T, int N>
constexpr inline ray<T, N> make_ray(
    const vec<T, N>& o, const vec<T, N>& d, T eps = default_ray_epsf) {
    return {o, d, eps, maxt<T>()};
}
template <typename T, int N>
constexpr inline ray<T, N> make_segment(
    const vec<T, N>& p1, const vec<T, N>& p2, T eps = default_ray_epsf) {
    return {p1, normalize(p2 - p1), eps, length(p2 - p1) - 2 * eps};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// TRANSFORMS
// -----------------------------------------------------------------------------
namespace yocto {

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
constexpr inline pair<vec<T, 3>, T> rotation_axisangle(const vec<T, 4>& quat) {
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

}  // namespace yocto

// -----------------------------------------------------------------------------
// RANDOM NUMBER GENERATION
// -----------------------------------------------------------------------------
namespace yocto {

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

}  // namespace yocto

// -----------------------------------------------------------------------------
// MONETACARLO SAMPLING FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

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

}  // namespace yocto

// -----------------------------------------------------------------------------
// PERLIN NOISE FUNCTION
// -----------------------------------------------------------------------------
namespace yocto {

// Compute the revised Perlin noise function. Wrap provides a wrapping noise
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

}  // namespace yocto

// -----------------------------------------------------------------------------
// ANIMATION UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

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

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR MATRICES
// -----------------------------------------------------------------------------
namespace yocto {

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

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR UI UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

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

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PERLIN NOISE
// -----------------------------------------------------------------------------
namespace yocto {

// clang-format off
inline float stb__perlin_lerp(float a, float b, float t)
{
   return a + (b-a) * t;
}

inline int stb__perlin_fastfloor(float a)
{
    int ai = (int) a;
    return (a < ai) ? ai-1 : ai;
}

// different grad function from Perlin's, but easy to modify to match reference
inline float stb__perlin_grad(int hash, float x, float y, float z)
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
inline float perlin_turbulence_noise(const vec3f& p, float lacunarity, float gain,
    int octaves, const vec3i& wrap) {
    return stb_perlin_turbulence_noise3(
        p.x, p.y, p.z, lacunarity, gain, octaves, wrap.x, wrap.y, wrap.z);
}

}  // namespace yocto

#endif
