//
// YOCTO_MATH: a collection of vector math functions and simple containers
// used to implement YOCTO. Features include
// - static length float vectors, with specialization for 2, 3, 4 length
// - static length matrices, with specialization for 2x2, 3x3, 4x4
// - affine and rigid transforms
// - linear algebra operations and transforms for fixed length matrices/vecs
// - axis aligned bounding boxes
// - rays
// - random number generation via PCG32
// - a few hash functions
// - timer (depends on C++11 chrono)
//
// While we tested this library in the implementation of our other ones, we
// consider this code incomplete and remommend to use a more complete math
// library. So use it at your own peril.
//
// We developed our own library since we felt that all existing ones are either
// complete, but unreadable or with lots of dependencies, or just as incomplete
// and untested as ours.
//

//
// COMPILATION:
//
// This library can only be used as a header only library in C++ since it uses
// templates for its basic types. To use STL containers in YOCTO, instead of
// the built in ones, #define YGL_USESTL before including this file or any
// other YOCTO file that dependes on this.
//

//
// HISTORY:
// - v 0.5: simplification of constructors, raname bbox -> bbox
// - v 0.4: overall type simplification
// - v 0.3: internal C++ refactoring
// - v 0.2: use of STL containers; removal of yocto containers
// - v 0.1: C++ only implementation
// - v 0.0: initial release in C99
//

//
// ACKNOLEDGEMENTS
//
// This library includes code from the PCG random number generator,
// boost hash_combine, Pixar multijittered sampling, code from "Real-Time
// Collision Detection" by Christer Ericson and public domain code from
// - https://github.com/sgorsten/linalg
// - https://gist.github.com/badboy/6267743
//
// For design ideas of a modern math library, it was very helpful to look at
// - https://github.com/sgorsten/linalg
//

//
// LICENSE:
//
// Copyright (c) 2016 Fabio Pellacini
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
//  LICENSE of included software
//
// This code also includes a small exerpt from http://www.pcg-random.org/
// licensed as follows
// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)
//

#ifndef _YMATH_H_
#define _YMATH_H_

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <functional>
#include <initializer_list>
#include <limits>
#include <type_traits>
#include <unordered_map>
#include <vector>

// -----------------------------------------------------------------------------
// CONSTANTS
// -----------------------------------------------------------------------------

namespace ym {

//
// Using directives
//
using std::fmod;
using std::sin;
using std::cos;
using std::tan;
using std::asin;
using std::acos;
using std::atan;
using std::atan2;
using std::sqrt;
using std::pow;
using std::abs;

// types from the standard library for consistency
using size_t = std::size_t;

// standard containers
using std::vector;
using std::array;
using std::initializer_list;
using std::function;
using std::pair;
using std::unordered_map;

// pi (float)
const float pif = 3.14159265f;
// pi (double)
const double pi = 3.1415926535897932384626433832795;

// shortcat for max values
const auto maxf = std::numeric_limits<float>::max();
const auto maxi = std::numeric_limits<int>::max();

// -----------------------------------------------------------------------------
// BASIC MATH FUNCTIONS
// -----------------------------------------------------------------------------

// Safe minimum value.
template <typename T>
constexpr inline T min(const T& x, const T& y) {
    return (x < y) ? x : y;
}

// Safe maximum value.
template <typename T>
constexpr inline T max(const T& x, const T& y) {
    return (x > y) ? x : y;
}

// Swap values.
template <typename T>
constexpr inline void swap(T& a, T& b) {
    T c = a;
    a = b;
    b = c;
}

// Clamp a value between a minimum and a maximum.
template <typename T>
constexpr inline T clamp(const T& x, const T& min_, const T& max_) {
    return min(max(x, min_), max_);
}

// Linear interpolation.
template <typename T, typename T1>
constexpr inline T lerp(const T& a, const T& b, T1 t) {
    return a * (1 - t) + b * t;
}

// Integer power of two
constexpr inline int pow2(int x) { return 1 << x; }

// -----------------------------------------------------------------------------
// VECTORS
// -----------------------------------------------------------------------------

//
// Vector of element of compile time length with default initializer,
// constant initialization and initialization from a C array. Data access
// via operator[].
//
template <typename T, size_t N>
struct vec : array<T, N> {
    // default constructor
    vec() : array<T, N>{} {}

    // member constructor
    template <typename... Args,
              typename = std::enable_if_t<sizeof...(Args) == N>>
    constexpr vec(Args... vv) : array<T, N>{T(vv)...} {}

    // list constructor
    constexpr vec(initializer_list<T> vv) {
        assert(N == vv.size());
        auto i = 0;
        for (auto&& e : vv) (*this)[i++] = e;
    }

    // copy constructor form array
    constexpr vec(const array<T, N>& vv) : array<T, N>{vv} {}
};

//
// Typedef for standard vectors.
//

using vec1f = vec<float, 1>;
using vec2f = vec<float, 2>;
using vec3f = vec<float, 3>;
using vec4f = vec<float, 4>;

using vec1i = vec<int, 1>;
using vec2i = vec<int, 2>;
using vec3i = vec<int, 3>;
using vec4i = vec<int, 4>;

using vec2b = vec<unsigned char, 2>;
using vec3b = vec<unsigned char, 3>;
using vec4b = vec<unsigned char, 4>;

//
// Vector Constants.
//

const auto zero2f = vec2f();
const auto zero3f = vec3f();
const auto zero4f = vec4f();

const auto zero2i = vec2i();
const auto zero3i = vec3i();
const auto zero4i = vec4i();

//
// Component-wise arithmentic.
//

template <typename T, size_t N>
inline vec<T, N> operator+(const vec<T, N>& a) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = +a[i];
    return c;
}

template <typename T, size_t N>
inline vec<T, N> operator-(const vec<T, N>& a) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = -a[i];
    return c;
}

template <typename T, size_t N>
inline vec<T, N> operator+(const vec<T, N>& a, const vec<T, N>& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a[i] + b[i];
    return c;
}

template <typename T, size_t N, typename T1>
inline vec<T, N> operator+(const vec<T, N>& a, const T1& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a[i] + b;
    return c;
}

template <typename T, size_t N, typename T1>
inline vec<T, N> operator+(const T1& a, const vec<T, N>& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a + b[i];
    return c;
}

template <typename T, size_t N>
inline vec<T, N> operator-(const vec<T, N>& a, const vec<T, N>& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a[i] - b[i];
    return c;
}

template <typename T, size_t N, typename T1>
inline vec<T, N> operator-(const vec<T, N>& a, const T1& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a[i] - b;
    return c;
}

template <typename T, size_t N, typename T1>
inline vec<T, N> operator-(const T1& a, const vec<T, N>& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a - b[i];
    return c;
}

template <typename T, size_t N>
inline vec<T, N> operator*(const vec<T, N>& a, const vec<T, N>& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a[i] * b[i];
    return c;
}

template <typename T, size_t N, typename T1>
inline vec<T, N> operator*(const vec<T, N>& a, const T1& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a[i] * b;
    return c;
}

template <typename T1, typename T, size_t N>
inline vec<T, N> operator*(const T1& a, const vec<T, N>& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a * b[i];
    return c;
}

template <typename T, size_t N>
inline vec<T, N> operator/(const vec<T, N>& a, const vec<T, N>& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a[i] / b[i];
    return c;
}

template <typename T, size_t N, typename T1>
inline vec<T, N> operator/(const vec<T, N>& a, const T1 b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a[i] / b;
    return c;
}

template <typename T1, typename T, size_t N>
inline vec<T, N> operator/(const T1& a, const vec<T, N>& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a / b[i];
    return c;
}

//
// Component-wise assignment arithmentic.
//

template <typename T, size_t N>
inline vec<T, N>& operator+=(vec<T, N>& a, const vec<T, N>& b) {
    return a = a + b;
}

template <typename T, size_t N>
inline vec<T, N>& operator-=(vec<T, N>& a, const vec<T, N>& b) {
    return a = a - b;
}

template <typename T, size_t N>
inline vec<T, N>& operator*=(vec<T, N>& a, const vec<T, N>& b) {
    return a = a * b;
}

template <typename T, size_t N, typename T1>
inline vec<T, N>& operator*=(vec<T, N>& a, const T1& b) {
    return a = a * b;
}

template <typename T, size_t N>
inline vec<T, N>& operator/=(vec<T, N>& a, const vec<T, N>& b) {
    return a = a / b;
}

template <typename T, size_t N, typename T1>
inline vec<T, N>& operator/=(vec<T, N>& a, const T1 b) {
    return a = a / b;
}

//
// Vector operations
//
template <typename T, size_t N>
inline T dot(const vec<T, N>& a, const vec<T, N>& b) {
    auto c = T(0);
    for (auto i = 0; i < N; i++) c += a[i] * b[i];
    return c;
}

template <typename T, size_t N>
inline T length(const vec<T, N>& a) {
    return sqrt(dot(a, a));
}

template <typename T, size_t N>
inline T lengthsqr(const vec<T, N>& a) {
    return dot(a, a);
}

template <typename T, size_t N>
inline vec<T, N> normalize(const vec<T, N>& a) {
    auto l = length(a);
    if (l == 0) return a;
    return a * (1 / l);
}

template <typename T, size_t N>
inline T dist(const vec<T, N>& a, const vec<T, N>& b) {
    return length(a - b);
}

template <typename T, size_t N>
inline T distsqr(const vec<T, N>& a, const vec<T, N>& b) {
    return lengthsqr(a - b);
}

template <typename T>
inline T cross(const vec<T, 2>& a, const vec<T, 2>& b) {
    return a[0] * b[1] - a[1] * b[0];
}

template <typename T>
inline vec<T, 3> cross(const vec<T, 3>& a, const vec<T, 3>& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]};
}

template <typename T, size_t N>
T uangle(const vec<T, N>& a, const vec<T, N>& b) {
    auto d = dot(a, b);
    return d > 1 ? 0 : acos(d < -1 ? -1 : d);
}

template <typename T, size_t N>
T angle(const vec<T, N>& a, const vec<T, N>& b) {
    return uangle(normalize(a), normalize(b));
}

template <typename T, size_t N>
inline vec<T, N> nlerp(const vec<T, N>& a, const vec<T, N>& b, T t) {
    return normalize(lerp(a, b, t));
}

template <typename T, size_t N>
inline vec<T, N> slerp(const vec<T, N>& a, const vec<T, N>& b, T t) {
    auto th = uangle(a, b);
    return th == 0 ? a
                   : a * (sin(th * (1 - t)) / sin(th)) +
                         b * (sin(th * t) / sin(th));
}

// http://lolengine.net/blog/2013/09/21/picking-orthogonal-vector-combing-coconuts)
template <typename T>
inline vec<T, 3> orthogonal(const vec<T, 3>& v) {
    return abs(v[0]) > abs(v[2]) ? vec<T, 3>{-v[1], v[0], 0}
                                 : vec<T, 3>{0, -v[2], v[1]};
}

template <typename T>
inline vec<T, 3> orthonormalize(const vec<T, 3>& a, const vec<T, 3>& b) {
    return normalize(a - b * ym_dot(a, b));
}

//
// Aggregate operations.
//
template <typename T, size_t M>
inline T sum(const vec<T, M>& a) {
    auto s = T(0);
    for (auto i = 0; i < M; i++) s += a[i];
    return s;
}

template <typename T, size_t M>
inline T mean(const vec<T, M>& a) {
    return sum(a) / M;
}

//
// Component-wise min/max and clamping.
//

template <typename T, size_t N>
inline vec<T, N> min(const vec<T, N>& a, const vec<T, N>& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = min(a[i], b[i]);
    return c;
}

template <typename T, size_t N>
inline vec<T, N> max(const vec<T, N>& a, const vec<T, N>& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = max(a[i], b[i]);
    return c;
}

template <typename T, size_t N>
inline vec<T, N> min(const vec<T, N>& a, const T& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = min(a[i], b);
    return c;
}

template <typename T, size_t N>
inline vec<T, N> max(const vec<T, N>& a, const T& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = max(a[i], b);
    return c;
}

template <typename T, size_t N>
inline vec<T, N> clamp(const vec<T, N>& x, const T& min, const T& max) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = clamp(x[i], min, max);
    return c;
}

template <typename T, size_t N>
inline vec<T, N> clamp(const vec<T, N>& x, const vec<T, N>& min,
                       const vec<T, N>& max) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = clamp(x[i], min[i], max[i]);
    return c;
}

template <typename T, size_t N, typename T1>
inline vec<T, N> clamplen(const vec<T, N> x, T1 max) {
    auto l = length(x);
    return (l > (T)max) ? x * (T)max / l : x;
}

//
// Element min/max.
//
template <typename T, size_t N>
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
template <typename T, size_t N>
inline int max_element(const vec<T, N>& a) {
    auto v = -std::numeric_limits<T>::max();
    auto pos = -1;
    for (auto i = 0; i < N; i++) {
        if (v < a[i]) {
            v = a[i];
            pos = i;
        }
    }
    return pos;
}

// -----------------------------------------------------------------------------
// MATRICES
// -----------------------------------------------------------------------------

//
// Matrix of element of compile time dimensions stored in column major format.
// Allows for default initializer, initalization with columns and to/from
// smaller matrices.
// Access to columns via operator [] and raw data pointer data().
//
template <typename T, size_t N, size_t M>
struct mat : array<vec<T, N>, M> {
    using V = vec<T, N>;

    constexpr mat() : array<vec<T, N>, M>{} {}

    // member constructor
    template <typename... Args,
              typename = std::enable_if_t<sizeof...(Args) == M>>
    constexpr mat(Args... args) : array<vec<T, N>, M>{V(args)...} {}

    // list constructor
    constexpr mat(initializer_list<V> vv) {
        assert(M == vv.size());
        auto i = 0;
        for (auto&& e : vv) (*this)[i++] = e;
    }

    // conversions
    constexpr mat(const array<array<float, N>, M>& v) {
        *this = *(mat<T, N, M>*)&v;
    }

    constexpr operator array<array<T, N>, M>() const {
        return *(array<array<T, N>, M>*)this;
    }

    constexpr mat(const array<float, N * M>& v) { *this = *(mat<T, N, M>*)&v; }

    constexpr operator array<T, N * M>() const {
        return *(array<T, N * M>*)this;
    }
};

//
// Typedef for standard matrices.
//

using mat2f = mat<float, 2, 2>;
using mat3f = mat<float, 3, 3>;
using mat4f = mat<float, 4, 4>;

//
// Matrix Constants.
//

const auto identity_mat2f = mat2f{{1, 0}, {0, 1}};
const auto identity_mat3f = mat3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
const auto identity_mat4f =
    mat4f{{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};

//
// Components-wise arithmetic
//

template <typename T, size_t N, size_t M>
inline mat<T, M, N> operator-(const mat<T, N, M>& a) {
    mat<T, N, M> c;
    for (auto i = 0; i < M; i++) c[i] = -a[i];
    return c;
}

template <typename T, size_t N, size_t M>
inline mat<T, M, N> operator+(const mat<T, N, M>& a, const mat<T, N, M>& b) {
    mat<T, N, M> c;
    for (auto i = 0; i < M; i++) c[i] = a[i] + b[i];
    return c;
}

//
// Vector/Matrix multiplies
//

template <typename T, size_t N, size_t M>
inline mat<T, M, N> operator*(const mat<T, N, M>& a, T b) {
    mat<T, N, M> c;
    for (auto i = 0; i < M; i++) c[i] = a[i] * b;
    return c;
}

template <typename T, size_t N, size_t M>
inline mat<T, M, N> operator/(const mat<T, N, M>& a, T b) {
    mat<T, N, M> c;
    for (auto i = 0; i < M; i++) c[i] = a[i] / b;
    return c;
}

template <typename T, size_t N, size_t M>
inline vec<T, N> operator*(const mat<T, N, M>& a, const vec<T, M>& b) {
    vec<T, N> c;
    for (auto j = 0; j < M; j++) c += a[j] * b[j];
    return c;
}

template <typename T, size_t N, size_t M>
inline vec<T, M> operator*(const vec<T, N>& a, const mat<T, N, M>& b) {
    vec<T, M> c;
    for (auto j = 0; j < M; j++) c[j] = dot(a, b[j]);
    return c;
}

template <typename T, size_t N, size_t M, size_t K>
inline mat<T, N, M> operator*(const mat<T, N, K>& a, const mat<T, K, M>& b) {
    mat<T, N, M> c;
    for (auto j = 0; j < M; j++) c[j] = a * b[j];
    return c;
}

template <typename T, size_t N, size_t M>
inline mat<T, M, N>& operator+=(mat<T, N, M>& a, const mat<T, N, M>& b) {
    return a = a + b;
}

template <typename T, size_t N, size_t M>
inline mat<T, M, N>& operator*=(mat<T, N, M>& a, const mat<T, N, M>& b) {
    return a = a * b;
}

template <typename T, size_t N, size_t M>
inline mat<T, M, N>& operator*=(mat<T, N, M>& a, const T& b) {
    return a = a * b;
}

template <typename T, size_t N, size_t M>
inline mat<T, M, N>& operator/=(mat<T, N, M>& a, const T& b) {
    return a = a / b;
}

//
// Linear Algebra operations.
//

template <typename T, size_t N>
vec<T, N> diagonal(const mat<T, N, N>& a) {
    vec<T, N> d;
    for (auto i = 0; i < N; i++) d[i] = a[i][i];
    return d;
}

template <typename T, size_t N, size_t M>
inline mat<T, M, N> transpose(const mat<T, N, M>& a) {
    mat<T, M, N> c;
    for (auto j = 0; j < M; j++) {
        for (auto i = 0; i < N; i++) {
            c[i][j] = a[j][i];
        }
    }
    return c;
}

template <typename T>
inline mat<T, 2, 2> adjugate(const mat<T, 2, 2>& a) {
    return {{a[1][1], -a[0][1]}, {-a[1][0], a[0][0]}};
}

template <typename T>
inline mat<T, 3, 3> adjugate(const mat<T, 3, 3>& a) {
    return {{a[1][1] * a[2][2] - a[2][1] * a[1][2],
             a[2][1] * a[0][2] - a[0][1] * a[2][2],
             a[0][1] * a[1][2] - a[1][1] * a[0][2]},
            {a[1][2] * a[2][0] - a[2][2] * a[1][0],
             a[2][2] * a[0][0] - a[0][2] * a[2][0],
             a[0][2] * a[1][0] - a[1][2] * a[0][0]},
            {a[1][0] * a[2][1] - a[2][0] * a[1][1],
             a[2][0] * a[0][1] - a[0][0] * a[2][1],
             a[0][0] * a[1][1] - a[1][0] * a[0][1]}};
}

template <typename T>
inline mat<T, 4, 4> adjugate(const mat<T, 4, 4>& a) {
    return {{a[1][1] * a[2][2] * a[3][3] + a[3][1] * a[1][2] * a[2][3] +
                 a[2][1] * a[3][2] * a[1][3] - a[1][1] * a[3][2] * a[2][3] -
                 a[2][1] * a[1][2] * a[3][3] - a[3][1] * a[2][2] * a[1][3],
             a[0][1] * a[3][2] * a[2][3] + a[2][1] * a[0][2] * a[3][3] +
                 a[3][1] * a[2][2] * a[0][3] - a[3][1] * a[0][2] * a[2][3] -
                 a[2][1] * a[3][2] * a[0][3] - a[0][1] * a[2][2] * a[3][3],
             a[0][1] * a[1][2] * a[3][3] + a[3][1] * a[0][2] * a[1][3] +
                 a[1][1] * a[3][2] * a[0][3] - a[0][1] * a[3][2] * a[1][3] -
                 a[1][1] * a[0][2] * a[3][3] - a[3][1] * a[1][2] * a[0][3],
             a[0][1] * a[2][2] * a[1][3] + a[1][1] * a[0][2] * a[2][3] +
                 a[2][1] * a[1][2] * a[0][3] - a[0][1] * a[1][2] * a[2][3] -
                 a[2][1] * a[0][2] * a[1][3] - a[1][1] * a[2][2] * a[0][3]},
            {a[1][2] * a[3][3] * a[2][0] + a[2][2] * a[1][3] * a[3][0] +
                 a[3][2] * a[2][3] * a[1][0] - a[1][2] * a[2][3] * a[3][0] -
                 a[3][2] * a[1][3] * a[2][0] - a[2][2] * a[3][3] * a[1][0],
             a[0][2] * a[2][3] * a[3][0] + a[3][2] * a[0][3] * a[2][0] +
                 a[2][2] * a[3][3] * a[0][0] - a[0][2] * a[3][3] * a[2][0] -
                 a[2][2] * a[0][3] * a[3][0] - a[3][2] * a[2][3] * a[0][0],
             a[0][2] * a[3][3] * a[1][0] + a[1][2] * a[0][3] * a[3][0] +
                 a[3][2] * a[1][3] * a[0][0] - a[0][2] * a[1][3] * a[3][0] -
                 a[3][2] * a[0][3] * a[1][0] - a[1][2] * a[3][3] * a[0][0],
             a[0][2] * a[1][3] * a[2][0] + a[2][2] * a[0][3] * a[1][0] +
                 a[1][2] * a[2][3] * a[0][0] - a[0][2] * a[2][3] * a[1][0] -
                 a[1][2] * a[0][3] * a[2][0] - a[2][2] * a[1][3] * a[0][0]},
            {a[1][3] * a[2][0] * a[3][1] + a[3][3] * a[1][0] * a[2][1] +
                 a[2][3] * a[3][0] * a[1][1] - a[1][3] * a[3][0] * a[2][1] -
                 a[2][3] * a[1][0] * a[3][1] - a[3][3] * a[2][0] * a[1][1],
             a[0][3] * a[3][0] * a[2][1] + a[2][3] * a[0][0] * a[3][1] +
                 a[3][3] * a[2][0] * a[0][1] - a[0][3] * a[2][0] * a[3][1] -
                 a[3][3] * a[0][0] * a[2][1] - a[2][3] * a[3][0] * a[0][1],
             a[0][3] * a[1][0] * a[3][1] + a[3][3] * a[0][0] * a[1][1] +
                 a[1][3] * a[3][0] * a[0][1] - a[0][3] * a[3][0] * a[1][1] -
                 a[1][3] * a[0][0] * a[3][1] - a[3][3] * a[1][0] * a[0][1],
             a[0][3] * a[2][0] * a[1][1] + a[1][3] * a[0][0] * a[2][1] +
                 a[2][3] * a[1][0] * a[0][1] - a[0][3] * a[1][0] * a[2][1] -
                 a[2][3] * a[0][0] * a[1][1] - a[1][3] * a[2][0] * a[0][1]},
            {a[1][0] * a[3][1] * a[2][2] + a[2][0] * a[1][1] * a[3][2] +
                 a[3][0] * a[2][1] * a[1][2] - a[1][0] * a[2][1] * a[3][2] -
                 a[3][0] * a[1][1] * a[2][2] - a[2][0] * a[3][1] * a[1][2],
             a[0][0] * a[2][1] * a[3][2] + a[3][0] * a[0][1] * a[2][2] +
                 a[2][0] * a[3][1] * a[0][2] - a[0][0] * a[3][1] * a[2][2] -
                 a[2][0] * a[0][1] * a[3][2] - a[3][0] * a[2][1] * a[0][2],
             a[0][0] * a[3][1] * a[1][2] + a[1][0] * a[0][1] * a[3][2] +
                 a[3][0] * a[1][1] * a[0][2] - a[0][0] * a[1][1] * a[3][2] -
                 a[3][0] * a[0][1] * a[1][2] - a[1][0] * a[3][1] * a[0][2],
             a[0][0] * a[1][1] * a[2][2] + a[2][0] * a[0][1] * a[1][2] +
                 a[1][0] * a[2][1] * a[0][2] - a[0][0] * a[2][1] * a[1][2] -
                 a[1][0] * a[0][1] * a[2][2] - a[2][0] * a[1][1] * a[0][2]}};
}

template <typename T>
inline T determinant(const mat<T, 2, 2>& a) {
    return a[0][0] * a[1][1] - a[0][1] * a[1][0];
}

template <typename T>
inline T determinant(const mat<T, 3, 3>& a) {
    return a[0][0] * (a[1][1] * a[2][2] - a[2][1] * a[1][2]) +
           a[0][1] * (a[1][2] * a[2][0] - a[2][2] * a[1][0]) +
           a[0][2] * (a[1][0] * a[2][1] - a[2][0] * a[1][1]);
}

template <typename T>
inline T determinant(const mat<T, 4, 4>& a) {
    return a[0][0] *
               (a[1][1] * a[2][2] * a[3][3] + a[3][1] * a[1][2] * a[2][3] +
                a[2][1] * a[3][2] * a[1][3] - a[1][1] * a[3][2] * a[2][3] -
                a[2][1] * a[1][2] * a[3][3] - a[3][1] * a[2][2] * a[1][3]) +
           a[0][1] *
               (a[1][2] * a[3][3] * a[2][0] + a[2][2] * a[1][3] * a[3][0] +
                a[3][2] * a[2][3] * a[1][0] - a[1][2] * a[2][3] * a[3][0] -
                a[3][2] * a[1][3] * a[2][0] - a[2][2] * a[3][3] * a[1][0]) +
           a[0][2] *
               (a[1][3] * a[2][0] * a[3][1] + a[3][3] * a[1][0] * a[2][1] +
                a[2][3] * a[3][0] * a[1][1] - a[1][3] * a[3][0] * a[2][1] -
                a[2][3] * a[1][0] * a[3][1] - a[3][3] * a[2][0] * a[1][1]) +
           a[0][3] *
               (a[1][0] * a[3][1] * a[2][2] + a[2][0] * a[1][1] * a[3][2] +
                a[3][0] * a[2][1] * a[1][2] - a[1][0] * a[2][1] * a[3][2] -
                a[3][0] * a[1][1] * a[2][2] - a[2][0] * a[3][1] * a[1][2]);
}

template <typename T, size_t N>
inline mat<T, N, N> inverse(const mat<T, N, N>& a) {
    return adjugate(a) / determinant(a);
}

// -----------------------------------------------------------------------------
// RIGID TRANSFORMS
// -----------------------------------------------------------------------------

//
// Rigid transforms stored as a linear tranform NxN-matrix m and a translation
// N-vector t. This is a special case of an affine matrix.
// Provides converion to/from (N+1)x(N+1) matrices. The matrix
// columns are the axis of the coordinate system, while the translation is
// its origin.
//
template <typename T, size_t N>
struct frame : array<vec<T, N>, N + 1> {
    using V = vec<T, N>;
    using M = mat<T, N, N>;

    constexpr frame() : array<vec<T, N>, N + 1>{} {}

    constexpr frame(const M& m, const V& t) {
        for (auto i = 0; i < N; i++) (*this)[i] = m[i];
        (*this)[N] = t;
    }

    // member constructor
    template <typename... Args,
              typename = std::enable_if_t<sizeof...(Args) == N + 1>>
    constexpr frame(Args... args) : array<vec<T, N>, N + 1>{V(args)...} {}

    // list constructor
    constexpr frame(initializer_list<V> vv) {
        assert(N + 1 == vv.size());
        auto i = 0;
        for (auto&& e : vv) (*this)[i++] = e;
    }

    // view as linear + translation
    constexpr const V& o() const { return (*this)[N]; }
    constexpr V& o() { return (*this)[N]; }
    constexpr const M& m() const { return *(M*)this; }
    constexpr M& m() { return *(M*)this; }
};

//
// Typedef for standard rigid transforms / frames.
//

using frame2f = frame<float, 2>;
using frame3f = frame<float, 3>;

//
// Rigid Transforms/Frames Constants.
//

const auto identity_frame2f = frame2f{{1, 0}, {0, 1}, {0, 0}};
const auto identity_frame3f =
    frame3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};

//
// Conversions
//
template <typename T, size_t N>
inline mat<T, N + 1, N + 1> to_mat(const frame<T, N>& a) {
    auto m = mat<T, N + 1, N + 1>();
    for (auto j = 0; j < N; j++) {
        (vec<T, N>&)m[j] = a[j];
        m[j][N] = 0;
    }
    (vec<T, N>&)m[N] = a[N];
    m[N][N] = 1;
    return m;
}

template <typename T, size_t N>
inline frame<T, N - 1> to_frame(const mat<T, N, N>& a) {
    auto f = frame<T, N - 1>();
    for (auto j = 0; j < N; j++) {
        for (auto i = 0; i < N - 1; i++) {
            f[j][i] = a[j][i];
        }
    }
    return f;
}

//
// Matrix operations for transforms.
//

template <typename T, size_t N>
inline frame<T, N> operator*(const frame<T, N>& a, const frame<T, N>& b) {
    return {a.m() * b.m(), a.m() * b.o() + a.o()};
}

template <typename T, size_t N>
inline frame<T, N> inverse(const frame<T, N>& a) {
    auto minv = transpose(a.m());
    return {minv, -(minv * a.o())};
}

// -----------------------------------------------------------------------------
// QUATERNIONS
// -----------------------------------------------------------------------------

//
// Quaternions implemented as a vec4. Data access via operator[].
// Quaterions are xi + yj + zk + w
//
template <typename T, size_t N>
struct quat : array<T, N> {
    // default constructor
    constexpr quat() : array<T, 4>{} {}

    // member constructor
    template <typename... Args,
              typename = std::enable_if_t<sizeof...(Args) == N>>
    constexpr quat(Args... vv) : array<T, 4>{T(vv)...} {}

    // list constructor
    constexpr quat(initializer_list<T> vv) {
        assert(N == vv.size());
        auto i = 0;
        for (auto&& e : vv) (*this)[i++] = e;
    }

    // axis angle access
    constexpr T angle() const { return acos((*this)[2]) * 2; }
    constexpr vec<T, 3> axis() {
        return normalize(vec<T, 3>((*this)[0], (*this)[1], (*this)[2]));
    }
};

//
// Conversions
//
template <typename T>
inline mat<T, 4, 4> to_mat(const quat<T, 4>& v) {
    return {{v[3] * v[3] + v[0] * v[0] - v[1] * v[1] - v[2] * v[2],
             (v[0] * v[1] + v[2] * v[3]) * 2, (v[2] * v[0] - v[1] * v[3]) * 2},
            {(v[0] * v[1] - v[2] * v[3]) * 2,
             v[3] * v[3] - v[0] * v[0] + v[1] * v[1] - v[2] * v[2],
             (v[1] * v[2] + v[0] * v[3]) * 2},
            {(v[2] * v[0] + v[1] * v[3]) * 2, (v[1] * v[2] - v[0] * v[3]) * 2,
             v[3] * v[3] - v[0] * v[0] - v[1] * v[1] + v[2] * v[2]}};
}

//
// Typedef for quaterions.
//

using quat4f = quat<float, 4>;

//
// Support for quaternion algebra using 4D vectors.
//

template <typename T>
constexpr quat<T, 4> conjugate(const quat<T, 4>& v) {
    return {-v[0], -v[1], -v[2], v[3]};
}

template <typename T>
quat<T, 4> inverse(const quat<T, 4>& v) {
    return qconj(v) / lengthsqr(vec<T, 4>(v));
}

template <typename T>
constexpr quat<T, 4> operator*(const quat<T, 4>& a, const quat<T, 4>& b) {
    return {a[0] * b[3] + a[3] * b[0] + a[1] * b[3] - a[2] * b[1],
            a[1] * b[3] + a[3] * b[1] + a[2] * b[0] - a[0] * b[2],
            a[2] * b[3] + a[3] * b[2] + a[0] * b[1] - a[1] * b[0],
            a[3] * b[3] - a[0] * b[0] - a[1] * b[1] - a[2] * b[2]};
}

template <typename T>
quat<T, 4> nlerp(const quat<T, 4>& a, const quat<T, 4>& b, T t) {
    return nlerp(
        vec<T, 4>(a),
        dot(vec<T, 4>(a), vec<T, 4>(b)) < 0 ? -vec<T, 4>(b) : vec<T, 4>(b), t);
}

template <typename T>
quat<T, 4> slerp(const quat<T, 4>& a, const quat<T, 4>& b, T t) {
    return slerp(
        vec<T, 4>(a),
        dot(vec<T, 4>(a), vec<T, 4>(b)) < 0 ? -vec<T, 4>(b) : vec<T, 4>(b), t);
}

// -----------------------------------------------------------------------------
// AXIS ALIGNED BOUNDING BOXES
// -----------------------------------------------------------------------------

//
// Axis aligned bounding box.
//
template <typename T, size_t N>
struct bbox : array<vec<T, N>, 2> {
    constexpr bbox() {
        for (auto i = 0; i < N; i++) {
            (*this)[0][i] = std::numeric_limits<T>::max();
            (*this)[1][i] = std::numeric_limits<T>::lowest();
        }
    }

    constexpr bbox(const vec<T, N>& m, const vec<T, N>& M)
        : array<vec<T, N>, 2>{m, M} {}

    constexpr bbox(const initializer_list<vec<T, N>>& v) : bbox() {
        for (auto&& vv : v) expand(vv);
    }

    constexpr vec<T, N> diagonal() const { return (*this)[1] - (*this)[0]; }
    constexpr vec<T, N> center() const { return ((*this)[1] + (*this)[0]) / 2; }

    constexpr void expand(const vec<T, N>& b) {
        (*this)[0] = min((*this)[0], b);
        (*this)[1] = max((*this)[1], b);
    }
};

//
// Axis aligned bounding box usings.
//

using bbox2f = bbox<float, 2>;
using bbox3f = bbox<float, 3>;

//
// Axis aligned bounding box constants.
//

const auto invalid_bbox2f = bbox2f();
const auto invalid_bbox3f = bbox3f();

//
// Axis aligned bounding box operations.
//

template <typename T, size_t N>
inline bbox<T, N> expand(const bbox<T, N>& a, const vec<T, N>& b) {
    return {min(a[0], b), max(a[1], b)};
}

template <typename T, size_t N>
inline bbox<T, N> expand(const bbox<T, N>& a, const bbox<T, N>& b) {
    return {min(a[0], b[0]), max(a[1], b[1])};
}

template <typename T, size_t N>
inline bbox<T, N> operator+(const bbox<T, N>& a, const T& b) {
    return expand(a, b);
}

template <typename T, size_t N>
inline bbox<T, N> operator+(const bbox<T, N>& a, const bbox<T, N>& b) {
    return expand(a, b);
}

template <typename T, size_t N>
inline bbox<T, N>& operator+=(bbox<T, N>& a, const vec<T, N>& b) {
    return a = expand(a, b);
}

template <typename T, size_t N>
inline bbox<T, N>& operator+=(bbox<T, N>& a, const bbox<T, N>& b) {
    return a = expand(a, b);
}

template <typename T, size_t N>
inline vec<T, N> center(const bbox<T, N>& a) {
    return (a[0] + a[1]) / 2;
}

template <typename T, size_t N>
inline vec<T, N> diagonal(const bbox<T, N>& a) {
    return a[1] - a[0];
}

// -----------------------------------------------------------------------------
// RAYS
// -----------------------------------------------------------------------------

//
// Rays with origin, direction and min/max t value.
//
template <typename T, size_t N>
struct ray {
    vec<T, N> o;
    vec<T, N> d;
    T tmin;
    T tmax;

    ray() : o(), d(0, 0, 1), tmin(0), tmax(std::numeric_limits<T>::max()) {}
    ray(const vec<T, N>& o, const vec<T, N>& d, T tmin = 0,
        T tmax = std::numeric_limits<T>::max())
        : o(o), d(d), tmin(tmin), tmax(tmax) {}

    vec<T, N> eval(T t) const { return o + t * d; }
};

//
// Rays usings.
//

using ray2f = ray<float, 2>;
using ray3f = ray<float, 3>;

// -----------------------------------------------------------------------------
// TRANSFORMS
// -----------------------------------------------------------------------------

//
// Tranform operations.
//

template <typename T, size_t N>
inline vec<T, N> transform_point(const mat<T, N + 1, N + 1>& a,
                                 const vec<T, N>& b) {
    // make it generic
    auto vb = vec<T, N + 1>();
    (vec<T, N>&)vb = b;
    vb[N] = 1;
    auto tvb = a * vb;
    return *(vec<T, N>*)(&tvb) / tvb[N];
}

template <typename T, size_t N>
inline vec<T, N> transform_vector(const mat<T, N + 1, N + 1>& a,
                                  const vec<T, N>& b) {
    // make it generic
    auto vb = vec<T, N + 1>();
    (vec<T, N>&)vb = b;
    vb[N] = 0;
    auto tvb = a * vb;
    return *(vec<T, N>*)(&tvb);
}

template <typename T, size_t N>
inline vec<T, N> transform_direction(const mat<T, N + 1, N + 1>& a,
                                     const vec<T, N>& b) {
    return normalize(transform_vector(a, b));
}

template <typename T, size_t N>
inline vec<T, N> transform_point(const frame<T, N>& a, const vec<T, N>& b) {
    return a.m() * b + a.o();
}

template <typename T, size_t N>
inline vec<T, N> transform_vector(const frame<T, N>& a, const vec<T, N>& b) {
    return a.m() * b;
}

template <typename T, size_t N>
inline vec<T, N> transform_direction(const frame<T, N>& a, const vec<T, N>& b) {
    return a.m() * b;
}

template <typename T, size_t N>
inline frame<T, N> transform_frame(const frame<T, N>& a, const frame<T, N>& b) {
    return {a.m() * b.m(), a.m() * b.o() + a.o()};
}

template <typename T, size_t N>
inline ray<T, N> transform_ray(const frame<T, N>& a, const ray<T, N>& b) {
    return {transform_point(a, b.o), transform_direction(a, b.d), b.tmin,
            b.tmax};
}

template <typename T, size_t N>
inline vec<T, N> transform_point_inverse(const frame<T, N>& a,
                                         const vec<T, N>& b) {
    return (b - a.o()) * a.m();
}

template <typename T, size_t N>
inline vec<T, N> transform_vector_inverse(const frame<T, N>& a,
                                          const vec<T, N>& b) {
    return b * a.m();
}

template <typename T, size_t N>
inline vec<T, N> transform_direction_inverse(const frame<T, N>& a,
                                             const vec<T, N>& b) {
    return b * a.m();
}

template <typename T, size_t N>
inline ray<T, N> transform_ray_inverse(const frame<T, N>& a,
                                       const ray<T, N>& b) {
    return {transform_point_inverse(a, b.o),
            transform_direction_inverse(a, b.d), b.tmin, b.tmax};
}

template <typename T, size_t N>
inline ray<T, N> transform_ray(const mat<T, N + 1, N + 1>& a,
                               const ray<T, N>& b) {
    return {transform_point(a, b.o), transform_direction(a, b.d), b.tmin,
            b.tmax};
}

#if 0
template <typename T>
inline bbox<T, 3> transform_bbox(const frame<T, 3>& a, const bbox<T, 3>& b) {
    vec<T, 3> corners[8] = {
        {b[0][0], b[0][1], b[0][2]}, {b[0][0], b[0][1], b[1][2]},
        {b[0][0], b[1][1], b[0][2]}, {b[0][0], b[1][1], b[1][2]},
        {b[1][0], b[0][1], b[0][2]}, {b[1][0], b[0][1], b[1][2]},
        {b[1][0], b[1][1], b[0][2]}, {b[1][0], b[1][1], b[1][2]},
    };
    auto xformed = bbox<T, 3>();
    for (auto j = 0; j < 8; j++) xformed += transform_point(a, corners[j]);
    return xformed;
}
#else
// Code from Real-time Collision Detection by Christer Ericson Sect. 4.2.6
// Transform AABB a by the matrix m and translation t,
// find maximum extents, and store result into AABB b.
template <typename T>
inline bbox<T, 3> transform_bbox(const frame<T, 3>& a, const bbox<T, 3>& b) {
    // start by adding in translation
    auto c = bbox<T, 3>(a.o(), a.o());
    // compute extent
    auto e = a.m() * b[0], f = a.m() * b[1];
    // for all three axes
    for (auto i = 0; i < 3; i++) {
        // form extent by summing smaller and larger terms respectively
        for (auto j = 0; j < 3; j++) {
            auto e = a.m()[j][i] * b[0][j];
            auto f = a.m()[j][i] * b[1][j];
            if (e < f) {
                c[0][i] += e;
                c[1][i] += f;
            } else {
                c[0][i] += f;
                c[1][i] += e;
            }
        }
    }
    return c;
}
#endif

template <typename T>
inline bbox<T, 3> transform_bbox_inverse(const frame<T, 3>& a,
                                         const bbox<T, 3>& b) {
    return transform_bbox(inverse(a), b);
}

template <typename T>
inline bbox<T, 3> transform_bbox(const mat<T, 4, 4>& a, const bbox<T, 3>& b) {
    vec<T, 3> corners[8] = {
        {b[0][0], b[0][1], b[0][2]}, {b[0][0], b[0][1], b[1][2]},
        {b[0][0], b[1][1], b[0][2]}, {b[0][0], b[1][1], b[1][2]},
        {b[1][0], b[0][1], b[0][2]}, {b[1][0], b[0][1], b[1][2]},
        {b[1][0], b[1][1], b[0][2]}, {b[1][0], b[1][1], b[1][2]},
    };
    auto xformed = bbox<T, 3>();
    for (auto j = 0; j < 8; j++) xformed += transform_point(a, corners[j]);
    return xformed;
}

//
// Rotation matrix crom axis/angle
//
template <typename T>
inline mat<T, 3, 3> rotation_mat3(const vec<T, 3>& axis, T angle) {
    auto s = sin(angle), c = cos(angle);
    auto vv = normalize(axis);
    return {{c + (1 - c) * vv[0] * vv[0], (1 - c) * vv[0] * vv[1] + s * vv[2],
             (1 - c) * vv[0] * vv[2] - s * vv[1]},
            {(1 - c) * vv[0] * vv[1] - s * vv[2], c + (1 - c) * vv[1] * vv[1],
             (1 - c) * vv[1] * vv[2] + s * vv[0]},
            {(1 - c) * vv[0] * vv[2] + s * vv[1],
             (1 - c) * vv[1] * vv[2] - s * vv[0], c + (1 - c) * vv[2] * vv[2]}};
}

//
// Translation transform
//
template <typename T>
inline frame<T, 3> translation_frame3(const vec<T, 3>& a) {
    return {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, a};
}

template <typename T>
inline mat<T, 4, 4> translation_mat4(const vec<T, 3>& a) {
    return (mat<T, 4, 4>)translation_frame3(a);
}

//
// Scaling transform (in this case the frame is broken and used only
// as
// affine)
//
template <typename T>
inline frame<T, 3> scaling_frame3(const vec<T, 3>& a) {
    return {{a[0], 0, 0}, {0, a[1], 0}, {0, 0, a[2]}, {0, 0, 0}};
}

template <typename T>
inline mat<T, 4, 4> scaling_mat4(const vec<T, 3>& a) {
    return (mat<T, 4, 4>)scaling_frame3(a);
}

//
// Rotation transform
//
template <typename T>
inline frame<T, 3> rotation_frame3(const vec<T, 3>& axis, T angle) {
    return frame<T, 3>{rotation_mat3(axis, angle), {0, 0, 0}};
}

template <typename T>
inline mat<T, 4, 4> rotation_mat4(const vec<T, 3>& axis, T angle) {
    return (mat<T, 4, 4>)rotation_frame3(axis, angle);
}

//
// Lookat tranform
//
template <typename T>
inline frame<T, 3> lookat_frame3(const vec<T, 3>& eye, const vec<T, 3>& center,
                                 const vec<T, 3>& up) {
    auto w = normalize(eye - center);
    auto u = normalize(cross(up, w));
    auto v = normalize(cross(w, u));
    return {u, v, w, eye};
}

template <typename T>
inline mat<T, 4, 4> lookat_mat4(const vec<T, 3>& eye, const vec<T, 3>& center,
                                const vec<T, 3>& up) {
    return to_mat(lookat_frame3(eye, center, up));
}

//
// Frame/transform for origin and z.
//
template <typename T>
inline frame<T, 3> make_frame3(const vec<T, 3>& o, const vec<T, 3>& z_) {
    auto z = normalize(z_);
    auto x = normalize(orthogonal(z));
    auto y = normalize(cross(z, x));
    return {x, y, z, o};
}

//
// OpenGL perspective frustum matrix
//
template <typename T>
inline mat<T, 4, 4> frustum_mat4(T l, T r, T b, T t, T n, T f) {
    return {{2 * n / (r - l), 0, 0, 0},
            {0, 2 * n / (t - b), 0, 0},
            {(r + l) / (r - l), (t + b) / (t - b), -(f + n) / (f - n), -1},
            {0, 0, -2 * f * n / (f - n), 0}};
}

//
// OpenGL orthographic matrix
//
template <typename T>
inline mat<T, 4, 4> ortho_mat4(T l, T r, T b, T t, T n, T f) {
    return {{2 / (r - l), 0, 0, 0},
            {0, 2 / (t - b), 0, 0},
            {0, 0, -2 / (f - n), 0},
            {-(r + l) / (r - l), -(t + b) / (t - b), -2 / (f - n),
             -(f + n) / (f - n), 1}};
}

//
// OpenGL orthographic 2D matrix
//
template <typename T>
inline mat<T, 4, 4> ortho2d_mat4(T l, T r, T b, T t) {
    return ortho_mat4(l, r, b, t, -1, 1);
}

//
// OpenGL perspective matrix
//
template <typename T>
inline mat<T, 4, 4> perspective_mat4(T fovy, T aspect, T near, T far) {
    auto y = near * tan(fovy / 2);
    auto x = y * aspect;
    return frustum_mat4<T>(-x, x, -y, y, near, far);
}

// -----------------------------------------------------------------------------
// GEOMETRY UTILITIES
// -----------------------------------------------------------------------------

template <typename T>
inline vec<T, 3> triangle_normal(const vec<T, 3>& v0, const vec<T, 3>& v1,
                                 const vec<T, 3>& v2) {
    return normalize(cross(v1 - v0, v2 - v0));
}

template <typename T>
inline T triangle_area(const vec<T, 3>& v0, const vec<T, 3>& v1,
                       const vec<T, 3>& v2) {
    return length(cross(v1 - v0, v2 - v0)) / 2;
}

template <typename T>
inline T tetrahedron_volume(const vec<T, 3>& v0, const vec<T, 3>& v1,
                            const vec<T, 3>& v2, const vec<T, 3>& v3) {
    return dot(cross(v1 - v0, v2 - v0), v3 - v0) / 6;
}

//
// Baricentric interpolation
//
template <typename T, typename T1>
inline T blerp(const T& a, const T& b, const T& c, const T1& w) {
    return a * w[0] + b * w[1] + c * w[2];
}

// -----------------------------------------------------------------------------
// UI UTILITIES
// -----------------------------------------------------------------------------

//
// Turntable for UI navigation from a from/to/up parametrization of
// the
// camera.
//
template <typename T>
inline void turntable(vec<T, 3>& from, vec<T, 3>& to, vec<T, 3>& up,
                      const vec<T, 2>& rotate, T dolly, const vec<T, 2>& pan) {
    // rotate if necessary
    if (rotate[0] || rotate[1]) {
        auto z = ym_normalize(*to - *from);
        auto lz = ym_dist(*to, *from);
        auto phi = atan2(z[2], z[0]) + rotate[0];
        auto theta = acos(z[1]) + rotate[1];
        theta = max(T(0.001), min(theta, T(pi - 0.001)));
        auto nz = vec<T, 3>{sin(theta) * cos(phi) * lz, cos(theta) * lz,
                            sin(theta) * sin(phi) * lz};
        *from = *to - nz;
    }

    // dolly if necessary
    if (dolly) {
        auto z = ym_normalize(*to - *from);
        auto lz = ym_max(T(0.001), ym_dist(*to, *from) * (1 + dolly));
        z *= lz;
        *from = *to - z;
    }

    // pan if necessary
    if (pan[0] || pan[1]) {
        auto z = ym_normalize(*to - *from);
        auto x = ym_normalize(ym_cross(*up, z));
        auto y = ym_normalize(ym_cross(z, x));
        auto t = vec<T, 3>{pan[0] * x[0] + pan[1] * y[0],
                           pan[0] * x[1] + pan[1] * y[1],
                           pan[0] * x[2] + pan[1] * y[2]};
        *from += t;
        *to += t;
    }
}

//
// Turntable for UI navigation for a frame/distance parametrization
// of the
// camera.
//
template <typename T>
inline void turntable(frame<T, 3>& frame, float& focus, const vec<T, 2>& rotate,
                      T dolly, const vec<T, 2>& pan) {
    // rotate if necessary
    if (rotate[0] || rotate[1]) {
        auto phi = atan2(frame[2][2], frame[2][0]) + rotate[0];
        auto theta = acos(frame[2][1]) + rotate[1];
        theta = max(T(0.001), min(theta, T(pi - 0.001)));
        auto new_z =
            vec<T, 3>(sin(theta) * cos(phi), cos(theta), sin(theta) * sin(phi));
        auto new_center = frame.o() - frame[2] * focus;
        auto new_o = new_center + new_z * focus;
        frame = lookat_frame3(new_o, new_center, {0, 1, 0});
        focus = dist(new_o, new_center);
    }

    // pan if necessary
    if (dolly) {
        auto c = frame.o() - frame[2] * focus;
        focus = max(focus + dolly, T(0.001));
        frame.o() = c + frame[2] * focus;
    }

    // pan if necessary
    if (pan[0] || pan[1]) {
        frame.o() += frame[0] * pan[0] + frame[1] * pan[1];
    }
}

// -----------------------------------------------------------------------------
// RANDOM NUMBER GENERATION
// -----------------------------------------------------------------------------

//
// PCG random numbers. A family of random number generators that
// supports
// multiple sequences. In our code, we allocate one sequence for
// each
// sample.
// PCG32 from http://www.pcg-random.org/
//

//
// Random number state (PCG32)
//
struct rng_pcg32 {
    uint64_t state, inc;
};

//
// Next random number
//
inline uint32_t rng_next(rng_pcg32& rng) {
    uint64_t oldstate = rng.state;
    rng.state = oldstate * 6364136223846793005ull + (rng.inc | 1u);
    uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-((int32_t)rot)) & 31));
}

//
// Init a random number generator with a state state from the
// sequence seq.
//
inline void rng_init(rng_pcg32& rng, uint64_t state, uint64_t seq) {
    rng.state = 0U;
    rng.inc = (seq << 1u) | 1u;
    rng_next(rng);
    rng.state += state;
    rng_next(rng);
}

//
// Next random float in [0,1).
//
inline float rng_nextf(rng_pcg32& rng) {
    return (float)ldexp(rng_next(rng), -32);
}

//
// Next random float in [0,1)x[0,1).
//
inline vec2f rng_next2f(rng_pcg32& rng) {
    return {rng_nextf(rng), rng_nextf(rng)};
}

// -----------------------------------------------------------------------------
// HASHING
// -----------------------------------------------------------------------------

//
// Computes the i-th term of a permutation of l values keyed by p.
// From Correlated Multi-Jittered Sampling by Kensler @ Pixar
//
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

//
// Computes a float value by hashing i with a key p.
// From Correlated Multi-Jittered Sampling by Kensler @ Pixar
//
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

//
// 64 bit integer hash. Public domain code.
//
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

//
// 64-to-32 bit integer hash. Public domain code.
//
inline uint32_t hash_uint64_32(uint64_t a) {
    a = (~a) + (a << 18);  // a = (a << 18) - a - 1;
    a ^= (a >> 31);
    a *= 21;  // a = (a + (a << 2)) + (a << 4);
    a ^= (a >> 11);
    a += (a << 6);
    a ^= (a >> 22);
    return (uint32_t)a;
}

//
// A function to combine 64 bit hashes with semantics as
// boost::hash_combine
//
inline size_t hash_combine(size_t a, size_t b) {
    return a ^ (b + 0x9e3779b9 + (a << 6) + (a >> 2));
}

//
// Hash a vector for std
//
template <typename T, size_t N>
inline size_t hash_vec(const vec<T, N>& v) {
    std::hash<T> Th;
    size_t h = 0;
    for (auto i = 0; i < N; i++) {
        h ^= (Th(v[i]) + 0x9e3779b9 + (h << 6) + (h >> 2));
    }
    return h;
}

// -----------------------------------------------------------------------------
// VIEW CONTAINERS
// -----------------------------------------------------------------------------

//
// An array_view is a non-owining reference to an array with an API
// similar
// to a vector/array containers, but without reallocation.
// This is inspired, but significantly simpler than
// gsl::span https://github.com/Microsoft/GSL or array_view.
//
template <typename T>
struct array_view {
    // constructors
    constexpr array_view() noexcept : _num(0), _data(nullptr) {}
    constexpr array_view(size_t num, T* data) noexcept
        : _num(num), _data(data) {}
    constexpr array_view(T* begin, T* end) noexcept
        : _num((end - begin)), _data(begin) {}
    constexpr array_view(const array_view<std::add_const_t<T>>& av)
        : _num(av._num), _data(av._data) {}
    constexpr array_view(vector<T>& av) : _num(av.size()), _data(av.data()) {}
    constexpr array_view(const vector<T>& av)
        : _num(av.size()), _data((T*)av.data()) {}

    // size
    constexpr size_t size() const noexcept { return _num; }
    constexpr bool empty() const noexcept { return _num == 0; }

    // raw data access
    constexpr T* data() noexcept { return _data; }
    constexpr const T* data() const noexcept { return _data; }

    // iterators
    constexpr T* begin() noexcept { return _data; }
    constexpr T* end() noexcept { return _data + _num; }
    constexpr const T* begin() const noexcept { return _data; }
    constexpr const T* end() const noexcept { return _data + _num; }

    // elements access
    constexpr T& operator[](int i) noexcept { return _data[i]; }
    constexpr const T& operator[](int i) const noexcept { return _data[i]; }

    constexpr T& at(int i) { return _data[i]; }
    constexpr const T& at(int i) const { return _data[i]; }

    constexpr T& front() { return _data[0]; }
    constexpr const T& front() const { return _data[0]; }
    constexpr T& back() { return _data[_num - 1]; }
    constexpr const T& back() const { return _data[_num - 1]; }

   private:
    size_t _num;
    T* _data;
};

// -----------------------------------------------------------------------------
// IMAGES
// -----------------------------------------------------------------------------

//
// An image_view is a non-owining reference to an array that allows
// to
// access it
// as a row-major image, but without reallocation.
// This is inspired, but significantly simpler than
// gsl::multi_span https://github.com/Microsoft/GSL or array_view.
//
template <typename T>
struct image_view {
    // constructors
    constexpr image_view() noexcept : _size(0, 0), _data(nullptr) {}
    constexpr image_view(const vec2i& size, T* data) noexcept
        : _size(size), _data(data) {}

    // size
    constexpr vec2i size() const noexcept { return _size; }
    constexpr bool empty() const noexcept {
        return _size[0] == 0 || _size[1] == 0;
    }

    // raw data access
    constexpr T* data() noexcept { return _data; }
    constexpr const T* data() const noexcept { return _data; }

    // iterators
    constexpr T* begin() noexcept { return _data; }
    constexpr T* end() noexcept { return _data + _size[0] * _size[1]; }
    constexpr const T* begin() const noexcept { return _data; }
    constexpr const T* end() const noexcept {
        return _data + _size[0] * _size[1];
    }

    // elements access
    constexpr T& operator[](const vec<int, 2>& ij) noexcept {
        return _data[ij[1] * _size[0] + ij[0]];
    }
    constexpr const T& operator[](const vec<int, 2>& ij) const noexcept {
        return _data[ij[1] * _size[0] + ij[0]];
    }

    constexpr T& at(const vec<int, 2>& ij) {
        return _data[ij[1] * _size[0] + ij[0]];
    }
    constexpr const T& at(const vec<int, 2>& ij) const {
        return _data[ij[1] * _size[0] + ij[0]];
    }

   private:
    vec2i _size;
    T* _data;
};

//
// Generic image with pixels stored as a vector. For operations use
// image_views.
//
template <typename T>
struct image {
    // constructors
    constexpr image() : _size(0, 0), _data() {}
    constexpr image(const vec2i& size, const T& v = T()) noexcept
        : _size(size), _data(size[0] * size[1], v) {}
    constexpr image(const vec2i& size, const T* v) noexcept
        : _size(size), _data(v, v + size[0] * size[1]) {}

    // access via image_view
    constexpr operator const image_view<T>() const {
        return image_view<T>(_size, const_cast<T*>(_data.data()));
    }
    constexpr operator image_view<T>() {
        return image_view<T>(_size, _data.data());
    }

    // size
    constexpr vec2i size() const noexcept { return _size; }
    constexpr bool empty() const noexcept {
        return _size[0] == 0 || _size[0] == 0;
    }

    // modify
    constexpr void resize(const vec2i& size) {
        _size = size;
        _data.resize(size[0] * size[1]);
    }

    constexpr void clear() {
        _size = {0, 0};
        _data.clear();
    }

    // raw data access
    constexpr T* data() noexcept { return _data.data(); }
    constexpr const T* data() const noexcept { return _data.data(); }

    // iterators
    constexpr T* begin() noexcept { return _data; }
    constexpr T* end() noexcept { return _data + _size[0] * _size[1]; }
    constexpr const T* begin() const noexcept { return _data; }
    constexpr const T* end() const noexcept {
        return _data + _size[0] * _size[1];
    }

    // elements access
    constexpr T& operator[](const vec<int, 2>& ij) noexcept {
        return _data[ij[1] * _size[0] + ij[0]];
    }
    constexpr const T& operator[](const vec<int, 2>& ij) const noexcept {
        return _data[ij[1] * _size[0] + ij[0]];
    }

    constexpr T& at(const vec<int, 2>& ij) {
        return _data[ij[1] * _size[0] + ij[0]];
    }
    constexpr const T& at(const vec<int, 2>& ij) const {
        return _data[ij[1] * _size[0] + ij[0]];
    }

   private:
    vec2i _size;
    vector<T> _data;
};

//
// Creates an image with 4 components
//
template <typename T>
inline image<vec<T, 4>> make_image4(int width, int height, int ncomp,
                                    const T* data, T alpha) {
    assert(ncomp >= 1 && ncomp <= 4);
    auto img = image<vec<T, 4>>{{width, height}};
    for (auto j = 0; j < height; j++) {
        for (auto i = 0; i < width; i++) {
            auto d = data + (j * width + i) * ncomp;
            switch (ncomp) {
                case 1: img[{i, j}] = {d[0], d[0], d[0], alpha}; break;
                case 2: img[{i, j}] = {d[0], d[1], 0, alpha}; break;
                case 3: img[{i, j}] = {d[0], d[1], d[2], alpha}; break;
                case 4: img[{i, j}] = {d[0], d[1], d[2], d[3]}; break;
            }
        }
    }
    return img;
}

//
// Conversion from srgb.
//
inline vec3f srgb_to_linear(const vec3f& srgb) {
    return {pow(srgb[0], 2.2f), pow(srgb[1], 2.2f), pow(srgb[2], 2.2f)};
}
inline vec4f srgb_to_linear(const vec4f& srgb) {
    return {pow(srgb[0], 2.2f), pow(srgb[1], 2.2f), pow(srgb[2], 2.2f),
            srgb[3]};
}
inline vec3f srgb_to_linear(const vec3b& srgb) {
    return srgb_to_linear(vec3f(srgb[0], srgb[1], srgb[2]) / 255.0f);
}
inline vec4f srgb_to_linear(const vec4b& srgb) {
    return srgb_to_linear(vec4f(srgb[0], srgb[1], srgb[2], srgb[3]) / 255.0f);
}

//
// Exposure/gamma correction
//
inline void exposure_gamma(const image_view<vec4f>& hdr, image_view<vec4f> ldr,
                           float exposure, float gamma) {
    assert(hdr.size() == ldr.size());
    auto s = pow(2.0f, exposure);
    for (auto j = 0; j < hdr.size()[1]; j++) {
        for (auto i = 0; i < hdr.size()[0]; i++) {
            auto v = hdr[{i, j}];
            v = {pow(s * v[0], 1 / gamma), pow(s * v[1], 1 / gamma),
                 pow(s * v[2], 1 / gamma), v[3]};
            ldr[{i, j}] = {
                (float)clamp(v[0], 0.0f, 1.0f), (float)clamp(v[1], 0.0f, 1.0f),
                (float)clamp(v[2], 0.0f, 1.0f), (float)clamp(v[3], 0.0f, 1.0f),
            };
        }
    }
}

inline void exposure_gamma(const image_view<vec4f>& hdr, image_view<vec4b> ldr,
                           float exposure, float gamma,
                           const vec2i& xy = {0, 0},
                           const vec2i& wh = {-1, -1}) {
    assert(hdr.size() == ldr.size());
    auto s = pow(2.0f, exposure);
    auto wh_ = vec2i{(wh[0] < 0) ? hdr.size()[0] - xy[0] : wh[0],
                     (wh[1] < 0) ? hdr.size()[1] - xy[1] : wh[1]};
    for (auto j = xy[1]; j < xy[1] + wh_[1]; j++) {
        for (auto i = xy[0]; i < xy[0] + wh_[0]; i++) {
            auto v = hdr[{i, j}];
            v = {pow(s * v[0], 1 / gamma), pow(s * v[1], 1 / gamma),
                 pow(s * v[2], 1 / gamma), v[3]};
            ldr[{i, j}] = {
                (unsigned char)(clamp(v[0], 0.0f, 1.0f) * 255),
                (unsigned char)(clamp(v[1], 0.0f, 1.0f) * 255),
                (unsigned char)(clamp(v[2], 0.0f, 1.0f) * 255),
                (unsigned char)(clamp(v[3], 0.0f, 1.0f) * 255),
            };
        }
    }
}

// -----------------------------------------------------------------------------
// TIMER
// -----------------------------------------------------------------------------

struct timer {
    timer(bool autostart = true) {
        if (autostart) start();
    }

    void start() {
        _start = std::chrono::steady_clock::now();
        _started = true;
    }

    void stop() {
        _end = std::chrono::steady_clock::now();
        _started = false;
    }

    double elapsed() {
        if (_started) stop();
        std::chrono::duration<double> diff = (_end - _start);
        return diff.count();
    }

   private:
    bool _started = false;
    std::chrono::time_point<std::chrono::steady_clock> _start, _end;
};

}  // namespace

#endif
