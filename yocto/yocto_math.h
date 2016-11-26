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
// complete, but unreadable of with tons of dependencies, or just as incomplete
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
// boost hash_combine, Pixar multijittered sampling
// and public domain code from
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

//
//
//

#ifndef _YMATH_H_
#define _YMATH_H_

#include <algorithm>
#include <array>
#include <cassert>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <type_traits>
#include <vector>

// -----------------------------------------------------------------------------
// CONSTANTS
// -----------------------------------------------------------------------------

namespace ym {

//
// pi
//
const float pif = 3.14159265f;
const double pi = 3.1415926535897932384626433832795;

//
// shortcut for numeric limits
//
const float max_float = std::numeric_limits<float>::max();

// -----------------------------------------------------------------------------
// BASIC MATH FUNCTIONS
// -----------------------------------------------------------------------------

//
// min/max
//
template <typename T>
inline T min(const T& x, const T& y) {
    return (x < y) ? x : y;
}
template <typename T>
inline T max(const T& x, const T& y) {
    return (x > y) ? x : y;
}

//
// swap utility
//
template <typename T>
inline void swap(T& a, T& b) {
    T c = a;
    a = b;
    b = c;
}

//
// Clamp a value between a min and max values
//
template <typename T>
inline T clamp(const T& x, const T& min_, const T& max_) {
    return min(max(x, min_), max_);
}
template <typename T, typename T1>
inline T clamp(const T& x, const T1& min_, const T1& max_) {
    return min(max(x, (T)min_), (T)max_);
}

//
// Check if a value if finite (will then be specialized on vectors)
//
template <typename T>
inline bool isfinite(const T& v) {
    return std::isfinite(v);
}

//
// Integer power of two
//
inline int pow2(int x) { return 1 << x; }

//
// Linear interpolation.
//
template <typename T, typename T1>
inline T lerp(const T& a, const T& b, T1 t) {
    return a * (1 - t) + b * t;
}

// -----------------------------------------------------------------------------
// LINEAR ALGEBRA TYPES (VECTORS and MATRICES)
// -----------------------------------------------------------------------------

//
// Vector of element of compile time length with default initializer,
// constant initialization and initialization from a C array. Data access
// via operator[].
//
template <typename T, int N>
struct vec {
    // default constructor
    vec() {
        for (auto i = 0; i < N; i++) (*this)[i] = T();
    }

    // member constructor
    template <typename... Args,
              typename = std::enable_if_t<sizeof...(Args) == N>>
    constexpr vec(Args... vv) : _v{T(vv)...} {}

    // list constructor
    constexpr vec(std::initializer_list<T> vv) {
        assert(N == vv.size());
        auto i = 0;
        for (auto&& e : vv) (*this)[i++] = e;
    }

    // one-member constructor
    constexpr explicit vec(T x) {
        for (auto i = 0; i < N; i++) (*this)[i] = x;
    }

    // constructor from shorter vector
    constexpr vec(const vec<T, N - 1>& va, const T& vb) {
        for (auto i = 0; i < N - 1; i++) (*this)[i] = va[i];
        (*this)[N - 1] = vb;
    }

    // conversion to shorter vector (ignore last value)
    explicit operator vec<T, N - 1>() {
        vec<T, N - 1> vv;
        for (auto i = 0; i < N - 1; i++) vv[i] = (*this)[i];
        return vv;
    }

    // conversion from different types
    template <typename T1>
    constexpr explicit vec(const vec<T1, N>& vv) {
        for (auto i = 0; i < N; i++) (*this)[i] = T(vv[i]);
    }

    // convert from STL
    constexpr vec(const std::array<T, N>& vv) {
        for (auto i = 0; i < N; i++) (*this)[i] = vv[i];
    }

    // convert to STL
    constexpr operator std::array<T, N>() const {
        std::array<T, N> vv;
        for (auto i = 0; i < N; i++) vv[i] = (*this)[i];
        return vv;
    }

    // element access
    constexpr const T& operator[](int i) const { return _v[i]; }
    constexpr T& operator[](int i) { return _v[i]; }

    // raw data access
    constexpr T* data() { return (T*)this; }
    constexpr const T* data() const { return (T*)this; }

   private:
    // elements
    T _v[N];
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

const auto one2f = vec2f(1, 1);
const auto one3f = vec3f(1, 1, 1);
const auto one4f = vec4f(1, 1, 1, 1);

const auto x2f = vec2f(1, 0);
const auto y2f = vec2f(0, 1);

const auto x3f = vec3f(1, 0, 0);
const auto y3f = vec3f(0, 1, 0);
const auto z3f = vec3f(0, 0, 1);

//
// Matrix of element of compile time dimensions stored in column major format.
// Allows for default initializer, initalization with columns and to/from
// smaller matrices.
// Access to columns via operator [] and raw data pointer data().
//
template <typename T, int N, int M>
struct mat {
    using V = vec<T, N>;

    constexpr mat() {
        for (int j = 0; j < M; j++) (*this)[j] = V();
    }

    // member constructor
    template <typename... Args,
              typename = std::enable_if_t<sizeof...(Args) == M>>
    constexpr mat(Args... args) : _v{V(args)...} {}

    // list constructor
    constexpr mat(std::initializer_list<V> vv) {
        assert(M == vv.size());
        auto i = 0;
        for (auto&& e : vv) (*this)[i++] = e;
    }

    // convert from STL
    constexpr mat(const std::array<V, M>& vv) {
        for (auto j = 0; j < M; j++) (*this)[j] = vv[j];
    }
    constexpr mat(const std::array<std::array<T, N>, M>& vv) {
        for (auto j = 0; j < M; j++) (*this)[j] = vv[j];
    }
    constexpr mat(const std::array<T, N * M>& vv) {
        for (auto j = 0; j < M; j++) {
            for (auto i = 0; i < N; i++) {
                (*this)[j][i] = vv[j * N + i];
            }
        }
    }

    // convert to STL
    constexpr operator std::array<V, M>() {
        std::array<V, M> vv;
        for (auto j = 0; j < M; j++) vv[j] = (*this)[j];
        return vv;
    }
    constexpr operator std::array<std::array<T, N>, M>() {
        std::array<std::array<T, N>, M> vv;
        for (auto j = 0; j < M; j++) vv[j] = (*this)[j];
        return vv;
    }
    constexpr operator std::array<T, N * M>() {
        std::array<T, N * M> vv;
        for (auto j = 0; j < M; j++)
            for (auto i = 0; i < N; i++) vv[j * N + i] = (*this)[j][i];
        return vv;
    }

    // element access
    constexpr const V& operator[](int i) const { return _v[i]; }
    constexpr V& operator[](int i) { return _v[i]; }

    // raw element access
    constexpr V* data() { return (V*)this; }
    constexpr const V* data() const { return (V*)this; }

   private:
    V _v[M];
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
// Components-wise comparison
//

template <typename T, int N>
inline bool isfinite(const vec<T, N>& a) {
    for (auto i = 0; i < N; i++)
        if (!isfinite(a[i])) return false;
    return true;
}

template <typename T, int N>
inline bool operator==(const vec<T, N>& a, const vec<T, N>& b) {
    for (auto i = 0; i < N; i++)
        if (a[i] != b[i]) return false;
    return true;
}

template <typename T, int N>
inline bool operator!=(const vec<T, N>& a, const vec<T, N>& b) {
    return !operator==(a, b);
}

template <typename T, int N>
inline bool operator<(const vec<T, N>& a, const vec<T, N>& b) {
    for (auto i = 0; i < N; i++)
        if (a[i] >= b[i]) return false;
    return true;
}

//
// Component-wise arithmentic.
//

template <typename T, int N>
inline vec<T, N> operator-(const vec<T, N>& a) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = -a[i];
    return c;
}

template <typename T, int N>
inline vec<T, N> operator+(const vec<T, N>& a, const vec<T, N>& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a[i] + b[i];
    return c;
}

template <typename T, int N>
inline vec<T, N> operator-(const vec<T, N>& a, const vec<T, N>& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a[i] - b[i];
    return c;
}

template <typename T, int N>
inline vec<T, N> operator*(const vec<T, N>& a, const vec<T, N>& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a[i] * b[i];
    return c;
}

template <typename T, int N, typename T1>
inline vec<T, N> operator*(const vec<T, N>& a, const T1& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a[i] * b;
    return c;
}

template <typename T1, typename T, int N>
inline vec<T, N> operator*(const T1& a, const vec<T, N>& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a * b[i];
    return c;
}

template <typename T, int N>
inline vec<T, N> operator/(const vec<T, N>& a, const vec<T, N>& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a[i] / b[i];
    return c;
}

template <typename T, int N, typename T1>
inline vec<T, N> operator/(const vec<T, N>& a, const T1 b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a[i] / b;
    return c;
}

template <typename T1, typename T, int N>
inline vec<T, N> operator/(const T1& a, const vec<T, N>& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a / b[i];
    return c;
}

//
// Component-wise assignment arithmentic.
//

template <typename T, int N>
inline vec<T, N>& operator+=(vec<T, N>& a, const vec<T, N>& b) {
    return a = a + b;
}

template <typename T, int N>
inline vec<T, N>& operator-=(vec<T, N>& a, const vec<T, N>& b) {
    return a = a - b;
}

template <typename T, int N>
inline vec<T, N>& operator*=(vec<T, N>& a, const vec<T, N>& b) {
    return a = a * b;
}

template <typename T, int N, typename T1>
inline vec<T, N>& operator*=(vec<T, N>& a, const T1& b) {
    return a = a * b;
}

template <typename T, int N>
inline vec<T, N>& operator/=(vec<T, N>& a, const vec<T, N>& b) {
    return a = a / b;
}

template <typename T, int N, typename T1>
inline vec<T, N>& operator/=(vec<T, N>& a, const T1 b) {
    return a = a / b;
}

//
// Component-wise min/max and clamping.
//

template <typename T, int N>
inline vec<T, N> min(const vec<T, N>& a, const vec<T, N>& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = min(a[i], b[i]);
    return c;
}

template <typename T, int N>
inline vec<T, N> max(const vec<T, N>& a, const vec<T, N>& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = max(a[i], b[i]);
    return c;
}

template <typename T, int N>
inline vec<T, N> min(const vec<T, N>& a, const T& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = min(a[i], b);
    return c;
}

template <typename T, int N>
inline vec<T, N> max(const vec<T, N>& a, const T& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = max(a[i], b);
    return c;
}

template <typename T, int N>
inline vec<T, N> clamp(const vec<T, N>& x, const T& min, const T& max) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = clamp(x[i], min, max);
    return c;
}

template <typename T, int N>
inline vec<T, N> clamp(const vec<T, N>& x, const vec<T, N>& min,
                       const vec<T, N>& max) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = clamp(x[i], min[i], max[i]);
    return c;
}

template <typename T, int N>
inline vec<T, N> clamplen(const vec<T, N> x, float max) {
    auto l = length(x);
    return (l > max) ? x * max / l : x;
}

//
// Element min/max.
//
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
template <typename T, int N>
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

//
// Component-wise functions.
//
template <typename T, int N>
inline vec<T, N> sqrt(const vec<T, N>& a) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = std::sqrt(a[i]);
    return c;
}

template <typename T, int N>
inline vec<T, N> pow(const vec<T, N>& a, const vec<T, N>& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = std::pow(a[i], b[i]);
    return c;
}

template <typename T, int N>
inline vec<T, N> pow(const vec<T, N>& a, const T& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = std::pow(a[i], b);
    return c;
}

template <typename T, int N>
inline vec<T, N> fmod(const vec<T, N>& a, const vec<T, N>& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = std::fmod(a[i], b[i]);
    return c;
}

template <typename T, int N>
inline vec<T, N> round(const vec<T, N>& a) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = std::round(a[i]);
    return c;
}

//
// Vector operations
//
template <typename T, int N>
inline T dot(const vec<T, N>& a, const vec<T, N>& b) {
    auto c = T(0);
    for (auto i = 0; i < N; i++) c += a[i] * b[i];
    return c;
}

template <typename T, int N>
inline T length(const vec<T, N>& a) {
    return std::sqrt(dot(a, a));
}

template <typename T, int N>
inline T lengthsqr(const vec<T, N>& a) {
    return dot(a, a);
}

template <typename T, int N>
inline vec<T, N> normalize(const vec<T, N>& a) {
    auto l = length(a);
    if (l == 0) return a;
    return a * (1 / l);
}

template <typename T, int N>
inline T dist(const vec<T, N>& a, const vec<T, N>& b) {
    return length(a - b);
}

template <typename T, int N>
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

template <typename T, int N>
T uangle(const vec<T, N>& a, const vec<T, N>& b) {
    auto d = dot(a, b);
    return d > 1 ? 0 : std::acos(d < -1 ? -1 : d);
}

template <typename T, int N>
T angle(const vec<T, N>& a, const vec<T, N>& b) {
    return uangle(normalize(a), normalize(b));
}

template <typename T, int N>
inline vec<T, N> nlerp(const vec<T, N>& a, const vec<T, N>& b, T t) {
    return normalize(lerp(a, b, t));
}

template <typename T, int N>
inline vec<T, N> slerp(const vec<T, N>& a, const vec<T, N>& b, T t) {
    auto th = uangle(a, b);
    return th == 0 ? a
                   : a * (std::sin(th * (1 - t)) / std::sin(th)) +
                         b * (std::sin(th * t) / std::sin(th));
}

// http://lolengine.net/blog/2013/09/21/picking-orthogonal-vector-combing-coconuts)
template <typename T>
inline vec<T, 3> orthogonal(const vec<T, 3>& v) {
    return std::abs(v[0]) > std::abs(v[2]) ? vec<T, 3>{-v[1], v[0], 0}
                                           : vec<T, 3>{0, -v[2], v[1]};
}

template <typename T>
inline vec<T, 3> orthonormalize(const vec<T, 3>& a, const vec<T, 3>& b) {
    return ym_normalize(a - b * ym_dot(a, b));
}

//
// Aggregate operations.
//
template <typename T, int M>
inline T sum(const vec<T, M>& a) {
    auto s = T(0);
    for (auto i = 0; i < M; i++) s += a[i];
    return s;
}

template <typename T, int M>
inline T mean(const vec<T, M>& a) {
    return sum(a) / M;
}

//
// Components-wise comparison
//

template <typename T, int N, int M>
inline bool operator==(const mat<T, N, M>& a, const mat<T, N, M>& b) {
    for (auto i = 0; i < M; i++)
        if (!(a[i] == b[i])) return false;
    return true;
}

//
// Components-wise arithmetic
//

template <typename T, int N, int M>
inline mat<T, M, N> operator-(const mat<T, N, M>& a) {
    mat<T, N, M> c;
    for (auto i = 0; i < M; i++) c[i] = -a[i];
    return c;
}

template <typename T, int N, int M>
inline mat<T, M, N> operator+(const mat<T, N, M>& a, const mat<T, N, M>& b) {
    mat<T, N, M> c;
    for (auto i = 0; i < M; i++) c[i] = a[i] + b[i];
    return c;
}

//
// Vector/Matrix multiplies
//

template <typename T, int N, int M>
inline mat<T, M, N> operator*(const mat<T, N, M>& a, T b) {
    mat<T, N, M> c;
    for (auto i = 0; i < M; i++) c[i] = a[i] * b;
    return c;
}

template <typename T, int N, int M>
inline mat<T, M, N> operator/(const mat<T, N, M>& a, T b) {
    mat<T, N, M> c;
    for (auto i = 0; i < M; i++) c[i] = a[i] / b;
    return c;
}

template <typename T, int N, int M>
inline vec<T, N> operator*(const mat<T, N, M>& a, const vec<T, M>& b) {
    vec<T, N> c;
    for (auto j = 0; j < M; j++) c += a[j] * b[j];
    return c;
}

template <typename T, int N, int M, int K>
inline mat<T, N, M> operator*(const mat<T, N, K>& a, const mat<T, K, M>& b) {
    mat<T, N, M> c;
    for (auto j = 0; j < M; j++) c[j] = a * b[j];
    return c;
}

//
// Linear Algebra operations.
//

template <typename T, int N>
vec<T, N> diagonal(const mat<T, N, N>& a) {
    vec<T, N> d;
    for (auto i = 0; i < N; i++) d[i] = a[i][i];
    return d;
}

template <typename T, int N, int M>
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

template <typename T, int N>
inline mat<T, N, N> inverse(const mat<T, N, N>& a) {
    return adjugate(a) / determinant(a);
}

//
// Support for iteration
//

template <typename T, int N>
T* begin(vec<T, N>& a) {
    return a;
}

template <typename T, int N>
const T* begin(const vec<T, N>& a) {
    return a;
}

template <typename T, int N>
T* end(vec<T, N>& a) {
    return a + N;
}

template <typename T, int N>
const T* end(const vec<T, N>& a) {
    return a + N;
}

template <typename T, int N, int M>
vec<T, M>* begin(mat<T, M, N>& a) {
    return a;
}

template <typename T, int N, int M>
const vec<T, N>* begin(const mat<T, M, N>& a) {
    return a;
}

template <typename T, int N, int M>
vec<T, N>* end(mat<T, M, N>& a) {
    return a + M;
}

template <typename T, int N, int M>
const vec<T, N>* end(const mat<T, M, N>& a) {
    return a + M;
}

};  // namespace

// -----------------------------------------------------------------------------
// SPECIALIZATIONS FOR STD
// -----------------------------------------------------------------------------

namespace std {

//
// Min/max type extents
//

template <typename T, int N>
struct std::numeric_limits<ym::vec<T, N>> {
    using V = ym::vec<T, N>;
    static constexpr bool is_specialized = true;
    static constexpr bool is_signed = std::numeric_limits<T>::is_signed;
    static constexpr bool is_integer = std::numeric_limits<T>::is_integer;
    static constexpr bool is_exact = std::numeric_limits<T>::is_exact;
    static constexpr bool has_infinity = std::numeric_limits<T>::has_infinity;
    constexpr static V lowest() { return V(std::numeric_limits<T>::lowest()); }
    constexpr static V max() { return V(std::numeric_limits<T>::max()); }
    constexpr static V min() { return V(std::numeric_limits<T>::min()); }
    constexpr static V epsilon() {
        return V(std::numeric_limits<T>::epsilon());
    }
    constexpr static V round_arror() {
        return V(std::numeric_limits<T>::round_arror());
    }
    constexpr static V infinity() {
        return V(std::numeric_limits<T>::infinity());
    }
    constexpr static V quite_NaN() {
        return V(std::numeric_limits<T>::quite_NaN());
    }
    constexpr static V signaling_NaN() {
        return V(std::numeric_limits<T>::signaling_NaN());
    }
    constexpr static V denorm_min() {
        return V(std::numeric_limits<T>::denorm_min());
    }
};

//
// A hash function for vecs
//
template <typename T, int N>
struct hash<ym::vec<T, N>> {
    using V = ym::vec<T, N>;
    std::hash<T> Th;
    constexpr size_t operator()(const V& v) const {
        size_t h = 0;
        for (auto i = 0; i < N; i++) {
            // embads hash_combine below
            h ^= (Th(v[i]) + 0x9e3779b9 + (h << 6) + (h >> 2));
        }
        return h;
    }
};

}  // namepace

namespace ym {

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
template <typename T, int N>
struct frame {
    using V = vec<T, N>;
    using M = mat<T, N, N>;

    constexpr frame() {
        for (auto i = 0; i < N + 1; i++) (*this)[i] = V();
    }

    constexpr frame(const M& m, const V& t) {
        for (auto i = 0; i < N; i++) (*this)[i] = m[i];
        (*this)[N] = t;
    }

    // member constructor
    template <typename... Args,
              typename = std::enable_if_t<sizeof...(Args) == N + 1>>
    constexpr frame(Args... args) : _v{V(args)...} {}

    // list constructor
    constexpr frame(std::initializer_list<V> vv) {
        assert(N + 1 == vv.size());
        auto i = 0;
        for (auto&& e : vv) (*this)[i++] = e;
    }

    constexpr explicit frame(const mat<T, N + 1, N + 1>& mat) {
        for (auto j = 0; j < N + 1; j++)
            for (auto i = 0; i < N; i++) (*this)[j][i] = mat[j][i];
    }

    constexpr explicit operator mat<T, N + 1, N + 1>() const {
        mat<T, N + 1, N + 1> m;
        for (auto j = 0; j < N; j++) {
            m[j] = {(*this)[j], 0};
        }
        m[N] = {(*this)[N], 1};
        return m;
    }

    // element access
    constexpr const V& operator[](int i) const { return _v[i]; }
    constexpr V& operator[](int i) { return _v[i]; }

    // view as linear + translation
    constexpr const V& o() const {
        static_assert(N > 1, "length");
        return (*this)[N];
    }
    constexpr V& o() {
        static_assert(N > 1, "length");
        return (*this)[N];
    }
    constexpr const M& m() const {
        static_assert(N > 1, "length");
        return *(M*)this;
    }
    constexpr M& m() {
        static_assert(N > 1, "length");
        return *(M*)this;
    }

    // raw element access
    constexpr V* data() { return (V*)this; }
    constexpr const V* data() const { return (V*)this; }

   private:
    V _v[N + 1];
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
// Component-wise comparisons
//

template <typename T, int N>
inline bool operator==(const frame<T, N>& a, const frame<T, N>& b) {
    return a.m() == b.m() && a.o() == b.o();
}

//
// Matrix operations for transforms.
//

template <typename T, int N>
inline frame<T, N> operator*(const frame<T, N>& a, const frame<T, N>& b) {
    return {a.m() * b.m(), a.m() * b.o() + a.o()};
}

template <typename T, int N>
inline frame<T, N> inverse(const frame<T, N>& a) {
    auto minv = transpose(a.m());
    return {minv, -(minv * a.o())};
}

template <typename T, int N>
vec<T, N>* begin(frame<T, N>& a) {
    return a;
}

template <typename T, int N>
const vec<T, N>* begin(const frame<T, N>& a) {
    return a;
}

template <typename T, int N>
vec<T, N>* end(frame<T, N>& a) {
    return a + N + 1;
}

template <typename T, int N>
const vec<T, N>* end(const frame<T, N>& a) {
    return a + N + 1;
}

// -----------------------------------------------------------------------------
// QUATERNIONS
// -----------------------------------------------------------------------------

//
// Quaternions implemented as a vec4. Data access via operator[].
// Quaterions are xi + yj + zk + w
//
template <typename T, int N>
struct quat : vec<T, 4> {
    // elements
    T v[N];

    // default constructor
    constexpr quat() {
        for (auto i = 0; i < N; i++) v[i] = T();
    }

    // member constructor
    template <typename... Args,
              typename = std::enable_if_t<sizeof...(Args) == N>>
    constexpr quat(Args... vv) : v{T(vv)...} {}

    // list constructor
    constexpr quat(std::initializer_list<T> vv) {
        assert(N == vv.size());
        auto i = 0;
        for (auto&& e : vv) v[i++] = e;
    }

    // conversions
    operator mat<T, 3, 3>() {
        return {
            {v[3] * v[3] + v[0] * v[0] - v[1] * v[1] - v[2] * v[2],
             (v[0] * v[1] + v[2] * v[3]) * 2, (v[2] * v[0] - v[1] * v[3]) * 2},
            {(v[0] * v[1] - v[2] * v[3]) * 2,
             v[3] * v[3] - v[0] * v[0] + v[1] * v[1] - v[2] * v[2],
             (v[1] * v[2] + v[0] * v[3]) * 2},
            {(v[2] * v[0] + v[1] * v[3]) * 2, (v[1] * v[2] - v[0] * v[3]) * 2,
             v[3] * v[3] - v[0] * v[0] - v[1] * v[1] + v[2] * v[2]}};
    }

    // axis angle access
    constexpr T angle() const { return std::acos(v[2]) * 2; }
    constexpr vec<T, 3> axis() {
        return normalize(vec<T, 3>(v[0], v[1], v[2]));
    }

    // raw data access
    constexpr T* data() { return v; }
    constexpr const T* data() const { return v; }
};

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
template <typename T>
struct range {
    T min, max;

    constexpr range()
        : min(std::numeric_limits<T>::max()),
          max(std::numeric_limits<T>::lowest()) {}

    constexpr range(const T& m, const T& M) : min(m), max(M) {}

    constexpr T size() const { return max - min; }
    constexpr T center() const { return (max + min) / 2; }

    constexpr T& operator[](int i) { return (&min)[i]; }
    constexpr const T& operator[](int i) const { return (&min)[i]; }
};

//
// Axis aligned bounding box usings.
//

using range1f = range<float>;
using range2f = range<vec2f>;
using range3f = range<vec3f>;

using range1i = range<int>;
using range2i = range<vec2i>;
using range3i = range<vec3i>;

//
// Axis aligned bounding box constants.
//

const auto invalid_range2f = range2f();
const auto invalid_range3f = range3f();

//
// Axis aligned bounding box operations.
//

template <typename T>
inline range<T> make_range(const std::initializer_list<T>& v) {
    auto r = range<T>();
    for (auto&& vv : v) r += vv;
    return r;
}

template <typename T>
inline range<T> rexpand(const range<T>& a, const T& b) {
    return {min(a.min, b), max(a.max, b)};
}

template <typename T>
inline range<T> rexpand(const range<T>& a, const range<T>& b) {
    return {min(a.min, b.min), max(a.max, b.max)};
}

template <typename T>
inline range<T> operator+(const range<T>& a, const T& b) {
    return rexpand(a, b);
}

template <typename T>
inline range<T> operator+(const range<T>& a, const range<T>& b) {
    return rexpand(a, b);
}

template <typename T>
inline range<T>& operator+=(range<T>& a, const T& b) {
    return a = rexpand(a, b);
}

template <typename T>
inline range<T>& operator+=(range<T>& a, const range<T>& b) {
    return a = rexpand(a, b);
}

template <typename T>
inline T rcenter(const range<T>& a) {
    return (a.min + a.max) / 2;
}

template <typename T>
inline T rsize(const range<T>& a) {
    return a.max - a.min;
}

// -----------------------------------------------------------------------------
// RAYS
// -----------------------------------------------------------------------------

//
// Rays with origin, direction and min/max t value.
//
template <typename T, int N>
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

template <typename T, int N>
inline vec<T, N> transform_point(const mat<T, N + 1, N + 1>& a,
                                 const vec<T, N>& b) {
    // make it generic
    auto vb = vec<T, N + 1>{b, 1};
    auto tvb = a * vb;
    return vec<T, N>(tvb) / tvb[3];
}

template <typename T, int N>
inline vec<T, N> transform_vector(const mat<T, N + 1, N + 1>& a,
                                  const vec<T, N>& b) {
    // make it generic
    auto vb = vec<T, N + 1>{b, 0};
    auto tvb = a * vb;
    return vec<T, N>(tvb);
}

template <typename T, int N>
inline vec<T, N> transform_direction(const mat<T, N + 1, N + 1>& a,
                                     const vec<T, N>& b) {
    return normalize(transform_vector(a, b));
}

template <typename T, int N>
inline vec<T, N> transform_point(const frame<T, N>& a, const vec<T, N>& b) {
    return a.m() * b + a.o();
}

template <typename T, int N>
inline vec<T, N> transform_vector(const frame<T, N>& a, const vec<T, N>& b) {
    return a.m() * b;
}

template <typename T, int N>
inline vec<T, N> transform_direction(const frame<T, N>& a, const vec<T, N>& b) {
    return a.m() * b;
}

template <typename T, int N>
inline frame<T, N> transform_frame(const frame<T, N>& a, const frame<T, N>& b) {
    return {a.m() * b.m(), a.m() * b.o() + a.o()};
}

template <typename T, int N>
inline ray<T, N> transform_ray(const frame<T, N>& a, const ray<T, N>& b) {
    return {transform_point(a, b.o), transform_direction(a, b.d), b.tmin,
            b.tmax};
}

template <typename T, int N>
inline ray<T, N> transform_ray(const mat<T, N + 1, N + 1>& a,
                               const ray<T, N>& b) {
    return {transform_point(a, b.o), transform_direction(a, b.d), b.tmin,
            b.tmax};
}

template <typename T>
inline range<vec<T, 3>> transform_bbox(const frame<T, 3>& xform,
                                       const range<vec<T, 3>>& bbox) {
    vec<T, 3> corners[8] = {
        {bbox.min[0], bbox.min[1], bbox.min[2]},
        {bbox.min[0], bbox.min[1], bbox.max[2]},
        {bbox.min[0], bbox.max[1], bbox.min[2]},
        {bbox.min[0], bbox.max[1], bbox.max[2]},
        {bbox.max[0], bbox.min[1], bbox.min[2]},
        {bbox.max[0], bbox.min[1], bbox.max[2]},
        {bbox.max[0], bbox.max[1], bbox.min[2]},
        {bbox.max[0], bbox.max[1], bbox.max[2]},
    };
    auto xformed = range<vec<T, 3>>();
    for (auto j = 0; j < 8; j++) xformed += transform_point(xform, corners[j]);
    return xformed;
}

template <typename T>
inline range<vec<T, 3>> transform_bbox(const mat<T, 4, 4>& xform,
                                       const range<vec<T, 3>>& bbox) {
    vec<T, 3> corners[8] = {
        {bbox.min[0], bbox.min[1], bbox.min[2]},
        {bbox.min[0], bbox.min[1], bbox.max[2]},
        {bbox.min[0], bbox.max[1], bbox.min[2]},
        {bbox.min[0], bbox.max[1], bbox.max[2]},
        {bbox.max[0], bbox.min[1], bbox.min[2]},
        {bbox.max[0], bbox.min[1], bbox.max[2]},
        {bbox.max[0], bbox.max[1], bbox.min[2]},
        {bbox.max[0], bbox.max[1], bbox.max[2]},
    };
    auto xformed = range<vec<T, 3>>();
    for (auto j = 0; j < 8; j++) xformed += transform_point(xform, corners[j]);
    return xformed;
}

//
// Rotation matrix crom axis/angle
//
template <typename T>
inline mat<T, 3, 3> rotation_mat3(const vec<T, 3>& axis, T angle) {
    auto s = std::sin(angle), c = std::cos(angle);
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
// Scaling transform (in this case the frame is broken and used only as affine)
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
    return (mat<T, 4, 4>)lookat_frame3(eye, center, up);
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
    auto y = near * std::tan(fovy / 2);
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
// Turntable for UI navigation from a from/to/up parametrization of the camera.
//
template <typename T>
inline void turntable(vec<T, 3>& from, vec<T, 3>& to, vec<T, 3>& up,
                      const vec<T, 2>& rotate, T dolly, const vec<T, 2>& pan) {
    // rotate if necessary
    if (rotate[0] || rotate[1]) {
        auto z = ym_normalize(*to - *from);
        auto lz = ym_dist(*to, *from);
        auto phi = std::atan2(z[2], z[0]) + rotate[0];
        auto theta = std::acos(z[1]) + rotate[1];
        theta = max(T(0.001), min(theta, T(pi - 0.001)));
        auto nz = vec<T, 3>{std::sin(theta) * std::cos(phi) * lz,
                            std::cos(theta) * lz,
                            std::sin(theta) * std::sin(phi) * lz};
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
// Turntable for UI navigation for a frame/distance parametrization of the
// camera.
//
template <typename T>
inline void turntable(frame<T, 3>& frame, float& focus, const vec<T, 2>& rotate,
                      T dolly, const vec<T, 2>& pan) {
    // rotate if necessary
    if (rotate[0] || rotate[1]) {
        auto phi = std::atan2(frame[2][2], frame[2][0]) + rotate[0];
        auto theta = std::acos(frame[2][1]) + rotate[1];
        theta = max(T(0.001), min(theta, T(pi - 0.001)));
        auto new_z = vec<T, 3>(std::sin(theta) * std::cos(phi), std::cos(theta),
                               std::sin(theta) * std::sin(phi));
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
// PCG random numbers. A family of random number generators that supports
// multiple sequences. In our code, we allocate one sequence for each sample.
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
// Init a random number generator with a state state from the sequence seq.
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
// A function to combine 64 bit hashes with semantics as boost::hash_combine
//
inline size_t hash_combine(size_t a, size_t b) {
    return a ^ (b + 0x9e3779b9 + (a << 6) + (a >> 2));
}

}  // namespace

// -----------------------------------------------------------------------------
// VIEW CONTAINERS
// -----------------------------------------------------------------------------

namespace ym {

//
// An array_view is a non-owining reference to an array with an API similar
// to a std::vector/std::array containers, but without reallocation.
// This is inspired, but significantly simpler than
// gsl::span https://github.com/Microsoft/GSL or std::array_view.
//
template <typename T>
struct array_view {
    // constructors
    constexpr array_view() noexcept : _num(0), _data(nullptr) {}
    constexpr array_view(int num, T* data) noexcept : _num(num), _data(data) {}
    constexpr array_view(T* begin, T* end) noexcept
        : _num((int)(end - begin)), _data(begin) {}
    constexpr array_view(const array_view<std::add_const_t<T>>& av)
        : _num(av._num), _data(av._data) {}
    constexpr array_view(std::vector<T>& av)
        : _num((int)av.size()), _data(av.data()) {}
    constexpr array_view(const std::vector<T>& av)
        : _num((int)av.size()), _data((T*)av.data()) {}

    // size
    constexpr int size() const noexcept { return _num; }
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
    int _num;
    T* _data;
};

}  // namesapce

// -----------------------------------------------------------------------------
// CONTAINER OPERATIONS
// -----------------------------------------------------------------------------

#include <array>
#include <string>
#include <vector>

//
// Vector append operations using operators (Python-like).
//

template <typename T>
inline std::vector<T>& operator+=(std::vector<T>& a, const T& b) {
    a.push_back(b);
    return a;
}

template <typename T>
inline std::vector<T>& operator+=(std::vector<T>& a, const std::vector<T>& b) {
    for (int i = 0; i < b.size(); i++) a.push_back(b[i]);
    return a;
}

template <typename T>
inline std::vector<T> operator+(const std::vector<T>& a,
                                const std::vector<T>& b) {
    std::vector<T> c = a;
    c += b;
    return c;
}

template <typename T>
inline std::vector<T> operator+(const std::vector<T>& a, const T& b) {
    std::vector<T> c = a;
    c += b;
    return c;
}

// -----------------------------------------------------------------------------
// IMAGES
// -----------------------------------------------------------------------------

namespace ym {

//
// An image_view is a non-owining reference to an array that allows to access it
// as a row-major image, but without reallocation.
// This is inspired, but significantly simpler than
// gsl::multi_span https://github.com/Microsoft/GSL or std::array_view.
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
// Generic image with pixels stored as a vector. For operations use image_views.
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
        return image_view<T>(_size, _data.data());
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
    std::vector<T> _data;
};

//
// Conversion from srgb.
//
inline vec3f srgb_to_linear(const vec3f& srgb) { return pow(srgb, 1 / 2.2f); }
inline vec4f srgb_to_linear(const vec4f& srgb) {
    auto rgb = pow(vec3f{srgb[0], srgb[1], srgb[2]}, 2.2f);
    return vec4f(rgb[0], rgb[1], rgb[2], srgb[3]);
}
inline vec3f srgb_to_linear(const vec3b& srgb) {
    return srgb_to_linear(vec3f(srgb) / 255.0f);
}
inline vec4f srgb_to_linear(const vec4b& srgb) {
    return srgb_to_linear(vec4f(srgb) / 255.0f);
}

}  // namespace

// -----------------------------------------------------------------------------
// TIMER
// -----------------------------------------------------------------------------

namespace ym {

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
