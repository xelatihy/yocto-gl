///
/// YOCTO_MATH: a collection of vector math functions and simple containers
/// used to implement YOCTO. Features include
/// - static length float vectors, with specialization for 2, 3, 4 length
/// - static length matrices, with specialization for 2x2, 3x3, 4x4
/// - affine and rigid transforms
/// - linear algebra operations and transforms for fixed length matrices/vecs
/// - axis aligned bounding boxes
/// - rays
/// - random number generation via PCG32
/// - a few hash functions
/// - timer (depends on C++11 chrono)
///
/// While we tested this library in the implementation of our other ones, we
/// consider this code incomplete and remommend to use a more complete math
/// library. So use it at your own peril.
///
/// We developed our own library since we felt that all existing ones are either
/// complete, but unreadable or with lots of dependencies, or just as incomplete
/// and untested as ours.
///
///
/// COMPILATION:
///
/// This library can only be used as a header only library in C++ since it uses
/// templates for its basic types. Some templated types and functions use
/// specialization for easier access and simpler compilation. Specialization
/// can be disabled by defining YM_NO_SPECIALIZATION.
///
///
/// HISTORY:
/// - v 0.10: specialize some type and functions
/// - v 0.9: bbox containment tests
/// - v 0.8: remove std:array as base class for better control
/// - v 0.7: doxygen comments
/// - v 0.6: uniformed internal names
/// - v 0.5: simplification of constructors, raname bbox -> bbox
/// - v 0.4: overall type simplification
/// - v 0.3: internal C++ refactoring
/// - v 0.2: use of STL containers; removal of yocto containers
/// - v 0.1: C++ only implementation
/// - v 0.0: initial release in C99
///
///
///
/// ACKNOLEDGEMENTS
///
/// This library includes code from the PCG random number generator,
/// boost hash_combine, Pixar multijittered sampling, code from "Real-Time
/// Collision Detection" by Christer Ericson and public domain code from
/// - https://github.com/sgorsten/linalg
/// - https://gist.github.com/badboy/6267743
///
/// For design ideas of a modern math library, it was very helpful to look at
/// - https://github.com/sgorsten/linalg
///
namespace ym {}

//
// LICENSE:
//
// Copyright (c) 2016 -- 2017 Fabio Pellacini
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
#include <vector>

// HACK to avoid compilation with MSVC2015 without dirtying code
#ifdef _WIN32
#define constexpr
#endif

namespace ym {

// -----------------------------------------------------------------------------
// BASIC TYPEDEFS
// -----------------------------------------------------------------------------

/// @name basic typedefs
/// @{

/// convenient typedef for bytes
using byte = unsigned char;

/// @}

// -----------------------------------------------------------------------------
// MATH CONSTANTS
// -----------------------------------------------------------------------------

/// @name main constants
/// @{

/// pi (float)
const float pif = 3.14159265f;
/// pi (double)
const double pi = 3.1415926535897932384626433832795;

/// shortcat for float max value
constexpr const auto flt_max = std::numeric_limits<float>::max();
/// shortcat for float min value
constexpr const auto flt_min = std::numeric_limits<float>::lowest();
/// shortcat for int max value
constexpr const auto int_max = std::numeric_limits<int>::max();
/// shortcat for int min value
constexpr const auto int_min = std::numeric_limits<int>::min();

/// @}

// -----------------------------------------------------------------------------
// BASIC MATH FUNCTIONS
// -----------------------------------------------------------------------------

/// @name basic math functions
/// @{

/// Safe minimum value.
template <typename T>
constexpr inline T min(T x, T y) {
    return (x < y) ? x : y;
}

/// Safe maximum value.
template <typename T>
constexpr inline T max(T x, T y) {
    return (x > y) ? x : y;
}

/// Clamp a value between a minimum and a maximum.
template <typename T>
constexpr inline T clamp(T x, T min_, T max_) {
    return min(max(x, min_), max_);
}

/// Linear interpolation.
template <typename T>
constexpr inline T lerp(T a, T b, T t) {
    return a * (1 - t) + b * t;
}

/// Integer power of two
constexpr inline int pow2(int x) { return 1 << x; }

/// @}

// -----------------------------------------------------------------------------
// VECTORS
// -----------------------------------------------------------------------------

/// @name vectors
/// @{

///
/// Vector of elements of compile time dimension with default initializer,
/// and conversione to/from std::array. Element access via operator[]. Supports
/// all std::array operations.
///
template <typename T, int N>
struct vec {
    /// default constructor
    constexpr vec() {
        for (auto i = 0; i < N; i++) v[i] = 0;
    }

    /// list constructor
    constexpr vec(const std::initializer_list<T>& vv) {
        assert(N == vv.size());
        auto i = 0;
        for (auto&& e : vv) v[i++] = e;
    }

    /// copy constructor from std::array
    constexpr vec(const std::array<T, N>& vv) {
        for (auto i = 0; i < N; i++) v[i] = vv[i];
    }

    /// conversion to std::array
    constexpr operator std::array<T, N>() const {
        auto c = std::array<T, N>();
        for (auto i = 0; i < N; i++) c[i] = v[i];
        return c;
    }

    /// element access
    constexpr T& operator[](int i) { return v[i]; }
    /// element access
    constexpr const T& operator[](int i) const { return v[i]; }

    /// element data
    T v[N];
};

#ifndef YM_NO_SPECIALIZATION

///
/// Specialization of vectors for 1 component and float coordinates.
///
template <>
struct vec<float, 1> {
    /// size
    constexpr static const int N = 1;
    /// type
    using T = float;

    /// default constructor
    constexpr vec() : x{0} {}
    /// element constructor
    constexpr vec(T x) : x{x} {}

    /// copy constructor from std::array
    constexpr vec(const std::array<T, N>& vv) : x{vv[0]} {}
    /// conversion to std::array
    constexpr operator std::array<T, N>() const { return {x}; }

    /// element access
    constexpr T& operator[](int i) { return v[i]; }
    /// element access
    constexpr const T& operator[](int i) const { return v[i]; }

    /// element data
    union {
        struct {
            T v[N];
        };
        struct {
            T x;
        };
    };
};

///
/// Specialization of vectors for 2 components and float coordinates.
///
template <>
struct vec<float, 2> {
    /// size
    constexpr static const int N = 2;
    /// type
    using T = float;

    /// default constructor
    constexpr vec() : x{0}, y{0} {}
    /// element constructor
    constexpr vec(T x, T y) : x{x}, y{y} {}

    /// copy constructor from std::array
    constexpr vec(const std::array<T, N>& vv) : x{vv[0]}, y{vv[1]} {}
    /// conversion to std::array
    constexpr operator std::array<T, N>() const { return {x, y}; }

    /// element access
    constexpr T& operator[](int i) { return v[i]; }
    /// element access
    constexpr const T& operator[](int i) const { return v[i]; }

    /// element data
    union {
        struct {
            T v[N];
        };
        struct {
            T x, y;
        };
    };
};

///
/// Specialization of vectors for 3 components and float coordinates.
///
template <>
struct vec<float, 3> {
    /// size
    constexpr static const int N = 3;
    /// type
    using T = float;

    /// default constructor
    constexpr vec() : x{0}, y{0}, z{0} {}
    /// element constructor
    constexpr vec(T x, T y, T z) : x{x}, y{y}, z{z} {}

    /// copy constructor from std::array
    constexpr vec(const std::array<T, N>& vv) : x{vv[0]}, y{vv[1]}, z{vv[2]} {}
    /// conversion to std::array
    constexpr operator std::array<T, N>() const { return {x, y, z}; }

    /// element access
    constexpr T& operator[](int i) { return v[i]; }
    /// element access
    constexpr const T& operator[](int i) const { return v[i]; }

    /// element data
    union {
        struct {
            T v[N];
        };
        struct {
            T x, y, z;
        };
    };
};

///
/// Specialization of vectors for 4 components and float coordinates.
///
template <>
struct vec<float, 4> {
    /// size
    constexpr static const int N = 4;
    /// type
    using T = float;

    /// default constructor
    constexpr vec() : x{0}, y{0}, z{0}, w{0} {}
    /// element constructor
    constexpr vec(T x, T y, T z, T w) : x{x}, y{y}, z{z}, w{w} {}

    /// copy constructor from std::array
    constexpr vec(const std::array<T, N>& vv)
        : x{vv[0]}, y{vv[1]}, z{vv[2]}, w{vv[3]} {}
    /// conversion to std::array
    constexpr operator std::array<T, N>() const { return {x, y, z, w}; }

    /// element access
    constexpr T& operator[](int i) { return v[i]; }
    /// element access
    constexpr const T& operator[](int i) const { return v[i]; }

    /// element data
    union {
        struct {
            T v[N];
        };
        struct {
            T x, y, z, w;
        };
    };
};

#endif

/// @}

/// @name vector typedefs

/// 1-dimensional float vector
using vec1f = vec<float, 1>;
/// 2-dimensional float vector
using vec2f = vec<float, 2>;
/// 3-dimensional float vector
using vec3f = vec<float, 3>;
/// 4-dimensional float vector
using vec4f = vec<float, 4>;

/// 1-dimensional int vector
using vec1i = vec<int, 1>;
/// 2-dimensional int vector
using vec2i = vec<int, 2>;
/// 3-dimensional int vector
using vec3i = vec<int, 3>;
/// 4-dimensional int vector
using vec4i = vec<int, 4>;

/// 1-dimensional byte vector
using vec1b = vec<byte, 1>;
/// 2-dimensional byte vector
using vec2b = vec<byte, 2>;
/// 3-dimensional byte vector
using vec3b = vec<byte, 3>;
/// 4-dimensional byte vector
using vec4b = vec<byte, 4>;

/// @}

/// @name vector initializations
/// @{

/// Initialize a zero vector.
template <typename T, int N>
constexpr inline vec<T, N> zero_vec() {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = 0;
    return c;
}

/// Sepcialization of Initialize a zero vector.
template <>
constexpr inline vec3f zero_vec() {
    return {0, 0, 0};
}

/// @}

/// @name vector constants
/// @{

/// 1-dimensional float zero vector
const auto zero1f = zero_vec<float, 1>();
/// 2-dimensional float zero vector
const auto zero2f = zero_vec<float, 2>();
/// 3-dimensional float zero vector
const auto zero3f = zero_vec<float, 3>();
/// 4-dimensional float zero vector
const auto zero4f = zero_vec<float, 4>();

/// 1-dimensional int zero vector
const auto zero1i = zero_vec<int, 1>();
/// 2-dimensional int zero vector
const auto zero2i = zero_vec<int, 2>();
/// 3-dimensional int zero vector
const auto zero3i = zero_vec<int, 3>();
/// 4-dimensional int zero vector
const auto zero4i = zero_vec<int, 4>();

/// @}

/// @name vector iteration functions
/// @{

/// iteration support
template <typename T, int N>
constexpr inline T* begin(vec<T, N>& a) {
    return a.v;
}

/// iteration support
template <typename T, int N>
constexpr inline const T* begin(const vec<T, N>& a) {
    return a.v;
}

/// iteration support
template <typename T, int N>
constexpr inline T* end(vec<T, N>& a) {
    return a.v + N;
}

/// iteration support
template <typename T, int N>
constexpr inline const T* end(const vec<T, N>& a) {
    return a.v + N;
}

/// @}

/// @name vector logical operations
/// @{

/// vector operator ==
template <typename T, int N>
constexpr inline bool operator==(const vec<T, N>& a, const vec<T, N>& b) {
    for (auto i = 0; i < N; i++)
        if (a[i] != b[i]) return false;
    return true;
}

/// vector operator !=
template <typename T, int N>
constexpr inline bool operator!=(const vec<T, N>& a, const vec<T, N>& b) {
    return !(a == b);
}

/// @}

#ifndef YM_NO_SPECIALIZATION

/// @name specialization of vector logical operations
/// @{

/// vector operator ==
template <>
constexpr inline bool operator==(const vec2f& a, const vec2f& b) {
    return a.x == b.x && a.y == b.y;
}

/// vector operator !=
template <>
constexpr inline bool operator!=(const vec2f& a, const vec2f& b) {
    return a.x != b.x || a.y != b.y;
}

/// vector operator ==
template <>
constexpr inline bool operator==(const vec3f& a, const vec3f& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}

/// vector operator !=
template <>
constexpr inline bool operator!=(const vec3f& a, const vec3f& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z;
}

/// vector operator ==
template <>
constexpr inline bool operator==(const vec4f& a, const vec4f& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}

/// vector operator !=
template <>
constexpr inline bool operator!=(const vec4f& a, const vec4f& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w;
}

/// @}

#endif

/// @name vector arithmetic operations
/// @{

/// vector operator +
template <typename T, int N>
constexpr inline vec<T, N> operator+(const vec<T, N>& a) {
    return a;
}

/// vector operator -
template <typename T, int N>
constexpr inline vec<T, N> operator-(const vec<T, N>& a) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = -a[i];
    return c;
}

/// vector operator +
template <typename T, int N>
constexpr inline vec<T, N> operator+(const vec<T, N>& a, const vec<T, N>& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a[i] + b[i];
    return c;
}

/// vector operator -
template <typename T, int N>
constexpr inline vec<T, N> operator-(const vec<T, N>& a, const vec<T, N>& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a[i] - b[i];
    return c;
}

/// vector operator *
template <typename T, int N>
constexpr inline vec<T, N> operator*(const vec<T, N>& a, const vec<T, N>& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a[i] * b[i];
    return c;
}

/// vector operator *
template <typename T, int N>
constexpr inline vec<T, N> operator*(const vec<T, N>& a, const T b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a[i] * b;
    return c;
}

/// vector operator *
template <typename T, int N>
constexpr inline vec<T, N> operator*(const T a, const vec<T, N>& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a * b[i];
    return c;
}

/// vector operator /
template <typename T, int N>
constexpr inline vec<T, N> operator/(const vec<T, N>& a, const vec<T, N>& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a[i] / b[i];
    return c;
}

/// vector operator /
template <typename T, int N>
constexpr inline vec<T, N> operator/(const vec<T, N>& a, const T b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a[i] / b;
    return c;
}

/// vector operator /
template <typename T, int N>
constexpr inline vec<T, N> operator/(const T a, const vec<T, N>& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a / b[i];
    return c;
}

/// @}

#ifndef YM_NO_SPECIALIZATION

/// @name specialization of vector arithmetic operations
/// @{

/// vector operator +
template <>
constexpr inline vec2f operator+(const vec2f& a) {
    return a;
}

/// vector operator -
template <>
constexpr inline vec2f operator-(const vec2f& a) {
    return {-a.x, -a.y};
}

/// vector operator +
template <>
constexpr inline vec2f operator+(const vec2f& a, const vec2f& b) {
    return {a.x + b.x, a.y + b.y};
}

/// vector operator -
template <>
constexpr inline vec2f operator-(const vec2f& a, const vec2f& b) {
    return {a.x - b.x, a.y - b.y};
}

/// vector operator *
template <>
constexpr inline vec2f operator*(const vec2f& a, const vec2f& b) {
    return {a.x * b.x, a.y * b.y};
}

/// vector operator *
template <>
constexpr inline vec2f operator*(const vec2f& a, const float b) {
    return {a.x * b, a.y * b};
}

/// vector operator *
template <>
constexpr inline vec2f operator*(const float a, const vec2f& b) {
    return {a * b.x, a * b.y};
}

/// vector operator /
template <>
constexpr inline vec2f operator/(const vec2f& a, const vec2f& b) {
    return {a.x / b.x, a.y / b.y};
}

/// vector operator /
template <>
constexpr inline vec2f operator/(const vec2f& a, const float b) {
    return {a.x / b, a.y / b};
}

/// vector operator /
template <>
constexpr inline vec2f operator/(const float a, const vec2f& b) {
    return {a / b.x, a / b.y};
}

/// vector operator +
template <>
constexpr inline vec3f operator+(const vec3f& a) {
    return a;
}

/// vector operator -
template <>
constexpr inline vec3f operator-(const vec3f& a) {
    return {-a.x, -a.y, -a.z};
}

/// vector operator +
template <>
constexpr inline vec3f operator+(const vec3f& a, const vec3f& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}

/// vector operator -
template <>
constexpr inline vec3f operator-(const vec3f& a, const vec3f& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}

/// vector operator *
template <>
constexpr inline vec3f operator*(const vec3f& a, const vec3f& b) {
    return {a.x * b.x, a.y * b.y, a.z * b.z};
}

/// vector operator *
template <>
constexpr inline vec3f operator*(const vec3f& a, const float b) {
    return {a.x * b, a.y * b, a.z * b};
}

/// vector operator *
template <>
constexpr inline vec3f operator*(const float a, const vec3f& b) {
    return {a * b.x, a * b.y, a * b.z};
}

/// vector operator /
template <>
constexpr inline vec3f operator/(const vec3f& a, const vec3f& b) {
    return {a.x / b.x, a.y / b.y, a.z / b.z};
}

/// vector operator /
template <>
constexpr inline vec3f operator/(const vec3f& a, const float b) {
    return {a.x / b, a.y / b, a.z / b};
}

/// vector operator /
template <>
constexpr inline vec3f operator/(const float a, const vec3f& b) {
    return {a / b.x, a / b.y, a / b.z};
}

/// vector operator +
template <>
constexpr inline vec4f operator+(const vec4f& a) {
    return a;
}

/// vector operator -
template <>
constexpr inline vec4f operator-(const vec4f& a) {
    return {-a.x, -a.y, -a.z, -a.w};
}

/// vector operator +
template <>
constexpr inline vec4f operator+(const vec4f& a, const vec4f& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}

/// vector operator -
template <>
constexpr inline vec4f operator-(const vec4f& a, const vec4f& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w};
}

/// vector operator *
template <>
constexpr inline vec4f operator*(const vec4f& a, const vec4f& b) {
    return {a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w};
}

/// vector operator *
template <>
constexpr inline vec4f operator*(const vec4f& a, const float b) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
}

/// vector operator *
template <>
constexpr inline vec4f operator*(const float a, const vec4f& b) {
    return {a * b.x, a * b.y, a * b.z, a * b.w};
}

/// vector operator /
template <>
constexpr inline vec4f operator/(const vec4f& a, const vec4f& b) {
    return {a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w};
}

/// vector operator /
template <>
constexpr inline vec4f operator/(const vec4f& a, const float b) {
    return {a.x / b, a.y / b, a.z / b, a.w / b};
}

/// vector operator /
template <>
constexpr inline vec4f operator/(const float a, const vec4f& b) {
    return {a / b.x, a / b.y, a / b.z, a / b.w};
}

/// @}

#endif

/// @name vector arithmetic assignment operations
/// @{

/// vector operator +=
template <typename T, int N>
constexpr inline vec<T, N>& operator+=(vec<T, N>& a, const vec<T, N>& b) {
    return a = a + b;
}

/// vector operator -=
template <typename T, int N>
constexpr inline vec<T, N>& operator-=(vec<T, N>& a, const vec<T, N>& b) {
    return a = a - b;
}

/// vector operator *=
template <typename T, int N>
constexpr inline vec<T, N>& operator*=(vec<T, N>& a, const vec<T, N>& b) {
    return a = a * b;
}

/// vector operator *=
template <typename T, int N>
constexpr inline vec<T, N>& operator*=(vec<T, N>& a, const T b) {
    return a = a * b;
}

/// vector operator /=
template <typename T, int N>
constexpr inline vec<T, N>& operator/=(vec<T, N>& a, const vec<T, N>& b) {
    return a = a / b;
}

/// vector operator /=
template <typename T, int N>
constexpr inline vec<T, N>& operator/=(vec<T, N>& a, const T b) {
    return a = a / b;
}

/// @}

/// @name vector product operations
/// @{

/// vector dot product
template <typename T, int N>
constexpr inline T dot(const vec<T, N>& a, const vec<T, N>& b) {
    auto c = T(0);
    for (auto i = 0; i < N; i++) c += a[i] * b[i];
    return c;
}

/// vector cross product (2d)
template <typename T>
constexpr inline T cross(const vec<T, 2>& a, const vec<T, 2>& b) {
    return a[0] * b[1] - a[1] * b[0];
}

/// vector cross product (3d)
template <typename T>
constexpr inline vec<T, 3> cross(const vec<T, 3>& a, const vec<T, 3>& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]};
}

/// @}

/// @name specialization of vector product operations
/// @{

/// vector dot product
template <>
constexpr inline float dot(const vec2f& a, const vec2f& b) {
    return a.x * b.x + a.y * b.y;
}

/// vector dot product
template <>
constexpr inline float dot(const vec3f& a, const vec3f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/// vector dot product
template <>
constexpr inline float dot(const vec4f& a, const vec4f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

/// vector cross product (2d)
template <>
constexpr inline float cross(const vec2f& a, const vec2f& b) {
    return a.x * b.y - a.y * b.x;
}

/// vector cross product (3d)
template <>
constexpr inline vec3f cross(const vec3f& a, const vec3f& b) {
    return {
        a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

/// @}

/// @name vector geometric operations
/// @{

/// vector length
template <typename T, int N>
constexpr inline T length(const vec<T, N>& a) {
    return std::sqrt(dot(a, a));
}

/// vector length squared
template <typename T, int N>
constexpr inline T lengthsqr(const vec<T, N>& a) {
    return dot(a, a);
}

/// vector normalization
template <typename T, int N>
constexpr inline vec<T, N> normalize(const vec<T, N>& a) {
    auto l = length(a);
    if (l == 0) return a;
    return a * (1 / l);
}

/// point distance
template <typename T, int N>
constexpr inline T dist(const vec<T, N>& a, const vec<T, N>& b) {
    return length(a - b);
}

/// point distance squared
template <typename T, int N>
constexpr inline T distsqr(const vec<T, N>& a, const vec<T, N>& b) {
    return lengthsqr(a - b);
}

/// angle between normalized vectors
template <typename T, int N>
constexpr inline T uangle(const vec<T, N>& a, const vec<T, N>& b) {
    auto d = dot(a, b);
    return d > 1 ? 0 : std::acos(d < -1 ? -1 : d);
}

/// angle between vectors
template <typename T, int N>
constexpr inline T angle(const vec<T, N>& a, const vec<T, N>& b) {
    return uangle(normalize(a), normalize(b));
}

/// vector linear interpolation
template <typename T, int N>
constexpr inline vec<T, N> lerp(const vec<T, N>& a, const vec<T, N>& b, T t) {
    return a * (1 - t) + b * t;
}

/// vector normalized linear interpolation
template <typename T, int N>
constexpr inline vec<T, N> nlerp(const vec<T, N>& a, const vec<T, N>& b, T t) {
    return normalize(lerp(a, b, t));
}

/// vector spherical linear interpolation (vectors have to be normalized)
template <typename T, int N>
constexpr inline vec<T, N> slerp(const vec<T, N>& a, const vec<T, N>& b, T t) {
    auto th = uangle(a, b);
    return th == 0 ? a :
                     a * (std::sin(th * (1 - t)) / std::sin(th)) +
                         b * (std::sin(th * t) / std::sin(th));
}

/// orthogonal vector
// http://lolengine.net/blog/2013/09/21/picking-orthogonal-vector-combing-coconuts)
template <typename T>
constexpr inline vec<T, 3> orthogonal(const vec<T, 3>& v) {
    return std::abs(v[0]) > std::abs(v[2]) ? vec<T, 3>{-v[1], v[0], 0} :
                                             vec<T, 3>{0, -v[2], v[1]};
}

/// orthonormalize two vectors
template <typename T>
constexpr inline vec<T, 3> orthonormalize(
    const vec<T, 3>& a, const vec<T, 3>& b) {
    return normalize(a - b * dot(a, b));
}

/// vector component-wise clamp
template <typename T, int N>
constexpr inline vec<T, N> clamp(
    const vec<T, N>& x, const T& min, const T& max) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = clamp(x[i], min, max);
    return c;
}

/// vector component-wise clamp
template <typename T, int N>
constexpr inline vec<T, N> clamp(
    const vec<T, N>& x, const vec<T, N>& min, const vec<T, N>& max) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = clamp(x[i], min[i], max[i]);
    return c;
}

/// clamp the length of a vector
template <typename T, int N, typename T1>
constexpr inline vec<T, N> clamplen(const vec<T, N> x, T1 max) {
    auto l = length(x);
    return (l > (T)max) ? x * (T)max / l : x;
}

/// @}

/// @name vector element comparison operations
/// @{

/// index of the min vector element
template <typename T, int N>
constexpr inline int min_element_idx(const vec<T, N>& a) {
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

/// index of the max vector element
template <typename T, int N>
constexpr inline int max_element_idx(const vec<T, N>& a) {
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

/// index of the min vector element
template <typename T, int N>
constexpr inline T min_element_val(const vec<T, N>& a) {
    auto v = std::numeric_limits<T>::max();
    for (auto i = 0; i < N; i++) {
        if (v > a[i]) v = a[i];
    }
    return v;
}

/// index of the max vector element
template <typename T, int N>
constexpr inline T max_element_val(const vec<T, N>& a) {
    auto v = -std::numeric_limits<T>::max();
    for (auto i = 0; i < N; i++) {
        if (v < a[i]) v = a[i];
    }
    return v;
}

// @}

// -----------------------------------------------------------------------------
// MATRICES
// -----------------------------------------------------------------------------

/// @name matrices
/// @{

///
/// Matrix of elements of compile time dimensions, stored in column major
/// format, with default initializer, and conversions to/from std::array.
/// Colums access via operator[]. Supports all std::array operations.
///
template <typename T, int N, int M>
struct mat {
    /// column data type
    using V = vec<T, N>;

    /// default constructor
    constexpr mat() {
        for (auto j = 0; j < M; j++) v[j] = V{};
    }

    /// list constructor
    constexpr mat(const std::initializer_list<V>& vv) {
        assert(M == vv.size());
        auto i = 0;
        for (auto&& e : vv) v[i++] = e;
    }

    /// conversion from std::array
    constexpr mat(const std::array<std::array<float, N>, M>& v) {
        *this = *(mat<T, N, M>*)&v;
    }

    /// conversion to std::array
    constexpr operator std::array<std::array<T, N>, M>() const {
        return *(std::array<std::array<T, N>, M>*)this;
    }

    /// conversion from flattened std::array
    constexpr mat(const std::array<float, N * M>& v) {
        *this = *(mat<T, N, M>*)&v;
    }

    /// conversion to flattened std::array
    constexpr operator std::array<T, N * M>() const {
        return *(std::array<T, N * M>*)this;
    }

    /// element access
    constexpr V& operator[](int i) { return v[i]; }
    /// element access
    constexpr const V& operator[](int i) const { return v[i]; }

    /// element data
    V v[M];
};

#ifndef YM_NO_SPECIALIZATION

///
/// Specialization for 2x2 float matrices.
///
template <>
struct mat<float, 2, 2> {
    /// size
    constexpr static const int N = 2, M = 2;
    /// type
    using T = float;

    /// column data type
    using V = vec<T, N>;

    /// default constructor
    constexpr mat() : x{0, 0}, y{0, 0} {}

    /// list constructor
    constexpr mat(const V& x, const V& y) : x(x), y(y) {}

    /// conversion from std::array
    constexpr mat(const std::array<std::array<T, N>, M>& vv)
        : x{vv[0]}, y{vv[1]} {}

    /// conversion to std::array
    constexpr operator std::array<std::array<T, N>, M>() const {
        return {x, y};
    }

    /// conversion from flattened std::array
    constexpr mat(const std::array<T, N * M>& vv)
        : x{vv[0], vv[1]}, y{vv[2], vv[3]} {}

    /// conversion to flattened std::array
    constexpr operator std::array<T, N * M>() const {
        return {x.x, x.y, y.x, y.y};
    }

    /// element access
    constexpr V& operator[](int i) { return v[i]; }
    /// element access
    constexpr const V& operator[](int i) const { return v[i]; }

    /// element data
    union {
        struct {
            V v[M];
        };
        struct {
            V x, y;
        };
    };
};

///
/// Specialization for 3x3 float matrices.
///
template <>
struct mat<float, 3, 3> {
    /// size
    constexpr static const int N = 3, M = 3;
    /// type
    using T = float;

    /// column data type
    using V = vec<T, N>;

    /// default constructor
    constexpr mat() : x{0, 0, 0}, y{0, 0, 0}, z{0, 0, 0} {}

    /// list constructor
    constexpr mat(const V& x, const V& y, const V& z) : x(x), y(y), z(z) {}

    /// conversion from std::array
    constexpr mat(const std::array<std::array<T, N>, M>& vv)
        : x{vv[0]}, y{vv[1]}, z{vv[2]} {}

    /// conversion to std::array
    constexpr operator std::array<std::array<T, N>, M>() const {
        return {x, y, z};
    }

    /// conversion from flattened std::array
    constexpr mat(const std::array<T, N * M>& vv)
        : x{vv[0], vv[1], vv[2]}
        , y{vv[3], vv[4], vv[5]}
        , z{vv[6], vv[7], vv[8]} {}

    /// conversion to flattened std::array
    constexpr operator std::array<T, N * M>() const {
        return {x.x, x.y, x.z, y.x, y.y, y.z, z.x, z.y, z.z};
    }

    /// element access
    constexpr V& operator[](int i) { return v[i]; }
    /// element access
    constexpr const V& operator[](int i) const { return v[i]; }

    /// element data
    union {
        struct {
            V v[M];
        };
        struct {
            V x, y, z;
        };
    };
};

///
/// Specialization for 4x4 float matrices.
///
template <>
struct mat<float, 4, 4> {
    /// size
    constexpr static const int N = 4, M = 4;
    /// type
    using T = float;

    /// column data type
    using V = vec<T, N>;

    /// default constructor
    constexpr mat()
        : x{0, 0, 0, 0}, y{0, 0, 0, 0}, z{0, 0, 0, 0}, w{0, 0, 0, 0} {}

    /// list constructor
    constexpr mat(const V& x, const V& y, const V& z, const V& w)
        : x(x), y(y), z(z), w(w) {}

    /// conversion from std::array
    constexpr mat(const std::array<std::array<T, N>, M>& vv)
        : x{vv[0]}, y{vv[1]}, z{vv[2]}, w{vv[3]} {}

    /// conversion to std::array
    constexpr operator std::array<std::array<T, N>, M>() const {
        return {x, y, z, w};
    }

    /// conversion from flattened std::array
    constexpr mat(const std::array<T, N * M>& vv)
        : x{vv[0], vv[1], vv[2], vv[3]}
        , y{vv[4], vv[5], vv[6], vv[7]}
        , z{vv[8], vv[9], vv[10], vv[11]}
        , w{vv[12], vv[13], vv[14], vv[15]} {}

    /// conversion to flattened std::array
    constexpr operator std::array<T, N * M>() const {
        return {x.x, x.y, x.z, x.w, y.x, y.y, y.z, y.w, z.x, z.y, z.z, z.w, w.x,
            w.y, w.z, w.w};
    }

    /// element access
    constexpr V& operator[](int i) { return v[i]; }
    /// element access
    constexpr const V& operator[](int i) const { return v[i]; }

    /// element data
    union {
        struct {
            V v[M];
        };
        struct {
            V x, y, z, w;
        };
    };
};

#endif

/// @}

/// @name matrix typedefs
/// @{

/// 1-dimensional float matrix
using mat1f = mat<float, 1, 1>;
/// 2-dimensional float matrix
using mat2f = mat<float, 2, 2>;
/// 3-dimensional float matrix
using mat3f = mat<float, 3, 3>;
/// 4-dimensional float matrix
using mat4f = mat<float, 4, 4>;

/// @}

/// @name matrix initializations

/// Initialize an identity matrix.
template <typename T, int N>
constexpr inline mat<T, N, N> identity_mat() {
    mat<T, N, N> c;
    for (auto j = 0; j < N; j++)
        for (auto i = 0; i < N; i++) c[j][i] = (i == j) ? 1 : 0;
    return c;
}

#ifndef YM_NO_SPECIALIZATION

/// Specialization for Initialize an identity matrix.
template <>
constexpr inline mat3f identity_mat() {
    return {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
}

/// Specialization for Initialize an identity matrix.
template <>
constexpr inline mat4f identity_mat() {
    return {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
}

#endif

/// @}

/// @name matrix constants
/// @{

/// 1-dimensional float identity matrix
const auto identity_mat1f = identity_mat<float, 1>();
/// 2-dimensional float identity matrix
const auto identity_mat2f = identity_mat<float, 2>();
/// 3-dimensional float identity matrix
const auto identity_mat3f = identity_mat<float, 3>();
/// 4-dimensional float identity matrix
const auto identity_mat4f = identity_mat<float, 4>();

/// @}

/// @name matrix iteration functions
/// @{

/// iteration support
template <typename T, int N, int M>
constexpr inline vec<T, N>* begin(mat<T, N, M>& a) {
    return a.v;
}

/// iteration support
template <typename T, int N, int M>
constexpr inline const vec<T, N>* begin(const mat<T, N, M>& a) {
    return a.v;
}

/// iteration support
template <typename T, int N, int M>
constexpr inline vec<T, N>* end(mat<T, N, M>& a) {
    return a.v + M;
}

/// iteration support
template <typename T, int N, int M>
constexpr inline const vec<T, N>* end(const mat<T, N, M>& a) {
    return a.v + M;
}

/// @}

/// @name matrix comparison operations
/// @{

/// vector operator ==
template <typename T, int N, int M>
constexpr inline bool operator==(const mat<T, N, M>& a, const mat<T, N, M>& b) {
    for (auto i = 0; i < M; i++)
        if (a[i] != b[i]) return false;
    return true;
}

/// vector operator !=
template <typename T, int N, int M>
constexpr inline bool operator!=(const mat<T, N, M>& a, const mat<T, N, M>& b) {
    return !(a == b);
}

/// @}

/// @name matrix arithmetic operations
/// @{

/// matrix operator -
template <typename T, int N, int M>
constexpr inline mat<T, M, N> operator-(const mat<T, N, M>& a) {
    mat<T, N, M> c;
    for (auto i = 0; i < M; i++) c[i] = -a[i];
    return c;
}

/// matrix operator +
template <typename T, int N, int M>
constexpr inline mat<T, M, N> operator+(
    const mat<T, N, M>& a, const mat<T, N, M>& b) {
    mat<T, N, M> c;
    for (auto i = 0; i < M; i++) c[i] = a[i] + b[i];
    return c;
}

/// matrix scalar multiply
template <typename T, int N, int M>
constexpr inline mat<T, M, N> operator*(const mat<T, N, M>& a, T b) {
    mat<T, N, M> c;
    for (auto i = 0; i < M; i++) c[i] = a[i] * b;
    return c;
}

/// matrix scalar division
template <typename T, int N, int M>
constexpr inline mat<T, M, N> operator/(const mat<T, N, M>& a, T b) {
    mat<T, N, M> c;
    for (auto i = 0; i < M; i++) c[i] = a[i] / b;
    return c;
}

/// matrix-vector right multiply
template <typename T, int N, int M>
constexpr inline vec<T, N> operator*(
    const mat<T, N, M>& a, const vec<T, M>& b) {
    vec<T, N> c = zero_vec<T, N>();
    for (auto j = 0; j < M; j++) c += a[j] * b[j];
    return c;
}

/// matrix-vector left multiply
template <typename T, int N, int M>
constexpr inline vec<T, M> operator*(
    const vec<T, N>& a, const mat<T, N, M>& b) {
    vec<T, M> c;
    for (auto j = 0; j < M; j++) c[j] = dot(a, b[j]);
    return c;
}

/// matrix-matrix multiply
template <typename T, int N, int M, int K>
constexpr inline mat<T, N, M> operator*(
    const mat<T, N, K>& a, const mat<T, K, M>& b) {
    mat<T, N, M> c;
    for (auto j = 0; j < M; j++) c[j] = a * b[j];
    return c;
}

/// @}

#ifndef YM_NO_SPECIALIZATION

/// @name specialization of matrix arithmetic operations
/// @{

/// matrix-vector right multiply
template <>
constexpr inline vec2f operator*(const mat2f& a, const vec2f& b) {
    return a.x * b.x + a.y * b.y;
}

/// matrix-vector left multiply
template <>
constexpr inline vec2f operator*(const vec2f& a, const mat2f& b) {
    return {dot(a, b.x), dot(a, b.y)};
}

/// matrix-matrix multiply
template <>
constexpr inline mat2f operator*(const mat2f& a, const mat2f& b) {
    return {a * b.x, a * b.y};
}

/// matrix-vector right multiply
template <>
constexpr inline vec3f operator*(const mat3f& a, const vec3f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

/// matrix-vector left multiply
template <>
constexpr inline vec3f operator*(const vec3f& a, const mat3f& b) {
    return {dot(a, b.x), dot(a, b.y), dot(a, b.z)};
}

/// matrix-matrix multiply
template <>
constexpr inline mat3f operator*(const mat3f& a, const mat3f& b) {
    return {a * b.x, a * b.y, a * b.z};
}

/// matrix-vector right multiply
template <>
constexpr inline vec4f operator*(const mat4f& a, const vec4f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}

/// matrix-vector left multiply
template <>
constexpr inline vec4f operator*(const vec4f& a, const mat4f& b) {
    return {dot(a, b.x), dot(a, b.y), dot(a, b.z), dot(a, b.w)};
}

/// matrix-matrix multiply
template <>
constexpr inline mat4f operator*(const mat4f& a, const mat4f& b) {
    return {a * b.x, a * b.y, a * b.z, a * b.w};
}

/// @}

#endif

/// @name matrix arithmetic assignment operations
/// @{

/// matrix sum assignment
template <typename T, int N, int M>
constexpr inline mat<T, M, N>& operator+=(
    mat<T, N, M>& a, const mat<T, N, M>& b) {
    return a = a + b;
}

/// matrix-matrix multiply assignment
template <typename T, int N, int M>
constexpr inline mat<T, M, N>& operator*=(
    mat<T, N, M>& a, const mat<T, N, M>& b) {
    return a = a * b;
}

/// matrix scaling assignment
template <typename T, int N, int M>
constexpr inline mat<T, M, N>& operator*=(mat<T, N, M>& a, const T& b) {
    return a = a * b;
}

/// matrix scaling assignment
template <typename T, int N, int M>
constexpr inline mat<T, M, N>& operator/=(mat<T, N, M>& a, const T& b) {
    return a = a / b;
}

/// @}

/// @name matrix algebra operations
/// @{

/// matrix diagonal
template <typename T, int N>
constexpr vec<T, N> mat_diagonal(const mat<T, N, N>& a) {
    vec<T, N> d;
    for (auto i = 0; i < N; i++) d[i] = a[i][i];
    return d;
}

/// matrix transpose
template <typename T, int N, int M>
constexpr inline mat<T, M, N> transpose(const mat<T, N, M>& a) {
    mat<T, M, N> c;
    for (auto j = 0; j < M; j++) {
        for (auto i = 0; i < N; i++) { c[i][j] = a[j][i]; }
    }
    return c;
}

/// matrix adjugate (2x2)
template <typename T>
constexpr inline mat<T, 2, 2> adjugate(const mat<T, 2, 2>& a) {
    return {{a[1][1], -a[0][1]}, {-a[1][0], a[0][0]}};
}

/// matrix adjugate (3x3)
template <typename T>
constexpr inline mat<T, 3, 3> adjugate(const mat<T, 3, 3>& a) {
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

/// matrix adjugate (4x4)
template <typename T>
constexpr inline mat<T, 4, 4> adjugate(const mat<T, 4, 4>& a) {
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

/// matrix determinant (2x2)
template <typename T>
constexpr inline T determinant(const mat<T, 2, 2>& a) {
    return a[0][0] * a[1][1] - a[0][1] * a[1][0];
}

/// matrix determinant (3x3)
template <typename T>
constexpr inline T determinant(const mat<T, 3, 3>& a) {
    return a[0][0] * (a[1][1] * a[2][2] - a[2][1] * a[1][2]) +
           a[0][1] * (a[1][2] * a[2][0] - a[2][2] * a[1][0]) +
           a[0][2] * (a[1][0] * a[2][1] - a[2][0] * a[1][1]);
}

/// matrix determinant (4x4)
template <typename T>
constexpr inline T determinant(const mat<T, 4, 4>& a) {
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

/// matrix inverse (uses adjugate and determinant)
template <typename T, int N>
constexpr inline mat<T, N, N> inverse(const mat<T, N, N>& a) {
    return adjugate(a) / determinant(a);
}

/// @}

// -----------------------------------------------------------------------------
// RIGID BODY TRANSFORMS/FRAMES
// -----------------------------------------------------------------------------

/// @name frames
/// @{

///
/// Rigid transforms stored as a column-major affine matrix Nx(N+1).
/// In memory, this representation is equivalent to storing an NxN rotation
/// followed by a Nx1 translation. Viewed this way, the representation allows
/// also to retrive the axis of the coordinate frame as the first N column and
/// the translation as the N+1 column.
/// Colums access via operator[]. Access rotation and position with pos() and
/// rot().
///
template <typename T, int N>
struct frame {
    /// column data type
    using V = vec<T, N>;
    /// rotation data type
    using M = mat<T, N, N>;

    /// default constructor
    constexpr frame() {
        for (auto i = 0; i < N + 1; i++) v[i] = V{};
    }

    /// element constructor
    constexpr frame(const std::initializer_list<vec<T, N>>& vv) {
        assert(N + 1 == vv.size());
        auto i = 0;
        for (auto&& e : vv) v[i++] = e;
    }

    /// element constructor
    constexpr frame(const M& m, const V& t) {
        for (auto i = 0; i < N; i++) v[i] = m[i];
        v[N + 1] = t;
    }

    /// conversion from std::array
    constexpr frame(const std::array<std::array<float, N>, N + 1>& v) {
        *this = *(frame<T, N>*)&v;
    }

    /// conversion to std::array
    constexpr operator std::array<std::array<T, N>, N + 1>() const {
        return *(std::array<std::array<T, N>, N + 1>*)this;
    }

    /// conversion from flattened std::array
    constexpr frame(const std::array<float, N*(N + 1)>& v) {
        *this = *(frame<T, N>*)&v;
    }

    /// conversion to flattened std::array
    constexpr operator std::array<T, N*(N + 1)>() const {
        return *(std::array<T, N*(N + 1)>*)this;
    }

    /// element access
    constexpr V& operator[](int i) { return v[i]; }
    /// element access
    constexpr const V& operator[](int i) const { return v[i]; }

    /// access position
    constexpr V& pos() { return v[N]; }
    /// access position
    constexpr const V& pos() const { return v[N]; }

    /// access rotation
    constexpr M& rot() { return *(M*)v; }
    /// access rotation
    constexpr const M& rot() const { return *(M*)v; }

    /// element data
    V v[N + 1];
};

#ifndef YM_NO_SPECIALIZATION

///
/// Specialization for 3D float frames.
///
template <>
struct frame<float, 2> {
    /// size
    constexpr static const int N = 2;
    /// type
    using T = float;

    /// column data type
    using V = vec<T, N>;
    /// rotation data type
    using M = mat<T, N, N>;

    /// default constructor
    constexpr frame() : x{0, 0}, y{0, 0}, o{0, 0} {}

    /// element constructor
    constexpr frame(const V& x, const V& y, const V& o) : x(x), y(y), o(o) {}

    /// element constructor
    constexpr frame(const M& m, const V& t) : m(m), t(t) {}

    /// conversion from std::array
    constexpr frame(const std::array<std::array<float, N>, N + 1>& vv)
        : x(vv[0]), y(vv[1]), o(vv[2]) {}

    /// conversion to std::array
    constexpr operator std::array<std::array<T, N>, N + 1>() const {
        return {x, y, o};
    }

    /// conversion from flattened std::array
    constexpr frame(const std::array<float, N*(N + 1)>& vv)
        : x{vv[0], vv[1]}, y{vv[2], vv[3]}, o{vv[4], vv[5]} {}

    /// conversion to flattened std::array
    constexpr operator std::array<T, N*(N + 1)>() const {
        return {x.x, x.y, y.x, y.y, o.x, o.y};
    }

    /// element access
    constexpr V& operator[](int i) { return v[i]; }
    /// element access
    constexpr const V& operator[](int i) const { return v[i]; }

    /// access position
    constexpr V& pos() { return v[N]; }
    /// access position
    constexpr const V& pos() const { return v[N]; }

    /// access rotation
    constexpr M& rot() { return *(M*)v; }
    /// access rotation
    constexpr const M& rot() const { return *(M*)v; }

    /// element data
    union {
        struct {
            V v[N + 1];
        };
        struct {
            V x, y, o;
        };
        struct {
            M m;
            V t;
        };
    };
};

///
/// Specialization for 3D float frames.
///
template <>
struct frame<float, 3> {
    /// size
    constexpr static const int N = 3;
    /// type
    using T = float;

    /// column data type
    using V = vec<T, N>;
    /// rotation data type
    using M = mat<T, N, N>;

    /// default constructor
    constexpr frame() : x{0, 0, 0}, y{0, 0, 0}, z{0, 0, 0}, o{0, 0, 0} {}

    /// element constructor
    constexpr frame(const V& x, const V& y, const V& z, const V& o)
        : x(x), y(y), z(z), o(o) {}

    /// element constructor
    constexpr frame(const M& m, const V& t) : m(m), t(t) {}

    /// conversion from std::array
    constexpr frame(const std::array<std::array<float, N>, N + 1>& vv)
        : x(vv[0]), y(vv[1]), z(vv[2]), o(vv[3]) {}

    /// conversion to std::array
    constexpr operator std::array<std::array<T, N>, N + 1>() const {
        return {x, y, z, o};
    }

    /// conversion from flattened std::array
    constexpr frame(const std::array<float, N*(N + 1)>& vv)
        : x{vv[0], vv[1], vv[2]}
        , y{vv[3], vv[4], vv[5]}
        , z{vv[6], vv[7], vv[8]}
        , o{vv[9], vv[10], vv[11]} {}

    /// conversion to flattened std::array
    constexpr operator std::array<T, N*(N + 1)>() const {
        return {x.x, x.y, x.z, y.x, y.y, y.z, z.x, z.y, z.z, o.x, o.y, o.z};
    }

    /// element access
    constexpr V& operator[](int i) { return v[i]; }
    /// element access
    constexpr const V& operator[](int i) const { return v[i]; }

    /// access position
    constexpr V& pos() { return v[N]; }
    /// access position
    constexpr const V& pos() const { return v[N]; }

    /// access rotation
    constexpr M& rot() { return *(M*)v; }
    /// access rotation
    constexpr const M& rot() const { return *(M*)v; }

    /// element data
    union {
        struct {
            V v[N + 1];
        };
        struct {
            V x, y, z, o;
        };
        struct {
            M m;
            V t;
        };
    };
};

#endif

/// @}

/// @name frame typedefs
/// @{

/// 1-dimensional float frame
using frame1f = frame<float, 1>;
/// 2-dimensional float frame
using frame2f = frame<float, 2>;
/// 3-dimensional float frame
using frame3f = frame<float, 3>;
/// 4-dimensional float frame
using frame4f = frame<float, 4>;

/// @}

/// @name frame initializations
/// @{

/// Initialize an identity frame.
template <typename T, int N>
constexpr inline frame<T, N> identity_frame() {
    frame<T, N> c;
    for (auto j = 0; j < N; j++)
        for (auto i = 0; i < N; i++) c[j][i] = (i == j) ? 1 : 0;
    for (auto i = 0; i < N; i++) c[N][i] = 0;
    return c;
}

// initializes a frame3 from origin and z.
template <typename T>
constexpr inline frame<T, 3> make_frame3_fromz(
    const vec<T, 3>& o, const vec<T, 3>& z_) {
    auto z = normalize(z_);
    auto x = normalize(orthogonal(z));
    auto y = normalize(cross(z, x));
    return {x, y, z, o};
}

// initializes a frame3 from origin, z and x.
template <typename T>
constexpr inline frame<T, 3> make_frame3_fromzx(
    const vec<T, 3>& o, const vec<T, 3>& z_, const vec<T, 3>& x_) {
    auto z = normalize(z_);
    auto x = orthonormalize(x_, z);
    auto y = normalize(cross(z, x));
    return {x, y, z, o};
}

/// @}

#ifndef YM_NO_SPECIALIZATION

/// @name specialization for frame initializations
/// @{

/// Initialize an identity frame.
template <>
constexpr inline frame2f identity_frame() {
    return {{1, 0}, {0, 1}, {0, 0}};
}

/// Initialize an identity frame.
template <>
constexpr inline frame3f identity_frame() {
    return {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};
}

// initializes a frame3 from origin and z.
template <>
constexpr inline frame3f make_frame3_fromz(const vec3f& o, const vec3f& z_) {
    auto z = normalize(z_);
    auto x = normalize(orthogonal(z));
    auto y = normalize(cross(z, x));
    return {x, y, z, o};
}

// initializes a frame3 from origin, z and x.
template <typename T>
constexpr inline frame3f make_frame3_fromzx(
    const vec3f& o, const vec3f& z_, const vec3f& x_) {
    auto z = normalize(z_);
    auto x = orthonormalize(x_, z);
    auto y = normalize(cross(z, x));
    return {x, y, z, o};
}

/// @}

#endif

/// @name frame constants
/// @{

/// 1-dimensional float identity frame
const auto identity_frame1f = identity_frame<float, 1>();
/// 2-dimensional float identity frame
const auto identity_frame2f = identity_frame<float, 2>();
/// 3-dimensional float identity frame
const auto identity_frame3f = identity_frame<float, 3>();
/// 4-dimensional float identity frame
const auto identity_frame4f = identity_frame<float, 4>();

/// @}

/// @name frame element access functions
/// @{

/// frame position const access
template <typename T, int N>
constexpr inline const vec<T, N>& pos(const frame<T, N>& f) {
    return f[N];
}

/// frame rotation const access
template <typename T, int N>
constexpr inline const mat<T, N, N>& rot(const frame<T, N>& f) {
    return *(mat<T, N, N>*)&f;
}

/// frame position reference
template <typename T, int N>
constexpr inline vec<T, N>& pos(frame<T, N>& f) {
    return f[N];
}

/// frame rotation reference
template <typename T, int N>
constexpr inline mat<T, N, N>& rot(frame<T, N>& f) {
    return *(mat<T, N, N>*)&f;
}

/// @}

/// @name frame iteration functions
/// @{

/// iteration support
template <typename T, int N>
constexpr inline vec<T, N>* begin(frame<T, N>& a) {
    return a.v;
}

/// iteration support
template <typename T, int N>
constexpr inline const vec<T, N>* begin(const frame<T, N>& a) {
    return a.v;
}

/// iteration support
template <typename T, int N>
constexpr inline vec<T, N>* end(frame<T, N>& a) {
    return a.v + N + 1;
}

/// iteration support
template <typename T, int N>
constexpr inline const vec<T, N>* end(const frame<T, N>& a) {
    return a.v + N + 1;
}

/// @}

/// @name frame conversion operations
/// @{

/// frame to matrix conversion
template <typename T, int N>
constexpr inline mat<T, N + 1, N + 1> to_mat(const frame<T, N>& a) {
    auto m = mat<T, N + 1, N + 1>();
    for (auto j = 0; j < N; j++) {
        (vec<T, N>&)m[j] = a[j];
        m[j][N] = 0;
    }
    (vec<T, N>&)m[N] = a[N];
    m[N][N] = 1;
    return m;
}

/// matrix to frame conversion
template <typename T, int N>
constexpr inline frame<T, N - 1> to_frame(const mat<T, N, N>& a) {
    auto f = frame<T, N - 1>();
    for (auto j = 0; j < N; j++) {
        for (auto i = 0; i < N - 1; i++) { f[j][i] = a[j][i]; }
    }
    return f;
}

/// @}

/// @name frame comparison operations
/// @{

/// vector operator ==
template <typename T, int N>
constexpr inline bool operator==(const frame<T, N>& a, const frame<T, N>& b) {
    for (auto i = 0; i < N + 1; i++)
        if (a[i] != b[i]) return false;
    return true;
}

/// vector operator !=
template <typename T, int N>
constexpr inline bool operator!=(const frame<T, N>& a, const frame<T, N>& b) {
    return !(a == b);
}

/// @}

/// @name frame algebra operations
/// @{

/// frame composition (equivalent to affine matrix multiply)
template <typename T, int N>
constexpr inline frame<T, N> operator*(
    const frame<T, N>& a, const frame<T, N>& b) {
    return make_frame(
        ym::rot(a) * ym::rot(b), ym::rot(a) * ym::pos(b) + ym::pos(a));
}

/// frame inverse (equivalent to rigid affine inverse)
template <typename T, int N>
constexpr inline frame<T, N> inverse(const frame<T, N>& a) {
    auto minv = transpose(ym::rot(a));
    return make_frame(minv, -(minv * ym::pos(a)));
}

/// @}

#ifndef YM_NO_SPECIALIZATION

/// @name specialization of frame algebra operations
/// @{

/// frame composition (equivalent to affine matrix multiply)
template <>
constexpr inline frame3f operator*(const frame3f& a, const frame3f& b) {
    return {a.m * b.m, a.m * b.t + a.t};
}

/// frame inverse (equivalent to rigid affine inverse)
template <>
constexpr inline frame3f inverse(const frame3f& a) {
    auto minv = transpose(a.m);
    return {minv, -(minv * a.t)};
}

/// @}

#endif

// -----------------------------------------------------------------------------
// QUATERNIONS
// -----------------------------------------------------------------------------

/// @name quaterions
/// @{

///
/// Quaternions implemented as a vec<T,4>. Data access via operator[].
/// Quaterions are xi + yj + zk + w.
///
template <typename T>
struct quat {
    /// default constructor
    constexpr quat() : std::array<T, 4>{} {}

    // list constructor
    constexpr quat(const T& x, const T& y, const T& z, const T& w)
        : v{x, y, z, w} {}

    /// element access
    constexpr T& operator[](int i) { return v[i]; }
    constexpr const T& operator[](int i) const { return v[i]; }

    /// iteration
    constexpr T* begin() { return v; };
    /// iteration
    constexpr const T* begin() const { return v; };
    /// iteration
    constexpr T* end() { return v + 4; };
    /// iteration
    constexpr const T* end() const { return v + 4; };

   private:
    T v[4];
};

/// @}

/// @name quaterion typedefs
/// @{

/// float quaterion
using quat4f = quat<float>;

/// @}

/// @name quaterion element access
/// @{

/// quaternion angle const access
template <typename T>
constexpr inline T angle(const quat<T>& a) {
    return std::acos(a[3]) * 2;
}

/// quaternion axis const access
template <typename T>
constexpr inline vec<T, 3> axis(const quat<T>& a) {
    return {a[0], a[1], a[2]};
}

/// @}

/// @name quaterion conversion operations
/// @{

/// quaterion to matrix conversion
template <typename T>
constexpr inline mat<T, 4, 4> quat_to_mat(const quat<T>& v) {
    return {
        {v[3] * v[3] + v[0] * v[0] - v[1] * v[1] - v[2] * v[2],
            (v[0] * v[1] + v[2] * v[3]) * 2, (v[2] * v[0] - v[1] * v[3]) * 2},
        {(v[0] * v[1] - v[2] * v[3]) * 2,
            v[3] * v[3] - v[0] * v[0] + v[1] * v[1] - v[2] * v[2],
            (v[1] * v[2] + v[0] * v[3]) * 2},
        {(v[2] * v[0] + v[1] * v[3]) * 2, (v[1] * v[2] - v[0] * v[3]) * 2,
            v[3] * v[3] - v[0] * v[0] - v[1] * v[1] + v[2] * v[2]}};
}

/// @}

/// @name quaterion algebra operations
/// @{

/// quaterion multiply
template <typename T>
constexpr quat<T> operator*(const quat<T>& a, const quat<T>& b) {
    return {a[0] * b[3] + a[3] * b[0] + a[1] * b[3] - a[2] * b[1],
        a[1] * b[3] + a[3] * b[1] + a[2] * b[0] - a[0] * b[2],
        a[2] * b[3] + a[3] * b[2] + a[0] * b[1] - a[1] * b[0],
        a[3] * b[3] - a[0] * b[0] - a[1] * b[1] - a[2] * b[2]};
}

/// quaterion conjugate
template <typename T>
constexpr quat<T> conjugate(const quat<T>& v) {
    return {-v[0], -v[1], -v[2], v[3]};
}

/// quaterion inverse
template <typename T>
constexpr quat<T> inverse(const quat<T>& v) {
    return qconj(v) / lengthsqr(vec<T, 4>(v));
}

/// quaterion normalized linear interpolation
template <typename T>
constexpr quat<T> nlerp(const quat<T>& a, const quat<T>& b, T t) {
    return nlerp(vec<T, 4>(a),
        dot(vec<T, 4>(a), vec<T, 4>(b)) < 0 ? -vec<T, 4>(b) : vec<T, 4>(b), t);
}

/// quaterion spherical linear interpolation
template <typename T>
constexpr quat<T> slerp(const quat<T>& a, const quat<T>& b, T t) {
    return slerp(vec<T, 4>(a),
        dot(vec<T, 4>(a), vec<T, 4>(b)) < 0 ? -vec<T, 4>(b) : vec<T, 4>(b), t);
}

/// @}

// -----------------------------------------------------------------------------
// AXIS ALIGNED BOUNDING BOXES
// -----------------------------------------------------------------------------

/// @name axis-aligned bounding boxes
/// @{

///
/// Axis aligned bounding box represented as a min/max vector pair.
/// Access min/max with operator[].
///
template <typename T, int N>
struct bbox {
    /// column data type
    using V = vec<T, N>;

    /// initializes an invalid bbox
    constexpr bbox() {
        for (auto i = 0; i < N; i++) {
            v[0][i] = std::numeric_limits<T>::max();
            v[1][i] = std::numeric_limits<T>::lowest();
        }
    }

    /// list constructor
    constexpr bbox(const vec<T, N>& m, const vec<T, N>& M) : v{m, M} {}

    /// element access
    constexpr V& operator[](int i) { return v[i]; }
    constexpr const V& operator[](int i) const { return v[i]; }

    /// element data
    V v[2];
};

#ifndef YM_NO_SPECIALIZATION

///
/// Specialization for float 3D bounding boxes.
///
template <>
struct bbox<float, 1> {
    /// size
    constexpr static const int N = 1;
    /// type
    using T = float;

    /// column data type
    using V = vec<T, N>;

    /// initializes an invalid bbox
    constexpr bbox() : min{flt_max}, max{flt_min} {}
    /// list constructor
    constexpr bbox(const vec<T, N>& m, const vec<T, N>& M) : min{m}, max{M} {}

    /// element access
    constexpr V& operator[](int i) { return v[i]; }
    /// element access
    constexpr const V& operator[](int i) const { return v[i]; }

    /// element data
    union {
        struct {
            V v[2];
        };
        struct {
            V min, max;
        };
    };
};

///
/// Specialization for float 3D bounding boxes.
///
template <>
struct bbox<float, 2> {
    /// size
    constexpr static const int N = 2;
    /// type
    using T = float;

    /// column data type
    using V = vec<T, N>;

    /// initializes an invalid bbox
    constexpr bbox() : min{flt_max, flt_max}, max{flt_min, flt_min} {}
    /// list constructor
    constexpr bbox(const vec<T, N>& m, const vec<T, N>& M) : min{m}, max{M} {}

    /// element access
    constexpr V& operator[](int i) { return v[i]; }
    /// element access
    constexpr const V& operator[](int i) const { return v[i]; }

    /// element data
    union {
        struct {
            V v[2];
        };
        struct {
            V min, max;
        };
    };
};

///
/// Specialization for float 3D bounding boxes.
///
template <>
struct bbox<float, 3> {
    /// size
    constexpr static const int N = 3;
    /// type
    using T = float;

    /// column data type
    using V = vec<T, N>;

    /// initializes an invalid bbox
    constexpr bbox()
        : min{flt_max, flt_max, flt_max}, max{flt_min, flt_min, flt_min} {}
    /// list constructor
    constexpr bbox(const vec<T, N>& m, const vec<T, N>& M) : min{m}, max{M} {}

    /// element access
    constexpr V& operator[](int i) { return v[i]; }
    /// element access
    constexpr const V& operator[](int i) const { return v[i]; }

    /// element data
    union {
        struct {
            V v[2];
        };
        struct {
            V min, max;
        };
    };
};

///
/// Specialization for float 3D bounding boxes.
///
template <>
struct bbox<float, 4> {
    /// size
    constexpr static const int N = 4;
    /// type
    using T = float;

    /// column data type
    using V = vec<T, N>;

    /// initializes an invalid bbox
    constexpr bbox()
        : min{flt_max, flt_max, flt_max, flt_max}
        , max{flt_min, flt_min, flt_min, flt_min} {}
    /// list constructor
    constexpr bbox(const vec<T, N>& m, const vec<T, N>& M) : min{m}, max{M} {}

    /// element access
    constexpr V& operator[](int i) { return v[i]; }
    /// element access
    constexpr const V& operator[](int i) const { return v[i]; }

    /// element data
    union {
        struct {
            V v[2];
        };
        struct {
            V min, max;
        };
    };
};

#endif

/// @}

/// @name axis-aligned bounding box typedefs
/// @{

/// 1-dimensional float bbox
using bbox1f = bbox<float, 1>;
/// 2-dimensional float bbox
using bbox2f = bbox<float, 2>;
/// 3-dimensional float bbox
using bbox3f = bbox<float, 3>;
/// 4-dimensional float bbox
using bbox4f = bbox<float, 4>;

/// @}

/// @name axis-aligned bounding box initializations
/// @{

/// initializes an empty bbox
template <typename T, int N>
constexpr inline bbox<T, N> invalid_bbox() {
    auto a = bbox<T, N>();
    for (auto i = 0; i < N; i++) {
        a[0][i] = std::numeric_limits<T>::max();
        a[1][i] = std::numeric_limits<T>::lowest();
    }
    return a;
}

/// initialize a bonding box from a list of points
template <typename T, int N>
constexpr inline bbox<T, N> make_bbox(
    const std::initializer_list<vec<T, N>>& v) {
    auto a = invalid_bbox<T, N>();
    for (auto&& vv : v) {
        for (auto i = 0; i < N; i++) {
            a[0][i] = min(a[0][i], vv[i]);
            a[1][i] = max(a[1][i], vv[i]);
        }
    }
    return a;
}

/// @}

/// @name axis-aligned bounding box constants
/// @{

/// 1-dimensional float empty bbox
const auto invalid_bbox1f = bbox1f();
/// 2-dimensional float empty bbox
const auto invalid_bbox2f = bbox2f();
/// 3-dimensional float empty bbox
const auto invalid_bbox3f = bbox3f();
/// 4-dimensional float empty bbox
const auto invalid_bbox4f = bbox4f();

/// @}

/// @name axis-aligned bounding element access
/// @{

/// computes the center of a bbox
template <typename T, int N>
constexpr inline vec<T, N> center(const bbox<T, N>& a) {
    return (a[0] + a[1]) / (T)2;
}

/// computes the diagonal of a bbox
template <typename T, int N>
constexpr inline vec<T, N> diagonal(const bbox<T, N>& a) {
    return a[1] - a[0];
}

/// @}

/// @name frame iteration functions
/// @{

/// iteration support
template <typename T, int N>
constexpr inline vec<T, N>* begin(bbox<T, N>& a) {
    return a.v;
}

/// iteration support
template <typename T, int N>
constexpr inline const vec<T, N>* begin(const bbox<T, N>& a) {
    return a.v;
}

/// iteration support
template <typename T, int N>
constexpr inline vec<T, N>* end(bbox<T, N>& a) {
    return a.v + 2;
}

/// iteration support
template <typename T, int N>
constexpr inline const vec<T, N>* end(const bbox<T, N>& a) {
    return a.v + 2;
}

/// @}

/// @name axis-aligned bounding box operations
/// @{

/// expands a bounding box with a point
template <typename T, int N>
constexpr inline bbox<T, N> expand(const bbox<T, N>& a, const vec<T, N>& b) {
    bbox<T, N> c;
    for (auto i = 0; i < N; i++) {
        c[0][i] = min(a[0][i], b[i]);
        c[1][i] = max(a[1][i], b[i]);
    }
    return c;
}

/// expands a bounding box with a bounding box
template <typename T, int N>
constexpr inline bbox<T, N> expand(const bbox<T, N>& a, const bbox<T, N>& b) {
    bbox<T, N> c;
    for (auto i = 0; i < N; i++) {
        c[0][i] = min(a[0][i], b[0][i]);
        c[1][i] = max(a[1][i], b[1][i]);
    }
    return c;
}

/// check if a bounding box contains a point
template <typename T, int N>
constexpr inline bool contains(const bbox<T, N>& a, const vec<T, N>& b) {
    for (auto i = 0; i < N; i++) {
        if (a[0][i] > b[i] || a[1][i] < b[i]) return false;
    }
    return true;
}

/// check if a bounding box contains a bounding box
template <typename T, int N>
constexpr inline bool contains(const bbox<T, N>& a, const bbox<T, N>& b) {
    for (auto i = 0; i < N; i++) {
        if (a[0][i] > b[1][i] || a[1][i] < b[0][i]) return false;
    }
    return true;
}

/// @}

#ifndef YM_NO_SPECIALIZATION

/// @name specialization of axis-aligned bounding box operations
/// @{

/// expands a bounding box with a point
template <>
constexpr inline bbox3f expand(const bbox3f& a, const vec3f& b) {
    return {{min(a.min.x, b.x), min(a.min.y, b.y), min(a.min.z, b.z)},
        {max(a.max.x, b.x), max(a.max.y, b.y), max(a.max.z, b.z)}};
}

/// expands a bounding box with a bounding box
template <>
constexpr inline bbox3f expand(const bbox3f& a, const bbox3f& b) {
    return {
        {min(a.min.x, b.min.x), min(a.min.y, b.min.y), min(a.min.z, b.min.z)},
        {max(a.max.x, b.max.x), max(a.max.y, b.max.y), max(a.max.z, b.max.z)}};
}

/// check if a bounding box contains a point
template <>
constexpr inline bool contains(const bbox3f& a, const vec3f& b) {
    if (a.min.x > b.x || a.max.x < b.x) return false;
    if (a.min.y > b.y || a.max.y < b.y) return false;
    if (a.min.z > b.z || a.max.z < b.z) return false;
    return true;
}

/// check if a bounding box contains a bounding box
template <>
constexpr inline bool contains(const bbox3f& a, const bbox3f& b) {
    if (a.min.x > b.max.x || a.max.x < b.min.x) return false;
    if (a.min.y > b.max.y || a.max.y < b.min.y) return false;
    if (a.min.z > b.max.z || a.max.z < b.min.z) return false;
    return true;
}

/// @}

#endif

/// @name axis-aligned bounding box operations
/// @{

/// same as expand()
template <typename T, int N>
constexpr inline bbox<T, N> operator+(const bbox<T, N>& a, const T& b) {
    return expand(a, b);
}

/// same as expand()
template <typename T, int N>
constexpr inline bbox<T, N> operator+(
    const bbox<T, N>& a, const bbox<T, N>& b) {
    return expand(a, b);
}

/// assign to expand()
template <typename T, int N>
constexpr inline bbox<T, N>& operator+=(bbox<T, N>& a, const vec<T, N>& b) {
    return a = expand(a, b);
}

/// assign to expand()
template <typename T, int N>
constexpr inline bbox<T, N>& operator+=(bbox<T, N>& a, const bbox<T, N>& b) {
    return a = expand(a, b);
}

/// @}

// -----------------------------------------------------------------------------
// RAYS
// -----------------------------------------------------------------------------

/// @name rays
/// @{

///
/// Rays with origin, direction and min/max t value.
///
template <typename T, int N>
struct ray {
    /// origin
    vec<T, N> o;
    /// direction
    vec<T, N> d;
    /// minimum distance
    T tmin;
    /// maximum distance
    T tmax;

    /// default constructor
    constexpr ray()
        : o(), d(0, 0, 1), tmin(0), tmax(std::numeric_limits<T>::max()) {}
    /// initializes a ray from its elements
    constexpr ray(const vec<T, N>& o, const vec<T, N>& d, T tmin = 0,
        T tmax = std::numeric_limits<T>::max())
        : o(o), d(d), tmin(tmin), tmax(tmax) {}
};

#ifndef YM_NO_SPECIALIZATION

///
/// Sepcialization for 3D float rays.
///
template <>
struct ray<float, 3> {
    /// size
    constexpr static const int N = 3;
    /// type
    using T = float;

    /// origin
    vec<T, N> o;
    /// direction
    vec<T, N> d;
    /// minimum distance
    T tmin;
    /// maximum distance
    T tmax;

    /// default constructor
    constexpr ray() : o{0, 0, 0}, d{0, 0, 1}, tmin{0}, tmax{flt_max} {}
    /// initializes a ray from its elements
    constexpr ray(
        const vec<T, N>& o, const vec<T, N>& d, T tmin = 0, T tmax = flt_max)
        : o(o), d(d), tmin(tmin), tmax(tmax) {}
};

#endif

/// @}

/// @name ray typedefs
/// @{

/// 1-dimensional float ray
using ray1f = ray<float, 1>;
/// 2-dimensional float ray
using ray2f = ray<float, 2>;
/// 3-dimensional float ray
using ray3f = ray<float, 3>;
/// 4-dimensional float ray
using ray4f = ray<float, 4>;

/// @}

/// @name ray operations
/// @{

/// evalutes the position along the ray
template <typename T, int N>
constexpr inline vec<T, N> eval(const ray<T, N>& ray, T t) {
    return ray.o + t * ray.d;
}

/// @}

// -----------------------------------------------------------------------------
// TRANSFORMS
// -----------------------------------------------------------------------------

/// @name transform operations
/// @{

/// transforms a point by a matrix
template <typename T, int N>
constexpr inline vec<T, N> transform_point(
    const mat<T, N + 1, N + 1>& a, const vec<T, N>& b) {
    // make it generic
    auto vb = vec<T, N + 1>();
    (vec<T, N>&)vb = b;
    vb[N] = 1;
    auto tvb = a * vb;
    return *(vec<T, N>*)(&tvb) / tvb[N];
}

/// transforms a vector by a matrix
template <typename T, int N>
constexpr inline vec<T, N> transform_vector(
    const mat<T, N + 1, N + 1>& a, const vec<T, N>& b) {
    // make it generic
    auto vb = vec<T, N + 1>();
    (vec<T, N>&)vb = b;
    vb[N] = 0;
    auto tvb = a * vb;
    return *(vec<T, N>*)(&tvb);
}

/// transforms a direction by a matrix
template <typename T, int N>
constexpr inline vec<T, N> transform_direction(
    const mat<T, N + 1, N + 1>& a, const vec<T, N>& b) {
    return normalize(transform_vector(a, b));
}

/// transforms a point by a frame (rigid affine transform)
template <typename T, int N>
constexpr inline vec<T, N> transform_point(
    const frame<T, N>& a, const vec<T, N>& b) {
    return ym::rot(a) * b + ym::pos(a);
}

/// transforms a vector by a frame (rigid affine transform)
template <typename T, int N>
constexpr inline vec<T, N> transform_vector(
    const frame<T, N>& a, const vec<T, N>& b) {
    return ym::rot(a) * b;
}

/// transforms a direction by a frame (rigid affine transform)
template <typename T, int N>
constexpr inline vec<T, N> transform_direction(
    const frame<T, N>& a, const vec<T, N>& b) {
    return ym::rot(a) * b;
}

/// transforms a frame by a frame (rigid affine transform)
template <typename T, int N>
constexpr inline frame<T, N> transform_frame(
    const frame<T, N>& a, const frame<T, N>& b) {
    return {ym::rot(a) * ym::rot(b), ym::pos(a) * ym::pos(b) + ym::pos(a)};
}

/// inverse transforms a point by a frame (rigid affine transform)
template <typename T, int N>
constexpr inline vec<T, N> transform_point_inverse(
    const frame<T, N>& a, const vec<T, N>& b) {
    return (b - ym::pos(a)) * ym::rot(a);
}

/// inverse transforms a vector by a frame (rigid affine transform)
template <typename T, int N>
constexpr inline vec<T, N> transform_vector_inverse(
    const frame<T, N>& a, const vec<T, N>& b) {
    return b * ym::rot(a);
}

/// inverse transforms a direction by a frame (rigid affine transform)
template <typename T, int N>
constexpr inline vec<T, N> transform_direction_inverse(
    const frame<T, N>& a, const vec<T, N>& b) {
    return b * ym::rot(a);
}

/// @}

#ifndef YM_NO_SPECIALIZATION

/// @name specialization of transform operations
/// @{

/// transforms a point by a matrix
template <>
constexpr inline vec3f transform_point(const mat4f& a, const vec3f& b) {
    auto vb = vec4f{b.x, b.y, b.z, 1};
    auto tvb = a * vb;
    return vec3f{tvb.x, tvb.y, tvb.z} / tvb.w;
}

/// transforms a vector by a matrix
template <>
constexpr inline vec3f transform_vector(const mat4f& a, const vec3f& b) {
    auto vb = vec4f{b.x, b.y, b.z, 0};
    auto tvb = a * vb;
    return vec3f{tvb.x, tvb.y, tvb.z} / tvb.w;
}

/// transforms a direction by a matrix
template <>
constexpr inline vec3f transform_direction(const mat4f& a, const vec3f& b) {
    return normalize(transform_vector(a, b));
}

/// transforms a point by a frame (rigid affine transform)
template <>
constexpr inline vec3f transform_point(const frame3f& a, const vec3f& b) {
    return a.m * b + a.t;
}

/// transforms a vector by a frame (rigid affine transform)
template <>
constexpr inline vec3f transform_vector(const frame3f& a, const vec3f& b) {
    return a.m * b;
}

/// transforms a direction by a frame (rigid affine transform)
template <>
constexpr inline vec3f transform_direction(const frame3f& a, const vec3f& b) {
    return a.m * b;
}

/// transforms a frame by a frame (rigid affine transform)
template <>
constexpr inline frame3f transform_frame(const frame3f& a, const frame3f& b) {
    return {a.m * b.m, a.m * b.t + a.t};
}

/// inverse transforms a point by a frame (rigid affine transform)
template <>
constexpr inline vec3f transform_point_inverse(
    const frame3f& a, const vec3f& b) {
    return (b - a.t) * a.m;
}

/// inverse transforms a vector by a frame (rigid affine transform)
template <>
constexpr inline vec3f transform_vector_inverse(
    const frame3f& a, const vec3f& b) {
    return b * a.m;
}

/// inverse transforms a direction by a frame (rigid affine transform)
template <>
constexpr inline vec3f transform_direction_inverse(
    const frame3f& a, const vec3f& b) {
    return b * a.m;
}

/// @}

#endif

/// @name specialization of transform operations
/// @{

/// transforms a ray by a matrix
template <typename T, int N>
constexpr inline ray<T, N> transform_ray(
    const mat<T, N + 1, N + 1>& a, const ray<T, N>& b) {
    return {
        transform_point(a, b.o), transform_direction(a, b.d), b.tmin, b.tmax};
}

/// transforms a bbox by a matrix
template <typename T>
constexpr inline bbox<T, 3> transform_bbox(
    const mat<T, 4, 4>& a, const bbox<T, 3>& b) {
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

/// transforms a ray by a frame (rigid affine transform)
template <typename T, int N>
constexpr inline ray<T, N> transform_ray(
    const frame<T, N>& a, const ray<T, N>& b) {
    return {
        transform_point(a, b.o), transform_direction(a, b.d), b.tmin, b.tmax};
}

/// transforms a bbox by a frame (rigid affine transform)
template <typename T>
constexpr inline bbox<T, 3> transform_bbox(
    const frame<T, 3>& a, const bbox<T, 3>& b) {
#if 0
    vec<T, 3> corners[8] = {
        {b[0][0], b[0][1], b[0][2]}, {b[0][0], b[0][1], b[1][2]},
        {b[0][0], b[1][1], b[0][2]}, {b[0][0], b[1][1], b[1][2]},
        {b[1][0], b[0][1], b[0][2]}, {b[1][0], b[0][1], b[1][2]},
        {b[1][0], b[1][1], b[0][2]}, {b[1][0], b[1][1], b[1][2]},
    };
    auto xformed = bbox<T, 3>();
    for (auto j = 0; j < 8; j++) xformed += transform_point(a, corners[j]);
    return xformed;
#else
    // Code from Real-time Collision Detection by Christer Ericson Sect. 4.2.6
    // Transform AABB a by the matrix m and translation t,
    // find maximum extents, and store result into AABB b.
    // start by adding in translation
    auto c = bbox<T, 3>{ym::pos(a), ym::pos(a)};
    // for all three axes
    for (auto i = 0; i < 3; i++) {
        // form extent by summing smaller and larger terms respectively
        for (auto j = 0; j < 3; j++) {
            auto e = ym::rot(a)[j][i] * b[0][j];
            auto f = ym::rot(a)[j][i] * b[1][j];
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
#endif
}

/// inverse transforms a ray by a frame (rigid affine transform)
template <typename T, int N>
constexpr inline ray<T, N> transform_ray_inverse(
    const frame<T, N>& a, const ray<T, N>& b) {
    return {transform_point_inverse(a, b.o),
        transform_direction_inverse(a, b.d), b.tmin, b.tmax};
}

/// inverse transforms a bbox by a frame (rigid affine transform)
template <typename T>
constexpr inline bbox<T, 3> transform_bbox_inverse(
    const frame<T, 3>& a, const bbox<T, 3>& b) {
    return transform_bbox(inverse(a), b);
}

/// @}

/// @name transform construction operations
/// @{

/// rotation matrix from axis-angle
template <typename T>
constexpr inline mat<T, 3, 3> rotation_mat3(const vec<T, 3>& axis, T angle) {
    auto s = std::sin(angle), c = std::cos(angle);
    auto vv = normalize(axis);
    return {{c + (1 - c) * vv[0] * vv[0], (1 - c) * vv[0] * vv[1] + s * vv[2],
                (1 - c) * vv[0] * vv[2] - s * vv[1]},
        {(1 - c) * vv[0] * vv[1] - s * vv[2], c + (1 - c) * vv[1] * vv[1],
            (1 - c) * vv[1] * vv[2] + s * vv[0]},
        {(1 - c) * vv[0] * vv[2] + s * vv[1],
            (1 - c) * vv[1] * vv[2] - s * vv[0], c + (1 - c) * vv[2] * vv[2]}};
}

/// translation frame
template <typename T>
constexpr inline frame<T, 3> translation_frame3(const vec<T, 3>& a) {
    return {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, a};
}

/// translation matrix
template <typename T>
constexpr inline mat<T, 4, 4> translation_mat4(const vec<T, 3>& a) {
    return (mat<T, 4, 4>)translation_frame3(a);
}

/// scaling frame (this is not rigid and it is only here for symmatry of
/// API)
template <typename T>
constexpr inline frame<T, 3> scaling_frame3(const vec<T, 3>& a) {
    return {{a[0], 0, 0}, {0, a[1], 0}, {0, 0, a[2]}, {0, 0, 0}};
}

/// scaling matrix
template <typename T>
constexpr inline mat<T, 4, 4> scaling_mat4(const vec<T, 3>& a) {
    return (mat<T, 4, 4>)scaling_frame3(a);
}

/// rotation frame
template <typename T>
constexpr inline frame<T, 3> rotation_frame3(const vec<T, 3>& axis, T angle) {
    return {rotation_mat3(axis, angle), {0, 0, 0}};
}

/// rotation matrix
template <typename T>
constexpr inline mat<T, 4, 4> rotation_mat4(const vec<T, 3>& axis, T angle) {
    return (mat<T, 4, 4>)rotation_frame3(axis, angle);
}

/// OpenGL lookat frame
template <typename T>
constexpr inline frame<T, 3> lookat_frame3(
    const vec<T, 3>& eye, const vec<T, 3>& center, const vec<T, 3>& up) {
    auto w = normalize(eye - center);
    auto u = normalize(cross(up, w));
    auto v = normalize(cross(w, u));
    return {u, v, w, eye};
}

/// OpenGL lookat matrix
template <typename T>
constexpr inline mat<T, 4, 4> lookat_mat4(
    const vec<T, 3>& eye, const vec<T, 3>& center, const vec<T, 3>& up) {
    return to_mat(lookat_frame3(eye, center, up));
}

/// OpenGL frustum matrix
template <typename T>
constexpr inline mat<T, 4, 4> frustum_mat4(T l, T r, T b, T t, T n, T f) {
    return {{2 * n / (r - l), 0, 0, 0}, {0, 2 * n / (t - b), 0, 0},
        {(r + l) / (r - l), (t + b) / (t - b), -(f + n) / (f - n), -1},
        {0, 0, -2 * f * n / (f - n), 0}};
}

/// OpenGL orthographic matrix
template <typename T>
constexpr inline mat<T, 4, 4> ortho_mat4(T l, T r, T b, T t, T n, T f) {
    return {{2 / (r - l), 0, 0, 0}, {0, 2 / (t - b), 0, 0},
        {0, 0, -2 / (f - n), 0},
        {-(r + l) / (r - l), -(t + b) / (t - b), -(f + n) / (f - n), 1}};
}

/// OpenGL orthographic 2D matrix
template <typename T>
constexpr inline mat<T, 4, 4> ortho2d_mat4(T l, T r, T b, T t) {
    return ortho_mat4(l, r, b, t, -1, 1);
}

/// OpenGL perspective matrix
template <typename T>
constexpr inline mat<T, 4, 4> perspective_mat4(
    T fovy, T aspect, T near, T far) {
    auto y = near * std::tan(fovy / 2);
    auto x = y * aspect;
    return frustum_mat4<T>(-x, x, -y, y, near, far);
}

/// @}

// -----------------------------------------------------------------------------
// GEOMETRY UTILITIES
// -----------------------------------------------------------------------------

/// @name geometry utlities
/// @{

/// triangle normal
template <typename T>
constexpr inline vec<T, 3> triangle_normal(
    const vec<T, 3>& v0, const vec<T, 3>& v1, const vec<T, 3>& v2) {
    return normalize(cross(v1 - v0, v2 - v0));
}

/// triangle area
template <typename T>
constexpr inline T triangle_area(
    const vec<T, 3>& v0, const vec<T, 3>& v1, const vec<T, 3>& v2) {
    return length(cross(v1 - v0, v2 - v0)) / 2;
}

/// tetrahedron volume
template <typename T>
constexpr inline T tetrahedron_volume(const vec<T, 3>& v0, const vec<T, 3>& v1,
    const vec<T, 3>& v2, const vec<T, 3>& v3) {
    return dot(cross(v1 - v0, v2 - v0), v3 - v0) / 6;
}

/// triangle baricentric interpolation
template <typename T, typename T1>
constexpr inline T blerp(const T& a, const T& b, const T& c, const T1& w) {
    return a * w[0] + b * w[1] + c * w[2];
}

/// @}

// -----------------------------------------------------------------------------
// UI UTILITIES
// -----------------------------------------------------------------------------

/// @name ui utilities
/// @{

/// Turntable for UI navigation from a from/to/up parametrization of the camera.
template <typename T>
constexpr inline void turntable(vec<T, 3>& from, vec<T, 3>& to, vec<T, 3>& up,
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
            pan[0] * x[1] + pan[1] * y[1], pan[0] * x[2] + pan[1] * y[2]};
        *from += t;
        *to += t;
    }
}

/// Turntable for UI navigation for a frame/distance parametrization of the
/// camera.
template <typename T>
constexpr inline void turntable(frame<T, 3>& frame, float& focus,
    const vec<T, 2>& rotate, T dolly, const vec<T, 2>& pan) {
    // rotate if necessary
    if (rotate[0] || rotate[1]) {
        auto phi = std::atan2(frame[2][2], frame[2][0]) + rotate[0];
        auto theta = std::acos(frame[2][1]) + rotate[1];
        theta = max(T(0.001), min(theta, T(pi - 0.001)));
        auto new_z = vec<T, 3>{std::sin(theta) * std::cos(phi), std::cos(theta),
            std::sin(theta) * std::sin(phi)};
        auto new_center = pos(frame) - frame[2] * focus;
        auto new_o = new_center + new_z * focus;
        frame = lookat_frame3(new_o, new_center, {0, 1, 0});
        focus = dist(new_o, new_center);
    }

    // pan if necessary
    if (dolly) {
        auto c = pos(frame) - frame[2] * focus;
        focus = max(focus + dolly, T(0.001));
        pos(frame) = c + frame[2] * focus;
    }

    // pan if necessary
    if (pan[0] || pan[1]) {
        pos(frame) += frame[0] * pan[0] + frame[1] * pan[1];
    }
}

/// @}

// -----------------------------------------------------------------------------
// RANDOM NUMBER GENERATION
// -----------------------------------------------------------------------------

/// @name random numbers
/// @{

///
/// PCG random numbers. A family of random number generators that supports
/// multiple sequences. In our code, we allocate one sequence for each sample.
/// PCG32 from http://www.pcg-random.org/
///
struct rng_pcg32 {
    uint64_t state, inc;
};

/// Next random number
constexpr inline uint32_t next(rng_pcg32* rng) {
    uint64_t oldstate = rng->state;
    rng->state = oldstate * 6364136223846793005ull + (rng->inc | 1u);
    uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
    uint32_t rot = oldstate >> 59u;
    return (xorshifted >> rot) | (xorshifted << ((-((int32_t)rot)) & 31));
}

/// Init a random number generator with a state state from the sequence seq.
constexpr inline void init(rng_pcg32* rng, uint64_t state, uint64_t seq) {
    rng->state = 0U;
    rng->inc = (seq << 1u) | 1u;
    next(rng);
    rng->state += state;
    next(rng);
}

/// Next random float in [0,1).
inline float next1f(rng_pcg32* rng) { return (float)ldexp(next(rng), -32); }

/// Next random float in [0,1)x[0,1).
inline vec2f next2f(rng_pcg32* rng) { return {next1f(rng), next1f(rng)}; }

/// @}

// -----------------------------------------------------------------------------
// HASHING
// -----------------------------------------------------------------------------

/// @name hashing
/// @{

/// Computes the i-th term of a permutation of l values keyed by p.
/// From Correlated Multi-Jittered Sampling by Kensler @ Pixar
constexpr inline uint32_t hash_permute(uint32_t i, uint32_t n, uint32_t key) {
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

/// Computes a float value by hashing i with a key p.
/// From Correlated Multi-Jittered Sampling by Kensler @ Pixar
constexpr inline float hash_randfloat(uint32_t i, uint32_t key) {
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

/// 64 bit integer hash. Public domain code.
constexpr inline uint64_t hash_uint64(uint64_t a) {
    a = (~a) + (a << 21);  // a = (a << 21) - a - 1;
    a ^= (a >> 24);
    a += (a << 3) + (a << 8);  // a * 265
    a ^= (a >> 14);
    a += (a << 2) + (a << 4);  // a * 21
    a ^= (a >> 28);
    a += (a << 31);
    return a;
}

/// 64-to-32 bit integer hash. Public domain code.
constexpr inline uint32_t hash_uint64_32(uint64_t a) {
    a = (~a) + (a << 18);  // a = (a << 18) - a - 1;
    a ^= (a >> 31);
    a *= 21;  // a = (a + (a << 2)) + (a << 4);
    a ^= (a >> 11);
    a += (a << 6);
    a ^= (a >> 22);
    return (uint32_t)a;
}

/// Combines two 64 bit hashes as in boost::hash_combine
constexpr inline int hash_combine(int a, int b) {
    return a ^ (b + 0x9e3779b9 + (a << 6) + (a >> 2));
}

/// Hash a vector with hash_combine() and std::hash
template <typename T, int N>
constexpr inline int hash_vec(const vec<T, N>& v) {
    std::hash<T> Th;
    int h = 0;
    for (auto i = 0; i < N; i++) {
        h ^= (Th(v[i]) + 0x9e3779b9 + (h << 6) + (h >> 2));
    }
    return h;
}

/// @}

// -----------------------------------------------------------------------------
// VIEW CONTAINERS
// -----------------------------------------------------------------------------

/// @name view containers
/// @{

///
/// An array_view is a non-owining reference to an array with an API
/// similar
/// to a vector/array containers, but without reallocation.
/// This is inspired, but significantly simpler than
/// gsl::span https://github.com/Microsoft/GSL or array_view.
///
template <typename T>
struct array_view {
    // constructors
    constexpr array_view() noexcept : _num(0), _data(nullptr) {}
    constexpr array_view(int num, T* data) noexcept
        : _num((data) ? num : 0), _data(data) {}
    template <typename T1>
    explicit constexpr array_view(int num, T1* data) noexcept
        : _num((data) ? num : 0), _data(data) {}
    constexpr array_view(T* begin, T* end) noexcept
        : _num((end - begin)), _data(begin) {}
    constexpr array_view(const array_view<std::add_const_t<T>>& av)
        : _num(av.size()), _data(av.data()) {}
    template <typename T1>
    constexpr array_view(std::vector<T1>& av)
        : _num(av.size()), _data(av.data()) {}
    template <typename T1>
    constexpr array_view(const std::vector<T1>& av)
        : _num(av.size()), _data(av.data()) {}

    // size
    constexpr int size() const noexcept { return _num; }
    constexpr bool empty() const noexcept { return _num == 0; }
    constexpr operator bool() const noexcept { return !empty(); }

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

///
/// An image_view is a non-owining reference to an array that allows
/// to access it as a row-major image. This is inspired, but significantly
/// simpler than gsl::multi_span https://github.com/Microsoft/GSL or
/// array_view.
///
template <typename T>
struct image_view {
    // constructors
    constexpr image_view() noexcept : _size{0, 0}, _data(nullptr) {}
    constexpr image_view(const vec2i& size, T* data) noexcept
        : _size((data) ? size : vec2i{0, 0}), _data(data) {}

    // size
    constexpr vec2i size() const noexcept { return _size; }
    constexpr bool empty() const noexcept {
        return _size[0] == 0 || _size[1] == 0;
    }
    constexpr operator bool() const noexcept { return !empty(); }

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

/// @}

// -----------------------------------------------------------------------------
// IMAGE OPERATIONS
// -----------------------------------------------------------------------------

/// @name image operations
/// @{

/// Lookup an image value from a generic image
template <typename T>
constexpr inline vec<T, 4> image_lookup(
    int width, int height, int ncomp, const T* img, int x, int y, T alpha = 0) {
    auto v = img + (y * width + x) * ncomp;
    switch (ncomp) {
        case 1: return {v[0], 0, 0, alpha};
        case 2: return {v[0], v[1], 0, alpha};
        case 3: return {v[0], v[1], v[2], alpha};
        case 4: return {v[0], v[1], v[2], v[3]};
        default: assert(false); return {0, 0, 0, 0};
    }
}

/// Set an image value for a generic image
template <typename T>
constexpr inline void image_set(int width, int height, int ncomp, T* img, int x,
    int y, const vec<T, 4>& vv) {
    auto v = img + (y * width + x) * ncomp;
    switch (ncomp) {
        case 1: v[0] = vv[0]; break;
        case 2:
            v[0] = vv[0];
            v[1] = vv[1];
            break;
        case 3:
            v[0] = vv[0];
            v[1] = vv[1];
            v[2] = vv[2];
            break;
        case 4:
            v[0] = vv[0];
            v[1] = vv[1];
            v[2] = vv[2];
            v[3] = vv[3];
            break;
        default: assert(false);
    }
}

/// Conversion from srgb.
inline vec4f srgb_to_linear(const vec4b& srgb) {
    return {std::pow((float)srgb[0] / 255.0f, 2.2f),
        std::pow((float)srgb[1] / 255.0f, 2.2f),
        std::pow((float)srgb[2] / 255.0f, 2.2f), (float)srgb[3] / 255.0f};
}

/// Conversion to srgb.
inline vec4b linear_to_srgb(const vec4f& srgb) {
    auto v = vec4f{std::pow(srgb[0], 1 / 2.2f), std::pow(srgb[1], 1 / 2.2f),
        std::pow(srgb[2], 1 / 2.2f), srgb[3]};
    return {(unsigned char)(clamp(v[0], 0.0f, 1.0f) * 255),
        (unsigned char)(clamp(v[1], 0.0f, 1.0f) * 255),
        (unsigned char)(clamp(v[2], 0.0f, 1.0f) * 255),
        (unsigned char)(clamp(v[3], 0.0f, 1.0f) * 255)};
}

/// Conversion from clamped bytes.
inline vec4f byte_to_linear(const vec4b& v) {
    return {(float)v[0] / 255.0f, (float)v[1] / 255.0f, (float)v[2] / 255.0f,
        (float)v[3] / 255.0f};
}

/// Conversion to clamped bytes.
inline vec4b linear_to_byte(const vec4f& v) {
    return {(unsigned char)(clamp(v[0], 0.0f, 1.0f) * 255),
        (unsigned char)(clamp(v[1], 0.0f, 1.0f) * 255),
        (unsigned char)(clamp(v[2], 0.0f, 1.0f) * 255),
        (unsigned char)(clamp(v[3], 0.0f, 1.0f) * 255)};
}

/// Exposure/gamma correction
inline void exposure_gamma(int width, int height, int ncomp, const float* hdr,
    float* ldr, float exposure, float gamma, bool clamped) {
    auto s = std::pow(2.0f, exposure);
    for (auto j = 0; j < height; j++) {
        for (auto i = 0; i < width; i++) {
            auto v = image_lookup(width, height, ncomp, hdr, i, j);
            v = {std::pow(s * v[0], 1 / gamma), std::pow(s * v[1], 1 / gamma),
                std::pow(s * v[2], 1 / gamma), v[3]};
            if (clamped) v = clamp(v, 0.0f, 1.0f);
            image_set(width, height, ncomp, ldr, i, j, v);
        }
    }
}

/// Exposure/gamma correction
inline void exposure_gamma(int width, int height, int ncomp, const float* hdr,
    byte* ldr, float exposure, float gamma, bool srgb_output) {
    auto s = std::pow(2.0f, exposure);
    for (auto j = 0; j < height; j++) {
        for (auto i = 0; i < width; i++) {
            auto v = image_lookup(width, height, ncomp, hdr, i, j);
            v = {std::pow(s * v[0], 1 / gamma), std::pow(s * v[1], 1 / gamma),
                std::pow(s * v[2], 1 / gamma), v[3]};
            if (srgb_output)
                image_set(width, height, ncomp, ldr, i, j, linear_to_srgb(v));
            else
                image_set(width, height, ncomp, ldr, i, j, linear_to_byte(v));
        }
    }
}

/// Exposure/gamma correction
inline void exposure_gamma(int width, int height, int ncomp, const float* hdr,
    byte* ldr, float exposure, float gamma, bool srgb_output, int x, int y,
    int w, int h) {
    auto s = std::pow(2.0f, exposure);
    for (auto j = y; j < y + h; j++) {
        for (auto i = x; i < x + w; i++) {
            auto v = image_lookup(width, height, ncomp, hdr, i, j);
            v = {std::pow(s * v[0], 1 / gamma), std::pow(s * v[1], 1 / gamma),
                std::pow(s * v[2], 1 / gamma), v[3]};
            if (srgb_output)
                image_set(width, height, ncomp, ldr, i, j, linear_to_srgb(v));
            else
                image_set(width, height, ncomp, ldr, i, j, linear_to_byte(v));
        }
    }
}

/// linear to srgb correction
inline void linear_to_srgb(
    int width, int height, int ncomp, const float* hdr, byte* ldr) {
    for (auto j = 0; j < height; j++) {
        for (auto i = 0; i < width; i++) {
            auto v = image_lookup(width, height, ncomp, hdr, i, j);
            image_set(width, height, ncomp, ldr, i, j, linear_to_srgb(v));
        }
    }
}

/// linear to byte conversion
inline void linear_to_byte(
    int width, int height, int ncomp, const float* hdr, byte* ldr) {
    for (auto j = 0; j < height; j++) {
        for (auto i = 0; i < width; i++) {
            auto v = image_lookup(width, height, ncomp, hdr, i, j);
            image_set(width, height, ncomp, ldr, i, j, linear_to_byte(v));
        }
    }
}

/// @}

// -----------------------------------------------------------------------------
// TIMER
// -----------------------------------------------------------------------------

/// @name timer
/// @{

/// A simple wrapper for std::chrono.
struct timer {
    /// initialize a timer and start it if necessary
    timer(bool autostart = true) {
        if (autostart) start();
    }

    /// start a timer
    void start() {
        _start = std::chrono::steady_clock::now();
        _started = true;
    }

    /// stops a timer
    void stop() {
        _end = std::chrono::steady_clock::now();
        _started = false;
    }

    /// elapsed time
    double elapsed() {
        if (_started) stop();
        std::chrono::duration<double> diff = (_end - _start);
        return diff.count();
    }

   private:
    bool _started = false;
    std::chrono::time_point<std::chrono::steady_clock> _start, _end;
};

/// @}

}  // namespace

// HACK to avoid compilation with MSVC2015 without dirtying code
#ifdef _WIN32
#undef constexpr
#endif

#endif
