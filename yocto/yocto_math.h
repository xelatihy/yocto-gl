///
/// # Yocto/Math
///
/// A collection of vector math functions and simple containers
/// used to implement YOCTO. Features include
///
/// - a few convenience math functions
/// - static length float vectors, with specialization for 2, 3, 4 length
/// - static length matrices, with specialization for 2x2, 3x3, 4x4
/// - static length rigid transforms (frames), specialized for 2d and 3d space
/// - linear algebra operations and transforms for fixed length matrices/vecs
/// - axis aligned bounding boxes
/// - rays
/// - ray-primitive intersection
/// - point-primitive distance and overlap tests
/// - normal amd tangent computation for meshes and lines
/// - generation of tesselated meshes
/// - random number generation via PCG32
/// - a few hash functions
/// - trivial image data structue and a few image operations
///
/// We developed our own library since we felt that all existing ones are either
/// complete, but unreadable or with lots of dependencies, or just as incomplete
/// and untested as ours.
///
/// This library has no dependencies.
/// Some templated types and functions use specialization for easier access
/// and faster compilation. Specialization can be disabled by defining
/// YM_NO_SPECIALIZATION. Specialization is only supported fully on clang.
///
/// This library includes code from the PCG random number generator,
/// boost hash_combine, Pixar multijittered sampling, code from "Real-Time
/// Collision Detection" by Christer Ericson and public domain code from
/// - https://github.com/sgorsten/linalg
/// - https://gist.github.com/badboy/6267743
///
///
/// ## History
///
/// - v 0.15: enable specialization always
/// - v 0.14: move timer to Yocto/Utils
/// - v 0.13: more shape functions
/// - v 0.12: documentation update
/// - v 0.11: added more matrix and quaternion operations
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
#include <cassert>
#include <cmath>
#include <cstdint>
#include <functional>
#include <initializer_list>
#include <limits>
#include <vector>

// HACK to avoid compilation with MSVC2015 and C++11 without dirtying code
#if defined(_WIN32) || __cplusplus < 201402L
#define constexpr
#endif

///
/// Math types and utlities for 3D graphics and imaging
///
namespace ym {

// -----------------------------------------------------------------------------
// BASIC TYPEDEFS
// -----------------------------------------------------------------------------

/// convenient typedef for bytes
using byte = unsigned char;

/// convenient typedef for bytes
using uint = unsigned int;

// -----------------------------------------------------------------------------
// MATH CONSTANTS
// -----------------------------------------------------------------------------

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

// -----------------------------------------------------------------------------
// BASIC MATH FUNCTIONS
// -----------------------------------------------------------------------------

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

/// Safe float to byte conversion
constexpr inline byte float_to_byte(float x) {
    return max(0, min(int(x * 256), 255));
}

/// Safe byte to float conversion
constexpr inline float byte_to_float(byte x) { return (float)x / 255.0f; }

// -----------------------------------------------------------------------------
// FUNCTIONS BROUGHT INTO NAMESPACE
// -----------------------------------------------------------------------------

/// sqrt
using std::sqrt;
/// pow
using std::pow;
/// sin
using std::sin;
/// cos
using std::cos;
/// tan
using std::tan;
/// asin
using std::asin;
/// acos
using std::acos;
/// atan2
using std::atan2;
/// abs
using std::abs;

// -----------------------------------------------------------------------------
// VECTORS
// -----------------------------------------------------------------------------

///
/// Vector of elements of compile time dimension with default initializer.
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

    /// element access
    constexpr T& operator[](int i) { return v[i]; }
    /// element access
    constexpr const T& operator[](int i) const { return v[i]; }

    /// data access
    constexpr T* data() { return v; }
    /// data access
    constexpr const T* data() const { return v; }

    /// element data
    T v[N];
};

#ifndef YM_NO_SPECIALIZATION

///
/// Specialization of vectors for 1 component and float coordinates.
///
template <typename T>
struct vec<T, 1> {
    /// size
    constexpr static const int N = 1;

    /// default constructor
    constexpr vec() : x{0} {}
    /// element constructor
    constexpr vec(T x) : x{x} {}

    /// element access
    constexpr T& operator[](int i) { return (&x)[i]; }
    /// element access
    constexpr const T& operator[](int i) const { return (&x)[i]; }

    /// data access
    constexpr T* data() { return &x; }
    /// data access
    constexpr const T* data() const { return &x; }

    /// element data
    T x;
};

///
/// Specialization of vectors for 2 components and float coordinates.
///
template <typename T>
struct vec<T, 2> {
    /// size
    constexpr static const int N = 2;

    /// default constructor
    constexpr vec() : x{0}, y{0} {}
    /// element constructor
    constexpr vec(T x, T y) : x{x}, y{y} {}

    /// element access
    constexpr T& operator[](int i) { return (&x)[i]; }
    /// element access
    constexpr const T& operator[](int i) const { return (&x)[i]; }

    /// data access
    constexpr T* data() { return &x; }
    /// data access
    constexpr const T* data() const { return &x; }

    /// element data
    T x;
    /// element data
    T y;
};

///
/// Specialization of vectors for 3 components and float coordinates.
///
template <typename T>
struct vec<T, 3> {
    /// size
    constexpr static const int N = 3;

    /// default constructor
    constexpr vec() : x{0}, y{0}, z{0} {}
    /// element constructor
    constexpr vec(T x, T y, T z) : x{x}, y{y}, z{z} {}

    /// element access
    constexpr T& operator[](int i) { return (&x)[i]; }
    /// element access
    constexpr const T& operator[](int i) const { return (&x)[i]; }

    /// data access
    constexpr T* data() { return &x; }
    /// data access
    constexpr const T* data() const { return &x; }

    /// element data
    T x;
    /// element data
    T y;
    /// element data
    T z;
};

///
/// Specialization of vectors for 4 components and float coordinates.
///
template <typename T>
struct vec<T, 4> {
    /// size
    constexpr static const int N = 4;

    /// default constructor
    constexpr vec() : x{0}, y{0}, z{0}, w{0} {}
    /// element constructor
    constexpr vec(T x, T y, T z, T w) : x{x}, y{y}, z{z}, w{w} {}

    /// element access
    constexpr T& operator[](int i) { return (&x)[i]; }
    /// element access
    constexpr const T& operator[](int i) const { return (&x)[i]; }

    /// data access
    constexpr T* data() { return &x; }
    /// data access
    constexpr const T* data() const { return &x; }

    /// element data
    T x;
    /// element data
    T y;
    /// element data
    T z;
    /// element data
    T w;
};

#endif

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

/// iteration support
template <typename T, int N>
constexpr inline T* begin(vec<T, N>& a) {
    return &a[0];
}

/// iteration support
template <typename T, int N>
constexpr inline const T* begin(const vec<T, N>& a) {
    return &a[0];
}

/// iteration support
template <typename T, int N>
constexpr inline T* end(vec<T, N>& a) {
    return &a[0] + N;
}

/// iteration support
template <typename T, int N>
constexpr inline const T* end(const vec<T, N>& a) {
    return &a[0] + N;
}

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

/// vector operator < (lexicographic order - useful for std::map)
template <typename T, int N>
constexpr inline bool operator<(const vec<T, N>& a, const vec<T, N>& b) {
    for (auto i = 0; i < N; i++) {
        if (a[i] < b[i]) return true;
        if (a[i] > b[i]) return false;
    }
    return false;
}

#ifndef YM_NO_SPECIALIZATION

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

#endif

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

/// vector operator +
template <typename T, int N>
constexpr inline vec<T, N> operator+(const vec<T, N>& a, const T b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a[i] + b;
    return c;
}

/// vector operator -
template <typename T, int N>
constexpr inline vec<T, N> operator-(const vec<T, N>& a, const T b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a[i] - b;
    return c;
}

/// vector operator +
template <typename T, int N>
constexpr inline vec<T, N> operator+(T a, const vec<T, N>& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a + b[i];
    return c;
}

/// vector operator -
template <typename T, int N>
constexpr inline vec<T, N> operator-(T a, const vec<T, N>& b) {
    vec<T, N> c;
    for (auto i = 0; i < N; i++) c[i] = a - b[i];
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

#ifndef YM_NO_SPECIALIZATION

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

#endif

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

#ifndef YM_NO_SPECIALIZATION

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

#endif

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
    return th == 0 ?
               a :
               a * (sin(th * (1 - t)) / sin(th)) + b * (sin(th * t) / sin(th));
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

/// Element-wise pow
template <typename T, int N>
constexpr inline vec<T, N> pow(const vec<T, N>& a, const T b) {
    auto c = vec<T, N>();
    for (auto i = 0; i < N; i++) c[i] = pow(a[i], b);
    return c;
}

/// Element-wise conversion
template <int N>
constexpr inline vec<byte, N> float_to_byte(const vec<float, N>& a) {
    auto c = vec<byte, N>();
    for (auto i = 0; i < N; i++) c[i] = float_to_byte(a[i]);
    return c;
}

/// Element-wise conversion
template <int N>
constexpr inline vec<float, N> byte_to_float(const vec<byte, N>& a) {
    auto c = vec<float, N>();
    for (auto i = 0; i < N; i++) c[i] = byte_to_float(a[i]);
    return c;
}

// -----------------------------------------------------------------------------
// MATRICES
// -----------------------------------------------------------------------------

///
/// Matrix of elements of compile time dimensions, stored in column major
/// format, with default initializer.
/// Colums access via operator[].
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

    /// element access
    constexpr V& operator[](int i) { return v[i]; }
    /// element access
    constexpr const V& operator[](int i) const { return v[i]; }

    /// data access
    constexpr V* data() { return v; }
    /// data access
    constexpr const V* data() const { return v; }

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

    /// element access
    constexpr V& operator[](int i) { return (&x)[i]; }
    /// element access
    constexpr const V& operator[](int i) const { return (&x)[i]; }

    /// data access
    constexpr V* data() { return &x; }
    /// data access
    constexpr const V* data() const { return &x; }

    /// element data
    V x;
    /// element data
    V y;
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

    /// element access
    constexpr V& operator[](int i) { return (&x)[i]; }
    /// element access
    constexpr const V& operator[](int i) const { return (&x)[i]; }

    /// data access
    constexpr V* data() { return &x; }
    /// data access
    constexpr const V* data() const { return &x; }

    /// element data
    V x;
    /// element data
    V y;
    /// element data
    V z;
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

    /// element access
    constexpr V& operator[](int i) { return (&x)[i]; }
    /// element access
    constexpr const V& operator[](int i) const { return (&x)[i]; }

    /// data access
    constexpr V* data() { return &x; }
    /// data access
    constexpr const V* data() const { return &x; }

    /// element data
    V x;
    /// element data
    V y;
    /// element data
    V z;
    /// element data
    V w;
};

#endif

/// 1-dimensional float matrix
using mat1f = mat<float, 1, 1>;
/// 2-dimensional float matrix
using mat2f = mat<float, 2, 2>;
/// 3-dimensional float matrix
using mat3f = mat<float, 3, 3>;
/// 4-dimensional float matrix
using mat4f = mat<float, 4, 4>;

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

/// 1-dimensional float identity matrix
const auto identity_mat1f = identity_mat<float, 1>();
/// 2-dimensional float identity matrix
const auto identity_mat2f = identity_mat<float, 2>();
/// 3-dimensional float identity matrix
const auto identity_mat3f = identity_mat<float, 3>();
/// 4-dimensional float identity matrix
const auto identity_mat4f = identity_mat<float, 4>();

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

#ifndef YM_NO_SPECIALIZATION

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

#endif

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

// -----------------------------------------------------------------------------
// RIGID BODY TRANSFORMS/FRAMES
// -----------------------------------------------------------------------------

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
    constexpr frame(const M& m, const V& t) : x(m.x), y(m.y), o(t) {}

    /// element access
    constexpr V& operator[](int i) { return (&x)[i]; }
    /// element access
    constexpr const V& operator[](int i) const { return (&x)[i]; }

    /// access position
    constexpr V& pos() { return o; }
    /// access position
    constexpr const V& pos() const { return o; }

    /// access rotation
    constexpr M& rot() { return *(M*)(&x); }
    /// access rotation
    constexpr const M& rot() const { return *(M*)(&x); }

    /// element data
    V x;
    /// element data
    V y;
    /// element data
    V o;
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
    constexpr frame(const M& m, const V& t) : x(m.x), y(m.y), z(m.z), o(t) {}

    /// element access
    constexpr V& operator[](int i) { return (&x)[i]; }
    /// element access
    constexpr const V& operator[](int i) const { return (&x)[i]; }

    /// access position
    constexpr V& pos() { return o; }
    /// access position
    constexpr const V& pos() const { return o; }

    /// access rotation
    constexpr M& rot() { return *(M*)(&x); }
    /// access rotation
    constexpr const M& rot() const { return *(M*)(&x); }

    /// element data
    V x;
    /// element data
    V y;
    /// element data
    V z;
    /// element data
    V o;
};

#endif

/// 1-dimensional float frame
using frame1f = frame<float, 1>;
/// 2-dimensional float frame
using frame2f = frame<float, 2>;
/// 3-dimensional float frame
using frame3f = frame<float, 3>;
/// 4-dimensional float frame
using frame4f = frame<float, 4>;

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

#ifndef YM_NO_SPECIALIZATION

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

#endif

/// 1-dimensional float identity frame
const auto identity_frame1f = identity_frame<float, 1>();
/// 2-dimensional float identity frame
const auto identity_frame2f = identity_frame<float, 2>();
/// 3-dimensional float identity frame
const auto identity_frame3f = identity_frame<float, 3>();
/// 4-dimensional float identity frame
const auto identity_frame4f = identity_frame<float, 4>();

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

/// frame composition (equivalent to affine matrix multiply)
template <typename T, int N>
constexpr inline frame<T, N> operator*(
    const frame<T, N>& a, const frame<T, N>& b) {
    return {rot(a) * rot(b), rot(a) * pos(b) + pos(a)};
}

/// frame inverse (equivalent to rigid affine inverse)
template <typename T, int N>
constexpr inline frame<T, N> inverse(const frame<T, N>& a) {
    auto minv = transpose(rot(a));
    return {minv, -(minv * pos(a))};
}

#ifndef YM_NO_SPECIALIZATION

/// frame composition (equivalent to affine matrix multiply)
template <>
constexpr inline frame3f operator*(const frame3f& a, const frame3f& b) {
    return {a.rot() * b.rot(), a.rot() * b.pos() + a.pos()};
}

/// frame inverse (equivalent to rigid affine inverse)
template <>
constexpr inline frame3f inverse(const frame3f& a) {
    auto minv = transpose(a.rot());
    return {minv, -(minv * a.pos())};
}

#endif

// -----------------------------------------------------------------------------
// QUATERNIONS
// -----------------------------------------------------------------------------

///
/// Quaternion placeholder. Only helpful in the specialization.
///
template <typename T, int N>
struct quat;

///
/// Quaternions implemented as a vec<T,4>. Data access via operator[].
/// Quaterions are xi + yj + zk + w.
///
template <typename T>
struct quat<T, 4> {
    /// size
    constexpr static const int N = 4;

    /// default constructor
    constexpr quat() : x{0}, y{0}, z{0}, w{1} {}

    // list constructor
    constexpr quat(const T& x, const T& y, const T& z, const T& w)
        : x{x}, y{y}, z{z}, w{w} {}

    /// conversion from vec
    constexpr explicit quat(const vec<T, N>& vv)
        : x{vv.x}, y{vv.y}, z{vv.z}, w{vv.w} {}
    /// conversion to vec
    constexpr explicit operator vec<T, N>() const { return {x, y, z, w}; }

    /// element access
    constexpr T& operator[](int i) { return (&x)[i]; }
    /// element access
    constexpr const T& operator[](int i) const { return (&x)[i]; }

    /// data
    T x;
    /// data
    T y;
    /// data
    T z;
    /// data
    T w;
};

/// float quaterion
using quat4f = quat<float, 4>;

/// float identity quaterion
const auto identity_quat4f = quat<float, 4>{0, 0, 0, 1};

/// vector operator ==
template <typename T, int N>
constexpr inline bool operator==(const quat<T, N>& a, const quat<T, N>& b) {
    for (auto i = 0; i < N; i++)
        if (a[i] != b[i]) return false;
    return true;
}

/// vector operator !=
template <typename T, int N>
constexpr inline bool operator!=(const quat<T, N>& a, const quat<T, N>& b) {
    return !(a == b);
}

/// quaterion multiply
template <typename T>
constexpr quat<T, 4> operator*(const quat<T, 4>& a, const quat<T, 4>& b) {
    return {a[0] * b[3] + a[3] * b[0] + a[1] * b[3] - a[2] * b[1],
        a[1] * b[3] + a[3] * b[1] + a[2] * b[0] - a[0] * b[2],
        a[2] * b[3] + a[3] * b[2] + a[0] * b[1] - a[1] * b[0],
        a[3] * b[3] - a[0] * b[0] - a[1] * b[1] - a[2] * b[2]};
}

/// quaterion conjugate
template <typename T>
constexpr quat<T, 4> conjugate(const quat<T, 4>& v) {
    return {-v[0], -v[1], -v[2], v[3]};
}

/// quaterion inverse
template <typename T>
constexpr quat<T, 4> inverse(const quat<T, 4>& v) {
    return qconj(v) / lengthsqr(vec<T, 4>(v));
}

/// quaterion inverse
template <typename T>
constexpr quat<T, 4> normalize(const quat<T, 4>& v) {
    auto l = length(vec<T, 4>{v.x, v.y, v.z, v.w});
    if (!l) return {0, 0, 0, 1};
    return {v.x / l, v.y / l, v.z / l, v.w / l};
}

/// quaterion normalized linear interpolation
template <typename T>
constexpr quat<T, 4> nlerp(const quat<T, 4>& a, const quat<T, 4>& b, T t) {
    return nlerp(vec<T, 4>(a),
        dot(vec<T, 4>(a), vec<T, 4>(b)) < 0 ? -vec<T, 4>(b) : vec<T, 4>(b), t);
}

/// quaterion spherical linear interpolation
template <typename T>
constexpr quat<T, 4> slerp(const quat<T, 4>& a, const quat<T, 4>& b, T t) {
    auto a_ = vec<float, 4>{a.x, a.y, a.z, a.w};
    auto b_ = vec<float, 4>{b.x, b.y, b.z, b.w};
    return slerp(a_, dot(a_, b_) < 0 ? -b_ : b_, t);
}

// -----------------------------------------------------------------------------
// AXIS ALIGNED BOUNDING BOXES
// -----------------------------------------------------------------------------

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
            min[i] = std::numeric_limits<T>::max();
            max[i] = std::numeric_limits<T>::lowest();
        }
    }

    /// list constructor
    constexpr bbox(const vec<T, N>& m, const vec<T, N>& M) : min{m}, max{M} {}

    /// element access
    constexpr V& operator[](int i) { return (&min)[i]; }
    /// element access
    constexpr const V& operator[](int i) const { return (&min)[i]; }

    /// element data
    V min;
    /// element data
    V max;
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
    constexpr V& operator[](int i) { return (&min)[i]; }
    /// element access
    constexpr const V& operator[](int i) const { return (&min)[i]; }

    /// element data
    V min;
    /// element data
    V max;
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
    constexpr V& operator[](int i) { return (&min)[i]; }
    /// element access
    constexpr const V& operator[](int i) const { return (&min)[i]; }

    /// element data
    V min;
    /// element data
    V max;
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
    constexpr V& operator[](int i) { return (&min)[i]; }
    /// element access
    constexpr const V& operator[](int i) const { return (&min)[i]; }

    /// element data
    V min;
    /// element data
    V max;
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
    constexpr V& operator[](int i) { return (&min)[i]; }
    /// element access
    constexpr const V& operator[](int i) const { return (&min)[i]; }

    /// element data
    V min;
    /// element data
    V max;
};

#endif

/// 1-dimensional float bbox
using bbox1f = bbox<float, 1>;
/// 2-dimensional float bbox
using bbox2f = bbox<float, 2>;
/// 3-dimensional float bbox
using bbox3f = bbox<float, 3>;
/// 4-dimensional float bbox
using bbox4f = bbox<float, 4>;

/// initializes an empty bbox
template <typename T, int N>
constexpr inline bbox<T, N> invalid_bbox() {
    auto a = bbox<T, N>();
    for (auto i = 0; i < N; i++) {
        a.min[i] = std::numeric_limits<T>::max();
        a.max[i] = std::numeric_limits<T>::lowest();
    }
    return a;
}

/// initialize a bonding box from a list of points
template <typename T, int N>
constexpr inline bbox<T, N> make_bbox(int count, const vec<T, N>* v) {
    auto a = invalid_bbox<T, N>();
    for (auto j = 0; j < count; j++) {
        auto&& vv = v[j];
        for (auto i = 0; i < N; i++) {
            a.min[i] = min(a.min[i], vv[i]);
            a.max[i] = max(a.max[i], vv[i]);
        }
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
            a.min[i] = min(a.min[i], vv[i]);
            a.max[i] = max(a.max[i], vv[i]);
        }
    }
    return a;
}

/// 1-dimensional float empty bbox
const auto invalid_bbox1f = bbox1f();
/// 2-dimensional float empty bbox
const auto invalid_bbox2f = bbox2f();
/// 3-dimensional float empty bbox
const auto invalid_bbox3f = bbox3f();
/// 4-dimensional float empty bbox
const auto invalid_bbox4f = bbox4f();

/// computes the center of a bbox
template <typename T, int N>
constexpr inline vec<T, N> center(const bbox<T, N>& a) {
    return (a.min + a.max) / (T)2;
}

/// computes the diagonal of a bbox
template <typename T, int N>
constexpr inline vec<T, N> diagonal(const bbox<T, N>& a) {
    return a.max - a.min;
}

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

#ifndef YM_NO_SPECIALIZATION

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

#endif

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

// -----------------------------------------------------------------------------
// RAYS
// -----------------------------------------------------------------------------

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

/// 1-dimensional float ray
using ray1f = ray<float, 1>;
/// 2-dimensional float ray
using ray2f = ray<float, 2>;
/// 3-dimensional float ray
using ray3f = ray<float, 3>;
/// 4-dimensional float ray
using ray4f = ray<float, 4>;

/// evalutes the position along the ray
template <typename T, int N>
constexpr inline vec<T, N> eval(const ray<T, N>& ray, T t) {
    return ray.o + t * ray.d;
}

// -----------------------------------------------------------------------------
// TRANSFORMS
// -----------------------------------------------------------------------------

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
    return rot(a) * b + pos(a);
}

/// transforms a vector by a frame (rigid affine transform)
template <typename T, int N>
constexpr inline vec<T, N> transform_vector(
    const frame<T, N>& a, const vec<T, N>& b) {
    return rot(a) * b;
}

/// transforms a direction by a frame (rigid affine transform)
template <typename T, int N>
constexpr inline vec<T, N> transform_direction(
    const frame<T, N>& a, const vec<T, N>& b) {
    return rot(a) * b;
}

/// transforms a frame by a frame (rigid affine transform)
template <typename T, int N>
constexpr inline frame<T, N> transform_frame(
    const frame<T, N>& a, const frame<T, N>& b) {
    return {rot(a) * rot(b), pos(a) * pos(b) + pos(a)};
}

/// inverse transforms a point by a frame (rigid affine transform)
template <typename T, int N>
constexpr inline vec<T, N> transform_point_inverse(
    const frame<T, N>& a, const vec<T, N>& b) {
    return (b - pos(a)) * rot(a);
}

/// inverse transforms a vector by a frame (rigid affine transform)
template <typename T, int N>
constexpr inline vec<T, N> transform_vector_inverse(
    const frame<T, N>& a, const vec<T, N>& b) {
    return b * rot(a);
}

/// inverse transforms a direction by a frame (rigid affine transform)
template <typename T, int N>
constexpr inline vec<T, N> transform_direction_inverse(
    const frame<T, N>& a, const vec<T, N>& b) {
    return b * rot(a);
}

#ifndef YM_NO_SPECIALIZATION

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
    return vec3f{tvb.x, tvb.y, tvb.z};
}

/// transforms a direction by a matrix
template <>
constexpr inline vec3f transform_direction(const mat4f& a, const vec3f& b) {
    return normalize(transform_vector(a, b));
}

/// transforms a point by a frame (rigid affine transform)
template <>
constexpr inline vec3f transform_point(const frame3f& a, const vec3f& b) {
    return a.rot() * b + a.pos();
}

/// transforms a vector by a frame (rigid affine transform)
template <>
constexpr inline vec3f transform_vector(const frame3f& a, const vec3f& b) {
    return a.rot() * b;
}

/// transforms a direction by a frame (rigid affine transform)
template <>
constexpr inline vec3f transform_direction(const frame3f& a, const vec3f& b) {
    return a.rot() * b;
}

/// transforms a frame by a frame (rigid affine transform)
template <>
constexpr inline frame3f transform_frame(const frame3f& a, const frame3f& b) {
    return {a.rot() * b.rot(), a.rot() * b.pos() + a.pos()};
}

/// inverse transforms a point by a frame (rigid affine transform)
template <>
constexpr inline vec3f transform_point_inverse(
    const frame3f& a, const vec3f& b) {
    return (b - a.pos()) * a.rot();
}

/// inverse transforms a vector by a frame (rigid affine transform)
template <>
constexpr inline vec3f transform_vector_inverse(
    const frame3f& a, const vec3f& b) {
    return b * a.rot();
}

/// inverse transforms a direction by a frame (rigid affine transform)
template <>
constexpr inline vec3f transform_direction_inverse(
    const frame3f& a, const vec3f& b) {
    return b * a.rot();
}

#endif

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
    auto c = bbox<T, 3>{pos(a), pos(a)};
    // for all three axes
    for (auto i = 0; i < 3; i++) {
        // form extent by summing smaller and larger terms respectively
        for (auto j = 0; j < 3; j++) {
            auto e = rot(a)[j][i] * b[0][j];
            auto f = rot(a)[j][i] * b[1][j];
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

/// rotation matrix from axis-angle
template <typename T>
constexpr inline mat<T, 3, 3> rotation_mat3(const vec<T, 3>& axis, T angle) {
    auto s = sin(angle), c = cos(angle);
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
    return to_mat(translation_frame3(a));
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
    return to_mat(scaling_frame3(a));
}

/// rotation frame
template <typename T>
constexpr inline frame<T, 3> rotation_frame3(const vec<T, 3>& axis, T angle) {
    return {rotation_mat3(axis, angle), {0, 0, 0}};
}

/// rotation matrix
template <typename T>
constexpr inline mat<T, 4, 4> rotation_mat4(const mat<T, 3, 3>& rot) {
    return mat<T, 4, 4>{{rot.x.x, rot.x.y, rot.x.z, 0},
        {rot.y.x, rot.y.y, rot.y.z, 0}, {rot.z.x, rot.z.y, rot.z.z, 0},
        {0, 0, 0, 1}};
}

/// rotation matrix
template <typename T>
constexpr inline mat<T, 4, 4> rotation_mat4(const vec<T, 3>& axis, T angle) {
    return rotation_mat4(rotation_frame3(axis, angle));
}

/// quaternion axis-angle conversion
template <typename T>
constexpr inline vec<T, 4> rotation_axisangle4(const quat<T, 4>& a) {
    auto axis = normalize(vec<T, 3>{a.x, a.y, a.z});
    auto angle = std::acos(a.w) * 2;
    return {axis.x, axis.y, axis.z, angle};
}

/// axis-angle to quaternion
template <typename T>
constexpr inline quat<T, 4> rotation_quat4(const vec<T, 4>& axis_angle) {
    auto axis = vec<T, 3>{axis_angle.x, axis_angle.y, axis_angle.z};
    auto len = lenght(axis);
    auto angle = std::atan2(len, axis_angle.w);
    if (len)
        axis /= len;
    else
        axis = {0, 0, 1};
    return {axis.x, axis.y, axis.z, angle};
}

/// quaterion to matrix conversion
template <typename T>
constexpr inline mat<T, 3, 3> rotation_mat3(const quat<T, 4>& v) {
    return {
        {v[3] * v[3] + v[0] * v[0] - v[1] * v[1] - v[2] * v[2],
            (v[0] * v[1] + v[2] * v[3]) * 2, (v[2] * v[0] - v[1] * v[3]) * 2},
        {(v[0] * v[1] - v[2] * v[3]) * 2,
            v[3] * v[3] - v[0] * v[0] + v[1] * v[1] - v[2] * v[2],
            (v[1] * v[2] + v[0] * v[3]) * 2},
        {(v[2] * v[0] + v[1] * v[3]) * 2, (v[1] * v[2] - v[0] * v[3]) * 2,
            v[3] * v[3] - v[0] * v[0] - v[1] * v[1] + v[2] * v[2]}};
}

/// rotation matrix
template <typename T>
constexpr inline mat<T, 4, 4> rotation_mat4(const quat<T, 4>& v) {
    return rotation_mat4(rotation_mat3(v));
}

// http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm
/// matrix to quaternion
template <typename T>
constexpr inline quat<T, 4> rotation_quat4(const mat<T, 3, 3>& m_) {
    auto q = quat<T, 4>();
    auto m = transpose(m_);
#if 1
    auto trace = m[0][0] + m[1][1] + m[2][2];
    if (trace > 0) {
        float s = (T)0.5 / std::sqrt(trace + 1);
        q.w = (T)0.25 / s;
        q.x = (m[2][1] - m[1][2]) * s;
        q.y = (m[0][2] - m[2][0]) * s;
        q.z = (m[1][0] - m[0][1]) * s;
    } else {
        if (m[0][0] > m[1][1] && m[0][0] > m[2][2]) {
            float s = 2 * std::sqrt(max((T)0, 1 + m[0][0] - m[1][1] - m[2][2]));
            q.w = (m[2][1] - m[1][2]) / s;
            q.x = (T)0.25 * s;
            q.y = (m[0][1] + m[1][0]) / s;
            q.z = (m[0][2] + m[2][0]) / s;
        } else if (m[1][1] > m[2][2]) {
            float s = 2 * std::sqrt(max((T)0, 1 + m[1][1] - m[0][0] - m[2][2]));
            q.w = (m[0][2] - m[2][0]) / s;
            q.x = (m[0][1] + m[1][0]) / s;
            q.y = (T)0.25 * s;
            q.z = (m[1][2] + m[2][1]) / s;
        } else {
            float s = 2 * std::sqrt(max((T)0, 1 + m[2][2] - m[0][0] - m[1][1]));
            q.w = (m[1][0] - m[0][1]) / s;
            q.x = (m[0][2] + m[2][0]) / s;
            q.y = (m[1][2] + m[2][1]) / s;
            q.z = (T)0.25 * s;
        }
    }

#else
    q.w = std::sqrt(max(0, 1 + m[0][0] + m[1][1] + m[2][2])) / 2;
    q.x = std::sqrt(max(0, 1 + m[0][0] - m[1][1] - m[2][2])) / 2;
    q.y = std::sqrt(max(0, 1 - m[0][0] + m[1][1] - m[2][2])) / 2;
    q.z = std::sqrt(max(0, 1 - m[0][0] - m[1][1] + m[2][2])) / 2;
    Q.x = std::copysign(q.x, m[2][1] - m[1][2]);
    Q.y = std::copysign(q.y, m[0][2] - m[2][0]);
    Q.z = std::copysign(q.z, m[1][0] - m[0][1]);
#endif

    return q;
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
constexpr inline mat<T, 4, 4> ortho2d_mat4(T left, T right, T bottom, T top) {
    return ortho_mat4(left, right, bottom, top, -1, 1);
}

/// OpenGL/GLTF orthographic matrix
template <typename T>
constexpr inline mat<T, 4, 4> ortho_mat4(T xmag, T ymag, T near, T far) {
    return {{1 / xmag, 0, 0, 0}, {0, 1 / ymag, 0, 0},
        {0, 0, 2 / (near - far), 0}, {0, 0, (far + near) / (near - far), 1}};
}

/// OpenGL/GLTF perspective matrix
template <typename T>
constexpr inline mat<T, 4, 4> perspective_mat4(
    T fovy, T aspect, T near, T far) {
    auto tg = std::tan(fovy / 2);
    return {{1 / (aspect * tg), 0, 0, 0}, {0, 1 / tg, 0, 0},
        {0, 0, (far + near) / (near - far), -1},
        {0, 0, 2 * far * near / (near - far), 0}};
}

/// OpenGL/GLTF infinite perspective matrix
template <typename T>
constexpr inline mat<T, 4, 4> perspective_mat4(T fovy, T aspect, T near) {
    auto tg = std::tan(fovy / 2);
    return {{1 / (aspect * tg), 0, 0, 0}, {0, 1 / tg, 0, 0}, {0, 0, -1, -1},
        {0, 0, 2 * near, 0}};
}

/// Decompose an affine matrix into translation, rotation, scale.
/// Assumes there is no shear and the matrix is affine.
template <typename T>
constexpr inline void decompose_mat4(const mat<T, 4, 4>& m,
    vec<T, 3>& translation, mat<T, 3, 3>& rotation, vec<T, 3>& scale) {
    translation = {m.w.x, m.w.y, m.w.z};
    rotation.x = {m.x.x, m.x.y, m.x.z};
    rotation.y = {m.y.x, m.y.y, m.y.z};
    rotation.z = {m.z.x, m.z.y, m.z.z};
    scale = {length(rotation.x), length(rotation.y), length(rotation.z)};
    rotation = {
        normalize(rotation.x), normalize(rotation.y), normalize(rotation.z)};
}

/// Decompose an affine matrix into translation, rotation, scale.
/// Assumes there is no shear and the matrix is affine.
template <typename T>
constexpr inline void decompose_mat4(const mat<T, 4, 4>& m,
    vec<T, 3>& translation, quat<T, 4>& rotation, vec<T, 3>& scale) {
    auto rot_matrix = mat<T, 3, 3>();
    decompose_mat4(m, translation, rot_matrix, scale);
    rotation = to_quat4(rotation_mat4(rot_matrix));
}

/// Decompose an affine matrix into translation, rotation, scale.
/// Assumes there is no shear and the matrix is affine.
template <typename T>
constexpr inline mat<T, 4, 4> compose_mat4(const vec<T, 3>& translation,
    const mat<T, 3, 3>& rotation, const vec<T, 3>& scale) {
    return translation_mat4(translation) * scaling_mat4(scale) *
           rotation_mat4(rotation);
}

/// Decompose an affine matrix into translation, rotation, scale.
/// Assumes there is no shear and the matrix is affine.
template <typename T>
constexpr inline mat<T, 4, 4> compose_mat4(const vec<T, 3>& translation,
    const quat<T, 4>& rotation, const vec<T, 3>& scale) {
    return translation_mat4(translation) * scaling_mat4(scale) *
           rotation_mat4(rotation);
}

// -----------------------------------------------------------------------------
// GEOMETRY UTILITIES
// -----------------------------------------------------------------------------

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

//
// Triangle tangent and bitangent from uv (not othornormalized with themselfves
// not the normal). Follows the definition in
// http://www.terathon.com/code/tangent.html and
// https://gist.github.com/aras-p/2843984
template <typename T>
constexpr inline std::pair<vec<T, 3>, vec<T, 3>> triangle_tangents_fromuv(
    const vec<T, 3>& v0, const vec<T, 3>& v1, const vec<T, 3>& v2,
    const vec<T, 2>& uv0, const vec<T, 2>& uv1, const vec<T, 2>& uv2) {
    // normal points up from texture space
    auto p = v1 - v0;
    auto q = v2 - v0;
    auto s = vec<T, 2>{uv1[0] - uv0[0], uv2[0] - uv0[0]};
    auto t = vec<T, 2>{uv1[1] - uv0[1], uv2[1] - uv0[1]};
    auto div = s[0] * t[1] - s[1] * t[0];

    if (div != 0) {
        auto tu = vec<T, 3>{t[1] * p[0] - t[0] * q[0],
                      t[1] * p[1] - t[0] * q[1], t[1] * p[2] - t[0] * q[2]} /
                  div;
        auto tv = vec<T, 3>{s[0] * q[0] - s[1] * p[0],
                      s[0] * q[1] - s[1] * p[1], s[0] * q[2] - s[1] * p[2]} /
                  div;
        return {tu, tv};
    } else {
        return {{1, 0, 0}, {0, 1, 0}};
    }
}

/// triangle baricentric interpolation
template <typename T, typename T1>
constexpr inline T blerp(const T& a, const T& b, const T& c, const T1& w) {
    return a * w[0] + b * w[1] + c * w[2];
}

///
/// Compute smoothed tangents (for lines).
///
/// Parameters:
/// - nverts/pos: array pf vertex positions
/// - npoints/points: array of point indices
/// - nlines/lines: array of point indices
/// - ntriangles/triangles: array of point indices
/// - weighted: whether to use area weighting (typically true)
///
/// Out Parameters:
/// - tang: preallocated array of computed normals
///
inline void compute_tangents(int nlines, const vec2i* lines, int nverts,
    const vec3f* pos, vec3f* tang, bool weighted = true) {
    // clear tangents
    for (auto i = 0; i < nverts; i++) tang[i] = zero3f;

    // handle lines
    for (auto i = 0; i < nlines; i++) {
        auto line = lines[i];
        auto n = pos[line.y] - pos[line.x];
        if (!weighted) n = normalize(n);
        tang[line.x] += n;
        tang[line.y] += n;
    }

    // normalize result
    for (auto i = 0; i < nverts; i++) tang[i] = normalize(tang[i]);
}

///
/// Compute smoothed tangents.
///
/// Parameters:
/// - nverts/pos: array pf vertex positions
/// - npoints/points: array of point indices
/// - nlines/lines: array of point indices
/// - ntriangles/triangles: array of point indices
/// - weighted: whether to use area weighting (typically true)
///
/// Out Parameters:
/// - tang: array of computed tangents
///
inline void compute_tangents(const std::vector<vec2i>& lines,
    const std::vector<vec3f>& pos, std::vector<vec3f>& tang,
    bool weighted = true) {
    tang.resize(pos.size());
    compute_tangents(lines.size(), lines.data(), pos.size(), pos.data(),
        tang.data(), weighted);
}

///
/// Compute smoothed normals.
///
/// Parameters:
/// - nverts/pos: array pf vertex positions
/// - npoints/points: array of point indices
/// - nlines/lines: array of point indices
/// - ntriangles/triangles: array of point indices
/// - weighted: whether to use area weighting (typically true)
///
/// Out Parameters:
/// - norm: preallocated array of computed normals
///
inline void compute_normals(int ntriangles, const vec3i* triangles, int nverts,
    const vec3f* pos, vec3f* norm, bool weighted = true) {
    // clear normals
    for (auto i = 0; i < nverts; i++) norm[i] = zero3f;

    // handle triangles
    for (auto i = 0; i < ntriangles; i++) {
        auto triangle = triangles[i];
        auto n = cross(pos[triangle.y] - pos[triangle.x],
            pos[triangle.z] - pos[triangle.x]);
        if (!weighted) n = normalize(n);
        norm[triangle.x] += n;
        norm[triangle.y] += n;
        norm[triangle.z] += n;
    }

    // normalize result
    for (auto i = 0; i < nverts; i++) norm[i] = normalize(norm[i]);
}

///
/// Compute smoothed normals.
///
/// Parameters:
/// - nverts/pos: array pf vertex positions
/// - npoints/points: array of point indices
/// - nlines/lines: array of point indices
/// - ntriangles/triangles: array of point indices
/// - weighted: whether to use area weighting (typically true)
///
/// Out Parameters:
/// - norm: array of computed normals
///
inline void compute_normals(const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, std::vector<vec3f>& norm,
    bool weighted = true) {
    norm.resize(pos.size());
    compute_normals(triangles.size(), triangles.data(), pos.size(), pos.data(),
        norm.data(), weighted);
}

///
/// Compute tangent frame for triangle mesh. Tangent space is defined by
/// a four component vector. The first three components are the tangent
/// with respect to the U texcoord. The fourth component is the sign of the
/// tangent wrt the V texcoord. Tangent frame is useful in normal mapping.
///
/// Parameters:
/// - nverts/pos: array pf vertex positions
/// - ntriangles/triangles: array of point indices
/// - weighted: whether to use area weighting (typically true)
///
/// Out Parameters:
/// - tangsp: preallocated array of computed tangent space
///
inline void compute_tangent_frame(int ntriangles, const vec3i* triangles,
    int nverts, const vec3f* pos, const vec3f* norm, const vec2f* texcoord,
    vec4f* tangsp, bool weighted = true) {
    auto tangu = std::vector<vec3f>(nverts, {0, 0, 0});
    auto tangv = std::vector<vec3f>(nverts, {0, 0, 0});

    for (auto i = 0; i < ntriangles; i++) {
        auto t = triangles[i];
        auto tutv = triangle_tangents_fromuv(pos[t.x], pos[t.y], pos[t.z],
            texcoord[t.x], texcoord[t.y], texcoord[t.z]);
        if (!weighted) tutv = {normalize(tutv.first), normalize(tutv.second)};
        tangu[t.x] += tutv.first;
        tangu[t.y] += tutv.first;
        tangu[t.z] += tutv.first;
        tangv[t.x] += tutv.second;
        tangv[t.y] += tutv.second;
        tangv[t.z] += tutv.second;
    }

    for (auto i = 0; i < nverts; i++) {
        tangu[i] = normalize(tangu[i]);
        tangv[i] = normalize(tangv[i]);
    }

    for (auto i = 0; i < nverts; i++) {
        tangu[i] = orthonormalize(tangu[i], norm[i]);
        auto s = (dot(cross(norm[i], tangu[i]), tangv[i]) < 0) ? -1.0f : 1.0f;
        tangsp[i] = {tangu[i][0], tangu[i][1], tangu[i][2], s};
    }
}

///
/// Compute tangent frame for triangle mesh. Tangent space is defined by
/// a four component vector. The first three components are the tangent
/// with respect to the U texcoord. The fourth component is the sign of the
/// tangent wrt the V texcoord. Tangent frame is useful in normal mapping.
///
/// Parameters:
/// - nverts/pos: array pf vertex positions
/// - ntriangles/triangles: array of point indices
/// - weighted: whether to use area weighting (typically true)
///
/// Out Parameters:
/// - tangsp: array of computed tangent space
///
inline void compute_tangent_frame(const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, const std::vector<vec3f>& norm,
    const std::vector<vec2f>& texcoord, std::vector<vec4f>& tangsp,
    bool weighted = true) {
    tangsp.resize(tangsp.size());
    compute_tangent_frame(triangles.size(), triangles.data(), pos.size(),
        pos.data(), norm.data(), texcoord.data(), tangsp.data(), weighted);
}

/// Apply skinning
inline void compute_skinning(int nverts, const vec3f* pos, const vec3f* norm,
    const vec4f* weights, const vec4i* joints, const mat4f* xforms,
    vec3f* skinned_pos, vec3f* skinned_norm) {
    for (auto i = 0; i < nverts; i++) {
        skinned_pos[i] =
            transform_point(xforms[joints[i].x], pos[i]) * weights[i].x +
            transform_point(xforms[joints[i].y], pos[i]) * weights[i].y +
            transform_point(xforms[joints[i].z], pos[i]) * weights[i].z +
            transform_point(xforms[joints[i].w], pos[i]) * weights[i].w;
    }
    for (auto i = 0; i < nverts; i++) {
        skinned_norm[i] = normalize(
            transform_direction(xforms[joints[i].x], norm[i]) * weights[i].x +
            transform_direction(xforms[joints[i].y], norm[i]) * weights[i].y +
            transform_direction(xforms[joints[i].z], norm[i]) * weights[i].z +
            transform_direction(xforms[joints[i].w], norm[i]) * weights[i].w);
    }
}

/// Apply skinning
inline void compute_skinning(int nverts, const vec3f* pos, const vec3f* norm,
    const vec4f* weights, const vec4i* joints, const frame3f* xforms,
    vec3f* skinned_pos, vec3f* skinned_norm) {
    for (auto i = 0; i < nverts; i++) {
        skinned_pos[i] =
            transform_point(xforms[joints[i].x], pos[i]) * weights[i].x +
            transform_point(xforms[joints[i].y], pos[i]) * weights[i].y +
            transform_point(xforms[joints[i].z], pos[i]) * weights[i].z +
            transform_point(xforms[joints[i].w], pos[i]) * weights[i].w;
    }
    for (auto i = 0; i < nverts; i++) {
        skinned_norm[i] = normalize(
            transform_direction(xforms[joints[i].x], norm[i]) * weights[i].x +
            transform_direction(xforms[joints[i].y], norm[i]) * weights[i].y +
            transform_direction(xforms[joints[i].z], norm[i]) * weights[i].z +
            transform_direction(xforms[joints[i].w], norm[i]) * weights[i].w);
    }
}

/// Apply skinning
inline void compute_skinning(const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec4f>& weights,
    const std::vector<vec4i>& joints, const std::vector<mat4f>& xforms,
    std::vector<vec3f>& skinned_pos, std::vector<vec3f>& skinned_norm) {
    skinned_pos.resize(pos.size());
    skinned_norm.resize(norm.size());
    compute_skinning(pos.size(), pos.data(), norm.data(), weights.data(),
        joints.data(), xforms.data(), skinned_pos.data(), skinned_norm.data());
}

/// Apply skinning
inline void compute_skinning(const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec4f>& weights,
    const std::vector<vec4i>& joints, const std::vector<frame3f>& xforms,
    std::vector<vec3f>& skinned_pos, std::vector<vec3f>& skinned_norm) {
    skinned_pos.resize(pos.size());
    skinned_norm.resize(norm.size());
    compute_skinning(pos.size(), pos.data(), norm.data(), weights.data(),
        joints.data(), xforms.data(), skinned_pos.data(), skinned_norm.data());
}

/// Apply skinning as specified in Khronos glTF
inline void compute_matrix_skinning(int nverts, const vec3f* pos,
    const vec3f* norm, const vec4f* weights, const vec4i* joints,
    const mat4f* xforms, vec3f* skinned_pos, vec3f* skinned_norm) {
    for (auto i = 0; i < nverts; i++) {
        auto xform = xforms[joints[i].x] * weights[i].x +
                     xforms[joints[i].y] * weights[i].y +
                     xforms[joints[i].z] * weights[i].z +
                     xforms[joints[i].w] * weights[i].w;
        skinned_pos[i] = transform_point(xform, pos[i]);
        skinned_norm[i] = normalize(transform_direction(xform, norm[i]));
    }
}

/// Apply skinning as specified in Khronos glTF
inline void compute_matrix_skinning(const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec4f>& weights,
    const std::vector<vec4i>& joints, const std::vector<mat4f>& xforms,
    std::vector<vec3f>& skinned_pos, std::vector<vec3f>& skinned_norm) {
    skinned_pos.resize(pos.size());
    skinned_norm.resize(norm.size());
    compute_matrix_skinning(pos.size(), pos.data(), norm.data(), weights.data(),
        joints.data(), xforms.data(), skinned_pos.data(), skinned_norm.data());
}

///
/// Generate a parametric surface with callbacks.
///
/// Parameters:
/// - usteps: subdivisions in u
/// - vsteps: subdivisions in v
/// - pos_fn: pos callbacks (vec2f -> vec3f)
/// - norm_fn: norm callbacks (vec2f -> vec3f)
/// - texcoord_fn: texcoord callbacks (vec2f -> vec2f)
///
/// Out Parameters:
/// - triangles: element array
/// - pos/norm/texcoord: vertex position/normal/texcoords
///
template <typename PosFunc, typename NormFunc, typename TexcoordFunc>
inline void make_triangles(int usteps, int vsteps,
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    const PosFunc& pos_fn, const NormFunc& norm_fn,
    const TexcoordFunc& texcoord_fn) {
    auto vid = [usteps](int i, int j) { return j * (usteps + 1) + i; };
    pos.resize((usteps + 1) * (vsteps + 1));
    norm.resize((usteps + 1) * (vsteps + 1));
    texcoord.resize((usteps + 1) * (vsteps + 1));
    for (auto j = 0; j <= vsteps; j++) {
        for (auto i = 0; i <= usteps; i++) {
            auto uv = vec2f{i / (float)usteps, j / (float)vsteps};
            pos[vid(i, j)] = pos_fn(uv);
            norm[vid(i, j)] = norm_fn(uv);
            texcoord[vid(i, j)] = texcoord_fn(uv);
        }
    }

    triangles.resize(usteps * vsteps * 2);
    for (auto j = 0; j < vsteps; j++) {
        for (auto i = 0; i < usteps; i++) {
            auto& f1 = triangles[(j * usteps + i) * 2 + 0];
            auto& f2 = triangles[(j * usteps + i) * 2 + 1];
            if ((i + j) % 2) {
                f1 = {vid(i, j), vid(i + 1, j), vid(i + 1, j + 1)};
                f2 = {vid(i + 1, j + 1), vid(i, j + 1), vid(i, j)};
            } else {
                f1 = {vid(i, j), vid(i + 1, j), vid(i, j + 1)};
                f2 = {vid(i + 1, j + 1), vid(i, j + 1), vid(i + 1, j)};
            }
        }
    }
}

///
/// Generate parametric lines with callbacks.
///
/// Parameters:
/// - usteps: subdivisions in u
/// - num: number of lines
/// - pos_fn: pos callbacks (vec2f -> vec3f)
/// - tang_fn: tangent callbacks (vec2f -> vec3f)
/// - texcoord_fn: texcoord callbacks (vec2f -> vec2f)
/// - radius_fn: radius callbacks (vec2f -> float)
///
/// Out Parameters:
/// - lines: element array
/// - pos/tang/texcoord/radius: vertex position/tangent/texcoords/radius
///
template <typename PosFunc, typename TangFunc, typename TexcoordFunc,
    typename RadiusFunc>
inline void make_lines(int usteps, int num, std::vector<vec2i>& lines,
    std::vector<vec3f>& pos, std::vector<vec3f>& tang,
    std::vector<vec2f>& texcoord, std::vector<float>& radius,
    const PosFunc& pos_fn, const TangFunc& tang_fn,
    const TexcoordFunc& texcoord_fn, const RadiusFunc& radius_fn) {
    auto vid = [usteps](int i, int j) { return j * (usteps + 1) + i; };
    pos.resize((usteps + 1) * num);
    tang.resize((usteps + 1) * num);
    texcoord.resize((usteps + 1) * num);
    radius.resize((usteps + 1) * num);
    for (auto j = 0; j < num; j++) {
        for (auto i = 0; i <= usteps; i++) {
            auto uv = vec2f{i / (float)usteps, j / (float)(num - 1)};
            pos[vid(i, j)] = pos_fn(uv);
            tang[vid(i, j)] = tang_fn(uv);
            texcoord[vid(i, j)] = texcoord_fn(uv);
            radius[vid(i, j)] = radius_fn(uv);
        }
    }

    lines.resize(usteps * num);
    for (int j = 0; j < num; j++) {
        for (int i = 0; i < usteps; i++) {
            lines[j * usteps + i] = {vid(i, j), vid(i + 1, j)};
        }
    }
}

///
/// Generate a parametric point set. Mostly here for completeness.
///
/// Parameters:
/// - num: number of points
/// - pos_fn: pos callbacks (float -> vec3f)
/// - norm_fn: norm callbacks (float -> vec3f)
/// - texcoord_fn: texcoord callbacks (float -> vec2f)
/// - radius_fn: radius callbacks (float -> float)
///
/// Out Parameters:
/// - points: element array
/// - pos/norm/texcoord/radius: vertex position/normal/texcoords/radius
///
template <typename PosFunc, typename NormFunc, typename TexcoordFunc,
    typename RadiusFunc>
inline void make_points(int num, std::vector<int>& points,
    std::vector<vec3f>& pos, std::vector<vec3f>& norm,
    std::vector<vec2f>& texcoord, std::vector<float>& radius,
    const PosFunc& pos_fn, const NormFunc& norm_fn,
    const TexcoordFunc& texcoord_fn, const RadiusFunc& radius_fn) {
    pos.resize(num);
    norm.resize(num);
    texcoord.resize(num);
    radius.resize(num);
    for (auto i = 0; i < num; i++) {
        auto u = i / (float)i;
        pos[i] = pos_fn(u);
        norm[i] = norm_fn(u);
        texcoord[i] = texcoord_fn(u);
        radius[i] = radius_fn(u);
    }

    points.resize(num);
    for (auto i = 0; i < num; i++) points[i] = i;
}

///
/// Merge a triangle mesh into another.
///
inline void merge_triangles(std::vector<ym::vec3i>& triangles,
    std::vector<ym::vec3f>& pos, std::vector<ym::vec3f>& norm,
    std::vector<ym::vec2f>& texcoord, const std::vector<ym::vec3i>& mtriangles,
    const std::vector<ym::vec3f>& mpos, const std::vector<ym::vec3f>& mnorm,
    const std::vector<ym::vec2f>& mtexcoord) {
    auto o = (int)pos.size();
    for (auto t : mtriangles)
        triangles.push_back({t[0] + o, t[1] + o, t[2] + o});
    for (auto p : mpos) pos.push_back(p);
    for (auto n : mnorm) norm.push_back(n);
    for (auto t : mtexcoord) texcoord.push_back(t);
}

// -----------------------------------------------------------------------------
// STANDARD SHAPES
// -----------------------------------------------------------------------------

///
/// Make a sphere.
///
inline void make_uvsphere(int usteps, int vsteps, std::vector<vec3i>& triangles,
    std::vector<vec3f>& pos, std::vector<vec3f>& norm,
    std::vector<vec2f>& texcoord) {
    return make_triangles(usteps, vsteps, triangles, pos, norm, texcoord,
        [](const vec2f& uv) {
            auto a = vec2f{2 * pif * uv[0], pif * (1 - uv[1])};
            return vec3f{
                cos(a[0]) * sin(a[1]), sin(a[0]) * sin(a[1]), cos(a[1])};
        },
        [](const vec2f& uv) {
            auto a = vec2f{2 * pif * uv[0], pif * (1 - uv[1])};
            return vec3f{
                cos(a[0]) * sin(a[1]), sin(a[0]) * sin(a[1]), cos(a[1])};
        },
        [](const vec2f& uv) { return uv; });
}

///
/// Make a sphere.
///
inline void make_uvhemisphere(int usteps, int vsteps,
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord) {
    return ym::make_triangles(usteps, vsteps, triangles, pos, norm, texcoord,
        [](const vec2f& uv) {
            auto a = vec2f{2 * pif * uv[0], pif * 0.5f * (1 - uv[1])};
            return vec3f{
                cos(a[0]) * sin(a[1]), sin(a[0]) * sin(a[1]), cos(a[1])};
        },
        [](const vec2f& uv) {
            auto a = vec2f{2 * pif * uv[0], pif * 0.5f * (1 - uv[1])};
            return vec3f{
                cos(a[0]) * sin(a[1]), sin(a[0]) * sin(a[1]), cos(a[1])};
        },
        [](const vec2f& uv) { return uv; });
}

///
/// Make an inside-out sphere.
///
inline void make_uvflippedsphere(int usteps, int vsteps,
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord) {
    return ym::make_triangles(usteps, vsteps, triangles, pos, norm, texcoord,
        [](const vec2f& uv) {
            auto a = vec2f{2 * pif * uv[0], pif * uv[1]};
            return vec3f{
                cos(a[0]) * sin(a[1]), sin(a[0]) * sin(a[1]), cos(a[1])};
        },
        [](const vec2f& uv) {
            auto a = vec2f{2 * pif * uv[0], pif * uv[1]};
            return vec3f{
                -cos(a[0]) * sin(a[1]), -sin(a[0]) * sin(a[1]), -cos(a[1])};
        },
        [](const vec2f& uv) {
            return vec2f{uv.x, 1 - uv.y};
        });
}

///
/// Make an inside-out hemisphere
///
inline void make_uvflippedhemisphere(int usteps, int vsteps,
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord) {
    return ym::make_triangles(usteps, vsteps, triangles, pos, norm, texcoord,
        [](const vec2f& uv) {
            auto a = vec2f{2 * pif * uv[0], pif * (0.5f + 0.5f * uv[1])};
            return vec3f{
                cos(a[0]) * sin(a[1]), sin(a[0]) * sin(a[1]), cos(a[1])};
        },
        [](const vec2f& uv) {
            auto a = vec2f{2 * pif * uv[0], pif * uv[1]};
            return vec3f{
                -cos(a[0]) * sin(a[1]), -sin(a[0]) * sin(a[1]), -cos(a[1])};
        },
        [](const vec2f& uv) {
            return vec2f{uv.x, 1 - uv.y};
        });
}

///
/// Make a quad.
///
inline void make_uvquad(int usteps, int vsteps, std::vector<vec3i>& triangles,
    std::vector<vec3f>& pos, std::vector<vec3f>& norm,
    std::vector<vec2f>& texcoord) {
    return make_triangles(usteps, vsteps, triangles, pos, norm, texcoord,
        [](const vec2f& uv) {
            return vec3f{(-1 + uv[0] * 2), (-1 + uv[1] * 2), 0};
        },
        [](const vec2f& uv) {
            return vec3f{0, 0, 1};
        },
        [](const vec2f& uv) { return uv; });
}

///
/// Make a quad.
///
inline void make_uvcube(int usteps, int vsteps, std::vector<vec3i>& triangles,
    std::vector<vec3f>& pos, std::vector<vec3f>& norm,
    std::vector<vec2f>& texcoord) {
    ym::frame3f frames[6] = {
        ym::frame3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 1}},
        ym::frame3f{{-1, 0, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, -1}},
        ym::frame3f{{-1, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 0}},
        ym::frame3f{{1, 0, 0}, {0, 0, 1}, {0, -1, 0}, {0, -1, 0}},
        ym::frame3f{{0, 1, 0}, {0, 0, 1}, {1, 0, 0}, {1, 0, 0}},
        ym::frame3f{{0, -1, 0}, {0, 0, 1}, {-1, 0, 0}, {-1, 0, 0}}};
    std::vector<vec3f> quad_pos, quad_norm;
    std::vector<vec2f> quad_texcoord;
    std::vector<ym::vec3i> quad_triangles;
    make_uvquad(
        usteps, vsteps, quad_triangles, quad_pos, quad_norm, quad_texcoord);
    for (auto i = 0; i < 6; i++) {
        pos.insert(pos.end(), quad_pos.begin(), quad_pos.end());
        norm.insert(norm.end(), quad_norm.begin(), quad_norm.end());
        texcoord.insert(
            texcoord.end(), quad_texcoord.begin(), quad_texcoord.end());
        triangles.insert(
            triangles.end(), quad_triangles.begin(), quad_triangles.end());
    }
    auto quad_faces = quad_triangles.size(), quad_verts = quad_pos.size();
    for (auto i = 0; i < 6; i++) {
        for (auto j = quad_verts * i; j < quad_verts * (i + 1); j++)
            pos[j] = ym::transform_point(frames[i], pos[j]);
        for (auto j = quad_verts * i; j < quad_verts * (i + 1); j++)
            norm[j] = ym::transform_direction(frames[i], norm[j]);
        for (auto j = quad_faces * i; j < quad_faces * (i + 1); j++) {
            triangles[j].x += quad_verts * i;
            triangles[j].y += quad_verts * i;
            triangles[j].z += quad_verts * i;
        }
    }
}

///
/// Make a quad.
///
inline void make_uvspherecube(int usteps, int vsteps,
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord) {
    make_uvcube(usteps, vsteps, triangles, pos, norm, texcoord);
    for (auto i = 0; i < pos.size(); i++) {
        pos[i] = normalize(pos[i]);
        norm[i] = normalize(pos[i]);
    }
}

///
/// Make a quad.
///
inline void make_uvspherizedcube(int usteps, int vsteps, float radius,
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord) {
    ym::make_uvcube(usteps, vsteps, triangles, pos, norm, texcoord);
    for (auto i = 0; i < pos.size(); i++) {
        norm[i] = ym::normalize(pos[i]);
        pos[i] *= 1 - radius;
        pos[i] += norm[i] * radius;
    }
    ym::compute_normals(triangles, pos, norm, true);
}

///
/// Make a quad.
///
inline void make_uvflipcapsphere(int usteps, int vsteps, float radius,
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord) {
    ym::make_uvsphere(usteps, vsteps, triangles, pos, norm, texcoord);
    for (auto i = 0; i < pos.size(); i++) {
        if (pos[i][2] > radius) {
            pos[i][2] = 2 * radius - pos[i][2];
            norm[i][0] = -norm[i][0];
            norm[i][1] = -norm[i][1];
        } else if (pos[i][2] < -radius) {
            pos[i][2] = -2 * radius - pos[i][2];
            norm[i][0] = -norm[i][0];
            norm[i][1] = -norm[i][1];
        }
    }
}

///
/// Make a quad.
///
inline void make_uvcutsphere(int usteps, int vsteps, float radius,
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord) {
    return make_triangles(usteps, vsteps, triangles, pos, norm, texcoord,
        [radius](const vec2f& uv) {
            auto p = 1 - std::acos(radius) / pif;
            auto a = vec2f{2 * pif * uv[0], pif * (1 - p * uv[1])};
            return vec3f{std::cos(a[0]) * std::sin(a[1]),
                std::sin(a[0]) * std::sin(a[1]), std::cos(a[1])};
        },
        [radius](const vec2f& uv) {
            auto p = 1 - std::acos(radius) / pif;
            auto a = vec2f{2 * pif * uv[0], pif * (1 - p * uv[1])};
            return vec3f{std::cos(a[0]) * std::sin(a[1]),
                std::sin(a[0]) * std::sin(a[1]), std::cos(a[1])};
        },
        [](const vec2f& uv) { return uv; });
}

///
/// Make a quad.
///
inline void make_uvflippedcutsphere(int usteps, int vsteps, float radius,
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord) {
    return make_triangles(usteps, vsteps, triangles, pos, norm, texcoord,
        [radius](const vec2f& uv) {
            auto p = 1 - acos(radius) / pif;
            auto a = vec2f{2 * pif * uv[0], pif * ((1 - p) + p * uv[1])};
            return vec3f{
                cos(a[0]) * sin(a[1]), sin(a[0]) * sin(a[1]), cos(a[1])};
        },
        [radius](const vec2f& uv) {
            auto p = 1 - acos(radius) / pif;
            auto a = vec2f{2 * pif * uv[0], pif * ((1 - p) + p * uv[1])};
            return vec3f{
                -cos(a[0]) * sin(a[1]), -sin(a[0]) * sin(a[1]), -cos(a[1])};
        },
        [](const vec2f& uv) {
            return vec2f{uv[0], (1 - uv[1])};
        });
}

// -----------------------------------------------------------------------------
// RAY-PRIMITIVE INTERSECTION FUNCTIONS
// -----------------------------------------------------------------------------

///
/// Intersect a ray with a point (approximate)
///
/// Parameters:
/// - ray: ray origin and direction, parameter min, max range
/// - p: point position
/// - r: point radius
///
/// Out Parameters:
/// - ray_t: ray parameter at the intersection point
/// - euv: primitive uv ( {0,0} for points )
///
/// Returns:
/// - whether the intersection occurred
///
/// Iplementation Notes:
/// - out Parameters and only writtent o if an intersection occurs
/// - algorithm finds the closest point on the ray segment to the point and
///    test their distance with the point radius
/// - based on http://geomalgorithms.com/a02-lines.html.
///
inline bool intersect_point(
    const ray3f& ray, const vec3f& p, float r, float& ray_t) {
    // find parameter for line-point minimum distance
    auto w = p - ray.o;
    auto t = dot(w, ray.d) / dot(ray.d, ray.d);

    // exit if not within bounds
    if (t < ray.tmin || t > ray.tmax) return false;

    // test for line-point distance vs point radius
    auto rp = eval(ray, t);
    auto prp = p - rp;
    if (dot(prp, prp) > r * r) return false;

    // intersection occurred: set params and exit
    ray_t = t;

    return true;
}

///
/// Intersect a ray with a line
///
/// Parameters:
/// - ray: ray origin and direction, parameter min, max range
/// - v0, v1: line segment points
/// - r0, r1: line segment radia
///
/// Out Parameters:
/// - ray_t: ray parameter at the intersection point
/// - euv: euv[0] is the line parameter at the intersection ( euv[1] is zero )
///
/// Returns:
/// - whether the intersection occurred
///
/// Notes:
/// - out Parameters and only writtent o if an intersection occurs
/// - algorithm find the closest points on line and ray segment and test
///   their distance with the line radius at that location
/// - based on http://geomalgorithms.com/a05-intersect-1.html
/// - based on http://geomalgorithms.com/a07-distance.html#
///     dist3D_Segment_to_Segment
///
inline bool intersect_line(const ray3f& ray, const vec3f& v0, const vec3f& v1,
    float r0, float r1, float& ray_t, vec2f& euv) {
    // setup intersection params
    auto u = ray.d;
    auto v = v1 - v0;
    auto w = ray.o - v0;

    // compute values to solve a linear system
    auto a = dot(u, u);
    auto b = dot(u, v);
    auto c = dot(v, v);
    auto d = dot(u, w);
    auto e = dot(v, w);
    auto det = a * c - b * b;

    // check determinant and exit if lines are parallel
    // (could use EPSILONS if desired)
    if (det == 0) return false;

    // compute Parameters on both ray and segment
    auto t = (b * e - c * d) / det;
    auto s = (a * e - b * d) / det;

    // exit if not within bounds
    if (t < ray.tmin || t > ray.tmax) return false;

    // clamp segment param to segment corners
    s = clamp(s, (float)0, (float)1);

    // compute segment-segment distance on the closest points
    auto p0 = eval(ray, t);
    auto p1 = eval(ray3f{v0, v1 - v0}, s);
    auto p01 = p0 - p1;

    // check with the line radius at the same point
    auto r = r0 * (1 - s) + r1 * s;
    if (dot(p01, p01) > r * r) return false;

    // intersection occurred: set params and exit
    ray_t = t;
    euv = {1 - s, s};

    return true;
}

///
/// Intersect a ray with a triangle
///
/// Parameters:
/// - ray: ray origin and direction, parameter min, max range
/// - v0, v1, v2: triangle vertices
///
/// Out Parameters:
/// - ray_t: ray parameter at the intersection point
/// - euv: baricentric coordinates of the intersection
///
/// Returns:
/// - whether the intersection occurred
///
/// Notes:
/// - out Parameters and only writtent o if an intersection occurs
/// - algorithm based on Muller-Trombone intersection test
///
inline bool intersect_triangle(const ray3f& ray, const vec3f& v0,
    const vec3f& v1, const vec3f& v2, float& ray_t, vec3f& euv) {
    // compute triangle edges
    auto edge1 = v1 - v0;
    auto edge2 = v2 - v0;

    // compute determinant to solve a linear system
    auto pvec = cross(ray.d, edge2);
    auto det = dot(edge1, pvec);

    // check determinant and exit if triangle and ray are parallel
    // (could use EPSILONS if desired)
    if (det == 0) return false;
    auto inv_det = 1.0f / det;

    // compute and check first bricentric coordinated
    auto tvec = ray.o - v0;
    auto u = dot(tvec, pvec) * inv_det;
    if (u < 0 || u > 1) return false;

    // compute and check second bricentric coordinated
    auto qvec = cross(tvec, edge1);
    auto v = dot(ray.d, qvec) * inv_det;
    if (v < 0 || u + v > 1) return false;

    // compute and check ray parameter
    auto t = dot(edge2, qvec) * inv_det;
    if (t < ray.tmin || t > ray.tmax) return false;

    // intersection occurred: set params and exit
    ray_t = t;
    euv = {1 - u - v, u, v};

    return true;
}

///
/// Intersect a ray with a tetrahedron. Note that we consider only intersection
/// wiht the tetrahedra surface and discount intersction with the interior.
///
/// Parameters:
/// - ray: ray to intersect with
/// - v0, v1, v2: triangle vertices
///
/// Out Parameters:
/// - ray_t: ray parameter at the intersection point
/// - euv: baricentric coordinates of the intersection
///
/// Returns:
/// - whether the intersection occurred
///
/// TODO: check order
/// TODO: uv
///
inline bool intersect_tetrahedron(const ray3f& ray_, const vec3f& v0,
    const vec3f& v1, const vec3f& v2, const vec3f& v3, float& ray_t,
    vec4f& euv) {
    // check intersction for each face
    auto hit = false;
    auto ray = ray_;
    auto tuv = zero3f;
    if (intersect_triangle(ray, v0, v1, v2, ray_t, tuv)) {
        hit = true;
        ray.tmax = ray_t;
    }
    if (intersect_triangle(ray, v0, v1, v3, ray_t, tuv)) {
        hit = true;
        ray.tmax = ray_t;
    }
    if (intersect_triangle(ray, v0, v2, v3, ray_t, tuv)) {
        hit = true;
        ray.tmax = ray_t;
    }
    if (intersect_triangle(ray, v1, v2, v3, ray_t, tuv)) {
        hit = true;
        ray.tmax = ray_t;
    }

    return hit;
}

///
/// Intersect a ray with a axis-aligned bounding box
///
/// Parameters:
/// - ray: ray to intersect with
/// - bbox: bounding box min/max bounds
///
/// Returns:
/// - whether the intersection occurred
///
inline bool intersect_check_bbox(const ray3f& ray, const bbox3f& bbox) {
    // set up convenient pointers for looping over axes
    auto tmin = ray.tmin, tmax = ray.tmax;

    // for each axis, clip intersection against the bounding planes
    for (int i = 0; i < 3; i++) {
        // determine intersection ranges
        auto invd = 1.0f / ray.d[i];
        auto t0 = (bbox[0][i] - ray.o[i]) * invd;
        auto t1 = (bbox[1][i] - ray.o[i]) * invd;
        // flip based on range directions
        if (invd < 0.0f) {
            float a = t0;
            t0 = t1;
            t1 = a;
        }
        // clip intersection
        tmin = t0 > tmin ? t0 : tmin;
        tmax = t1 < tmax ? t1 : tmax;
        // if intersection is empty, exit
        if (tmin > tmax) return false;
    }

    // passed all planes, then intersection occurred
    return true;
}

///
/// Min/max used in BVH traversal. Copied here since the traversal code relies
/// on the specific behaviour wrt NaNs.
///
template <typename T>
static inline const T& _safemin(const T& a, const T& b) {
    return (a < b) ? a : b;
}
///
/// Min/max used in BVH traversal. Copied here since the traversal code relies
/// on the specific behaviour wrt NaNs.
///
template <typename T>
static inline const T& _safemax(const T& a, const T& b) {
    return (a > b) ? a : b;
}

///
/// Intersect a ray with a axis-aligned bounding box
///
/// Parameters:
/// - ray_o, ray_d: ray origin and direction
/// - ray_tmin, ray_tmax: ray parameter min, max range
/// - ray_dinv: ray inverse direction
/// - ray_dsign: ray direction sign
/// - bbox_min, bbox_max: bounding box min/max bounds
///
/// Returns:
/// - whether the intersection occurred
///
/// Implementation Notes:
/// - based on "Robust BVH Ray Traversal" by T. Ize published at
/// http://jcgt.org/published/0002/02/02/paper.pdf
///
inline bool intersect_check_bbox(const ray3f& ray, const vec3f& ray_dinv,
    const vec3i& ray_dsign, const bbox3f& bbox) {
    auto txmin = (bbox[ray_dsign[0]][0] - ray.o[0]) * ray_dinv[0];
    auto txmax = (bbox[1 - ray_dsign[0]][0] - ray.o[0]) * ray_dinv[0];
    auto tymin = (bbox[ray_dsign[1]][1] - ray.o[1]) * ray_dinv[1];
    auto tymax = (bbox[1 - ray_dsign[1]][1] - ray.o[1]) * ray_dinv[1];
    auto tzmin = (bbox[ray_dsign[2]][2] - ray.o[2]) * ray_dinv[2];
    auto tzmax = (bbox[1 - ray_dsign[2]][2] - ray.o[2]) * ray_dinv[2];
    auto tmin = _safemax(tzmin, _safemax(tymin, _safemax(txmin, ray.tmin)));
    auto tmax = _safemin(tzmax, _safemin(tymax, _safemin(txmax, ray.tmax)));
    tmax *= 1.00000024f;  // for double: 1.0000000000000004
    return tmin <= tmax;
}

// -----------------------------------------------------------------------------
// POINT-PRIMITIVE DISTANCE FUNCTIONS
// -----------------------------------------------------------------------------

// TODO: documentation
inline bool overlap_point(
    const vec3f& pos, float dist_max, const vec3f& p, float r, float& dist) {
    auto d2 = distsqr(pos, p);
    if (d2 > (dist_max + r) * (dist_max + r)) return false;
    dist = sqrt(d2);
    return true;
}

// TODO: documentation
inline float closestuv_line(
    const vec3f& pos, const vec3f& v0, const vec3f& v1) {
    auto ab = v1 - v0;
    auto d = dot(ab, ab);
    // Project c onto ab, computing parameterized position d(t) = a + t*(b – a)
    auto u = dot(pos - v0, ab) / d;
    u = clamp(u, (float)0, (float)1);
    return u;
}

// TODO: documentation
inline bool overlap_line(const vec3f& pos, float dist_max, const vec3f& v0,
    const vec3f& v1, float r0, float r1, float& dist, vec2f& euv) {
    auto u = closestuv_line(pos, v0, v1);
    // Compute projected position from the clamped t d = a + t * ab;
    auto p = lerp(v0, v1, u);
    auto r = lerp(r0, r1, u);
    auto d2 = distsqr(pos, p);
    // check distance
    if (d2 > (dist_max + r) * (dist_max + r)) return false;
    // done
    dist = sqrt(d2);
    euv = {1 - u, u};
    return true;
}

// TODO: documentation
// this is a complicated test -> I probably prefer to use a sequence of test
// (triangle body, and 3 edges)
inline vec2f closestuv_triangle(
    const vec3f& pos, const vec3f& v0, const vec3f& v1, const vec3f& v2) {
    auto ab = v1 - v0;
    auto ac = v2 - v0;
    auto ap = pos - v0;

    auto d1 = dot(ab, ap);
    auto d2 = dot(ac, ap);

    // corner and edge cases
    if (d1 <= 0 && d2 <= 0) return vec2f{0, 0};

    auto bp = pos - v1;
    auto d3 = dot(ab, bp);
    auto d4 = dot(ac, bp);
    if (d3 >= 0 && d4 <= d3) return vec2f{1, 0};

    auto vc = d1 * d4 - d3 * d2;
    if ((vc <= 0) && (d1 >= 0) && (d3 <= 0)) return vec2f{d1 / (d1 - d3), 0};

    auto cp = pos - v2;
    auto d5 = dot(ab, cp);
    auto d6 = dot(ac, cp);
    if (d6 >= 0 && d5 <= d6) return vec2f{0, 1};

    auto vb = d5 * d2 - d1 * d6;
    if ((vb <= 0) && (d2 >= 0) && (d6 <= 0)) return vec2f{0, d2 / (d2 - d6)};

    auto va = d3 * d6 - d5 * d4;
    if ((va <= 0) && (d4 - d3 >= 0) && (d5 - d6 >= 0)) {
        auto w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return vec2f{1 - w, w};
    }

    // face case
    auto denom = 1 / (va + vb + vc);
    auto v = vb * denom;
    auto w = vc * denom;
    return vec2f{v, w};
}

// TODO: documentation
inline bool overlap_triangle(const vec3f& pos, float dist_max, const vec3f& v0,
    const vec3f& v1, const vec3f& v2, float r0, float r1, float r2, float& dist,
    vec3f& euv) {
    auto uv = closestuv_triangle(pos, v0, v1, v2);
    auto p = blerp(v0, v1, v2, vec3f{1 - uv[0] - uv[1], uv[0], uv[1]});
    auto r = blerp(r0, r1, r2, vec3f{1 - uv[0] - uv[1], uv[0], uv[1]});
    auto dd = distsqr(p, pos);
    if (dd > (dist_max + r) * (dist_max + r)) return false;
    dist = sqrt(dd);
    euv = {1 - uv[0] - uv[1], uv[0], uv[1]};
    return true;
}

// TODO: documentation
inline bool overlap_tetrahedron(const vec3f& pos, const vec3f& v0,
    const vec3f& v1, const vec3f& v2, const vec3f& v3, vec4f& euv) {
    auto vol = dot(v3 - v0, cross(v3 - v1, v3 - v0));
    if (vol == 0) return false;
    auto u = dot(v3 - v0, cross(v3 - v1, v3 - v0)) / vol;
    if (u < 0 || u > 1) return false;
    auto v = dot(v3 - v0, cross(v3 - v1, v3 - v0)) / vol;
    if (v < 0 || v > 1 || u + v > 1) return false;
    auto w = dot(v3 - v0, cross(v3 - v1, v3 - v0)) / vol;
    if (w < 0 || w > 1 || u + v + w > 1) return false;
    euv = {u, v, w, 1 - u - v - w};
    return true;
}

// TODO: documentation
inline bool overlap_tetrahedron(const vec3f& pos, float dist_max,
    const vec3f& v0, const vec3f& v1, const vec3f& v2, const vec3f& v3,
    float r0, float r1, float r2, float r3, float& dist, vec4f& euv) {
    // check interior
    if (overlap_tetrahedron(pos, v0, v1, v2, v3, euv)) {
        dist = 0;
        return true;
    }

    // check faces
    auto hit = false;
    auto tuv = zero3f;
    if (overlap_triangle(pos, dist_max, v0, v1, v2, r0, r1, r2, dist, tuv)) {
        hit = true;
        dist_max = dist;
    }
    if (overlap_triangle(pos, dist_max, v0, v1, v3, r0, r1, r3, dist, tuv)) {
        hit = true;
        dist_max = dist;
    }
    if (overlap_triangle(pos, dist_max, v0, v2, v3, r0, r2, r3, dist, tuv)) {
        hit = true;
        dist_max = dist;
    }
    if (overlap_triangle(pos, dist_max, v1, v2, v3, r1, r2, r3, dist, tuv)) {
        hit = true;
        dist_max = dist;
    }

    return hit;
}

// TODO: documentation
inline bool distance_check_bbox(
    const vec3f& pos, float dist_max, const bbox3f& bbox) {
    // computing distance
    auto dd = 0.0f;

    // For each axis count any excess distance outside box extents
    for (int i = 0; i < 3; i++) {
        auto v = pos[i];
        if (v < bbox[0][i]) dd += (bbox[0][i] - v) * (bbox[0][i] - v);
        if (v > bbox[1][i]) dd += (v - bbox[1][i]) * (v - bbox[1][i]);
    }

    // check distance
    return dd < dist_max * dist_max;
}

// TODO: doc
inline bool overlap_bbox(const bbox3f& bbox1, const bbox3f& bbox2) {
    if (bbox1[1][0] < bbox2[0][0] || bbox1[0][0] > bbox2[1][0]) return false;
    if (bbox1[1][1] < bbox2[0][1] || bbox1[0][1] > bbox2[1][1]) return false;
    if (bbox1[1][2] < bbox2[0][2] || bbox1[0][2] > bbox2[1][2]) return false;
    return true;
}

// TODO: doc
// from "Real-Time Collision Detection" by Christer Ericson, Sect. 4.4.1
inline bool overlap_bbox(const bbox3f& bbox1, const bbox3f& bbox2,
    const frame3f& frame1, const frame3f& frame2) {
    // epsilon
    const auto epsilon = 1e-5f;

    // compute centered frames and extents
    auto cframe1 = frame3f{rot(frame1), transform_point(frame1, center(bbox1))};
    auto cframe2 = frame3f{rot(frame2), transform_point(frame2, center(bbox2))};
    auto ext1 = diagonal(bbox1) / 2.0f, ext2 = diagonal(bbox2) / 2.0f;

    // compute frame from 2 to 1
    auto cframe2to1 = inverse(cframe1) * cframe2;

    // split frame components and move to row-major
    auto rot = transpose(ym::rot(cframe2to1));
    auto t = pos(cframe2to1);

    // Compute common subexpressions. Add in an epsilon term to
    // counteract arithmetic errors when two edges are parallel and
    // their cross product is (near) null (see text for details)
    mat3f absrot;
    auto parallel_axis = false;
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            absrot[i][j] = abs(rot[i][j]) + epsilon;
            if (absrot[i][j] > 1) parallel_axis = true;
        }

    // Test axes L = A0, L = A1, L = A2
    for (int i = 0; i < 3; i++) {
        if (abs(t[i]) > ext1[i] + dot(ext2, absrot[i])) return false;
    }

    // Test axes L = B0, L = B1, L = B2
    for (int i = 0; i < 3; i++) {
        auto ra = ext1[0] * absrot[0][i] + ext1[1] * absrot[1][i] +
                  ext1[2] * absrot[2][i];
        auto rb = ext2[i];
        if (abs(t[0] * rot[0][i] + t[1] * rot[1][i] + t[2] * rot[2][i]) >
            ra + rb)
            return false;
    }

    // if axis were nearly parallel, can exit here
    if (parallel_axis) return true;

    // Test axis L = A0 x B0
    auto ra = ext1[1] * absrot[2][0] + ext1[2] * absrot[1][0];
    auto rb = ext2[1] * absrot[0][2] + ext2[2] * absrot[0][1];
    if (abs(t[2] * rot[1][0] - t[1] * rot[2][0]) > ra + rb) return false;

    // Test axis L = A0 x B1
    ra = ext1[1] * absrot[2][1] + ext1[2] * absrot[1][1];
    rb = ext2[0] * absrot[0][2] + ext2[2] * absrot[0][0];
    if (abs(t[2] * rot[1][1] - t[1] * rot[2][1]) > ra + rb) return false;

    // Test axis L = A0 x B2
    ra = ext1[1] * absrot[2][2] + ext1[2] * absrot[1][2];
    rb = ext2[0] * absrot[0][1] + ext2[1] * absrot[0][0];
    if (abs(t[2] * rot[1][2] - t[1] * rot[2][2]) > ra + rb) return false;

    // Test axis L = A1 x B0
    ra = ext1[0] * absrot[2][0] + ext1[2] * absrot[0][0];
    rb = ext2[1] * absrot[1][2] + ext2[2] * absrot[1][1];
    if (abs(t[0] * rot[2][0] - t[2] * rot[0][0]) > ra + rb) return false;

    // Test axis L = A1 x B1
    ra = ext1[0] * absrot[2][1] + ext1[2] * absrot[0][1];
    rb = ext2[0] * absrot[1][2] + ext2[2] * absrot[1][0];
    if (abs(t[0] * rot[2][1] - t[2] * rot[0][1]) > ra + rb) return false;

    // Test axis L = A1 x B2
    ra = ext1[0] * absrot[2][2] + ext1[2] * absrot[0][2];
    rb = ext2[0] * absrot[1][1] + ext2[1] * absrot[1][0];
    if (abs(t[0] * rot[2][2] - t[2] * rot[0][2]) > ra + rb) return false;

    // Test axis L = A2 x B0
    ra = ext1[0] * absrot[1][0] + ext1[1] * absrot[0][0];
    rb = ext2[1] * absrot[2][2] + ext2[2] * absrot[2][1];
    if (abs(t[1] * rot[0][0] - t[0] * rot[1][0]) > ra + rb) return false;

    // Test axis L = A2 x B1
    ra = ext1[0] * absrot[1][1] + ext1[1] * absrot[0][1];
    rb = ext2[0] * absrot[2][2] + ext2[2] * absrot[2][0];
    if (abs(t[1] * rot[0][1] - t[0] * rot[1][1]) > ra + rb) return false;

    // Test axis L = A2 x B2
    ra = ext1[0] * absrot[1][2] + ext1[1] * absrot[0][2];
    rb = ext2[0] * absrot[2][1] + ext2[1] * absrot[2][0];
    if (abs(t[1] * rot[0][2] - t[0] * rot[1][2]) > ra + rb) return false;

    // Since no separating axis is found, the OBBs must be intersecting
    return true;
}

// this is only a conservative test!
// TODO: rename to something more clear
// TODO: doc
inline bool overlap_bbox_conservative(const bbox3f& bbox1, const bbox3f& bbox2,
    const frame3f& frame1, const frame3f& frame2) {
    return overlap_bbox(bbox1, transform_bbox(inverse(frame1) * frame2, bbox2));
}

// -----------------------------------------------------------------------------
// PRIMITIVE BBOX FUNCTIONS
// -----------------------------------------------------------------------------

///
/// Point bounds
///
inline ym::bbox3f point_bbox(const vec3f& p, float r = 0) {
    return ym::bbox3f{p - vec3f{r, r, r}, p + vec3f{r, r, r}};
}

///
/// Line bounds
///
inline ym::bbox3f line_bbox(
    const vec3f& v0, const vec3f& v1, float r0 = 0, float r1 = 0) {
    return ym::make_bbox({v0 - vec3f{r0, r0, r0}, v0 + vec3f{r0, r0, r0},
        v1 - vec3f{r1, r1, r1}, v1 + vec3f{r1, r1, r1}});
}

///
/// Triangle bounds
///
inline ym::bbox3f triangle_bbox(
    const vec3f& v0, const vec3f& v1, const vec3f& v2) {
    return ym::make_bbox({v0, v1, v2});
}

///
/// Tetrahedron bounds
///
inline ym::bbox3f tetrahedron_bbox(
    const vec3f& v0, const vec3f& v1, const vec3f& v2, const vec3f& v3) {
    return ym::make_bbox({v0, v1, v2, v3});
}

// -----------------------------------------------------------------------------
// UI UTILITIES
// -----------------------------------------------------------------------------

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
        auto z = normalize(*to - *from);
        auto lz = max(T(0.001), dist(*to, *from) * (1 + dolly));
        z *= lz;
        *from = *to - z;
    }

    // pan if necessary
    if (pan[0] || pan[1]) {
        auto z = normalize(*to - *from);
        auto x = normalize(cross(*up, z));
        auto y = normalize(cross(z, x));
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
        auto phi = atan2(frame[2][2], frame[2][0]) + rotate[0];
        auto theta = acos(frame[2][1]) + rotate[1];
        theta = max(T(0.001), min(theta, T(pi - 0.001)));
        auto new_z =
            vec<T, 3>{sin(theta) * cos(phi), cos(theta), sin(theta) * sin(phi)};
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

// -----------------------------------------------------------------------------
// RANDOM NUMBER GENERATION
// -----------------------------------------------------------------------------

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

// -----------------------------------------------------------------------------
// HASHING
// -----------------------------------------------------------------------------

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

// -----------------------------------------------------------------------------
// IMAGE CONTAINERS
// -----------------------------------------------------------------------------

///
/// Image of a specified type
///
template <typename T>
struct image {
    /// empty image constructor
    constexpr image() : _w{0}, _h{0}, _d{} {}
    /// image constructor
    constexpr image(int w, int h, const T& v = {})
        : _w{w}, _h{h}, _d{w * h, v} {}

    /// width
    int width() const { return _w; }
    /// height
    int height() const { return _h; }
    /// size
    vec2i size() const { return {_w, _h}; }

    /// reallocate memory
    void resize(int w, int h, const T& v = {}) {
        _w = w;
        _h = h;
        _d.resize(_w * _h);
    }
    /// reallocate memory
    void assign(int w, int h, const T& v) {
        _w = w;
        _h = h;
        _d.assign(_w * _h, v);
    }

    /// element access
    T& operator[](const vec2i& ij) { return _d[ij.y * _w + ij.x]; }
    /// element access
    const T& operator[](const vec2i& ij) const { return _d[ij.y * _w + ij.x]; }
    /// element access
    T& at(const vec2i& ij) { return _d.at(ij.y * _w + ij.x); }
    /// element access
    const T& at(const vec2i& ij) const { return _d.at(ij.y * _w + ij.x); }

    /// data access
    T* data() { return _d.data(); }
    /// data access
    const T* data() const { return _d.data(); }

   private:
    int _w, _h;
    std::vector<T> _d;
};

// -----------------------------------------------------------------------------
// IMAGE OPERATIONS
// -----------------------------------------------------------------------------

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

/// Approximate conversion from srgb.
inline vec3f srgb_to_linear(const vec3b& srgb) {
    return pow(byte_to_float(srgb), 2.2f);
}

/// Conversion from srgb.
inline vec4f srgb_to_linear(const vec4b& srgb) {
    return {pow(byte_to_float(srgb[0]), 2.2f),
        std::pow(byte_to_float(srgb[1]), 2.2f),
        std::pow(byte_to_float(srgb[2]), 2.2f), byte_to_float(srgb[3])};
}

//
// Tone mapping configurations
//
enum struct tonemap_type { none = 0, srgb, gamma, filmic };

#if 1
///
/// Tone map with a fitted filmic curve.
///
/// Implementation from
/// https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
///
inline vec3f tonemap_filmic(const vec3f& hdr) {
    // rescale
    auto x = hdr * 2.05f;
    // fitted values
    float a = 2.51f, b = 0.03f, c = 2.43f, d = 0.59f, e = 0.14f;
    auto y = ((x * (a * x + b)) / (x * (c * x + d) + e));
    return pow(clamp(y, 0.0f, 1.0f), 1 / 2.2f);
}
#else
inline float tonemap_filmic(float x) {
    auto y =
        (x * (x * (x * (x * 2708.7142 + 6801.1525) + 1079.5474) + 1.1614649) -
            0.00004139375) /
        (x * (x * (x * (x * 983.38937 + 4132.0662) + 2881.6522) + 128.35911) +
            1.0);
    return (float)std::max(y, 0.0);
}
#endif

///
/// Tone mapping HDR to LDR images.
///
inline void tonemap_image(int width, int height, int ncomp, const float* hdr,
    byte* ldr, tonemap_type tm, float exposure, float gamma) {
    if (ncomp < 3 || ncomp > 4)
        throw std::invalid_argument("tonemap supports 3-4 channels only");
    auto scale = pow(2.0f, exposure);
    for (auto j = 0; j < height; j++) {
        for (auto i = 0; i < width; i++) {
            auto h_ptr = hdr + (j * width + i) * ncomp;
            auto l_ptr = ldr + (j * width + i) * ncomp;
            auto h = *(vec3f*)h_ptr;
            h *= scale;
            switch (tm) {
                case tonemap_type::none: break;
                case tonemap_type::srgb: h = pow(h, 1 / 2.2f); break;
                case tonemap_type::gamma: h = pow(h, 1 / gamma); break;
                case tonemap_type::filmic: h = tonemap_filmic(h); break;
            }
            *(vec3b*)l_ptr = float_to_byte(h);
            if (ncomp == 4) l_ptr[3] = float_to_byte(h_ptr[3]);
        }
    }
}

}  // namespace ym

// HACK to avoid compilation with MSVC2015 without dirtying code
#ifdef constexpr
#undef constexpr
#endif

#endif
