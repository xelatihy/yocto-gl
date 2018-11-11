//
// # Yocto/Math: Tiny library for math support in graphics applications.
//
// Yocto/Math provides the basic math primitives used in grahics, including
// small-sized vectors and matrixes, frames, bounding boxes and transforms.
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
// `translation_mat()` or `make_translation_frame()` respectively, etc.
// For rotation we support axis-angle and quaternions, with slerp.
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

#ifndef _YOCTO_MATH_H_
#define _YOCTO_MATH_H_

#ifndef YOCTO_SPECIALIZATIONS
#define YOCTO_SPECIALIZATIONS 1
#endif

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <initializer_list>
#include <limits>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
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
using std::initializer_list;
using std::pair;

using byte = unsigned char;
using uint = unsigned int;

}  // namespace yocto

// -----------------------------------------------------------------------------
// MATH CONSTANTS AND FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

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
constexpr const T type_min = std::numeric_limits<T>::lowest();
template <class T>
constexpr const T type_max = std::numeric_limits<T>::max();
template <class T>
constexpr const T    type_epsilon  = std::numeric_limits<T>::epsilon();
constexpr const auto float_max     = type_max<float>;
constexpr const auto float_epsilon = type_epsilon<float>;

template <class T>
constexpr const T    pi  = (T)3.14159265358979323846;
constexpr const auto pif = 3.14159265f;

}  // namespace yocto

// -----------------------------------------------------------------------------
// VECTORS
// -----------------------------------------------------------------------------
namespace yocto {

// Small-sized vectors
template <typename T, int N>
struct vec;

template <typename T>
struct vec<T, 1> {
    T elements[1] = {0};

    constexpr vec() : elements{0} {}
    constexpr vec(const T& x) : elements{x} {}

    constexpr T&       operator[](int idx) { return elements[idx]; }
    constexpr const T& operator[](int idx) const { return elements[idx]; }
};
template <typename T>
struct vec<T, 2> {
    T elements[2] = {0};

    constexpr vec() : elements{0, 0} {}
    constexpr vec(const T& x, const T& y) : elements{x, y} {}

    constexpr T&       operator[](int idx) { return elements[idx]; }
    constexpr const T& operator[](int idx) const { return elements[idx]; }
};
template <typename T>
struct vec<T, 3> {
    T elements[3] = {0};

    constexpr vec() : elements{0, 0, 0} {}
    constexpr vec(const T& x, const T& y, const T& z) : elements{x, y, z} {}

    constexpr T&       operator[](int idx) { return elements[idx]; }
    constexpr const T& operator[](int idx) const { return elements[idx]; }
};
template <typename T>
struct vec<T, 4> {
    T elements[4] = {0};

    constexpr vec() : elements{0, 0, 0, 0} {}
    constexpr vec(const T& x, const T& y, const T& z, const T& w)
        : elements{x, y, z, w} {}

    constexpr T&       operator[](int idx) { return elements[idx]; }
    constexpr const T& operator[](int idx) const { return elements[idx]; }
};

// Vector constants
template <typename T, int N>
constexpr inline vec<T, N> make_zero_vec() {
    auto v = vec<T, N>{};
    for (auto i = 0; i < N; i++) v[i] = 0;
    return v;
}
template <typename T, int N>
constexpr inline vec<T, N> make_one_vec() {
    auto v = vec<T, N>{};
    for (auto i = 0; i < N; i++) v[i] = 1;
    return v;
}
template <typename T, int N>
constexpr inline vec<T, N> make_uniform_vec(T value) {
    auto v = vec<T, N>{};
    for (auto i = 0; i < N; i++) v[i] = value;
    return v;
}
template <typename T, int N>
constexpr const vec<T, N> zero_vec = make_zero_vec<T, N>();
template <typename T, int N>
constexpr const vec<T, N> one_vec = make_one_vec<T, N>();

// Type aliases.
using vec1f = vec<float, 1>;
using vec2f = vec<float, 2>;
using vec3f = vec<float, 3>;
using vec4f = vec<float, 4>;
using vec1i = vec<int, 1>;
using vec2i = vec<int, 2>;
using vec3i = vec<int, 3>;
using vec4i = vec<int, 4>;
using vec4b = vec<byte, 4>;

// Zero vector constants.
constexpr const auto zero_vec1f = zero_vec<float, 1>;
constexpr const auto zero_vec2f = zero_vec<float, 2>;
constexpr const auto zero_vec3f = zero_vec<float, 3>;
constexpr const auto zero_vec4f = zero_vec<float, 4>;
constexpr const auto zero_vec1i = zero_vec<int, 1>;
constexpr const auto zero_vec2i = zero_vec<int, 2>;
constexpr const auto zero_vec3i = zero_vec<int, 3>;
constexpr const auto zero_vec4i = zero_vec<int, 4>;
constexpr const auto zero_vec4b = zero_vec<byte, 4>;

// Element access
template <typename T, int N>
constexpr inline vec<T, N - 1> make_shorter_vec(const vec<T, N>& value) {
    auto smaller = vec<T, N - 1>{};
    for (auto i = 0; i < N - 1; i++) smaller[i] = value[i];
    return smaller;
}
template <typename T, int N>
constexpr inline vec<T, N + 1> make_longer_vec(const vec<T, N>& value, T last) {
    auto longer = vec<T, N + 1>{};
    for (auto i = 0; i < N; i++) longer[i] = value[i];
    longer[N] = last;
    return longer;
}

// Vector comparison operations.
template <typename T, int N>
constexpr inline bool operator==(const vec<T, N>& a, const vec<T, N>& b) {
    for (auto i = 0; i < N; i++)
        if (a[i] != b[i]) return false;
    return true;
}
template <typename T, int N>
constexpr inline bool operator!=(const vec<T, N>& a, const vec<T, N>& b) {
    for (auto i = 0; i < N; i++)
        if (a[i] != b[i]) return true;
    return false;
}
template <typename T, int N, typename T1>
constexpr inline bool operator==(const vec<T, N>& a, T1 b) {
    for (auto i = 0; i < N; i++)
        if (a[i] != b) return false;
    return true;
}
template <typename T, int N, typename T1>
constexpr inline bool operator!=(const vec<T, N>& a, T1 b) {
    for (auto i = 0; i < N; i++)
        if (a[i] == b) return false;
    return true;
}

// Vector operations.
template <typename T, int N>
constexpr inline vec<T, N> operator+(const vec<T, N>& a) {
    return a;
}
template <typename T, int N>
constexpr inline vec<T, N> operator-(const vec<T, N>& a) {
    auto c = vec<T, N>{};
    for (auto i = 0; i < N; i++) c[i] = -a[i];
    return c;
}
template <typename T, int N>
constexpr inline vec<T, N> operator+(const vec<T, N>& a, const vec<T, N>& b) {
    auto c = vec<T, N>{};
    for (auto i = 0; i < N; i++) c[i] = a[i] + b[i];
    return c;
}
template <typename T, int N, typename T1>
constexpr inline vec<T, N> operator+(const vec<T, N>& a, T1 b) {
    auto c = vec<T, N>{};
    for (auto i = 0; i < N; i++) c[i] = a[i] + b;
    return c;
}
template <typename T, int N, typename T1>
constexpr inline vec<T, N> operator+(T1 a, const vec<T, N>& b) {
    auto c = vec<T, N>{};
    for (auto i = 0; i < N; i++) c[i] = a + b[i];
    return c;
}
template <typename T, int N>
constexpr inline vec<T, N> operator-(const vec<T, N>& a, const vec<T, N>& b) {
    auto c = vec<T, N>{};
    for (auto i = 0; i < N; i++) c[i] = a[i] - b[i];
    return c;
}
template <typename T, int N, typename T1>
constexpr inline vec<T, N> operator-(const vec<T, N>& a, T1 b) {
    auto c = vec<T, N>{};
    for (auto i = 0; i < N; i++) c[i] = a[i] - b;
    return c;
}
template <typename T, int N, typename T1>
constexpr inline vec<T, N> operator-(T1 a, const vec<T, N>& b) {
    auto c = vec<T, N>{};
    for (auto i = 0; i < N; i++) c[i] = a - b[i];
    return c;
}
template <typename T, int N>
constexpr inline vec<T, N> operator*(const vec<T, N>& a, const vec<T, N>& b) {
    auto c = vec<T, N>{};
    for (auto i = 0; i < N; i++) c[i] = a[i] * b[i];
    return c;
}
template <typename T, int N, typename T1>
constexpr inline vec<T, N> operator*(const vec<T, N>& a, T1 b) {
    auto c = vec<T, N>{};
    for (auto i = 0; i < N; i++) c[i] = a[i] * b;
    return c;
}
template <typename T, int N, typename T1>
constexpr inline vec<T, N> operator*(T1 a, const vec<T, N>& b) {
    auto c = vec<T, N>{};
    for (auto i = 0; i < N; i++) c[i] = a * b[i];
    return c;
}
template <typename T, int N>
constexpr inline vec<T, N> operator/(const vec<T, N>& a, const vec<T, N>& b) {
    auto c = vec<T, N>{};
    for (auto i = 0; i < N; i++) c[i] = a[i] / b[i];
    return c;
}
template <typename T, int N, typename T1>
constexpr inline vec<T, N> operator/(const vec<T, N>& a, T1 b) {
    auto c = vec<T, N>{};
    for (auto i = 0; i < N; i++) c[i] = a[i] / b;
    return c;
}
template <typename T, int N, typename T1>
constexpr inline vec<T, N> operator/(T1 a, const vec<T, N>& b) {
    auto c = vec<T, N>{};
    for (auto i = 0; i < N; i++) c[i] = a / b[i];
    return c;
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
template <typename T, int N>
constexpr inline T dot(const vec<T, N>& a, const vec<T, N>& b) {
    auto c = T{0};
    for (auto i = 0; i < N; i++) c += a[i] * b[i];
    return c;
}
template <typename T>
constexpr inline T cross(const vec<T, 2>& a, const vec<T, 2>& b) {
    return a[0] * b[1] - a[1] * b[0];
}
template <typename T>
constexpr inline vec<T, 3> cross(const vec<T, 3>& a, const vec<T, 3>& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0]};
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
    return fabs(v[0]) > fabs(v[2]) ? vec<T, 3>{-v[1], v[0], 0} :
                                     vec<T, 3>{0, -v[2], v[1]};
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
template <typename T, int N, typename T1, typename T2>
constexpr inline vec<T, N> clamp(const vec<T, N>& value, T1 min, T2 max) {
    auto clamped = T{};
    for (auto i = 0; i < N; i++) clamped[i] = clamp(value[i], min, max);
    return clamped;
}
template <typename T, int N>
constexpr inline T max(const vec<T, N>& a) {
    auto c = T{a[0]};
    for (auto i = 1; i < N; i++) c = max(c, a[i]);
    return c;
}
template <typename T, int N>
constexpr inline T min(const vec<T, N>& a) {
    auto c = T{a[0]};
    for (auto i = 1; i < N; i++) c = min(c, a[i]);
    return c;
}
template <typename T, int N>
constexpr inline T mean(const vec<T, N>& a) {
    auto c = T{0};
    for (auto i = 0; i < N; i++) c += a[i];
    return c / N;
}

// Quaternion operatons represented as xi + yj + zk + w
const auto identity_quat4f = vec4f{0, 0, 0, 1};
template <typename T>
constexpr inline vec<T, 4> quat_mul(const vec<T, 4>& a, float b) {
    return {a[0] * b, a[1] * b, a[2] * b, a[3] * b};
}
template <typename T>
constexpr inline vec<T, 4> quat_mul(const vec<T, 4>& a, const vec<T, 4>& b) {
    return {a[0] * b[3] + a[3] * b[0] + a[1] * b[3] - a[2] * b[1],
        a[1] * b[3] + a[3] * b[1] + a[2] * b[0] - a[0] * b[2],
        a[2] * b[3] + a[3] * b[2] + a[0] * b[1] - a[1] * b[0],
        a[3] * b[3] - a[0] * b[0] - a[1] * b[1] - a[2] * b[2]};
}
template <typename T>
constexpr inline vec<T, 4> quat_conjugate(const vec<T, 4>& a) {
    return {-a[0], -a[1], -a[2], a[3]};
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

// Small fixed-size matrices stored in column major format.
template <typename T, int N, int M>
struct mat;

template <typename T, int N>
struct mat<T, N, 1> {
    vec<T, N> columns[1] = {};

    constexpr mat() : columns{} {}
    constexpr mat(const vec<T, N>& c0) : columns{c0} {}

    constexpr vec<T, N>&       operator[](int idx) { return columns[idx]; }
    constexpr const vec<T, N>& operator[](int idx) const {
        return columns[idx];
    }
};
template <typename T, int N>
struct mat<T, N, 2> {
    vec<T, N> columns[2] = {};

    constexpr mat() : columns{} {}
    constexpr mat(const vec<T, N>& c0, const vec<T, N>& c1) : columns{c0, c1} {}

    constexpr vec<T, N>&       operator[](int idx) { return columns[idx]; }
    constexpr const vec<T, N>& operator[](int idx) const {
        return columns[idx];
    }
};
template <typename T, int N>
struct mat<T, N, 3> {
    vec<T, N> columns[3] = {};

    constexpr mat() : columns{} {}
    constexpr mat(const vec<T, N>& c0, const vec<T, N>& c1, const vec<T, N>& c2)
        : columns{c0, c1, c2} {}

    constexpr vec<T, N>&       operator[](int idx) { return columns[idx]; }
    constexpr const vec<T, N>& operator[](int idx) const {
        return columns[idx];
    }
};
template <typename T, int N>
struct mat<T, N, 4> {
    vec<T, N> columns[4] = {};

    constexpr mat() : columns{} {}
    constexpr mat(const vec<T, N>& c0, const vec<T, N>& c1, const vec<T, N>& c2,
        const vec<T, N>& c3)
        : columns{c0, c1, c2, c3} {}

    constexpr vec<T, N>&       operator[](int idx) { return columns[idx]; }
    constexpr const vec<T, N>& operator[](int idx) const {
        return columns[idx];
    }
};

// Matrix contants
template <typename T, int N>
constexpr inline mat<T, N, N> make_identity_mat() {
    auto c = mat<T, N, N>{};
    for (auto j = 0; j < N; j++)
        for (auto i = 0; i < N; i++) c[j][i] = j == i ? 1 : 0;
    return c;
}
template <typename T, int N>
constexpr inline mat<T, N, N> make_diagonal_mat(const vec<T, N>& diagonal) {
    auto c = mat<T, N, N>{};
    for (auto j = 0; j < N; j++)
        for (auto i = 0; i < N; i++) c[j][i] = j == i ? diagonal[i] : 0;
    return c;
}
template <typename T, int N>
constexpr const mat<T, N, N> identity_mat = make_identity_mat<T, N>();

// Type aliases.
using mat1f = mat<float, 1, 1>;
using mat2f = mat<float, 2, 2>;
using mat3f = mat<float, 3, 3>;
using mat4f = mat<float, 4, 4>;

// Identity matrices constants.
constexpr const auto identity_mat1f = identity_mat<float, 1>;
constexpr const auto identity_mat2f = identity_mat<float, 2>;
constexpr const auto identity_mat3f = identity_mat<float, 3>;
constexpr const auto identity_mat4f = identity_mat<float, 4>;

// Matrix comparisons.
template <typename T, int N, int M>
constexpr inline bool operator==(const mat<T, N, M>& a, const mat<T, N, M>& b) {
    for (auto j = 0; j < M; j++)
        if (a[j] != b[j]) return false;
    return true;
}
template <typename T, int N, int M>
constexpr inline bool operator!=(const mat<T, N, 1>& a, const mat<T, N, 1>& b) {
    for (auto j = 0; j < M; j++)
        if (a[j] != b[j]) return true;
    return false;
}

// Matrix operations.
template <typename T, int N, int M>
constexpr inline mat<T, N, M> operator+(
    const mat<T, N, M>& a, const mat<T, N, M>& b) {
    auto c = mat<T, N, M>{};
    for (auto j = 0; j < M; j++) c[j] = a[j] + b[j];
    return c;
}
template <typename T, int N, int M>
constexpr inline mat<T, N, M> operator*(const mat<T, N, M>& a, T b) {
    auto c = mat<T, N, M>{};
    for (auto j = 0; j < M; j++) c[j] = a[j] * b;
    return c;
}
template <typename T, int N, int M>
constexpr inline vec<T, N> operator*(const mat<T, N, M>& a, const vec<T, M>& b) {
    auto c = vec<T, N>{};
    for (auto j = 0; j < M; j++) c += a[j] * b[j];
    return c;
}
template <typename T, int N, int M>
constexpr inline vec<T, M> operator*(const vec<T, N>& a, const mat<T, N, M>& b) {
    auto c = vec<T, M>{};
    for (auto j = 0; j < M; j++) c[j] = dot(a, b[j]);
    return c;
}
template <typename T, int N, int M, int K>
constexpr inline mat<T, N, M> operator*(
    const mat<T, N, K>& a, const mat<T, K, M>& b) {
    auto c = mat<T, N, M>{};
    for (auto j = 0; j < M; j++) c[j] = a * b[j];
    return c;
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
template <typename T, int N>
constexpr inline vec<T, N> diagonal(const mat<T, N, N>& a) {
    auto c = vec<T, N>{};
    for (auto j = 0; j < N; j++) c[j] = a[j][j];
    return c;
}
template <typename T, int N, int M>
constexpr inline mat<T, M, N> transpose(const mat<T, N, M>& a) {
    auto c = mat<T, M, N>{};
    for (auto j = 0; j < M; j++)
        for (auto i = 0; i < N; i++) c[i][j] = a[j][i];
    return c;
}

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
struct frame {
    mat<T, N, N> axes   = {};
    vec<T, N>    origin = {};

    constexpr frame() : axes{identity_mat<T, N>}, origin{zero_vec<T, N>} {}
    constexpr frame(const mat<T, N, N>& axes_, const vec<T, N>& origin_)
        : axes{axes_}, origin{origin_} {}

    constexpr vec<T, N>& operator[](int idx) {
        return idx < N ? axes[idx] : origin;
    }
    constexpr const vec<T, N>& operator[](int idx) const {
        return idx < N ? axes[idx] : origin;
    }
};

// Frame contants
template <typename T, int N>
constexpr inline frame<T, N> make_identity_frame() {
    return {make_identity_mat<T, N>(), make_zero_vec<T, N>()};
}
template <typename T, int N>
constexpr const frame<T, N> identity_frame = make_identity_frame<T, N>();

// Type aliases.
using frame2f = frame<float, 2>;
using frame3f = frame<float, 3>;

// Indentity frames.
constexpr const auto identity_frame2f = identity_frame<float, 2>;
constexpr const auto identity_frame3f = identity_frame<float, 3>;

// Frame construction from axis.
template <typename T>
constexpr inline frame<T, 3> make_frame_fromz(
    const vec<T, 3>& o, const vec<T, 3>& v) {
    auto z = normalize(v);
    auto x = normalize(orthogonal(z));
    auto y = normalize(cross(z, x));
    return {{x, y, z}, o};
}
template <typename T>
constexpr inline frame<T, 3> make_frame_fromzx(
    const vec<T, 3>& o, const vec<T, 3>& z_, const vec<T, 3>& x_) {
    auto z = normalize(z_);
    auto x = orthonormalize(x_, z);
    auto y = normalize(cross(z, x));
    return {{x, y, z}, o};
}

// Frame to matrix conversion.
template <typename T, int N>
constexpr inline mat<T, N + 1, N + 1> frame_to_mat(const frame<T, N>& a) {
    auto c = mat<T, N + 1, N + 1>{};
    for (auto j = 0; j < N + 1; j++)
        c[j] = make_longer_vec(a[j], j == N ? (T)1 : (T)0);
    return c;
}
template <typename T, int N>
constexpr inline frame<T, N - 1> mat_to_frame(const mat<T, N, N>& a) {
    auto c = frame<T, N - 1>{};
    for (auto j = 0; j < N; j++) c[j] = make_shorter_vec(a[j]);
    return c;
}

// Frame comparisons.
template <typename T, int N>
constexpr inline bool operator==(const frame<T, N>& a, const frame<T, N>& b) {
    return a.axes == b.axes && a.origin == b.origin;
}
template <typename T, int N>
constexpr inline bool operator!=(const frame<T, N>& a, const frame<T, N>& b) {
    return a.axes != b.axes || a.origin != b.origin;
}

// Frame composition, equivalent to affine matrix product.
template <typename T, int N>
constexpr inline frame<T, N> operator*(
    const frame<T, N>& a, const frame<T, N>& b) {
    return {a.axes * b.axes, a.axes * b.origin + a.origin};
}
// Frame inverse, equivalent to rigid affine inverse.
template <typename T, int N>
constexpr inline frame<T, N> inverse(const frame<T, N>& a, bool is_rigid = true) {
    auto axes_inv = (is_rigid) ? transpose(a.axes) : inverse(a.axes);
    return {axes_inv, -(axes_inv * a.origin)};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// AXIS ALIGNED BOUNDING BOXES
// -----------------------------------------------------------------------------
namespace yocto {

// Describes a range of values in N dimensions.
template <typename T, int N>
struct bbox {
    using V = vec<T, N>;

    vec<T, N> min = make_uniform_vec<T, N>(type_max<T>);
    vec<T, N> max = make_uniform_vec<T, N>(type_min<T>);

    constexpr bbox()
        : min{make_uniform_vec<T, N>(type_max<T>)}
        , max{make_uniform_vec<T, N>(type_min<T>)} {}
    constexpr bbox(const vec<T, N>& min_, const vec<T, N>& max_)
        : min{min_}, max{max_} {}
    constexpr bbox(const bbox&) = default;

    constexpr vec<T, N>& operator[](int idx) { return idx == 0 ? min : max; }
    constexpr const vec<T, N>& operator[](int idx) const {
        return idx == 0 ? min : max;
    }
};

// Bbox constants
template <typename T, int N>
constexpr inline bbox<T, N> make_invalid_bbox() {
    return {make_uniform_vec<T, N>(type_max<T>),
        make_uniform_vec<T, N>(type_min<T>)};
}
template <typename T, int N>
constexpr const bbox<T, N> invalid_bbox = make_invalid_bbox<T, N>();

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
constexpr const auto invalid_bbox1f = invalid_bbox<float, 1>;
constexpr const auto invalid_bbox2f = invalid_bbox<float, 2>;
constexpr const auto invalid_bbox3f = invalid_bbox<float, 3>;
constexpr const auto invalid_bbox4f = invalid_bbox<float, 4>;
constexpr const auto invalid_bbox1i = invalid_bbox<int, 1>;
constexpr const auto invalid_bbox2i = invalid_bbox<int, 2>;
constexpr const auto invalid_bbox3i = invalid_bbox<int, 3>;
constexpr const auto invalid_bbox4i = invalid_bbox<int, 4>;

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
template <typename T, int N>
constexpr inline bool operator==(const bbox<T, N>& a, const bbox<T, N>& b) {
    return a.min == b.min && a.max == b.max;
}
template <typename T, int N>
constexpr inline bool operator!=(const bbox<T, N>& a, const bbox<T, N>& b) {
    return a.min != b.min || a.max != b.max;
}

// Bounding box expansions with points and other boxes.
template <typename T>
constexpr inline bbox<T, 1>& operator+=(bbox<T, 1>& a, T b) {
    a.min[0] = {min(a.min[0], b)};
    a.max[0] = {max(a.max[0], b)};
    return a;
}
// Bounding box expansions with points and other boxes.
template <typename T, int N>
constexpr inline bbox<T, N>& operator+=(bbox<T, N>& a, const vec<T, N>& b) {
    for (auto i = 0; i < N; i++) {
        a.min[i] = min(a.min[i], b[i]);
        a.max[i] = max(a.max[i], b[i]);
    }
    return a;
}
template <typename T, int N>
constexpr inline bbox<T, N>& operator+=(bbox<T, N>& a, const bbox<T, N>& b) {
    for (auto i = 0; i < N; i++) {
        a.min[i] = min(a.min[i], b.min[i]);
        a.max[i] = max(a.max[i], b.max[i]);
    }
    return a;
}

// Primitive bounds.
template <typename T>
constexpr inline bbox<T, 3> point_bounds(const vec<T, 3>& p, T r = 0) {
    auto bounds = bbox<T, 3>{};
    bounds += p - vec<T, 3>{r, r, r};
    bounds += p + vec<T, 3>{r, r, r};
    return bounds;
}
template <typename T>
constexpr inline bbox<T, 3> line_bounds(
    const vec<T, 3>& p0, const vec<T, 3>& p1, T r0 = 0, T r1 = 0) {
    auto bounds = bbox<T, 3>{};
    bounds += p0 - vec<T, 3>{r0, r0, r0};
    bounds += p0 + vec<T, 3>{r0, r0, r0};
    bounds += p1 - vec<T, 3>{r1, r1, r1};
    bounds += p1 + vec<T, 3>{r1, r1, r1};
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
struct ray {
    vec<T, N> origin    = {0, 0};
    vec<T, N> direction = {0, 1};
    T         tmin      = 0;
    T         tmax      = type_max<T>;

    constexpr ray() : origin{}, direction{}, tmin{0}, tmax{max<T>()} {}
    constexpr ray(const vec<T, N>& origin_, const vec<T, N>& direction_,
        const T& tmin_, const T& tmax_)
        : origin{origin_}, direction{direction_}, tmin{tmin_}, tmax{tmax_} {}
    constexpr ray(const ray&) = default;
};

// Type aliases.
using ray2f = ray<float, 2>;
using ray3f = ray<float, 3>;

// Default ray epsilon
const auto default_ray_eps = 1e-4;

// Construct a ray from direction or segments using a default epsilon.
template <typename T, int N>
constexpr inline ray<T, N> make_ray(const vec<T, N>& origin,
    const vec<T, N>& direction, T eps = (T)default_ray_eps) {
    return {origin, direction, eps, type_max<T>};
}
template <typename T, int N>
constexpr inline ray<T, N> make_segment(
    const vec<T, N>& p1, const vec<T, N>& p2, T eps = (T)default_ray_eps) {
    return {p1, normalize(p2 - p1), eps, length(p2 - p1) - 2 * eps};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// TRANSFORMS
// -----------------------------------------------------------------------------
namespace yocto {

// Transforms points, vectors and directions by matrices.
template <typename T, int N>
constexpr inline vec<T, N> transform_point(
    const mat<T, N + 1, N + 1>& a, const vec<T, N>& b) {
    auto tvb = a * make_longer_vec(b, (T)1);
    return make_shorter_vec(tvb) / tvb[N];
}
template <typename T, int N>
constexpr inline vec<T, N> transform_vector(
    const mat<T, N + 1, N + 1>& a, const vec<T, N>& b) {
    auto tvb = a * make_longer_vec(b, (T)0);
    return make_shorter_vec(tvb) / tvb[N];
}
template <typename T, int N>
constexpr inline vec<T, N> transform_vector(
    const mat<T, N, N>& a, const vec<T, N>& b) {
    return a * b;
}
template <typename T, int N>
constexpr inline vec<T, N> transform_direction(
    const mat<T, N + 1, N + 1>& a, const vec<T, N>& b) {
    return normalize(transform_vector(a, b));
}

// Transforms points, vectors and directions by frames.
template <typename T, int N>
constexpr inline vec<T, N> transform_point(
    const frame<T, N>& a, const vec<T, N>& b) {
    return a.axes * b + a.origin;
}
template <typename T, int N>
constexpr inline vec<T, N> transform_vector(
    const frame<T, N>& a, const vec<T, N>& b) {
    return a.axes * b;
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
    return {transform_point(a, b.origin), transform_vector(a, b.direction),
        b.tmin, b.tmax};
}
template <typename T, int N>
constexpr inline ray<T, N> transform_ray(
    const mat<T, N + 1, N + 1>& a, const ray<T, N>& b) {
    return {transform_point(a, b.origin), transform_vector(a, b.direction),
        b.tmin, b.tmax};
}

template <typename T>
constexpr inline bbox<T, 2> transform_bbox(
    const frame<T, 2>& a, const bbox<T, 2>& b) {
    auto corners = {vec<T, 2>{b.min[0], b.min[1]}, vec<T, 2>{b.min[0], b.max[1]},
        vec<T, 2>{b.max[0], b.min[1]}, vec<T, 2>{b.max[0], b.max[1]}};
    auto xformed = bbox<T, 3>();
    for (auto& corner : corners) xformed += transform_point(a, corner);
    return xformed;
}
template <typename T>
constexpr inline bbox<T, 2> transform_bbox(
    const mat<T, 3, 3>& a, const bbox<T, 2>& b) {
    auto corners = {vec<T, 2>{b.min[0], b.min[1]}, vec<T, 2>{b.min[0], b.max[1]},
        vec<T, 2>{b.max[0], b.min[1]}, vec<T, 2>{b.max[0], b.max[1]}};
    auto xformed = bbox<T, 3>();
    for (auto& corner : corners) xformed += transform_point(a, corner);
    return xformed;
}

template <typename T>
constexpr inline bbox<T, 3> transform_bbox(
    const frame<T, 3>& a, const bbox<T, 3>& b) {
    auto corners = {vec<T, 3>{b.min[0], b.min[1], b.min[2]},
        vec<T, 3>{b.min[0], b.min[1], b.max[2]},
        vec<T, 3>{b.min[0], b.max[1], b.min[2]},
        vec<T, 3>{b.min[0], b.max[1], b.max[2]},
        vec<T, 3>{b.max[0], b.min[1], b.min[2]},
        vec<T, 3>{b.max[0], b.min[1], b.max[2]},
        vec<T, 3>{b.max[0], b.max[1], b.min[2]},
        vec<T, 3>{b.max[0], b.max[1], b.max[2]}};
    auto xformed = bbox<T, 3>();
    for (auto& corner : corners) xformed += transform_point(a, corner);
    return xformed;
}
template <typename T>
constexpr inline bbox<T, 3> transform_bbox(
    const mat<T, 4, 4>& a, const bbox<T, 3>& b) {
    auto corners = {vec<T, 3>{b.min[0], b.min[1], b.min[2]},
        vec<T, 3>{b.min[0], b.min[1], b.max[2]},
        vec<T, 3>{b.min[0], b.max[1], b.min[2]},
        vec<T, 3>{b.min[0], b.max[1], b.max[2]},
        vec<T, 3>{b.max[0], b.min[1], b.min[2]},
        vec<T, 3>{b.max[0], b.min[1], b.max[2]},
        vec<T, 3>{b.max[0], b.max[1], b.min[2]},
        vec<T, 3>{b.max[0], b.max[1], b.max[2]}};
    auto xformed = bbox<T, 3>();
    for (auto& corner : corners) xformed += transform_point(a, corner);
    return xformed;
}

// Inverse transforms by frames, assuming they are rigid transforms.
template <typename T, int N>
constexpr inline vec<T, N> transform_point_inverse(
    const frame<T, N>& a, const vec<T, N>& b) {
    return (b - a.origin) * a.axes;
}
template <typename T, int N>
constexpr inline vec<T, N> transform_vector_inverse(
    const frame<T, N>& a, const vec<T, N>& b) {
    return b * a.axes;
}
template <typename T, int N>
constexpr inline vec<T, N> transform_direction_inverse(
    const frame<T, N>& a, const vec<T, N>& b) {
    return normalize(transform_vector_inverse(a, b));
}
template <typename T, int N>
constexpr inline ray<T, N> transform_ray_inverse(
    const frame<T, N>& a, const ray<T, N>& b) {
    return {transform_point_inverse(a, b.origin),
        transform_direction_inverse(a, b.direction), b.tmin, b.tmax};
}
template <typename T, int N>
constexpr inline bbox<T, N> transform_bbox_inverse(
    const frame<T, N>& a, const bbox<T, N>& b) {
    return transform_bbox(inverse(a), b);
}

// Translation, scaling and rotations transforms.
template <typename T, int N>
constexpr inline frame<T, N> make_translation_frame(const vec<T, N>& a) {
    return {identity_mat<T, N>, a};
}
template <typename T, int N>
constexpr inline frame<T, N> make_scaling_frame(const vec<T, N>& a) {
    return {make_diagonal_mat<T, N>(a), {}};
}
template <typename T>
constexpr inline frame<T, 3> make_rotation_frame(const vec<T, 3>& axis, T angle) {
    auto s = sin(angle), c = cos(angle);
    auto vv = normalize(axis);
    return {{{c + (1 - c) * vv[0] * vv[0], (1 - c) * vv[0] * vv[1] + s * vv[2],
                 (1 - c) * vv[0] * vv[2] - s * vv[1]},
                {(1 - c) * vv[0] * vv[1] - s * vv[2], c + (1 - c) * vv[1] * vv[1],
                    (1 - c) * vv[1] * vv[2] + s * vv[0]},
                {(1 - c) * vv[0] * vv[2] + s * vv[1],
                    (1 - c) * vv[1] * vv[2] - s * vv[0],
                    c + (1 - c) * vv[2] * vv[2]}},
        {0, 0, 0}};
}
template <typename T>
constexpr inline frame<T, 3> make_rotation_frame(const vec<T, 4>& quat) {
    auto v = quat;
    return {
        {{v[3] * v[3] + v[0] * v[0] - v[1] * v[1] - v[2] * v[2],
             (v[0] * v[1] + v[2] * v[3]) * 2, (v[2] * v[0] - v[1] * v[3]) * 2},
            {(v[0] * v[1] - v[2] * v[3]) * 2,
                v[3] * v[3] - v[0] * v[0] + v[1] * v[1] - v[2] * v[2],
                (v[1] * v[2] + v[0] * v[3]) * 2},
            {(v[2] * v[0] + v[1] * v[3]) * 2, (v[1] * v[2] - v[0] * v[3]) * 2,
                v[3] * v[3] - v[0] * v[0] - v[1] * v[1] + v[2] * v[2]}},
        {0, 0, 0}};
}
template <typename T>
constexpr inline frame<T, 3> make_rotation_frame(const mat<T, 3, 3>& rot) {
    return {rot, rot, rot, {}};
}

// Lookat frame. Z-axis can be inverted with inv_xz.
template <typename T>
constexpr inline frame<T, 3> make_lookat_frame(const vec<T, 3>& eye,
    const vec<T, 3>& center, const vec<T, 3>& up, bool inv_xz = false) {
    auto w = normalize(eye - center);
    auto u = normalize(cross(up, w));
    auto v = normalize(cross(w, u));
    if (inv_xz) {
        w = -w;
        u = -u;
    }
    return {{u, v, w}, eye};
}

// OpenGL frustum, ortho and perspecgive matrices.
template <typename T>
constexpr inline mat<T, 4, 4> make_frustum_mat(T l, T r, T b, T t, T n, T f) {
    return {{2 * n / (r - l), 0, 0, 0}, {0, 2 * n / (t - b), 0, 0},
        {(r + l) / (r - l), (t + b) / (t - b), -(f + n) / (f - n), -1},
        {0, 0, -2 * f * n / (f - n), 0}};
}
template <typename T>
constexpr inline mat<T, 4, 4> make_orthographic_mat(T l, T r, T b, T t, T n, T f) {
    return {{2 / (r - l), 0, 0, 0}, {0, 2 / (t - b), 0, 0},
        {0, 0, -2 / (f - n), 0},
        {-(r + l) / (r - l), -(t + b) / (t - b), -(f + n) / (f - n), 1}};
}
template <typename T>
constexpr inline mat<T, 4, 4> make_orthographic2d_mat(
    T left, T right, T bottom, T top) {
    return make_orthographic_mat(left, right, bottom, top, -1, 1);
}
template <typename T>
constexpr inline mat<T, 4, 4> make_orthographic_mat(
    T xmag, T ymag, T near, T far) {
    return {{1 / xmag, 0, 0, 0}, {0, 1 / ymag, 0, 0},
        {0, 0, 2 / (near - far), 0}, {0, 0, (far + near) / (near - far), 1}};
}
template <typename T>
constexpr inline mat<T, 4, 4> make_perspective_mat(
    T fovy, T aspect, T near, T far) {
    auto tg = tan(fovy / 2);
    return {{1 / (aspect * tg), 0, 0, 0}, {0, 1 / tg, 0, 0},
        {0, 0, (far + near) / (near - far), -1},
        {0, 0, 2 * far * near / (near - far), 0}};
}
template <typename T>
constexpr inline mat<T, 4, 4> make_perspective_mat(T fovy, T aspect, T near) {
    auto tg = tan(fovy / 2);
    return {{1 / (aspect * tg), 0, 0, 0}, {0, 1 / tg, 0, 0}, {0, 0, -1, -1},
        {0, 0, 2 * near, 0}};
}

// Rotation conversions.
template <typename T>
constexpr inline pair<vec<T, 3>, T> make_rotation_axisangle(
    const vec<T, 4>& quat) {
    return {normalize(vec<T, 3>{quat[0], quat[1], quat[2]}), 2 * acos(quat[3])};
}
template <typename T>
constexpr inline vec<T, 4> make_rotation_quat(const vec<T, 3>& axis, T angle) {
    auto len = length(axis);
    if (!len) return {0, 0, 0, 1};
    return vec4f{sin(angle / 2) * axis[0] / len, sin(angle / 2) * axis[1] / len,
        sin(angle / 2) * axis[2] / len, cos(angle / 2)};
}
template <typename T>
constexpr inline vec<T, 4> make_rotation_quat(const vec<T, 4>& axisangle) {
    return make_rotation_quat(
        vec<T, 3>{axisangle[0], axisangle[1], axisangle[2]}, axisangle[3]);
}

// Turntable and FPS Camera navigation.
template <typename T>
inline void update_camera_turntable(vec<T, 3>& from, vec<T, 3>& to,
    vec<T, 3>& up, const vec<T, 2>& rotate, T dolly, const vec<T, 2>& pan);
template <typename T>
inline void update_camera_turntable(frame<T, 3>& frame, T& focus,
    const vec<T, 2>& rotate, T dolly, const vec<T, 2>& pan);
template <typename T>
inline void update_camera_firstperson(
    frame<T, 3>& frame, const vec<T, 3>& transl, const vec<T, 2>& rotate);

// Computes the image uv coordinates corresponding to the view parameters.
// Returns negative coordinates if out of the image.
template <typename T>
inline vec<int, 2> get_image_coords(const vec<T, 2>& mouse_pos,
    const vec<T, 2>& center, T scale, const vec<int, 2>& txt_size) {
    auto xyf = (mouse_pos - center) / scale;
    return {(int)round(xyf[0] + txt_size[0] / 2),
        (int)round(xyf[1] + txt_size[1] / 2)};
}

// Center image and autofit.
inline void update_image_view(vec2f& center, float& scale,
    const vec2i& image_size, const vec2i& window_size, bool zoom_to_fit) {
    if (zoom_to_fit) {
        scale  = min(window_size[0] / (float)image_size[0],
            window_size[1] / (float)image_size[1]);
        center = {(float)window_size[0] / 2, (float)window_size[1] / 2};
    } else {
        if (window_size[0] >= image_size[0] * scale)
            center[0] = window_size[0] / 2;
        if (window_size[1] >= image_size[1] * scale)
            center[1] = window_size[1] / 2;
    }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR MATRICES
// -----------------------------------------------------------------------------
namespace yocto {

// Matrix diagonals and transposes.
template <typename T>
constexpr inline mat<T, 1, 1> transpose(const mat<T, 1, 1>& a) {
    return {{a[0]}};
}
template <typename T>
constexpr inline mat<T, 2, 2> transpose(const mat<T, 2, 2>& a) {
    return {{a[0][0], a[1][0]}, {a[0][1], a[1][1]}};
}
template <typename T>
constexpr inline mat<T, 3, 3> transpose(const mat<T, 3, 3>& a) {
    return {
        {a[0][0], a[1][0], a[2][0]},
        {a[0][1], a[1][1], a[2][1]},
        {a[0][2], a[1][2], a[2][2]},
    };
}
template <typename T>
constexpr inline mat<T, 4, 4> transpose(const mat<T, 4, 4>& a) {
    return {
        {a[0][0], a[1][0], a[2][0], a[3][0]},
        {a[0][1], a[1][1], a[2][1], a[3][1]},
        {a[0][2], a[1][2], a[2][2], a[3][2]},
        {a[0][3], a[1][3], a[2][3], a[3][3]},
    };
}

// Matrix adjugates, determinant and inverses.
template <typename T>
constexpr inline mat<T, 1, 1> adjugate(const mat<T, 1, 1>& a) {
    return {{a[0]}};
}
template <typename T>
constexpr inline mat<T, 2, 2> adjugate(const mat<T, 2, 2>& a) {
    return {{a[1][1], -a[0][1]}, {-a[1][0], a[0][0]}};
}
template <typename T>
constexpr inline mat<T, 3, 3> adjugate(const mat<T, 3, 3>& a) {
    return {
        {
            a[1][1] * a[2][2] - a[2][1] * a[1][2],
            a[2][1] * a[0][2] - a[0][1] * a[2][2],
            a[0][1] * a[1][2] - a[1][1] * a[0][2],
        },
        {
            a[1][2] * a[2][0] - a[2][2] * a[1][0],
            a[2][2] * a[0][0] - a[0][2] * a[2][0],
            a[0][2] * a[1][0] - a[1][2] * a[0][0],
        },
        {
            a[1][0] * a[2][1] - a[2][0] * a[1][1],
            a[2][0] * a[0][1] - a[0][0] * a[2][1],
            a[0][0] * a[1][1] - a[1][0] * a[0][1],
        },
    };
}
template <typename T>
constexpr inline mat<T, 4, 4> adjugate(const mat<T, 4, 4>& a) {
    return {
        {
            a[1][1] * a[2][2] * a[3][3] + a[3][1] * a[1][2] * a[2][3] +
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
                a[2][1] * a[0][2] * a[1][3] - a[1][1] * a[2][2] * a[0][3],
        },
        {
            a[1][2] * a[3][3] * a[2][0] + a[2][2] * a[1][3] * a[3][0] +
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
                a[1][2] * a[0][3] * a[2][0] - a[2][2] * a[1][3] * a[0][0],
        },
        {
            a[1][3] * a[2][0] * a[3][1] + a[3][3] * a[1][0] * a[2][1] +
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
                a[2][3] * a[0][0] * a[1][1] - a[1][3] * a[2][0] * a[0][1],
        },
        {
            a[1][0] * a[3][1] * a[2][2] + a[2][0] * a[1][1] * a[3][2] +
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
                a[1][0] * a[0][1] * a[2][2] - a[2][0] * a[1][1] * a[0][2],
        },
    };
}
template <typename T>
constexpr inline T determinant(const mat<T, 1, 1>& a) {
    return a[0];
}
template <typename T>
constexpr inline T determinant(const mat<T, 2, 2>& a) {
    return a[0][0] * a[1][1] - a[0][1] * a[1][0];
}
template <typename T>
constexpr inline T determinant(const mat<T, 3, 3>& a) {
    return a[0][0] * (a[1][1] * a[2][2] - a[2][1] * a[1][2]) +
           a[0][1] * (a[1][2] * a[2][0] - a[2][2] * a[1][0]) +
           a[0][2] * (a[1][0] * a[2][1] - a[2][0] * a[1][1]);
}
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

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR UI UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Turntable for UI navigation.
template <typename T>
inline void update_camera_turntable(vec<T, 3>& from, vec<T, 3>& to,
    vec<T, 3>& up, const vec<T, 2>& rotate, T dolly, const vec<T, 2>& pan) {
    // rotate if necessary
    if (rotate[0] || rotate[1]) {
        auto z     = normalize(to - from);
        auto lz    = length(to - from);
        auto phi   = atan2(z[2], z[0]) + rotate[0];
        auto theta = acos(z[1]) + rotate[1];
        theta      = clamp(theta, 0.001f, pif - 0.001f);
        auto nz    = vec<T, 3>{sin(theta) * cos(phi) * lz, cos(theta) * lz,
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
    if (pan[0] || pan[1]) {
        auto z = normalize(to - from);
        auto x = normalize(cross(up, z));
        auto y = normalize(cross(z, x));
        auto t = vec<T, 3>{pan[0] * x[0] + pan[1] * y[0],
            pan[0] * x[1] + pan[1] * y[1], pan[0] * x[2] + pan[1] * y[2]};
        from += t;
        to += t;
    }
}

// Turntable for UI navigation.
template <typename T>
inline void update_camera_turntable(frame<T, 3>& frame, T& focus,
    const vec<T, 2>& rotate, T dolly, const vec<T, 2>& pan) {
    // rotate if necessary
    if (rotate != zero_vec2f) {
        auto phi   = atan2(frame[2][2], frame[2][0]) + rotate[0];
        auto theta = acos(frame[2][1]) + rotate[1];
        theta      = clamp(theta, 0.001f, pif - 0.001f);
        auto new_z = vec<T, 3>{
            sin(theta) * cos(phi), cos(theta), sin(theta) * sin(phi)};
        auto new_center = frame.origin - frame.axes[2] * focus;
        auto new_o      = new_center + new_z * focus;
        frame           = make_lookat_frame(new_o, new_center, {0, 1, 0});
        focus           = length(new_o - new_center);
    }

    // pan if necessary
    if (dolly) {
        auto c       = frame.origin - frame.axes[2] * focus;
        focus        = max(focus * (1 + dolly), 0.001f);
        frame.origin = c + frame.axes[2] * focus;
    }

    // pan if necessary
    if (pan[0] || pan[1]) {
        frame.origin += frame.axes[0] * pan[0] + frame.axes[1] * pan[1];
    }
}

// FPS camera for UI navigation for a frame parametrization.
template <typename T>
inline void update_camera_firstperson(const frame<T, 3>& frame,
    const vec<T, 3>& transl, const vec<T, 2>& rotate) {
    // https://gamedev.stackexchange.com/questions/30644/how-to-keep-my-quaternion-using-fps-camera-from-tilting-and-messing-up
    auto y = vec<T, 3>{0, 1, 0};
    auto z = orthonormalize(frame.axes[2], y);
    auto x = cross(y, z);

    auto rot = make_rotation_frame(vec<T, 3>{1, 0, 0}, rotate[1]) *
               frame3f{frame[0], frame[1], frame[2], vec<T, 3>{0, 0, 0}} *
               make_rotation_frame(vec<T, 3>{0, 1, 0}, rotate[0]);
    auto pos = frame.origin + transl[0] * x + transl[1] * y + transl[2] * z;

    frame = {rot[0], rot[1], rot[2], pos};
}

}  // namespace yocto

#endif
