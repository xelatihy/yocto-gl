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
// in graphics. In particular, we support 1-4 dimensional vectors
// coordinates (`vec<T, 1>`, `vec<T, 2>`, `vec<T, 3>`, `vec<T, 4>`).
//
// We support 2-4 dimensional matrices (`mat<T, 2, 2>`, `mat<T, 3, 3>`,
// `mat<T, 4, 4>`) with matrix-matrix and matrix-vector products, transposes and
// inverses. Matrices are stored in column-major order and are accessed and
// constructed by column. The one dimensional version is for completeness only.
//
// To represent transformations, most of the library facilities prefer the use
// coordinate frames, aka rigid transforms, represented as `frame2f` and
// `frame3f`. The structure store three coordinate axes and the origin.
// This is equivalent to a rigid transform written as a column-major affine
// matrix. Transform operations are fater with this representation.
//
// We represent bounding boxes in 2-3 dimensions with `bbox2f`, `bbox3f`
// Each bounding box support construction from points and other bounding box. 
// We provide operations to compute bounds for points, lines, triangles and 
// quads.
//
// For both matrices and frames we support transform operations for points,
// vectors and directions (`transform_point()`, `transform_vector()`,
// `transform_direction()`). Transform matrices and frames can be
// constructed from basic translation, rotation and scaling, e.g. with
// `make_translation_mat()` or `make_translation_frame()` respectively, etc.
// For rotation we support axis-angle and quaternions, with slerp.
//
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2019 Fabio Pellacini
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

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <algorithm>
#include <array>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <functional>
#include <limits>
#include <memory>
#include <numeric>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

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
using std::log2;
using std::pow;
using std::round;
using std::sin;
using std::sqrt;
using std::swap;
using std::tan;
using std::to_string;

using std::array;
using std::exception;
using std::function;
using std::invalid_argument;
using std::make_shared;
using std::make_unique;
using std::numeric_limits;
using std::out_of_range;
using std::pair;
using std::runtime_error;
using std::shared_ptr;
using std::string;
using std::string_view;
using std::tuple;
using std::unique_ptr;
using std::unordered_map;
using std::unordered_set;
using std::vector;
using namespace std::literals::string_literals;
using namespace std::literals::string_view_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// MATH CONSTANTS AND FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

using byte = unsigned char;
using uint = unsigned int;

constexpr double pi  = 3.14159265358979323846;
constexpr float  pif = (float)pi;

template <typename T>
constexpr T type_max = numeric_limits<T>::max();
template <typename T>
constexpr T type_min = numeric_limits<T>::lowest();

constexpr auto int_max   = type_max<int>;
constexpr auto int_min   = type_min<int>;
constexpr auto float_max = type_max<float>;
constexpr auto float_min = type_min<float>;

template <typename T>
constexpr T min(T x, T y) {
    return (x < y) ? x : y;
}
template <typename T>
constexpr T max(T x, T y) {
    return (x > y) ? x : y;
}
template <typename T>
constexpr T clamp(T x, T min_, T max_) {
    return min(max(x, min_), max_);
}
template <typename T>
constexpr T clamp01(T x) {
    return min(max(x, (T)0), (T)1);
}
template <typename T, typename T1>
constexpr T lerp(const T& a, const T& b, const T1& u) {
    return a * (1 - u) + b * u;
}
template <typename T, typename T1>
constexpr T bilerp(const T& c00, const T& c10, const T& c11, const T& c01,
    const T1& u, const T1& v) {
    return c00 * (1 - u) * (1 - v) + c10 * u * (1 - v) + c01 * (1 - u) * v +
           c11 * u * v;
}
template <typename T, typename T1>
constexpr T bias(const T& a, const T1& bias) {
    return a / ((1 / bias - 2) * (1 - a) + 1);
}
template <typename T, typename T1>
constexpr T gain(const T& a, const T1& gain) {
    if (a < (T)0.5) {
        return bias(a * 2, gain) / 2;
    } else {
        return bias(a * 2 - 1, 1 - gain) / 2 + (T)0.5;
    }
}
constexpr int pow2(int x) { return 1 << x; }
template <typename T>
inline T radians(T x) {
    return x * (T)pi / 180;
}
template <typename T>
inline T degrees(T x) {
    return x * 180 / (T)pi;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// VECTORS
// -----------------------------------------------------------------------------
namespace yocto {

// Small size vectors.
template <typename T, int N>
struct vec;

template <typename T>
struct vec<T, 1> {
    T x;

    constexpr vec() : x{0} {}
    constexpr vec(T x) : x{x} {}

    constexpr T&       operator[](int i) { return (&x)[i]; }
    constexpr const T& operator[](int i) const { return (&x)[i]; }
};

template <typename T>
struct vec<T, 2> {
    T x, y;

    constexpr vec() : x{0}, y{0} {}
    constexpr vec(T x, T y) : x{x}, y{y} {}
    constexpr vec(const vec<T, 1>& v, T y) : x{v.x}, y{y} {}
    constexpr explicit vec(T v) : x{v}, y{v} {}

    constexpr T&       operator[](int i) { return (&x)[i]; }
    constexpr const T& operator[](int i) const { return (&x)[i]; }
};

template <typename T>
struct vec<T, 3> {
    T x, y, z;

    constexpr vec() : x{0}, y{0}, z{0} {}
    constexpr vec(T x, T y, T z) : x{x}, y{y}, z{z} {}
    constexpr vec(const vec<T, 2>& v, T z) : x{v.x}, y{v.y}, z{z} {}
    constexpr explicit vec(T v) : x{v}, y{v}, z{v} {}

    constexpr T&       operator[](int i) { return (&x)[i]; }
    constexpr const T& operator[](int i) const { return (&x)[i]; }
};

template <typename T>
struct vec<T, 4> {
    T x, y, z, w;

    constexpr vec() : x{0}, y{0}, z{0}, w{0} {}
    constexpr vec(T x, T y, T z, T w) : x{x}, y{y}, z{z}, w{w} {}
    constexpr vec(const vec<T, 3>& v, T w) : x{v.x}, y{v.y}, z{v.z}, w{w} {}
    constexpr explicit vec(T v) : x{v}, y{v}, z{v}, w{v} {}

    constexpr T&       operator[](int i) { return (&x)[i]; }
    constexpr const T& operator[](int i) const { return (&x)[i]; }
};

// Typedefs
using vec1f = vec<float, 1>;
using vec2f = vec<float, 2>;
using vec3f = vec<float, 3>;
using vec4f = vec<float, 4>;
using vec1i = vec<int, 1>;
using vec2i = vec<int, 2>;
using vec3i = vec<int, 3>;
using vec4i = vec<int, 4>;
using vec1b = vec<byte, 1>;
using vec2b = vec<byte, 2>;
using vec3b = vec<byte, 3>;
using vec4b = vec<byte, 4>;
using vec1d = vec<double, 1>;
using vec2d = vec<double, 2>;
using vec3d = vec<double, 3>;
using vec4d = vec<double, 4>;

// Zero vector constants.
template <typename T, int N>
constexpr auto zero   = vec<T, N>{};
constexpr auto zero1f = vec1f{0};
constexpr auto zero2f = vec2f{0, 0};
constexpr auto zero3f = vec3f{0, 0, 0};
constexpr auto zero4f = vec4f{0, 0, 0, 0};
constexpr auto zero1i = vec1i{0};
constexpr auto zero2i = vec2i{0, 0};
constexpr auto zero3i = vec3i{0, 0, 0};
constexpr auto zero4i = vec4i{0, 0, 0, 0};
constexpr auto zero1b = vec1b{0};
constexpr auto zero2b = vec2b{0, 0};
constexpr auto zero3b = vec3b{0, 0, 0};
constexpr auto zero4b = vec4b{0, 0, 0, 0};

// Vector comparison operations.
template <typename T, int N>
constexpr bool operator==(const vec<T, N>& a, const vec<T, N>& b) {
    if constexpr (N == 1) {
        return a.x == b.x;
    } else if constexpr (N == 2) {
        return a.x == b.x && a.y == b.y;
    } else if constexpr (N == 3) {
        return a.x == b.x && a.y == b.y && a.z == b.z;
    } else if constexpr (N == 4) {
        return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
    } else {
        for (auto i = 0; i < N; i++)
            if (a[i] != b[i]) return false;
        return true;
    }
}
template <typename T, int N>
constexpr bool operator!=(const vec<T, N>& a, const vec<T, N>& b) {
    if constexpr (N == 1) {
        return a.x != b.x;
    } else if constexpr (N == 2) {
        return a.x != b.x || a.y != b.y;
    } else if constexpr (N == 3) {
        return a.x != b.x || a.y != b.y || a.z != b.z;
    } else if constexpr (N == 4) {
        return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w;
    } else {
        return !(a == b);
    }
}

// Element access
template <typename T, int N>
constexpr vec<T, 3>& xyz(vec<T, N>& a) {
    if constexpr (N == 3) {
        return a;
    } else if constexpr (N == 4) {
        return (vec<T, 3>&)a;
    } else {
        throw runtime_error("operation not supported");
    }
}
template <typename T, int N>
constexpr const vec<T, 3>& xyz(const vec<T, N>& a) {
    if constexpr (N == 3) {
        return a;
    } else if constexpr (N == 4) {
        return (const vec<T, 3>&)a;
    } else {
        throw runtime_error("operation not supported");
    }
}

// Vector operations.
template <typename T, int N>
constexpr vec<T, N> operator+(const vec<T, N>& a) {
    return a;
}
template <typename T, int N>
constexpr vec<T, N> operator-(const vec<T, N>& a) {
    if constexpr (N == 1) {
        return {-a.x};
    } else if constexpr (N == 2) {
        return {-a.x, -a.y};
    } else if constexpr (N == 3) {
        return {-a.x, -a.y, -a.z};
    } else if constexpr (N == 4) {
        return {-a.x, -a.y, -a.z, -a.w};
    } else {
        auto c = vec<T, N>{};
        for (auto i = 0; i < N; i++) c[i] = -a[i];
        return c;
    }
}
template <typename T, int N>
constexpr vec<T, N> operator+(const vec<T, N>& a, const vec<T, N>& b) {
    if constexpr (N == 1) {
        return {a.x + b.x};
    } else if constexpr (N == 2) {
        return {a.x + b.x, a.y + b.y};
    } else if constexpr (N == 3) {
        return {a.x + b.x, a.y + b.y, a.z + b.z};
    } else if constexpr (N == 4) {
        return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
    } else {
        auto c = vec<T, N>{};
        for (auto i = 0; i < N; i++) c[i] = a[i] + b[i];
        return c;
    }
}
template <typename T, int N, typename T1>
constexpr vec<T, N> operator+(const vec<T, N>& a, T1 b) {
    if constexpr (N == 1) {
        return {a.x + b};
    } else if constexpr (N == 2) {
        return {a.x + b, a.y + b};
    } else if constexpr (N == 3) {
        return {a.x + b, a.y + b, a.z + b};
    } else if constexpr (N == 4) {
        return {a.x + b, a.y + b, a.z + b, a.w + b};
    } else {
        auto c = vec<T, N>{};
        for (auto i = 0; i < N; i++) c[i] = a[i] + b;
        return c;
    }
}
template <typename T, int N, typename T1>
constexpr vec<T, N> operator+(T1 a, const vec<T, N>& b) {
    if constexpr (N == 1) {
        return {a + b.x};
    } else if constexpr (N == 2) {
        return {a + b.x, a + b.y};
    } else if constexpr (N == 3) {
        return {a + b.x, a + b.y, a + b.z};
    } else if constexpr (N == 4) {
        return {a + b.x, a + b.y, a + b.z, a + b.w};
    } else {
        auto c = vec<T, N>{};
        for (auto i = 0; i < N; i++) c[i] = a + b[i];
        return c;
    }
}
template <typename T, int N>
constexpr vec<T, N> operator-(const vec<T, N>& a, const vec<T, N>& b) {
    if constexpr (N == 1) {
        return {a.x - b.x};
    } else if constexpr (N == 2) {
        return {a.x - b.x, a.y - b.y};
    } else if constexpr (N == 3) {
        return {a.x - b.x, a.y - b.y, a.z - b.z};
    } else if constexpr (N == 4) {
        return {a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w};
    } else {
        auto c = vec<T, N>{};
        for (auto i = 0; i < N; i++) c[i] = a[i] - b[i];
        return c;
    }
}
template <typename T, int N, typename T1>
constexpr vec<T, N> operator-(const vec<T, N>& a, T1 b) {
    if constexpr (N == 1) {
        return {a.x - b};
    } else if constexpr (N == 2) {
        return {a.x - b, a.y - b};
    } else if constexpr (N == 3) {
        return {a.x - b, a.y - b, a.z - b};
    } else if constexpr (N == 4) {
        return {a.x - b, a.y - b, a.z - b, a.w - b};
    } else {
        auto c = vec<T, N>{};
        for (auto i = 0; i < N; i++) c[i] = a[i] - b;
        return c;
    }
}
template <typename T, int N, typename T1>
constexpr vec<T, N> operator-(T1 a, const vec<T, N>& b) {
    if constexpr (N == 1) {
        return {a - b.x};
    } else if constexpr (N == 2) {
        return {a - b.x, a - b.y};
    } else if constexpr (N == 3) {
        return {a - b.x, a - b.y, a - b.z};
    } else if constexpr (N == 4) {
        return {a - b.x, a - b.y, a - b.z, a - b.w};
    } else {
        auto c = vec<T, N>{};
        for (auto i = 0; i < N; i++) c[i] = a - b[i];
        return c;
    }
}
template <typename T, int N>
constexpr vec<T, N> operator*(const vec<T, N>& a, const vec<T, N>& b) {
    if constexpr (N == 1) {
        return {a.x * b.x};
    } else if constexpr (N == 2) {
        return {a.x * b.x, a.y * b.y};
    } else if constexpr (N == 3) {
        return {a.x * b.x, a.y * b.y, a.z * b.z};
    } else if constexpr (N == 4) {
        return {a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w};
    } else {
        auto c = vec<T, N>{};
        for (auto i = 0; i < N; i++) c[i] = a[i] * b[i];
        return c;
    }
}
template <typename T, int N, typename T1>
constexpr vec<T, N> operator*(const vec<T, N>& a, T1 b) {
    if constexpr (N == 1) {
        return {a.x * b};
    } else if constexpr (N == 2) {
        return {a.x * b, a.y * b};
    } else if constexpr (N == 3) {
        return {a.x * b, a.y * b, a.z * b};
    } else if constexpr (N == 4) {
        return {a.x * b, a.y * b, a.z * b, a.w * b};
    } else {
        auto c = vec<T, N>{};
        for (auto i = 0; i < N; i++) c[i] = a[i] * b;
        return c;
    }
}
template <typename T, int N, typename T1>
constexpr vec<T, N> operator*(T1 a, const vec<T, N>& b) {
    if constexpr (N == 1) {
        return {a * b.x};
    } else if constexpr (N == 2) {
        return {a * b.x, a * b.y};
    } else if constexpr (N == 3) {
        return {a * b.x, a * b.y, a * b.z};
    } else if constexpr (N == 4) {
        return {a * b.x, a * b.y, a * b.z, a * b.w};
    } else {
        auto c = vec<T, N>{};
        for (auto i = 0; i < N; i++) c[i] = a * b[i];
        return c;
    }
}
template <typename T, int N>
constexpr vec<T, N> operator/(const vec<T, N>& a, const vec<T, N>& b) {
    if constexpr (N == 1) {
        return {a.x / b.x};
    } else if constexpr (N == 2) {
        return {a.x / b.x, a.y / b.y};
    } else if constexpr (N == 3) {
        return {a.x / b.x, a.y / b.y, a.z / b.z};
    } else if constexpr (N == 4) {
        return {a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w};
    } else {
        auto c = vec<T, N>{};
        for (auto i = 0; i < N; i++) c[i] = a[i] / b[i];
        return c;
    }
}
template <typename T, int N, typename T1>
constexpr vec<T, N> operator/(const vec<T, N>& a, T1 b) {
    if constexpr (N == 1) {
        return {a.x / b};
    } else if constexpr (N == 2) {
        return {a.x / b, a.y / b};
    } else if constexpr (N == 3) {
        return {a.x / b, a.y / b, a.z / b};
    } else if constexpr (N == 4) {
        return {a.x / b, a.y / b, a.z / b, a.w / b};
    } else {
        auto c = vec<T, N>{};
        for (auto i = 0; i < N; i++) c[i] = a[i] / b;
        return c;
    }
}
template <typename T, int N, typename T1>
constexpr vec<T, N> operator/(T1 a, const vec<T, N>& b) {
    if constexpr (N == 1) {
        return {a / b.x};
    } else if constexpr (N == 2) {
        return {a / b.x, a / b.y};
    } else if constexpr (N == 3) {
        return {a / b.x, a / b.y, a / b.z};
    } else if constexpr (N == 4) {
        return {a / b.x, a / b.y, a / b.z, a / b.w};
    } else {
        auto c = vec<T, N>{};
        for (auto i = 0; i < N; i++) c[i] = a / b[i];
        return c;
    }
}

// Vector assignments
template <typename T, int N>
constexpr vec<T, N>& operator+=(vec<T, N>& a, const vec<T, N>& b) {
    return a = a + b;
}
template <typename T, int N, typename T1>
constexpr vec<T, N>& operator+=(vec<T, N>& a, T1 b) {
    return a = a + b;
}
template <typename T, int N>
constexpr vec<T, N>& operator-=(vec<T, N>& a, const vec<T, N>& b) {
    return a = a - b;
}
template <typename T, int N, typename T1>
constexpr vec<T, N>& operator-=(vec<T, N>& a, T1 b) {
    return a = a - b;
}
template <typename T, int N>
constexpr vec<T, N>& operator*=(vec<T, N>& a, const vec<T, N>& b) {
    return a = a * b;
}
template <typename T, int N, typename T1>
constexpr vec<T, N>& operator*=(vec<T, N>& a, T1 b) {
    return a = a * b;
}
template <typename T, int N>
constexpr vec<T, N>& operator/=(vec<T, N>& a, const vec<T, N>& b) {
    return a = a / b;
}
template <typename T, int N, typename T1>
constexpr vec<T, N>& operator/=(vec<T, N>& a, T1 b) {
    return a = a / b;
}

// Vector products and lengths.
template <typename T, int N>
constexpr T dot(const vec<T, N>& a, const vec<T, N>& b) {
    if constexpr (N == 1) {
        return a.x * b.x;
    } else if constexpr (N == 2) {
        return a.x * b.x + a.y * b.y;
    } else if constexpr (N == 3) {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    } else if constexpr (N == 4) {
        return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
    } else {
        auto c = T{0};
        for (auto i = 0; i < N; i++) c += a[i] * b[i];
        return c;
    }
}
template <typename T>
constexpr T cross(const vec<T, 2>& a, const vec<T, 2>& b) {
    return a.x * b.y - a.y * b.x;
}
template <typename T>
constexpr vec<T, 3> cross(const vec<T, 3>& a, const vec<T, 3>& b) {
    return {
        a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

template <typename T, int N>
constexpr T length(const vec<T, N>& a) {
    return sqrt(dot(a, a));
}
template <typename T, int N>
constexpr vec<T, N> normalize(const vec<T, N>& a) {
    auto l = length(a);
    return (l != 0) ? a / l : a;
}
template <typename T, int N>
constexpr T distance(const vec<T, N>& a, const vec<T, N>& b) {
    return length(a - b);
}
template <typename T, int N>
constexpr T distance_squared(const vec<T, N>& a, const vec<T, N>& b) {
    return dot(a - b, a - b);
}

template <typename T>
constexpr T angle(const vec<T, 3>& a, const vec<T, 3>& b) {
    return acos(clamp(dot(normalize(a), normalize(b)), (T)-1, (T)1));
}
template <typename T, typename T1>
constexpr vec<T, 4> slerp(const vec<T, 4>& a, const vec<T, 4>& b, T1 u) {
    // https://en.wikipedia.org/wiki/Slerp
    auto an = normalize(a), bn = normalize(b);
    auto d = dot(an, bn);
    if (d < 0) {
        bn = -bn;
        d  = -d;
    }
    if (d > (T)0.9995) return normalize(an + u * (bn - an));
    auto th = acos(clamp(d, (T)-1, (T)1));
    if (!th) return an;
    return an * (sin(th * (1 - u)) / sin(th)) + bn * (sin(th * u) / sin(th));
}

// Orthogonal vectors.
template <typename T>
constexpr vec<T, 3> orthogonal(const vec<T, 3>& v) {
    // http://lolengine.net/blog/2013/09/21/picking-orthogonal-vector-combing-coconuts)
    return abs(v.x) > abs(v.z) ? vec<T, 3>{-v.y, v.x, 0}
                               : vec<T, 3>{0, -v.z, v.y};
}
template <typename T>
constexpr vec<T, 3> orthonormalize(const vec<T, 3>& a, const vec<T, 3>& b) {
    return normalize(a - b * dot(a, b));
}

// Reflected and refracted vector.
template <typename T>
constexpr vec<T, 3> reflect(const vec<T, 3>& w, const vec<T, 3>& n) {
    return -w + 2 * dot(n, w) * n;
}
template <typename T, typename T1>
constexpr vec<T, 3> refract(const vec<T, 3>& w, const vec<T, 3>& n, T1 eta) {
    // auto k = 1.0 - eta * eta * (1.0 - dot(n, w) * dot(n, w));
    auto k = 1 - eta * eta * max((T)0, 1 - dot(n, w) * dot(n, w));
    if (k < 0) return {0, 0, 0};  // tir
    return -w * eta + (eta * dot(n, w) - sqrt(k)) * n;
}

// Max element and clamp.
template <typename T, int N>
constexpr vec<T, N> max(const vec<T, N>& a, T b) {
    if constexpr (N == 1) {
        return {max(a.x, b)};
    } else if constexpr (N == 2) {
        return {max(a.x, b), max(a.y, b)};
    } else if constexpr (N == 3) {
        return {max(a.x, b), max(a.y, b), max(a.z, b)};
    } else if constexpr (N == 4) {
        return {max(a.x, b), max(a.y, b), max(a.z, b), max(a.w, b)};
    } else {
        auto c = vec<T, N>{};
        for (auto i = 0; i < N; i++) c[i] = max(a[i], b);
        return c;
    }
}
template <typename T, int N>
constexpr vec<T, N> min(const vec<T, N>& a, T b) {
    if constexpr (N == 1) {
        return {min(a.x, b)};
    } else if constexpr (N == 2) {
        return {min(a.x, b), min(a.y, b)};
    } else if constexpr (N == 3) {
        return {min(a.x, b), min(a.y, b), min(a.z, b)};
    } else if constexpr (N == 4) {
        return {min(a.x, b), min(a.y, b), min(a.z, b), min(a.w, b)};
    } else {
        auto c = vec<T, N>{};
        for (auto i = 0; i < N; i++) c[i] = min(a[i], b);
        return c;
    }
}
template <typename T, int N>
constexpr vec<T, N> max(const vec<T, N>& a, const vec<T, N>& b) {
    if constexpr (N == 1) {
        return {max(a.x, b.x)};
    } else if constexpr (N == 2) {
        return {max(a.x, b.x), max(a.y, b.y)};
    } else if constexpr (N == 3) {
        return {max(a.x, b.x), max(a.y, b.y), max(a.z, b.z)};
    } else if constexpr (N == 4) {
        return {max(a.x, b.x), max(a.y, b.y), max(a.z, b.z), max(a.w, b.w)};
    } else {
        auto c = vec<T, N>{};
        for (auto i = 0; i < N; i++) c[i] = max(a[i], b[i]);
        return c;
    }
}
template <typename T, int N>
constexpr vec<T, N> min(const vec<T, N>& a, const vec<T, N>& b) {
    if constexpr (N == 1) {
        return {min(a.x, b.x)};
    } else if constexpr (N == 2) {
        return {min(a.x, b.x), min(a.y, b.y)};
    } else if constexpr (N == 3) {
        return {min(a.x, b.x), min(a.y, b.y), min(a.z, b.z)};
    } else if constexpr (N == 4) {
        return {min(a.x, b.x), min(a.y, b.y), min(a.z, b.z), min(a.w, b.w)};
    } else {
        auto c = vec<T, N>{};
        for (auto i = 0; i < N; i++) c[i] = min(a[i], b[i]);
        return c;
    }
}
template <typename T, int N, typename T1, typename T2>
constexpr vec<T, N> clamp(const vec<T, N>& x, T1 min, T2 max) {
    if constexpr (N == 1) {
        return {clamp(x.x, min, max)};
    } else if constexpr (N == 2) {
        return {clamp(x.x, min, max), clamp(x.y, min, max)};
    } else if constexpr (N == 3) {
        return {
            clamp(x.x, min, max), clamp(x.y, min, max), clamp(x.z, min, max)};
    } else if constexpr (N == 4) {
        return {clamp(x.x, min, max), clamp(x.y, min, max),
            clamp(x.z, min, max), clamp(x.w, min, max)};
    } else {
        auto c = vec<T, N>{};
        for (auto i = 0; i < N; i++) c[i] = clamp(x[i], min, max);
        return c;
    }
}
template <typename T, int N>
constexpr vec<T, N> clamp01(const vec<T, N>& x) {
    if constexpr (N == 1) {
        return {clamp01(x.x)};
    } else if constexpr (N == 2) {
        return {clamp01(x.x), clamp01(x.y)};
    } else if constexpr (N == 3) {
        return {clamp01(x.x), clamp01(x.y), clamp01(x.z)};
    } else if constexpr (N == 4) {
        return {clamp01(x.x), clamp01(x.y), clamp01(x.z), clamp01(x.w)};
    } else {
        auto c = vec<T, N>{};
        for (auto i = 0; i < N; i++) c[i] = clamp01(x[i]);
        return c;
    }
}

template <typename T, int N>
constexpr T max(const vec<T, N>& a) {
    if constexpr (N == 1) {
        return a.x;
    } else if constexpr (N == 2) {
        return max(a.x, a.y);
    } else if constexpr (N == 3) {
        return max(max(a.x, a.y), a.z);
    } else if constexpr (N == 4) {
        return max(max(max(a.x, a.y), a.z), a.w);
    } else {
        auto m = a[0];
        for (auto i = 1; i < N; i++) m = max(a[i], m);
        return m;
    }
}
template <typename T, int N>
constexpr T min(const vec<T, N>& a) {
    if constexpr (N == 1) {
        return a.x;
    } else if constexpr (N == 2) {
        return min(a.x, a.y);
    } else if constexpr (N == 3) {
        return min(min(a.x, a.y), a.z);
    } else if constexpr (N == 4) {
        return min(min(min(a.x, a.y), a.z), a.w);
    } else {
        auto m = a[0];
        for (auto i = 1; i < N; i++) m = min(a[i], m);
        return m;
    }
}
template <typename T, int N>
constexpr T sum(const vec<T, N>& a) {
    if constexpr (N == 1) {
        return a.x;
    } else if constexpr (N == 2) {
        return a.x + a.y;
    } else if constexpr (N == 3) {
        return a.x + a.y + a.z;
    } else if constexpr (N == 4) {
        return a.x + a.y + a.z + a.w;
    } else {
        auto m = a[0];
        for (auto i = 1; i < N; i++) m += a[i];
        return m;
    }
}
template <typename T, int N>
constexpr T mean(const vec<T, N>& a) {
    return sum(a) / N;
}
template <typename T, int N>
constexpr int max_element(const vec<T, N>& a) {
    auto pos = 0;
    for (auto i = 1; i < N; i++) {
        if (a[i] > a[pos]) pos = i;
    }
    return pos;
}
template <typename T, int N>
constexpr int min_element(const vec<T, N>& a) {
    auto pos = 0;
    for (auto i = 1; i < N; i++) {
        if (a[i] < a[pos]) pos = i;
    }
    return pos;
}

// Apply a unary function to all vector elements
template <typename T, int N, typename Func>
constexpr vec<T, N> apply(const Func& func, const vec<T, N>& a) {
    if constexpr (N == 1) {
        return {func(a.x)};
    } else if constexpr (N == 2) {
        return {func(a.x), func(a.y)};
    } else if constexpr (N == 3) {
        return {func(a.x), func(a.y), func(a.z)};
    } else if constexpr (N == 4) {
        return {func(a.x), func(a.y), func(a.z), func(a.w)};
    } else {
        auto c = vec<T, N>{};
        for (auto i = 0; i < N; i++) c[i] = func(a[i]);
        return c;
    }
}
// Apply a binary function to all vector elements
template <typename T, int N, typename Func>
constexpr vec<T, N> apply(const Func& func, const vec<T, N>& a, T b) {
    if constexpr (N == 1) {
        return {func(a.x, b)};
    } else if constexpr (N == 2) {
        return {func(a.x, b), func(a.y, b)};
    } else if constexpr (N == 3) {
        return {func(a.x, b), func(a.y, b), func(a.z, b)};
    } else if constexpr (N == 4) {
        return {func(a.x, b), func(a.y, b), func(a.z, b), func(a.w, b)};
    } else {
        auto c = vec<T, N>{};
        for (auto i = 0; i < N; i++) c[i] = func(a[i], b);
        return c;
    }
}
// Apply a binary function to all vector elements
template <typename T, int N, typename Func>
constexpr vec<T, N> apply(
    const Func& func, const vec<T, N>& a, const vec<T, N>& b) {
    if constexpr (N == 1) {
        return {func(a.x, b.x)};
    } else if constexpr (N == 2) {
        return {func(a.x, b.x), func(a.y, b.x)};
    } else if constexpr (N == 3) {
        return {func(a.x, b.x), func(a.y, b.y), func(a.z, b.z)};
    } else if constexpr (N == 4) {
        return {func(a.x, b.x), func(a.y, b.y), func(a.z, b.z), func(a.w, b.w)};
    } else {
        auto c = vec<T, N>{};
        for (auto i = 0; i < N; i++) c[i] = func(a[i], b);
        return c;
    }
}

// Functions applied to vector elements
template <typename T, int N>
constexpr vec<T, N> abs(const vec<T, N>& a) {
    return apply([](const T& a) { return abs(a); }, a);
};
template <typename T, int N>
constexpr vec<T, N> sqrt(const vec<T, N>& a) {
    return apply([](const T& a) { return sqrt(a); }, a);
};
template <typename T, int N, typename T1>
constexpr vec<T, N> pow(const vec<T, N>& a, const T1& b) {
    return apply([](const T& a, const T& b) { return pow(a, b); }, a, b);
};
template <typename T, int N>
constexpr vec<T, N> exp(const vec<T, N>& a) {
    return apply([](const T& a) { return exp(a); }, a);
};
template <typename T, int N>
constexpr vec<T, N> log(const vec<T, N>& a) {
    return apply([](const T& a) { return log(a); }, a);
};
template <typename T, int N>
constexpr vec<T, N> exp2(const vec<T, N>& a) {
    return apply([](const T& a) { return exp2(a); }, a);
};
template <typename T, int N>
constexpr vec<T, N> log2(const vec<T, N>& a) {
    return apply([](const T& a) { return log2(a); }, a);
};
template <typename T, int N, typename T1>
constexpr vec<T, N> gain(const vec<T, N>& a, const T1& b) {
    return apply([](const T& a, const T& b) { return gain(a, b); }, a, b);
};

template <typename T, int N>
constexpr bool isfinite(const vec<T, N>& a) {
    if constexpr (N == 1) {
        return isfinite(a.x);
    } else if constexpr (N == 2) {
        return isfinite(a.x) && isfinite(a.y);
    } else if constexpr (N == 3) {
        return isfinite(a.x) && isfinite(a.y) && isfinite(a.z);
    } else if constexpr (N == 4) {
        return isfinite(a.x) && isfinite(a.y) && isfinite(a.z) && isfinite(a.w);
    } else {
        throw invalid_argument("not implemented");
    }
};

// Quaternion operatons represented as xi + yj + zk + w
// const auto identity_quat4f = vec<float, 4>{0, 0, 0, 1};
template <typename T, typename T1>
constexpr vec<T, 4> quat_mul(const vec<T, 4>& a, T1 b) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
}
template <typename T>
constexpr vec<T, 4> quat_mul(const vec<T, 4>& a, const vec<T, 4>& b) {
    return {a.x * b.w + a.w * b.x + a.y * b.w - a.z * b.y,
        a.y * b.w + a.w * b.y + a.z * b.x - a.x * b.z,
        a.z * b.w + a.w * b.z + a.x * b.y - a.y * b.x,
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z};
}
template <typename T>
constexpr vec<T, 4> quat_conjugate(const vec<T, 4>& a) {
    return {-a.x, -a.y, -a.z, a.w};
}
template <typename T>
constexpr vec<T, 4> quat_inverse(const vec<T, 4>& a) {
    return quat_conjugate(a) / dot(a, a);
}

}  // namespace yocto

namespace std {

// Hash functor for vector for use with unordered_map
template <typename T, int N>
struct hash<yocto::vec<T, N>> {
    static constexpr std::hash<T> hasher = std::hash<T>();

    size_t operator()(const yocto::vec<T, N>& v) const {
        auto h = (size_t)0;
        for (auto i = 0; i < N; i++)
            h ^= hasher((&v.x)[i]) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

}  // namespace std

// -----------------------------------------------------------------------------
// MATRICES
// -----------------------------------------------------------------------------
namespace yocto {

// Small Fixed-size matrices stored in column major format.
template <typename T, int N, int M>
struct mat;

// Small Fixed-size matrices stored in column major format.
template <typename T, int N>
struct mat<T, N, 1> {
    vec<T, N> x;

    constexpr mat() : x{} {}
    constexpr mat(const vec<T, N>& x) : x{x} {}

    constexpr vec<T, N>&       operator[](int i) { return (&x)[i]; }
    constexpr const vec<T, N>& operator[](int i) const { return (&x)[i]; }
};

// Small Fixed-size matrices stored in column major format.
template <typename T, int N>
struct mat<T, N, 2> {
    vec<T, N> x, y;

    constexpr mat() : x{}, y{} {}
    constexpr mat(const vec<T, N>& x, const vec<T, N>& y) : x{x}, y{y} {}

    constexpr vec<T, N>&       operator[](int i) { return (&x)[i]; }
    constexpr const vec<T, N>& operator[](int i) const { return (&x)[i]; }
};

// Small Fixed-size matrices stored in column major format.
template <typename T, int N>
struct mat<T, N, 3> {
    vec<T, N> x, y, z;

    constexpr mat() : x{}, y{}, z{} {}
    constexpr mat(const vec<T, N>& x, const vec<T, N>& y, const vec<T, N>& z)
        : x{x}, y{y}, z{z} {}

    constexpr vec<T, N>&       operator[](int i) { return (&x)[i]; }
    constexpr const vec<T, N>& operator[](int i) const { return (&x)[i]; }
};

// Small Fixed-size matrices stored in column major format.
template <typename T, int N>
struct mat<T, N, 4> {
    vec<T, N> x, y, z, w;

    constexpr mat() : x{}, y{}, z{}, w{} {}
    constexpr mat(const vec<T, N>& x, const vec<T, N>& y, const vec<T, N>& z,
        const vec<T, N>& w)
        : x{x}, y{y}, z{z}, w{w} {}

    constexpr vec<T, N>&       operator[](int i) { return (&x)[i]; }
    constexpr const vec<T, N>& operator[](int i) const { return (&x)[i]; }
};

// Typedefs
using mat2f = mat<float, 2, 2>;
using mat3f = mat<float, 3, 3>;
using mat4f = mat<float, 4, 4>;
using mat2d = mat<double, 2, 2>;
using mat3d = mat<double, 3, 3>;
using mat4d = mat<double, 4, 4>;

// Identity matrix
template <typename T, int N>
constexpr mat<T, N, N> make_identity_mat();

// Identity matrices constants.
constexpr auto identity_mat2f = mat2f{{1, 0}, {0, 1}};
constexpr auto identity_mat3f = mat3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
constexpr auto identity_mat4f = mat4f{
    {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};

// Identity matrix
template <typename T, int N>
constexpr mat<T, N, N> make_identity_mat() {
    if constexpr (N == 1) {
        return {{1}};
    } else if constexpr (N == 2) {
        return {{1, 0}, {0, 1}};
    } else if constexpr (N == 3) {
        return {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
    } else if constexpr (N == 4) {
        return {{1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};
    } else {
        auto m = mat<T, N, N>{};
        for (auto i = 0; i < N; i++) m[i][i] = 1;
        return m;
    }
}

// Matrix comparisons.
template <typename T, int N, int M>
constexpr bool operator==(const mat<T, N, M>& a, const mat<T, N, M>& b) {
    if constexpr (M == 1) {
        return a.x == b.x;
    } else if constexpr (M == 2) {
        return a.x == b.x && a.y == b.y;
    } else if constexpr (M == 3) {
        return a.x == b.x && a.y == b.y && a.z == b.z;
    } else if constexpr (M == 4) {
        return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
    } else {
        for (auto i = 0; i < M; i++)
            if (a[i] != b[i]) return false;
        return true;
    }
}
template <typename T, int N, int M>
constexpr bool operator!=(const mat<T, N, M>& a, const mat<T, N, M>& b) {
    return !(a == b);
}

// Matrix operations.
template <typename T, int N, int M>
constexpr mat<T, N, M> operator+(const mat<T, N, M>& a, const mat<T, N, M>& b) {
    if constexpr (M == 1) {
        return {a.x + b.x};
    } else if constexpr (M == 2) {
        return {a.x + b.x, a.y + b.y};
    } else if constexpr (M == 3) {
        return {a.x + b.x, a.y + b.y, a.z + b.z};
    } else if constexpr (M == 4) {
        return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
    } else {
        auto c = mat<T, N, M>{};
        for (auto i = 0; i < M; i++) c[i] = a[i] + b[i];
        return c;
    }
}
template <typename T, int N, int M, typename T1>
constexpr mat<T, N, M> operator*(const mat<T, N, M>& a, T1 b) {
    if constexpr (M == 1) {
        return {a.x * b};
    } else if constexpr (M == 2) {
        return {a.x * b, a.y * b};
    } else if constexpr (M == 3) {
        return {a.x * b, a.y * b, a.z * b};
    } else if constexpr (M == 4) {
        return {a.x * b, a.y * b, a.z * b, a.w * b};
    } else {
        auto c = mat<T, N, M>{};
        for (auto i = 0; i < M; i++) c[i] = a[i] * b;
        return c;
    }
}
template <typename T, int N, int M>
constexpr vec<T, N> operator*(const mat<T, N, M>& a, const vec<T, M>& b) {
    if constexpr (M == 1) {
        return a.x * b.x;
    } else if constexpr (M == 2) {
        return a.x * b.x + a.y * b.y;
    } else if constexpr (M == 3) {
        return a.x * b.x + a.y * b.y + a.z * b.z;
    } else if constexpr (M == 4) {
        return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
    } else {
        auto c = vec<T, N>{};
        for (auto i = 0; i < M; i++) c += a[i] + b[i];
        return c;
    }
}
template <typename T, int N, int M>
constexpr vec<T, M> operator*(const vec<T, N>& a, const mat<T, N, M>& b) {
    if constexpr (M == 1) {
        return {dot(a, b.x)};
    } else if constexpr (M == 2) {
        return {dot(a, b.x), dot(a, b.y)};
    } else if constexpr (M == 3) {
        return {dot(a, b.x), dot(a, b.y), dot(a, b.z)};
    } else if constexpr (M == 4) {
        return {dot(a, b.x), dot(a, b.y), dot(a, b.z), dot(a, b.w)};
    } else {
        auto c = vec<T, M>{};
        for (auto i = 0; i < M; i++) c[i] = dot(a, b[i]);
        return c;
    }
}
template <typename T, int N, int M, int K>
constexpr mat<T, N, M> operator*(const mat<T, N, K>& a, const mat<T, K, M>& b) {
    if constexpr (M == 1) {
        return {dot(a, b.x)};
    } else if constexpr (M == 2) {
        return {a * b.x, a * b.y};
    } else if constexpr (M == 3) {
        return {a * b.x, a * b.y, a * b.z};
    } else if constexpr (M == 4) {
        return {a * b.x, a * b.y, a * b.z, a * b.w};
    } else {
        auto c = mat<T, N, M>{};
        for (auto i = 0; i < M; i++) c[i] = a * b[i];
        return c;
    }
}

// Matrix assignments.
template <typename T, int N, int M>
constexpr mat<T, N, M>& operator+=(mat<T, N, M>& a, const mat<T, N, M>& b) {
    return a = a + b;
}
template <typename T, int N>
constexpr mat<T, N, N>& operator*=(mat<T, N, N>& a, const mat<T, N, N>& b) {
    return a = a * b;
}
template <typename T, int N, int M, typename T1>
constexpr mat<T, N, M>& operator*=(mat<T, N, M>& a, T1 b) {
    return a = a * b;
}

// Matrix diagonals and transposes.
template <typename T, int N>
constexpr vec<T, N> diagonal(const mat<T, N, N>& a) {
    if constexpr (N == 1) {
        return {a.x.x};
    } else if constexpr (N == 2) {
        return {a.x.x, a.y.y};
    } else if constexpr (N == 3) {
        return {a.x.x, a.y.y, a.z.z};
    } else if constexpr (N == 4) {
        return {a.x.x, a.y.y, a.z.z, a.w.w};
    } else {
        auto c = vec<T, N>{};
        for (auto i = 0; i < N; i++) c[i] += a[i][i];
        return c;
    }
}
template <typename T, int N, int M>
constexpr mat<T, M, N> transpose(const mat<T, N, M>& a) {
    if constexpr (N == 1 && M == 1) {
        return a;
    } else if constexpr (N == 2 && M == 2) {
        return {{a.x.x, a.y.x}, {a.x.y, a.y.y}};
    } else if constexpr (N == 3 && M == 3) {
        return {
            {a.x.x, a.y.x, a.z.x},
            {a.x.y, a.y.y, a.z.y},
            {a.x.z, a.y.z, a.z.z},
        };
    } else if constexpr (N == 4 && M == 4) {
        return {
            {a.x.x, a.y.x, a.z.x, a.w.x},
            {a.x.y, a.y.y, a.z.y, a.w.y},
            {a.x.z, a.y.z, a.z.z, a.w.z},
            {a.x.w, a.y.w, a.z.w, a.w.w},
        };
    } else {
        auto c = mat<T, M, N>{};
        for (auto i = 0; i < M; i++)
            for (auto j = 0; j < N; j++) c[j][i] = a[i][j];
    }
}

// Matrix adjoints, determinants and inverses.
template <typename T>
constexpr T determinant(const mat<T, 2, 2>& a) {
    return cross(a.x, a.y);
}
template <typename T>
constexpr T determinant(const mat<T, 3, 3>& a) {
    return dot(a.x, cross(a.y, a.z));
}
template <typename T>
constexpr mat<T, 2, 2> adjoint(const mat<T, 2, 2>& a) {
    return {{a.y.y, -a.x.y}, {-a.y.x, a.x.x}};
}
template <typename T>
constexpr mat<T, 3, 3> adjoint(const mat<T, 3, 3>& a) {
    return transpose(
        mat<T, 3, 3>{cross(a.y, a.z), cross(a.z, a.x), cross(a.x, a.y)});
}
template <typename T>
constexpr mat<T, 2, 2> inverse(const mat<T, 2, 2>& a) {
    return adjoint(a) * (1 / determinant(a));
}
template <typename T>
constexpr mat<T, 3, 3> inverse(const mat<T, 3, 3>& a) {
    return adjoint(a) * (1 / determinant(a));
}

// Constructs a basis from a direction
template <typename T>
constexpr mat<T, 3, 3> make_basis_fromz(const vec<T, 3>& v) {
    // https://graphics.pixar.com/library/OrthonormalB/paper.pdf
    auto z    = normalize(v);
    auto sign = copysignf(1.0f, z.z);
    auto a    = -1.0f / (sign + z.z);
    auto b    = z.x * z.y * a;
    auto x    = vec<T, 3>{1.0f + sign * z.x * z.x * a, sign * b, -sign * z.x};
    auto y    = vec<T, 3>{b, sign + z.y * z.y * a, -z.y};
    return {x, y, z};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// RIGID BODY TRANSFORMS/FRAMES
// -----------------------------------------------------------------------------
namespace yocto {

// Rigid frames stored as a column-major affine transform matrix.
struct frame2f {
    vec2f x = {1, 0};
    vec2f y = {0, 1};
    vec2f o = {0, 0};

    constexpr frame2f() : x{}, y{}, o{} {}
    constexpr frame2f(const vec2f& x, const vec2f& y, const vec2f& o)
        : x{x}, y{y}, o{o} {}
    constexpr explicit frame2f(const vec2f& o) : x{1, 0}, y{0, 1}, o{o} {}
    constexpr frame2f(const mat2f& m, const vec2f& t)
        : x{m.x}, y{m.y}, o{t} {}
    constexpr explicit frame2f(const mat3f& m)
        : x{m.x.x, m.x.y}, y{m.y.x, m.y.y}, o{m.z.x, m.z.y} {}
    constexpr operator mat3f() const { return {{x, 0}, {y, 0}, {o, 1}}; }

    constexpr vec2f&       operator[](int i) { return (&x)[i]; }
    constexpr const vec2f& operator[](int i) const { return (&x)[i]; }

    constexpr mat2f&       m() { return *(mat2f*)&x; }
    constexpr const mat2f& m() const { return *(mat2f*)&x; }
    constexpr vec2f&          t() { return o; }
    constexpr const vec2f&    t() const { return o; }
};

// Rigid frames stored as a column-major affine transform matrix.
struct frame3f {
    vec3f x = {1,0,0};
    vec3f y = {0,1,0};
    vec3f z = {0,0,1};
    vec3f o = {0,0,0};

    constexpr frame3f() : x{}, y{}, z{}, o{} {}
    constexpr frame3f(const vec3f& x, const vec3f& y, const vec3f& z,
        const vec3f& o)
        : x{x}, y{y}, z{z}, o{o} {}
    constexpr explicit frame3f(const vec3f& o)
        : x{1, 0, 0}, y{0, 1, 0}, z{0, 0, 1}, o{o} {}
    constexpr frame3f(const mat3f& m, const vec3f& t)
        : x{m.x}, y{m.y}, z{m.z}, o{t} {}
    constexpr explicit frame3f(const mat4f& m)
        : x{m.x.x, m.x.y, m.x.z}
        , y{m.y.x, m.y.y, m.y.z}
        , z{m.z.x, m.z.y, m.z.z}
        , o{m.w.x, m.w.y, m.w.z} {}
    constexpr operator mat4f() const {
        return {{x, 0}, {y, 0}, {z, 0}, {o, 1}};
    }

    constexpr vec3f&       operator[](int i) { return (&x)[i]; }
    constexpr const vec3f& operator[](int i) const { return (&x)[i]; }

    constexpr mat3f&       m() { return *(mat3f*)&x; }
    constexpr const mat3f& m() const { return *(mat3f*)&x; }
    constexpr vec3f&          t() { return o; }
    constexpr const vec3f&    t() const { return o; }
};

// Indentity frames.
constexpr auto identity_frame2f = frame2f{{1, 0}, {0, 1}, {0, 0}};
constexpr auto identity_frame3f = frame3f{
    {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};

// Frame properties
inline const mat2f& linear_component(const frame2f& a) {
    return (const mat2f&)a;
}

// Frame comparisons.
constexpr bool operator==(const frame2f& a, const frame2f& b) {
        return a.x == b.x && a.y == b.y && a.o == b.o;
}
constexpr bool operator!=(const frame2f& a, const frame2f& b) {
    return !(a == b);
}

// Frame composition, equivalent to affine matrix product.
constexpr frame2f operator*(const frame2f& a, const frame2f& b) {
    return {a.m() * b.m(), a.m() * b.o + a.o};
}
constexpr frame2f& operator*=(frame2f& a, const frame2f& b) {
    return a = a * b;
}

// Frame inverse, equivalent to rigid affine inverse.
constexpr frame2f inverse(const frame2f& a, bool non_rigid = false) {
    if (non_rigid) {
        auto minv = inverse(a.m());
        return {minv, -(minv * a.o)};
    } else {
        auto minv = transpose(a.m());
        return {minv, -(minv * a.o)};
    }
}

// Frame properties
inline const mat3f& linear_component(const frame3f& a) {
    return (const mat3f&)a;
}

// Frame comparisons.
constexpr bool operator==(const frame3f& a, const frame3f& b) {
        return a.x == b.x && a.y == b.y && a.z == b.z && a.o == b.o;
}
constexpr bool operator!=(const frame3f& a, const frame3f& b) {
    return !(a == b);
}

// Frame composition, equivalent to affine matrix product.
constexpr frame3f operator*(const frame3f& a, const frame3f& b) {
    return {a.m() * b.m(), a.m() * b.o + a.o};
}
constexpr frame3f& operator*=(frame3f& a, const frame3f& b) {
    return a = a * b;
}

// Frame inverse, equivalent to rigid affine inverse.
constexpr frame3f inverse(const frame3f& a, bool non_rigid = false) {
    if (non_rigid) {
        auto minv = inverse(a.m());
        return {minv, -(minv * a.o)};
    } else {
        auto minv = transpose(a.m());
        return {minv, -(minv * a.o)};
    }
}

// Frame construction from axis.
inline frame3f make_frame_fromz(const vec3f& o, const vec3f& v) {
    // https://graphics.pixar.com/library/OrthonormalB/paper.pdf
    auto z    = normalize(v);
    auto sign = copysignf(1.0f, z.z);
    auto a    = -1.0f / (sign + z.z);
    auto b    = z.x * z.y * a;
    auto x    = vec3f{1.0f + sign * z.x * z.x * a, sign * b, -sign * z.x};
    auto y    = vec3f{b, sign + z.y * z.y * a, -z.y};
    return {x, y, z, o};
}
inline frame3f make_frame_fromzx(
    const vec3f& o, const vec3f& z_, const vec3f& x_) {
    auto z = normalize(z_);
    auto x = orthonormalize(x_, z);
    auto y = normalize(cross(z, x));
    return {x, y, z, o};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// QUATERNIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Quaternions to represent rotations
struct quat4f {
    float x, y, z, w;

    // constructors
    constexpr quat4f() : x{0}, y{0}, z{0}, w{1} {}
    constexpr quat4f(float x, float y, float z, float w)
        : x{x}, y{y}, z{z}, w{w} {}
};

// Constants
constexpr auto identity_quat4f = quat4f{0, 0, 0, 1};

// Quaternion operatons
inline quat4f operator+(const quat4f& a, const quat4f& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}
inline quat4f operator*(const quat4f& a, float b) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
}
inline quat4f operator/(const quat4f& a, float b) {
    return {a.x / b, a.y / b, a.z / b, a.w / b};
}
inline quat4f operator*(const quat4f& a, const quat4f& b) {
    return {a.x * b.w + a.w * b.x + a.y * b.w - a.z * b.y,
        a.y * b.w + a.w * b.y + a.z * b.x - a.x * b.z,
        a.z * b.w + a.w * b.z + a.x * b.y - a.y * b.x,
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z};
}

// Quaterion operations
inline float dot(const quat4f& a, const quat4f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}
inline float  length(const quat4f& a) { return sqrt(dot(a, a)); }
inline quat4f normalize(const quat4f& a) {
    auto l = length(a);
    return (l != 0) ? a / l : a;
}
inline quat4f conjugate(const quat4f& a) { return {-a.x, -a.y, -a.z, a.w}; }
inline quat4f inverse(const quat4f& a) { return conjugate(a) / dot(a, a); }
inline float  uangle(const quat4f& a, const quat4f& b) {
    auto d = dot(a, b);
    return d > 1 ? 0 : std::acos(d < -1 ? -1 : d);
}
inline quat4f lerp(const quat4f& a, const quat4f& b, float t) {
    return a * (1 - t) + b * t;
}
inline quat4f nlerp(const quat4f& a, const quat4f& b, float t) {
    return normalize(lerp(a, b, t));
}
inline quat4f slerp(const quat4f& a, const quat4f& b, float t) {
    auto th = uangle(a, b);
    return th == 0 ? a
                   : a * (sin(th * (1 - t)) / sin(th)) +
                         b * (sin(th * t) / sin(th));
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// AXIS ALIGNED BOUNDING BOXES
// -----------------------------------------------------------------------------
namespace yocto {

// Axis aligned bounding box represented as a min/max vector pairs.
struct bbox2f {
    vec2f min = {float_max, float_max};
    vec2f max = {float_min, float_min};

    constexpr bbox2f() {}
    constexpr bbox2f(const vec2f& min, const vec2f& max)
        : min{min}, max{max} {}

    constexpr vec2f&       operator[](int i) { return (&min)[i]; }
    constexpr const vec2f& operator[](int i) const { return (&min)[i]; }
};

// Axis aligned bounding box represented as a min/max vector pairs.
struct bbox3f {
    vec3f min = {float_max, float_max, float_max};
    vec3f max = {float_min, float_min, float_min};

    constexpr bbox3f() {}
    constexpr bbox3f(const vec3f& min, const vec3f& max)
        : min{min}, max{max} {}

    constexpr vec3f&       operator[](int i) { return (&min)[i]; }
    constexpr const vec3f& operator[](int i) const { return (&min)[i]; }
};

// Empty bbox constant.
constexpr auto invalid_bbox2f = bbox2f{};
constexpr auto invalid_bbox3f = bbox3f{};

// Bounding box properties
constexpr vec2f bbox_center(const bbox2f& a) {
    return (a.min + a.max) / 2;
}
constexpr vec2f bbox_size(const bbox2f& a) {
    return a.max - a.min;
}

// Bounding box comparisons.
constexpr bool operator==(const bbox2f& a, const bbox2f& b) {
    return a.min == b.min && a.max == b.max;
}
constexpr bool operator!=(const bbox2f& a, const bbox2f& b) {
    return a.min != b.min || a.max != b.max;
}

// Bounding box expansions with points and other boxes.
constexpr bbox2f operator+(const bbox2f& a, const vec2f& b) {
    return {min(a.min, b), max(a.max, b)};
}
constexpr bbox2f operator+(const bbox2f& a, const bbox2f& b) {
    return {min(a.min, b.min), max(a.max, b.max)};
}
constexpr bbox2f& operator+=(bbox2f& a, const vec2f& b) {
    return a = a + b;
}
constexpr bbox2f& operator+=(bbox2f& a, const bbox2f& b) {
    return a = a + b;
}

// Bounding box properties
constexpr vec3f bbox_center(const bbox3f& a) {
    return (a.min + a.max) / 2;
}
constexpr vec3f bbox_size(const bbox3f& a) {
    return a.max - a.min;
}

// Bounding box comparisons.
constexpr bool operator==(const bbox3f& a, const bbox3f& b) {
    return a.min == b.min && a.max == b.max;
}
constexpr bool operator!=(const bbox3f& a, const bbox3f& b) {
    return a.min != b.min || a.max != b.max;
}

// Bounding box expansions with points and other boxes.
constexpr bbox3f operator+(const bbox3f& a, const vec3f& b) {
    return {min(a.min, b), max(a.max, b)};
}
constexpr bbox3f operator+(const bbox3f& a, const bbox3f& b) {
    return {min(a.min, b.min), max(a.max, b.max)};
}
constexpr bbox3f& operator+=(bbox3f& a, const vec3f& b) {
    return a = a + b;
}
constexpr bbox3f& operator+=(bbox3f& a, const bbox3f& b) {
    return a = a + b;
}

// Primitive bounds.
constexpr bbox3f point_bounds(const vec3f& p) {
    auto a = bbox3f{};
    a += p;
    return a;
}
constexpr bbox3f point_bounds(const vec3f& p, float r) {
    auto a = bbox3f{};
    a += p - r;
    a += p + r;
    return a;
}
constexpr bbox3f line_bounds(const vec3f& p0, const vec3f& p1) {
    auto a = bbox3f{};
    a += p0;
    a += p1;
    return a;
}
constexpr bbox3f line_bounds(
    const vec3f& p0, const vec3f& p1, float r0, float r1) {
    auto a = bbox3f{};
    a += p0 - r0;
    a += p0 + r0;
    a += p1 - r1;
    a += p1 + r1;
    return a;
}
constexpr bbox3f triangle_bounds(
    const vec3f& p0, const vec3f& p1, const vec3f& p2) {
    auto a = bbox3f{};
    a += p0;
    a += p1;
    a += p2;
    return a;
}
constexpr bbox3f quad_bounds(const vec3f& p0, const vec3f& p1,
    const vec3f& p2, const vec3f& p3) {
    auto a = bbox3f{};
    a += p0;
    a += p1;
    a += p2;
    a += p3;
    return a;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// RAYS
// -----------------------------------------------------------------------------
namespace yocto {

// Ray esplison
constexpr auto ray_eps = 1e-4f;

struct ray2f {
    vec2f o = {0, 0};
    vec2f d = {0, 1};
    float tmin = ray_eps;
    float tmax = float_max;

    constexpr ray2f() {}
    constexpr ray2f(const vec2f& o, const vec2f& d, float tmin = ray_eps,
        float tmax = float_max)
        : o{o}, d{d}, tmin{tmin}, tmax{tmax} {}
};

// Rays with origin, direction and min/max t value.
struct ray3f {
    vec3f o = {0, 0, 0};
    vec3f d = {0, 0, 1};
    float tmin = ray_eps;
    float tmax = float_max;

    constexpr ray3f() {}
    constexpr ray3f(const vec3f& o, const vec3f& d, float tmin = ray_eps,
        float tmax = float_max)
        : o{o}, d{d}, tmin{tmin}, tmax{tmax} {}
};

}  // namespace yocto

// -----------------------------------------------------------------------------
// TRANSFORMS
// -----------------------------------------------------------------------------
namespace yocto {

// Transforms points, vectors and directions by matrices.
inline vec2f transform_point(const mat3f& a, const vec2f& b) {
    auto tvb = a * vec3f{b.x, b.y, 1};
    return vec2f{tvb.x, tvb.y} / tvb.z;
}
inline vec2f transform_vector(const mat3f& a, const vec2f& b) {
    auto tvb = a * vec3f{b.x, b.y, 0};
    return vec2f{tvb.x, tvb.y} / tvb.z;
}
inline vec2f transform_direction(const mat3f& a, const vec2f& b) {
    return normalize(transform_vector(a, b));
}
inline vec2f transform_normal(const mat3f& a, const vec2f& b) {
    return normalize(transform_vector(transpose(inverse(a)), b));
}
inline vec2f transform_vector(const mat2f& a, const vec2f& b) { return a * b; }
inline vec2f transform_direction(const mat2f& a, const vec2f& b) {
    return normalize(transform_vector(a, b));
}
inline vec2f transform_normal(const mat2f& a, const vec2f& b) {
    return normalize(transform_vector(transpose(inverse(a)), b));
}

inline vec3f transform_point(const mat4f& a, const vec3f& b) {
    auto tvb = a * vec4f{b.x, b.y, b.z, 1};
    return vec3f{tvb.x, tvb.y, tvb.z} / tvb.w;
}
inline vec3f transform_vector(const mat4f& a, const vec3f& b) {
    auto tvb = a * vec4f{b.x, b.y, b.z, 0};
    return vec3f{tvb.x, tvb.y, tvb.z};
}
inline vec3f transform_direction(const mat4f& a, const vec3f& b) {
    return normalize(transform_vector(a, b));
}
inline vec3f transform_vector(const mat3f& a, const vec3f& b) { return a * b; }
inline vec3f transform_direction(const mat3f& a, const vec3f& b) {
    return normalize(transform_vector(a, b));
}
inline vec3f transform_normal(const mat3f& a, const vec3f& b) {
    return normalize(transform_vector(transpose(inverse(a)), b));
}

// Transforms points, vectors and directions by frames.
inline vec2f transform_point(const frame2f& a, const vec2f& b) {
    return a.x * b.x + a.y * b.y + a.o;
}
inline vec2f transform_vector(const frame2f& a, const vec2f& b) {
    return a.x * b.x + a.y * b.y;
}
inline vec2f transform_direction(const frame2f& a, const vec2f& b) {
    return normalize(transform_vector(a, b));
}
inline vec2f transform_normal(
    const frame2f& a, const vec2f& b, bool non_rigid = false) {
    if (non_rigid) {
        return transform_normal(a.m(), b);
    } else {
        return normalize(transform_vector(a, b));
    }
}

// Transforms points, vectors and directions by frames.
inline vec3f transform_point(const frame3f& a, const vec3f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.o;
}
inline vec3f transform_vector(const frame3f& a, const vec3f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline vec3f transform_direction(const frame3f& a, const vec3f& b) {
    return normalize(transform_vector(a, b));
}
inline vec3f transform_normal(
    const frame3f& a, const vec3f& b, bool non_rigid = false) {
    if (non_rigid) {
        return transform_normal(a.m(), b);
    } else {
        return normalize(transform_vector(a, b));
    }
}

// Transforms rays and bounding boxes by matrices.
inline ray3f transform_ray(const mat4f& a, const ray3f& b) {
    return {transform_point(a, b.o), transform_vector(a, b.d), b.tmin, b.tmax};
}
inline ray3f transform_ray(const frame3f& a, const ray3f& b) {
    return {transform_point(a, b.o), transform_vector(a, b.d), b.tmin, b.tmax};
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

// Translation, scaling and rotations transforms.
inline frame3f make_translation_frame(const vec3f& a) {
    return {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, a};
}
inline frame3f make_scaling_frame(const vec3f& a) {
    return {{a.x, 0, 0}, {0, a.y, 0}, {0, 0, a.z}, {0, 0, 0}};
}
inline frame3f make_rotation_frame(const vec3f& axis, float angle) {
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
inline frame3f make_rotation_frame(const vec4f& quat) {
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
inline frame3f make_rotation_frame(const quat4f& quat) {
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
inline frame3f make_rotation_frame(const mat3f& rot) {
    return {rot.x, rot.y, rot.z, {0, 0, 0}};
}

// Lookat frame. Z-axis can be inverted with inv_xz.
inline frame3f make_lookat_frame(const vec3f& eye, const vec3f& center,
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
inline mat4f make_frustum_mat(
    float l, float r, float b, float t, float n, float f) {
    return {{2 * n / (r - l), 0, 0, 0}, {0, 2 * n / (t - b), 0, 0},
        {(r + l) / (r - l), (t + b) / (t - b), -(f + n) / (f - n), -1},
        {0, 0, -2 * f * n / (f - n), 0}};
}
inline mat4f make_ortho_mat(
    float l, float r, float b, float t, float n, float f) {
    return {{2 / (r - l), 0, 0, 0}, {0, 2 / (t - b), 0, 0},
        {0, 0, -2 / (f - n), 0},
        {-(r + l) / (r - l), -(t + b) / (t - b), -(f + n) / (f - n), 1}};
}
inline mat4f make_ortho2d_mat(
    float left, float right, float bottom, float top) {
    return make_ortho_mat(left, right, bottom, top, -1, 1);
}
inline mat4f make_ortho_mat(float xmag, float ymag, float near, float far) {
    return {{1 / xmag, 0, 0, 0}, {0, 1 / ymag, 0, 0},
        {0, 0, 2 / (near - far), 0}, {0, 0, (far + near) / (near - far), 1}};
}
inline mat4f make_perspective_mat(
    float fovy, float aspect, float near, float far) {
    auto tg = tan(fovy / 2);
    return {{1 / (aspect * tg), 0, 0, 0}, {0, 1 / tg, 0, 0},
        {0, 0, (far + near) / (near - far), -1},
        {0, 0, 2 * far * near / (near - far), 0}};
}
inline mat4f make_perspective_mat(float fovy, float aspect, float near) {
    auto tg = tan(fovy / 2);
    return {{1 / (aspect * tg), 0, 0, 0}, {0, 1 / tg, 0, 0}, {0, 0, -1, -1},
        {0, 0, 2 * near, 0}};
}

// Rotation conversions.
inline pair<vec3f, float> make_rotation_axisangle(const vec4f& quat) {
    return {normalize(vec3f{quat.x, quat.y, quat.z}), 2 * acos(quat.w)};
}
inline vec4f make_rotation_quat(const vec3f& axis, float angle) {
    auto len = length(axis);
    if (!len) return {0, 0, 0, 1};
    return vec4f{sin(angle / 2) * axis.x / len, sin(angle / 2) * axis.y / len,
        sin(angle / 2) * axis.z / len, cos(angle / 2)};
}
inline vec4f make_rotation_quat(const vec4f& axisangle) {
    return make_rotation_quat(
        vec3f{axisangle.x, axisangle.y, axisangle.z}, axisangle.w);
}

// Computes the image uv coordinates corresponding to the view parameters.
// Returns negative coordinates if out of the image.
inline vec2i get_image_coords(const vec2f& mouse_pos, const vec2f& center,
    float scale, const vec2i& txt_size) {
    auto xyf = (mouse_pos - center) / scale;
    return vec2i{(int)round(xyf.x + txt_size.x / 2.0f),
        (int)round(xyf.y + txt_size.y / 2.0f)};
}

// Center image and autofit.
inline void update_image_view(vec2f& center, float& scale, const vec2i& imsize,
    const vec2i& winsize, bool zoom_to_fit) {
    if (zoom_to_fit) {
        scale  = min(winsize.x / (float)imsize.x, winsize.y / (float)imsize.y);
        center = {(float)winsize.x / 2, (float)winsize.y / 2};
    } else {
        if (winsize.x >= imsize.x * scale) center.x = winsize.x / 2;
        if (winsize.y >= imsize.y * scale) center.y = winsize.y / 2;
    }
}

// Turntable for UI navigation.
inline void update_camera_turntable(vec3f& from, vec3f& to, vec3f& up,
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
inline void update_camera_turntable(frame3f& frame, float& focus,
    const vec2f& rotate, float dolly, const vec2f& pan) {
    // rotate if necessary
    if (rotate != zero2f) {
        auto phi   = atan2(frame.z.z, frame.z.x) + rotate.x;
        auto theta = acos(frame.z.y) + rotate.y;
        theta      = clamp(theta, 0.001f, pif - 0.001f);
        auto new_z = vec3f{
            sin(theta) * cos(phi), cos(theta), sin(theta) * sin(phi)};
        auto new_center = frame.o - frame.z * focus;
        auto new_o      = new_center + new_z * focus;
        frame           = make_lookat_frame(new_o, new_center, {0, 1, 0});
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
inline void update_camera_first_person(
    frame3f& frame, const vec3f& transl, const vec2f& rotate) {
    // https://gamedev.stackexchange.com/questions/30644/how-to-keep-my-quaternion-using-fps-camera-from-tilting-and-messing-up
    auto y = vec3f{0, 1, 0};
    auto z = orthonormalize(frame.z, y);
    auto x = cross(y, z);

    auto rot = make_rotation_frame(vec3f{1, 0, 0}, rotate.y) *
               yocto::frame3f{frame.x, frame.y, frame.z, vec3f{0, 0, 0}} *
               make_rotation_frame(vec3f{0, 1, 0}, rotate.x);
    auto pos = frame.o + transl.x * x + transl.y * y + transl.z * z;

    frame = {rot.x, rot.y, rot.z, pos};
}

}  // namespace yocto

#endif
