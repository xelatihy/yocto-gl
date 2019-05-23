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
// coordinates in float and int coordinates (`vec1f`, `vec2f`, `vec3f`, `vec4f`,
// `vec1i`, `vec2i`, `vec3i`, `vec4i`).
//
// We support 2-4 dimensional matrices (`mat2f`, `mat3f`, `mat4f`) with
// matrix-matrix and matrix-vector products, transposes and inverses.
// Matrices are stored in column-major order and are accessed and
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
#include <memory>
#include <numeric>
#include <string>
#include <string_view>
#include <tuple>
#include <unordered_map>
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

using std::array;
using std::exception;
using std::function;
using std::invalid_argument;
using std::make_shared;
using std::make_unique;
using std::out_of_range;
using std::pair;
using std::runtime_error;
using std::shared_ptr;
using std::string;
using std::string_view;
using std::tuple;
using std::unique_ptr;
using std::unordered_map;
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

constexpr auto int_max   = INT_MAX;
constexpr auto int_min   = INT_MIN;
constexpr auto float_max = FLT_MAX;
constexpr auto float_min = -FLT_MAX;

constexpr float min(float x, float y) { return (x < y) ? x : y; }
constexpr float max(float x, float y) { return (x > y) ? x : y; }
constexpr float clamp(float x, float min_, float max_) {
    return min(max(x, min_), max_);
}
constexpr float clamp01(float x) { return min(max(x, 0.0f), 1.0f); }
constexpr float lerp(float a, float b, float u) { return a * (1 - u) + b * u; }
constexpr float bias(float a, float bias) {
    return a / ((1 / bias - 2) * (1 - a) + 1);
}
constexpr float gain(float a, float gain) {
    if (a < 0.5f) {
        return bias(a * 2, gain) / 2;
    } else {
        return bias(a * 2 - 1, 1 - gain) / 2 + 0.5f;
    }
}
inline float radians(float x) { return x * pif / 180; }
inline float degrees(float x) { return x * 180 / pif; }

constexpr int min(int x, int y) { return (x < y) ? x : y; }
constexpr int max(int x, int y) { return (x > y) ? x : y; }
constexpr int clamp(int x, int min_, int max_) {
    return min(max(x, min_), max_);
}
constexpr int pow2(int x) { return 1 << x; }

}  // namespace yocto

// -----------------------------------------------------------------------------
// VECTORS
// -----------------------------------------------------------------------------
namespace yocto {

struct vec2f {
    float x = 0;
    float y = 0;

    constexpr vec2f() {}
    constexpr vec2f(float x, float y) : x{x}, y{y} {}
    constexpr explicit vec2f(float v) : x{v}, y{v} {}

    constexpr float&       operator[](int i) { return (&x)[i]; }
    constexpr const float& operator[](int i) const { return (&x)[i]; }
};

struct vec3f {
    float x = 0;
    float y = 0;
    float z = 0;

    constexpr vec3f() {}
    constexpr vec3f(float x, float y, float z) : x{x}, y{y}, z{z} {}
    constexpr vec3f(const vec2f& v, float z) : x{v.x}, y{v.y}, z{z} {}
    constexpr explicit vec3f(float v) : x{v}, y{v}, z{v} {}

    constexpr float&       operator[](int i) { return (&x)[i]; }
    constexpr const float& operator[](int i) const { return (&x)[i]; }
};

struct vec4f {
    float x = 0;
    float y = 0;
    float z = 0;
    float w = 0;

    constexpr vec4f() {}
    constexpr vec4f(float x, float y, float z, float w)
        : x{x}, y{y}, z{z}, w{w} {}
    constexpr vec4f(const vec3f& v, float w) : x{v.x}, y{v.y}, z{v.z}, w{w} {}
    constexpr explicit vec4f(float v) : x{v}, y{v}, z{v}, w{v} {}

    constexpr float&       operator[](int i) { return (&x)[i]; }
    constexpr const float& operator[](int i) const { return (&x)[i]; }
};

// Zero vector constants.
constexpr auto zero2f = vec2f{0, 0};
constexpr auto zero3f = vec3f{0, 0, 0};
constexpr auto zero4f = vec4f{0, 0, 0, 0};

// Element access
inline vec3f&       xyz(vec4f& a) { return (vec3f&)a; }
inline const vec3f& xyz(const vec4f& a) { return (const vec3f&)a; }

// Vector comparison operations.
constexpr bool operator==(const vec2f& a, const vec2f& b) {
    return a.x == b.x && a.y == b.y;
}
constexpr bool operator!=(const vec2f& a, const vec2f& b) {
    return a.x != b.x || a.y != b.y;
}

// Vector operations.
constexpr vec2f operator+(const vec2f& a) { return a; }
constexpr vec2f operator-(const vec2f& a) { return {-a.x, -a.y}; }
constexpr vec2f operator+(const vec2f& a, const vec2f& b) {
    return {a.x + b.x, a.y + b.y};
}
constexpr vec2f operator+(const vec2f& a, float b) {
    return {a.x + b, a.y + b};
}
constexpr vec2f operator+(float a, const vec2f& b) {
    return {a + b.x, a + b.y};
}
constexpr vec2f operator-(const vec2f& a, const vec2f& b) {
    return {a.x - b.x, a.y - b.y};
}
constexpr vec2f operator-(const vec2f& a, float b) {
    return {a.x - b, a.y - b};
}
constexpr vec2f operator-(float a, const vec2f& b) {
    return {a - b.x, a - b.y};
}
constexpr vec2f operator*(const vec2f& a, const vec2f& b) {
    return {a.x * b.x, a.y * b.y};
}
constexpr vec2f operator*(const vec2f& a, float b) {
    return {a.x * b, a.y * b};
}
constexpr vec2f operator*(float a, const vec2f& b) {
    return {a * b.x, a * b.y};
}
constexpr vec2f operator/(const vec2f& a, const vec2f& b) {
    return {a.x / b.x, a.y / b.y};
}
constexpr vec2f operator/(const vec2f& a, float b) {
    return {a.x / b, a.y / b};
}
constexpr vec2f operator/(float a, const vec2f& b) {
    return {a / b.x, a / b.y};
}

// Vector assignments
constexpr vec2f& operator+=(vec2f& a, const vec2f& b) { return a = a + b; }
constexpr vec2f& operator+=(vec2f& a, float b) { return a = a + b; }
constexpr vec2f& operator-=(vec2f& a, const vec2f& b) { return a = a - b; }
constexpr vec2f& operator-=(vec2f& a, float b) { return a = a - b; }
constexpr vec2f& operator*=(vec2f& a, const vec2f& b) { return a = a * b; }
constexpr vec2f& operator*=(vec2f& a, float b) { return a = a * b; }
constexpr vec2f& operator/=(vec2f& a, const vec2f& b) { return a = a / b; }
constexpr vec2f& operator/=(vec2f& a, float b) { return a = a / b; }

// Vector products and lengths.
constexpr float dot(const vec2f& a, const vec2f& b) {
    return a.x * b.x + a.y * b.y;
}
constexpr float cross(const vec2f& a, const vec2f& b) {
    return a.x * b.y - a.y * b.x;
}

inline float length(const vec2f& a) { return sqrt(dot(a, a)); }
inline vec2f normalize(const vec2f& a) {
    auto l = length(a);
    return (l != 0) ? a / l : a;
}
inline float distance(const vec2f& a, const vec2f& b) { return length(a - b); }
inline float distance_squared(const vec2f& a, const vec2f& b) {
    return dot(a - b, a - b);
}

// Max element and clamp.
constexpr vec2f max(const vec2f& a, float b) {
    return {max(a.x, b), max(a.y, b)};
}
constexpr vec2f min(const vec2f& a, float b) {
    return {min(a.x, b), min(a.y, b)};
}
constexpr vec2f max(const vec2f& a, const vec2f& b) {
    return {max(a.x, b.x), max(a.y, b.y)};
}
constexpr vec2f min(const vec2f& a, const vec2f& b) {
    return {min(a.x, b.x), min(a.y, b.y)};
}
constexpr vec2f clamp(const vec2f& x, float min, float max) {
    return {clamp(x.x, min, max), clamp(x.y, min, max)};
}
constexpr vec2f clamp01(const vec2f& x) { return {clamp01(x.x), clamp01(x.y)}; }
inline vec2f    lerp(const vec2f& a, const vec2f& b, float u) {
    return a * (1 - u) + b * u;
}
inline vec2f lerp(const vec2f& a, const vec2f& b, const vec2f& u) {
    return a * (1 - u) + b * u;
}

constexpr float max(const vec2f& a) { return max(a.x, a.y); }
constexpr float min(const vec2f& a) { return min(a.x, a.y); }
constexpr float sum(const vec2f& a) { return a.x + a.y; }
constexpr float mean(const vec2f& a) { return sum(a) / 2; }
constexpr int   max_element(const vec2f& a) {
    auto pos = 0;
    for (auto i = 1; i < 2; i++) {
        if (a[i] > a[pos]) pos = i;
    }
    return pos;
}
constexpr int min_element(const vec2f& a) {
    auto pos = 0;
    for (auto i = 1; i < 2; i++) {
        if (a[i] < a[pos]) pos = i;
    }
    return pos;
}

// Functions applied to vector elements
inline vec2f abs(const vec2f& a) { return {abs(a.x), abs(a.y)}; };
inline vec2f sqrt(const vec2f& a) { return {sqrt(a.x), sqrt(a.y)}; };
inline vec2f exp(const vec2f& a) { return {exp(a.x), exp(a.y)}; };
inline vec2f log(const vec2f& a) { return {log(a.x), log(a.y)}; };
inline vec2f exp2(const vec2f& a) { return {exp2(a.x), exp2(a.y)}; };
inline vec2f log2(const vec2f& a) { return {log2(a.x), log2(a.y)}; };
inline bool isfinite(const vec2f& a) { return isfinite(a.x) && isfinite(a.y); };
inline vec2f pow(const vec2f& a, float b) {
    return {pow(a.x, b), pow(a.y, b)};
};
inline vec2f pow(const vec2f& a, const vec2f& b) {
    return {pow(a.x, b.x), pow(a.y, b.y)};
};
inline vec2f gain(const vec2f& a, float b) {
    return {gain(a.x, b), gain(a.y, b)};
};

// Vector comparison operations.
constexpr bool operator==(const vec3f& a, const vec3f& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}
constexpr bool operator!=(const vec3f& a, const vec3f& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z;
}

// Vector operations.
constexpr vec3f operator+(const vec3f& a) { return a; }
constexpr vec3f operator-(const vec3f& a) { return {-a.x, -a.y, -a.z}; }
constexpr vec3f operator+(const vec3f& a, const vec3f& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
constexpr vec3f operator+(const vec3f& a, float b) {
    return {a.x + b, a.y + b, a.z + b};
}
constexpr vec3f operator+(float a, const vec3f& b) {
    return {a + b.x, a + b.y, a + b.z};
}
constexpr vec3f operator-(const vec3f& a, const vec3f& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
constexpr vec3f operator-(const vec3f& a, float b) {
    return {a.x - b, a.y - b, a.z - b};
}
constexpr vec3f operator-(float a, const vec3f& b) {
    return {a - b.x, a - b.y, a - b.z};
}
constexpr vec3f operator*(const vec3f& a, const vec3f& b) {
    return {a.x * b.x, a.y * b.y, a.z * b.z};
}
constexpr vec3f operator*(const vec3f& a, float b) {
    return {a.x * b, a.y * b, a.z * b};
}
constexpr vec3f operator*(float a, const vec3f& b) {
    return {a * b.x, a * b.y, a * b.z};
}
constexpr vec3f operator/(const vec3f& a, const vec3f& b) {
    return {a.x / b.x, a.y / b.y, a.z / b.z};
}
constexpr vec3f operator/(const vec3f& a, float b) {
    return {a.x / b, a.y / b, a.z / b};
}
constexpr vec3f operator/(float a, const vec3f& b) {
    return {a / b.x, a / b.y, a / b.z};
}

// Vector assignments
constexpr vec3f& operator+=(vec3f& a, const vec3f& b) { return a = a + b; }
constexpr vec3f& operator+=(vec3f& a, float b) { return a = a + b; }
constexpr vec3f& operator-=(vec3f& a, const vec3f& b) { return a = a - b; }
constexpr vec3f& operator-=(vec3f& a, float b) { return a = a - b; }
constexpr vec3f& operator*=(vec3f& a, const vec3f& b) { return a = a * b; }
constexpr vec3f& operator*=(vec3f& a, float b) { return a = a * b; }
constexpr vec3f& operator/=(vec3f& a, const vec3f& b) { return a = a / b; }
constexpr vec3f& operator/=(vec3f& a, float b) { return a = a / b; }

// Vector products and lengths.
constexpr float dot(const vec3f& a, const vec3f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
constexpr vec3f cross(const vec3f& a, const vec3f& b) {
    return {
        a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

inline float length(const vec3f& a) { return sqrt(dot(a, a)); }
inline vec3f normalize(const vec3f& a) {
    auto l = length(a);
    return (l != 0) ? a / l : a;
}
inline float distance(const vec3f& a, const vec3f& b) { return length(a - b); }
inline float distance_squared(const vec3f& a, const vec3f& b) {
    return dot(a - b, a - b);
}

inline float angle(const vec3f& a, const vec3f& b) {
    return acos(clamp(dot(normalize(a), normalize(b)), (float)-1, (float)1));
}

// Orthogonal vectors.
inline vec3f orthogonal(const vec3f& v) {
    // http://lolengine.net/blog/2013/09/21/picking-orthogonal-vector-combing-coconuts)
    return abs(v.x) > abs(v.z) ? vec3f{-v.y, v.x, 0} : vec3f{0, -v.z, v.y};
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
    auto k = 1 - eta * eta * max((float)0, 1 - dot(n, w) * dot(n, w));
    if (k < 0) return {0, 0, 0};  // tir
    return -w * eta + (eta * dot(n, w) - sqrt(k)) * n;
}

// Max element and clamp.
constexpr vec3f max(const vec3f& a, float b) {
    return {max(a.x, b), max(a.y, b), max(a.z, b)};
}
constexpr vec3f min(const vec3f& a, float b) {
    return {min(a.x, b), min(a.y, b), min(a.z, b)};
}
constexpr vec3f max(const vec3f& a, const vec3f& b) {
    return {max(a.x, b.x), max(a.y, b.y), max(a.z, b.z)};
}
constexpr vec3f min(const vec3f& a, const vec3f& b) {
    return {min(a.x, b.x), min(a.y, b.y), min(a.z, b.z)};
}
constexpr vec3f clamp(const vec3f& x, float min, float max) {
    return {clamp(x.x, min, max), clamp(x.y, min, max), clamp(x.z, min, max)};
}
constexpr vec3f clamp01(const vec3f& x) {
    return {clamp01(x.x), clamp01(x.y), clamp01(x.z)};
}
inline vec3f lerp(const vec3f& a, const vec3f& b, float u) {
    return a * (1 - u) + b * u;
}
inline vec3f lerp(const vec3f& a, const vec3f& b, const vec3f& u) {
    return a * (1 - u) + b * u;
}

constexpr float max(const vec3f& a) { return max(max(a.x, a.y), a.z); }
constexpr float min(const vec3f& a) { return min(min(a.x, a.y), a.z); }
constexpr float sum(const vec3f& a) { return a.x + a.y + a.z; }
inline float    mean(const vec3f& a) { return sum(a) / 3; }
constexpr int   max_element(const vec3f& a) {
    auto pos = 0;
    for (auto i = 1; i < 3; i++) {
        if (a[i] > a[pos]) pos = i;
    }
    return pos;
}
constexpr int min_element(const vec3f& a) {
    auto pos = 0;
    for (auto i = 1; i < 3; i++) {
        if (a[i] < a[pos]) pos = i;
    }
    return pos;
}

// Functions applied to vector elements
inline vec3f abs(const vec3f& a) { return {abs(a.x), abs(a.y), abs(a.z)}; };
inline vec3f sqrt(const vec3f& a) { return {sqrt(a.x), sqrt(a.y), sqrt(a.z)}; };
inline vec3f exp(const vec3f& a) { return {exp(a.x), exp(a.y), exp(a.z)}; };
inline vec3f log(const vec3f& a) { return {log(a.x), log(a.y), log(a.z)}; };
inline vec3f exp2(const vec3f& a) { return {exp2(a.x), exp2(a.y), exp2(a.z)}; };
inline vec3f log2(const vec3f& a) { return {log2(a.x), log2(a.y), log2(a.z)}; };
inline vec3f pow(const vec3f& a, float b) {
    return {pow(a.x, b), pow(a.y, b), pow(a.z, b)};
};
inline vec3f pow(const vec3f& a, const vec3f& b) {
    return {pow(a.x, b.x), pow(a.y, b.y), pow(a.z, b.z)};
};
inline vec3f gain(const vec3f& a, float b) {
    return {gain(a.x, b), gain(a.y, b), gain(a.z, b)};
};
inline bool isfinite(const vec3f& a) {
    return isfinite(a.x) && isfinite(a.y) && isfinite(a.z);
};

// Vector comparison operations.
constexpr bool operator==(const vec4f& a, const vec4f& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}
constexpr bool operator!=(const vec4f& a, const vec4f& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w;
}

// Vector operations.
constexpr vec4f operator+(const vec4f& a) { return a; }
constexpr vec4f operator-(const vec4f& a) { return {-a.x, -a.y, -a.z, -a.w}; }
constexpr vec4f operator+(const vec4f& a, const vec4f& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}
constexpr vec4f operator+(const vec4f& a, float b) {
    return {a.x + b, a.y + b, a.z + b, a.w + b};
}
constexpr vec4f operator+(float a, const vec4f& b) {
    return {a + b.x, a + b.y, a + b.z, a + b.w};
}
constexpr vec4f operator-(const vec4f& a, const vec4f& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w};
}
constexpr vec4f operator-(const vec4f& a, float b) {
    return {a.x - b, a.y - b, a.z - b, a.w - b};
}
constexpr vec4f operator-(float a, const vec4f& b) {
    return {a - b.x, a - b.y, a - b.z, a - b.w};
}
constexpr vec4f operator*(const vec4f& a, const vec4f& b) {
    return {a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w};
}
constexpr vec4f operator*(const vec4f& a, float b) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
}
constexpr vec4f operator*(float a, const vec4f& b) {
    return {a * b.x, a * b.y, a * b.z, a * b.w};
}
constexpr vec4f operator/(const vec4f& a, const vec4f& b) {
    return {a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w};
}
constexpr vec4f operator/(const vec4f& a, float b) {
    return {a.x / b, a.y / b, a.z / b, a.w / b};
}
constexpr vec4f operator/(float a, const vec4f& b) {
    return {a / b.x, a / b.y, a / b.z, a / b.w};
}

// Vector assignments
constexpr vec4f& operator+=(vec4f& a, const vec4f& b) { return a = a + b; }
constexpr vec4f& operator+=(vec4f& a, float b) { return a = a + b; }
constexpr vec4f& operator-=(vec4f& a, const vec4f& b) { return a = a - b; }
constexpr vec4f& operator-=(vec4f& a, float b) { return a = a - b; }
constexpr vec4f& operator*=(vec4f& a, const vec4f& b) { return a = a * b; }
constexpr vec4f& operator*=(vec4f& a, float b) { return a = a * b; }
constexpr vec4f& operator/=(vec4f& a, const vec4f& b) { return a = a / b; }
constexpr vec4f& operator/=(vec4f& a, float b) { return a = a / b; }

// Vector products and lengths.
constexpr float dot(const vec4f& a, const vec4f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}
inline float length(const vec4f& a) { return sqrt(dot(a, a)); }
inline vec4f normalize(const vec4f& a) {
    auto l = length(a);
    return (l != 0) ? a / l : a;
}
inline float distance(const vec4f& a, const vec4f& b) { return length(a - b); }
inline float distance_squared(const vec4f& a, const vec4f& b) {
    return dot(a - b, a - b);
}

inline vec4f slerp(const vec4f& a, const vec4f& b, float u) {
    // https://en.wikipedia.org/wiki/Slerp
    auto an = normalize(a), bn = normalize(b);
    auto d = dot(an, bn);
    if (d < 0) {
        bn = -bn;
        d  = -d;
    }
    if (d > (float)0.9995) return normalize(an + u * (bn - an));
    auto th = acos(clamp(d, (float)-1, (float)1));
    if (!th) return an;
    return an * (sin(th * (1 - u)) / sin(th)) + bn * (sin(th * u) / sin(th));
}

// Max element and clamp.
constexpr vec4f max(const vec4f& a, float b) {
    return {max(a.x, b), max(a.y, b), max(a.z, b), max(a.w, b)};
}
constexpr vec4f min(const vec4f& a, float b) {
    return {min(a.x, b), min(a.y, b), min(a.z, b), min(a.w, b)};
}
constexpr vec4f max(const vec4f& a, const vec4f& b) {
    return {max(a.x, b.x), max(a.y, b.y), max(a.z, b.z), max(a.w, b.w)};
}
constexpr vec4f min(const vec4f& a, const vec4f& b) {
    return {min(a.x, b.x), min(a.y, b.y), min(a.z, b.z), min(a.w, b.w)};
}
constexpr vec4f clamp(const vec4f& x, float min, float max) {
    return {clamp(x.x, min, max), clamp(x.y, min, max), clamp(x.z, min, max),
        clamp(x.w, min, max)};
}
constexpr vec4f clamp01(const vec4f& x) {
    return {clamp01(x.x), clamp01(x.y), clamp01(x.z), clamp01(x.w)};
}
inline vec4f lerp(const vec4f& a, const vec4f& b, float u) {
    return a * (1 - u) + b * u;
}
inline vec4f lerp(const vec4f& a, const vec4f& b, const vec4f& u) {
    return a * (1 - u) + b * u;
}

constexpr float max(const vec4f& a) {
    return max(max(max(a.x, a.y), a.z), a.w);
}
constexpr float min(const vec4f& a) {
    return min(min(min(a.x, a.y), a.z), a.w);
}
constexpr float sum(const vec4f& a) { return a.x + a.y + a.z + a.w; }
constexpr float mean(const vec4f& a) { return sum(a) / 4; }
constexpr int   max_element(const vec4f& a) {
    auto pos = 0;
    for (auto i = 1; i < 4; i++) {
        if (a[i] > a[pos]) pos = i;
    }
    return pos;
}
constexpr int min_element(const vec4f& a) {
    auto pos = 0;
    for (auto i = 1; i < 4; i++) {
        if (a[i] < a[pos]) pos = i;
    }
    return pos;
}

// Functions applied to vector elements
inline vec4f abs(const vec4f& a) {
    return {abs(a.x), abs(a.y), abs(a.z), abs(a.w)};
};
inline vec4f sqrt(const vec4f& a) {
    return {sqrt(a.x), sqrt(a.y), sqrt(a.z), sqrt(a.w)};
};
inline vec4f exp(const vec4f& a) {
    return {exp(a.x), exp(a.y), exp(a.z), exp(a.w)};
};
inline vec4f log(const vec4f& a) {
    return {log(a.x), log(a.y), log(a.z), log(a.w)};
};
inline vec4f exp2(const vec4f& a) {
    return {exp2(a.x), exp2(a.y), exp2(a.z), exp2(a.w)};
};
inline vec4f log2(const vec4f& a) {
    return {log2(a.x), log2(a.y), log2(a.z), log2(a.w)};
};
inline vec4f pow(const vec4f& a, float b) {
    return {pow(a.x, b), pow(a.y, b), pow(a.z, b), pow(a.w, b)};
};
inline vec4f pow(const vec4f& a, const vec4f& b) {
    return {pow(a.x, b.x), pow(a.y, b.y), pow(a.z, b.z), pow(a.w, b.w)};
};
inline vec4f gain(const vec4f& a, float b) {
    return {gain(a.x, b), gain(a.y, b), gain(a.z, b), gain(a.w, b)};
};
inline bool isfinite(const vec4f& a) {
    return isfinite(a.x) && isfinite(a.y) && isfinite(a.z) && isfinite(a.w);
};

// Quaternion operatons represented as xi + yj + zk + w
// const auto identity_quat4f = vec4f{0, 0, 0, 1};
constexpr vec4f quat_mul(const vec4f& a, float b) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
}
constexpr vec4f quat_mul(const vec4f& a, const vec4f& b) {
    return {a.x * b.w + a.w * b.x + a.y * b.w - a.z * b.y,
        a.y * b.w + a.w * b.y + a.z * b.x - a.x * b.z,
        a.z * b.w + a.w * b.z + a.x * b.y - a.y * b.x,
        a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z};
}
constexpr vec4f quat_conjugate(const vec4f& a) {
    return {-a.x, -a.y, -a.z, a.w};
}
constexpr vec4f quat_inverse(const vec4f& a) {
    return quat_conjugate(a) / dot(a, a);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// INTEGER VECTORS
// -----------------------------------------------------------------------------
namespace yocto {

struct vec2i {
    int x = 0;
    int y = 0;

    constexpr vec2i() {}
    constexpr vec2i(int x, int y) : x{x}, y{y} {}
    constexpr explicit vec2i(int v) : x{v}, y{v} {}

    constexpr int&       operator[](int i) { return (&x)[i]; }
    constexpr const int& operator[](int i) const { return (&x)[i]; }
};

struct vec3i {
    int x = 0;
    int y = 0;
    int z = 0;

    constexpr vec3i() {}
    constexpr vec3i(int x, int y, int z) : x{x}, y{y}, z{z} {}
    constexpr vec3i(const vec2i& v, int z) : x{v.x}, y{v.y}, z{z} {}
    constexpr explicit vec3i(int v) : x{v}, y{v}, z{v} {}

    constexpr int&       operator[](int i) { return (&x)[i]; }
    constexpr const int& operator[](int i) const { return (&x)[i]; }
};

struct vec4i {
    int x = 0;
    int y = 0;
    int z = 0;
    int w = 0;

    constexpr vec4i() {}
    constexpr vec4i(int x, int y, int z, int w) : x{x}, y{y}, z{z}, w{w} {}
    constexpr vec4i(const vec3i& v, int w) : x{v.x}, y{v.y}, z{v.z}, w{w} {}
    constexpr explicit vec4i(int v) : x{v}, y{v}, z{v}, w{v} {}

    constexpr int&       operator[](int i) { return (&x)[i]; }
    constexpr const int& operator[](int i) const { return (&x)[i]; }
};

struct vec4b {
    byte x = 0;
    byte y = 0;
    byte z = 0;
    byte w = 0;

    constexpr vec4b() {}
    constexpr vec4b(byte x, byte y, byte z, byte w) : x{x}, y{y}, z{z}, w{w} {}
    constexpr explicit vec4b(byte v) : x{v}, y{v}, z{v}, w{v} {}

    constexpr byte&       operator[](int i) { return (&x)[i]; }
    constexpr const byte& operator[](int i) const { return (&x)[i]; }
};

// Zero vector constants.
constexpr auto zero2i = vec2i{0, 0};
constexpr auto zero3i = vec3i{0, 0, 0};
constexpr auto zero4i = vec4i{0, 0, 0, 0};
constexpr auto zero4b = vec4b{0, 0, 0, 0};

// Element access
inline vec3i&       xyz(vec4i& a) { return (vec3i&)a; }
inline const vec3i& xyz(const vec4i& a) { return (const vec3i&)a; }

// Vector comparison operations.
constexpr bool operator==(const vec2i& a, const vec2i& b) {
    return a.x == b.x && a.y == b.y;
}
constexpr bool operator!=(const vec2i& a, const vec2i& b) {
    return a.x != b.x || a.y != b.y;
}

// Vector operations.
constexpr vec2i operator+(const vec2i& a) { return a; }
constexpr vec2i operator-(const vec2i& a) { return {-a.x, -a.y}; }
constexpr vec2i operator+(const vec2i& a, const vec2i& b) {
    return {a.x + b.x, a.y + b.y};
}
constexpr vec2i operator+(const vec2i& a, int b) { return {a.x + b, a.y + b}; }
constexpr vec2i operator+(int a, const vec2i& b) { return {a + b.x, a + b.y}; }
constexpr vec2i operator-(const vec2i& a, const vec2i& b) {
    return {a.x - b.x, a.y - b.y};
}
constexpr vec2i operator-(const vec2i& a, int b) { return {a.x - b, a.y - b}; }
constexpr vec2i operator-(int a, const vec2i& b) { return {a - b.x, a - b.y}; }
constexpr vec2i operator*(const vec2i& a, const vec2i& b) {
    return {a.x * b.x, a.y * b.y};
}
constexpr vec2i operator*(const vec2i& a, int b) { return {a.x * b, a.y * b}; }
constexpr vec2i operator*(int a, const vec2i& b) { return {a * b.x, a * b.y}; }
constexpr vec2i operator/(const vec2i& a, const vec2i& b) {
    return {a.x / b.x, a.y / b.y};
}
constexpr vec2i operator/(const vec2i& a, int b) { return {a.x / b, a.y / b}; }
constexpr vec2i operator/(int a, const vec2i& b) { return {a / b.x, a / b.y}; }

// Vector assignments
constexpr vec2i& operator+=(vec2i& a, const vec2i& b) { return a = a + b; }
constexpr vec2i& operator+=(vec2i& a, int b) { return a = a + b; }
constexpr vec2i& operator-=(vec2i& a, const vec2i& b) { return a = a - b; }
constexpr vec2i& operator-=(vec2i& a, int b) { return a = a - b; }
constexpr vec2i& operator*=(vec2i& a, const vec2i& b) { return a = a * b; }
constexpr vec2i& operator*=(vec2i& a, int b) { return a = a * b; }
constexpr vec2i& operator/=(vec2i& a, const vec2i& b) { return a = a / b; }
constexpr vec2i& operator/=(vec2i& a, int b) { return a = a / b; }

// Max element and clamp.
constexpr vec2i max(const vec2i& a, int b) {
    return {max(a.x, b), max(a.y, b)};
}
constexpr vec2i min(const vec2i& a, int b) {
    return {min(a.x, b), min(a.y, b)};
}
constexpr vec2i max(const vec2i& a, const vec2i& b) {
    return {max(a.x, b.x), max(a.y, b.y)};
}
constexpr vec2i min(const vec2i& a, const vec2i& b) {
    return {min(a.x, b.x), min(a.y, b.y)};
}
constexpr vec2i clamp(const vec2i& x, int min, int max) {
    return {clamp(x.x, min, max), clamp(x.y, min, max)};
}

constexpr int max(const vec2i& a) { return max(a.x, a.y); }
constexpr int min(const vec2i& a) { return min(a.x, a.y); }
constexpr int sum(const vec2i& a) { return a.x + a.y; }
constexpr int max_element(const vec2i& a) {
    auto pos = 0;
    for (auto i = 1; i < 2; i++) {
        if (a[i] > a[pos]) pos = i;
    }
    return pos;
}
constexpr int min_element(const vec2i& a) {
    auto pos = 0;
    for (auto i = 1; i < 2; i++) {
        if (a[i] < a[pos]) pos = i;
    }
    return pos;
}

// Functions applied to vector elements
inline vec2i abs(const vec2i& a) { return {abs(a.x), abs(a.y)}; };

// Vector comparison operations.
constexpr bool operator==(const vec3i& a, const vec3i& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}
constexpr bool operator!=(const vec3i& a, const vec3i& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z;
}

// Vector operations.
constexpr vec3i operator+(const vec3i& a) { return a; }
constexpr vec3i operator-(const vec3i& a) { return {-a.x, -a.y, -a.z}; }
constexpr vec3i operator+(const vec3i& a, const vec3i& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
constexpr vec3i operator+(const vec3i& a, int b) {
    return {a.x + b, a.y + b, a.z + b};
}
constexpr vec3i operator+(int a, const vec3i& b) {
    return {a + b.x, a + b.y, a + b.z};
}
constexpr vec3i operator-(const vec3i& a, const vec3i& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
constexpr vec3i operator-(const vec3i& a, int b) {
    return {a.x - b, a.y - b, a.z - b};
}
constexpr vec3i operator-(int a, const vec3i& b) {
    return {a - b.x, a - b.y, a - b.z};
}
constexpr vec3i operator*(const vec3i& a, const vec3i& b) {
    return {a.x * b.x, a.y * b.y, a.z * b.z};
}
constexpr vec3i operator*(const vec3i& a, int b) {
    return {a.x * b, a.y * b, a.z * b};
}
constexpr vec3i operator*(int a, const vec3i& b) {
    return {a * b.x, a * b.y, a * b.z};
}
constexpr vec3i operator/(const vec3i& a, const vec3i& b) {
    return {a.x / b.x, a.y / b.y, a.z / b.z};
}
constexpr vec3i operator/(const vec3i& a, int b) {
    return {a.x / b, a.y / b, a.z / b};
}
constexpr vec3i operator/(int a, const vec3i& b) {
    return {a / b.x, a / b.y, a / b.z};
}

// Vector assignments
constexpr vec3i& operator+=(vec3i& a, const vec3i& b) { return a = a + b; }
constexpr vec3i& operator+=(vec3i& a, int b) { return a = a + b; }
constexpr vec3i& operator-=(vec3i& a, const vec3i& b) { return a = a - b; }
constexpr vec3i& operator-=(vec3i& a, int b) { return a = a - b; }
constexpr vec3i& operator*=(vec3i& a, const vec3i& b) { return a = a * b; }
constexpr vec3i& operator*=(vec3i& a, int b) { return a = a * b; }
constexpr vec3i& operator/=(vec3i& a, const vec3i& b) { return a = a / b; }
constexpr vec3i& operator/=(vec3i& a, int b) { return a = a / b; }

// Max element and clamp.
constexpr vec3i max(const vec3i& a, int b) {
    return {max(a.x, b), max(a.y, b), max(a.z, b)};
}
constexpr vec3i min(const vec3i& a, int b) {
    return {min(a.x, b), min(a.y, b), min(a.z, b)};
}
constexpr vec3i max(const vec3i& a, const vec3i& b) {
    return {max(a.x, b.x), max(a.y, b.y), max(a.z, b.z)};
}
constexpr vec3i min(const vec3i& a, const vec3i& b) {
    return {min(a.x, b.x), min(a.y, b.y), min(a.z, b.z)};
}
constexpr vec3i clamp(const vec3i& x, int min, int max) {
    return {clamp(x.x, min, max), clamp(x.y, min, max), clamp(x.z, min, max)};
}

constexpr int max(const vec3i& a) { return max(max(a.x, a.y), a.z); }
constexpr int min(const vec3i& a) { return min(min(a.x, a.y), a.z); }
constexpr int sum(const vec3i& a) { return a.x + a.y + a.z; }
constexpr int max_element(const vec3i& a) {
    auto pos = 0;
    for (auto i = 1; i < 3; i++) {
        if (a[i] > a[pos]) pos = i;
    }
    return pos;
}
constexpr int min_element(const vec3i& a) {
    auto pos = 0;
    for (auto i = 1; i < 3; i++) {
        if (a[i] < a[pos]) pos = i;
    }
    return pos;
}

// Functions applied to vector elements
inline vec3i abs(const vec3i& a) { return {abs(a.x), abs(a.y), abs(a.z)}; };

// Vector comparison operations.
constexpr bool operator==(const vec4i& a, const vec4i& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}
constexpr bool operator!=(const vec4i& a, const vec4i& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w;
}

// Vector operations.
constexpr vec4i operator+(const vec4i& a) { return a; }
constexpr vec4i operator-(const vec4i& a) { return {-a.x, -a.y, -a.z, -a.w}; }
constexpr vec4i operator+(const vec4i& a, const vec4i& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}
constexpr vec4i operator+(const vec4i& a, int b) {
    return {a.x + b, a.y + b, a.z + b, a.w + b};
}
constexpr vec4i operator+(int a, const vec4i& b) {
    return {a + b.x, a + b.y, a + b.z, a + b.w};
}
constexpr vec4i operator-(const vec4i& a, const vec4i& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w};
}
constexpr vec4i operator-(const vec4i& a, int b) {
    return {a.x - b, a.y - b, a.z - b, a.w - b};
}
constexpr vec4i operator-(int a, const vec4i& b) {
    return {a - b.x, a - b.y, a - b.z, a - b.w};
}
constexpr vec4i operator*(const vec4i& a, const vec4i& b) {
    return {a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w};
}
constexpr vec4i operator*(const vec4i& a, int b) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
}
constexpr vec4i operator*(int a, const vec4i& b) {
    return {a * b.x, a * b.y, a * b.z, a * b.w};
}
constexpr vec4i operator/(const vec4i& a, const vec4i& b) {
    return {a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w};
}
constexpr vec4i operator/(const vec4i& a, int b) {
    return {a.x / b, a.y / b, a.z / b, a.w / b};
}
constexpr vec4i operator/(int a, const vec4i& b) {
    return {a / b.x, a / b.y, a / b.z, a / b.w};
}

// Vector assignments
constexpr vec4i& operator+=(vec4i& a, const vec4i& b) { return a = a + b; }
constexpr vec4i& operator+=(vec4i& a, int b) { return a = a + b; }
constexpr vec4i& operator-=(vec4i& a, const vec4i& b) { return a = a - b; }
constexpr vec4i& operator-=(vec4i& a, int b) { return a = a - b; }
constexpr vec4i& operator*=(vec4i& a, const vec4i& b) { return a = a * b; }
constexpr vec4i& operator*=(vec4i& a, int b) { return a = a * b; }
constexpr vec4i& operator/=(vec4i& a, const vec4i& b) { return a = a / b; }
constexpr vec4i& operator/=(vec4i& a, int b) { return a = a / b; }

// Max element and clamp.
constexpr vec4i max(const vec4i& a, int b) {
    return {max(a.x, b), max(a.y, b), max(a.z, b), max(a.w, b)};
}
constexpr vec4i min(const vec4i& a, int b) {
    return {min(a.x, b), min(a.y, b), min(a.z, b), min(a.w, b)};
}
constexpr vec4i max(const vec4i& a, const vec4i& b) {
    return {max(a.x, b.x), max(a.y, b.y), max(a.z, b.z), max(a.w, b.w)};
}
constexpr vec4i min(const vec4i& a, const vec4i& b) {
    return {min(a.x, b.x), min(a.y, b.y), min(a.z, b.z), min(a.w, b.w)};
}
constexpr vec4i clamp(const vec4i& x, int min, int max) {
    return {clamp(x.x, min, max), clamp(x.y, min, max), clamp(x.z, min, max),
        clamp(x.w, min, max)};
}

constexpr int max(const vec4i& a) { return max(max(max(a.x, a.y), a.z), a.w); }
constexpr int min(const vec4i& a) { return min(min(min(a.x, a.y), a.z), a.w); }
constexpr int sum(const vec4i& a) { return a.x + a.y + a.z + a.w; }
constexpr int max_element(const vec4i& a) {
    auto pos = 0;
    for (auto i = 1; i < 4; i++) {
        if (a[i] > a[pos]) pos = i;
    }
    return pos;
}
constexpr int min_element(const vec4i& a) {
    auto pos = 0;
    for (auto i = 1; i < 4; i++) {
        if (a[i] < a[pos]) pos = i;
    }
    return pos;
}

// Functions applied to vector elements
inline vec4i abs(const vec4i& a) {
    return {abs(a.x), abs(a.y), abs(a.z), abs(a.w)};
};

}  // namespace yocto

namespace std {

// Hash functor for vector for use with unordered_map
template <typename T, size_t N>
struct hash<std::array<T, N>> {
    static constexpr std::hash<T> hasher = std::hash<T>();
    size_t                        operator()(const std::array<T, N>& v) const {
        auto h = (size_t)0;
        for (auto i = 0; i < N; i++)
            h ^= hasher(v[i]) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};
template <>
struct hash<yocto::vec2i> {
    size_t operator()(const yocto::vec2i& v) const {
        return hash<std::array<int, 2>>()((const std::array<int, 2>&)v);
    }
};
template <>
struct hash<yocto::vec3i> {
    size_t operator()(const yocto::vec3i& v) const {
        return hash<std::array<int, 3>>()((const std::array<int, 3>&)v);
    }
};
template <>
struct hash<yocto::vec4i> {
    size_t operator()(const yocto::vec4i& v) const {
        return hash<std::array<int, 4>>()((const std::array<int, 4>&)v);
    }
};

}  // namespace std

// -----------------------------------------------------------------------------
// MATRICES
// -----------------------------------------------------------------------------
namespace yocto {

// Small Fixed-size matrices stored in column major format.
struct mat2f {
    vec2f x = {0, 0};
    vec2f y = {0, 0};

    constexpr mat2f() {}
    constexpr mat2f(const vec2f& x, const vec2f& y) : x{x}, y{y} {}

    constexpr vec2f&       operator[](int i) { return (&x)[i]; }
    constexpr const vec2f& operator[](int i) const { return (&x)[i]; }
};

// Small Fixed-size matrices stored in column major format.
struct mat3f {
    vec3f x = {0, 0, 0};
    vec3f y = {0, 0, 0};
    vec3f z = {0, 0, 0};

    constexpr mat3f() {}
    constexpr mat3f(const vec3f& x, const vec3f& y, const vec3f& z)
        : x{x}, y{y}, z{z} {}

    constexpr vec3f&       operator[](int i) { return (&x)[i]; }
    constexpr const vec3f& operator[](int i) const { return (&x)[i]; }
};

// Small Fixed-size matrices stored in column major format.
struct mat4f {
    vec4f x = {0, 0, 0, 0};
    vec4f y = {0, 0, 0, 0};
    vec4f z = {0, 0, 0, 0};
    vec4f w = {0, 0, 0, 0};

    constexpr mat4f() {}
    constexpr mat4f(
        const vec4f& x, const vec4f& y, const vec4f& z, const vec4f& w)
        : x{x}, y{y}, z{z}, w{w} {}

    constexpr vec4f&       operator[](int i) { return (&x)[i]; }
    constexpr const vec4f& operator[](int i) const { return (&x)[i]; }
};

// Identity matrices constants.
constexpr auto identity_mat2f = mat2f{{1, 0}, {0, 1}};
constexpr auto identity_mat3f = mat3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
constexpr auto identity_mat4f = mat4f{
    {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};

// Matrix comparisons.
constexpr bool operator==(const mat2f& a, const mat2f& b) {
    return a.x == b.x && a.y == b.y;
}
constexpr bool operator!=(const mat2f& a, const mat2f& b) { return !(a == b); }

// Matrix operations.
constexpr mat2f operator+(const mat2f& a, const mat2f& b) {
    return {a.x + b.x, a.y + b.y};
}
constexpr mat2f operator*(const mat2f& a, float b) {
    return {a.x * b, a.y * b};
}
constexpr vec2f operator*(const mat2f& a, const vec2f& b) {
    return a.x * b.x + a.y * b.y;
}
constexpr vec2f operator*(const vec2f& a, const mat2f& b) {
    return {dot(a, b.x), dot(a, b.y)};
}
constexpr mat2f operator*(const mat2f& a, const mat2f& b) {
    return {a * b.x, a * b.y};
}

// Matrix assignments.
constexpr mat2f& operator+=(mat2f& a, const mat2f& b) { return a = a + b; }
constexpr mat2f& operator*=(mat2f& a, const mat2f& b) { return a = a * b; }
constexpr mat2f& operator*=(mat2f& a, float b) { return a = a * b; }

// Matrix diagonals and transposes.
constexpr vec2f diagonal(const mat2f& a) { return {a.x.x, a.y.y}; }
constexpr mat2f transpose(const mat2f& a) {
    return {{a.x.x, a.y.x}, {a.x.y, a.y.y}};
}

// Matrix adjoints, determinants and inverses.
constexpr float determinant(const mat2f& a) { return cross(a.x, a.y); }
constexpr mat2f adjoint(const mat2f& a) {
    return {{a.y.y, -a.x.y}, {-a.y.x, a.x.x}};
}
constexpr mat2f inverse(const mat2f& a) {
    return adjoint(a) * (1 / determinant(a));
}

// Matrix comparisons.
constexpr bool operator==(const mat3f& a, const mat3f& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}
constexpr bool operator!=(const mat3f& a, const mat3f& b) { return !(a == b); }

// Matrix operations.
constexpr mat3f operator+(const mat3f& a, const mat3f& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
constexpr mat3f operator*(const mat3f& a, float b) {
    return {a.x * b, a.y * b, a.z * b};
}
constexpr vec3f operator*(const mat3f& a, const vec3f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
constexpr vec3f operator*(const vec3f& a, const mat3f& b) {
    return {dot(a, b.x), dot(a, b.y), dot(a, b.z)};
}
constexpr mat3f operator*(const mat3f& a, const mat3f& b) {
    return {a * b.x, a * b.y, a * b.z};
}

// Matrix assignments.
constexpr mat3f& operator+=(mat3f& a, const mat3f& b) { return a = a + b; }
constexpr mat3f& operator*=(mat3f& a, const mat3f& b) { return a = a * b; }
constexpr mat3f& operator*=(mat3f& a, float b) { return a = a * b; }

// Matrix diagonals and transposes.
constexpr vec3f diagonal(const mat3f& a) { return {a.x.x, a.y.y, a.z.z}; }
constexpr mat3f transpose(const mat3f& a) {
    return {
        {a.x.x, a.y.x, a.z.x},
        {a.x.y, a.y.y, a.z.y},
        {a.x.z, a.y.z, a.z.z},
    };
}

// Matrix adjoints, determinants and inverses.
constexpr float determinant(const mat3f& a) {
    return dot(a.x, cross(a.y, a.z));
}
constexpr mat3f adjoint(const mat3f& a) {
    return transpose(mat3f{cross(a.y, a.z), cross(a.z, a.x), cross(a.x, a.y)});
}
constexpr mat3f inverse(const mat3f& a) {
    return adjoint(a) * (1 / determinant(a));
}

// Constructs a basis from a direction
inline mat3f make_basis_fromz(const vec3f& v) {
    // https://graphics.pixar.com/library/OrthonormalB/paper.pdf
    auto z    = normalize(v);
    auto sign = copysignf(1.0f, z.z);
    auto a    = -1.0f / (sign + z.z);
    auto b    = z.x * z.y * a;
    auto x    = vec3f{1.0f + sign * z.x * z.x * a, sign * b, -sign * z.x};
    auto y    = vec3f{b, sign + z.y * z.y * a, -z.y};
    return {x, y, z};
}

// Matrix comparisons.
constexpr bool operator==(const mat4f& a, const mat4f& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}
constexpr bool operator!=(const mat4f& a, const mat4f& b) { return !(a == b); }

// Matrix operations.
constexpr mat4f operator+(const mat4f& a, const mat4f& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}
constexpr mat4f operator*(const mat4f& a, float b) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
}
constexpr vec4f operator*(const mat4f& a, const vec4f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}
constexpr vec4f operator*(const vec4f& a, const mat4f& b) {
    return {dot(a, b.x), dot(a, b.y), dot(a, b.z), dot(a, b.w)};
}
constexpr mat4f operator*(const mat4f& a, const mat4f& b) {
    return {a * b.x, a * b.y, a * b.z, a * b.w};
}

// Matrix assignments.
constexpr mat4f& operator+=(mat4f& a, const mat4f& b) { return a = a + b; }
constexpr mat4f& operator*=(mat4f& a, const mat4f& b) { return a = a * b; }
constexpr mat4f& operator*=(mat4f& a, float b) { return a = a * b; }

// Matrix diagonals and transposes.
constexpr vec4f diagonal(const mat4f& a) {
    return {a.x.x, a.y.y, a.z.z, a.w.w};
}
constexpr mat4f transpose(const mat4f& a) {
    return {
        {a.x.x, a.y.x, a.z.x, a.w.x},
        {a.x.y, a.y.y, a.z.y, a.w.y},
        {a.x.z, a.y.z, a.z.z, a.w.z},
        {a.x.w, a.y.w, a.z.w, a.w.w},
    };
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
    constexpr frame2f(const mat2f& m, const vec2f& t) : x{m.x}, y{m.y}, o{t} {}
    constexpr explicit frame2f(const mat3f& m)
        : x{m.x.x, m.x.y}, y{m.y.x, m.y.y}, o{m.z.x, m.z.y} {}
    constexpr operator mat3f() const { return {{x, 0}, {y, 0}, {o, 1}}; }

    constexpr vec2f&       operator[](int i) { return (&x)[i]; }
    constexpr const vec2f& operator[](int i) const { return (&x)[i]; }

    constexpr mat2f&       m() { return *(mat2f*)&x; }
    constexpr const mat2f& m() const { return *(mat2f*)&x; }
    constexpr vec2f&       t() { return o; }
    constexpr const vec2f& t() const { return o; }
};

// Rigid frames stored as a column-major affine transform matrix.
struct frame3f {
    vec3f x = {1, 0, 0};
    vec3f y = {0, 1, 0};
    vec3f z = {0, 0, 1};
    vec3f o = {0, 0, 0};

    constexpr frame3f() : x{}, y{}, z{}, o{} {}
    constexpr frame3f(
        const vec3f& x, const vec3f& y, const vec3f& z, const vec3f& o)
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
    constexpr vec3f&       t() { return o; }
    constexpr const vec3f& t() const { return o; }
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
    constexpr bbox2f(const vec2f& min, const vec2f& max) : min{min}, max{max} {}

    constexpr vec2f&       operator[](int i) { return (&min)[i]; }
    constexpr const vec2f& operator[](int i) const { return (&min)[i]; }
};

// Axis aligned bounding box represented as a min/max vector pairs.
struct bbox3f {
    vec3f min = {float_max, float_max, float_max};
    vec3f max = {float_min, float_min, float_min};

    constexpr bbox3f() {}
    constexpr bbox3f(const vec3f& min, const vec3f& max) : min{min}, max{max} {}

    constexpr vec3f&       operator[](int i) { return (&min)[i]; }
    constexpr const vec3f& operator[](int i) const { return (&min)[i]; }
};

// Empty bbox constant.
constexpr auto invalid_bbox2f = bbox2f{};
constexpr auto invalid_bbox3f = bbox3f{};

// Bounding box properties
constexpr vec2f bbox_center(const bbox2f& a) { return (a.min + a.max) / 2; }
constexpr vec2f bbox_size(const bbox2f& a) { return a.max - a.min; }

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
constexpr bbox2f& operator+=(bbox2f& a, const vec2f& b) { return a = a + b; }
constexpr bbox2f& operator+=(bbox2f& a, const bbox2f& b) { return a = a + b; }

// Bounding box properties
constexpr vec3f bbox_center(const bbox3f& a) { return (a.min + a.max) / 2; }
constexpr vec3f bbox_size(const bbox3f& a) { return a.max - a.min; }

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
constexpr bbox3f& operator+=(bbox3f& a, const vec3f& b) { return a = a + b; }
constexpr bbox3f& operator+=(bbox3f& a, const bbox3f& b) { return a = a + b; }

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
constexpr bbox3f quad_bounds(
    const vec3f& p0, const vec3f& p1, const vec3f& p2, const vec3f& p3) {
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
    vec2f o    = {0, 0};
    vec2f d    = {0, 1};
    float tmin = ray_eps;
    float tmax = float_max;

    constexpr ray2f() {}
    constexpr ray2f(const vec2f& o, const vec2f& d, float tmin = ray_eps,
        float tmax = float_max)
        : o{o}, d{d}, tmin{tmin}, tmax{tmax} {}
};

// Rays with origin, direction and min/max t value.
struct ray3f {
    vec3f o    = {0, 0, 0};
    vec3f d    = {0, 0, 1};
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

// -----------------------------------------------------------------------------
// CONVERSION TO STRING
// -----------------------------------------------------------------------------
namespace yocto {

inline string to_string(const char* value, bool quoted = false) {
    return quoted ? "\""s + value + "\""s : value;
}
inline string to_string(const string& value, bool quoted = false) {
    return quoted ? "\"" + value + "\"" : value;
}

inline string to_string(bool value, bool alpha = false) {
    if (alpha) {
        return value ? "true" : "false";
    } else {
        return value ? "1" : "0";
    }
}

inline string to_string(int value) { return std::to_string(value); }
inline string to_string(float value) { return std::to_string(value); }
inline string to_string(double value) { return std::to_string(value); }
inline string to_string(uint64_t value) { return std::to_string(value); }
inline string to_string(size_t value) { return std::to_string(value); }

template <typename T>
inline string to_string(const T* values, int num, bool bracketed = false) {
    auto str = ""s;
    for (auto i = 0; i < num; i++) {
        if (i) str += " ";
        str += to_string(values[i]);
    }
    return str;
}

inline string to_string(const vec2f& value, bool bracketed = false) {
    return to_string(&value.x, 2, bracketed);
}
inline string to_string(const vec3f& value, bool bracketed = false) {
    return to_string(&value.x, 3, bracketed);
}
inline string to_string(const vec4f& value, bool bracketed = false) {
    return to_string(&value.x, 4, bracketed);
}

inline string to_string(const vec2i& value, bool bracketed = false) {
    return to_string(&value.x, 2, bracketed);
}
inline string to_string(const vec3i& value, bool bracketed = false) {
    return to_string(&value.x, 3, bracketed);
}
inline string to_string(const vec4i& value, bool bracketed = false) {
    return to_string(&value.x, 4, bracketed);
}

inline string to_string(const mat2f& value, bool bracketed = false) {
    return to_string(&value.x.x, 4, bracketed);
}
inline string to_string(const mat3f& value, bool bracketed = false) {
    return to_string(&value.x.x, 9, bracketed);
}
inline string to_string(const mat4f& value, bool bracketed = false) {
    return to_string(&value.x.x, 16, bracketed);
}

inline string to_string(const frame2f& value, bool bracketed = false) {
    return to_string(&value.x.x, 6, bracketed);
}
inline string to_string(const frame3f& value, bool bracketed = false) {
    return to_string(&value.x.x, 12, bracketed);
}

inline string to_string(const bbox2f& value, bool bracketed = false) {
    return to_string(&value.min.x, 4, bracketed);
}
inline string to_string(const bbox3f& value, bool bracketed = false) {
    return to_string(&value.min.x, 6, bracketed);
}

}  // namespace yocto

#endif
