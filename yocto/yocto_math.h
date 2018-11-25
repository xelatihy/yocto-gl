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
// in graphics. In particular, we support 1-4 dimensional vectors of float
// coordinates (`vec1f`, `vec2f`, `vec3f`, `vec4f`) and int coordinates
// (`vec1i`, `vec2i`, `vec3i`, `vec4i`).
//
// We support 2-4 dimensional float matrices (`mat2f`, `mat3f`, `mat4f`) with
// matrix-matrix and matrix-vector products, transposes and inverses.
// Matrices are stored in column-major order and are accessed and constructed
// by column. The one dimensional version is for completeness only.
//
// To represent transformations, most of the library facilities prefer the use
// coordinate frames, aka rigid transforms, represented as `frame2f` and
// `frame3f`. The structure store three coordinate axes and the origin.
// This is equivalent to a rigid transform written as a column-major affine
// matrix. Transform operations are better behaved with this representation.
//
// We represent ranges of values in 1-4 dimensions with `bbox1f`, `bbox2f`,
// `bbox3f`, `bbox4f`, and `bbox1i`, `bbox2i`, `bbox3i`, `bbox4i`. Each range
// support construction from points and other ranges.
// These can be used to represent generic ranges and axis-aligned bounding
// boxes, for which we define the aliases `bbox1f`, `bbox2f`, `bbox3f`,`bbox4f`.
// We provide operations to compute bounds for points, lines, triangles and
// quads.
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

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <algorithm>
#include <array>
#include <atomic>
#include <cfloat>
#include <climits>
#include <cmath>
#include <cstdint>
#include <functional>
#include <string>
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
using std::pow;
using std::round;
using std::sin;
using std::sqrt;
using std::swap;
using std::tan;

using std::array;
using std::atomic;
using std::function;
using std::get;
using std::ignore;
using std::pair;
using std::string;
using std::tie;
using std::tuple;
using std::unordered_map;
using std::vector;
using namespace std::string_literals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// MATH CONSTANTS AND FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

using byte = unsigned char;
using uint = unsigned int;

const auto pi  = 3.14159265358979323846;
const auto pif = 3.14159265f;

const auto int_max       = INT_MAX;
const auto int_min       = -INT_MAX;
const auto float_max     = FLT_MAX;
const auto float_min     = -FLT_MAX;
const auto float_epsilon = FLT_EPSILON;

inline int   min(int x, int y) { return (x < y) ? x : y; }
inline float min(float x, float y) { return (x < y) ? x : y; }
inline int   max(int x, int y) { return (x > y) ? x : y; }
inline float max(float x, float y) { return (x > y) ? x : y; }
inline int clamp(int x, int min_, int max_) { return min(max(x, min_), max_); }
inline float clamp(float x, float min_, float max_) {
    return min(max(x, min_), max_);
}
inline float lerp(float a, float b, float u) { return a * (1 - u) + b * u; }
inline int   pow2(int x) { return 1 << x; }

}  // namespace yocto

// -----------------------------------------------------------------------------
// VECTORS
// -----------------------------------------------------------------------------
namespace yocto {

// Small size vectors.
struct vec1f {
    float x = 0;
};
struct vec2f {
    float x = 0;
    float y = 0;
};
struct vec3f {
    float x = 0;
    float y = 0;
    float z = 0;
};
struct vec4f {
    float x = 0;
    float y = 0;
    float z = 0;
    float w = 0;
};

// Zero vector constants.
const auto zero2f = vec2f{0, 0};
const auto zero3f = vec3f{0, 0, 0};
const auto zero4f = vec4f{0, 0, 0, 0};

// Access component by index.
inline float  at(const vec2f& v, int i) { return *(&v.x + i); }
inline float& at(vec2f& v, int i) { return *(&v.x + i); }
inline float  at(const vec3f& v, int i) { return *(&v.x + i); }
inline float& at(vec3f& v, int i) { return *(&v.x + i); }
inline float  at(const vec4f& v, int i) { return *(&v.x + i); }
inline float& at(vec4f& v, int i) { return *(&v.x + i); }

// Access xyz component of a vec4 typically used for color operation.
inline vec3f& xyz(const vec4f& a) { return (vec3f&)a; }
inline vec3f& xyz(vec4f& a) { return (vec3f&)a; }

// Vector comparison operations.
inline bool operator==(const vec2f& a, const vec2f& b) {
    return a.x == b.x && a.y == b.y;
}
inline bool operator!=(const vec2f& a, const vec2f& b) {
    return a.x != b.x || a.y != b.y;
}
inline bool operator==(const vec3f& a, const vec3f& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}
inline bool operator!=(const vec3f& a, const vec3f& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z;
}
inline bool operator==(const vec4f& a, const vec4f& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}
inline bool operator!=(const vec4f& a, const vec4f& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w;
}

// Vector operations.
inline vec2f operator-(const vec2f& a) { return {-a.x, -a.y}; }
inline vec2f operator+(const vec2f& a, const vec2f& b) {
    return {a.x + b.x, a.y + b.y};
}
inline vec2f operator+(const vec2f& a, float b) { return {a.x + b, a.y + b}; }
inline vec2f operator+(float a, const vec2f& b) { return {a + b.x, a + b.y}; }
inline vec2f operator-(const vec2f& a, const vec2f& b) {
    return {a.x - b.x, a.y - b.y};
}
inline vec2f operator-(const vec2f& a, float b) { return {a.x - b, a.y - b}; }
inline vec2f operator-(float a, const vec2f& b) { return {a - b.x, a - b.y}; }
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
inline vec3f operator+(const vec3f& a, float b) {
    return {a.x + b, a.y + b, a.z + b};
}
inline vec3f operator+(float a, const vec3f& b) {
    return {a + b.x, a + b.y, a + b.z};
}
inline vec3f operator-(const vec3f& a, const vec3f& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
inline vec3f operator-(const vec3f& a, float b) {
    return {a.x - b, a.y - b, a.z - b};
}
inline vec3f operator-(float a, const vec3f& b) {
    return {a - b.x, a - b.y, a - b.z};
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
inline vec4f operator+(const vec4f& a, float b) {
    return {a.x + b, a.y + b, a.z + b, a.w + b};
}
inline vec4f operator+(float a, const vec4f& b) {
    return {a + b.x, a + b.y, a + b.z, a + b.w};
}
inline vec4f operator-(const vec4f& a, const vec4f& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w};
}
inline vec4f operator-(const vec4f& a, float b) {
    return {a.x - b, a.y - b, a.z - b, a.w - b};
}
inline vec4f operator-(float a, const vec4f& b) {
    return {a - b.x, a - b.y, a - b.z, a - b.w};
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
inline vec2f& operator+=(vec2f& a, float b) { return a = a + b; }
inline vec2f& operator-=(vec2f& a, const vec2f& b) { return a = a - b; }
inline vec2f& operator-=(vec2f& a, float b) { return a = a - b; }
inline vec2f& operator*=(vec2f& a, const vec2f& b) { return a = a * b; }
inline vec2f& operator*=(vec2f& a, float b) { return a = a * b; }
inline vec2f& operator/=(vec2f& a, const vec2f& b) { return a = a / b; }
inline vec2f& operator/=(vec2f& a, float b) { return a = a / b; }

// Vector assignments
inline vec3f& operator+=(vec3f& a, const vec3f& b) { return a = a + b; }
inline vec3f& operator+=(vec3f& a, float b) { return a = a + b; }
inline vec3f& operator-=(vec3f& a, const vec3f& b) { return a = a - b; }
inline vec3f& operator-=(vec3f& a, float b) { return a = a - b; }
inline vec3f& operator*=(vec3f& a, const vec3f& b) { return a = a * b; }
inline vec3f& operator*=(vec3f& a, float b) { return a = a * b; }
inline vec3f& operator/=(vec3f& a, const vec3f& b) { return a = a / b; }
inline vec3f& operator/=(vec3f& a, float b) { return a = a / b; }

// Vector assignments
inline vec4f& operator+=(vec4f& a, const vec4f& b) { return a = a + b; }
inline vec4f& operator+=(vec4f& a, float b) { return a = a + b; }
inline vec4f& operator-=(vec4f& a, const vec4f& b) { return a = a - b; }
inline vec4f& operator-=(vec4f& a, float b) { return a = a - b; }
inline vec4f& operator*=(vec4f& a, const vec4f& b) { return a = a * b; }
inline vec4f& operator*=(vec4f& a, float b) { return a = a * b; }
inline vec4f& operator/=(vec4f& a, const vec4f& b) { return a = a / b; }
inline vec4f& operator/=(vec4f& a, float b) { return a = a / b; }

// Vector products and lengths.
// Vector products and lengths.
inline float dot(const vec2f& a, const vec2f& b) {
    return a.x * b.x + a.y * b.y;
}
inline float length(const vec2f& a) { return sqrt(a.x * a.x + a.y * a.y); }
inline vec2f normalize(const vec2f& a) {
    auto l = length(a);
    return (l) ? a / l : a;
}
inline float cross(const vec2f& a, const vec2f& b) {
    return a.x * b.y - a.y * b.x;
}
inline float dot(const vec3f& a, const vec3f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline float length(const vec3f& a) {
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z);
}
inline vec3f normalize(const vec3f& a) {
    auto l = length(a);
    return (l) ? a / l : a;
}
inline vec3f cross(const vec3f& a, const vec3f& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline float dot(const vec4f& a, const vec4f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}
inline float length(const vec4f& a) {
    return sqrt(a.x * a.x + a.y * a.y + a.z * a.z + a.w * a.w);
}
inline vec4f normalize(const vec4f& a) {
    auto l = length(a);
    return (l) ? a / l : a;
}
inline float distance(const vec2f& a, const vec2f& b) { return length(a - b); }
inline float distance(const vec3f& a, const vec3f& b) { return length(a - b); }
inline float distance(const vec4f& a, const vec4f& b) { return length(a - b); }
inline float distance_squared(const vec2f& a, const vec2f& b) {
    return dot(a - b, a - b);
}
inline float distance_squared(const vec3f& a, const vec3f& b) {
    return dot(a - b, a - b);
}
inline float distance_squared(const vec4f& a, const vec4f& b) {
    return dot(a - b, a - b);
}

// Vector angles and slerps.
inline float angle(const vec3f& a, const vec3f& b) {
    return acos(clamp(dot(normalize(a), normalize(b)), -1.0f, 1.0f));
}
inline vec4f slerp(const vec4f& a, const vec4f& b, float u) {
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
    if (k < 0) return {0, 0, 0};  // tir
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
inline float mean(const vec2f& a) { return (a.x + a.y) / 2; }
inline float mean(const vec3f& a) { return (a.x + a.y + a.z) / 3; }
inline float mean(const vec4f& a) { return (a.x + a.y + a.z + a.w) / 4; }

// Apply a unary function to all vector elements
inline vec2f apply(float (*func)(float), const vec2f& a) {
    return {func(a.x), func(a.y)};
}
inline vec3f apply(float (*func)(float), const vec3f& a) {
    return {func(a.x), func(a.y), func(a.z)};
}
inline vec4f apply(float (*func)(float), const vec4f& a) {
    return {func(a.x), func(a.y), func(a.z), func(a.w)};
}
// Apply a binary function to all vector elements
inline vec2f apply(float (*func)(float, float), const vec2f& a, float b) {
    return {func(a.x, b), func(a.y, b)};
}
inline vec3f apply(float (*func)(float, float), const vec3f& a, float b) {
    return {func(a.x, b), func(a.y, b), func(a.z, b)};
}
inline vec4f apply(float (*func)(float, float), const vec4f& a, float b) {
    return {func(a.x, b), func(a.y, b), func(a.z, b), func(a.w, b)};
}

// Functions applied to vector elements
inline vec2f exp(const vec2f& a) { return apply(exp, a); };
inline vec3f exp(const vec3f& a) { return apply(exp, a); };
inline vec4f exp(const vec4f& a) { return apply(exp, a); };
inline vec2f pow(const vec2f& a, float b) { return apply(pow, a, b); };
inline vec3f pow(const vec3f& a, float b) { return apply(pow, a, b); };
inline vec4f pow(const vec4f& a, float b) { return apply(pow, a, b); };

inline bool isfinite(const vec2f& a) { return isfinite(a.x) && isfinite(a.x); };
inline bool isfinite(const vec3f& a) {
    return isfinite(a.x) && isfinite(a.x) && isfinite(a.z);
};
inline bool isfinite(const vec4f& a) {
    return isfinite(a.x) && isfinite(a.x) && isfinite(a.z) && isfinite(a.w);
};

// Quaternion operatons represented as xi + yj + zk + w
const auto   identity_quat4f = vec4f{0, 0, 0, 1};
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
    return quat_conjugate(a) / dot(a, a);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// INT VECTORS
// -----------------------------------------------------------------------------
namespace yocto {

// Small size vectors.
struct vec2i {
    int x = 0;
    int y = 0;
};
struct vec3i {
    int x = 0;
    int y = 0;
    int z = 0;
};
struct vec4i {
    int x = 0;
    int y = 0;
    int z = 0;
    int w = 0;
};
struct vec4b {
    byte x = 0;
    byte y = 0;
    byte z = 0;
    byte w = 0;
};

// Zero vector constants.
const auto zero2i = vec2i{0, 0};
const auto zero3i = vec3i{0, 0, 0};
const auto zero4i = vec4i{0, 0, 0, 0};
const auto zero4b = vec4b{0, 0, 0, 0};

// Access component by index.
inline int  at(const vec2i& v, int i) { return *(&v.x + i); }
inline int& at(vec2i& v, int i) { return *(&v.x + i); }
inline int  at(const vec3i& v, int i) { return *(&v.x + i); }
inline int& at(vec3i& v, int i) { return *(&v.x + i); }
inline int  at(const vec4i& v, int i) { return *(&v.x + i); }
inline int& at(vec4i& v, int i) { return *(&v.x + i); }

// Access xyz component of a vec4 typically used for color operation.
inline vec3i& xyz(const vec4i& a) { return (vec3i&)a; }
inline vec3i& xyz(vec4i& a) { return (vec3i&)a; }

// Vector comparison operations.
inline bool operator==(const vec2i& a, const vec2i& b) {
    return a.x == b.x && a.y == b.y;
}
inline bool operator!=(const vec2i& a, const vec2i& b) {
    return a.x != b.x || a.y != b.y;
}
inline bool operator==(const vec3i& a, const vec3i& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
}
inline bool operator!=(const vec3i& a, const vec3i& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z;
}
inline bool operator==(const vec4i& a, const vec4i& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}
inline bool operator!=(const vec4i& a, const vec4i& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w;
}
inline bool operator==(const vec4b& a, const vec4b& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}
inline bool operator!=(const vec4b& a, const vec4b& b) {
    return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w;
}

// Vector operations.
inline vec2i operator-(const vec2i& a) { return {-a.x, -a.y}; }
inline vec2i operator+(const vec2i& a, const vec2i& b) {
    return {a.x + b.x, a.y + b.y};
}
inline vec2i operator+(const vec2i& a, int b) { return {a.x + b, a.y + b}; }
inline vec2i operator+(int a, const vec2i& b) { return {a + b.x, a + b.y}; }
inline vec2i operator-(const vec2i& a, const vec2i& b) {
    return {a.x - b.x, a.y - b.y};
}
inline vec2i operator-(const vec2i& a, int b) { return {a.x - b, a.y - b}; }
inline vec2i operator-(int a, const vec2i& b) { return {a - b.x, a - b.y}; }
inline vec2i operator*(const vec2i& a, const vec2i& b) {
    return {a.x * b.x, a.y * b.y};
}
inline vec2i operator*(const vec2i& a, int b) { return {a.x * b, a.y * b}; }
inline vec2i operator*(int a, const vec2i& b) { return {a * b.x, a * b.y}; }
inline vec2i operator/(const vec2i& a, const vec2i& b) {
    return {a.x / b.x, a.y / b.y};
}
inline vec2i operator/(const vec2i& a, int b) { return {a.x / b, a.y / b}; }
inline vec2i operator/(int a, const vec2i& b) { return {a / b.x, a / b.y}; }

// Vector operations.
inline vec3i operator+(const vec3i& a) { return a; }
inline vec3i operator-(const vec3i& a) { return {-a.x, -a.y, -a.z}; }
inline vec3i operator+(const vec3i& a, const vec3i& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
}
inline vec3i operator+(const vec3i& a, int b) {
    return {a.x + b, a.y + b, a.z + b};
}
inline vec3i operator+(int a, const vec3i& b) {
    return {a + b.x, a + b.y, a + b.z};
}
inline vec3i operator-(const vec3i& a, const vec3i& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
}
inline vec3i operator-(const vec3i& a, int b) {
    return {a.x - b, a.y - b, a.z - b};
}
inline vec3i operator-(int a, const vec3i& b) {
    return {a - b.x, a - b.y, a - b.z};
}
inline vec3i operator*(const vec3i& a, const vec3i& b) {
    return {a.x * b.x, a.y * b.y, a.z * b.z};
}
inline vec3i operator*(const vec3i& a, int b) {
    return {a.x * b, a.y * b, a.z * b};
}
inline vec3i operator*(int a, const vec3i& b) {
    return {a * b.x, a * b.y, a * b.z};
}
inline vec3i operator/(const vec3i& a, const vec3i& b) {
    return {a.x / b.x, a.y / b.y, a.z / b.z};
}
inline vec3i operator/(const vec3i& a, int b) {
    return {a.x / b, a.y / b, a.z / b};
}
inline vec3i operator/(int a, const vec3i& b) {
    return {a / b.x, a / b.y, a / b.z};
}

// Vector operations.
inline vec4i operator-(const vec4i& a) { return {-a.x, -a.y, -a.z, -a.w}; }
inline vec4i operator+(const vec4i& a, const vec4i& b) {
    return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}
inline vec4i operator+(const vec4i& a, int b) {
    return {a.x + b, a.y + b, a.z + b, a.w + b};
}
inline vec4i operator+(int a, const vec4i& b) {
    return {a + b.x, a + b.y, a + b.z, a + b.w};
}
inline vec4i operator-(const vec4i& a, const vec4i& b) {
    return {a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w};
}
inline vec4i operator-(const vec4i& a, int b) {
    return {a.x - b, a.y - b, a.z - b, a.w - b};
}
inline vec4i operator-(int a, const vec4i& b) {
    return {a - b.x, a - b.y, a - b.z, a - b.w};
}
inline vec4i operator*(const vec4i& a, const vec4i& b) {
    return {a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w};
}
inline vec4i operator*(const vec4i& a, int b) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
}
inline vec4i operator*(int a, const vec4i& b) {
    return {a * b.x, a * b.y, a * b.z, a * b.w};
}
inline vec4i operator/(const vec4i& a, const vec4i& b) {
    return {a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w};
}
inline vec4i operator/(const vec4i& a, int b) {
    return {a.x / b, a.y / b, a.z / b, a.w / b};
}
inline vec4i operator/(int a, const vec4i& b) {
    return {a / b.x, a / b.y, a / b.z, a / b.w};
}

// Vector assignments
inline vec2i& operator+=(vec2i& a, const vec2i& b) { return a = a + b; }
inline vec2i& operator+=(vec2i& a, int b) { return a = a + b; }
inline vec2i& operator-=(vec2i& a, const vec2i& b) { return a = a - b; }
inline vec2i& operator-=(vec2i& a, int b) { return a = a - b; }
inline vec2i& operator*=(vec2i& a, const vec2i& b) { return a = a * b; }
inline vec2i& operator*=(vec2i& a, int b) { return a = a * b; }
inline vec2i& operator/=(vec2i& a, const vec2i& b) { return a = a / b; }
inline vec2i& operator/=(vec2i& a, int b) { return a = a / b; }

// Vector assignments
inline vec3i& operator+=(vec3i& a, const vec3i& b) { return a = a + b; }
inline vec3i& operator+=(vec3i& a, int b) { return a = a + b; }
inline vec3i& operator-=(vec3i& a, const vec3i& b) { return a = a - b; }
inline vec3i& operator-=(vec3i& a, int b) { return a = a - b; }
inline vec3i& operator*=(vec3i& a, const vec3i& b) { return a = a * b; }
inline vec3i& operator*=(vec3i& a, int b) { return a = a * b; }
inline vec3i& operator/=(vec3i& a, const vec3i& b) { return a = a / b; }
inline vec3i& operator/=(vec3i& a, int b) { return a = a / b; }

// Vector assignments
inline vec4i& operator+=(vec4i& a, const vec4i& b) { return a = a + b; }
inline vec4i& operator+=(vec4i& a, int b) { return a = a + b; }
inline vec4i& operator-=(vec4i& a, const vec4i& b) { return a = a - b; }
inline vec4i& operator-=(vec4i& a, int b) { return a = a - b; }
inline vec4i& operator*=(vec4i& a, const vec4i& b) { return a = a * b; }
inline vec4i& operator*=(vec4i& a, int b) { return a = a * b; }
inline vec4i& operator/=(vec4i& a, const vec4i& b) { return a = a / b; }
inline vec4i& operator/=(vec4i& a, int b) { return a = a / b; }

}  // namespace yocto

namespace std {

// Hash functor for vector for use with unordered_map
template <>
struct hash<yocto::vec2i> {
    size_t operator()(const yocto::vec2i& v) const {
        auto vh = hash<int>();
        auto h  = (size_t)0;
        for (auto i = 0; i < 2; i++)
            h ^= vh((&v.x)[i]) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};
// Hash functor for vector for use with unordered_map
template <>
struct hash<yocto::vec3i> {
    size_t operator()(const yocto::vec3i& v) const {
        auto vh = hash<int>();
        auto h  = (size_t)0;
        for (auto i = 0; i < 3; i++)
            h ^= vh((&v.x)[i]) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};
// Hash functor for vector for use with unordered_map
template <>
struct hash<yocto::vec4i> {
    size_t operator()(const yocto::vec4i& v) const {
        auto vh = hash<int>();
        auto h  = (size_t)0;
        for (auto i = 0; i < 4; i++)
            h ^= vh((&v.x)[i]) + 0x9e3779b9 + (h << 6) + (h >> 2);
        return h;
    }
};

}  // namespace std

// -----------------------------------------------------------------------------
// MATRICES
// -----------------------------------------------------------------------------
namespace yocto {

// Small Fixed-size square matrices stored in column major format.
struct mat2f {
    vec2f x = {1, 0};
    vec2f y = {0, 1};
};
struct mat3f {
    vec3f x = {1, 0, 0};
    vec3f y = {0, 1, 0};
    vec3f z = {0, 0, 1};
};
struct mat4f {
    vec4f x = {1, 0, 0, 0};
    vec4f y = {0, 1, 0, 0};
    vec4f z = {0, 0, 1, 0};
    vec4f w = {0, 0, 0, 1};
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

// Matrix assignments.
inline mat3f& operator+=(mat3f& a, const mat3f& b) { return a = a + b; }
inline mat3f& operator*=(mat3f& a, const mat3f& b) { return a = a * b; }
inline mat3f& operator*=(mat3f& a, float b) { return a = a * b; }

// Matrix assignments.
inline mat4f& operator+=(mat4f& a, const mat4f& b) { return a = a + b; }
inline mat4f& operator*=(mat4f& a, const mat4f& b) { return a = a * b; }
inline mat4f& operator*=(mat4f& a, float b) { return a = a * b; }

// Matrix diagonals and transposes.
inline vec2f diagonal(const mat2f& a) { return {a.x.x, a.y.y}; }
inline vec3f diagonal(const mat3f& a) { return {a.x.x, a.y.y, a.z.z}; }
inline vec4f diagonal(const mat4f& a) { return {a.x.x, a.y.y, a.z.z, a.w.w}; }
inline mat2f transpose(const mat2f& a);
inline mat3f transpose(const mat3f& a);
inline mat4f transpose(const mat4f& a);

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

// Constructs a basis from a direction
inline mat3f make_basis_fromz(const vec3f& v) {
    // https://graphics.pixar.com/library/OrthonormalB/paper.pdf
    auto z = normalize(v);
    auto sign = copysignf(1.0f, z.z);
    auto a = -1.0f / (sign + z.z);
    auto b = z.x * z.y * a;
    auto x = vec3f{1.0f + sign * z.x * z.x * a, sign * b, -sign * z.x};
    auto y = vec3f{b, sign + z.y * z.y * a, -z.y};
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
};
struct frame3f {
    vec3f x = {1, 0, 0};
    vec3f y = {0, 1, 0};
    vec3f z = {0, 0, 1};
    vec3f o = {0, 0, 0};
};

// Indentity frames.
const auto identity_frame2f = frame2f{{1, 0}, {0, 1}, {0, 0}};
const auto identity_frame3f = frame3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};

// Frame construction from axis.
inline frame3f make_frame_fromz(const vec3f& o, const vec3f& v) {
    // https://graphics.pixar.com/library/OrthonormalB/paper.pdf
    auto z = normalize(v);
    auto sign = copysignf(1.0f, z.z);
    auto a = -1.0f / (sign + z.z);
    auto b = z.x * z.y * a;
    auto x = vec3f{1.0f + sign * z.x * z.x * a, sign * b, -sign * z.x};
    auto y = vec3f{b, sign + z.y * z.y * a, -z.y};
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
    return {
        {a.x.x, a.x.y, a.x.z, 0},
        {a.y.x, a.y.y, a.y.z, 0},
        {a.z.x, a.z.y, a.z.z, 0},
        {a.o.x, a.o.y, a.o.z, 1},
    };
}
inline frame3f mat_to_frame(const mat4f& a) {
    return {
        {a.x.x, a.x.y, a.x.z},
        {a.y.x, a.y.y, a.y.z},
        {a.z.x, a.z.y, a.z.z},
        {a.w.x, a.w.y, a.w.z},
    };
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
inline frame2f inverse(const frame2f& a, bool is_rigid = true) {
    auto minv = (is_rigid) ? transpose(mat2f{a.x, a.y}) :
                             inverse(mat2f{a.x, a.y});
    return {minv.x, minv.y, -(minv * a.o)};
}
inline frame3f inverse(const frame3f& a, bool is_rigid = true) {
    auto minv = (is_rigid) ? transpose(mat3f{a.x, a.y, a.z}) :
                             inverse(mat3f{a.x, a.y, a.z});
    return {minv.x, minv.y, minv.z, -(minv * a.o)};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// AXIS ALIGNED BOUNDING BOXES
// -----------------------------------------------------------------------------
namespace yocto {

// Range of values in 1D.
struct bbox1f {
    float min = float_max;
    float max = float_min;
};

// Axis aligned bounding box represented as a min/max vector pairs.
struct bbox2f {
    vec2f min = {float_max, float_max};
    vec2f max = {float_min, float_min};
};

// Axis aligned bounding box represented as a min/max vector pairs.
struct bbox3f {
    vec3f min = {float_max, float_max, float_max};
    vec3f max = {float_min, float_min, float_min};
};

// Axis aligned bounding box represented as a min/max vector pairs.
struct bbox4f {
    vec4f min = {float_max, float_max, float_max, float_max};
    vec4f max = {float_min, float_min, float_min, float_min};
};

// Empty bbox constant.
const auto invalid_bbox1f = bbox1f();
const auto invalid_bbox2f = bbox2f();
const auto invalid_bbox3f = bbox3f();
const auto invalid_bbox4f = bbox4f();

// Bounding box values
inline float bbox_size(const bbox1f& a) { return a.max - a.min; }
inline vec2f bbox_size(const bbox2f& a) { return a.max - a.min; }
inline vec3f bbox_size(const bbox3f& a) { return a.max - a.min; }
inline vec4f bbox_size(const bbox4f& a) { return a.max - a.min; }
inline float bbox_center(const bbox1f& a) { return (a.max + a.min) / 2; }
inline vec2f bbox_center(const bbox2f& a) { return (a.max + a.min) / 2; }
inline vec3f bbox_center(const bbox3f& a) { return (a.max + a.min) / 2; }
inline vec4f bbox_center(const bbox4f& a) { return (a.max + a.min) / 2; }

// Bounding box comparisons.
inline bool operator==(const bbox1f& a, const bbox1f& b) {
    return a.min == b.min && a.max == b.max;
}
inline bool operator!=(const bbox1f& a, const bbox1f& b) {
    return a.min != b.min || a.max != b.max;
}
inline bool operator==(const bbox2f& a, const bbox2f& b) {
    return a.min == b.min && a.max == b.max;
}
inline bool operator!=(const bbox2f& a, const bbox2f& b) {
    return a.min != b.min || a.max != b.max;
}
inline bool operator==(const bbox3f& a, const bbox3f& b) {
    return a.min == b.min && a.max == b.max;
}
inline bool operator!=(const bbox3f& a, const bbox3f& b) {
    return a.min != b.min || a.max != b.max;
}
inline bool operator==(const bbox4f& a, const bbox4f& b) {
    return a.min == b.min && a.max == b.max;
}
inline bool operator!=(const bbox4f& a, const bbox4f& b) {
    return a.min != b.min || a.max != b.max;
}

// Bounding box expansions with points and other boxes.
inline bbox1f& operator+=(bbox1f& a, float b) {
    a.min = min(a.min, b);
    a.max = max(a.max, b);
    return a;
}
inline bbox1f& operator+=(bbox1f& a, const bbox1f& b) {
    a.min = min(a.min, b.min);
    a.max = max(a.max, b.max);
    return a;
}
// Bounding box expansions with points and other boxes.
inline bbox2f& operator+=(bbox2f& a, const vec2f& b) {
    a.min = {min(a.min.x, b.x), min(a.min.y, b.y)};
    a.max = {max(a.max.x, b.x), max(a.max.y, b.y)};
    return a;
}
inline bbox2f& operator+=(bbox2f& a, const bbox2f& b) {
    a.min = {min(a.min.x, b.min.x), min(a.min.y, b.min.y)};
    a.max = {max(a.max.x, b.max.x), max(a.max.y, b.max.y)};
    return a;
}
// Bounding box expansions with points and other boxes.
inline bbox3f& operator+=(bbox3f& a, const vec3f& b) {
    a.min = {min(a.min.x, b.x), min(a.min.y, b.y), min(a.min.z, b.z)};
    a.max = {max(a.max.x, b.x), max(a.max.y, b.y), max(a.max.z, b.z)};
    return a;
}
inline bbox3f& operator+=(bbox3f& a, const bbox3f& b) {
    a.min = {min(a.min.x, b.min.x), min(a.min.y, b.min.y), min(a.min.z, b.min.z)};
    a.max = {max(a.max.x, b.max.x), max(a.max.y, b.max.y), max(a.max.z, b.max.z)};
    return a;
}
// Bounding box expansions with points and other boxes.
inline bbox4f& operator+=(bbox4f& a, const vec4f& b) {
    a.min = {min(a.min.x, b.x), min(a.min.y, b.y), min(a.min.z, b.z),
        min(a.min.w, b.w)};
    a.max = {max(a.max.x, b.x), max(a.max.y, b.y), max(a.max.z, b.z),
        max(a.max.w, b.w)};
    return a;
}
inline bbox4f& operator+=(bbox4f& a, const bbox4f& b) {
    a.min = {min(a.min.x, b.min.x), min(a.min.y, b.min.y),
        min(a.min.z, b.min.z), min(a.min.w, b.min.w)};
    a.max = {max(a.max.x, b.max.x), max(a.max.y, b.max.y),
        max(a.max.z, b.max.z), max(a.max.w, b.max.w)};
    return a;
}

// Create bounding boxes from arrays
inline bbox1f make_bbox(const vector<float>& values) {
    auto bbox = bbox1f{};
    for (auto& value : values) bbox += value;
    return bbox;
}
inline bbox2f make_bbox(const vector<vec2f>& values) {
    auto bbox = bbox2f{};
    for (auto& value : values) bbox += value;
    return bbox;
}
inline bbox3f make_bbox(const vector<vec3f>& values) {
    auto bbox = bbox3f{};
    for (auto& value : values) bbox += value;
    return bbox;
}
inline bbox4f make_bbox(const vector<vec4f>& values) {
    auto bbox = bbox4f{};
    for (auto& value : values) bbox += value;
    return bbox;
}

// Primitive bounds.
inline bbox3f point_bounds(const vec3f& p, float r = 0) {
    auto bbox = bbox3f{};
    bbox += p - vec3f{r, r, r};
    bbox += p + vec3f{r, r, r};
    return bbox;
}
inline bbox3f line_bounds(
    const vec3f& p0, const vec3f& p1, float r0 = 0, float r1 = 0) {
    auto bbox = bbox3f{};
    bbox += p0 - vec3f{r0, r0, r0};
    bbox += p0 + vec3f{r0, r0, r0};
    bbox += p1 - vec3f{r1, r1, r1};
    bbox += p1 + vec3f{r1, r1, r1};
    return bbox;
}
inline bbox3f triangle_bounds(const vec3f& p0, const vec3f& p1, const vec3f& p2) {
    auto bbox = bbox3f{};
    bbox += p0;
    bbox += p1;
    bbox += p2;
    return bbox;
}
inline bbox3f quad_bounds(
    const vec3f& p0, const vec3f& p1, const vec3f& p2, const vec3f& p3) {
    auto bbox = bbox3f{};
    bbox += p0;
    bbox += p1;
    bbox += p2;
    bbox += p3;
    return bbox;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// AXIS ALIGNED BOUNDING BOXES
// -----------------------------------------------------------------------------
namespace yocto {

// Range of values in 1D.
struct bbox1i {
    float min = int_max;
    float max = int_min;
};

// Axis aligned bounding box represented as a min/max vector pairs.
struct bbox2i {
    vec2i min = {int_max, int_max};
    vec2i max = {int_min, int_min};
};

// Axis aligned bounding box represented as a min/max vector pairs.
struct bbox3i {
    vec3i min = {int_max, int_max, int_max};
    vec3i max = {int_min, int_min, int_min};
};

// Axis aligned bounding box represented as a min/max vector pairs.
struct bbox4i {
    vec4i min = {int_max, int_max, int_max, int_max};
    vec4i max = {int_min, int_min, int_min, int_min};
};

// Empty bbox constant.
const auto invalid_bbox1i = bbox1i();
const auto invalid_bbox2i = bbox2i();
const auto invalid_bbox3i = bbox3i();
const auto invalid_bbox4i = bbox4i();

// Bbox properties
inline int   bbox_size(const bbox1i& a) { return a.max - a.min; }
inline vec2i bbox_size(const bbox2i& a) { return a.max - a.min; }
inline vec3i bbox_size(const bbox3i& a) { return a.max - a.min; }
inline vec4i bbox_size(const bbox4i& a) { return a.max - a.min; }

// Bounding box comparisons.
inline bool operator==(const bbox1i& a, const bbox1i& b) {
    return a.min == b.min && a.max == b.max;
}
inline bool operator!=(const bbox1i& a, const bbox1i& b) {
    return a.min != b.min || a.max != b.max;
}
inline bool operator==(const bbox2i& a, const bbox2i& b) {
    return a.min == b.min && a.max == b.max;
}
inline bool operator!=(const bbox2i& a, const bbox2i& b) {
    return a.min != b.min || a.max != b.max;
}
inline bool operator==(const bbox3i& a, const bbox3i& b) {
    return a.min == b.min && a.max == b.max;
}
inline bool operator!=(const bbox3i& a, const bbox3i& b) {
    return a.min != b.min || a.max != b.max;
}
inline bool operator==(const bbox4i& a, const bbox4i& b) {
    return a.min == b.min && a.max == b.max;
}
inline bool operator!=(const bbox4i& a, const bbox4i& b) {
    return a.min != b.min || a.max != b.max;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// RAYS
// -----------------------------------------------------------------------------
namespace yocto {

// Rays with origin, direction and min/max t value.
struct ray2f {
    vec2f o    = {0, 0};
    vec2f d    = {0, 1};
    float tmin = 0;
    float tmax = float_max;
};

// Rays with origin, direction and min/max t value.
struct ray3f {
    vec3f o    = {0, 0, 0};
    vec3f d    = {0, 0, 1};
    float tmin = 0;
    float tmax = float_max;
};

// Ray esplison
const auto ray_eps = 1e-4f;

// Construct a ray from direction or segments using a default epsilon.
inline ray3f make_ray(const vec3f& o, const vec3f& d, float eps = ray_eps) {
    return {o, d, eps, float_max};
}
inline ray3f make_segment(const vec3f& p1, const vec3f& p2, float eps = ray_eps) {
    return {p1, normalize(p2 - p1), eps, length(p2 - p1) - 2 * eps};
}

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
inline vec3f transform_point(const frame3f& a, const vec3f& b) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.o;
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
inline vec2f transform_vector_inverse(const frame2f& a, const vec2f& b) {
    return {dot(b, a.x), dot(b, a.y)};
}
inline vec2f transform_direction_inverse(const frame2f& a, const vec2f& b) {
    return normalize(transform_vector_inverse(a, b));
}
inline vec3f transform_point_inverse(const frame3f& a, const vec3f& b) {
    return {dot(b - a.o, a.x), dot(b - a.o, a.y), dot(b - a.o, a.z)};
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
inline mat4f make_ortho_mat(float l, float r, float b, float t, float n, float f) {
    return {{2 / (r - l), 0, 0, 0}, {0, 2 / (t - b), 0, 0},
        {0, 0, -2 / (f - n), 0},
        {-(r + l) / (r - l), -(t + b) / (t - b), -(f + n) / (f - n), 1}};
}
inline mat4f make_ortho2d_mat(float left, float right, float bottom, float top) {
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
inline mat4f perspective_mat(float fovy, float aspect, float near) {
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

// Turntable and FPS Camera navigation.
inline void update_camera_turntable(vec3f& from, vec3f& to, vec3f& up,
    const vec2f& rotate, float dolly, const vec2f& pan);
inline void update_camera_turntable(frame3f& frame, float& focus,
    const vec2f& rotate, float dolly, const vec2f& pan);
inline void update_camera_fps(
    frame3f& frame, const vec3f& transl, const vec2f& rotate);

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

}  // namespace yocto

// ---------------------------------------------------------------------------//
//                                                                            //
//                             IMPLEMENTATION                                 //
//                                                                            //
// ---------------------------------------------------------------------------//

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF MATRICES
// -----------------------------------------------------------------------------
namespace yocto {

// Matrix diagonals and transposes.
inline mat2f transpose(const mat2f& a) {
    return {{a.x.x, a.y.x}, {a.x.y, a.y.y}};
}
inline mat3f transpose(const mat3f& a) {
    return {
        {a.x.x, a.y.x, a.z.x},
        {a.x.y, a.y.y, a.z.y},
        {a.x.z, a.y.z, a.z.z},
    };
}
inline mat4f transpose(const mat4f& a) {
    return {
        {a.x.x, a.y.x, a.z.x, a.w.x},
        {a.x.y, a.y.y, a.z.y, a.w.y},
        {a.x.z, a.y.z, a.z.z, a.w.z},
        {a.x.w, a.y.w, a.z.w, a.w.w},
    };
}

// Matrix adjugates, determinant and inverses.
inline mat2f adjugate(const mat2f& a) {
    return {{a.y.y, -a.x.y}, {-a.y.x, a.x.x}};
}
inline mat3f adjugate(const mat3f& a) {
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
inline mat4f adjugate(const mat4f& a) {
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

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF UI UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Turntable for UI navigation.
inline void update_camera_turntable(vec3f& from, vec3f& to, vec3f& up,
    const vec2f& rotate, float dolly, const vec2f& pan) {
    // rotate if necessary
    if (rotate.x || rotate.y) {
        auto z     = normalize(to - from);
        auto lz    = length(to - from);
        auto phi   = atan2(z.z, z.x) + rotate.x;
        auto theta = acos(z.y) + rotate.y;
        theta      = clamp(theta, 0.001f, pi - 0.001f);
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
        theta      = clamp(theta, 0.001f, pi - 0.001f);
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
    frame3f& frame, vec3f transl, vec2f rotate) {
    // https://gamedev.stackexchange.com/questions/30644/how-to-keep-my-quaternion-using-fps-camera-from-tilting-and-messing-up
    auto y = vec3f{0, 1, 0};
    auto z = orthonormalize(frame.z, y);
    auto x = cross(y, z);

    auto rot = make_rotation_frame(vec3f{1, 0, 0}, rotate.y) *
               frame3f{frame.x, frame.y, frame.z, vec3f{0, 0, 0}} *
               make_rotation_frame(vec3f{0, 1, 0}, rotate.x);
    auto pos = frame.o + transl.x * x + transl.y * y + transl.z * z;

    frame = {rot.x, rot.y, rot.z, pos};
}

}  // namespace yocto

#endif
