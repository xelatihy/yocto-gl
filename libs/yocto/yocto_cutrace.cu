//
// # Yocto/CuTrace: Path tracing on Cuda/Optix
//
// Yocto/CuTrace is a simple path tracer written on the Yocto/Scene model.
// Yocto/CuTrace is implemented in `yocto_cutrace.h`, `yocto_cutrace.cpp`,
// and `yocto_cutrace.cu`.
// This library includes a stand-alone implementaton of the PCG32 random number
// generator by M.E. O'Neill.
//
// THIS IS AN EXPERIMENTAL LIBRARY THAT IS NOT READY FOR PRIME TIME
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2021 Fabio Pellacini
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

#include <optix_device.h>
// do not flip it
#include <cuda_runtime.h>

// HACK TO ALLOW CUT&PASTING FROM YOCTO'S CODE
#define inline __forceinline__ __device__
#define static static __forceinline__ __device__
#define optix_shader extern "C" __global__
#define optix_constant extern "C" __constant__

// whether to use builtin compound types or yocto's ones
#define CUTRACE_BUILTIN_VECS 0

// -----------------------------------------------------------------------------
// SUBSTITUTES FOR STD TYPES
// -----------------------------------------------------------------------------
namespace yocto {

// pair
template <typename T1, typename T2>
struct pair {
  T1 first;
  T2 second;
};

}  // namespace yocto

// -----------------------------------------------------------------------------
// MATH TYPES
// -----------------------------------------------------------------------------
namespace yocto {

#if CUTRACE_BUILTIN_VECS

using vec2f = float2;
using vec3f = float3;
using vec4f = float4;
using vec2i = int2;
using vec3i = int3;
using vec4i = int4;

#else

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

using byte = unsigned char;

struct vec4b {
  byte x = 0;
  byte y = 0;
  byte z = 0;
  byte w = 0;
};

#endif

// Rigid frames stored as a column-major affine transform matrix.
struct frame2f {
  vec2f x = {1, 0};
  vec2f y = {0, 1};
  vec2f o = {0, 0};
};

// Rigid frames stored as a column-major affine transform matrix.
struct frame3f {
  vec3f x = {1, 0, 0};
  vec3f y = {0, 1, 0};
  vec3f z = {0, 0, 1};
  vec3f o = {0, 0, 0};
};

}  // namespace yocto

// -----------------------------------------------------------------------------
// MATH CONSTANTS AND FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

using byte   = unsigned char;
using uint   = unsigned int;
using ushort = unsigned short;

constexpr double pi  = 3.14159265358979323846;
constexpr float  pif = (float)pi;

inline float abs(float a) { return a < 0 ? -a : a; }
inline float min(float a, float b) { return (a < b) ? a : b; }
inline float max(float a, float b) { return (a > b) ? a : b; }
inline float clamp(float a, float min_, float max_) {
  return min(max(a, min_), max_);
}
inline float sign(float a) { return a < 0 ? -1.0f : 1.0f; }
inline float sqr(float a) { return a * a; }
inline float sqrt(float a) { return std::sqrt(a); }
inline float sin(float a) { return std::sin(a); }
inline float cos(float a) { return std::cos(a); }
inline float tan(float a) { return std::tan(a); }
inline float asin(float a) { return std::asin(a); }
inline float acos(float a) { return std::acos(a); }
inline float atan(float a) { return std::atan(a); }
inline float log(float a) { return std::log(a); }
inline float exp(float a) { return std::exp(a); }
inline float log2(float a) { return std::log2(a); }
inline float exp2(float a) { return std::exp2(a); }
inline float pow(float a, float b) { return std::pow(a, b); }
inline bool  isfinite(float a) { return std::isfinite(a); }
inline float atan2(float a, float b) { return std::atan2(a, b); }
inline float fmod(float a, float b) { return std::fmod(a, b); }
inline void  swap(float& a, float& b) { std::swap(a, b); }
inline float radians(float a) { return a * pif / 180; }
inline float degrees(float a) { return a * 180 / pif; }
inline float lerp(float a, float b, float u) { return a * (1 - u) + b * u; }
inline float step(float a, float u) { return u < a ? 0.0f : 1.0f; }
inline float smoothstep(float a, float b, float u) {
  auto t = clamp((u - a) / (b - a), 0.0f, 1.0f);
  return t * t * (3 - 2 * t);
}
inline float bias(float a, float bias) {
  return a / ((1 / bias - 2) * (1 - a) + 1);
}
inline float gain(float a, float gain) {
  return (a < 0.5f) ? bias(a * 2, gain) / 2
                    : bias(a * 2 - 1, 1 - gain) / 2 + 0.5f;
}

inline int  abs(int a) { return a < 0 ? -a : a; }
inline int  min(int a, int b) { return (a < b) ? a : b; }
inline int  max(int a, int b) { return (a > b) ? a : b; }
inline int  clamp(int a, int min_, int max_) { return min(max(a, min_), max_); }
inline int  sign(int a) { return a < 0 ? -1 : 1; }
inline int  pow2(int a) { return 1 << a; }
inline void swap(int& a, int& b) { std::swap(a, b); }

inline size_t min(size_t a, size_t b) { return (a < b) ? a : b; }
inline size_t max(size_t a, size_t b) { return (a > b) ? a : b; }

}  // namespace yocto

// -----------------------------------------------------------------------------
// VECTORS
// -----------------------------------------------------------------------------
namespace yocto {

// Element access
inline vec3f xyz(const vec4f& a) { return {a.x, a.y, a.z}; }

// Vector sequence operations.
inline int          size(const vec2f& a) { return 2; }
inline const float* begin(const vec2f& a) { return &a.x; }
inline const float* end(const vec2f& a) { return &a.x + 2; }
inline float*       begin(vec2f& a) { return &a.x; }
inline float*       end(vec2f& a) { return &a.x + 2; }
inline const float* data(const vec2f& a) { return &a.x; }
inline float*       data(vec2f& a) { return &a.x; }

// Vector comparison operations.
inline bool operator==(const vec2f& a, const vec2f& b) {
  return a.x == b.x && a.y == b.y;
}
inline bool operator!=(const vec2f& a, const vec2f& b) {
  return a.x != b.x || a.y != b.y;
}

// Vector operations.
inline vec2f operator+(const vec2f& a) { return a; }
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

// Vector assignments
inline vec2f& operator+=(vec2f& a, const vec2f& b) { return a = a + b; }
inline vec2f& operator+=(vec2f& a, float b) { return a = a + b; }
inline vec2f& operator-=(vec2f& a, const vec2f& b) { return a = a - b; }
inline vec2f& operator-=(vec2f& a, float b) { return a = a - b; }
inline vec2f& operator*=(vec2f& a, const vec2f& b) { return a = a * b; }
inline vec2f& operator*=(vec2f& a, float b) { return a = a * b; }
inline vec2f& operator/=(vec2f& a, const vec2f& b) { return a = a / b; }
inline vec2f& operator/=(vec2f& a, float b) { return a = a / b; }

// Vector products and lengths.
inline float dot(const vec2f& a, const vec2f& b) {
  return a.x * b.x + a.y * b.y;
}
inline float cross(const vec2f& a, const vec2f& b) {
  return a.x * b.y - a.y * b.x;
}

inline float length(const vec2f& a) { return sqrt(dot(a, a)); }
inline float length_squared(const vec2f& a) { return dot(a, a); }
inline vec2f normalize(const vec2f& a) {
  auto l = length(a);
  return (l != 0) ? a / l : a;
}
inline float distance(const vec2f& a, const vec2f& b) { return length(a - b); }
inline float distance_squared(const vec2f& a, const vec2f& b) {
  return dot(a - b, a - b);
}
inline float angle(const vec2f& a, const vec2f& b) {
  return acos(clamp(dot(normalize(a), normalize(b)), (float)-1, (float)1));
}

// Max element and clamp.
inline vec2f max(const vec2f& a, float b) { return {max(a.x, b), max(a.y, b)}; }
inline vec2f min(const vec2f& a, float b) { return {min(a.x, b), min(a.y, b)}; }
inline vec2f max(const vec2f& a, const vec2f& b) {
  return {max(a.x, b.x), max(a.y, b.y)};
}
inline vec2f min(const vec2f& a, const vec2f& b) {
  return {min(a.x, b.x), min(a.y, b.y)};
}
inline vec2f clamp(const vec2f& x, float min, float max) {
  return {clamp(x.x, min, max), clamp(x.y, min, max)};
}
inline vec2f lerp(const vec2f& a, const vec2f& b, float u) {
  return a * (1 - u) + b * u;
}
inline vec2f lerp(const vec2f& a, const vec2f& b, const vec2f& u) {
  return a * (1 - u) + b * u;
}

inline float max(const vec2f& a) { return max(a.x, a.y); }
inline float min(const vec2f& a) { return min(a.x, a.y); }
inline float sum(const vec2f& a) { return a.x + a.y; }
inline float mean(const vec2f& a) { return sum(a) / 2; }

// Functions applied to vector elements
inline vec2f abs(const vec2f& a) { return {abs(a.x), abs(a.y)}; }
inline vec2f sqr(const vec2f& a) { return {sqr(a.x), sqr(a.y)}; }
inline vec2f sqrt(const vec2f& a) { return {sqrt(a.x), sqrt(a.y)}; }
inline vec2f exp(const vec2f& a) { return {exp(a.x), exp(a.y)}; }
inline vec2f log(const vec2f& a) { return {log(a.x), log(a.y)}; }
inline vec2f exp2(const vec2f& a) { return {exp2(a.x), exp2(a.y)}; }
inline vec2f log2(const vec2f& a) { return {log2(a.x), log2(a.y)}; }
inline bool  isfinite(const vec2f& a) { return isfinite(a.x) && isfinite(a.y); }
inline vec2f pow(const vec2f& a, float b) { return {pow(a.x, b), pow(a.y, b)}; }
inline vec2f pow(const vec2f& a, const vec2f& b) {
  return {pow(a.x, b.x), pow(a.y, b.y)};
}
inline vec2f gain(const vec2f& a, float b) {
  return {gain(a.x, b), gain(a.y, b)};
}
inline void swap(vec2f& a, vec2f& b) { std::swap(a, b); }

// Vector sequence operations.
inline int          size(const vec3f& a) { return 3; }
inline const float* begin(const vec3f& a) { return &a.x; }
inline const float* end(const vec3f& a) { return &a.x + 3; }
inline float*       begin(vec3f& a) { return &a.x; }
inline float*       end(vec3f& a) { return &a.x + 3; }
inline const float* data(const vec3f& a) { return &a.x; }
inline float*       data(vec3f& a) { return &a.x; }

// Vector comparison operations.
inline bool operator==(const vec3f& a, const vec3f& b) {
  return a.x == b.x && a.y == b.y && a.z == b.z;
}
inline bool operator!=(const vec3f& a, const vec3f& b) {
  return a.x != b.x || a.y != b.y || a.z != b.z;
}

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

// Vector assignments
inline vec3f& operator+=(vec3f& a, const vec3f& b) { return a = a + b; }
inline vec3f& operator+=(vec3f& a, float b) { return a = a + b; }
inline vec3f& operator-=(vec3f& a, const vec3f& b) { return a = a - b; }
inline vec3f& operator-=(vec3f& a, float b) { return a = a - b; }
inline vec3f& operator*=(vec3f& a, const vec3f& b) { return a = a * b; }
inline vec3f& operator*=(vec3f& a, float b) { return a = a * b; }
inline vec3f& operator/=(vec3f& a, const vec3f& b) { return a = a / b; }
inline vec3f& operator/=(vec3f& a, float b) { return a = a / b; }

// Vector products and lengths.
inline float dot(const vec3f& a, const vec3f& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline vec3f cross(const vec3f& a, const vec3f& b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}

inline float length(const vec3f& a) { return sqrt(dot(a, a)); }
inline float length_squared(const vec3f& a) { return dot(a, a); }
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
inline vec3f refract(const vec3f& w, const vec3f& n, float inv_eta) {
  auto cosine = dot(n, w);
  auto k      = 1 + inv_eta * inv_eta * (cosine * cosine - 1);
  if (k < 0) return {0, 0, 0};  // tir
  return -w * inv_eta + (inv_eta * cosine - sqrt(k)) * n;
}

// Max element and clamp.
inline vec3f max(const vec3f& a, float b) {
  return {max(a.x, b), max(a.y, b), max(a.z, b)};
}
inline vec3f min(const vec3f& a, float b) {
  return {min(a.x, b), min(a.y, b), min(a.z, b)};
}
inline vec3f max(const vec3f& a, const vec3f& b) {
  return {max(a.x, b.x), max(a.y, b.y), max(a.z, b.z)};
}
inline vec3f min(const vec3f& a, const vec3f& b) {
  return {min(a.x, b.x), min(a.y, b.y), min(a.z, b.z)};
}
inline vec3f clamp(const vec3f& x, float min, float max) {
  return {clamp(x.x, min, max), clamp(x.y, min, max), clamp(x.z, min, max)};
}
inline vec3f lerp(const vec3f& a, const vec3f& b, float u) {
  return a * (1 - u) + b * u;
}
inline vec3f lerp(const vec3f& a, const vec3f& b, const vec3f& u) {
  return a * (1 - u) + b * u;
}

inline float max(const vec3f& a) { return max(max(a.x, a.y), a.z); }
inline float min(const vec3f& a) { return min(min(a.x, a.y), a.z); }
inline float sum(const vec3f& a) { return a.x + a.y + a.z; }
inline float mean(const vec3f& a) { return sum(a) / 3; }

// Functions applied to vector elements
inline vec3f abs(const vec3f& a) { return {abs(a.x), abs(a.y), abs(a.z)}; }
inline vec3f sqr(const vec3f& a) { return {sqr(a.x), sqr(a.y), sqr(a.z)}; }
inline vec3f sqrt(const vec3f& a) { return {sqrt(a.x), sqrt(a.y), sqrt(a.z)}; }
inline vec3f exp(const vec3f& a) { return {exp(a.x), exp(a.y), exp(a.z)}; }
inline vec3f log(const vec3f& a) { return {log(a.x), log(a.y), log(a.z)}; }
inline vec3f exp2(const vec3f& a) { return {exp2(a.x), exp2(a.y), exp2(a.z)}; }
inline vec3f log2(const vec3f& a) { return {log2(a.x), log2(a.y), log2(a.z)}; }
inline vec3f pow(const vec3f& a, float b) {
  return {pow(a.x, b), pow(a.y, b), pow(a.z, b)};
}
inline vec3f pow(const vec3f& a, const vec3f& b) {
  return {pow(a.x, b.x), pow(a.y, b.y), pow(a.z, b.z)};
}
inline vec3f gain(const vec3f& a, float b) {
  return {gain(a.x, b), gain(a.y, b), gain(a.z, b)};
}
inline bool isfinite(const vec3f& a) {
  return isfinite(a.x) && isfinite(a.y) && isfinite(a.z);
}
inline void swap(vec3f& a, vec3f& b) { std::swap(a, b); }

// Vector sequence operations.
inline int          size(const vec4f& a) { return 4; }
inline const float* begin(const vec4f& a) { return &a.x; }
inline const float* end(const vec4f& a) { return &a.x + 4; }
inline float*       begin(vec4f& a) { return &a.x; }
inline float*       end(vec4f& a) { return &a.x + 4; }
inline const float* data(const vec4f& a) { return &a.x; }
inline float*       data(vec4f& a) { return &a.x; }

// Vector comparison operations.
inline bool operator==(const vec4f& a, const vec4f& b) {
  return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}
inline bool operator!=(const vec4f& a, const vec4f& b) {
  return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w;
}

// Vector operations.
inline vec4f operator+(const vec4f& a) { return a; }
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
inline vec4f& operator+=(vec4f& a, const vec4f& b) { return a = a + b; }
inline vec4f& operator+=(vec4f& a, float b) { return a = a + b; }
inline vec4f& operator-=(vec4f& a, const vec4f& b) { return a = a - b; }
inline vec4f& operator-=(vec4f& a, float b) { return a = a - b; }
inline vec4f& operator*=(vec4f& a, const vec4f& b) { return a = a * b; }
inline vec4f& operator*=(vec4f& a, float b) { return a = a * b; }
inline vec4f& operator/=(vec4f& a, const vec4f& b) { return a = a / b; }
inline vec4f& operator/=(vec4f& a, float b) { return a = a / b; }

// Vector products and lengths.
inline float dot(const vec4f& a, const vec4f& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}
inline float length(const vec4f& a) { return sqrt(dot(a, a)); }
inline float length_squared(const vec4f& a) { return dot(a, a); }
inline vec4f normalize(const vec4f& a) {
  auto l = length(a);
  return (l != 0) ? a / l : a;
}
inline float distance(const vec4f& a, const vec4f& b) { return length(a - b); }
inline float distance_squared(const vec4f& a, const vec4f& b) {
  return dot(a - b, a - b);
}
inline float angle(const vec4f& a, const vec4f& b) {
  return acos(clamp(dot(normalize(a), normalize(b)), (float)-1, (float)1));
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
  if (th == 0) return an;
  return an * (sin(th * (1 - u)) / sin(th)) + bn * (sin(th * u) / sin(th));
}

// Max element and clamp.
inline vec4f max(const vec4f& a, float b) {
  return {max(a.x, b), max(a.y, b), max(a.z, b), max(a.w, b)};
}
inline vec4f min(const vec4f& a, float b) {
  return {min(a.x, b), min(a.y, b), min(a.z, b), min(a.w, b)};
}
inline vec4f max(const vec4f& a, const vec4f& b) {
  return {max(a.x, b.x), max(a.y, b.y), max(a.z, b.z), max(a.w, b.w)};
}
inline vec4f min(const vec4f& a, const vec4f& b) {
  return {min(a.x, b.x), min(a.y, b.y), min(a.z, b.z), min(a.w, b.w)};
}
inline vec4f clamp(const vec4f& x, float min, float max) {
  return {clamp(x.x, min, max), clamp(x.y, min, max), clamp(x.z, min, max),
      clamp(x.w, min, max)};
}
inline vec4f lerp(const vec4f& a, const vec4f& b, float u) {
  return a * (1 - u) + b * u;
}
inline vec4f lerp(const vec4f& a, const vec4f& b, const vec4f& u) {
  return a * (1 - u) + b * u;
}

inline float max(const vec4f& a) { return max(max(max(a.x, a.y), a.z), a.w); }
inline float min(const vec4f& a) { return min(min(min(a.x, a.y), a.z), a.w); }
inline float sum(const vec4f& a) { return a.x + a.y + a.z + a.w; }
inline float mean(const vec4f& a) { return sum(a) / 4; }

// Functions applied to vector elements
inline vec4f abs(const vec4f& a) {
  return {abs(a.x), abs(a.y), abs(a.z), abs(a.w)};
}
inline vec4f sqr(const vec4f& a) {
  return {sqr(a.x), sqr(a.y), sqr(a.z), sqr(a.w)};
}
inline vec4f sqrt(const vec4f& a) {
  return {sqrt(a.x), sqrt(a.y), sqrt(a.z), sqrt(a.w)};
}
inline vec4f exp(const vec4f& a) {
  return {exp(a.x), exp(a.y), exp(a.z), exp(a.w)};
}
inline vec4f log(const vec4f& a) {
  return {log(a.x), log(a.y), log(a.z), log(a.w)};
}
inline vec4f exp2(const vec4f& a) {
  return {exp2(a.x), exp2(a.y), exp2(a.z), exp2(a.w)};
}
inline vec4f log2(const vec4f& a) {
  return {log2(a.x), log2(a.y), log2(a.z), log2(a.w)};
}
inline vec4f pow(const vec4f& a, float b) {
  return {pow(a.x, b), pow(a.y, b), pow(a.z, b), pow(a.w, b)};
}
inline vec4f pow(const vec4f& a, const vec4f& b) {
  return {pow(a.x, b.x), pow(a.y, b.y), pow(a.z, b.z), pow(a.w, b.w)};
}
inline vec4f gain(const vec4f& a, float b) {
  return {gain(a.x, b), gain(a.y, b), gain(a.z, b), gain(a.w, b)};
}
inline bool isfinite(const vec4f& a) {
  return isfinite(a.x) && isfinite(a.y) && isfinite(a.z) && isfinite(a.w);
}
inline void swap(vec4f& a, vec4f& b) { std::swap(a, b); }

// Quaternion operatons represented as xi + yj + zk + w
// const auto identity_quat4f = vec4f{0, 0, 0, 1};
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
// TRANSFORMS
// -----------------------------------------------------------------------------
namespace yocto {

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
    const frame2f& a, const vec2f& b, bool non_rigid) {
  // if (non_rigid) {
  //  return transform_normal(rotation(a), b);
  //} else {
  return normalize(transform_vector(a, b));
  //}
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
  // if (non_rigid) {
  //   return transform_normal(rotation(a), b);
  // } else {
  return normalize(transform_vector(a, b));
  //}
}

// Translation, scaling and rotations transforms.
inline frame3f translation_frame(const vec3f& a) {
  return {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, a};
}
inline frame3f scaling_frame(const vec3f& a) {
  return {{a.x, 0, 0}, {0, a.y, 0}, {0, 0, a.z}, {0, 0, 0}};
}
inline frame3f rotation_frame(const vec3f& axis, float angle) {
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

// Lookat frame. Z-axis can be inverted with inv_xz.
inline frame3f lookat_frame(
    const vec3f& eye, const vec3f& center, const vec3f& up, bool inv_xz) {
  auto w = normalize(eye - center);
  auto u = normalize(cross(up, w));
  auto v = normalize(cross(w, u));
  if (inv_xz) {
    w = -w;
    u = -u;
  }
  return {u, v, w, eye};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// GEOMETRY TYPES
// -----------------------------------------------------------------------------
namespace yocto {

// Ray epsilon
constexpr auto ray_eps = 1e-4f;
constexpr auto ray_far = 1e20f;

struct ray2f {
  vec2f o    = {0, 0};
  vec2f d    = {0, 1};
  float tmin = ray_eps;
  float tmax = ray_far;
};

// Rays with origin, direction and min/max t value.
struct ray3f {
  vec3f o    = {0, 0, 0};
  vec3f d    = {0, 0, 1};
  float tmin = ray_eps;
  float tmax = ray_far;
};

}  // namespace yocto

// -----------------------------------------------------------------------------
// GEOMETRY UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Line properties.
inline vec3f line_tangent(const vec3f& p0, const vec3f& p1) {
  return normalize(p1 - p0);
}
inline float line_length(const vec3f& p0, const vec3f& p1) {
  return length(p1 - p0);
}

// Triangle properties.
inline vec3f triangle_normal(
    const vec3f& p0, const vec3f& p1, const vec3f& p2) {
  return normalize(cross(p1 - p0, p2 - p0));
}
inline float triangle_area(const vec3f& p0, const vec3f& p1, const vec3f& p2) {
  return length(cross(p1 - p0, p2 - p0)) / 2;
}

// Quad propeties.
inline vec3f quad_normal(
    const vec3f& p0, const vec3f& p1, const vec3f& p2, const vec3f& p3) {
  return normalize(triangle_normal(p0, p1, p3) + triangle_normal(p2, p3, p1));
}
inline float quad_area(
    const vec3f& p0, const vec3f& p1, const vec3f& p2, const vec3f& p3) {
  return triangle_area(p0, p1, p3) + triangle_area(p2, p3, p1);
}

// Interpolates values over a line parameterized from a to b by u. Same as lerp.
template <typename T>
inline T interpolate_line(const T& p0, const T& p1, float u) {
  return p0 * (1 - u) + p1 * u;
}
// Interpolates values over a triangle parameterized by u and v along the
// (p1-p0) and (p2-p0) directions. Same as barycentric interpolation.
template <typename T>
inline T interpolate_triangle(
    const T& p0, const T& p1, const T& p2, const vec2f& uv) {
  return p0 * (1 - uv.x - uv.y) + p1 * uv.x + p2 * uv.y;
}
// Interpolates values over a quad parameterized by u and v along the
// (p1-p0) and (p2-p1) directions. Same as bilinear interpolation.
template <typename T>
inline T interpolate_quad(
    const T& p0, const T& p1, const T& p2, const T& p3, const vec2f& uv) {
  if (uv.x + uv.y <= 1) {
    return interpolate_triangle(p0, p1, p3, uv);
  } else {
    return interpolate_triangle(p2, p3, p1, 1 - uv);
  }
}

// Interpolates values along a cubic Bezier segment parametrized by u.
template <typename T>
inline T interpolate_bezier(
    const T& p0, const T& p1, const T& p2, const T& p3, float u) {
  return p0 * (1 - u) * (1 - u) * (1 - u) + p1 * 3 * u * (1 - u) * (1 - u) +
         p2 * 3 * u * u * (1 - u) + p3 * u * u * u;
}
// Computes the derivative of a cubic Bezier segment parametrized by u.
template <typename T>
inline T interpolate_bezier_derivative(
    const T& p0, const T& p1, const T& p2, const T& p3, float u) {
  return (p1 - p0) * 3 * (1 - u) * (1 - u) + (p2 - p1) * 6 * u * (1 - u) +
         (p3 - p2) * 3 * u * u;
}

// Interpolated line properties.
inline vec3f line_point(const vec3f& p0, const vec3f& p1, float u) {
  return p0 * (1 - u) + p1 * u;
}
inline vec3f line_tangent(const vec3f& t0, const vec3f& t1, float u) {
  return normalize(t0 * (1 - u) + t1 * u);
}

// Interpolated triangle properties.
inline vec3f triangle_point(
    const vec3f& p0, const vec3f& p1, const vec3f& p2, const vec2f& uv) {
  return p0 * (1 - uv.x - uv.y) + p1 * uv.x + p2 * uv.y;
}
inline vec3f triangle_normal(
    const vec3f& n0, const vec3f& n1, const vec3f& n2, const vec2f& uv) {
  return normalize(n0 * (1 - uv.x - uv.y) + n1 * uv.x + n2 * uv.y);
}

// Interpolated quad properties.
inline vec3f quad_point(const vec3f& p0, const vec3f& p1, const vec3f& p2,
    const vec3f& p3, const vec2f& uv) {
  if (uv.x + uv.y <= 1) {
    return triangle_point(p0, p1, p3, uv);
  } else {
    return triangle_point(p2, p3, p1, 1 - uv);
  }
}
inline vec3f quad_normal(const vec3f& n0, const vec3f& n1, const vec3f& n2,
    const vec3f& n3, const vec2f& uv) {
  if (uv.x + uv.y <= 1) {
    return triangle_normal(n0, n1, n3, uv);
  } else {
    return triangle_normal(n2, n3, n1, 1 - uv);
  }
}

// Interpolated sphere properties.
inline vec3f sphere_point(const vec3f p, float r, const vec2f& uv) {
  return p + r * vec3f{cos(uv.x * 2 * pif) * sin(uv.y * pif),
                     sin(uv.x * 2 * pif) * sin(uv.y * pif), cos(uv.y * pif)};
}
inline vec3f sphere_normal(const vec3f p, float r, const vec2f& uv) {
  return normalize(vec3f{cos(uv.x * 2 * pif) * sin(uv.y * pif),
      sin(uv.x * 2 * pif) * sin(uv.y * pif), cos(uv.y * pif)});
}

// Triangle tangent and bitangent from uv
inline pair<vec3f, vec3f> triangle_tangents_fromuv(const vec3f& p0,
    const vec3f& p1, const vec3f& p2, const vec2f& uv0, const vec2f& uv1,
    const vec2f& uv2) {
  // Follows the definition in http://www.terathon.com/code/tangent.html and
  // https://gist.github.com/aras-p/2843984
  // normal points up from texture space
  auto p   = p1 - p0;
  auto q   = p2 - p0;
  auto s   = vec2f{uv1.x - uv0.x, uv2.x - uv0.x};
  auto t   = vec2f{uv1.y - uv0.y, uv2.y - uv0.y};
  auto div = s.x * t.y - s.y * t.x;

  if (div != 0) {
    auto tu = vec3f{t.y * p.x - t.x * q.x, t.y * p.y - t.x * q.y,
                  t.y * p.z - t.x * q.z} /
              div;
    auto tv = vec3f{s.x * q.x - s.y * p.x, s.x * q.y - s.y * p.y,
                  s.x * q.z - s.y * p.z} /
              div;
    return {tu, tv};
  } else {
    return {{1, 0, 0}, {0, 1, 0}};
  }
}

// Quad tangent and bitangent from uv.
inline pair<vec3f, vec3f> quad_tangents_fromuv(const vec3f& p0, const vec3f& p1,
    const vec3f& p2, const vec3f& p3, const vec2f& uv0, const vec2f& uv1,
    const vec2f& uv2, const vec2f& uv3, const vec2f& current_uv) {
  if (current_uv.x + current_uv.y <= 1) {
    return triangle_tangents_fromuv(p0, p1, p3, uv0, uv1, uv3);
  } else {
    return triangle_tangents_fromuv(p2, p3, p1, uv2, uv3, uv1);
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

  rng_state() = default;
  rng_state(uint64_t state, uint64_t inc);
};

// PCG random numbers from http://www.pcg-random.org/
inline rng_state::rng_state(uint64_t state, uint64_t inc)
    : state{state}, inc{inc} {}

// Next random number, used internally only.
inline uint32_t _advance_rng(rng_state& rng) {
  uint64_t oldstate = rng.state;
  rng.state         = oldstate * 6364136223846793005ULL + rng.inc;
  auto xorshifted   = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
  auto rot          = (uint32_t)(oldstate >> 59u);
  // return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
  return (xorshifted >> rot) | (xorshifted << ((~rot + 1u) & 31));
}

// Init a random number generator with a state state from the sequence seq.
inline rng_state make_rng(uint64_t seed, uint64_t seq) {
  auto rng  = rng_state();
  rng.state = 0U;
  rng.inc   = (seq << 1u) | 1u;
  _advance_rng(rng);
  rng.state += seed;
  _advance_rng(rng);
  return rng;
}

// Next random numbers: floats in [0,1), ints in [0,n).
inline int   rand1i(rng_state& rng, int n) { return _advance_rng(rng) % n; }
inline float rand1f(rng_state& rng) {
  union {
    uint32_t u;
    float    f;
  } x;
  x.u = (_advance_rng(rng) >> 9) | 0x3f800000u;
  return x.f - 1.0f;
  // alternate implementation
  // const static auto scale = (float)(1.0 / numeric_limits<uint32_t>::max());
  // return advance_rng(rng) * scale;
}
inline vec2f rand2f(rng_state& rng) {
  // force order of evaluation by using separate assignments.
  auto x = rand1f(rng);
  auto y = rand1f(rng);
  return {x, y};
}
inline vec3f rand3f(rng_state& rng) {
  // force order of evaluation by using separate assignments.
  auto x = rand1f(rng);
  auto y = rand1f(rng);
  auto z = rand1f(rng);
  return {x, y, z};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// CUDA HELPERS
// -----------------------------------------------------------------------------
namespace yocto {

template <typename T>
struct cubuffer {
  inline bool     empty() const { return _size == 0; }
  inline size_t   size() const { return _size; }
  inline T&       operator[](int idx) { return _data[idx]; }
  inline const T& operator[](int idx) const { return _data[idx]; }

  T*     _data = nullptr;
  size_t _size = 0;
};

inline void* unpackPointer(uint32_t i0, uint32_t i1) {
  const uint64_t uptr = static_cast<uint64_t>(i0) << 32 | i1;
  void*          ptr  = reinterpret_cast<void*>(uptr);
  return ptr;
}

inline void packPointer(void* ptr, uint32_t& i0, uint32_t& i1) {
  const uint64_t uptr = reinterpret_cast<uint64_t>(ptr);
  i0                  = uptr >> 32;
  i1                  = uptr & 0x00000000ffffffff;
}

template <typename T>
inline T* getPRD() {
  const uint32_t u0 = optixGetPayload_0();
  const uint32_t u1 = optixGetPayload_1();
  return reinterpret_cast<T*>(unpackPointer(u0, u1));
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// CUTRACE TYPES
// -----------------------------------------------------------------------------
namespace yocto {

constexpr int invalidid = -1;

struct cutrace_state {
  int                 width   = 0;
  int                 height  = 0;
  int                 samples = 0;
  cubuffer<vec4f>     image   = {};
  cubuffer<vec3f>     albedo  = {};
  cubuffer<vec3f>     normal  = {};
  cubuffer<int>       hits    = {};
  cubuffer<rng_state> rngs    = {};
  cubuffer<vec4f>     display = {};
};

struct cutrace_camera {
  frame3f frame        = {};
  float   lens         = {};
  float   film         = {};
  float   aspect       = {};
  float   focus        = {};
  float   aperture     = {};
  bool    orthographic = {};
};

struct cutrace_texture {
  cudaArray_t         array   = nullptr;
  cudaTextureObject_t texture = 0;
};

enum struct material_type {
  // clang-format off
  matte, glossy, reflective, transparent, refractive, subsurface, volumetric, 
  gltfpbr
  // clang-format on
};

struct cutrace_material {
  material_type type         = material_type::matte;
  vec3f         emission     = {0, 0, 0};
  vec3f         color        = {0, 0, 0};
  float         roughness    = 0;
  float         metallic     = 0;
  float         ior          = 1.5f;
  vec3f         scattering   = {0, 0, 0};
  float         scanisotropy = 0;
  float         trdepth      = 0.01f;
  float         opacity      = 1;

  int emission_tex   = invalidid;
  int color_tex      = invalidid;
  int roughness_tex  = invalidid;
  int scattering_tex = invalidid;
  int normal_tex     = invalidid;
};

struct cutrace_instance {
  frame3f frame    = {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};
  int     shape    = -1;
  int     material = -1;
};

struct cutrace_shape {
  cubuffer<vec3f> positions = {};
  cubuffer<vec3f> normals   = {};
  cubuffer<vec2f> texcoords = {};
  cubuffer<vec4f> colors    = {};
  cubuffer<vec3i> triangles = {};
};

struct cutrace_scene {
  cubuffer<cutrace_camera>   cameras   = {};
  cubuffer<cutrace_texture>  textures  = {};
  cubuffer<cutrace_material> materials = {};
  cubuffer<cutrace_shape>    shapes    = {};
  cubuffer<cutrace_instance> instances = {};
};

struct cutrace_lights {};

// Type of tracing algorithm
enum struct cutrace_sampler_type {
  path,        // path tracing
  pathdirect,  // path tracing with direct
  pathmis,     // path tracing with mis
  naive,       // naive path tracing
  eyelight,    // eyelight rendering
  eyelightao,  // eyelight with ambient occlusion
  furnace,     // furnace test
  falsecolor,  // false color rendering
};
// Type of false color visualization
enum struct cutrace_falsecolor_type {
  // clang-format off
  position, normal, frontfacing, gnormal, gfrontfacing, texcoord, mtype, color,
  emission, roughness, opacity, metallic, delta, instance, shape, material, 
  element, highlight
  // clang-format on
};

// Default trace seed
constexpr auto cutrace_default_seed = 961748941ull;

// params
struct cutrace_params {
  int                     camera         = 0;
  int                     resolution     = 1280;
  cutrace_sampler_type    sampler        = cutrace_sampler_type::path;
  cutrace_falsecolor_type falsecolor     = cutrace_falsecolor_type::color;
  int                     samples        = 512;
  int                     bounces        = 8;
  float                   clamp          = 10;
  bool                    nocaustics     = false;
  bool                    envhidden      = false;
  bool                    tentfilter     = false;
  uint64_t                seed           = cutrace_default_seed;
  bool                    embreebvh      = false;
  bool                    highqualitybvh = false;
  bool                    noparallel     = false;
  int                     pratio         = 8;
  float                   exposure       = 0;
  bool                    filmic         = false;
  bool                    denoise        = false;
  int                     batch          = 1;
};

using cutrace_bvh = OptixTraversableHandle;

struct cutrace_globals {
  cutrace_state          state  = {};
  cutrace_scene          scene  = {};
  OptixTraversableHandle bvh    = 0;
  cutrace_params         params = {};
};

// global data
optix_constant cutrace_globals globals;

}  // namespace yocto

// -----------------------------------------------------------------------------
// SCENE FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// compatibility aliases
using scene_data    = cutrace_scene;
using camera_data   = cutrace_camera;
using material_data = cutrace_material;
using texture_data  = cutrace_texture;
using instance_data = cutrace_instance;
using shape_data    = cutrace_shape;

// constant values
constexpr auto min_roughness = 0.03f * 0.03f;

// Evaluates an image at a point `uv`.
static vec4f eval_texture(const texture_data& texture, const vec2f& texcoord,
    bool as_linear, bool no_interpolation = false, bool clamp_to_edge = false) {
  auto fromTexture = tex2D<float4>(texture.texture, texcoord.x, texcoord.y);
  return {fromTexture.x, fromTexture.y, fromTexture.z, fromTexture.w};
}

// Helpers
static vec4f eval_texture(const scene_data& scene, int texture, const vec2f& uv,
    bool ldr_as_linear, bool no_interpolation = false,
    bool clamp_to_edge = false) {
  if (texture == invalidid) return {1, 1, 1, 1};
  return eval_texture(
      scene.textures[texture], uv, ldr_as_linear, no_interpolation);
}

// Material parameters evaluated at a point on the surface
struct material_point {
  material_type type         = material_type::gltfpbr;
  vec3f         emission     = {0, 0, 0};
  vec3f         color        = {0, 0, 0};
  float         opacity      = 1;
  float         roughness    = 0;
  float         metallic     = 0;
  float         ior          = 1;
  vec3f         density      = {0, 0, 0};
  vec3f         scattering   = {0, 0, 0};
  float         scanisotropy = 0;
  float         trdepth      = 0.01f;
};

// Evaluate material
static material_point eval_material(const scene_data& scene,
    const material_data& material, const vec2f& texcoord,
    const vec4f& color_shp) {
  // evaluate textures
  auto emission_tex = eval_texture(
      scene, material.emission_tex, texcoord, true);
  auto color_tex     = eval_texture(scene, material.color_tex, texcoord, true);
  auto roughness_tex = eval_texture(
      scene, material.roughness_tex, texcoord, false);
  auto scattering_tex = eval_texture(
      scene, material.scattering_tex, texcoord, true);

  // material point
  auto point         = material_point{};
  point.type         = material.type;
  point.emission     = material.emission * xyz(emission_tex);
  point.color        = material.color * xyz(color_tex) * xyz(color_shp);
  point.opacity      = material.opacity * color_tex.w * color_shp.w;
  point.metallic     = material.metallic * roughness_tex.z;
  point.roughness    = material.roughness * roughness_tex.y;
  point.roughness    = point.roughness * point.roughness;
  point.ior          = material.ior;
  point.scattering   = material.scattering * xyz(scattering_tex);
  point.scanisotropy = material.scanisotropy;
  point.trdepth      = material.trdepth;

  // volume density
  if (material.type == material_type::refractive ||
      material.type == material_type::volumetric ||
      material.type == material_type::subsurface) {
    point.density = -log(clamp(point.color, 0.0001f, 1.0f)) / point.trdepth;
  } else {
    point.density = {0, 0, 0};
  }

  // fix roughness
  if (point.type == material_type::matte ||
      point.type == material_type::gltfpbr ||
      point.type == material_type::glossy) {
    point.roughness = clamp(point.roughness, min_roughness, 1.0f);
  }

  return point;
}

// Eval position
static vec3f eval_position(const scene_data& scene,
    const instance_data& instance, int element, const vec2f& uv) {
  auto& shape = scene.shapes[instance.shape];
  if (!shape.triangles.empty()) {
    auto t = shape.triangles[element];
    return transform_point(
        instance.frame, interpolate_triangle(shape.positions[t.x],
                            shape.positions[t.y], shape.positions[t.z], uv));
  } else {
    return {0, 0, 0};
  }
}

// Shape element normal.
static vec3f eval_element_normal(
    const scene_data& scene, const instance_data& instance, int element) {
  auto& shape = scene.shapes[instance.shape];
  if (!shape.triangles.empty()) {
    auto t = shape.triangles[element];
    return transform_normal(
        instance.frame, triangle_normal(shape.positions[t.x],
                            shape.positions[t.y], shape.positions[t.z]));
  } else {
    return {0, 0, 0};
  }
}

// Eval normal
static vec3f eval_normal(const scene_data& scene, const instance_data& instance,
    int element, const vec2f& uv) {
  auto& shape = scene.shapes[instance.shape];
  if (shape.normals.empty())
    return eval_element_normal(scene, instance, element);
  if (!shape.triangles.empty()) {
    auto t = shape.triangles[element];
    return transform_normal(
        instance.frame, normalize(interpolate_triangle(shape.normals[t.x],
                            shape.normals[t.y], shape.normals[t.z], uv)));
  } else {
    return {0, 0, 0};
  }
}

// Eval texcoord
static vec2f eval_texcoord(const scene_data& scene,
    const instance_data& instance, int element, const vec2f& uv) {
  auto& shape = scene.shapes[instance.shape];
  if (shape.texcoords.empty()) return uv;
  if (!shape.triangles.empty()) {
    auto t = shape.triangles[element];
    return interpolate_triangle(
        shape.texcoords[t.x], shape.texcoords[t.y], shape.texcoords[t.z], uv);
  } else {
    return {0, 0};
  }
}

// Shape element normal.
static pair<vec3f, vec3f> eval_element_tangents(
    const scene_data& scene, const instance_data& instance, int element) {
  auto& shape = scene.shapes[instance.shape];
  if (!shape.triangles.empty() && !shape.texcoords.empty()) {
    auto t   = shape.triangles[element];
    auto tuv = triangle_tangents_fromuv(shape.positions[t.x],
        shape.positions[t.y], shape.positions[t.z], shape.texcoords[t.x],
        shape.texcoords[t.y], shape.texcoords[t.z]);
    return {transform_direction(instance.frame, tuv.first),
        transform_direction(instance.frame, tuv.second)};
  } else {
    return {};
  }
}

static vec3f eval_normalmap(const scene_data& scene,
    const instance_data& instance, int element, const vec2f& uv) {
  auto& shape    = scene.shapes[instance.shape];
  auto& material = scene.materials[instance.material];
  // apply normal mapping
  auto normal   = eval_normal(scene, instance, element, uv);
  auto texcoord = eval_texcoord(scene, instance, element, uv);
  if (material.normal_tex != invalidid && (!shape.triangles.empty())) {
    auto& normal_tex = scene.textures[material.normal_tex];
    auto  normalmap  = -1 + 2 * xyz(eval_texture(normal_tex, texcoord, false));
    auto  tuv        = eval_element_tangents(scene, instance, element);
    auto  frame      = frame3f{tuv.first, tuv.second, normal, {0, 0, 0}};
    frame.x          = orthonormalize(frame.x, frame.z);
    frame.y          = normalize(cross(frame.z, frame.x));
    auto flip_v      = dot(frame.y, tuv.second) < 0;
    normalmap.y *= flip_v ? 1 : -1;  // flip vertical axis
    normal = transform_normal(frame, normalmap);
  }
  return normal;
}

// Eval shading position
static vec3f eval_shading_position(const scene_data& scene,
    const instance_data& instance, int element, const vec2f& uv,
    const vec3f& outgoing) {
  auto& shape = scene.shapes[instance.shape];
  if (!shape.triangles.empty()) {
    return eval_position(scene, instance, element, uv);
  } else {
    return {0, 0, 0};
  }
}

// Eval shading normal
static vec3f eval_shading_normal(const scene_data& scene,
    const instance_data& instance, int element, const vec2f& uv,
    const vec3f& outgoing) {
  auto& shape    = scene.shapes[instance.shape];
  auto& material = scene.materials[instance.material];
  if (!shape.triangles.empty()) {
    auto normal = eval_normal(scene, instance, element, uv);
    if (material.normal_tex != invalidid) {
      normal = eval_normalmap(scene, instance, element, uv);
    }
    if (material.type == material_type::refractive) return normal;
    return dot(normal, outgoing) >= 0 ? normal : -normal;
  } else {
    return {0, 0, 0};
  }
}

// Eval color
static vec4f eval_color(const scene_data& scene, const instance_data& instance,
    int element, const vec2f& uv) {
  auto& shape = scene.shapes[instance.shape];
  if (shape.colors.empty()) return {1, 1, 1, 1};
  if (!shape.triangles.empty()) {
    auto t = shape.triangles[element];
    return interpolate_triangle(
        shape.colors[t.x], shape.colors[t.y], shape.colors[t.z], uv);
  } else {
    return {0, 0, 0, 0};
  }
}

// Evaluate material
static material_point eval_material(const scene_data& scene,
    const instance_data& instance, int element, const vec2f& uv) {
  auto& material = scene.materials[instance.material];
  auto  texcoord = eval_texcoord(scene, instance, element, uv);

  // evaluate textures
  auto emission_tex = eval_texture(
      scene, material.emission_tex, texcoord, true);
  auto color_shp     = eval_color(scene, instance, element, uv);
  auto color_tex     = eval_texture(scene, material.color_tex, texcoord, true);
  auto roughness_tex = eval_texture(
      scene, material.roughness_tex, texcoord, false);
  auto scattering_tex = eval_texture(
      scene, material.scattering_tex, texcoord, true);

  // material point
  auto point         = material_point{};
  point.type         = material.type;
  point.emission     = material.emission * xyz(emission_tex);
  point.color        = material.color * xyz(color_tex) * xyz(color_shp);
  point.opacity      = material.opacity * color_tex.w * color_shp.w;
  point.metallic     = material.metallic * roughness_tex.z;
  point.roughness    = material.roughness * roughness_tex.y;
  point.roughness    = point.roughness * point.roughness;
  point.ior          = material.ior;
  point.scattering   = material.scattering * xyz(scattering_tex);
  point.scanisotropy = material.scanisotropy;
  point.trdepth      = material.trdepth;

  // volume density
  if (material.type == material_type::refractive ||
      material.type == material_type::volumetric ||
      material.type == material_type::subsurface) {
    point.density = -log(clamp(point.color, 0.0001f, 1.0f)) / point.trdepth;
  } else {
    point.density = {0, 0, 0};
  }

  // fix roughness
  if (point.type == material_type::matte ||
      point.type == material_type::gltfpbr ||
      point.type == material_type::glossy) {
    point.roughness = clamp(point.roughness, min_roughness, 1.0f);
  } else if (material.type == material_type::volumetric) {
    point.roughness = 0;
  } else {
    if (point.roughness < min_roughness) point.roughness = 0;
  }

  return point;
}

// check if a material is a delta or volumetric
static bool is_delta(const material_data& material) {
  return (material.type == material_type::reflective &&
             material.roughness == 0) ||
         (material.type == material_type::refractive &&
             material.roughness == 0) ||
         (material.type == material_type::transparent &&
             material.roughness == 0) ||
         (material.type == material_type::volumetric);
}
static bool is_volumetric(const material_data& material) {
  return material.type == material_type::refractive ||
         material.type == material_type::volumetric ||
         material.type == material_type::subsurface;
}

// check if an instance is volumetric
static bool is_volumetric(
    const scene_data& scene, const instance_data& instance) {
  return is_volumetric(scene.materials[instance.material]);
}

// check if a brdf is a delta
static bool is_delta(const material_point& material) {
  return (material.type == material_type::reflective &&
             material.roughness == 0) ||
         (material.type == material_type::refractive &&
             material.roughness == 0) ||
         (material.type == material_type::transparent &&
             material.roughness == 0) ||
         (material.type == material_type::volumetric);
}
static bool has_volume(const material_point& material) {
  return material.type == material_type::refractive ||
         material.type == material_type::volumetric ||
         material.type == material_type::subsurface;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// RAY-SCENE INTERSECTION
// -----------------------------------------------------------------------------
namespace yocto {

// intersection result
struct scene_intersection {
  int   instance = -1;
  int   element  = -1;
  vec2f uv       = {0, 0};
  float distance = 0;
  bool  hit      = false;
  float _pad     = 0;
};

// closest hit
optix_shader void __closesthit__intersect_scene() {
  auto& intersection    = *getPRD<scene_intersection>();
  intersection.instance = optixGetInstanceIndex();
  intersection.element  = optixGetPrimitiveIndex();
  intersection.uv       = {
      optixGetTriangleBarycentrics().x, optixGetTriangleBarycentrics().y};
  intersection.distance = 0;
  intersection.hit      = true;
}

// anyhit shader
optix_shader void __anyhit__intersect_scene() {}

// miss shader
optix_shader void __miss__intersect_scene() {
  auto& intersection    = *getPRD<scene_intersection>();
  intersection.instance = 0;
  intersection.element  = 0;
  intersection.uv       = {0, 0};
  intersection.distance = 0;
  intersection.hit      = false;
}

// scene intersection via shaders
static scene_intersection intersect_scene(
    const cutrace_scene& scene, OptixTraversableHandle bvh, const ray3f& ray) {
  auto     intersection = scene_intersection{};
  uint32_t u0, u1;
  packPointer(&intersection, u0, u1);
  optixTrace(bvh, {ray.o.x, ray.o.y, ray.o.z}, {ray.d.x, ray.d.y, ray.d.z},
      ray.tmin, ray.tmax, 0.0f, OptixVisibilityMask(255),
      OPTIX_RAY_FLAG_DISABLE_ANYHIT, 0, 0, 0, u0, u1);
  return intersection;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// TRACE FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Convenience functions
[[maybe_unused]] static vec3f eval_position(
    const scene_data& scene, const scene_intersection& intersection) {
  return eval_position(scene, scene.instances[intersection.instance],
      intersection.element, intersection.uv);
}
[[maybe_unused]] static vec3f eval_normal(
    const scene_data& scene, const scene_intersection& intersection) {
  return eval_normal(scene, scene.instances[intersection.instance],
      intersection.element, intersection.uv);
}
[[maybe_unused]] static vec3f eval_element_normal(
    const scene_data& scene, const scene_intersection& intersection) {
  return eval_element_normal(
      scene, scene.instances[intersection.instance], intersection.element);
}
[[maybe_unused]] static vec3f eval_shading_position(const scene_data& scene,
    const scene_intersection& intersection, const vec3f& outgoing) {
  return eval_shading_position(scene, scene.instances[intersection.instance],
      intersection.element, intersection.uv, outgoing);
}
[[maybe_unused]] static vec3f eval_shading_normal(const scene_data& scene,
    const scene_intersection& intersection, const vec3f& outgoing) {
  return eval_shading_normal(scene, scene.instances[intersection.instance],
      intersection.element, intersection.uv, outgoing);
}
[[maybe_unused]] static vec2f eval_texcoord(
    const scene_data& scene, const scene_intersection& intersection) {
  return eval_texcoord(scene, scene.instances[intersection.instance],
      intersection.element, intersection.uv);
}
[[maybe_unused]] static material_point eval_material(
    const scene_data& scene, const scene_intersection& intersection) {
  return eval_material(scene, scene.instances[intersection.instance],
      intersection.element, intersection.uv);
}
[[maybe_unused]] static bool is_volumetric(
    const scene_data& scene, const scene_intersection& intersection) {
  return is_volumetric(scene, scene.instances[intersection.instance]);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// TRACE FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

static ray3f eval_camera(
    const cutrace_camera& camera, const vec2f& image_uv, const vec2f& lens_uv) {
  auto film = camera.aspect >= 1
                  ? vec2f{camera.film, camera.film / camera.aspect}
                  : vec2f{camera.film * camera.aspect, camera.film};
  auto q    = vec3f{
      film.x * (0.5f - image_uv.x), film.y * (image_uv.y - 0.5f), camera.lens};
  // ray direction through the lens center
  auto dc = -normalize(q);
  // point on the lens
  auto e = vec3f{
      lens_uv.x * camera.aperture / 2, lens_uv.y * camera.aperture / 2, 0};
  // point on the focus plane
  auto p = dc * camera.focus / abs(dc.z);
  // correct ray direction to account for camera focusing
  auto d = normalize(p - e);
  // done
  return ray3f{
      transform_point(camera.frame, e), transform_direction(camera.frame, d)};
}

static void trace_pixel(cutrace_state& state, const cutrace_scene& scene,
    const cutrace_bvh& bvh, const cutrace_lights& lights, int i, int j,
    const cutrace_params& params) {
  auto& camera = scene.cameras[params.camera];
  // auto  sampler = get_trace_sampler_func(params);
  auto idx = state.width * j + i;
  auto uv  = vec2f{(i + rand1f(state.rngs[idx])) / state.width,
      (j + rand1f(state.rngs[idx])) / state.height};
  auto ray = eval_camera(camera, uv, {0, 0});

  // trace ray
  auto intersection = intersect_scene(scene, bvh, ray);

  // and write to frame buffer ...
  if (intersection.hit) {
    auto& instance = scene.instances[intersection.instance];
    auto& shape    = scene.shapes[instance.shape];
    auto  material = eval_material(scene, intersection);

    auto color = material.color;

    state.image[i + j * state.width] += {color.x, color.y, color.z, 1};
  } else {
    state.image[i + j * state.width] += {0, 1, 0, 1};
  }
}

// raygen shader
optix_shader void __raygen__trace_pixel() {
  // pixel index
  auto ij  = optixGetLaunchIndex();
  auto idx = ij.y * globals.state.width + ij.x;

  // initialize state on first sample
  if (globals.state.samples == 0) {
    globals.state.image[idx] = {0, 0, 0, 0};
    globals.state.rngs[idx]  = make_rng(98273987, idx * 2 + 1);
  }

  // run shading
  auto nsamples = 16;
  for (auto sample = 0; sample < nsamples; sample++) {
    trace_pixel(globals.state, globals.scene, globals.bvh, cutrace_lights{},
        optixGetLaunchIndex().x, optixGetLaunchIndex().y, cutrace_params{});
  }

  // normalize output
  globals.state.image[idx] /= nsamples;
}

}  // namespace yocto