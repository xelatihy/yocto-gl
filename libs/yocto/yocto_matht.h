//
// # Yocto/MathT: Tiny library for math support in graphics applications.
//
// **This is an experimental implementation of Yocto/Math using templates.**
//
// Yocto/Math provides the basic math primitives used in grahics, including
// small-sized vectors and matrixes, frames, bounding boxes, transforms,
// color and geometry functions, random number generation, noise. 
//
//
// ## Small vectors, matrices and frames
//
// We provide common operations for small vectors and matrices typically used
// in graphics. In particular, we support 1-4 dimensional vectors `vec<T, N>`
// with aliases for float and int coordinates (`vec1f`, `vec2f`, `vec3f`, 
// `vec4f`, `vec1i`, `vec2i`, `vec3i`, `vec4i`).
//
// We support 1-4 dimensional matrices `mat<T, N>` with aliases for float 
// coordinates (`mat2f`, `mat3f`, `mat4f`) with
// matrix-matrix and matrix-vector products, transposes and inverses.
// Matrices are stored in column-major order and are accessed and
// constructed by column. The one dimensional version is for completeness only.
//
// To represent transformations, most of the library facilities prefer the use
// coordinate frames, aka rigid transforms, represented as `frame<T, N>` with 
// aliases for float coordinates (`frame2f`, `frame3f`). 
// The structure store three coordinate axes and the origin.
// This is equivalent to a rigid transform written as a column-major affine
// matrix. Transform operations are fater with this representation.
//
//
// ## Rays and bounding boxes
//
// We represent 2-3 dimensinal rays as `ray<T, N>` with float alises 
// (`ray2f`, `ray3f`). Each ray support initialization and evaluation.
//
// We represent 1-4 dimensional bounding boxes `bbox<T, N>` with float alises 
// (`bbox1f`, `bbox2f`, `bbox3f`, `bbox4f`).
// Each bounding box support construction from points and other bounding box.
// We provide operations to compute bounds for points, lines, triangles and
// quads.
//
//
// ## Transforms
//
// For both matrices and frames we support transform operations for points,
// vectors and directions (`transform_point()`, `transform_vector()`,
// `transform_direction()`). Transform matrices and frames can be
// constructed from basic translation, rotation and scaling, e.g. with
// `translation_mat()` or `translation_frame()` respectively, etc.
// For rotation we support axis-angle and quaternions, with slerp.
//
// TODO: better documentation
//
//
// ## Geometry functions
//
// The library supports basic geomtry functions such as computing
// line/triangle/quad normals and areas, picking points on triangles
// and the like. In these functions, triangles are parameterized with uv written
// w.r.t the (p1-p0) and (p2-p0) axis respectively. Quads are internally handled
// as pairs of two triangles (p0,p1,p3) and (p2,p3,p1), with the uv coordinates
// of the second triangle corrected as 1-u and 1-v to produce a quad
// parametrization where u and v go from 0 to 1. Degenerate quads with p2==p3
// represent triangles correctly, and this convention is used throught the
// library. This is equivalent to Intel's Embree.
//
// TODO: better documentation
//
//
// ## Color funtions
//
// This library support a small number of color operations helpful in writing
// graphics applications. In particular, we support color conversion to/from
// linear rgb, srgb, hsv, xyz, byte to flot color conversions and a few color
// manipulations like contrast and saturation.
//
// TODO: better documentation
//
//
// ## Random Number Generation
//
// This library supports generting random numbers using the PCG32 algorithm,
// that is a portable generator well suited for graphics applications.
//
// 1. initialize the random number generator with `make_rng()`
// 2. if necessary, you can reseed the rng with `seed_rng()`
// 3. generate random integers in an interval with `rand1i()`
// 4. generate random floats in the [0,1) range with `rand1f()`,
//    `rand2f()`, `rand3f()`, `rand4f()`
//
//
// ## Noise Functions
//
// We support generation of Perlin noise based on the stb libraries.
//
// 1. use `perlin_noise()` to generate Perlin noise with optional wrapping
// 2. use `perlin_ridge()`, `perlin_fbm()` and `perlin_turbulence()` for fractal
//    noises
//
//
// ## Shading functions
//
// We include a few functions to help writing shaders for path tracing.
//
// 1. use `fresnel_dielectric()` or `fresnel_conductor()` to evaluate the
//    fresnel term for dielectrics or conductors; use `fresnel_schlick()` for
//    the Schlick fresnel approximation
// 2. use `eta_to_reflectivity()` and `reflective_to_eta()` to convert eta to
//    reflectivity and vice-versa; use `eta_to_edgetint()` and
//    `edgetint_to_eta()`
// 3. use `microfacet_distribution()` and `microfacet_shadowing()` to evaluate
//    a microfacet distribution and its associated shadowing term
// 4. evaluate BRDF lobes with
//    - `eval_diffuse_reflection()`: diffuse brdf
//    - `eval_microfacet_reflection()`: specular brdf for dielectrics and metals
//    - `eval_microfacet_transmission()`: transmission brdf for thin dielectrics
//    - `eval_microfacet_refraction()`: refraction brdf for dielectrics
// 5. sample BRDF lobes with `sample_XXX()` using the above lobe names
// 6. compute the PDF for BRDF lobe sampling with `sample_XXX_pdf()` using the
//    above lobe names
//
//
// ## Monte Carlo helpers
//
// We include many method to generate random points and directions. These may be
// used in path tracing or procedural generation.
//
// 1. use `sample_XXX()` to warp random numbers in [0,1)^k domains to the
//   desired domain; in particular we support `sample_hemisphere()`,
//   `sample_sphere()`, `sample_hemisphere_cos()`,
//   `sample_hemisphere_cospower()`. `sample_disk()`. `sample_cylinder()`.
//   `sample_triangle()`, `sample_quad()`
// 2. use `sample_discrete()` to sample from a descreet distribution
// 3. use `sample_XXX_pdf()` to compute the PDF of the sampling functions
//
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2020 Fabio Pellacini
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
// -----------------------------------------------------------------------------
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
// -----------------------------------------------------------------------------
//
//

#ifndef _YOCTO_MATH_H_
#define _YOCTO_MATH_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <vector>

// -----------------------------------------------------------------------------
// MATH CONSTANTS AND FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto::math {

using byte   = unsigned char;
using uint   = unsigned int;
using ushort = unsigned short;

inline const double pi  = 3.14159265358979323846;
inline const float  pif = (float)pi;

inline const auto int_max = std::numeric_limits<int>::max();
inline const auto int_min = std::numeric_limits<int>::lowest();
inline const auto flt_max = std::numeric_limits<float>::max();
inline const auto flt_min = std::numeric_limits<float>::lowest();
inline const auto flt_eps = std::numeric_limits<float>::epsilon();

template <typename T>
inline const auto type_max = std::numeric_limits<T>::max();
template <typename T>
inline const auto type_min = std::numeric_limits<T>::lowest();
template <typename T>
inline const auto type_eps = std::numeric_limits<T>::epsilon();

template <typename T>
inline T abs(T a);
template <typename T>
inline T min(T a, T b);
template <typename T>
inline T max(T a, T b);
template <typename T>
inline T clamp(T a, T min, T max);
template <typename T>
inline T sign(T a);
template <typename T>
inline T sqrt(T a);
template <typename T>
inline T sin(T a);
template <typename T>
inline T cos(T a);
template <typename T>
inline T tan(T a);
template <typename T>
inline T asin(T a);
template <typename T>
inline T acos(T a);
template <typename T>
inline T atan(T a);
template <typename T>
inline T log(T a);
template <typename T>
inline T exp(T a);
template <typename T>
inline T log2(T a);
template <typename T>
inline T exp2(T a);
template <typename T>
inline T pow(T a, T b);
template <typename T>
inline T pow(T a, int b);
template <typename T>
inline T isfinite(T a);
template <typename T>
inline T atan2(T a, T b);
template <typename T>
inline T fmod(T a, T b);
template <typename T>
inline T radians(T a);
template <typename T>
inline T degrees(T a);
template <typename T>
inline T lerp(T a, T b, T u);
template <typename T>
inline void swap(T& a, T& b);
template <typename T>
inline T smoothstep(T a, T b, T u);
template <typename T>
inline T bias(T a, T bias);
template <typename T>
inline T gain(T a, T gain);
template <typename T>
inline int pow2(int a);

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// VECTORS
// -----------------------------------------------------------------------------
namespace yocto::math {

// Short vector with static size. Implementation uses specializations to support
// direct member access.
template <typename T, size_t N>
struct vec;

// One-dimensional vector.
template <typename T>
struct vec<T, 1> {
  T x = 0;

  vec();
  explicit vec(T x);
  explicit operator bool() const;

  T&       operator[](int i);
  const T& operator[](int i) const;
};

// Two-dimensional vector.
template <typename T>
struct vec<T, 2> {
  T x = 0;
  T y = 0;

  vec();
  vec(T x, T y);
  explicit vec(T v);
  explicit operator bool() const;

  T&       operator[](int i);
  const T& operator[](int i) const;
};

// Three-dimensional vector.
template <typename T>
struct vec<T, 3> {
  T x = 0;
  T y = 0;
  T z = 0;

  vec();
  vec(T x, T y, T z);
  vec(const vec<T, 2>& v, T z);
  explicit vec(T v);
  explicit operator bool() const;

  T&       operator[](int i);
  const T& operator[](int i) const;
};

// Four-dimensional vector.
template <typename T>
struct vec<T, 4> {
  T x = 0;
  T y = 0;
  T z = 0;
  T w = 0;

  vec();
  vec(T x, T y, T z, T w);
  vec(const vec<T, 3>& v, T w);
  explicit vec(T v);
  explicit operator bool() const;

  T&       operator[](int i);
  const T& operator[](int i) const;
};

// Type aliases
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

// Zero std::vector constants.
inline const auto zero1f = vec1f{0};
inline const auto zero2f = vec2f{0, 0};
inline const auto zero3f = vec3f{0, 0, 0};
inline const auto zero4f = vec4f{0, 0, 0, 0};
inline const auto zero1i = vec1i{0};
inline const auto zero2i = vec2i{0, 0};
inline const auto zero3i = vec3i{0, 0, 0};
inline const auto zero4i = vec4i{0, 0, 0, 0};
inline const auto zero3b = vec3b{0, 0, 0};
inline const auto zero4b = vec4b{0, 0, 0, 0};

// Element access
template <typename T>
inline vec<T, 3>& xyz(vec<T, 4>& a);
template <typename T>
inline const vec<T, 3>& xyz(const vec<T, 4>& a);

// Vector comparison operations.
template <typename T, size_t N>
inline bool operator==(const vec<T, N>& a, const vec<T, N>& b);
template <typename T, size_t N>
inline bool operator!=(const vec<T, N>& a, const vec<T, N>& b);

// Vector operations.
template <typename T, size_t N>
inline vec<T, N> operator+(const vec<T, N>& a);
template <typename T, size_t N>
inline vec<T, N> operator-(const vec2f& a);
template <typename T, size_t N>
inline vec<T, N> operator+(const vec<T, N>& a, const vec<T, N>& b);
template <typename T, size_t N, typename T1>
inline vec<T, N> operator+(const vec<T, N>& a, T1 b);
template <typename T, size_t N, typename T1>
inline vec<T, N> operator+(T1 a, const vec<T, N>& b);
template <typename T, size_t N>
inline vec<T, N> operator-(const vec<T, N>& a, const vec<T, N>& b);
template <typename T, size_t N, typename T1>
inline vec<T, N> operator-(const vec<T, N>& a, T1 b);
template <typename T, size_t N, typename T1>
inline vec<T, N> operator-(T1 a, const vec<T, N>& b);
template <typename T, size_t N>
inline vec<T, N> operator*(const vec<T, N>& a, const vec<T, N>& b);
template <typename T, size_t N, typename T1>
inline vec<T, N> operator*(const vec<T, N>& a, T1 b);
template <typename T, size_t N, typename T1>
inline vec<T, N> operator*(T1 a, const vec2f& b);
template <typename T, size_t N>
inline vec<T, N> operator/(const vec<T, N>& a, const vec<T, N>& b);
template <typename T, size_t N, typename T1>
inline vec<T, N> operator/(const vec<T, N>& a, T1 b);
template <typename T, size_t N, typename T1>
inline vec<T, N> operator/(T1 a, const vec<T, N>& b);

// Vector assignments
template <typename T, size_t N>
inline vec<T, N>& operator+=(vec<T, N>& a, const vec<T, N>& b);
template <typename T, size_t N, typename T1>
inline vec<T, N>& operator+=(vec<T, N>& a, T1 b);
template <typename T, size_t N>
inline vec<T, N>& operator-=(vec<T, N>& a, const vec<T, N>& b);
template <typename T, size_t N, typename T1>
inline vec<T, N>& operator-=(vec<T, N>& a, T1 b);
template <typename T, size_t N>
inline vec<T, N>& operator*=(vec<T, N>& a, const vec<T, N>& b);
template <typename T, size_t N, typename T1>
inline vec<T, N>& operator*=(vec<T, N>& a, T1 b);
template <typename T, size_t N>
inline vec<T, N>& operator/=(vec<T, N>& a, const vec<T, N>& b);
template <typename T, size_t N, typename T1>
inline vec<T, N>& operator/=(vec<T, N>& a, T1 b);

// Vector products and lengths.
template <typename T, size_t N>
inline T dot(const vec<T, N>& a, const vec<T, N>& b);
template <typename T, size_t N>
inline T length(const vec<T, N>& a);
template <typename T, size_t N>
inline vec<T, N> normalize(const vec<T, N>& a);
template <typename T, size_t N>
inline T distance(const vec<T, N>& a, const vec<T, N>& b);
template <typename T, size_t N>
inline T distance_squared(const vec<T, N>& a, const vec<T, N>& b);
template <typename T>
inline T cross(const vec<T, 2>& a, const vec<T, 2>& b);
template <typename T>
inline vec<T, 3> cross(const vec<T, 3>& a, const vec<T, 3>& b);
template <typename T>
inline T angle(const vec<T, 3>& a, const vec<T, 3>& b);

// Orthogonal vectors.
template <typename T>
inline vec<T, 3> orthogonal(const vec<T, 3>& v);
template <typename T>
inline vec<T, 3> orthonormalize(const vec<T, 3>& a, const vec<T, 3>& b);

// Reflected and refracted std::vector.
template <typename T>
inline vec<T, 3> reflect(const vec<T, 3>& w, const vec<T, 3>& n);
template <typename T>
inline vec<T, 3> refract(const vec<T, 3>& w, const vec<T, 3>& n, T inv_eta);

// Slerp
template <typename T, size_t N, typename T1>
inline vec<T, 4> slerp(const vec<T, 4>& a, const vec<T, 4>& b, T1 u);

// Max element and clamp.
template <typename T, size_t N, typename T1>
inline vec<T, N> max(const vec<T, N>& a, T1 b);
template <typename T, size_t N, typename T1>
inline vec<T, N> min(const vec<T, N>& a, T1 b);
template <typename T, size_t N>
inline vec<T, N> max(const vec<T, N>& a, const vec<T, N>& b);
template <typename T, size_t N>
inline vec<T, N> min(const vec<T, N>& a, const vec<T, N>& b);
template <typename T, size_t N, typename T1, typename T2>
inline vec<T, N> clamp(const vec<T, N>& x, T1 min, T2 max);
template <typename T, size_t N, typename T1>
inline vec<T, N> lerp(const vec<T, N>& a, const vec<T, N>& b, T1 u);
template <typename T, size_t N>
inline vec<T, N> lerp(
    const vec<T, N>& a, const vec<T, N>& b, const vec<T, N>& u);

template <typename T, size_t N>
inline T max(const vec<T, N>& a);
template <typename T, size_t N>
inline T min(const vec<T, N>& a);
template <typename T, size_t N>
inline T sum(const vec<T, N>& a);
template <typename T, size_t N>
inline T mean(const vec<T, N>& a);

// Functions applied to std::vector elements
template <typename T, size_t N>
inline vec<T, N> abs(const vec<T, N>& a);
template <typename T, size_t N>
inline vec<T, N> sqrt(const vec<T, N>& a);
template <typename T, size_t N>
inline vec<T, N> exp(const vec<T, N>& a);
template <typename T, size_t N>
inline vec<T, N> log(const vec<T, N>& a);
template <typename T, size_t N>
inline vec<T, N> exp2(const vec<T, N>& a);
template <typename T, size_t N>
inline vec<T, N> log2(const vec<T, N>& a);
template <typename T, size_t N, typename T1>
inline vec<T, N> pow(const vec<T, N>& a, T1 b);
template <typename T, size_t N>
inline vec<T, N> pow(const vec<T, N>& a, const vec<T, N>& b);
template <typename T, size_t N, typename T1>
inline vec<T, N> gain(const vec<T, N>& a, T1 b);
template <typename T, size_t N>
inline bool isfinite(const vec<T, N>& a);
template <typename T, size_t N>
inline void swap(vec<T, N>& a, vec<T, N>& b);

// Quaternion operatons represented as xi + yj + zk + w
// const auto identity_quat4f = vec4f{0, 0, 0, 1};
template <typename T, typename T1>
inline vec<T, 4> quat_mul(const vec<T, 4>& a, T1 b);
template <typename T>
inline vec<T, 4> quat_mul(const vec<T, 4>& a, const vec<T, 4>& b);
template <typename T>
inline vec<T, 4> quat_conjugate(const vec<T, 4>& a);
template <typename T>
inline vec<T, 4> quat_inverse(const vec<T, 4>& a);

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// VECTOR HASHING
// -----------------------------------------------------------------------------
namespace std {

// Hash functor for std::vector for use with hash_map
template <typename T, size_t N>
struct hash<yocto::math::vec<T, N>>;

}  // namespace std

// -----------------------------------------------------------------------------
// MATRICES
// -----------------------------------------------------------------------------
namespace yocto::math {

// Small Fixed-size matrices stored in column major format.
template <typename T, size_t M, size_t N>
struct mat;

// One-dimensional matrix stored in column major format.
template <typename T, size_t M>
struct mat<T, M, 1> {
  vec<T, M> x = vec<T, M>{0};

  mat();
  mat(const vec<T, M>& x);

  vec<T, M>&       operator[](int i);
  const vec<T, M>& operator[](int i) const;
};

// Two-dimensional matrix stored in column major format.
template <typename T, size_t M>
struct mat<T, M, 2> {
  vec<T, M> x = vec<T, M>{0};
  vec<T, M> y = vec<T, M>{0};

  mat();
  mat(const vec<T, M>& x, const vec<T, M>& y);

  vec<T, M>&       operator[](int i);
  const vec<T, M>& operator[](int i) const;
};

// Three-dimensional matrix stored in column major format.
template <typename T, size_t M>
struct mat<T, M, 3> {
  vec<T, M> x = vec<T, M>{0};
  vec<T, M> y = vec<T, M>{0};
  vec<T, M> z = vec<T, M>{0};

  mat();
  mat(const vec<T, M>& x, const vec<T, M>& y, const vec<T, M>& z);

  vec<T, M>&       operator[](int i);
  const vec<T, M>& operator[](int i) const;
};

// Four-dimensional matrix stored in column major format.
template <typename T, size_t M>
struct mat<T, M, 4> {
  vec<T, M> x = vec<T, M>{0};
  vec<T, M> y = vec<T, M>{0};
  vec<T, M> z = vec<T, M>{0};
  vec<T, M> w = vec<T, M>{0};

  mat();
  mat(const vec<T, M>& x, const vec<T, M>& y, const vec<T, M>& z,
      const vec<T, M>& w);

  vec<T, M>&       operator[](int i);
  const vec<T, M>& operator[](int i) const;
};

// Type aliases
using mat1f = mat<float, 1, 1>;
using mat2f = mat<float, 2, 2>;
using mat3f = mat<float, 3, 3>;
using mat4f = mat<float, 4, 4>;

// Identity matrices constants.
inline const auto identity2x2f = mat2f{{1, 0}, {0, 1}};
inline const auto identity3x3f = mat3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
inline const auto identity4x4f = mat4f{
    {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};

// Matrix comparisons.
template <typename T, size_t M, size_t N>
inline bool operator==(const mat<T, M, N>& a, const mat<T, M, N>& b);
template <typename T, size_t M, size_t N>
inline bool operator!=(const mat<T, M, N>& a, const mat<T, M, N>& b);

// Matrix operations.
template <typename T, size_t M, size_t N>
inline mat<T, M, N> operator+(const mat<T, M, N>& a, const mat<T, M, N>& b);
template <typename T, size_t M, size_t N, typename T1>
inline mat<T, M, N> operator*(const mat<T, M, N>& a, T1 b);
template <typename T, size_t M, size_t N>
inline vec<T, M> operator*(const mat<T, M, N>& a, const vec<T, N>& b);
template <typename T, size_t M, size_t N>
inline vec<T, N> operator*(const vec<T, M>& a, const mat<T, M, N>& b);
template <typename T, size_t M, size_t N, size_t K>
inline mat<T, M, K> operator*(const mat<T, M, N>& a, const mat<T, N, K>& b);

// Matrix assignments.
template <typename T, size_t M, size_t N>
inline mat<T, M, N>& operator+=(mat<T, M, N>& a, const mat<T, M, N>& b);
template <typename T, size_t M, size_t N>
inline mat2f& operator*=(mat<T, M, N>& a, const mat<T, M, N>& b);
template <typename T, size_t M, size_t N, typename T1>
inline mat<T, M, N>& operator*=(mat<T, M, N>& a, T1 b);

// Matrix diagonals and transposes.
template <typename T, size_t N>
inline vec<T, N> diagonal(const mat<T, N, N>& a);
template <typename T, size_t N>
inline mat<T, N, N> transpose(const mat<T, N, N>& a);

// Matrix adjoints, determinants and inverses.
template <typename T, size_t N>
inline T determinant(const mat<T, N, N>& a);
template <typename T, size_t N>
inline mat<T, N, N> adjoint(const mat<T, N, N>& a);
template <typename T, size_t N>
inline mat<T, N, N> inverse(const mat<T, N, N>& a);

// Constructs a basis from a direction
template <typename T>
inline mat<T, 3, 3> basis_fromz(const vec<T, 3>& v);

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// RIGID BODY TRANSFORMS/FRAMES
// -----------------------------------------------------------------------------
namespace yocto::math {

// Fixed-size rigid frames stored as a column-major affine transform matrix.
template <typename T, size_t N>
struct frame;

// One-dimensional frame stored in column major format.
template <typename T>
struct frame<T, 1> {
  vec<T, 1> x = {1};
  vec<T, 1> o = {0};

  frame();
  frame(const vec<T, 1>& x, const vec<T, 1>& o);
  explicit frame(const vec<T, 1>& o);
  frame(const mat<T, 1, 1>& m, const vec<T, 1>& t);
  explicit frame(const mat<T, 2, 2>& m);
  operator mat<T, 2, 2>() const;

  vec<T, 1>&       operator[](int i);
  const vec<T, 1>& operator[](int i) const;
};

// Two-dimensional frame stored in column major format.
template <typename T>
struct frame<T, 2> {
  vec<T, 2> x = {1, 0};
  vec<T, 2> y = {0, 1};
  vec<T, 2> o = {0, 0};

  frame();
  frame(const vec<T, 2>& x, const vec<T, 2>& y, const vec<T, 2>& o);
  explicit frame(const vec<T, 2>& o);
  frame(const mat<T, 2, 2>& m, const vec<T, 2>& t);
  explicit frame(const mat<T, 3, 3>& m);
  operator mat<T, 3, 3>() const;

  vec<T, 2>&       operator[](int i);
  const vec<T, 2>& operator[](int i) const;
};

// Rigid frames stored as a column-major affine transform matrix.
template <typename T>
struct frame<T, 3> {
  vec<T, 3> x = {1, 0, 0};
  vec<T, 3> y = {0, 1, 0};
  vec<T, 3> z = {0, 0, 1};
  vec<T, 3> o = {0, 0, 0};

  frame();
  frame(const vec<T, 3>& x, const vec<T, 3>& y, const vec<T, 3>& z,
      const vec<T, 3>& o);
  explicit frame(const vec<T, 3>& o);
  frame(const mat<T, 3, 3>& m, const vec<T, 3>& t);
  explicit frame(const mat<T, 4, 4>& m);
  operator mat<T, 4, 4>() const;

  vec<T, 3>&       operator[](int i);
  const vec<T, 3>& operator[](int i) const;
};

// Type aliases
using frame2f = frame<float, 2>;
using frame3f = frame<float, 3>;

// Indentity frames.
inline const auto identity2x3f = frame2f{{1, 0}, {0, 1}, {0, 0}};
inline const auto identity3x4f = frame3f{
    {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};

// Frame properties
template <typename T, size_t N>
inline const mat<T, N, N>& rotation(const frame<T, N>& a);

// Frame comparisons.
template <typename T, size_t N>
inline bool operator==(const frame<T, N>& a, const frame<T, N>& b);
template <typename T, size_t N>
inline bool operator!=(const frame<T, N>& a, const frame<T, N>& b);

// Frame composition, equivalent to affine matrix product.
template <typename T, size_t N>
inline frame<T, N> operator*(const frame<T, N>& a, const frame<T, N>& b);
template <typename T, size_t N>
inline frame<T, N>& operator*=(frame<T, N>& a, const frame<T, N>& b);

// Frame inverse, equivalent to rigid affine inverse.
template <typename T, size_t N>
inline frame<T, N> inverse(const frame<T, N>& a, bool non_rigid = false);

// Frame construction from axis.
template <typename T>
inline frame<T, 3> frame_fromz(const vec<T, 3>& o, const vec<T, 3>& v);
template <typename T>
inline frame<T, 3> frame_fromzx(
    const vec<T, 3>& o, const vec<T, 3>& z_, const vec<T, 3>& x_);

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// QUATERNIONS
// -----------------------------------------------------------------------------
namespace yocto::math {

// Quaternions to represent rotations
template <typename T, size_t N>
struct quat;

// Quaternions to represent rotations
template <typename T>
struct quat<T, 4> {
  T x = 0;
  T y = 0;
  T z = 0;
  T w = 0;

  // constructors
  quat();
  quat(T x, T y, T z, T w);
};

// Type alias
using quat4f = quat<float, 4>;

// Constants
inline const auto identity_quat4f = quat4f{0, 0, 0, 1};

// Quaternion operatons
template <typename T>
inline quat<T, 4> operator+(const quat<T, 4>& a, const quat<T, 4>& b);
template <typename T, typename T1>
inline quat<T, 4> operator*(const quat<T, 4>& a, T1 b);
template <typename T, typename T1>
inline quat<T, 4> operator/(const quat<T, 4>& a, T1 b);
template <typename T>
inline quat<T, 4> operator*(const quat<T, 4>& a, const quat<T, 4>& b);

// Quaterion operations
template <typename T>
inline T dot(const quat<T, 4>& a, const quat<T, 4>& b);
template <typename T>
inline T length(const quat<T, 4>& a);
template <typename T>
inline quat<T, 4> normalize(const quat<T, 4>& a);
template <typename T>
inline quat<T, 4> conjugate(const quat<T, 4>& a);
template <typename T>
inline quat<T, 4> inverse(const quat<T, 4>& a);
template <typename T>
inline T uangle(const quat<T, 4>& a, const quat<T, 4>& b);
template <typename T, typename T1>
inline quat<T, 4> lerp(const quat<T, 4>& a, const quat<T, 4>& b, T1 t);
template <typename T, typename T1>
inline quat<T, 4> nlerp(const quat<T, 4>& a, const quat<T, 4>& b, T1 t);
template <typename T, typename T1>
inline quat<T, 4> slerp(const quat<T, 4>& a, const quat<T, 4>& b, T1 t);

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// AXIS ALIGNED BOUNDING BOXES
// -----------------------------------------------------------------------------
namespace yocto::math {

// Axis aligned bounding box represented as a min/max vector pairs.
template <typename T, size_t N>
struct bbox {
  vec<T, N> min = vec<T, N>{type_max<T>};
  vec<T, N> max = vec<T, N>{type_min<T>};

  bbox();
  bbox(const vec<T, N>& min, const vec<T, N>& max);

  vec<T, N>&       operator[](int i);
  const vec<T, N>& operator[](int i) const;
};

// Type alias
using bbox1f = bbox<float, 1>;
using bbox2f = bbox<float, 2>;
using bbox3f = bbox<float, 3>;
using bbox4f = bbox<float, 4>;

// Empty bbox constant.
inline const auto invalidb1f = bbox1f{};
inline const auto invalidb2f = bbox2f{};
inline const auto invalidb3f = bbox3f{};
inline const auto invalidb4f = bbox4f{};

// Bounding box properties
template <typename T, size_t N>
inline vec<T, N> center(const bbox<T, N>& a);
template <typename T, size_t N>
inline vec<T, N> size(const bbox<T, N>& a);

// Bounding box comparisons.
template <typename T, size_t N>
inline bool operator==(const bbox<T, N>& a, const bbox<T, N>& b);
template <typename T, size_t N>
inline bool operator!=(const bbox<T, N>& a, const bbox<T, N>& b);

// Bounding box expansions with points and other boxes.
template <typename T, size_t N>
inline bbox<T, N> merge(const bbox<T, N>& a, const vec<T, N>& b);
template <typename T, size_t N>
inline bbox<T, N> merge(const bbox<T, N>& a, const bbox<T, N>& b);
template <typename T, size_t N>
inline void expand(bbox<T, N>& a, const vec<T, N>& b);
template <typename T, size_t N>
inline void expand(bbox<T, N>& a, const bbox<T, N>& b);

// Primitive bounds.
template <typename T, size_t N>
inline bbox<T, N> point_bounds(const vec<T, N>& p);
template <typename T, size_t N, typename T1>
inline bbox<T, N> point_bounds(const vec<T, N>& p, T1 r);
template <typename T, size_t N>
inline bbox<T, N> line_bounds(const vec<T, N>& p0, const vec<T, N>& p1);
template <typename T, size_t N, typename T1, typename T2>
inline bbox<T, N> line_bounds(
    const vec<T, N>& p0, const vec<T, N>& p1, T1 r0, T2 r1);
template <typename T, size_t N>
inline bbox<T, N> triangle_bounds(
    const vec<T, N>& p0, const vec<T, N>& p1, const vec<T, N>& p2);
template <typename T, size_t N>
inline bbox<T, N> quad_bounds(const vec<T, N>& p0, const vec<T, N>& p1,
    const vec<T, N>& p2, const vec<T, N>& p3);

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// RAYS
// -----------------------------------------------------------------------------
namespace yocto::math {

// Ray esplison
inline const auto ray_eps  = 1e-4;
inline const auto ray_epsd = 1e-4;
inline const auto ray_epsf = 1e-4f;

// Ray with origin, direction and min/max t value.
template <typename T, size_t N>
struct ray;

// Two-dimensional ray with origin, direction and min/max t value.
template <typename T>
struct ray<T, 2> {
  vec<T, 2> o    = {0, 0};
  vec<T, 2> d    = {0, 1};
  T         tmin = ray_eps;
  T         tmax = flt_max;

  ray();
  ray(const vec<T, 2>& o, const vec<T, 2>& d, T tmin = ray_eps,
      T tmax = flt_max);
};

// Two-dimensional ray with origin, direction and min/max t value.
template <typename T>
struct ray<T, 3> {
  vec<T, 3> o    = {0, 0, 0};
  vec<T, 3> d    = {0, 0, 1};
  T         tmin = ray_eps;
  T         tmax = flt_max;

  ray();
  ray(const vec<T, 3>& o, const vec<T, 3>& d, T tmin = ray_eps,
      T tmax = flt_max);
};

// Type alias
using ray2f = ray<float, 2>;
using ray3f = ray<float, 3>;

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// TRANSFORMS
// -----------------------------------------------------------------------------
namespace yocto::math {

// Transforms points, vectors and directions by matrices.
template <typename T>
inline vec<T, 2> transform_point(const mat<T, 3, 3>& a, const vec<T, 2>& b);
template <typename T>
inline vec<T, 2> transform_vector(const mat<T, 3, 3>& a, const vec<T, 2>& b);
template <typename T>
inline vec<T, 2> transform_direction(const mat<T, 3, 3>& a, const vec<T, 2>& b);
template <typename T>
inline vec<T, 2> transform_normal(const mat<T, 3, 3>& a, const vec<T, 2>& b);
template <typename T>
inline vec<T, 2> transform_vector(const mat<T, 2, 2>& a, const vec<T, 2>& b);
template <typename T>
inline vec<T, 2> transform_direction(const mat<T, 2, 2>& a, const vec<T, 2>& b);
template <typename T>
inline vec<T, 2> transform_normal(const mat<T, 2, 2>& a, const vec<T, 2>& b);

template <typename T>
inline vec<T, 3> transform_point(const mat<T, 4, 4>& a, const vec<T, 3>& b);
template <typename T>
inline vec<T, 3> transform_vector(const mat<T, 4, 4>& a, const vec<T, 3>& b);
template <typename T>
inline vec<T, 3> transform_direction(const mat<T, 4, 4>& a, const vec<T, 3>& b);
template <typename T>
inline vec<T, 3> transform_vector(const mat<T, 3, 3>& a, const vec<T, 3>& b);
template <typename T>
inline vec<T, 3> transform_direction(const mat<T, 3, 3>& a, const vec<T, 3>& b);
template <typename T>
inline vec<T, 3> transform_normal(const mat<T, 3, 3>& a, const vec<T, 3>& b);

// Transforms points, vectors and directions by frames.
template <typename T>
inline vec<T, 2> transform_point(const frame<T, 2>& a, const vec<T, 2>& b);
template <typename T>
inline vec<T, 2> transform_vector(const frame<T, 2>& a, const vec<T, 2>& b);
template <typename T>
inline vec<T, 2> transform_direction(const frame<T, 2>& a, const vec<T, 2>& b);
template <typename T>
inline vec<T, 2> transform_normal(
    const frame<T, 2>& a, const vec<T, 2>& b, bool non_rigid = false);

// Transforms points, vectors and directions by frames.
template <typename T>
inline vec<T, 3> transform_point(const frame<T, 3>& a, const vec<T, 3>& b);
template <typename T>
inline vec<T, 3> transform_vector(const frame<T, 3>& a, const vec<T, 3>& b);
template <typename T>
inline vec<T, 3> transform_direction(const frame<T, 3>& a, const vec<T, 3>& b);
template <typename T>
inline vec<T, 3> transform_normal(
    const frame<T, 3>& a, const vec<T, 3>& b, bool non_rigid = false);

// Transforms rays and bounding boxes by matrices.
template <typename T>
inline ray<T, 3> transform_ray(const mat<T, 4, 4>& a, const ray<T, 3>& b);
template <typename T>
inline ray<T, 3> transform_ray(const frame<T, 3>& a, const ray<T, 3>& b);
template <typename T>
inline bbox<T, 3> transform_bbox(const mat<T, 4, 4>& a, const bbox<T, 3>& b);
template <typename T>
inline bbox<T, 3> transform_bbox(const frame<T, 3>& a, const bbox<T, 3>& b);

// Translation, scaling and rotations transforms.
template <typename T>
inline frame<T, 3> translation_frame(const vec<T, 3>& a);
template <typename T>
inline frame<T, 3> scaling_frame(const vec<T, 3>& a);
template <typename T>
inline frame<T, 3> rotation_frame(const vec<T, 3>& axis, T angle);
template <typename T>
inline frame<T, 3> rotation_frame(const vec<T, 4>& quat);
template <typename T>
inline frame<T, 3> rotation_frame(const quat4f& quat);
template <typename T>
inline frame<T, 3> rotation_frame(const mat<T, 3, 3>& rot);

// Lookat frame. Z-axis can be inverted with inv_xz.
template <typename T>
inline frame<T, 3> lookat_frame(const vec<T, 3>& eye, const vec<T, 3>& center,
    const vec<T, 3>& up, bool inv_xz = false);

// OpenGL frustum, ortho and perspecgive matrices.
template <typename T>
inline mat<T, 4, 4> frustum_mat(T l, T r, T b, T t, T n, T f);
template <typename T>
inline mat<T, 4, 4> ortho_mat(T l, T r, T b, T t, T n, T f);
template <typename T>
inline mat<T, 4, 4> ortho2d_mat(T left, T right, T bottom, T top);
template <typename T>
inline mat<T, 4, 4> ortho_mat(T xmag, T ymag, T near, T far);
template <typename T>
inline mat<T, 4, 4> perspective_mat(T fovy, T aspect, T near, T far);
template <typename T>
inline mat<T, 4, 4> perspective_mat(T fovy, T aspect, T near);

// Rotation conversions.
template <typename T>
inline std::pair<vec<T, 3>, T> rotation_axisangle(const vec<T, 4>& quat);
template <typename T, typename T1>
inline vec<T, 4> rotation_quat(const vec<T, 3>& axis, T1 angle);
template <typename T>
inline vec<T, 4> rotation_quat(const vec<T, 4>& axisangle);

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// GEOMETRY UTILITIES
// -----------------------------------------------------------------------------
namespace yocto::math {

// Line properties.
template <typename T>
inline vec<T, 3> line_tangent(const vec<T, 3>& p0, const vec<T, 3>& p1);
template <typename T>
inline T line_length(const vec<T, 3>& p0, const vec<T, 3>& p1);

// Triangle properties.
template <typename T>
inline vec<T, 3> triangle_normal(
    const vec<T, 3>& p0, const vec<T, 3>& p1, const vec<T, 3>& p2);
template <typename T>
inline T triangle_area(
    const vec<T, 3>& p0, const vec<T, 3>& p1, const vec<T, 3>& p2);

// Quad propeties.
template <typename T>
inline vec<T, 3> quad_normal(const vec<T, 3>& p0, const vec<T, 3>& p1,
    const vec<T, 3>& p2, const vec<T, 3>& p3);
template <typename T>
inline T quad_area(const vec<T, 3>& p0, const vec<T, 3>& p1,
    const vec<T, 3>& p2, const vec<T, 3>& p3);

// Triangle tangent and bitangent from uv
template <typename T>
inline std::pair<vec<T, 3>, vec<T, 3>> triangle_tangents_fromuv(
    const vec<T, 3>& p0, const vec<T, 3>& p1, const vec<T, 3>& p2,
    const vec<T, 2>& uv0, const vec<T, 2>& uv1, const vec<T, 2>& uv2);

// Quad tangent and bitangent from uv. Note that we pass a current_uv since
// internally we may want to split the quad in two and we need to known where
// to do it. If not interested in the split, just pass zero2f here.
template <typename T>
inline std::pair<vec<T, 3>, vec<T, 3>> quad_tangents_fromuv(const vec<T, 3>& p0,
    const vec<T, 3>& p1, const vec<T, 3>& p2, const vec<T, 3>& p3,
    const vec<T, 2>& uv0, const vec<T, 2>& uv1, const vec<T, 2>& uv2,
    const vec<T, 2>& uv3, const vec<T, 2>& current_uv);

// Interpolates values over a line parameterized from a to b by u. Same as lerp.
template <typename T, typename T1>
inline T interpolate_line(const T& p0, const T& p1, T1 u);

// Interpolates values over a triangle parameterized by u and v along the
// (p1-p0) and (p2-p0) directions. Same as barycentric interpolation.
template <typename T, typename T1>
inline T interpolate_triangle(
    const T& p0, const T& p1, const T& p2, const vec<T1, 2>& uv);

// Interpolates values over a quad parameterized by u and v along the
// (p1-p0) and (p2-p1) directions. Same as bilinear interpolation.
template <typename T, typename T1>
inline T interpolate_quad(
    const T& p0, const T& p1, const T& p2, const T& p3, const vec<T1, 2>& uv);

// Interpolates values along a cubic Bezier segment parametrized by u.
template <typename T, typename T1>
inline T interpolate_bezier(
    const T& p0, const T& p1, const T& p2, const T& p3, T1 u);
// Computes the derivative of a cubic Bezier segment parametrized by u.

template <typename T, typename T1>
inline T interpolate_bezier_derivative(
    const T& p0, const T& p1, const T& p2, const T& p3, T1 u);

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// RAY-PRIMITIVE INTERSECTION FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto::math {

// Intersect a ray with a point (approximate)
template <typename T>
inline bool intersect_point(
    const ray<T, 3>& ray, const vec<T, 3>& p, T r, vec<T, 2>& uv, T& dist);

// Intersect a ray with a line
template <typename T>
inline bool intersect_line(const ray<T, 3>& ray, const vec<T, 3>& p0,
    const vec<T, 3>& p1, T r0, T r1, vec<T, 2>& uv, T& dist);

// Intersect a ray with a triangle
template <typename T>
inline bool intersect_triangle(const ray<T, 3>& ray, const vec<T, 3>& p0,
    const vec<T, 3>& p1, const vec<T, 3>& p2, vec<T, 2>& uv, T& dist);

// Intersect a ray with a quad.
template <typename T>
inline bool intersect_quad(const ray<T, 3>& ray, const vec<T, 3>& p0,
    const vec<T, 3>& p1, const vec<T, 3>& p2, const vec<T, 3>& p3,
    vec<T, 2>& uv, T& dist);

// Intersect a ray with a axis-aligned bounding box
template <typename T>
inline bool intersect_bbox(const ray<T, 3>& ray, const bbox<T, 3>& bbox);

// Intersect a ray with a axis-aligned bounding box
template <typename T>
inline bool intersect_bbox(
    const ray<T, 3>& ray, const vec<T, 3>& ray_dinv, const bbox<T, 3>& bbox);

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// POINT-PRIMITIVE DISTANCE FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto::math {

// Check if a point overlaps a position pos withint a maximum distance dist_max.
template <typename T>
inline bool overlap_point(const vec<T, 3>& pos, T dist_max, const vec<T, 3>& p,
    T r, vec<T, 2>& uv, T& dist);

// Compute the closest line uv to a give position pos.
template <typename T>
inline T closestuv_line(
    const vec<T, 3>& pos, const vec<T, 3>& p0, const vec<T, 3>& p1);

// Check if a line overlaps a position pos withint a maximum distance dist_max.
template <typename T>
inline bool overlap_line(const vec<T, 3>& pos, T dist_max, const vec<T, 3>& p0,
    const vec<T, 3>& p1, T r0, T r1, vec<T, 2>& uv, T& dist);

// Compute the closest triangle uv to a give position pos.
template <typename T>
inline vec<T, 2> closestuv_triangle(const vec<T, 3>& pos, const vec<T, 3>& p0,
    const vec<T, 3>& p1, const vec<T, 3>& p2);

// Check if a triangle overlaps a position pos withint a maximum distance
// dist_max.
template <typename T>
inline bool overlap_triangle(const vec<T, 3>& pos, T dist_max,
    const vec<T, 3>& p0, const vec<T, 3>& p1, const vec<T, 3>& p2, T r0, T r1,
    T r2, vec<T, 2>& uv, T& dist);

// Check if a quad overlaps a position pos withint a maximum distance dist_max.
template <typename T>
inline bool overlap_quad(const vec<T, 3>& pos, T dist_max, const vec<T, 3>& p0,
    const vec<T, 3>& p1, const vec<T, 3>& p2, const vec<T, 3>& p3, T r0, T r1,
    T r2, T r3, vec<T, 2>& uv, T& dist);

// Check if a bbox overlaps a position pos withint a maximum distance dist_max.
template <typename T>
inline bool distance_check_bbox(
    const vec<T, 3>& pos, T dist_max, const bbox<T, 3>& bbox);

// Check if two bboxe overlap.
template <typename T>
inline bool overlap_bbox(const bbox<T, 3>& bbox1, const bbox<T, 3>& bbox2);

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// COLOR OPERATIONS
// -----------------------------------------------------------------------------
namespace yocto::math {

// Conversion between flots and bytes
inline vec3b  float_to_byte(const vec3f& a);
inline vec3f  byte_to_float(const vec3b& a);
inline vec4b  float_to_byte(const vec4f& a);
inline vec4f  byte_to_float(const vec4b& a);
inline byte   float_to_byte(float a);
inline float  byte_to_float(byte a);
inline ushort float_to_ushort(float a);
inline float  ushort_to_float(ushort a);

// Conversion between reals in [0,1] and normalized ints [0, max_int]
template <typename I, typename T>
inline I real_to_nint(T a);
template <typename T, typename I>
inline T nint_to_real(I a);
template <typename I, typename T, size_t N>
inline vec<I, N> real_to_nint(const vec<T, N>& a);
template <typename T, typename I, size_t N>
inline vec<T, N> nint_to_real(const vec<I, N>& a);

// Luminance
template <typename T>
inline T luminance(const vec<T, 3>& a);

// sRGB non-linear curve
template <typename T>
inline T srgb_to_rgb(T srgb);
template <typename T>
inline T rgb_to_srgb(T rgb);
template <typename T>
inline vec<T, 3> srgb_to_rgb(const vec<T, 3>& srgb);
template <typename T>
inline vec<T, 3> rgb_to_srgb(const vec<T, 3>& rgb);
template <typename T>
inline vec<T, 4> srgb_to_rgb(const vec<T, 4>& srgb);
template <typename T>
inline vec<T, 4> rgb_to_srgb(const vec<T, 4>& rgb);

// Apply contrast. Grey should be 0.18 for linear and 0.5 for gamma.
template <typename T>
inline vec<T, 3> lincontrast(const vec<T, 3>& rgb, T contrast, T grey);
// Apply contrast in log2. Grey should be 0.18 for linear and 0.5 for gamma.
template <typename T>
inline vec<T, 3> logcontrast(const vec<T, 3>& rgb, T logcontrast, T grey);
// Apply an s-shaped contrast.
template <typename T>
inline vec<T, 3> contrast(const vec<T, 3>& rgb, T contrast);
// Apply saturation.
template <typename T>
inline vec<T, 3> saturate(const vec<T, 3>& rgb, T saturation,
    const vec<T, 3>& weights = vec<T, 3>{0.333333});

// Apply tone mapping
template <typename T>
inline vec<T, 3> tonemap(
    const vec<T, 3>& hdr, T exposure, bool filmic = false, bool srgb = true);
template <typename T>
inline vec<T, 4> tonemap(
    const vec<T, 4>& hdr, T exposure, bool filmic = false, bool srgb = true);

// Convert between CIE XYZ and RGB
template <typename T>
inline vec<T, 3> rgb_to_xyz(const vec<T, 3>& rgb);
template <typename T>
inline vec<T, 3> xyz_to_rgb(const vec<T, 3>& xyz);

// Convert between CIE XYZ and xyY
template <typename T>
inline vec<T, 3> xyz_to_xyY(const vec<T, 3>& xyz);
template <typename T>
inline vec<T, 3> xyY_to_xyz(const vec<T, 3>& xyY);

// Converts between HSV and RGB color spaces.
template <typename T>
inline vec<T, 3> hsv_to_rgb(const vec<T, 3>& hsv);
template <typename T>
inline vec<T, 3> rgb_to_hsv(const vec<T, 3>& rgb);

// Approximate color of blackbody radiation from wavelength in nm.
template <typename T>
inline vec<T, 3> blackbody_to_rgb(T temperature);

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// RANDOM NUMBER GENERATION
// -----------------------------------------------------------------------------
namespace yocto::math {

// PCG random numbers from http://www.pcg-random.org/
struct rng_state {
  uint64_t state = 0x853c49e6748fea9bULL;
  uint64_t inc   = 0xda3e39cb94b95bdbULL;

  rng_state();
  rng_state(uint64_t state, uint64_t inc);
};

// Init a random number generator with a state state from the sequence seq.
inline rng_state make_rng(uint64_t seed, uint64_t seq = 1);

// Next random numbers: floats in [0,1), ints in [0,n).
inline int   rand1i(rng_state& rng, int n);
inline vec2i rand2i(rng_state& rng, int n);
inline vec3i rand3i(rng_state& rng, int n);
inline vec4i rand4i(rng_state& rng, int n);
inline float rand1f(rng_state& rng);
inline vec2f rand2f(rng_state& rng);
inline vec3f rand3f(rng_state& rng);
inline vec4f rand4f(rng_state& rng);

// Shuffles a sequence of elements
template <typename T>
inline void shuffle(std::vector<T>& vals, rng_state& rng);

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// PERLIN NOISE FUNCTION
// -----------------------------------------------------------------------------
namespace yocto::math {

// Compute the revised Perlin noise function. Wrap provides a wrapping noise
// but must be power of two (wraps at 256 anyway). For octave based noise,
// good values are obtained with octaves=6 (numerber of noise calls),
// lacunarity=~2.0 (spacing between successive octaves: 2.0 for warpping
// output), gain=0.5 (relative weighting applied to each successive octave),
// offset=1.0 (used to invert the ridges).
template <typename T>
inline T perlin_noise(const vec<T, 3>& p, const vec<int, 3>& wrap = zero3i);
template <typename T>
inline T perlin_ridge(const vec<T, 3>& p, T lacunarity = 2, T gain = 0.5,
    int octaves = 6, T offset = 1, const vec<int, 3>& wrap = zero3i);
template <typename T>
inline T perlin_fbm(const vec<T, 3>& p, T lacunarity = 2, T gain = 0.5,
    int octaves = 6, const vec<int, 3>& wrap = zero3i);
template <typename T>
inline T perlin_turbulence(const vec<T, 3>& p, T lacunarity = 2, T gain = 0.5,
    int octaves = 6, const vec<int, 3>& wrap = zero3i);

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// SHADING FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto::math {

// Schlick approximation of the Fresnel term.
template <typename T>
inline vec<T, 3> fresnel_schlick(const vec<T, 3>& specular,
    const vec<T, 3>& normal, const vec<T, 3>& outgoing);
// Compute the fresnel term for dielectrics.
template <typename T>
inline T fresnel_dielectric(
    T eta, const vec<T, 3>& normal, const vec<T, 3>& outgoing);
// Compute the fresnel term for metals.
template <typename T>
inline vec<T, 3> fresnel_conductor(const vec<T, 3>& eta, const vec<T, 3>& etak,
    const vec<T, 3>& normal, const vec<T, 3>& outgoing);

// Convert eta to reflectivity
template <typename T>
inline vec<T, 3> eta_to_reflectivity(const vec<T, 3>& eta);
// Convert reflectivity to  eta.
template <typename T>
inline vec<T, 3> reflectivity_to_eta(const vec<T, 3>& reflectivity);
// Convert conductor eta to reflectivity.
template <typename T>
inline vec<T, 3> eta_to_reflectivity(
    const vec<T, 3>& eta, const vec<T, 3>& etak);
// Convert eta to edge tint parametrization.
template <typename T>
inline std::pair<vec<T, 3>, vec<T, 3>> eta_to_edgetint(
    const vec<T, 3>& eta, const vec<T, 3>& etak);
// Convert reflectivity and edge tint to eta.
template <typename T>
inline std::pair<vec<T, 3>, vec<T, 3>> edgetint_to_eta(
    const vec<T, 3>& reflectivity, const vec<T, 3>& edgetint);

// Evaluates the microfacet distribution.
template <typename T>
inline T microfacet_distribution(T roughness, const vec<T, 3>& normal,
    const vec<T, 3>& halfway, bool ggx = true);
// Evaluates the microfacet shadowing.
template <typename T>
inline T microfacet_shadowing(T roughness, const vec<T, 3>& normal,
    const vec<T, 3>& halfway, const vec<T, 3>& outgoing,
    const vec<T, 3>& incoming, bool ggx = true);

// Samples a microfacet distribution.
template <typename T>
inline vec<T, 3> sample_microfacet(
    T roughness, const vec<T, 3>& normal, const vec<T, 2>& rn, bool ggx = true);
// Pdf for microfacet distribution sampling.
template <typename T>
inline T sample_microfacet_pdf(T roughness, const vec<T, 3>& normal,
    const vec<T, 3>& halfway, bool ggx = true);

// Samples a microfacet distribution with the distribution of visible normals.
template <typename T>
inline vec<T, 3> sample_microfacet(T roughness, const vec<T, 3>& normal,
    const vec<T, 3>& outgoing, const vec<T, 2>& rn, bool ggx = true);
// Pdf for microfacet distribution sampling with the distribution of visible
// normals.
template <typename T>
inline T sample_microfacet_pdf(T roughness, const vec<T, 3>& normal,
    const vec<T, 3>& halfway, const vec<T, 3>& outgoing, bool ggx = true);

// Evaluates a diffuse BRDF lobe.
template <typename T>
inline vec<T, 3> eval_diffuse_reflection(const vec<T, 3>& normal,
    const vec<T, 3>& outgoing, const vec<T, 3>& incoming);
// Evaluates a specular BRDF lobe.
template <typename T>
inline vec<T, 3> eval_microfacet_reflection(T ior, T roughness,
    const vec<T, 3>& normal, const vec<T, 3>& outgoing,
    const vec<T, 3>& incoming);
// Evaluates a metal BRDF lobe.
template <typename T>
inline vec<T, 3> eval_microfacet_reflection(const vec<T, 3>& eta,
    const vec<T, 3>& etak, T roughness, const vec<T, 3>& normal,
    const vec<T, 3>& outgoing, const vec<T, 3>& incoming);
// Evaluates a transmission BRDF lobe.
template <typename T>
inline vec<T, 3> eval_microfacet_transmission(T ior, T roughness,
    const vec<T, 3>& normal, const vec<T, 3>& outgoing,
    const vec<T, 3>& incoming);
// Evaluates a refraction BRDF lobe.
template <typename T>
inline vec<T, 3> eval_microfacet_refraction(T ior, T roughness,
    const vec<T, 3>& normal, const vec<T, 3>& outgoing,
    const vec<T, 3>& incoming);

// Sample a diffuse BRDF lobe.
template <typename T>
inline vec<T, 3> sample_diffuse_reflection(
    const vec<T, 3>& normal, const vec<T, 3>& outgoing, const vec<T, 2>& rn);
// Sample a specular BRDF lobe.
template <typename T>
inline vec<T, 3> sample_microfacet_reflection(T ior, T roughness,
    const vec<T, 3>& normal, const vec<T, 3>& outgoing, const vec<T, 2>& rn);
// Sample a metal BRDF lobe.
template <typename T>
inline vec<T, 3> sample_microfacet_reflection(const vec<T, 3>& eta,
    const vec<T, 3>& etak, T roughness, const vec<T, 3>& normal,
    const vec<T, 3>& outgoing, const vec<T, 2>& rn);
// Sample a transmission BRDF lobe.
template <typename T>
inline vec<T, 3> sample_microfacet_transmission(T ior, T roughness,
    const vec<T, 3>& normal, const vec<T, 3>& outgoing, const vec<T, 2>& rn);
// Sample a refraction BRDF lobe.
template <typename T>
inline vec<T, 3> sample_microfacet_refraction(T ior, T roughness,
    const vec<T, 3>& normal, const vec<T, 3>& outgoing, T rnl,
    const vec<T, 2>& rn);

// Pdf for diffuse BRDF lobe sampling.
template <typename T>
inline T sample_diffuse_reflection_pdf(const vec<T, 3>& normal,
    const vec<T, 3>& outgoing, const vec<T, 3>& incoming);
// Pdf for specular BRDF lobe sampling.
template <typename T>
inline T sample_microfacet_reflection_pdf(T ior, T roughness,
    const vec<T, 3>& normal, const vec<T, 3>& outgoing,
    const vec<T, 3>& incoming);
// Pdf for metal BRDF lobe sampling.
template <typename T>
inline T sample_microfacet_reflection_pdf(const vec<T, 3>& eta,
    const vec<T, 3>& etak, T roughness, const vec<T, 3>& normal,
    const vec<T, 3>& outgoing, const vec<T, 3>& incoming);
// Pdf for transmission BRDF lobe sampling.
template <typename T>
inline T sample_microfacet_transmission_pdf(T ior, T roughness,
    const vec<T, 3>& normal, const vec<T, 3>& outgoing,
    const vec<T, 3>& incoming);
// Pdf for refraction BRDF lobe sampling.
template <typename T>
inline T sample_microfacet_refraction_pdf(T ior, T roughness,
    const vec<T, 3>& normal, const vec<T, 3>& outgoing,
    const vec<T, 3>& incoming);

// Evaluate a delta specular BRDF lobe.
template <typename T>
inline vec<T, 3> eval_delta_reflection(T ior, const vec<T, 3>& normal,
    const vec<T, 3>& outgoing, const vec<T, 3>& incoming);
// Evaluate a delta metal BRDF lobe.
template <typename T>
inline vec<T, 3> eval_delta_reflection(const vec<T, 3>& eta,
    const vec<T, 3>& etak, const vec<T, 3>& normal, const vec<T, 3>& outgoing,
    const vec<T, 3>& incoming);
// Evaluate a delta transmission BRDF lobe.
template <typename T>
inline vec<T, 3> eval_delta_transmission(T ior, const vec<T, 3>& normal,
    const vec<T, 3>& outgoing, const vec<T, 3>& incoming);
// Evaluate a delta refraction BRDF lobe.
template <typename T>
inline vec<T, 3> eval_delta_refraction(T ior, const vec<T, 3>& normal,
    const vec<T, 3>& outgoing, const vec<T, 3>& incoming);

// Sample a delta specular BRDF lobe.
template <typename T>
inline vec<T, 3> sample_delta_reflection(
    T ior, const vec<T, 3>& normal, const vec<T, 3>& outgoing);
// Sample a delta metal BRDF lobe.
template <typename T>
inline vec<T, 3> sample_delta_reflection(const vec<T, 3>& eta,
    const vec<T, 3>& etak, const vec<T, 3>& normal, const vec<T, 3>& outgoing);
// Sample a delta transmission BRDF lobe.
template <typename T>
inline vec<T, 3> sample_delta_transmission(
    T ior, const vec<T, 3>& normal, const vec<T, 3>& outgoing);
// Sample a delta refraction BRDF lobe.
template <typename T>
inline vec<T, 3> sample_delta_refraction(
    T ior, const vec<T, 3>& normal, const vec<T, 3>& outgoing, T rnl);

// Pdf for delta specular BRDF lobe sampling.
template <typename T>
inline T sample_delta_reflection_pdf(T ior, const vec<T, 3>& normal,
    const vec<T, 3>& outgoing, const vec<T, 3>& incoming);
// Pdf for delta metal BRDF lobe sampling.
template <typename T>
inline T sample_delta_reflection_pdf(const vec<T, 3>& eta,
    const vec<T, 3>& etak, const vec<T, 3>& normal, const vec<T, 3>& outgoing,
    const vec<T, 3>& incoming);
// Pdf for delta transmission BRDF lobe sampling.
template <typename T>
inline T sample_delta_transmission_pdf(T ior, const vec<T, 3>& normal,
    const vec<T, 3>& outgoing, const vec<T, 3>& incoming);
// Pdf for delta refraction BRDF lobe sampling.
template <typename T>
inline T sample_delta_refraction_pdf(T ior, const vec<T, 3>& normal,
    const vec<T, 3>& outgoing, const vec<T, 3>& incoming);

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// MONETACARLO SAMPLING FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto::math {

// Sample an hemispherical direction with uniform distribution.
template <typename T>
inline vec<T, 3> sample_hemisphere(const vec<T, 2>& ruv);
template <typename T>
inline T sample_hemisphere_pdf(const vec<T, 3>& direction);

// Sample an hemispherical direction with uniform distribution.
template <typename T>
inline vec<T, 3> sample_hemisphere(
    const vec<T, 3>& normal, const vec<T, 2>& ruv);
template <typename T>
inline T sample_hemisphere_pdf(
    const vec<T, 3>& normal, const vec<T, 3>& direction);

// Sample a spherical direction with uniform distribution.
template <typename T>
inline vec<T, 3> sample_sphere(const vec<T, 2>& ruv);
template <typename T>
inline T sample_sphere_pdf(const vec<T, 3>& w);

// Sample an hemispherical direction with cosine distribution.
template <typename T>
inline vec<T, 3> sample_hemisphere_cos(const vec<T, 2>& ruv);
template <typename T>
inline T sample_hemisphere_cos_pdf(const vec<T, 3>& direction);

// Sample an hemispherical direction with cosine distribution.
template <typename T>
inline vec<T, 3> sample_hemisphere_cos(
    const vec<T, 3>& normal, const vec<T, 2>& ruv);
template <typename T>
inline T sample_hemisphere_cos_pdf(
    const vec<T, 3>& normal, const vec<T, 3>& direction);

// Sample an hemispherical direction with cosine power distribution.
template <typename T>
inline vec<T, 3> sample_hemisphere_cospower(T exponent, const vec<T, 2>& ruv);
template <typename T>
inline T sample_hemisphere_cospower_pdf(T exponent, const vec<T, 3>& direction);

// Sample a point uniformly on a disk.
template <typename T>
inline vec<T, 2> sample_disk(const vec<T, 2>& ruv);
template <typename T>
inline T sample_disk_pdf(const vec<T, 2>& point);

// Sample a point uniformly on a cylinder, without caps.
template <typename T>
inline vec<T, 3> sample_cylinder(const vec<T, 2>& ruv);
template <typename T>
inline T sample_cylinder_pdf(const vec<T, 3>& point);

// Sample a point uniformly on a triangle returning the baricentric coordinates.
template <typename T>
inline vec<T, 2> sample_triangle(const vec<T, 2>& ruv);

// Sample a point uniformly on a triangle.
template <typename T>
inline vec<T, 3> sample_triangle(const vec<T, 3>& p0, const vec<T, 3>& p1,
    const vec<T, 3>& p2, const vec<T, 2>& ruv);
// Pdf for uniform triangle sampling, i.e. triangle area.
template <typename T>
inline T sample_triangle_pdf(
    const vec<T, 3>& p0, const vec<T, 3>& p1, const vec<T, 3>& p2);

// Sample an index with uniform distribution.
template <typename T>
inline int sample_uniform(int size, T r);
template <typename T = float>
inline T sample_uniform_pdf(int size);

// Sample an index with uniform distribution.
template <typename T>
inline T sample_uniform(const std::vector<T>& elements, T r);
template <typename T>
inline T sample_uniform_pdf(const std::vector<T>& elements);

// Sample a discrete distribution represented by its cdf.
template <typename T>
inline int sample_discrete(const std::vector<T>& cdf, T r);
// Pdf for uniform discrete distribution sampling.
template <typename T>
inline T sample_discrete_pdf(const std::vector<T>& cdf, int idx);

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// USER INTERFACE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto::math {

// Computes the image uv coordinates corresponding to the view parameters.
// Returns negative coordinates if out of the image.
template <typename T>
inline vec<int, 2> get_image_coords(const vec<T, 2>& mouse_pos,
    const vec<T, 2>& center, T scale, const vec<int, 2>& txt_size);

// Center image and autofit.
template <typename T>
inline void update_imview(vec<T, 2>& center, T& scale,
    const vec<int, 2>& imsize, const vec<int, 2>& winsize, bool zoom_to_fit);

// Turntable for UI navigation.
template <typename T>
inline void update_turntable(vec<T, 3>& from, vec<T, 3>& to, vec<T, 3>& up,
    const vec<T, 2>& rotate, T dolly, const vec<T, 2>& pan);

// Turntable for UI navigation.
template <typename T>
inline void update_turntable(frame<T, 3>& frame, T& focus,
    const vec<T, 2>& rotate, T dolly, const vec<T, 2>& pan);

// FPS camera for UI navigation for a frame parametrization.
template <typename T>
inline void update_fpscam(
    frame<T, 3>& frame, const vec<T, 3>& transl, const vec<T, 2>& rotate);

// Generate a ray from a camera
template <typename T>
inline ray<T, 3> camera_ray(const frame<T, 3>& frame, T lens,
    const vec<T, 2>& film, const vec<T, 2>& image_uv);

}  // namespace yocto::math

// -----------------------------------------------------------------------------
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF MATH CONSTANTS AND FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto::math {

template <typename T>
inline T abs(T a) {
  return a < 0 ? -a : a;
}
template <typename T>
inline T min(T a, T b) {
  return (a < b) ? a : b;
}
template <typename T>
inline T max(T a, T b) {
  return (a > b) ? a : b;
}
template <typename T>
inline T clamp(T a, T min_, T max_) {
  return min(max(a, min_), max_);
}
template <typename T>
inline T sign(T a) {
  return a < 0 ? -1 : 1;
}
template <typename T>
inline T sqrt(T a) {
  return std::sqrt(a);
}
template <typename T>
inline T sin(T a) {
  return std::sin(a);
}
template <typename T>
inline T cos(T a) {
  return std::cos(a);
}
template <typename T>
inline T tan(T a) {
  return std::tan(a);
}
template <typename T>
inline T asin(T a) {
  return std::asin(a);
}
template <typename T>
inline T acos(T a) {
  return std::acos(a);
}
template <typename T>
inline T atan(T a) {
  return std::atan(a);
}
template <typename T>
inline T log(T a) {
  return std::log(a);
}
template <typename T>
inline T exp(T a) {
  return std::exp(a);
}
template <typename T>
inline T log2(T a) {
  return std::log2(a);
}
template <typename T>
inline T exp2(T a) {
  return std::exp2(a);
}
template <typename T>
inline T pow(T a, T b) {
  return std::pow(a, b);
}
template <typename T>
inline T pow(T a, int b) {
  return std::pow(a, (T)b);
}
template <typename T>
inline T isfinite(T a) {
  return std::isfinite(a);
}
template <typename T>
inline T atan2(T a, T b) {
  return std::atan2(a, b);
}
template <typename T>
inline T fmod(T a, T b) {
  return std::fmod(a, b);
}
template <typename T>
inline void swap(T& a, T& b) {
  std::swap(a, b);
}
template <typename T>
inline T radians(T a) {
  return a * (T)pi / 180;
}
template <typename T>
inline T degrees(T a) {
  return a * 180 / (T)pi;
}
template <typename T>
inline T lerp(T a, T b, T u) {
  return a * (1 - u) + b * u;
}
template <typename T>
inline T step(T a, T u) {
  return u < a ? 0 : 1;
}
template <typename T>
inline T smoothstep(T a, T b, T u) {
  auto t = clamp((u - a) / (b - a), (T)0, (T)1);
  return t * t * (3 - 2 * t);
}
template <typename T>
inline T bias(T a, T bias) {
  return a / ((1 / bias - 2) * (1 - a) + 1);
}
template <typename T>
inline T gain(T a, T gain) {
  return (a < 0.5f) ? bias(a * 2, gain) / 2
                    : bias(a * 2 - 1, 1 - gain) / 2 + 0.5f;
}
template <typename T>
inline T pow2(T a) {
  return 1 << a;
}

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// VECTORS
// -----------------------------------------------------------------------------
namespace yocto::math {

// Vec2
template <typename T>
inline vec<T, 1>::vec() {}
template <typename T>
inline vec<T, 1>::vec(T v) : x{v} {}
template <typename T>
inline vec<T, 1>::operator bool() const {
  return x;
}

template <typename T>
inline T& vec<T, 1>::operator[](int i) {
  return (&x)[i];
}
template <typename T>
inline const T& vec<T, 1>::operator[](int i) const {
  return (&x)[i];
}

// Vec2
template <typename T>
inline vec<T, 2>::vec() {}
template <typename T>
inline vec<T, 2>::vec(T x, T y) : x{x}, y{y} {}
template <typename T>
inline vec<T, 2>::vec(T v) : x{v}, y{v} {}
template <typename T>
inline vec<T, 2>::operator bool() const {
  return x || y;
}

template <typename T>
inline T& vec<T, 2>::operator[](int i) {
  return (&x)[i];
}
template <typename T>
inline const T& vec<T, 2>::operator[](int i) const {
  return (&x)[i];
}

// Vec3
template <typename T>
inline vec<T, 3>::vec() {}
template <typename T>
inline vec<T, 3>::vec(T x, T y, T z) : x{x}, y{y}, z{z} {}
template <typename T>
inline vec<T, 3>::vec(const vec<T, 2>& v, T z) : x{v.x}, y{v.y}, z{z} {}
template <typename T>
inline vec<T, 3>::vec(T v) : x{v}, y{v}, z{v} {}
template <typename T>
inline vec<T, 3>::operator bool() const {
  return x || y || z;
}

template <typename T>
inline T& vec<T, 3>::operator[](int i) {
  return (&x)[i];
}
template <typename T>
inline const T& vec<T, 3>::operator[](int i) const {
  return (&x)[i];
}

// Vec4
template <typename T>
inline vec<T, 4>::vec() {}
template <typename T>
inline vec<T, 4>::vec(T x, T y, T z, T w) : x{x}, y{y}, z{z}, w{w} {}
template <typename T>
inline vec<T, 4>::vec(const vec<T, 3>& v, T w) : x{v.x}, y{v.y}, z{v.z}, w{w} {}
template <typename T>
inline vec<T, 4>::vec(T v) : x{v}, y{v}, z{v}, w{v} {}
template <typename T>
inline vec<T, 4>::operator bool() const {
  return x || y || z || w;
}

template <typename T>
inline T& vec<T, 4>::operator[](int i) {
  return (&x)[i];
}
template <typename T>
inline const T& vec<T, 4>::operator[](int i) const {
  return (&x)[i];
}

// Element access
template <typename T>
inline vec<T, 3>& xyz(vec<T, 4>& a) {
  return (vec<T, 3>&)a;
}
template <typename T>
inline const vec<T, 3>& xyz(const vec<T, 4>& a) {
  return (const vec<T, 3>&)a;
}

// Vector comparison operations.
template <typename T, size_t N>
inline bool operator==(const vec<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return a.x == b.x;
  } else if constexpr (N == 2) {
    return a.x == b.x || a.y == b.y;
  } else if constexpr (N == 3) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
  } else if constexpr (N == 4) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N>
inline bool operator!=(const vec<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return a.x != b.x;
  } else if constexpr (N == 2) {
    return a.x != b.x || a.y != b.y;
  } else if constexpr (N == 3) {
    return a.x != b.x || a.y != b.y || a.z != b.z;
  } else if constexpr (N == 4) {
    return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w;
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}

// Vector operations.
template <typename T, size_t N>
inline vec<T, N> operator+(const vec<T, N>& a) {
  return a;
}
template <typename T, size_t N>
inline vec<T, N> operator-(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return {-a.x};
  } else if constexpr (N == 2) {
    return {-a.x, -a.y};
  } else if constexpr (N == 3) {
    return {-a.x, -a.y, -a.z};
  } else if constexpr (N == 4) {
    return {-a.x, -a.y, -a.z, -a.w};
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N>
inline vec<T, N> operator+(const vec<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return {a.x + b.x};
  } else if constexpr (N == 2) {
    return {a.x + b.x, a.y + b.y};
  } else if constexpr (N == 3) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
  } else if constexpr (N == 4) {
    return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N, typename T1>
inline vec<T, N> operator+(const vec<T, N>& a, T1 b) {
  if constexpr (N == 1) {
    return {a.x + b};
  } else if constexpr (N == 2) {
    return {a.x + b, a.y + b};
  } else if constexpr (N == 3) {
    return {a.x + b, a.y + b, a.z + b};
  } else if constexpr (N == 4) {
    return {a.x + b, a.y + b, a.z + b, a.w + b};
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N, typename T1>
inline vec<T, N> operator+(T1 a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return {a + b.x};
  } else if constexpr (N == 2) {
    return {a + b.x, a + b.y};
  } else if constexpr (N == 3) {
    return {a + b.x, a + b.y, a + b.z};
  } else if constexpr (N == 4) {
    return {a + b.x, a + b.y, a + b.z, a + b.w};
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N>
inline vec<T, N> operator-(const vec<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return {a.x - b.x};
  } else if constexpr (N == 2) {
    return {a.x - b.x, a.y - b.y};
  } else if constexpr (N == 3) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
  } else if constexpr (N == 4) {
    return {a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w};
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N, typename T1>
inline vec<T, N> operator-(const vec<T, N>& a, T1 b) {
  if constexpr (N == 1) {
    return {a.x - b};
  } else if constexpr (N == 2) {
    return {a.x - b, a.y - b};
  } else if constexpr (N == 3) {
    return {a.x - b, a.y - b, a.z - b};
  } else if constexpr (N == 4) {
    return {a.x - b, a.y - b, a.z - b, a.w - b};
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N, typename T1>
inline vec<T, N> operator-(T1 a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return {a - b.x};
  } else if constexpr (N == 2) {
    return {a - b.x, a - b.y};
  } else if constexpr (N == 3) {
    return {a - b.x, a - b.y, a - b.z};
  } else if constexpr (N == 4) {
    return {a - b.x, a - b.y, a - b.z, a - b.w};
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N>
inline vec<T, N> operator*(const vec<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return {a.x * b.x};
  } else if constexpr (N == 2) {
    return {a.x * b.x, a.y * b.y};
  } else if constexpr (N == 3) {
    return {a.x * b.x, a.y * b.y, a.z * b.z};
  } else if constexpr (N == 4) {
    return {a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w};
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N, typename T1>
inline vec<T, N> operator*(const vec<T, N>& a, T1 b) {
  if constexpr (N == 1) {
    return {a.x * b};
  } else if constexpr (N == 2) {
    return {a.x * b, a.y * b};
  } else if constexpr (N == 3) {
    return {a.x * b, a.y * b, a.z * b};
  } else if constexpr (N == 4) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N, typename T1>
inline vec<T, N> operator*(T1 a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return {a * b.x};
  } else if constexpr (N == 2) {
    return {a * b.x, a * b.y};
  } else if constexpr (N == 3) {
    return {a * b.x, a * b.y, a * b.z};
  } else if constexpr (N == 4) {
    return {a * b.x, a * b.y, a * b.z, a * b.w};
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N>
inline vec<T, N> operator/(const vec<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return {a.x / b.x};
  } else if constexpr (N == 2) {
    return {a.x / b.x, a.y / b.y};
  } else if constexpr (N == 3) {
    return {a.x / b.x, a.y / b.y, a.z / b.z};
  } else if constexpr (N == 4) {
    return {a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w};
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N, typename T1>
inline vec<T, N> operator/(const vec<T, N>& a, T1 b) {
  if constexpr (N == 1) {
    return {a.x / b};
  } else if constexpr (N == 2) {
    return {a.x / b, a.y / b};
  } else if constexpr (N == 3) {
    return {a.x / b, a.y / b, a.z / b};
  } else if constexpr (N == 4) {
    return {a.x / b, a.y / b, a.z / b, a.w / b};
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N, typename T1>
inline vec<T, N> operator/(T1 a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return {a / b.x};
  } else if constexpr (N == 2) {
    return {a / b.x, a / b.y};
  } else if constexpr (N == 3) {
    return {a / b.x, a / b.y, a / b.z};
  } else if constexpr (N == 4) {
    return {a / b.x, a / b.y, a / b.z, a / b.w};
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}

// Vector assignments
template <typename T, size_t N>
inline vec<T, N>& operator+=(vec<T, N>& a, const vec<T, N>& b) {
  return a = a + b;
}
template <typename T, size_t N, typename T1>
inline vec<T, N>& operator+=(vec<T, N>& a, T1 b) {
  return a = a + b;
}
template <typename T, size_t N>
inline vec<T, N>& operator-=(vec<T, N>& a, const vec<T, N>& b) {
  return a = a - b;
}
template <typename T, size_t N, typename T1>
inline vec<T, N>& operator-=(vec<T, N>& a, T1 b) {
  return a = a - b;
}
template <typename T, size_t N>
inline vec<T, N>& operator*=(vec<T, N>& a, const vec<T, N>& b) {
  return a = a * b;
}
template <typename T, size_t N, typename T1>
inline vec<T, N>& operator*=(vec<T, N>& a, T1 b) {
  return a = a * b;
}
template <typename T, size_t N>
inline vec<T, N>& operator/=(vec<T, N>& a, const vec<T, N>& b) {
  return a = a / b;
}
template <typename T, size_t N, typename T1>
inline vec<T, N>& operator/=(vec<T, N>& a, T1 b) {
  return a = a / b;
}

// Vector products and lengths.
template <typename T, size_t N>
inline T dot(const vec<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return a.x * b.x;
  } else if constexpr (N == 2) {
    return a.x * b.x + a.y * b.y;
  } else if constexpr (N == 3) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
  } else if constexpr (N == 4) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N>
inline T length(const vec<T, N>& a) {
  return sqrt(dot(a, a));
}
template <typename T, size_t N>
inline vec<T, N> normalize(const vec<T, N>& a) {
  auto l = length(a);
  return (l != 0) ? a / l : a;
}
template <typename T, size_t N>
inline T distance(const vec<T, N>& a, const vec<T, N>& b) {
  return length(a - b);
}
template <typename T, size_t N>
inline T distance_squared(const vec<T, N>& a, const vec<T, N>& b) {
  return dot(a - b, a - b);
}
template <typename T>
inline T cross(const vec<T, 2>& a, const vec<T, 2>& b) {
  return a.x * b.y - a.y * b.x;
}
template <typename T>
inline vec<T, 3> cross(const vec<T, 3>& a, const vec<T, 3>& b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
template <typename T>
inline T angle(const vec<T, 3>& a, const vec<T, 3>& b) {
  return acos(clamp(dot(normalize(a), normalize(b)), (float)-1, (float)1));
}

// Orthogonal vectors.
template <typename T>
inline vec<T, 3> orthogonal(const vec<T, 3>& v) {
  // http://lolengine.net/blog/2013/09/21/picking-orthogonal-std::vector-combing-coconuts)
  return abs(v.x) > abs(v.z) ? vec<T, 3>{-v.y, v.x, 0}
                             : vec<T, 3>{0, -v.z, v.y};
}
template <typename T>
inline vec<T, 3> orthonormalize(const vec<T, 3>& a, const vec<T, 3>& b) {
  return normalize(a - b * dot(a, b));
}

// Reflected and refracted vector.
template <typename T>
inline vec<T, 3> reflect(const vec<T, 3>& w, const vec<T, 3>& n) {
  return -w + 2 * dot(n, w) * n;
}
template <typename T>
inline vec<T, 3> refract(const vec<T, 3>& w, const vec<T, 3>& n, T inv_eta) {
  auto cosine = dot(n, w);
  auto k      = 1 + inv_eta * inv_eta * (cosine * cosine - 1);
  if (k < 0) return {0, 0, 0};  // tir
  return -w * inv_eta + (inv_eta * cosine - sqrt(k)) * n;
}

// Slerp
template <typename T, size_t N>
inline vec<T, N> slerp(const vec<T, N>& a, const vec<T, N>& b, T u) {
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
template <typename T, size_t N, typename T1>
inline vec<T, N> max(const vec<T, N>& a, T1 b) {
  if constexpr (N == 1) {
    return {max(a.x, (T)b)};
  } else if constexpr (N == 2) {
    return {max(a.x, (T)b), max(a.y, (T)b)};
  } else if constexpr (N == 3) {
    return {max(a.x, (T)b), max(a.y, (T)b), max(a.z, (T)b)};
  } else if constexpr (N == 4) {
    return {max(a.x, (T)b), max(a.y, (T)b), max(a.z, (T)b), max(a.w, (T)b)};
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N, typename T1>
inline vec<T, N> min(const vec<T, N>& a, T1 b) {
  if constexpr (N == 1) {
    return {min(a.x, (T)b)};
  } else if constexpr (N == 2) {
    return {min(a.x, (T)b), min(a.y, (T)b)};
  } else if constexpr (N == 3) {
    return {min(a.x, (T)b), min(a.y, (T)b), min(a.z, (T)b)};
  } else if constexpr (N == 4) {
    return {min(a.x, (T)b), min(a.y, (T)b), min(a.z, (T)b), min(a.w, (T)b)};
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N>
inline vec<T, N> max(const vec<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return {max(a.x, b.x)};
  } else if constexpr (N == 2) {
    return {max(a.x, b.x), max(a.y, b.y)};
  } else if constexpr (N == 3) {
    return {max(a.x, b.x), max(a.y, b.y), max(a.z, b.z)};
  } else if constexpr (N == 4) {
    return {max(a.x, b.x), max(a.y, b.y), max(a.z, b.z), max(a.w, b.w)};
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N>
inline vec<T, N> min(const vec<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return {min(a.x, b.x)};
  } else if constexpr (N == 2) {
    return {min(a.x, b.x), min(a.y, b.y)};
  } else if constexpr (N == 3) {
    return {min(a.x, b.x), min(a.y, b.y), min(a.z, b.z)};
  } else if constexpr (N == 4) {
    return {min(a.x, b.x), min(a.y, b.y), min(a.z, b.z), min(a.w, b.w)};
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N, typename T1, typename T2>
inline vec<T, N> clamp(const vec<T, N>& a, T1 min, T2 max) {
  if constexpr (N == 1) {
    return {clamp(a.x, (T)min, (T)max)};
  } else if constexpr (N == 2) {
    return {clamp(a.x, (T)min, (T)max), clamp(a.y, (T)min, (T)max)};
  } else if constexpr (N == 3) {
    return {clamp(a.x, (T)min, (T)max), clamp(a.y, (T)min, (T)max),
        clamp(a.z, (T)min, (T)max)};
  } else if constexpr (N == 4) {
    return {clamp(a.x, (T)min, (T)max), clamp(a.y, (T)min, (T)max),
        clamp(a.z, (T)min, (T)max), clamp(a.w, (T)min, (T)max)};
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N>
inline vec<T, N> lerp(const vec<T, N>& a, const vec<T, N>& b, T u) {
  return a * (1 - u) + b * u;
}
template <typename T, size_t N>
inline vec<T, N> lerp(
    const vec<T, N>& a, const vec<T, N>& b, const vec<T, N>& u) {
  return a * (1 - u) + b * u;
}

template <typename T, size_t N>
inline T max(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return a.x;
  } else if constexpr (N == 2) {
    return max(a.x, a.y);
  } else if constexpr (N == 3) {
    return max(max(a.x, a.y), a.z);
  } else if constexpr (N == 4) {
    return max(max(max(a.x, a.y), a.z), a.w);
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N>
inline T min(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return a.x;
  } else if constexpr (N == 2) {
    return min(a.x, a.y);
  } else if constexpr (N == 3) {
    return min(min(a.x, a.y), a.z);
  } else if constexpr (N == 4) {
    return min(min(min(a.x, a.y), a.z), a.w);
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N>
inline T sum(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return a.x;
  } else if constexpr (N == 2) {
    return a.x + a.y;
  } else if constexpr (N == 3) {
    return a.x + a.y + a.z;
  } else if constexpr (N == 4) {
    return a.x + a.y + a.z + a.w;
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N>
inline T mean(const vec<T, N>& a) {
  return sum(a) / (T)N;
}

// Functions applied to std::vector elements
template <typename T, size_t N>
inline vec<T, N> abs(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return {abs(a.x)};
  } else if constexpr (N == 2) {
    return {abs(a.x), abs(a.y)};
  } else if constexpr (N == 3) {
    return {abs(a.x), abs(a.y), abs(a.z)};
  } else if constexpr (N == 4) {
    return {abs(a.x), abs(a.y), abs(a.z), abs(a.w)};
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N>
inline vec<T, N> sqrt(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return {sqrt(a.x)};
  } else if constexpr (N == 2) {
    return {sqrt(a.x), sqrt(a.y)};
  } else if constexpr (N == 3) {
    return {sqrt(a.x), sqrt(a.y), sqrt(a.z)};
  } else if constexpr (N == 4) {
    return {sqrt(a.x), sqrt(a.y), sqrt(a.z), sqrt(a.w)};
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N>
inline vec<T, N> exp(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return {exp(a.x)};
  } else if constexpr (N == 2) {
    return {exp(a.x), exp(a.y)};
  } else if constexpr (N == 3) {
    return {exp(a.x), exp(a.y), exp(a.z)};
  } else if constexpr (N == 4) {
    return {exp(a.x), exp(a.y), exp(a.z), exp(a.w)};
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N>
inline vec<T, N> log(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return {log(a.x)};
  } else if constexpr (N == 2) {
    return {log(a.x), log(a.y)};
  } else if constexpr (N == 3) {
    return {log(a.x), log(a.y), log(a.z)};
  } else if constexpr (N == 4) {
    return {log(a.x), log(a.y), log(a.z), log(a.w)};
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N>
inline vec<T, N> exp2(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return {exp2(a.x)};
  } else if constexpr (N == 2) {
    return {exp2(a.x), exp2(a.y)};
  } else if constexpr (N == 3) {
    return {exp2(a.x), exp2(a.y), exp2(a.z)};
  } else if constexpr (N == 4) {
    return {exp2(a.x), exp2(a.y), exp2(a.z), exp2(a.w)};
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N>
inline vec<T, N> log2(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return {log2(a.x)};
  } else if constexpr (N == 2) {
    return {log2(a.x), log2(a.y)};
  } else if constexpr (N == 3) {
    return {log2(a.x), log2(a.y), log2(a.z)};
  } else if constexpr (N == 4) {
    return {log2(a.x), log2(a.y), log2(a.z), log2(a.w)};
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N>
inline vec<T, N> pow(const vec<T, N>& a, T b) {
  if constexpr (N == 1) {
    return {pow(a.x, b)};
  } else if constexpr (N == 2) {
    return {pow(a.x, b), pow(a.y, b)};
  } else if constexpr (N == 3) {
    return {pow(a.x, b), pow(a.y, b), pow(a.z, b)};
  } else if constexpr (N == 4) {
    return {pow(a.x, b), pow(a.y, b), pow(a.z, b), pow(a.w, b)};
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N>
inline vec<T, N> pow(const vec<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return {pow(a.x, b.x)};
  } else if constexpr (N == 2) {
    return {pow(a.x, b.x), pow(a.y, b.y)};
  } else if constexpr (N == 3) {
    return {pow(a.x, b.x), pow(a.y, b.y), pow(a.z, b.z)};
  } else if constexpr (N == 4) {
    return {pow(a.x, b.x), pow(a.y, b.y), pow(a.z, b.z), pow(a.w, b.w)};
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N>
inline vec<T, N> gain(const vec<T, N>& a, T b) {
  if constexpr (N == 1) {
    return {gain(a.x, b)};
  } else if constexpr (N == 2) {
    return {gain(a.x, b), gain(a.y, b)};
  } else if constexpr (N == 3) {
    return {gain(a.x, b), gain(a.y, b), gain(a.z, b)};
  } else if constexpr (N == 4) {
    return {gain(a.x, b), gain(a.y, b), gain(a.z, b), gain(a.w, b)};
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N>
inline bool isfinite(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return isfinite(a.x);
  } else if constexpr (N == 2) {
    return isfinite(a.x) && isfinite(a.y);
  } else if constexpr (N == 3) {
    return isfinite(a.x) && isfinite(a.y) && isfinite(a.z);
  } else if constexpr (N == 4) {
    return isfinite(a.x) && isfinite(a.y) && isfinite(a.z) && isfinite(a.w);
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, size_t N>
inline void swap(vec<T, N>& a, vec<T, N>& b) {
  std::swap(a, b);
}

// Quaternion operatons represented as xi + yj + zk + w
// const auto identity_quat4f = vec4f{0, 0, 0, 1};
template <typename T, size_t N, typename T1>
inline vec<T, N> quat_mul(const vec<T, N>& a, T1 b) {
  return {a.x * b, a.y * b, a.z * b, a.w * b};
}
template <typename T, size_t N>
inline vec<T, N> quat_mul(const vec<T, N>& a, const vec<T, N>& b) {
  return {a.x * b.w + a.w * b.x + a.y * b.w - a.z * b.y,
      a.y * b.w + a.w * b.y + a.z * b.x - a.x * b.z,
      a.z * b.w + a.w * b.z + a.x * b.y - a.y * b.x,
      a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z};
}
template <typename T, size_t N>
inline vec<T, N> quat_conjugate(const vec<T, N>& a) {
  return {-a.x, -a.y, -a.z, a.w};
}
template <typename T, size_t N>
inline vec<T, N> quat_inverse(const vec<T, N>& a) {
  return quat_conjugate(a) / dot(a, a);
}

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// IMPLEMENRTATION OF VECTOR HASHING
// -----------------------------------------------------------------------------
namespace std {

// Hash functor for std::vector for use with hash_map
template <typename T, size_t N>
struct hash<yocto::math::vec<T, N>> {
  size_t operator()(const yocto::math::vec<T, N>& v) const {
    static const auto hasher = std::hash<int>();
    auto              h      = (size_t)0;
    if constexpr (N == 1) {
      h ^= hasher(v.x) + 0x9e3779b9 + (h << 6) + (h >> 2);
    } else if constexpr (N == 2) {
      h ^= hasher(v.x) + 0x9e3779b9 + (h << 6) + (h >> 2);
      h ^= hasher(v.y) + 0x9e3779b9 + (h << 6) + (h >> 2);
    } else if constexpr (N == 3) {
      h ^= hasher(v.x) + 0x9e3779b9 + (h << 6) + (h >> 2);
      h ^= hasher(v.y) + 0x9e3779b9 + (h << 6) + (h >> 2);
      h ^= hasher(v.z) + 0x9e3779b9 + (h << 6) + (h >> 2);
    } else if constexpr (N == 4) {
      h ^= hasher(v.x) + 0x9e3779b9 + (h << 6) + (h >> 2);
      h ^= hasher(v.y) + 0x9e3779b9 + (h << 6) + (h >> 2);
      h ^= hasher(v.z) + 0x9e3779b9 + (h << 6) + (h >> 2);
      h ^= hasher(v.w) + 0x9e3779b9 + (h << 6) + (h >> 2);
    } else {
      static_assert(N >= 1 && N <= 4, "vector size unsupported");
    }
    return h;
  }
};

}  // namespace std

// -----------------------------------------------------------------------------
// MATRICES
// -----------------------------------------------------------------------------
namespace yocto::math {

// Small Fixed-size matrices stored in column major format.
template <typename T, size_t M>
inline mat<T, M, 1>::mat() {}
template <typename T, size_t M>
inline mat<T, M, 1>::mat(const vec<T, M>& x) : x{x} {}

template <typename T, size_t M>
inline vec<T, M>& mat<T, M, 1>::operator[](int i) {
  return (&x)[i];
}
template <typename T, size_t M>
inline const vec<T, M>& mat<T, M, 1>::operator[](int i) const {
  return (&x)[i];
}

// Small Fixed-size matrices stored in column major format.
template <typename T, size_t M>
inline mat<T, M, 2>::mat() {}
template <typename T, size_t M>
inline mat<T, M, 2>::mat(const vec<T, M>& x, const vec<T, M>& y) : x{x}, y{y} {}

template <typename T, size_t M>
inline vec<T, M>& mat<T, M, 2>::operator[](int i) {
  return (&x)[i];
}
template <typename T, size_t M>
inline const vec<T, M>& mat<T, M, 2>::operator[](int i) const {
  return (&x)[i];
}

// Small Fixed-size matrices stored in column major format.
template <typename T, size_t M>
inline mat<T, M, 3>::mat() {}
template <typename T, size_t M>
inline mat<T, M, 3>::mat(
    const vec<T, M>& x, const vec<T, M>& y, const vec<T, M>& z)
    : x{x}, y{y}, z{z} {}

template <typename T, size_t M>
inline vec<T, M>& mat<T, M, 3>::operator[](int i) {
  return (&x)[i];
}
template <typename T, size_t M>
inline const vec<T, M>& mat<T, M, 3>::operator[](int i) const {
  return (&x)[i];
}

// Small Fixed-size matrices stored in column major format.
template <typename T, size_t M>
inline mat<T, M, 4>::mat() {}
template <typename T, size_t M>
inline mat<T, M, 4>::mat(const vec<T, M>& x, const vec<T, M>& y,
    const vec<T, M>& z, const vec<T, M>& w)
    : x{x}, y{y}, z{z}, w{w} {}

template <typename T, size_t M>
inline vec<T, M>& mat<T, M, 4>::operator[](int i) {
  return (&x)[i];
}
template <typename T, size_t M>
inline const vec<T, M>& mat<T, M, 4>::operator[](int i) const {
  return (&x)[i];
}

// Matrix comparisons.
template <typename T, size_t M, size_t N>
inline bool operator==(const mat<T, M, N>& a, const mat<T, M, N>& b) {
  if constexpr (N == 1) {
    return a.x == b.x;
  } else if constexpr (N == 2) {
    return a.x == b.x || a.y == b.y;
  } else if constexpr (N == 3) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
  } else if constexpr (N == 4) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
  } else {
    static_assert(N >= 0 || N <= 4, "matrix size unsupported");
  }
}
template <typename T, size_t M, size_t N>
inline bool operator!=(const mat<T, M, N>& a, const mat<T, M, N>& b) {
  if constexpr (N == 1) {
    return a.x != b.x;
  } else if constexpr (N == 2) {
    return a.x != b.x || a.y != b.y;
  } else if constexpr (N == 3) {
    return a.x != b.x || a.y != b.y || a.z != b.z;
  } else if constexpr (N == 4) {
    return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w;
  } else {
    static_assert(N >= 0 || N <= 4, "matrix size unsupported");
  }
}

// Matrix operations.
template <typename T, size_t M, size_t N>
inline mat<T, M, N> operator+(const mat<T, M, N>& a, const mat<T, M, N>& b) {
  if constexpr (N == 1) {
    return {a.x + b.x};
  } else if constexpr (N == 2) {
    return {a.x + b.x, a.y + b.y};
  } else if constexpr (N == 3) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
  } else if constexpr (N == 4) {
    return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
  } else {
    static_assert(N >= 0 || N <= 4, "matrix size unsupported");
  }
}
template <typename T, size_t M, size_t N, typename T1>
inline mat<T, M, N> operator*(const mat<T, M, N>& a, T1 b) {
  if constexpr (N == 1) {
    return {a.x * b};
  } else if constexpr (N == 2) {
    return {a.x * b, a.y * b};
  } else if constexpr (N == 3) {
    return {a.x * b, a.y * b, a.z * b};
  } else if constexpr (N == 4) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
  } else {
    static_assert(N >= 0 || N <= 4, "matrix size unsupported");
  }
}
template <typename T, size_t M, size_t N>
inline vec<T, M> operator*(const mat<T, M, N>& a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return a.x * b.x;
  } else if constexpr (N == 2) {
    return a.x * b.x + a.y * b.y;
  } else if constexpr (N == 3) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
  } else if constexpr (N == 4) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
  } else {
    static_assert(N >= 0 || N <= 4, "matrix size unsupported");
  }
}
template <typename T, size_t M, size_t N>
inline vec<T, N> operator*(const vec<T, M>& a, const mat<T, M, N>& b) {
  if constexpr (N == 1) {
    return {dot(a, b.x)};
  } else if constexpr (N == 2) {
    return {dot(a, b.x), dot(a, b.y)};
  } else if constexpr (N == 3) {
    return {dot(a, b.x), dot(a, b.y), dot(a, b.z)};
  } else if constexpr (N == 4) {
    return {dot(a, b.x), dot(a, b.y), dot(a, b.z), dot(a, b.w)};
  } else {
    static_assert(N >= 0 || N <= 4, "matrix size unsupported");
  }
}
template <typename T, size_t M, size_t N, size_t K>
inline mat<T, M, K> operator*(const mat<T, M, N>& a, const mat<T, N, K>& b) {
  if constexpr (N == 1) {
    return {a * b.x};
  } else if constexpr (N == 2) {
    return {a * b.x, a * b.y};
  } else if constexpr (N == 3) {
    return {a * b.x, a * b.y, a * b.z};
  } else if constexpr (N == 4) {
    return {a * b.x, a * b.y, a * b.z, a * b.w};
  } else {
    static_assert(N >= 0 || N <= 4, "matrix size unsupported");
  }
}

// Matrix assignments.
template <typename T, size_t M, size_t N>
inline mat<T, M, N>& operator+=(mat<T, M, N>& a, const mat<T, M, N>& b) {
  return a = a + b;
}
template <typename T, size_t M, size_t N>
inline mat<T, M, N>& operator*=(mat<T, M, N>& a, const mat<T, N, N>& b) {
  return a = a * b;
}
template <typename T, size_t M, size_t N, typename T1>
inline mat<T, M, N>& operator*=(mat<T, M, N>& a, T1 b) {
  return a = a * b;
}

// Matrix diagonals and transposes.
template <typename T, size_t N>
inline vec<T, N> diagonal(const mat<T, N, N>& a) {
  if constexpr (N == 1) {
    return {a.x.x};
  } else if constexpr (N == 2) {
    return {a.x.x, a.y.y};
  } else if constexpr (N == 3) {
    return {a.x.x, a.y.y, a.z.z};
  } else if constexpr (N == 4) {
    return {a.x.x, a.y.y, a.z.z, a.w.w};
  } else {
    static_assert(N >= 0 || N <= 4, "matrix size unsupported");
  }
}
template <typename T, size_t N>
inline mat<T, N, N> transpose(const mat<T, N, N>& a) {
  if constexpr (N == 1) {
    return {{a.x.x, a.y.x}};
  } else if constexpr (N == 2) {
    return {{a.x.x, a.y.x}, {a.x.y, a.y.y}};
  } else if constexpr (N == 3) {
    return {
        {a.x.x, a.y.x, a.z.x},
        {a.x.y, a.y.y, a.z.y},
        {a.x.z, a.y.z, a.z.z},
    };
  } else if constexpr (N == 4) {
    return {
        {a.x.x, a.y.x, a.z.x, a.w.x},
        {a.x.y, a.y.y, a.z.y, a.w.y},
        {a.x.z, a.y.z, a.z.z, a.w.z},
        {a.x.w, a.y.w, a.z.w, a.w.w},
    };
  } else {
    static_assert(N >= 0 || N <= 4, "matrix size unsupported");
  }
}

// Matrix adjoints, determinants and inverses.
template <typename T, size_t N>
inline T determinant(const mat<T, N, N>& a) {
  if constexpr (N == 1) {
    return a.x.x;
  } else if constexpr (N == 2) {
    return cross(a.x, a.y);
  } else if constexpr (N == 3) {
    return dot(a.x, cross(a.y, a.z));
  } else {
    static_assert(N >= 0 || N <= 3, "matrix size unsupported");
  }
}
template <typename T, size_t N>
inline mat<T, N, N> adjoint(const mat<T, N, N>& a) {
  if constexpr (N == 1) {
    return {a.x.x};
  } else if constexpr (N == 2) {
    return {{a.y.y, -a.x.y}, {-a.y.x, a.x.x}};
  } else if constexpr (N == 3) {
    return transpose(mat3f{cross(a.y, a.z), cross(a.z, a.x), cross(a.x, a.y)});
  } else {
    static_assert(N >= 0 || N <= 3, "matrix size unsupported");
  }
}
template <typename T, size_t N>
inline mat<T, N, N> inverse(const mat<T, N, N>& a) {
  return adjoint(a) * (1 / determinant(a));
}

// Constructs a basis from a direction
template <typename T>
inline mat<T, 3, 3> basis_fromz(const vec<T, 3>& v) {
  // https://graphics.pixar.com/library/OrthonormalB/paper.pdf
  auto z    = normalize(v);
  auto sign = copysignf((T)1, z.z);
  auto a    = -1 / (sign + z.z);
  auto b    = z.x * z.y * a;
  auto x    = vec<T, 3>{1 + sign * z.x * z.x * a, sign * b, -sign * z.x};
  auto y    = vec<T, 3>{b, sign + z.y * z.y * a, -z.y};
  return {x, y, z};
}

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// RIGID BODY TRANSFORMS/FRAMES
// -----------------------------------------------------------------------------
namespace yocto::math {

// Rigid frames stored as a column-major affine transform matrix.
template <typename T>
inline frame<T, 1>::frame() {}
template <typename T>
inline frame<T, 1>::frame(const vec<T, 1>& x, const vec<T, 1>& o)
    : x{x}, o{o} {}
template <typename T>
inline frame<T, 1>::frame(const vec<T, 1>& o) : x{1, 0}, o{o} {}
template <typename T>
inline frame<T, 1>::frame(const mat<T, 1, 1>& m, const vec<T, 1>& t)
    : x{m.x}, o{t} {}
template <typename T>
inline frame<T, 1>::frame(const mat<T, 2, 2>& m)
    : x{m.x.x, m.x.y}, o{m.z.x, m.z.y} {}
template <typename T>
inline frame<T, 1>::operator mat<T, 2, 2>() const {
  return {{x, 0}, {o, 1}};
}

template <typename T>
inline vec<T, 1>& frame<T, 1>::operator[](int i) {
  return (&x)[i];
}
template <typename T>
inline const vec<T, 1>& frame<T, 1>::operator[](int i) const {
  return (&x)[i];
}

// Rigid frames stored as a column-major affine transform matrix.
template <typename T>
inline frame<T, 2>::frame() {}
template <typename T>
inline frame<T, 2>::frame(
    const vec<T, 2>& x, const vec<T, 2>& y, const vec<T, 2>& o)
    : x{x}, y{y}, o{o} {}
template <typename T>
inline frame<T, 2>::frame(const vec<T, 2>& o) : x{1, 0}, y{0, 1}, o{o} {}
template <typename T>
inline frame<T, 2>::frame(const mat<T, 2, 2>& m, const vec<T, 2>& t)
    : x{m.x}, y{m.y}, o{t} {}
template <typename T>
inline frame<T, 2>::frame(const mat<T, 3, 3>& m)
    : x{m.x.x, m.x.y}, y{m.y.x, m.y.y}, o{m.z.x, m.z.y} {}
template <typename T>
inline frame<T, 2>::operator mat<T, 3, 3>() const {
  return {{x, 0}, {y, 0}, {o, 1}};
}

template <typename T>
inline vec<T, 2>& frame<T, 2>::operator[](int i) {
  return (&x)[i];
}
template <typename T>
inline const vec<T, 2>& frame<T, 2>::operator[](int i) const {
  return (&x)[i];
}

// Rigid frames stored as a column-major affine transform matrix.
template <typename T>
inline frame<T, 3>::frame() {}
template <typename T>
inline frame<T, 3>::frame(const vec<T, 3>& x, const vec<T, 3>& y,
    const vec<T, 3>& z, const vec<T, 3>& o)
    : x{x}, y{y}, z{z}, o{o} {}
template <typename T>
inline frame<T, 3>::frame(const vec<T, 3>& o)
    : x{1, 0, 0}, y{0, 1, 0}, z{0, 0, 1}, o{o} {}
template <typename T>
inline frame<T, 3>::frame(const mat<T, 3, 3>& m, const vec<T, 3>& t)
    : x{m.x}, y{m.y}, z{m.z}, o{t} {}
template <typename T>
inline frame<T, 3>::frame(const mat<T, 4, 4>& m)
    : x{m.x.x, m.x.y, m.x.z}
    , y{m.y.x, m.y.y, m.y.z}
    , z{m.z.x, m.z.y, m.z.z}
    , o{m.w.x, m.w.y, m.w.z} {}
template <typename T>
inline frame<T, 3>::operator mat<T, 4, 4>() const {
  return {{x, 0}, {y, 0}, {z, 0}, {o, 1}};
}

template <typename T>
inline vec<T, 3>& frame<T, 3>::operator[](int i) {
  return (&x)[i];
}
template <typename T>
inline const vec<T, 3>& frame<T, 3>::operator[](int i) const {
  return (&x)[i];
}

// Frame properties
template <typename T, size_t N>
inline const mat<T, N, N>& rotation(const frame<T, N>& a) {
  return (const mat<T, N, N>&)a;
}

// Frame comparisons.
template <typename T, size_t N>
inline bool operator==(const frame<T, N>& a, const frame<T, N>& b) {
  if constexpr (N == 1) {
    return a.x == b.x && a.o == b.o;
  } else if constexpr (N == 2) {
    return a.x == b.x && a.y == b.y && a.o == b.o;
  } else if constexpr (N == 3) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.o == b.o;
  } else {
    static_assert(N >= 0 || N <= 3, "frame size unsupported");
  }
}
template <typename T, size_t N>
inline bool operator!=(const frame<T, N>& a, const frame<T, N>& b) {
  if constexpr (N == 1) {
    return a.x != b.x || a.o != b.o;
  } else if constexpr (N == 2) {
    return a.x != b.x || a.y != b.y || a.o != b.o;
  } else if constexpr (N == 3) {
    return a.x != b.x || a.y != b.y || a.z != b.z || a.o != b.o;
  } else {
    static_assert(N >= 0 || N <= 3, "frame size unsupported");
  }
}

// Frame composition, equivalent to affine matrix product.
template <typename T, size_t N>
inline frame<T, N> operator*(const frame<T, N>& a, const frame<T, N>& b) {
  return {rotation(a) * rotation(b), rotation(a) * b.o + a.o};
}
template <typename T, size_t N>
inline frame<T, N>& operator*=(frame<T, N>& a, const frame<T, N>& b) {
  return a = a * b;
}

// Frame inverse, equivalent to rigid affine inverse.
template <typename T, size_t N>
inline frame<T, N> inverse(const frame<T, N>& a, bool non_rigid) {
  if (non_rigid) {
    auto minv = inverse(rotation(a));
    return {minv, -(minv * a.o)};
  } else {
    auto minv = transpose(rotation(a));
    return {minv, -(minv * a.o)};
  }
}

// Frame construction from axis.
template <typename T>
inline frame<T, 3> frame_fromz(const vec<T, 3>& o, const vec<T, 3>& v) {
  // https://graphics.pixar.com/library/OrthonormalB/paper.pdf
  auto z    = normalize(v);
  auto sign = copysignf((T)1, z.z);
  auto a    = -1 / (sign + z.z);
  auto b    = z.x * z.y * a;
  auto x    = vec<T, 3>{1 + sign * z.x * z.x * a, sign * b, -sign * z.x};
  auto y    = vec<T, 3>{b, sign + z.y * z.y * a, -z.y};
  return {x, y, z, o};
}
template <typename T>
inline frame<T, 3> frame_fromzx(
    const vec<T, 3>& o, const vec<T, 3>& z_, const vec<T, 3>& x_) {
  auto z = normalize(z_);
  auto x = orthonormalize(x_, z);
  auto y = normalize(cross(z, x));
  return {x, y, z, o};
}

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// QUATERNIONS
// -----------------------------------------------------------------------------
namespace yocto::math {

// Quaternions to represent rotations
template <typename T>
inline quat<T, 4>::quat() : x{0}, y{0}, z{0}, w{1} {}
template <typename T>
inline quat<T, 4>::quat(T x, T y, T z, T w) : x{x}, y{y}, z{z}, w{w} {}

// Quaternion operatons
template <typename T>
inline quat<T, 4> operator+(const quat<T, 4>& a, const quat<T, 4>& b) {
  return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}
template <typename T, typename T1>
inline quat<T, 4> operator*(const quat<T, 4>& a, T1 b) {
  return {a.x * b, a.y * b, a.z * b, a.w * b};
}
template <typename T, typename T1>
inline quat<T, 4> operator/(const quat<T, 4>& a, T1 b) {
  return {a.x / b, a.y / b, a.z / b, a.w / b};
}
template <typename T>
inline quat<T, 4> operator*(const quat<T, 4>& a, const quat<T, 4>& b) {
  return {a.x * b.w + a.w * b.x + a.y * b.w - a.z * b.y,
      a.y * b.w + a.w * b.y + a.z * b.x - a.x * b.z,
      a.z * b.w + a.w * b.z + a.x * b.y - a.y * b.x,
      a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z};
}

// Quaterion operations
template <typename T>
inline T dot(const quat<T, 4>& a, const quat<T, 4>& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}
template <typename T>
inline T length(const quat<T, 4>& a) {
  return sqrt(dot(a, a));
}
template <typename T>
inline quat<T, 4> normalize(const quat<T, 4>& a) {
  auto l = length(a);
  return (l != 0) ? a / l : a;
}
template <typename T>
inline quat<T, 4> conjugate(const quat<T, 4>& a) {
  return {-a.x, -a.y, -a.z, a.w};
}
template <typename T>
inline quat<T, 4> inverse(const quat<T, 4>& a) {
  return conjugate(a) / dot(a, a);
}
template <typename T>
inline float uangle(const quat<T, 4>& a, const quat<T, 4>& b) {
  auto d = dot(a, b);
  return d > 1 ? 0 : acos(d < -1 ? -1 : d);
}
template <typename T, typename T1>
inline quat<T, 4> lerp(const quat<T, 4>& a, const quat<T, 4>& b, T1 t) {
  return a * (1 - t) + b * t;
}
template <typename T, typename T1>
inline quat<T, 4> nlerp(const quat<T, 4>& a, const quat<T, 4>& b, T1 t) {
  return normalize(lerp(a, b, t));
}
template <typename T, typename T1>
inline quat<T, 4> slerp(const quat<T, 4>& a, const quat<T, 4>& b, T1 t) {
  auto th = uangle(a, b);
  return th == 0
             ? a
             : a * (sin(th * (1 - t)) / sin(th)) + b * (sin(th * t) / sin(th));
}

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// AXIS ALIGNED BOUNDING BOXES
// -----------------------------------------------------------------------------
namespace yocto::math {

// Axis aligned bounding box represented as a min/max std::vector pairs.
template <typename T, size_t N>
inline bbox<T, N>::bbox() {}
template <typename T, size_t N>
inline bbox<T, N>::bbox(const vec<T, N>& min, const vec<T, N>& max)
    : min{min}, max{max} {}

template <typename T, size_t N>
inline vec<T, N>& bbox<T, N>::operator[](int i) {
  return (&min)[i];
}
template <typename T, size_t N>
inline const vec<T, N>& bbox<T, N>::operator[](int i) const {
  return (&min)[i];
}

// Bounding box properties
template <typename T, size_t N>
inline vec<T, N> center(const bbox<T, N>& a) {
  return (a.min + a.max) / 2;
}
template <typename T, size_t N>
inline vec<T, N> size(const bbox<T, N>& a) {
  return a.max - a.min;
}

// Bounding box comparisons.
template <typename T, size_t N>
inline bool operator==(const bbox<T, N>& a, const bbox<T, N>& b) {
  return a.min == b.min && a.max == b.max;
}
template <typename T, size_t N>
inline bool operator!=(const bbox<T, N>& a, const bbox<T, N>& b) {
  return a.min != b.min || a.max != b.max;
}

// Bounding box expansions with points and other boxes.
template <typename T, size_t N>
inline bbox<T, N> merge(const bbox<T, N>& a, const vec<T, N>& b) {
  return {min(a.min, b), max(a.max, b)};
}
template <typename T, size_t N>
inline bbox<T, N> merge(const bbox<T, N>& a, const bbox<T, N>& b) {
  return {min(a.min, b.min), max(a.max, b.max)};
}
template <typename T, size_t N>
inline void expand(bbox<T, N>& a, const vec<T, N>& b) {
  a = merge(a, b);
}
template <typename T, size_t N>
inline void expand(bbox<T, N>& a, const bbox<T, N>& b) {
  a = merge(a, b);
}

// Primitive bounds.
template <typename T, size_t N>
inline bbox<T, N> point_bounds(const vec<T, N>& p) {
  return {p, p};
}
template <typename T, size_t N, typename T1>
inline bbox<T, N> point_bounds(const vec<T, N>& p, T1 r) {
  return {min(p - r, p + r), max(p - r, p + r)};
}
template <typename T, size_t N>
inline bbox<T, N> line_bounds(const vec<T, N>& p0, const vec<T, N>& p1) {
  return {min(p0, p1), max(p0, p1)};
}
template <typename T, size_t N, typename T1, typename T2>
inline bbox<T, N> line_bounds(
    const vec<T, N>& p0, const vec<T, N>& p1, T1 r0, T2 r1) {
  return {min(p0 - r0, p1 - r1), max(p0 + r0, p1 + r1)};
}
template <typename T, size_t N>
inline bbox<T, N> triangle_bounds(
    const vec<T, N>& p0, const vec<T, N>& p1, const vec<T, N>& p2) {
  return {min(p0, min(p1, p2)), max(p0, max(p1, p2))};
}
template <typename T, size_t N>
inline bbox<T, N> quad_bounds(const vec<T, N>& p0, const vec<T, N>& p1,
    const vec<T, N>& p2, const vec<T, N>& p3) {
  return {min(p0, min(p1, min(p2, p3))), max(p0, max(p1, max(p2, p3)))};
}

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// RAYS
// -----------------------------------------------------------------------------
namespace yocto::math {

// Rays with origin, direction and min/max t value.
template <typename T>
inline ray<T, 2>::ray() {}
template <typename T>
inline ray<T, 2>::ray(const vec<T, 2>& o, const vec<T, 2>& d, T tmin, T tmax)
    : o{o}, d{d}, tmin{tmin}, tmax{tmax} {}

// Rays with origin, direction and min/max t value.
template <typename T>
inline ray<T, 3>::ray() {}
template <typename T>
inline ray<T, 3>::ray(const vec<T, 3>& o, const vec<T, 3>& d, T tmin, T tmax)
    : o{o}, d{d}, tmin{tmin}, tmax{tmax} {}

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// TRANSFORMS
// -----------------------------------------------------------------------------
namespace yocto::math {

// Transforms points, vectors and directions by matrices.
template <typename T>
inline vec<T, 2> transform_point(const mat<T, 3, 3>& a, const vec<T, 2>& b) {
  auto tvb = a * vec<T, 3>{b.x, b.y, 1};
  return vec<T, 2>{tvb.x, tvb.y} / tvb.z;
}
template <typename T>
inline vec<T, 2> transform_vector(const mat<T, 3, 3>& a, const vec<T, 2>& b) {
  auto tvb = a * vec<T, 3>{b.x, b.y, 0};
  return vec<T, 2>{tvb.x, tvb.y} / tvb.z;
}
template <typename T>
inline vec<T, 2> transform_direction(
    const mat<T, 3, 3>& a, const vec<T, 2>& b) {
  return normalize(transform_vector(a, b));
}
template <typename T>
inline vec<T, 2> transform_normal(const mat<T, 3, 3>& a, const vec<T, 2>& b) {
  return normalize(transform_vector(transpose(inverse(a)), b));
}
template <typename T>
inline vec<T, 2> transform_vector(const mat<T, 2, 2>& a, const vec<T, 2>& b) {
  return a * b;
}
template <typename T>
inline vec<T, 2> transform_direction(
    const mat<T, 2, 2>& a, const vec<T, 2>& b) {
  return normalize(transform_vector(a, b));
}
template <typename T>
inline vec<T, 2> transform_normal(const mat<T, 2, 2>& a, const vec<T, 2>& b) {
  return normalize(transform_vector(transpose(inverse(a)), b));
}

template <typename T>
inline vec<T, 3> transform_point(const mat<T, 4, 4>& a, const vec<T, 3>& b) {
  auto tvb = a * vec<T, 4>{b.x, b.y, b.z, 1};
  return vec<T, 3>{tvb.x, tvb.y, tvb.z} / tvb.w;
}
template <typename T>
inline vec<T, 3> transform_vector(const mat<T, 4, 4>& a, const vec<T, 3>& b) {
  auto tvb = a * vec<T, 4>{b.x, b.y, b.z, 0};
  return vec<T, 3>{tvb.x, tvb.y, tvb.z};
}
template <typename T>
inline vec<T, 3> transform_direction(
    const mat<T, 4, 4>& a, const vec<T, 3>& b) {
  return normalize(transform_vector(a, b));
}
template <typename T>
inline vec<T, 3> transform_vector(const mat<T, 3, 3>& a, const vec<T, 3>& b) {
  return a * b;
}
template <typename T>
inline vec<T, 3> transform_direction(
    const mat<T, 3, 3>& a, const vec<T, 3>& b) {
  return normalize(transform_vector(a, b));
}
template <typename T>
inline vec<T, 3> transform_normal(const mat<T, 3, 3>& a, const vec<T, 3>& b) {
  return normalize(transform_vector(transpose(inverse(a)), b));
}

// Transforms points, vectors and directions by frames.
template <typename T>
inline vec<T, 2> transform_point(const frame<T, 2>& a, const vec<T, 2>& b) {
  return a.x * b.x + a.y * b.y + a.o;
}
template <typename T>
inline vec<T, 2> transform_vector(const frame<T, 2>& a, const vec<T, 2>& b) {
  return a.x * b.x + a.y * b.y;
}
template <typename T>
inline vec<T, 2> transform_direction(const frame<T, 2>& a, const vec<T, 2>& b) {
  return normalize(transform_vector(a, b));
}
template <typename T>
inline vec<T, 2> transform_normal(
    const frame<T, 2>& a, const vec<T, 2>& b, bool non_rigid) {
  if (non_rigid) {
    return transform_normal(rotation(a), b);
  } else {
    return normalize(transform_vector(a, b));
  }
}

// Transforms points, vectors and directions by frames.
template <typename T>
inline vec<T, 3> transform_point(const frame<T, 3>& a, const vec<T, 3>& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z + a.o;
}
template <typename T>
inline vec<T, 3> transform_vector(const frame<T, 3>& a, const vec<T, 3>& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}
template <typename T>
inline vec<T, 3> transform_direction(const frame<T, 3>& a, const vec<T, 3>& b) {
  return normalize(transform_vector(a, b));
}
template <typename T>
inline vec<T, 3> transform_normal(
    const frame<T, 3>& a, const vec<T, 3>& b, bool non_rigid) {
  if (non_rigid) {
    return transform_normal(rotation(a), b);
  } else {
    return normalize(transform_vector(a, b));
  }
}

// Transforms rays and bounding boxes by matrices.
template <typename T>
inline ray<T, 3> transform_ray(const mat<T, 4, 4>& a, const ray<T, 3>& b) {
  return {transform_point(a, b.o), transform_vector(a, b.d), b.tmin, b.tmax};
}
template <typename T>
inline ray<T, 3> transform_ray(const frame<T, 3>& a, const ray<T, 3>& b) {
  return {transform_point(a, b.o), transform_vector(a, b.d), b.tmin, b.tmax};
}
template <typename T>
inline bbox<T, 3> transform_bbox(const mat<T, 4, 4>& a, const bbox<T, 3>& b) {
  auto corners = {vec<T, 3>{b.min.x, b.min.y, b.min.z},
      vec<T, 3>{b.min.x, b.min.y, b.max.z},
      vec<T, 3>{b.min.x, b.max.y, b.min.z},
      vec<T, 3>{b.min.x, b.max.y, b.max.z},
      vec<T, 3>{b.max.x, b.min.y, b.min.z},
      vec<T, 3>{b.max.x, b.min.y, b.max.z},
      vec<T, 3>{b.max.x, b.max.y, b.min.z},
      vec<T, 3>{b.max.x, b.max.y, b.max.z}};
  auto xformed = bbox<T, 3>();
  for (auto& corner : corners)
    xformed = merge(xformed, transform_point(a, corner));
  return xformed;
}
template <typename T>
inline bbox<T, 3> transform_bbox(const frame<T, 3>& a, const bbox<T, 3>& b) {
  auto corners = {vec<T, 3>{b.min.x, b.min.y, b.min.z},
      vec<T, 3>{b.min.x, b.min.y, b.max.z},
      vec<T, 3>{b.min.x, b.max.y, b.min.z},
      vec<T, 3>{b.min.x, b.max.y, b.max.z},
      vec<T, 3>{b.max.x, b.min.y, b.min.z},
      vec<T, 3>{b.max.x, b.min.y, b.max.z},
      vec<T, 3>{b.max.x, b.max.y, b.min.z},
      vec<T, 3>{b.max.x, b.max.y, b.max.z}};
  auto xformed = bbox<T, 3>();
  for (auto& corner : corners)
    xformed = merge(xformed, transform_point(a, corner));
  return xformed;
}

// Translation, scaling and rotations transforms.
template <typename T>
inline frame<T, 3> translation_frame(const vec<T, 3>& a) {
  return {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, a};
}
template <typename T>
inline frame<T, 3> scaling_frame(const vec<T, 3>& a) {
  return {{a.x, 0, 0}, {0, a.y, 0}, {0, 0, a.z}, {0, 0, 0}};
}
template <typename T>
inline frame<T, 3> rotation_frame(const vec<T, 3>& axis, T angle) {
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
inline frame<T, 3> rotation_frame(const vec<T, 4>& quat) {
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
inline frame<T, 3> rotation_frame(const quat<T, 4>& quat) {
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
inline frame<T, 3> rotation_frame(const mat<T, 3, 3>& rot) {
  return {rot.x, rot.y, rot.z, {0, 0, 0}};
}

// Lookat frame. Z-axis can be inverted with inv_xz.
template <typename T>
inline frame<T, 3> lookat_frame(const vec<T, 3>& eye, const vec<T, 3>& center,
    const vec<T, 3>& up, bool inv_xz) {
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
inline mat<T, 4, 4> frustum_mat(T l, T r, T b, T t, T n, T f) {
  return {{2 * n / (r - l), 0, 0, 0}, {0, 2 * n / (t - b), 0, 0},
      {(r + l) / (r - l), (t + b) / (t - b), -(f + n) / (f - n), -1},
      {0, 0, -2 * f * n / (f - n), 0}};
}
template <typename T>
inline mat<T, 4, 4> ortho_mat(T l, T r, T b, T t, T n, T f) {
  return {{2 / (r - l), 0, 0, 0}, {0, 2 / (t - b), 0, 0},
      {0, 0, -2 / (f - n), 0},
      {-(r + l) / (r - l), -(t + b) / (t - b), -(f + n) / (f - n), 1}};
}
template <typename T>
inline mat<T, 4, 4> ortho2d_mat(T left, T right, T bottom, T top) {
  return ortho_mat(left, right, bottom, top, -1, 1);
}
template <typename T>
inline mat<T, 4, 4> ortho_mat(T xmag, T ymag, T near, T far) {
  return {{1 / xmag, 0, 0, 0}, {0, 1 / ymag, 0, 0}, {0, 0, 2 / (near - far), 0},
      {0, 0, (far + near) / (near - far), 1}};
}
template <typename T>
inline mat<T, 4, 4> perspective_mat(T fovy, T aspect, T near, T far) {
  auto tg = tan(fovy / 2);
  return {{1 / (aspect * tg), 0, 0, 0}, {0, 1 / tg, 0, 0},
      {0, 0, (far + near) / (near - far), -1},
      {0, 0, 2 * far * near / (near - far), 0}};
}
template <typename T>
inline mat<T, 4, 4> perspective_mat(T fovy, T aspect, T near) {
  auto tg = tan(fovy / 2);
  return {{1 / (aspect * tg), 0, 0, 0}, {0, 1 / tg, 0, 0}, {0, 0, -1, -1},
      {0, 0, 2 * near, 0}};
}

// Rotation conversions.
template <typename T>
inline std::pair<vec<T, 3>, T> rotation_axisangle(const vec<T, 4>& quat) {
  return {normalize(vec<T, 3>{quat.x, quat.y, quat.z}), 2 * acos(quat.w)};
}
template <typename T, typename T1>
inline vec<T, 4> rotation_quat(const vec<T, 3>& axis, T1 angle) {
  auto len = length(axis);
  if (!len) return {0, 0, 0, 1};
  return vec<T, 4>{sin(angle / 2) * axis.x / len, sin(angle / 2) * axis.y / len,
      sin(angle / 2) * axis.z / len, cos(angle / 2)};
}
template <typename T>
inline vec<T, 4> rotation_quat(const vec<T, 4>& axisangle) {
  return rotation_quat(
      vec<T, 3>{axisangle.x, axisangle.y, axisangle.z}, axisangle.w);
}

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// GEOMETRY UTILITIES
// -----------------------------------------------------------------------------
namespace yocto::math {

// Line properties.
template <typename T>
inline vec<T, 3> line_tangent(const vec<T, 3>& p0, const vec<T, 3>& p1) {
  return normalize(p1 - p0);
}
template <typename T>
inline T line_length(const vec<T, 3>& p0, const vec<T, 3>& p1) {
  return length(p1 - p0);
}

// Triangle properties.
template <typename T>
inline vec<T, 3> triangle_normal(
    const vec<T, 3>& p0, const vec<T, 3>& p1, const vec<T, 3>& p2) {
  return normalize(cross(p1 - p0, p2 - p0));
}
template <typename T>
inline T triangle_area(
    const vec<T, 3>& p0, const vec<T, 3>& p1, const vec<T, 3>& p2) {
  return length(cross(p1 - p0, p2 - p0)) / 2;
}

// Quad propeties.
template <typename T>
inline vec<T, 3> quad_normal(const vec<T, 3>& p0, const vec<T, 3>& p1,
    const vec<T, 3>& p2, const vec<T, 3>& p3) {
  return normalize(triangle_normal(p0, p1, p3) + triangle_normal(p2, p3, p1));
}
template <typename T>
inline T quad_area(const vec<T, 3>& p0, const vec<T, 3>& p1,
    const vec<T, 3>& p2, const vec<T, 3>& p3) {
  return triangle_area(p0, p1, p3) + triangle_area(p2, p3, p1);
}

// Interpolates values over a line parameterized from a to b by u. Same as lerp.
template <typename T, typename T1>
inline T interpolate_line(const T& p0, const T& p1, T1 u) {
  return p0 * (1 - u) + p1 * u;
}
// Interpolates values over a triangle parameterized by u and v along the
// (p1-p0) and (p2-p0) directions. Same as barycentric interpolation.
template <typename T, typename T1>
inline T interpolate_triangle(
    const T& p0, const T& p1, const T& p2, const vec<T1, 2>& uv) {
  return p0 * (1 - uv.x - uv.y) + p1 * uv.x + p2 * uv.y;
}
// Interpolates values over a quad parameterized by u and v along the
// (p1-p0) and (p2-p1) directions. Same as bilinear interpolation.
template <typename T, typename T1>
inline T interpolate_quad(
    const T& p0, const T& p1, const T& p2, const T& p3, const vec<T1, 2>& uv) {
  if (uv.x + uv.y <= 1) {
    return interpolate_triangle(p0, p1, p3, uv);
  } else {
    return interpolate_triangle(p2, p3, p1, 1 - uv);
  }
}

// Interpolates values along a cubic Bezier segment parametrized by u.
template <typename T, typename T1>
inline T interpolate_bezier(
    const T& p0, const T& p1, const T& p2, const T& p3, T1 u) {
  return p0 * (1 - u) * (1 - u) * (1 - u) + p1 * 3 * u * (1 - u) * (1 - u) +
         p2 * 3 * u * u * (1 - u) + p3 * u * u * u;
}
// Computes the derivative of a cubic Bezier segment parametrized by u.
template <typename T, typename T1>
inline T interpolate_bezier_derivative(
    const T& p0, const T& p1, const T& p2, const T& p3, T1 u) {
  return (p1 - p0) * 3 * (1 - u) * (1 - u) + (p2 - p1) * 6 * u * (1 - u) +
         (p3 - p2) * 3 * u * u;
}

// Triangle tangent and bitangent from uv
template <typename T>
inline std::pair<vec<T, 3>, vec<T, 3>> triangle_tangents_fromuv(
    const vec<T, 3>& p0, const vec<T, 3>& p1, const vec<T, 3>& p2,
    const vec<T, 2>& uv0, const vec<T, 2>& uv1, const vec<T, 2>& uv2) {
  // Follows the definition in http://www.terathon.com/code/tangent.html and
  // https://gist.github.com/aras-p/2843984
  // normal points up from texture space
  auto p   = p1 - p0;
  auto q   = p2 - p0;
  auto s   = vec<T, 2>{uv1.x - uv0.x, uv2.x - uv0.x};
  auto t   = vec<T, 2>{uv1.y - uv0.y, uv2.y - uv0.y};
  auto div = s.x * t.y - s.y * t.x;

  if (div != 0) {
    auto tu = vec<T, 3>{t.y * p.x - t.x * q.x, t.y * p.y - t.x * q.y,
                  t.y * p.z - t.x * q.z} /
              div;
    auto tv = vec<T, 3>{s.x * q.x - s.y * p.x, s.x * q.y - s.y * p.y,
                  s.x * q.z - s.y * p.z} /
              div;
    return {tu, tv};
  } else {
    return {{1, 0, 0}, {0, 1, 0}};
  }
}

// Quad tangent and bitangent from uv.
template <typename T>
inline std::pair<vec<T, 3>, vec<T, 3>> quad_tangents_fromuv(const vec<T, 3>& p0,
    const vec<T, 3>& p1, const vec<T, 3>& p2, const vec<T, 3>& p3,
    const vec<T, 2>& uv0, const vec<T, 2>& uv1, const vec<T, 2>& uv2,
    const vec<T, 2>& uv3, const vec<T, 2>& current_uv) {
  if (current_uv.x + current_uv.y <= 1) {
    return triangle_tangents_fromuv(p0, p1, p3, uv0, uv1, uv3);
  } else {
    return triangle_tangents_fromuv(p2, p3, p1, uv2, uv3, uv1);
  }
}

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// IMPLEMENRTATION OF RAY-PRIMITIVE INTERSECTION FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto::math {

// Intersect a ray with a point (approximate)
template <typename T>
inline bool intersect_point(
    const ray<T, 3>& ray, const vec<T, 3>& p, T r, vec<T, 2>& uv, T& dist) {
  // find parameter for line-point minimum distance
  auto w = p - ray.o;
  auto t = dot(w, ray.d) / dot(ray.d, ray.d);

  // exit if not within bounds
  if (t < ray.tmin || t > ray.tmax) return false;

  // test for line-point distance vs point radius
  auto rp  = ray.o + ray.d * t;
  auto prp = p - rp;
  if (dot(prp, prp) > r * r) return false;

  // intersection occurred: set params and exit
  uv   = {0, 0};
  dist = t;
  return true;
}

// Intersect a ray with a line
template <typename T>
inline bool intersect_line(const ray<T, 3>& ray, const vec<T, 3>& p0,
    const vec<T, 3>& p1, T r0, T r1, vec<T, 2>& uv, T& dist) {
  // setup intersection params
  auto u = ray.d;
  auto v = p1 - p0;
  auto w = ray.o - p0;

  // compute values to solve a linear system
  auto a   = dot(u, u);
  auto b   = dot(u, v);
  auto c   = dot(v, v);
  auto d   = dot(u, w);
  auto e   = dot(v, w);
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
  s = clamp(s, (T)0, (T)1);

  // compute segment-segment distance on the closest points
  auto pr  = ray.o + ray.d * t;
  auto pl  = p0 + (p1 - p0) * s;
  auto prl = pr - pl;

  // check with the line radius at the same point
  auto d2 = dot(prl, prl);
  auto r  = r0 * (1 - s) + r1 * s;
  if (d2 > r * r) return {};

  // intersection occurred: set params and exit
  uv   = {s, sqrt(d2) / r};
  dist = t;
  return true;
}

// Intersect a ray with a triangle
template <typename T>
inline bool intersect_triangle(const ray<T, 3>& ray, const vec<T, 3>& p0,
    const vec<T, 3>& p1, const vec<T, 3>& p2, vec<T, 2>& uv, T& dist) {
  // compute triangle edges
  auto edge1 = p1 - p0;
  auto edge2 = p2 - p0;

  // compute determinant to solve a linear system
  auto pvec = cross(ray.d, edge2);
  auto det  = dot(edge1, pvec);

  // check determinant and exit if triangle and ray are parallel
  // (could use EPSILONS if desired)
  if (det == 0) return false;
  auto inv_det = 1 / det;

  // compute and check first bricentric coordinated
  auto tvec = ray.o - p0;
  auto u    = dot(tvec, pvec) * inv_det;
  if (u < 0 || u > 1) return false;

  // compute and check second bricentric coordinated
  auto qvec = cross(tvec, edge1);
  auto v    = dot(ray.d, qvec) * inv_det;
  if (v < 0 || u + v > 1) return false;

  // compute and check ray parameter
  auto t = dot(edge2, qvec) * inv_det;
  if (t < ray.tmin || t > ray.tmax) return false;

  // intersection occurred: set params and exit
  uv   = {u, v};
  dist = t;
  return true;
}

// Intersect a ray with a quad.
template <typename T>
inline bool intersect_quad(const ray<T, 3>& ray, const vec<T, 3>& p0,
    const vec<T, 3>& p1, const vec<T, 3>& p2, const vec<T, 3>& p3,
    vec<T, 2>& uv, T& dist) {
  if (p2 == p3) {
    return intersect_triangle(ray, p0, p1, p3, uv, dist);
  }
  auto hit  = false;
  auto tray = ray;
  if (intersect_triangle(tray, p0, p1, p3, uv, dist)) {
    hit       = true;
    tray.tmax = dist;
  }
  if (intersect_triangle(tray, p2, p3, p1, uv, dist)) {
    hit       = true;
    uv        = 1 - uv;
    tray.tmax = dist;
  }
  return hit;
}

// Intersect a ray with a axis-aligned bounding box
template <typename T>
inline bool intersect_bbox(const ray<T, 3>& ray, const bbox<T, 3>& bbox) {
  // determine intersection ranges
  auto invd = 1 / ray.d;
  auto t0   = (bbox.min - ray.o) * invd;
  auto t1   = (bbox.max - ray.o) * invd;
  // flip based on range directions
  if (invd.x < 0.0f) swap(t0.x, t1.x);
  if (invd.y < 0.0f) swap(t0.y, t1.y);
  if (invd.z < 0.0f) swap(t0.z, t1.z);
  auto tmin = max(t0.z, max(t0.y, max(t0.x, ray.tmin)));
  auto tmax = min(t1.z, min(t1.y, min(t1.x, ray.tmax)));
  if constexpr (std::is_same_v<T, float>) {
    tmax *= 1.00000024f;  // for double: 1.0000000000000004
  }
  if constexpr (std::is_same_v<T, double>) {
    tmax *= 1.0000000000000004;
  }
  return tmin <= tmax;
}

// Intersect a ray with a axis-aligned bounding box
template <typename T>
inline bool intersect_bbox(
    const ray<T, 3>& ray, const vec<T, 3>& ray_dinv, const bbox<T, 3>& bbox) {
  auto it_min = (bbox.min - ray.o) * ray_dinv;
  auto it_max = (bbox.max - ray.o) * ray_dinv;
  auto tmin   = min(it_min, it_max);
  auto tmax   = max(it_min, it_max);
  auto t0     = max(max(tmin), ray.tmin);
  auto t1     = min(min(tmax), ray.tmax);
  if constexpr (std::is_same_v<T, float>) {
    t1 *= 1.00000024f;  // for double: 1.0000000000000004
  }
  if constexpr (std::is_same_v<T, double>) {
    t1 *= 1.0000000000000004;
  }
  return t0 <= t1;
}

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF POINT-PRIMITIVE DISTANCE FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto::math {

// Check if a point overlaps a position pos withint a maximum distance dist_max.
template <typename T>
inline bool overlap_point(const vec<T, 3>& pos, T dist_max, const vec<T, 3>& p,
    T r, vec<T, 2>& uv, T& dist) {
  auto d2 = dot(pos - p, pos - p);
  if (d2 > (dist_max + r) * (dist_max + r)) return false;
  uv   = {0, 0};
  dist = sqrt(d2);
  return true;
}

// Compute the closest line uv to a give position pos.
template <typename T>
inline T closestuv_line(
    const vec<T, 3>& pos, const vec<T, 3>& p0, const vec<T, 3>& p1) {
  auto ab = p1 - p0;
  auto d  = dot(ab, ab);
  // Project c onto ab, computing parameterized position d(t) = a + t*(b –
  // a)
  auto u = dot(pos - p0, ab) / d;
  u      = clamp(u, (float)0, (float)1);
  return u;
}

// Check if a line overlaps a position pos withint a maximum distance dist_max.
template <typename T>
inline bool overlap_line(const vec<T, 3>& pos, T dist_max, const vec<T, 3>& p0,
    const vec<T, 3>& p1, T r0, T r1, vec<T, 2>& uv, T& dist) {
  auto u = closestuv_line(pos, p0, p1);
  // Compute projected position from the clamped t d = a + t * ab;
  auto p  = p0 + (p1 - p0) * u;
  auto r  = r0 + (r1 - r0) * u;
  auto d2 = dot(pos - p, pos - p);
  // check distance
  if (d2 > (dist_max + r) * (dist_max + r)) return false;
  // done
  uv   = {u, 0};
  dist = sqrt(d2);
  return true;
}

// Compute the closest triangle uv to a give position pos.
template <typename T>
inline vec<T, 2> closestuv_triangle(const vec<T, 3>& pos, const vec<T, 3>& p0,
    const vec<T, 3>& p1, const vec<T, 3>& p2) {
  // this is a complicated test -> I probably "--"+prefix to use a sequence of
  // test (triangle body, and 3 edges)
  auto ab = p1 - p0;
  auto ac = p2 - p0;
  auto ap = pos - p0;

  auto d1 = dot(ab, ap);
  auto d2 = dot(ac, ap);

  // corner and edge cases
  if (d1 <= 0 && d2 <= 0) return {0, 0};

  auto bp = pos - p1;
  auto d3 = dot(ab, bp);
  auto d4 = dot(ac, bp);
  if (d3 >= 0 && d4 <= d3) return {1, 0};

  auto vc = d1 * d4 - d3 * d2;
  if ((vc <= 0) && (d1 >= 0) && (d3 <= 0)) return {d1 / (d1 - d3), 0};

  auto cp = pos - p2;
  auto d5 = dot(ab, cp);
  auto d6 = dot(ac, cp);
  if (d6 >= 0 && d5 <= d6) return {0, 1};

  auto vb = d5 * d2 - d1 * d6;
  if ((vb <= 0) && (d2 >= 0) && (d6 <= 0)) return {0, d2 / (d2 - d6)};

  auto va = d3 * d6 - d5 * d4;
  if ((va <= 0) && (d4 - d3 >= 0) && (d5 - d6 >= 0)) {
    auto w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
    return {1 - w, w};
  }

  // face case
  auto denom = 1 / (va + vb + vc);
  auto u     = vb * denom;
  auto v     = vc * denom;
  return {u, v};
}

// Check if a triangle overlaps a position pos withint a maximum distance
// dist_max.
template <typename T>
inline bool overlap_triangle(const vec<T, 3>& pos, T dist_max,
    const vec<T, 3>& p0, const vec<T, 3>& p1, const vec<T, 3>& p2, T r0, T r1,
    T r2, vec<T, 2>& uv, T& dist) {
  auto cuv = closestuv_triangle(pos, p0, p1, p2);
  auto p   = p0 * (1 - cuv.x - cuv.y) + p1 * cuv.x + p2 * cuv.y;
  auto r   = r0 * (1 - cuv.x - cuv.y) + r1 * cuv.x + r2 * cuv.y;
  auto dd  = dot(p - pos, p - pos);
  if (dd > (dist_max + r) * (dist_max + r)) return false;
  uv   = cuv;
  dist = sqrt(dd);
  return true;
}

// Check if a quad overlaps a position pos withint a maximum distance dist_max.
template <typename T>
inline bool overlap_quad(const vec<T, 3>& pos, T dist_max, const vec<T, 3>& p0,
    const vec<T, 3>& p1, const vec<T, 3>& p2, const vec<T, 3>& p3, T r0, T r1,
    T r2, T r3, vec<T, 2>& uv, T& dist) {
  if (p2 == p3) {
    return overlap_triangle(pos, dist_max, p0, p1, p3, r0, r1, r2, uv, dist);
  }
  auto hit = false;
  if (overlap_triangle(pos, dist_max, p0, p1, p3, r0, r1, r2, uv, dist)) {
    hit      = true;
    dist_max = dist;
  }
  if (!overlap_triangle(pos, dist_max, p2, p3, p1, r2, r3, r1, uv, dist)) {
    hit = true;
    uv  = 1 - uv;
    // dist_max = dist;
  }
  return hit;
}

// Check if a bbox overlaps a position pos withint a maximum distance dist_max.
template <typename T>
inline bool distance_check_bbox(
    const vec<T, 3>& pos, T dist_max, const bbox<T, 3>& bbox) {
  // computing distance
  auto dd = 0.0f;

  // For each axis count any excess distance outside box extents
  if (pos.x < bbox.min.x) dd += (bbox.min.x - pos.x) * (bbox.min.x - pos.x);
  if (pos.x > bbox.max.x) dd += (pos.x - bbox.max.x) * (pos.x - bbox.max.x);
  if (pos.y < bbox.min.y) dd += (bbox.min.y - pos.y) * (bbox.min.y - pos.y);
  if (pos.y > bbox.max.y) dd += (pos.y - bbox.max.y) * (pos.y - bbox.max.y);
  if (pos.z < bbox.min.z) dd += (bbox.min.z - pos.z) * (bbox.min.z - pos.z);
  if (pos.z > bbox.max.z) dd += (pos.z - bbox.max.z) * (pos.z - bbox.max.z);

  // check distance
  return dd < dist_max * dist_max;
}

// Check if two bboxe overlap.
template <typename T>
inline bool overlap_bbox(const bbox<T, 3>& bbox1, const bbox<T, 3>& bbox2) {
  if (bbox1.max.x < bbox2.min.x || bbox1.min.x > bbox2.max.x) return false;
  if (bbox1.max.y < bbox2.min.y || bbox1.min.y > bbox2.max.y) return false;
  if (bbox1.max.z < bbox2.min.z || bbox1.min.z > bbox2.max.z) return false;
  return true;
}

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR COLOR CONVERSION UTILITIES
// -----------------------------------------------------------------------------
namespace yocto::math {

// Conversion between flots and bytes
inline vec3b float_to_byte(const vec3f& a) {
  return {(byte)clamp(int(a.x * 256), 0, 255),
      (byte)clamp(int(a.y * 256), 0, 255), (byte)clamp(int(a.z * 256), 0, 255)};
}
inline vec3f byte_to_float(const vec3b& a) {
  return {a.x / 255.0f, a.y / 255.0f, a.z / 255.0f};
}
inline vec4b float_to_byte(const vec4f& a) {
  return {(byte)clamp(int(a.x * 256), 0, 255),
      (byte)clamp(int(a.y * 256), 0, 255), (byte)clamp(int(a.z * 256), 0, 255),
      (byte)clamp(int(a.w * 256), 0, 255)};
}
inline vec4f byte_to_float(const vec4b& a) {
  return {a.x / 255.0f, a.y / 255.0f, a.z / 255.0f, a.w / 255.0f};
}
inline byte float_to_byte(float a) { return (byte)clamp(int(a * 256), 0, 255); }
inline float  byte_to_float(byte a) { return a / 255.0f; }
inline ushort float_to_ushort(float a) {
  return (ushort)clamp(int(a * 65536), 0, 65535);
}
inline float ushort_to_float(ushort a) { return a / 65535.0f; }

// Conversion between reals in [0,1] and normalized ints [0, max_int]
template <typename I, typename T>
inline I real_to_nint(T a) {
  return clamp(I(a * ((T)type_max<I> + (T)1)), (I)0, type_max<I>);
}
template <typename T, typename I>
inline T nint_to_real(I a) {
  return a / (T)type_max<I>;
}
template <typename I, typename T, size_t N>
inline vec<I, N> real_to_nint(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return {real_to_nint<I, T>(a.x)};
  } else if constexpr (N == 2) {
    return {real_to_nint<I, T>(a.x), real_to_nint<I, T>(a.y)};
  } else if constexpr (N == 3) {
    return {real_to_nint<I, T>(a.x), real_to_nint<I, T>(a.y),
        real_to_nint<I, T>(a.z)};
  } else if constexpr (N == 4) {
    return {real_to_nint<I, T>(a.x), real_to_nint<I, T>(a.y),
        real_to_nint<I, T>(a.z), real_to_nint<I, T>(a.w)};
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}
template <typename T, typename I, size_t N>
inline vec<T, N> nint_to_real(const vec<I, N>& a) {
  if constexpr (N == 1) {
    return {nint_to_real<T, I>(a.x)};
  } else if constexpr (N == 2) {
    return {nint_to_real<T, I>(a.x), nint_to_real<T, I>(a.y)};
  } else if constexpr (N == 3) {
    return {nint_to_real<T, I>(a.x), nint_to_real<T, I>(a.y),
        nint_to_real<T, I>(a.z)};
  } else if constexpr (N == 4) {
    return {nint_to_real<T, I>(a.x), nint_to_real<T, I>(a.y),
        nint_to_real<T, I>(a.z), nint_to_real<T, I>(a.w)};
  } else {
    static_assert(N >= 0 || N <= 4, "vector size unsupported");
  }
}

// Luminance
template <typename T>
inline T luminance(const vec<T, 3>& a) {
  return (0.2126f * a.x + 0.7152f * a.y + 0.0722f * a.z);
}

// sRGB non-linear curve
template <typename T>
inline T srgb_to_rgb(T srgb) {
  return (srgb <= (T)0.04045)
             ? srgb / (T)12.92
             : pow((srgb + (T)0.055) / ((T)1.0 + (T)0.055), (T)2.4);
}
template <typename T>
inline T rgb_to_srgb(T rgb) {
  return (rgb <= (T)0.0031308)
             ? (T)12.92 * rgb
             : (1 + (T)0.055) * pow(rgb, 1 / (T)2.4) - (T)0.055;
}
template <typename T>
inline vec<T, 3> srgb_to_rgb(const vec<T, 3>& srgb) {
  return {srgb_to_rgb(srgb.x), srgb_to_rgb(srgb.y), srgb_to_rgb(srgb.z)};
}
template <typename T>
inline vec<T, 4> srgb_to_rgb(const vec<T, 4>& srgb) {
  return {
      srgb_to_rgb(srgb.x), srgb_to_rgb(srgb.y), srgb_to_rgb(srgb.z), srgb.w};
}
template <typename T>
inline vec<T, 3> rgb_to_srgb(const vec<T, 3>& rgb) {
  return {rgb_to_srgb(rgb.x), rgb_to_srgb(rgb.y), rgb_to_srgb(rgb.z)};
}
template <typename T>
inline vec<T, 4> rgb_to_srgb(const vec<T, 4>& rgb) {
  return {rgb_to_srgb(rgb.x), rgb_to_srgb(rgb.y), rgb_to_srgb(rgb.z), rgb.w};
}

// Apply contrast. Grey should be 0.18 for linear and 0.5 for gamma.
template <typename T>
inline vec<T, 3> lincontrast(const vec<T, 3>& rgb, T contrast, T grey) {
  return max(grey + (rgb - grey) * (contrast * 2), 0);
}
// Apply contrast in log2. Grey should be 0.18 for linear and 0.5 for gamma.
template <typename T>
inline vec<T, 3> logcontrast(const vec<T, 3>& rgb, T logcontrast, T grey) {
  auto epsilon  = (T)0.0001;
  auto log_grey = log2(grey);
  auto log_ldr  = log2(rgb + epsilon);
  auto adjusted = log_grey + (log_ldr - log_grey) * (logcontrast * 2);
  return max(exp2(adjusted) - epsilon, 0);
}
// Apply an s-shaped contrast.
template <typename T>
inline vec<T, 3> contrast(const vec<T, 3>& rgb, T contrast) {
  return gain(rgb, 1 - contrast);
}
// Apply saturation.
template <typename T>
inline vec<T, 3> saturate(
    const vec<T, 3>& rgb, T saturation, const vec<T, 3>& weights) {
  auto grey = dot(weights, rgb);
  return max(grey + (rgb - grey) * (saturation * 2), 0);
}

// Filmic tonemapping
template <typename T>
inline vec<T, 3> tonemap_filmic(
    const vec<T, 3>& hdr_, bool accurate_fit = false) {
  if (!accurate_fit) {
    // https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
    auto hdr = hdr_ * (T)0.6;  // brings it back to ACES range
    auto ldr = (hdr * hdr * (T)2.51 + hdr * (T)0.03) /
               (hdr * hdr * (T)2.43 + hdr * (T)0.59 + (T)0.14);
    return max(ldr, 0);
  } else {
    // https://github.com/TheRealMJP/BakingLab/blob/master/BakingLab/ACES.hlsl
    // sRGB => XYZ => D65_2_D60 => AP1 => RRT_SAT
    static const auto ACESInputMat = transpose(mat3f{
        {0.59719, 0.35458, 0.04823},
        {0.07600, 0.90834, 0.01566},
        {0.02840, 0.13383, 0.83777},
    });
    // ODT_SAT => XYZ => D60_2_D65 => sRGB
    static const auto ACESOutputMat = transpose(mat3f{
        {1.60475, -0.53108, -0.07367},
        {-0.10208, 1.10813, -0.00605},
        {-0.00327, -0.07276, 1.07602},
    });
    // RRT => ODT
    auto RRTAndODTFit = [](const vec<T, 3>& v) -> vec<T, 3> {
      return (v * v + v * (T)0.0245786 - (T)0.000090537) /
             (v * v * (T)0.983729 + v * (T)0.4329510 + (T)0.238081);
    };

    auto ldr = ACESOutputMat * RRTAndODTFit(ACESInputMat * hdr_);
    return max(ldr, 0);
  }
}

template <typename T>
inline vec<T, 3> tonemap(
    const vec<T, 3>& hdr, T exposure, bool filmic, bool srgb) {
  auto rgb = hdr;
  if (exposure != 0) rgb *= exp2(exposure);
  if (filmic) rgb = tonemap_filmic(rgb);
  if (srgb) rgb = rgb_to_srgb(rgb);
  return rgb;
}
template <typename T>
inline vec<T, 4> tonemap(
    const vec<T, 4>& hdr, T exposure, bool filmic, bool srgb) {
  return {tonemap(xyz(hdr), exposure, filmic, srgb), hdr.w};
}

// Convert between CIE XYZ and RGB
template <typename T>
inline vec<T, 3> rgb_to_xyz(const vec<T, 3>& rgb) {
  // https://en.wikipedia.org/wiki/SRGB
  static const auto mat_ = mat<T, 3, 3>{
      {0.4124, 0.2126, 0.0193},
      {0.3576, 0.7152, 0.1192},
      {0.1805, 0.0722, 0.9504},
  };
  return mat_ * rgb;
}
template <typename T>
inline vec<T, 3> xyz_to_rgb(const vec<T, 3>& xyz) {
  // https://en.wikipedia.org/wiki/SRGB
  static const auto mat_ = mat<T, 3, 3>{
      {+3.2406, -0.9689, +0.0557},
      {-1.5372, +1.8758, -0.2040},
      {-0.4986, +0.0415, +1.0570},
  };
  return mat_ * xyz;
}

// Convert between CIE XYZ and xyY
template <typename T>
inline vec<T, 3> xyz_to_xyY(const vec<T, 3>& xyz) {
  if (xyz == vec<T, 3>{0}) return {0};
  return {
      xyz.x / (xyz.x + xyz.y + xyz.z), xyz.y / (xyz.x + xyz.y + xyz.z), xyz.y};
}
template <typename T>
inline vec<T, 3> xyY_to_xyz(const vec<T, 3>& xyY) {
  if (xyY.y == 0) return vec<T, 3>{0};
  return {xyY.x * xyY.z / xyY.y, xyY.z, (1 - xyY.x - xyY.y) * xyY.z / xyY.y};
}

// Convert HSV to RGB
template <typename T>
inline vec<T, 3> hsv_to_rgb(const vec<T, 3>& hsv) {
  // from Imgui.cpp
  auto h = hsv.x, s = hsv.y, v = hsv.z;
  if (hsv.y == 0) return {v, v, v};

  h       = fmod(h, (T)1.0) / ((T)60.0 / (T)360.0);
  int   i = (int)h;
  float f = h - (float)i;
  float p = v * (1 - s);
  float q = v * (1 - s * f);
  float t = v * (1 - s * (1 - f));

  switch (i) {
    case 0: return {v, t, p};
    case 1: return {q, v, p};
    case 2: return {p, v, t};
    case 3: return {p, q, v};
    case 4: return {t, p, v};
    case 5: return {v, p, q};
    default: return {v, p, q};
  }
}

template <typename T>
inline vec<T, 3> rgb_to_hsv(const vec<T, 3>& rgb) {
  // from Imgui.cpp
  auto r = rgb.x, g = rgb.y, b = rgb.z;
  auto K = 0.f;
  if (g < b) {
    swap(g, b);
    K = -1;
  }
  if (r < g) {
    swap(r, g);
    K = -2 / 6.0f - K;
  }

  auto chroma = r - (g < b ? g : b);
  return {
      abs(K + (g - b) / (6 * chroma + (T)1e-20)), chroma / (r + (T)1e-20), r};
}

// Approximate color of blackbody radiation from wavelength in nm.
template <typename T>
inline vec<T, 3> blackbody_to_rgb(T temperature) {
  // https://github.com/neilbartlett/color-temperature
  auto rgb = vec<T, 3>{0};
  if ((temperature / 100) < 66) {
    rgb.x = 255;
  } else {
    // a + b x + c Log[x] /.
    // {a -> 351.97690566805693`,
    // b -> 0.114206453784165`,
    // c -> -40.25366309332127
    // x -> (kelvin/100) - 55}
    rgb.x = (temperature / 100) - 55;
    rgb.x = (T)351.97690566805693 + (T)0.114206453784165 * rgb.x -
            (T)40.25366309332127 * log(rgb.x);
    if (rgb.x < 0) rgb.x = 0;
    if (rgb.x > 255) rgb.x = 255;
  }

  if ((temperature / 100) < 66) {
    // a + b x + c Log[x] /.
    // {a -> -155.25485562709179`,
    // b -> -0.44596950469579133`,
    // c -> 104.49216199393888`,
    // x -> (kelvin/100) - 2}
    rgb.y = (temperature / 100) - 2;
    rgb.y = (T)-155.25485562709179 - (T)0.44596950469579133 * rgb.y +
            (T)104.49216199393888 * log(rgb.y);
    if (rgb.y < 0) rgb.y = 0;
    if (rgb.y > 255) rgb.y = 255;
  } else {
    // a + b x + c Log[x] /.
    // {a -> 325.4494125711974`,
    // b -> 0.07943456536662342`,
    // c -> -28.0852963507957`,
    // x -> (kelvin/100) - 50}
    rgb.y = (temperature / 100) - 50;
    rgb.y = (T)325.4494125711974 + (T)0.07943456536662342 * rgb.y -
            (T)28.0852963507957 * log(rgb.y);
    if (rgb.y < 0) rgb.y = 0;
    if (rgb.y > 255) rgb.y = 255;
  }

  if ((temperature / 100) >= 66) {
    rgb.z = 255;
  } else {
    if ((temperature / 100) <= 20) {
      rgb.z = 0;
    } else {
      // a + b x + c Log[x] /.
      // {a -> -254.76935184120902`,
      // b -> 0.8274096064007395`,
      // c -> 115.67994401066147`,
      // x -> kelvin/100 - 10}
      rgb.z = (temperature / 100) - 10;
      rgb.z = (T)-254.76935184120902 + (T)0.8274096064007395 * rgb.z +
              (T)115.67994401066147 * log(rgb.z);
      if (rgb.z < 0) rgb.z = 0;
      if (rgb.z > 255) rgb.z = 255;
    }
  }

  return srgb_to_rgb(rgb / 255);
}

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// IMPLEMENTATION RANDOM NUMBER GENERATION
// -----------------------------------------------------------------------------
namespace yocto::math {

// PCG random numbers from http://www.pcg-random.org/
inline rng_state::rng_state()
    : state{0x853c49e6748fea9bULL}, inc{0xda3e39cb94b95bdbULL} {}
inline rng_state::rng_state(uint64_t state, uint64_t inc)
    : state{state}, inc{inc} {}

// Next random number, used internally only.
inline uint32_t _advance_rng(rng_state& rng) {
  uint64_t oldstate   = rng.state;
  rng.state           = oldstate * 6364136223846793005ULL + rng.inc;
  uint32_t xorshifted = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
  uint32_t rot        = (uint32_t)(oldstate >> 59u);
  return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
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
inline vec2i rand2i(rng_state& rng, int n) {
  // force order of evaluation by using separate assignments.
  auto x = rand1i(rng, n);
  auto y = rand1i(rng, n);
  return {x, y};
}
inline vec3i rand3i(rng_state& rng, int n) {
  // force order of evaluation by using separate assignments.
  auto x = rand1i(rng, n);
  auto y = rand1i(rng, n);
  auto z = rand1i(rng, n);
  return {x, y, z};
}
inline vec4i rand4i(rng_state& rng, int n) {
  // force order of evaluation by using separate assignments.
  auto x = rand1i(rng, n);
  auto y = rand1i(rng, n);
  auto z = rand1i(rng, n);
  auto w = rand1i(rng, n);
  return {x, y, z, w};
}
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
inline vec4f rand4f(rng_state& rng) {
  // force order of evaluation by using separate assignments.
  auto x = rand1f(rng);
  auto y = rand1f(rng);
  auto z = rand1f(rng);
  auto w = rand1f(rng);
  return {x, y, z, w};
}

// Shuffles a sequence of elements
template <typename T>
inline void shuffle(std::vector<T>& vals, rng_state& rng) {
  // https://en.wikipedia.org/wiki/Fisher–Yates_shuffle
  for (auto i = (int)vals.size() - 1; i > 0; i--) {
    auto j = rand1i(rng, i + 1);
    std::swap(vals[j], vals[i]);
  }
}

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PERLIN NOISE
// -----------------------------------------------------------------------------
namespace yocto::math {

// clang-format off
inline float _stb__perlin_lerp(float a, float b, float t)
{
   return a + (b-a) * t;
}

inline int _stb__perlin_fastfloor(float a)
{
    int ai = (int) a;
    return (a < ai) ? ai-1 : ai;
}

// different grad function from Perlin's, but easy to modify to match reference
inline float _stb__perlin_grad(int hash, float x, float y, float z)
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

inline float _stb_perlin_noise3(float x, float y, float z, int x_wrap, int y_wrap, int z_wrap)
{
    // not same permutation table as Perlin's reference to avoid copyright issues;
    // Perlin's table can be found at http://mrl.nyu.edu/~perlin/noise/
    // @OPTIMIZE: should this be unsigned char instead of int for cache?
    static unsigned char _stb__perlin_randtab[512] =
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
   int px = _stb__perlin_fastfloor(x);
   int py = _stb__perlin_fastfloor(y);
   int pz = _stb__perlin_fastfloor(z);
   int x0 = px & x_mask, x1 = (px+1) & x_mask;
   int y0 = py & y_mask, y1 = (py+1) & y_mask;
   int z0 = pz & z_mask, z1 = (pz+1) & z_mask;
   int r0,r1, r00,r01,r10,r11;

   #define _stb__perlin_ease(a)   (((a*6-15)*a + 10) * a * a * a)

   x -= px; u = _stb__perlin_ease(x);
   y -= py; v = _stb__perlin_ease(y);
   z -= pz; w = _stb__perlin_ease(z);

   r0 = _stb__perlin_randtab[x0];
   r1 = _stb__perlin_randtab[x1];

   r00 = _stb__perlin_randtab[r0+y0];
   r01 = _stb__perlin_randtab[r0+y1];
   r10 = _stb__perlin_randtab[r1+y0];
   r11 = _stb__perlin_randtab[r1+y1];

   n000 = _stb__perlin_grad(_stb__perlin_randtab[r00+z0], x  , y  , z   );
   n001 = _stb__perlin_grad(_stb__perlin_randtab[r00+z1], x  , y  , z-1 );
   n010 = _stb__perlin_grad(_stb__perlin_randtab[r01+z0], x  , y-1, z   );
   n011 = _stb__perlin_grad(_stb__perlin_randtab[r01+z1], x  , y-1, z-1 );
   n100 = _stb__perlin_grad(_stb__perlin_randtab[r10+z0], x-1, y  , z   );
   n101 = _stb__perlin_grad(_stb__perlin_randtab[r10+z1], x-1, y  , z-1 );
   n110 = _stb__perlin_grad(_stb__perlin_randtab[r11+z0], x-1, y-1, z   );
   n111 = _stb__perlin_grad(_stb__perlin_randtab[r11+z1], x-1, y-1, z-1 );

   n00 = _stb__perlin_lerp(n000,n001,w);
   n01 = _stb__perlin_lerp(n010,n011,w);
   n10 = _stb__perlin_lerp(n100,n101,w);
   n11 = _stb__perlin_lerp(n110,n111,w);

   n0 = _stb__perlin_lerp(n00,n01,v);
   n1 = _stb__perlin_lerp(n10,n11,v);

   return _stb__perlin_lerp(n0,n1,u);
}

inline float _stb_perlin_ridge_noise3(float x, float y, float z,float lacunarity, float gain, float offset, int octaves,int x_wrap, int y_wrap, int z_wrap)
{
   int i;
   float frequency = 1.0f;
   float prev = 1.0f;
   float amplitude = 0.5f;
   float sum = 0.0f;

   for (i = 0; i < octaves; i++) {
      float r = (float)(_stb_perlin_noise3(x*frequency,y*frequency,z*frequency,x_wrap,y_wrap,z_wrap));
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

inline float _stb_perlin_fbm_noise3(float x, float y, float z,float lacunarity, float gain, int octaves,int x_wrap, int y_wrap, int z_wrap)
{
   int i;
   float frequency = 1.0f;
   float amplitude = 1.0f;
   float sum = 0.0f;

   for (i = 0; i < octaves; i++) {
      sum += _stb_perlin_noise3(x*frequency,y*frequency,z*frequency,x_wrap,y_wrap,z_wrap)*amplitude;
      frequency *= lacunarity;
      amplitude *= gain;
   }
   return sum;
}

inline float _stb_perlin_turbulence_noise3(float x, float y, float z, float lacunarity, float gain, int octaves,int x_wrap, int y_wrap, int z_wrap)
{
   int i;
   float frequency = 1.0f;
   float amplitude = 1.0f;
   float sum = 0.0f;

   for (i = 0; i < octaves; i++) {
      float r = _stb_perlin_noise3(x*frequency,y*frequency,z*frequency,x_wrap,y_wrap,z_wrap)*amplitude;
      r = r<0 ? -r : r; // fabs()
      sum += r;
      frequency *= lacunarity;
      amplitude *= gain;
   }
   return sum;
}
// clang-format on

// adapeted  stb_perlin.h
template <typename T>
inline T perlin_noise(const vec<T, 3>& p, const vec<int, 3>& wrap) {
  return _stb_perlin_noise3(
      (float)p.x, (float)p.y, (float)p.z, wrap.x, wrap.y, wrap.z);
}

// adapeted  stb_perlin.h
template <typename T>
inline T perlin_ridge(const vec<T, 3>& p, T lacunarity, T gain, int octaves,
    T offset, const vec<int, 3>& wrap) {
  return _stb_perlin_ridge_noise3((float)p.x, (float)p.y, (float)p.z,
      (float)lacunarity, (float)gain, offset, octaves, wrap.x, wrap.y, wrap.z);
}

// adapeted  stb_perlin.h
template <typename T>
inline T perlin_fbm(const vec<T, 3>& p, T lacunarity, T gain, int octaves,
    const vec<int, 3>& wrap) {
  return _stb_perlin_fbm_noise3((float)p.x, (float)p.y, (float)p.z,
      (float)lacunarity, (float)gain, octaves, wrap.x, wrap.y, wrap.z);
}

// adapeted  stb_perlin.h
template <typename T>
inline T perlin_turbulence(const vec<T, 3>& p, T lacunarity, T gain,
    int octaves, const vec<int, 3>& wrap) {
  return _stb_perlin_turbulence_noise3((float)p.x, (float)p.y, (float)p.z,
      (float)lacunarity, (float)gain, octaves, wrap.x, wrap.y, wrap.z);
}

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHADING FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto::math {

// Schlick approximation of the Fresnel term
template <typename T>
inline vec<T, 3> fresnel_schlick(const vec<T, 3>& specular,
    const vec<T, 3>& normal, const vec<T, 3>& outgoing) {
  if (specular == vec<T, 3>{0}) return vec<T, 3>{0};
  auto cosine = dot(normal, outgoing);
  return specular + (1 - specular) * pow(clamp(1 - abs(cosine), (T)0, (T)1), 5);
}

// Compute the fresnel term for dielectrics.
template <typename T>
inline T fresnel_dielectric(
    T eta, const vec<T, 3>& normal, const vec<T, 3>& outgoing) {
  // Implementation from
  // https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
  auto cosw = abs(dot(normal, outgoing));

  auto sin2 = 1 - cosw * cosw;
  auto eta2 = eta * eta;

  auto cos2t = 1 - sin2 / eta2;
  if (cos2t < 0) return 1;  // tir

  auto t0 = sqrt(cos2t);
  auto t1 = eta * t0;
  auto t2 = eta * cosw;

  auto rs = (cosw - t1) / (cosw + t1);
  auto rp = (t0 - t2) / (t0 + t2);

  return (rs * rs + rp * rp) / 2;
}

// Compute the fresnel term for metals.
template <typename T>
inline vec<T, 3> fresnel_conductor(const vec<T, 3>& eta, const vec<T, 3>& etak,
    const vec<T, 3>& normal, const vec<T, 3>& outgoing) {
  // Implementation from
  // https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
  auto cosw = dot(normal, outgoing);
  if (cosw <= 0) return vec<T, 3>{0};

  cosw       = clamp(cosw, (T)-1, (T)1);
  auto cos2  = cosw * cosw;
  auto sin2  = clamp(1 - cos2, (T)0, (T)1);
  auto eta2  = eta * eta;
  auto etak2 = etak * etak;

  auto t0       = eta2 - etak2 - sin2;
  auto a2plusb2 = sqrt(t0 * t0 + 4 * eta2 * etak2);
  auto t1       = a2plusb2 + cos2;
  auto a        = sqrt((a2plusb2 + t0) / 2);
  auto t2       = 2 * a * cosw;
  auto rs       = (t1 - t2) / (t1 + t2);

  auto t3 = cos2 * a2plusb2 + sin2 * sin2;
  auto t4 = t2 * sin2;
  auto rp = rs * (t3 - t4) / (t3 + t4);

  return (rp + rs) / 2;
}

// Convert eta to reflectivity
template <typename T>
inline vec<T, 3> eta_to_reflectivity(const vec<T, 3>& eta) {
  return ((eta - 1) * (eta - 1)) / ((eta + 1) * (eta + 1));
}
// Convert reflectivity to  eta.
template <typename T>
inline vec<T, 3> reflectivity_to_eta(const vec<T, 3>& reflectivity_) {
  auto reflectivity = clamp(reflectivity_, 0.0f, 0.99f);
  return (1 + sqrt(reflectivity)) / (1 - sqrt(reflectivity));
}
// Convert conductor eta to reflectivity
template <typename T>
inline vec<T, 3> eta_to_reflectivity(
    const vec<T, 3>& eta, const vec<T, 3>& etak) {
  return ((eta - 1) * (eta - 1) + etak * etak) /
         ((eta + 1) * (eta + 1) + etak * etak);
}
// Convert eta to edge tint parametrization
template <typename T>
inline std::pair<vec<T, 3>, vec<T, 3>> eta_to_edgetint(
    const vec<T, 3>& eta, const vec<T, 3>& etak) {
  auto reflectivity = eta_to_reflectivity(eta, etak);
  auto numer        = (1 + sqrt(reflectivity)) / (1 - sqrt(reflectivity)) - eta;
  auto denom        = (1 + sqrt(reflectivity)) / (1 - sqrt(reflectivity)) -
               (1 - reflectivity) / (1 + reflectivity);
  auto edgetint = numer / denom;
  return {reflectivity, edgetint};
}
// Convert reflectivity and edge tint to eta.
template <typename T>
inline std::pair<vec<T, 3>, vec<T, 3>> edgetint_to_eta(
    const vec<T, 3>& reflectivity, const vec<T, 3>& edgetint) {
  auto r = clamp(reflectivity, 0.0f, 0.99f);
  auto g = edgetint;

  auto r_sqrt = sqrt(r);
  auto n_min  = (1 - r) / (1 + r);
  auto n_max  = (1 + r_sqrt) / (1 - r_sqrt);

  auto n  = lerp(n_max, n_min, g);
  auto k2 = ((n + 1) * (n + 1) * r - (n - 1) * (n - 1)) / (1 - r);
  k2      = max(k2, 0.0f);
  auto k  = sqrt(k2);
  return {n, k};
}

// Evaluate microfacet distribution
template <typename T>
inline T microfacet_distribution(
    T roughness, const vec<T, 3>& normal, const vec<T, 3>& halfway, bool ggx) {
  // https://google.github.io/filament/Filament.html#materialsystem/specularbrdf
  // http://graphicrants.blogspot.com/2013/08/specular-brdf-reference.html
  auto cosine = dot(normal, halfway);
  if (cosine <= 0) return 0;
  auto roughness2 = roughness * roughness;
  auto cosine2    = cosine * cosine;
  if (ggx) {
    return roughness2 / ((T)pi * (cosine2 * roughness2 + 1 - cosine2) *
                            (cosine2 * roughness2 + 1 - cosine2));
  } else {
    return exp((cosine2 - 1) / (roughness2 * cosine2)) /
           ((T)pi * roughness2 * cosine2 * cosine2);
  }
}

// Evaluate the microfacet shadowing1
template <typename T>
inline T microfacet_shadowing1(T roughness, const vec<T, 3>& normal,
    const vec<T, 3>& halfway, const vec<T, 3>& direction, bool ggx) {
  // https://google.github.io/filament/Filament.html#materialsystem/specularbrdf
  // http://graphicrants.blogspot.com/2013/08/specular-brdf-reference.html
  auto cosine  = dot(normal, direction);
  auto cosineh = dot(halfway, direction);
  if (cosine * cosineh <= 0) return 0;
  auto roughness2 = roughness * roughness;
  auto cosine2    = cosine * cosine;
  if (ggx) {
    return 2 * abs(cosine) /
           (abs(cosine) + sqrt(cosine2 - roughness2 * cosine2 + roughness2));
  } else {
    auto ci = abs(cosine) / (roughness * sqrt(1 - cosine2));
    return ci < (T)1.6 ? ((T)3.535 * ci + (T)2.181 * ci * ci) /
                             ((T)1.0 + (T)2.276 * ci + (T)2.577 * ci * ci)
                       : (T)1.0;
  }
}

// Evaluate microfacet shadowing
template <typename T>
inline T microfacet_shadowing(T roughness, const vec<T, 3>& normal,
    const vec<T, 3>& halfway, const vec<T, 3>& outgoing,
    const vec<T, 3>& incoming, bool ggx) {
  return microfacet_shadowing1(roughness, normal, halfway, outgoing, ggx) *
         microfacet_shadowing1(roughness, normal, halfway, incoming, ggx);
}

// Sample a microfacet ditribution.
template <typename T>
inline vec<T, 3> sample_microfacet(
    T roughness, const vec<T, 3>& normal, const vec<T, 2>& rn, bool ggx) {
  auto phi   = 2 * (T)pi * rn.x;
  auto theta = (T)0.0;
  if (ggx) {
    theta = atan(roughness * sqrt(rn.y / (1 - rn.y)));
  } else {
    auto roughness2 = roughness * roughness;
    theta           = atan(sqrt(-roughness2 * log(1 - rn.y)));
  }
  auto local_half_vector = vec<T, 3>{
      cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)};
  return transform_direction(basis_fromz(normal), local_half_vector);
}

// Pdf for microfacet distribution sampling.
template <typename T>
inline T sample_microfacet_pdf(
    T roughness, const vec<T, 3>& normal, const vec<T, 3>& halfway, bool ggx) {
  auto cosine = dot(normal, halfway);
  if (cosine < 0) return 0;
  return microfacet_distribution(roughness, normal, halfway, ggx) * cosine;
}

// Sample a microfacet ditribution with the distribution of visible normals.
template <typename T>
inline vec<T, 3> sample_microfacet(T roughness, const vec<T, 3>& normal,
    const vec<T, 3>& outgoing, const vec<T, 2>& rn, bool ggx) {
  // http://jcgt.org/published/0007/04/01/
  if (ggx) {
    // move to local coordinate system
    auto basis   = basis_fromz(normal);
    auto Ve      = transform_direction(transpose(basis), outgoing);
    auto alpha_x = roughness, alpha_y = roughness;
    // Section 3.2: transforming the view direction to the hemisphere
    // configuration
    auto Vh = normalize(vec<T, 3>{alpha_x * Ve.x, alpha_y * Ve.y, Ve.z});
    // Section 4.1: orthonormal basis (with special case if cross product is
    // zero)
    auto lensq = Vh.x * Vh.x + Vh.y * Vh.y;
    auto T1    = lensq > 0 ? vec<T, 3>{-Vh.y, Vh.x, 0} * (1 / sqrt(lensq))
                        : vec<T, 3>{1, 0, 0};
    auto T2 = cross(Vh, T1);
    // Section 4.2: parameterization of the projected area
    auto r   = sqrt(rn.y);
    auto phi = 2 * (T)pi * rn.x;
    auto t1  = r * cos(phi);
    auto t2  = r * sin(phi);
    auto s   = 0.5f * (1 + Vh.z);
    t2       = (1 - s) * sqrt(1 - t1 * t1) + s * t2;
    // Section 4.3: reprojection onto hemisphere
    auto Nh = t1 * T1 + t2 * T2 + sqrt(max(1 - t1 * t1 - t2 * t2, (T)0)) * Vh;
    // Section 3.4: transforming the normal back to the ellipsoid configuration
    auto Ne = normalize(
        vec<T, 3>{alpha_x * Nh.x, alpha_y * Nh.y, max(Nh.z, (T)0)});
    // move to world coordinate
    auto local_halfway = Ne;
    return transform_direction(basis, local_halfway);
  } else {
    throw std::invalid_argument{"not implemented yet"};
  }
}

// Pdf for microfacet distribution sampling with the distribution of visible
// normals.
template <typename T>
inline T sample_microfacet_pdf(T roughness, const vec<T, 3>& normal,
    const vec<T, 3>& halfway, const vec<T, 3>& outgoing, bool ggx) {
  // http://jcgt.org/published/0007/04/01/
  if (dot(normal, halfway) < 0) return 0;
  if (dot(halfway, outgoing) < 0) return 0;
  return microfacet_distribution(roughness, normal, halfway, ggx) *
         microfacet_shadowing1(roughness, normal, halfway, outgoing, ggx) *
         max(0.0f, dot(halfway, outgoing)) / abs(dot(normal, outgoing));
}

// Evaluate a diffuse BRDF lobe.
template <typename T>
inline vec<T, 3> eval_diffuse_reflection(const vec<T, 3>& normal,
    const vec<T, 3>& outgoing, const vec<T, 3>& incoming) {
  if (dot(normal, incoming) <= 0 || dot(normal, outgoing) <= 0)
    return vec<T, 3>{0};
  return vec<T, 3>{1} / (T)pi * dot(normal, incoming);
}

// Evaluate a specular BRDF lobe.
template <typename T>
inline vec<T, 3> eval_microfacet_reflection(T ior, T roughness,
    const vec<T, 3>& normal, const vec<T, 3>& outgoing,
    const vec<T, 3>& incoming) {
  if (dot(normal, incoming) <= 0 || dot(normal, outgoing) <= 0)
    return vec<T, 3>{0};
  auto halfway = normalize(incoming + outgoing);
  auto F       = fresnel_dielectric(ior, halfway, incoming);
  auto D       = microfacet_distribution(roughness, normal, halfway);
  auto G = microfacet_shadowing(roughness, normal, halfway, outgoing, incoming);
  return vec<T, 3>{1} * F * D * G /
         (4 * dot(normal, outgoing) * dot(normal, incoming)) *
         dot(normal, incoming);
}

// Evaluate a metal BRDF lobe.
template <typename T>
inline vec<T, 3> eval_microfacet_reflection(const vec<T, 3>& eta,
    const vec<T, 3>& etak, T roughness, const vec<T, 3>& normal,
    const vec<T, 3>& outgoing, const vec<T, 3>& incoming) {
  if (dot(normal, incoming) <= 0 || dot(normal, outgoing) <= 0)
    return vec<T, 3>{0};
  auto halfway = normalize(incoming + outgoing);
  auto F       = fresnel_conductor(eta, etak, halfway, incoming);
  auto D       = microfacet_distribution(roughness, normal, halfway);
  auto G = microfacet_shadowing(roughness, normal, halfway, outgoing, incoming);
  return F * D * G / (4 * dot(normal, outgoing) * dot(normal, incoming)) *
         dot(normal, incoming);
}

// Evaluate a transmission BRDF lobe.
template <typename T>
inline vec<T, 3> eval_microfacet_transmission(T ior, T roughness,
    const vec<T, 3>& normal, const vec<T, 3>& outgoing,
    const vec<T, 3>& incoming) {
  if (dot(normal, incoming) >= 0 || dot(normal, outgoing) <= 0)
    return vec<T, 3>{0};
  auto reflected = reflect(-incoming, normal);
  auto halfway   = normalize(reflected + outgoing);
  // auto F       = fresnel_schlick(
  //     point.reflectance, abs(dot(halfway, outgoing)), entering);
  auto D = microfacet_distribution(roughness, normal, halfway);
  auto G = microfacet_shadowing(
      roughness, normal, halfway, outgoing, reflected);
  return vec<T, 3>{1} * D * G /
         (4 * dot(normal, outgoing) * dot(normal, reflected)) *
         (dot(normal, reflected));
}

// Evaluate a refraction BRDF lobe.
template <typename T>
inline vec<T, 3> eval_microfacet_refraction(T ior, T roughness,
    const vec<T, 3>& normal, const vec<T, 3>& outgoing,
    const vec<T, 3>& incoming) {
  auto entering  = dot(normal, outgoing) >= 0;
  auto up_normal = entering ? normal : -normal;
  auto rel_ior   = entering ? ior : (1 / ior);
  if (dot(normal, incoming) * dot(normal, outgoing) >= 0) {
    auto halfway = normalize(incoming + outgoing);
    auto F       = fresnel_dielectric(rel_ior, halfway, outgoing);
    auto D       = microfacet_distribution(roughness, up_normal, halfway);
    auto G       = microfacet_shadowing(
        roughness, up_normal, halfway, outgoing, incoming);
    return vec<T, 3>{1} * F * D * G /
           abs(4 * dot(normal, outgoing) * dot(normal, incoming)) *
           abs(dot(normal, incoming));
  } else {
    auto halfway = -normalize(rel_ior * incoming + outgoing) *
                   (entering ? 1 : -1);
    auto F = fresnel_dielectric(rel_ior, halfway, outgoing);
    auto D = microfacet_distribution(roughness, up_normal, halfway);
    auto G = microfacet_shadowing(
        roughness, up_normal, halfway, outgoing, incoming);
    // [Walter 2007] equation 21
    return vec<T, 3>{1} *
           abs((dot(outgoing, halfway) * dot(incoming, halfway)) /
               (dot(outgoing, normal) * dot(incoming, normal))) *
           (1 - F) * D * G /
           pow(rel_ior * dot(halfway, incoming) + dot(halfway, outgoing), 2) *
           abs(dot(normal, incoming));
  }
}

// Sample a diffuse BRDF lobe.
template <typename T>
inline vec<T, 3> sample_diffuse_reflection(
    const vec<T, 3>& normal, const vec<T, 3>& outgoing, const vec<T, 2>& rn) {
  if (dot(normal, outgoing) <= 0) return vec<T, 3>{0};
  return sample_hemisphere_cos(normal, rn);
}

// Sample a specular BRDF lobe.
template <typename T>
inline vec<T, 3> sample_microfacet_reflection(T ior, T roughness,
    const vec<T, 3>& normal, const vec<T, 3>& outgoing, const vec<T, 2>& rn) {
  if (dot(normal, outgoing) <= 0) return vec<T, 3>{0};
  // auto halfway = sample_microfacet(roughness, normal, outgoing, rn);
  auto halfway = sample_microfacet(roughness, normal, rn);
  return reflect(outgoing, halfway);
}

// Sample a metal BRDF lobe.
template <typename T>
inline vec<T, 3> sample_microfacet_reflection(const vec<T, 3>& eta,
    const vec<T, 3>& etak, T roughness, const vec<T, 3>& normal,
    const vec<T, 3>& outgoing, const vec<T, 2>& rn) {
  if (dot(normal, outgoing) <= 0) return vec<T, 3>{0};
  // auto halfway = sample_microfacet(roughness, normal, outgoing, rn);
  auto halfway = sample_microfacet(roughness, normal, rn);
  return reflect(outgoing, halfway);
}

// Sample a transmission BRDF lobe.
template <typename T>
inline vec<T, 3> sample_microfacet_transmission(T ior, T roughness,
    const vec<T, 3>& normal, const vec<T, 3>& outgoing, const vec<T, 2>& rn) {
  if (dot(normal, outgoing) <= 0) return vec<T, 3>{0};
  auto halfway = sample_microfacet(roughness, normal, rn);
  // auto halfway   = sample_microfacet(roughness, normal, outgoing, rn);
  auto reflected = reflect(outgoing, halfway);
  return -reflect(reflected, normal);
}

// Sample a refraction BRDF lobe.
template <typename T>
inline vec<T, 3> sample_microfacet_refraction(T ior, T roughness,
    const vec<T, 3>& normal, const vec<T, 3>& outgoing, T rnl,
    const vec<T, 2>& rn) {
  auto entering  = dot(normal, outgoing) >= 0;
  auto up_normal = entering ? normal : -normal;
  // auto halfway   = sample_microfacet(roughness, up_normal, outgoing, rn);
  auto halfway = sample_microfacet(roughness, up_normal, rn);
  if (rnl < fresnel_dielectric(entering ? ior : (1 / ior), halfway, outgoing)) {
    return reflect(outgoing, halfway);
  } else {
    return refract(outgoing, halfway, entering ? (1 / ior) : ior);
  }
}

// Pdf for diffuse BRDF lobe sampling.
template <typename T>
inline T sample_diffuse_reflection_pdf(const vec<T, 3>& normal,
    const vec<T, 3>& outgoing, const vec<T, 3>& incoming) {
  if (dot(normal, incoming) <= 0 || dot(normal, outgoing) <= 0) return 0;
  return sample_hemisphere_cos_pdf(normal, incoming);
}

// Pdf for specular BRDF lobe sampling.
template <typename T>
inline T sample_microfacet_reflection_pdf(T ior, T roughness,
    const vec<T, 3>& normal, const vec<T, 3>& outgoing,
    const vec<T, 3>& incoming) {
  if (dot(normal, incoming) <= 0 || dot(normal, outgoing) <= 0) return 0;
  auto halfway = normalize(outgoing + incoming);
  // return sample_microfacet_pdf(roughness, normal, halfway, outgoing) /
  return sample_microfacet_pdf(roughness, normal, halfway) /
         (4 * abs(dot(outgoing, halfway)));
}

// Pdf for metal BRDF lobe sampling.
template <typename T>
inline T sample_microfacet_reflection_pdf(const vec<T, 3>& eta,
    const vec<T, 3>& etak, T roughness, const vec<T, 3>& normal,
    const vec<T, 3>& outgoing, const vec<T, 3>& incoming) {
  if (dot(normal, incoming) <= 0 || dot(normal, outgoing) <= 0) return 0;
  auto halfway = normalize(outgoing + incoming);
  // return sample_microfacet_pdf(roughness, normal, halfway, outgoing) /
  return sample_microfacet_pdf(roughness, normal, halfway) /
         (4 * abs(dot(outgoing, halfway)));
}

// Pdf for transmission BRDF lobe sampling.
template <typename T>
inline T sample_microfacet_transmission_pdf(T ior, T roughness,
    const vec<T, 3>& normal, const vec<T, 3>& outgoing,
    const vec<T, 3>& incoming) {
  if (dot(normal, incoming) >= 0 || dot(normal, outgoing) <= 0) return 0;
  auto reflected = reflect(-incoming, normal);
  auto halfway   = normalize(reflected + outgoing);
  // auto d         = sample_microfacet_pdf(roughness, normal, halfway,
  // outgoing);
  auto d = sample_microfacet_pdf(roughness, normal, halfway);
  return d / (4 * abs(dot(outgoing, halfway)));
}

// Pdf for refraction BRDF lobe sampling.
template <typename T>
inline T sample_microfacet_refraction_pdf(T ior, T roughness,
    const vec<T, 3>& normal, const vec<T, 3>& outgoing,
    const vec<T, 3>& incoming) {
  auto entering  = dot(normal, outgoing) >= 0;
  auto up_normal = entering ? normal : -normal;
  auto rel_ior   = entering ? ior : (1 / ior);
  if (dot(normal, incoming) * dot(normal, outgoing) >= 0) {
    auto halfway = normalize(incoming + outgoing);
    return fresnel_dielectric(rel_ior, halfway, outgoing) *
           //  sample_microfacet_pdf(roughness, up_normal, halfway, outgoing) /
           sample_microfacet_pdf(roughness, up_normal, halfway) /
           (4 * abs(dot(outgoing, halfway)));
  } else {
    auto halfway = -normalize(rel_ior * incoming + outgoing) *
                   (entering ? 1 : -1);
    // [Walter 2007] equation 17
    return (1 - fresnel_dielectric(rel_ior, halfway, outgoing)) *
           //  sample_microfacet_pdf(roughness, up_normal, halfway, outgoing) *
           sample_microfacet_pdf(roughness, up_normal, halfway) *
           abs(dot(halfway, outgoing)) /
           pow(rel_ior * dot(halfway, incoming) + dot(halfway, outgoing), 2);
  }
}

// Evaluate a delta specular BRDF lobe.
template <typename T>
inline vec<T, 3> eval_delta_reflection(T ior, const vec<T, 3>& normal,
    const vec<T, 3>& outgoing, const vec<T, 3>& incoming) {
  if (dot(normal, incoming) <= 0 || dot(normal, outgoing) <= 0)
    return vec<T, 3>{0};
  return vec<T, 3>{1} * fresnel_dielectric(ior, normal, outgoing);
}

// Evaluate a delta metal BRDF lobe.
template <typename T>
inline vec<T, 3> eval_delta_reflection(const vec<T, 3>& eta,
    const vec<T, 3>& etak, const vec<T, 3>& normal, const vec<T, 3>& outgoing,
    const vec<T, 3>& incoming) {
  if (dot(normal, incoming) <= 0 || dot(normal, outgoing) <= 0)
    return vec<T, 3>{0};
  return fresnel_conductor(eta, etak, normal, outgoing);
}

// Evaluate a delta transmission BRDF lobe.
template <typename T>
inline vec<T, 3> eval_delta_transmission(T ior, const vec<T, 3>& normal,
    const vec<T, 3>& outgoing, const vec<T, 3>& incoming) {
  if (dot(normal, incoming) >= 0 || dot(normal, outgoing) <= 0)
    return vec<T, 3>{0};
  return vec<T, 3>{1};
}

// Evaluate a delta refraction BRDF lobe.
template <typename T>
inline vec<T, 3> eval_delta_refraction(T ior, const vec<T, 3>& normal,
    const vec<T, 3>& outgoing, const vec<T, 3>& incoming) {
  auto entering  = dot(normal, outgoing) >= 0;
  auto up_normal = entering ? normal : -normal;
  auto rel_ior   = entering ? ior : (1 / ior);
  if (dot(normal, incoming) * dot(normal, outgoing) >= 0) {
    return vec<T, 3>{1} * fresnel_dielectric(rel_ior, up_normal, outgoing);
  } else {
    return vec<T, 3>{1} *
           (1 - fresnel_dielectric(rel_ior, up_normal, outgoing));
  }
}

// Sample a delta specular BRDF lobe.
template <typename T>
inline vec<T, 3> sample_delta_reflection(
    T ior, const vec<T, 3>& normal, const vec<T, 3>& outgoing) {
  if (dot(normal, outgoing) <= 0) return vec<T, 3>{0};
  return reflect(outgoing, normal);
}

// Sample a delta metal BRDF lobe.
template <typename T>
inline vec<T, 3> sample_delta_reflection(const vec<T, 3>& eta,
    const vec<T, 3>& etak, const vec<T, 3>& normal, const vec<T, 3>& outgoing) {
  if (dot(normal, outgoing) <= 0) return vec<T, 3>{0};
  return reflect(outgoing, normal);
}

// Sample a delta transmission BRDF lobe.
template <typename T>
inline vec<T, 3> sample_delta_transmission(
    T ior, const vec<T, 3>& normal, const vec<T, 3>& outgoing) {
  if (dot(normal, outgoing) <= 0) return vec<T, 3>{0};
  return -outgoing;
}

// Sample a delta refraction BRDF lobe.
template <typename T>
inline vec<T, 3> sample_delta_refraction(
    T ior, const vec<T, 3>& normal, const vec<T, 3>& outgoing, T rnl) {
  auto entering  = dot(normal, outgoing) >= 0;
  auto up_normal = entering ? normal : -normal;
  auto rel_ior   = entering ? ior : (1 / ior);
  if (rnl < fresnel_dielectric(rel_ior, up_normal, outgoing)) {
    return reflect(outgoing, up_normal);
  } else {
    return refract(outgoing, up_normal, 1 / rel_ior);
  }
}

// Pdf for delta specular BRDF lobe sampling.
template <typename T>
inline T sample_delta_reflection_pdf(T ior, const vec<T, 3>& normal,
    const vec<T, 3>& outgoing, const vec<T, 3>& incoming) {
  if (dot(normal, incoming) <= 0 || dot(normal, outgoing) <= 0) return 0;
  return 1;
}

// Pdf for delta metal BRDF lobe sampling.
template <typename T>
inline T sample_delta_reflection_pdf(const vec<T, 3>& eta,
    const vec<T, 3>& etak, const vec<T, 3>& normal, const vec<T, 3>& outgoing,
    const vec<T, 3>& incoming) {
  if (dot(normal, incoming) <= 0 || dot(normal, outgoing) <= 0) return 0;
  return 1;
}

// Pdf for delta transmission BRDF lobe sampling.
template <typename T>
inline T sample_delta_transmission_pdf(T ior, const vec<T, 3>& normal,
    const vec<T, 3>& outgoing, const vec<T, 3>& incoming) {
  if (dot(normal, incoming) >= 0 || dot(normal, outgoing) <= 0) return 0;
  return 1;
}

// Pdf for delta refraction BRDF lobe sampling.
template <typename T>
inline T sample_delta_refraction_pdf(T ior, const vec<T, 3>& normal,
    const vec<T, 3>& outgoing, const vec<T, 3>& incoming) {
  auto entering  = dot(normal, outgoing) >= 0;
  auto up_normal = entering ? normal : -normal;
  auto rel_ior   = entering ? ior : (1 / ior);
  if (dot(normal, incoming) * dot(normal, outgoing) >= 0) {
    return fresnel_dielectric(rel_ior, up_normal, outgoing);
  } else {
    return (1 - fresnel_dielectric(rel_ior, up_normal, outgoing));
  }
}

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF MONETACARLO SAMPLING FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto::math {

// Sample an hemispherical direction with uniform distribution.
template <typename T>
inline vec<T, 3> sample_hemisphere(const vec<T, 2>& ruv) {
  auto z   = ruv.y;
  auto r   = sqrt(clamp(1 - z * z, 0, 1));
  auto phi = 2 * (T)pi * ruv.x;
  return {r * cos(phi), r * sin(phi), z};
}
template <typename T>
inline T sample_hemisphere_pdf(const vec<T, 3>& direction) {
  return (direction.z <= 0) ? 0 : 1 / (2 * (T)pi);
}

// Sample an hemispherical direction with uniform distribution.
template <typename T>
inline vec<T, 3> sample_hemisphere(
    const vec<T, 3>& normal, const vec<T, 2>& ruv) {
  auto z               = ruv.y;
  auto r               = sqrt(clamp(1 - z * z, (T)0, (T)1));
  auto phi             = 2 * (T)pi * ruv.x;
  auto local_direction = vec<T, 3>{r * cos(phi), r * sin(phi), z};
  return transform_direction(basis_fromz(normal), local_direction);
}
template <typename T>
inline T sample_hemisphere_pdf(
    const vec<T, 3>& normal, const vec<T, 3>& direction) {
  return (dot(normal, direction) <= 0) ? 0 : 1 / (2 * (T)pi);
}

// Sample a spherical direction with uniform distribution.
template <typename T>
inline vec<T, 3> sample_sphere(const vec<T, 2>& ruv) {
  auto z   = 2 * ruv.y - 1;
  auto r   = sqrt(clamp(1 - z * z, (T)0, (T)1));
  auto phi = 2 * (T)pi * ruv.x;
  return {r * cos(phi), r * sin(phi), z};
}
template <typename T>
inline T sample_sphere_pdf(const vec<T, 3>& w) {
  return 1 / (4 * (T)pi);
}

// Sample an hemispherical direction with cosine distribution.
template <typename T>
inline vec<T, 3> sample_hemisphere_cos(const vec<T, 2>& ruv) {
  auto z   = sqrt(ruv.y);
  auto r   = sqrt(1 - z * z);
  auto phi = 2 * (T)pi * ruv.x;
  return {r * cos(phi), r * sin(phi), z};
}
template <typename T>
inline T sample_hemisphere_cos_pdf(const vec<T, 3>& direction) {
  return (direction.z <= 0) ? 0 : direction.z / (T)pi;
}

// Sample an hemispherical direction with cosine distribution.
template <typename T>
inline vec<T, 3> sample_hemisphere_cos(
    const vec<T, 3>& normal, const vec<T, 2>& ruv) {
  auto z               = sqrt(ruv.y);
  auto r               = sqrt(1 - z * z);
  auto phi             = 2 * (T)pi * ruv.x;
  auto local_direction = vec<T, 3>{r * cos(phi), r * sin(phi), z};
  return transform_direction(basis_fromz(normal), local_direction);
}
template <typename T>
inline T sample_hemisphere_cos_pdf(
    const vec<T, 3>& normal, const vec<T, 3>& direction) {
  auto cosw = dot(normal, direction);
  return (cosw <= 0) ? 0 : cosw / (T)pi;
}

// Sample an hemispherical direction with cosine power distribution.
template <typename T>
inline vec<T, 3> sample_hemisphere_cospower(T exponent, const vec<T, 2>& ruv) {
  auto z   = pow(ruv.y, 1 / (exponent + 1));
  auto r   = sqrt(1 - z * z);
  auto phi = 2 * (T)pi * ruv.x;
  return {r * cos(phi), r * sin(phi), z};
}
template <typename T>
inline T sample_hemisphere_cospower_pdf(
    T exponent, const vec<T, 3>& direction) {
  return (direction.z <= 0)
             ? 0
             : pow(direction.z, exponent) * (exponent + 1) / (2 * (T)pi);
}

// Sample a point uniformly on a disk.
template <typename T>
inline vec<T, 2> sample_disk(const vec<T, 2>& ruv) {
  auto r   = sqrt(ruv.y);
  auto phi = 2 * (T)pi * ruv.x;
  return {cos(phi) * r, sin(phi) * r};
}
template <typename T>
inline T sample_disk_pdf(const vec<T, 2>& point) {
  return 1 / (T)pi;
}

// Sample a point uniformly on a cylinder, without caps.
template <typename T>
inline vec<T, 3> sample_cylinder(const vec<T, 2>& ruv) {
  auto phi = 2 * (T)pi * ruv.x;
  return {sin(phi), cos(phi), ruv.y * 2 - 1};
}
template <typename T>
inline T sample_cylinder_pdf(const vec<T, 3>& point) {
  return 1 / (T)pi;
}

// Sample a point uniformly on a triangle returning the baricentric coordinates.
template <typename T>
inline vec<T, 2> sample_triangle(const vec<T, 2>& ruv) {
  return {1 - sqrt(ruv.x), ruv.y * sqrt(ruv.x)};
}

// Sample a point uniformly on a triangle.
template <typename T>
inline vec<T, 3> sample_triangle(const vec<T, 3>& p0, const vec<T, 3>& p1,
    const vec<T, 3>& p2, const vec<T, 2>& ruv) {
  auto uv = sample_triangle(ruv);
  return p0 * (1 - uv.x - uv.y) + p1 * uv.x + p2 * uv.y;
}
// Pdf for uniform triangle sampling, i.e. triangle area.
template <typename T>
inline T sample_triangle_pdf(
    const vec<T, 3>& p0, const vec<T, 3>& p1, const vec<T, 3>& p2) {
  return 2 / length(cross(p1 - p0, p2 - p0));
}

// Sample an index with uniform distribution.
template <typename T>
inline int sample_uniform(int size, T r) {
  return clamp((int)(r * size), 0, size - 1);
}
template <typename T>
inline T sample_uniform_pdf(int size) {
  return (T)1 / (T)size;
}

// Sample an index with uniform distribution.
template <typename T>
inline T sample_uniform(const std::vector<T>& elements, T r) {
  if (elements.empty()) return {};
  auto size = (int)elements.size();
  return elements[clamp((int)(r * size), 0, size - 1)];
}
template <typename T>
inline T sample_uniform_pdf(const std::vector<T>& elements) {
  if (elements.empty()) return 0;
  return (T)1 / (int)elements.size();
}

// Sample a discrete distribution represented by its cdf.
template <typename T>
inline int sample_discrete(const std::vector<T>& cdf, T r) {
  r        = clamp(r * cdf.back(), (T)0, cdf.back() - (T)0.00001);
  auto idx = (int)(std::upper_bound(cdf.data(), cdf.data() + cdf.size(), r) -
                   cdf.data());
  return clamp(idx, 0, (int)cdf.size() - 1);
}
// Pdf for uniform discrete distribution sampling.
template <typename T>
inline T sample_discrete_pdf(const std::vector<T>& cdf, int idx) {
  if (idx == 0) return cdf.at(0);
  return cdf.at(idx) - cdf.at(idx - 1);
}

}  // namespace yocto::math

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF USER INTERFACE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto::math {

// Computes the image uv coordinates corresponding to the view parameters.
// Returns negative coordinates if out of the image.
template <typename T>
inline vec<int, 2> get_image_coords(const vec<T, 2>& mouse_pos,
    const vec<T, 2>& center, T scale, const vec<int, 2>& txt_size) {
  auto xyf = (mouse_pos - center) / scale;
  return vec2i{(int)round(xyf.x + (float)txt_size.x / 2),
      (int)round(xyf.y + (float)txt_size.y / 2)};
}

// Center image and autofit.
template <typename T>
inline void update_imview(vec<T, 2>& center, T& scale,
    const vec<int, 2>& imsize, const vec<int, 2>& winsize, bool zoom_to_fit) {
  if (zoom_to_fit) {
    scale  = min(winsize.x / (float)imsize.x, winsize.y / (float)imsize.y);
    center = {(float)winsize.x / 2, (float)winsize.y / 2};
  } else {
    if (winsize.x >= imsize.x * scale) center.x = winsize.x / 2;
    if (winsize.y >= imsize.y * scale) center.y = winsize.y / 2;
  }
}

// Turntable for UI navigation.
template <typename T>
inline void update_turntable(vec<T, 3>& from, vec<T, 3>& to, vec<T, 3>& up,
    const vec<T, 2>& rotate, T dolly, const vec<T, 2>& pan) {
  // rotate if necessary
  if (rotate.x || rotate.y) {
    auto z     = normalize(to - from);
    auto lz    = length(to - from);
    auto phi   = atan2(z.z, z.x) + rotate.x;
    auto theta = acos(z.y) + rotate.y;
    theta      = clamp(theta, (T)0.001, (T)pi - (T)0.001);
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
template <typename T>
inline void update_turntable(frame<T, 3>& frame, T& focus,
    const vec<T, 2>& rotate, T dolly, const vec<T, 2>& pan) {
  // rotate if necessary
  if (rotate != zero2f) {
    auto phi   = atan2(frame.z.z, frame.z.x) + rotate.x;
    auto theta = acos(frame.z.y) + rotate.y;
    theta      = clamp(theta, (T)0.001, (T)pi - (T)0.001);
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
    focus   = max(focus * (1 + dolly), (T)0.001);
    frame.o = c + frame.z * focus;
  }

  // pan if necessary
  if (pan.x || pan.y) {
    frame.o += frame.x * pan.x + frame.y * pan.y;
  }
}

// FPS camera for UI navigation for a frame parametrization.
template <typename T>
inline void update_fpscam(
    frame<T, 3>& frame, const vec<T, 3>& transl, const vec<T, 2>& rotate) {
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

// Generate a ray from a camera
template <typename T>
inline ray<T, 3> camera_ray(const frame<T, 3>& frame, T lens,
    const vec<T, 2>& film, const vec<T, 2>& image_uv) {
  auto e = vec3f{0};
  auto q = vec3f{
      film.x * (0.5f - image_uv.x), film.y * (image_uv.y - 0.5f), lens};
  auto q1  = -q;
  auto d   = normalize(q1 - e);
  auto ray = ray3f{transform_point(frame, e), transform_direction(frame, d)};
  return ray;
}

}  // namespace yocto::math

#endif
