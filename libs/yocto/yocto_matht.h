//
// # Yocto/Math: Math types
//
// Yocto/Math defines the basic math primitives used in graphics, including
// small-sized vectors, matrices, frames, quaternions, rays, bounding boxes
// and their transforms. Yocto/Math is implemented in `yocto_math.h`.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2022 Fabio Pellacini
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

#ifndef _YOCTO_MATH_H_
#define _YOCTO_MATH_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <utility>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::pair;

}  // namespace yocto

// -----------------------------------------------------------------------------
// MATH CONSTANTS AND FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

using byte   = unsigned char;
using uint   = unsigned int;
using ushort = unsigned short;

inline const double pi  = 3.14159265358979323846;
inline const float  pif = (float)pi;

constexpr auto int_max = std::numeric_limits<int>::max();
constexpr auto int_min = std::numeric_limits<int>::lowest();
constexpr auto flt_max = std::numeric_limits<float>::max();
constexpr auto flt_min = std::numeric_limits<float>::lowest();
constexpr auto flt_eps = std::numeric_limits<float>::epsilon();

template <typename T>
constexpr auto num_max = std::numeric_limits<T>::max();
template <typename T>
constexpr auto num_min = std::numeric_limits<T>::lowest();
template <typename T>
constexpr auto num_eps = std::numeric_limits<T>::epsilon();

using std::swap;

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
inline T sqr(T a);
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
inline bool isfinite(T a);
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
// template <typename T>
// inline void swap(T& a, T& b);
template <typename T>
inline T smoothstep(T a, T b, T u);
template <typename T>
inline T bias(T a, T bias);
template <typename T>
inline T gain(T a, T gain);
template <typename I>
inline int pow2(int a);

}  // namespace yocto

// -----------------------------------------------------------------------------
// VECTORS
// -----------------------------------------------------------------------------
namespace yocto {

template <typename T, int N>
struct vec;

template <typename T>
struct vec<T, 1> {
  T x = 0;

  T&       operator[](int i);
  const T& operator[](int i) const;
};

template <typename T>
struct vec<T, 2> {
  T x = 0;
  T y = 0;

  T&       operator[](int i);
  const T& operator[](int i) const;
};

template <typename T>
struct vec<T, 3> {
  T x = 0;
  T y = 0;
  T z = 0;

  T&       operator[](int i);
  const T& operator[](int i) const;
};

template <typename T>
struct vec<T, 4> {
  T x = 0;
  T y = 0;
  T z = 0;
  T w = 0;

  T&       operator[](int i);
  const T& operator[](int i) const;
};

// Vector aliases
using vec1f = vec<float, 1>;
using vec2f = vec<float, 2>;
using vec3f = vec<float, 3>;
using vec4f = vec<float, 4>;
using vec1i = vec<int, 1>;
using vec2i = vec<int, 2>;
using vec3i = vec<int, 3>;
using vec4i = vec<int, 4>;
using vec4b = vec<byte, 4>;

// Element access
template <typename T>
inline vec<T, 3> xyz(const vec<T, 4>& a);

// Vector sequence operations.
template <typename T, int N>
inline int size(const vec<T, N>& a);
template <typename T, int N>
inline const T* begin(const vec<T, N>& a);
template <typename T, int N>
inline const T* end(const vec<T, N>& a);
template <typename T, int N>
inline T* begin(vec<T, N>& a);
template <typename T, int N>
inline T* end(vec<T, N>& a);
template <typename T, int N>
inline const T* data(const vec<T, N>& a);
template <typename T, int N>
inline T* data(vec<T, N>& a);

// Vector comparison operations.
template <typename T, int N>
inline bool operator==(const vec<T, N>& a, const vec<T, N>& b);
template <typename T, int N>
inline bool operator!=(const vec<T, N>& a, const vec<T, N>& b);

// Vector operations.
template <typename T, int N>
inline vec<T, N> operator+(const vec<T, N>& a);
template <typename T, int N>
inline vec<T, N> operator-(const vec<T, N>& a);
template <typename T, int N>
inline vec<T, N> operator+(const vec<T, N>& a, const vec<T, N>& b);
template <typename T, int N, typename T1>
inline vec<T, N> operator+(const vec<T, N>& a, T1 b);
template <typename T, int N, typename T1>
inline vec<T, N> operator+(T1 a, const vec<T, N>& b);
template <typename T, int N>
inline vec<T, N> operator-(const vec<T, N>& a, const vec<T, N>& b);
template <typename T, int N, typename T1>
inline vec<T, N> operator-(const vec<T, N>& a, T1 b);
template <typename T, int N, typename T1>
inline vec<T, N> operator-(T1 a, const vec<T, N>& b);
template <typename T, int N>
inline vec<T, N> operator*(const vec<T, N>& a, const vec<T, N>& b);
template <typename T, int N, typename T1>
inline vec<T, N> operator*(const vec<T, N>& a, T1 b);
template <typename T, int N, typename T1>
inline vec<T, N> operator*(T1 a, const vec<T, N>& b);
template <typename T, int N>
inline vec<T, N> operator/(const vec<T, N>& a, const vec<T, N>& b);
template <typename T, int N, typename T1>
inline vec<T, N> operator/(const vec<T, N>& a, T1 b);
template <typename T, int N, typename T1>
inline vec<T, N> operator/(T1 a, const vec<T, N>& b);

// Vector assignments
template <typename T, int N>
inline vec<T, N>& operator+=(vec<T, N>& a, const vec<T, N>& b);
template <typename T, int N, typename T1>
inline vec<T, N>& operator+=(vec<T, N>& a, T1 b);
template <typename T, int N>
inline vec<T, N>& operator-=(vec<T, N>& a, const vec<T, N>& b);
template <typename T, int N, typename T1>
inline vec<T, N>& operator-=(vec<T, N>& a, T1 b);
template <typename T, int N>
inline vec<T, N>& operator*=(vec<T, N>& a, const vec<T, N>& b);
template <typename T, int N, typename T1>
inline vec<T, N>& operator*=(vec<T, N>& a, T1 b);
template <typename T, int N>
inline vec<T, N>& operator/=(vec<T, N>& a, const vec<T, N>& b);
template <typename T, int N, typename T1>
inline vec<T, N>& operator/=(vec<T, N>& a, T1 b);

// Vector products and lengths.
template <typename T, int N>
inline T dot(const vec<T, N>& a, const vec<T, N>& b);
template <typename T>
inline T cross(const vec<T, 2>& a, const vec<T, 2>& b);
template <typename T>
inline vec<T, 3> cross(const vec<T, 3>& a, const vec<T, 3>& b);

template <typename T, int N>
inline T length(const vec<T, N>& a);
template <typename T, int N>
inline T length_squared(const vec<T, N>& a);
template <typename T, int N>
inline vec<T, N> normalize(const vec<T, N>& a);
template <typename T, int N>
inline T distance(const vec<T, N>& a, const vec<T, N>& b);
template <typename T, int N>
inline T distance_squared(const vec<T, N>& a, const vec<T, N>& b);
template <typename T, int N>
inline T angle(const vec<T, N>& a, const vec<T, N>& b);

// Orthogonal vectors.
template <typename T>
inline vec<T, 3> orthogonal(const vec<T, 3>& v);
template <typename T>
inline vec<T, 3> orthonormalize(const vec<T, 3>& a, const vec<T, 3>& b);

// Reflected and refracted vector.
template <typename T>
inline vec<T, 3> reflect(const vec<T, 3>& w, const vec<T, 3>& n);
template <typename T>
inline vec<T, 3> refract(const vec<T, 3>& w, const vec<T, 3>& n, T inv_eta);

// Slerp
template <typename T>
inline vec<T, 4> slerp(const vec<T, 4>& a, const vec<T, 4>& b, T u);

// Max element and clamp.
template <typename T, int N>
inline vec<T, N> max(const vec<T, N>& a, T b);
template <typename T, int N>
inline vec<T, N> min(const vec<T, N>& a, T b);
template <typename T, int N>
inline vec<T, N> max(const vec<T, N>& a, const vec<T, N>& b);
template <typename T, int N>
inline vec<T, N> min(const vec<T, N>& a, const vec<T, N>& b);
template <typename T, int N>
inline vec<T, N> clamp(const vec<T, N>& x, T min, T max);
template <typename T, int N>
inline vec<T, N> lerp(const vec<T, N>& a, const vec<T, N>& b, T u);
template <typename T, int N>
inline vec<T, N> lerp(
    const vec<T, N>& a, const vec<T, N>& b, const vec<T, N>& u);

template <typename T, int N>
inline T max(const vec<T, N>& a);
template <typename T, int N>
inline T min(const vec<T, N>& a);
template <typename T, int N>
inline T sum(const vec<T, N>& a);
template <typename T, int N>
inline T mean(const vec<T, N>& a);

// Functions applied to vector elements
template <typename T, int N>
inline vec<T, N> abs(const vec<T, N>& a);
template <typename T, int N>
inline vec<T, N> sqr(const vec<T, N>& a);
template <typename T, int N>
inline vec<T, N> sqrt(const vec<T, N>& a);
template <typename T, int N>
inline vec<T, N> exp(const vec<T, N>& a);
template <typename T, int N>
inline vec<T, N> log(const vec<T, N>& a);
template <typename T, int N>
inline vec<T, N> exp2(const vec<T, N>& a);
template <typename T, int N>
inline vec<T, N> log2(const vec<T, N>& a);
template <typename T, int N>
inline bool isfinite(const vec<T, N>& a);
template <typename T, int N>
inline vec<T, N> pow(const vec<T, N>& a, T b);
template <typename T, int N>
inline vec<T, N> pow(const vec<T, N>& a, const vec<T, N>& b);
template <typename T, int N>
inline vec<T, N> gain(const vec<T, N>& a, T b);
// template <typename T, int N>
// inline void swap(vec<T, N>& a, vec<T, N>& b);

// Quaternion operatons represented as xi + yj + zk + w
// const auto identity_quat4f = vec4f{0, 0, 0, 1};
template <typename T>
inline vec<T, 4> quat_mul(const vec<T, 4>& a, T b);
template <typename T>
inline vec<T, 4> quat_mul(const vec<T, 4>& a, const vec<T, 4>& b);
template <typename T>
inline vec<T, 4> quat_conjugate(const vec<T, 4>& a);
template <typename T>
inline vec<T, 4> quat_inverse(const vec<T, 4>& a);

}  // namespace yocto

// -----------------------------------------------------------------------------
// MATRICES
// -----------------------------------------------------------------------------
namespace yocto {

// Small Fixed-size matrices stored in column major format.
template <typename T, int N>
struct mat;

// Small Fixed-size matrices stored in column major format.
template <typename T>
struct mat<T, 1> {
  vec<T, 1> x = {1};

  vec<T, 1>&       operator[](int i);
  const vec<T, 1>& operator[](int i) const;
};

// Small Fixed-size matrices stored in column major format.
template <typename T>
struct mat<T, 2> {
  vec<T, 2> x = {1, 0};
  vec<T, 2> y = {0, 1};

  vec<T, 2>&       operator[](int i);
  const vec<T, 2>& operator[](int i) const;
};

// Small Fixed-size matrices stored in column major format.
template <typename T>
struct mat<T, 3> {
  vec<T, 3> x = {1, 0, 0};
  vec<T, 3> y = {0, 1, 0};
  vec<T, 3> z = {0, 0, 1};

  vec<T, 3>&       operator[](int i);
  const vec<T, 3>& operator[](int i) const;
};

// Small Fixed-size matrices stored in column major format.
template <typename T>
struct mat<T, 4> {
  vec<T, 4> x = {1, 0, 0, 0};
  vec<T, 4> y = {0, 1, 0, 0};
  vec<T, 4> z = {0, 0, 1, 0};
  vec<T, 4> w = {0, 0, 0, 1};

  vec<T, 4>&       operator[](int i);
  const vec<T, 4>& operator[](int i) const;
};

// Matrix aliases
using mat1f = mat<float, 1>;
using mat2f = mat<float, 2>;
using mat3f = mat<float, 3>;
using mat4f = mat<float, 4>;

// Matrix comparisons.
template <typename T, int N>
inline bool operator==(const mat<T, N>& a, const mat<T, N>& b);
template <typename T, int N>
inline bool operator!=(const mat<T, N>& a, const mat<T, N>& b);

// Matrix operations.
template <typename T, int N>
inline mat<T, N> operator+(const mat<T, N>& a, const mat<T, N>& b);
template <typename T, int N, typename T1>
inline mat<T, N> operator*(const mat<T, N>& a, T1 b);
template <typename T, int N>
inline vec<T, N> operator*(const mat<T, N>& a, const vec<T, N>& b);
template <typename T, int N>
inline vec<T, N> operator*(const vec<T, N>& a, const mat<T, N>& b);
template <typename T, int N>
inline mat<T, N> operator*(const mat<T, N>& a, const mat<T, N>& b);

// Matrix assignments.
template <typename T, int N>
inline mat<T, N>& operator+=(mat<T, N>& a, const mat<T, N>& b);
template <typename T, int N>
inline mat<T, N>& operator*=(mat<T, N>& a, const mat<T, N>& b);
template <typename T, int N, typename T1>
inline mat<T, N>& operator*=(mat<T, N>& a, T1 b);

// Matrix diagonals and transposes.
template <typename T, int N>
inline vec<T, N> diagonal(const mat<T, N>& a);
template <typename T, int N>
inline mat<T, N> transpose(const mat<T, N>& a);

// Matrix adjoints, determinants and inverses.
template <typename T, int N>
inline T determinant(const mat<T, N>& a);
template <typename T, int N>
inline mat<T, N> adjoint(const mat<T, N>& a);
template <typename T, int N>
inline mat<T, N> inverse(const mat<T, N>& a);

// Constructs a basis from a direction
template <typename T>
inline mat<T, 3> basis_fromz(const vec<T, 3>& v);

}  // namespace yocto

// -----------------------------------------------------------------------------
// RIGID BODY TRANSFORMS/FRAMES
// -----------------------------------------------------------------------------
namespace yocto {

// Rigid frames stored as a column-major affine transform matrix.
template <typename T, int N>
struct frame;

// Rigid frames stored as a column-major affine transform matrix.
template <typename T>
struct frame<T, 2> {
  vec<T, 2> x = {1, 0};
  vec<T, 2> y = {0, 1};
  vec<T, 2> o = {0, 0};

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

  vec<T, 3>&       operator[](int i);
  const vec<T, 3>& operator[](int i) const;
};

// Frame aliases
using frame2f = frame<float, 2>;
using frame3f = frame<float, 3>;

// Frame properties
template <typename T, int N>
inline mat<T, N> rotation(const frame<T, N>& a);
template <typename T, int N>
inline vec<T, N> translation(const frame<T, N>& a);

// Frame construction
template <typename T, int N>
inline frame<T, N> make_frame(const mat<T, 2>& m, const vec<T, 2>& t);

// Conversion between frame and mat
template <typename T, int N>
inline mat<T, N + 1> frame_to_mat(const frame<T, N>& f);
template <typename T, int N>
inline frame<T, N - 1> mat_to_frame(const mat<T, N>& ma);

// Frame comparisons.
template <typename T, int N>
inline bool operator==(const frame<T, N>& a, const frame<T, N>& b);
template <typename T, int N>
inline bool operator!=(const frame<T, N>& a, const frame<T, N>& b);

// Frame composition, equivalent to affine matrix product.
template <typename T, int N>
inline frame<T, N> operator*(const frame<T, N>& a, const frame<T, N>& b);
template <typename T, int N>
inline frame<T, N>& operator*=(frame<T, N>& a, const frame<T, N>& b);

// Frame inverse, equivalent to rigid affine inverse.
template <typename T, int N>
inline frame<T, N> inverse(const frame<T, N>& a, bool non_rigid = false);

// Frame construction from axis.
template <typename T>
inline frame<T, 3> frame_fromz(const vec<T, 3>& o, const vec<T, 3>& v);
template <typename T, int N>
inline frame<T, 3> frame_fromzx(
    const vec<T, 3>& o, const vec<T, 3>& z_, const vec<T, 3>& x_);

}  // namespace yocto

// -----------------------------------------------------------------------------
// QUATERNIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Quaternions to represent rotations
template <typename T, int N>
struct quat;

// Quaternions to represent rotations
template <typename T>
struct quat<T, 4> {
  T x = 0;
  T y = 0;
  T z = 0;
  T w = 1;
};

// Quaternion aliases
using quat4f = quat<float, 4>;

// Quaternion operatons
template <typename T>
inline quat<T, 4> operator+(const quat<T, 4>& a, const quat<T, 4>& b);
template <typename T>
inline quat<T, 4> operator*(const quat<T, 4>& a, T b);
template <typename T>
inline quat<T, 4> operator/(const quat<T, 4>& a, T b);
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
template <typename T>
inline quat<T, 4> lerp(const quat<T, 4>& a, const quat<T, 4>& b, T t);
template <typename T>
inline quat<T, 4> nlerp(const quat<T, 4>& a, const quat<T, 4>& b, T t);
template <typename T>
inline quat<T, 4> slerp(const quat<T, 4>& a, const quat<T, 4>& b, T t);

}  // namespace yocto

// -----------------------------------------------------------------------------
// TRANSFORMS
// -----------------------------------------------------------------------------
namespace yocto {

// Transforms points, vectors and directions by matrices.
template <typename T, int N>
inline vec<T, N> transform_point(const mat<T, N + 1>& a, const vec<T, N>& b);
template <typename T, int N>
inline vec<T, N> transform_vector(const mat<T, N + 1>& a, const vec<T, N>& b);
template <typename T, int N>
inline vec<T, N> transform_direction(
    const mat<T, N + 1>& a, const vec<T, N>& b);
template <typename T, int N>
inline vec<T, N> transform_normal(const mat<T, N + 1>& a, const vec<T, N>& b);
template <typename T, int N>
inline vec<T, N> transform_vector(const mat<T, N>& a, const vec<T, N>& b);
template <typename T, int N>
inline vec<T, N> transform_direction(const mat<T, N>& a, const vec<T, N>& b);
template <typename T, int N>
inline vec<T, N> transform_normal(const mat<T, N>& a, const vec<T, N>& b);

// Transforms points, vectors and directions by frames.
template <typename T, int N>
inline vec<T, N> transform_point(const frame<T, N>& a, const vec<T, N>& b);
template <typename T, int N>
inline vec<T, N> transform_vector(const frame<T, N>& a, const vec<T, N>& b);
template <typename T, int N>
inline vec<T, N> transform_direction(const frame<T, N>& a, const vec<T, N>& b);
template <typename T, int N>
inline vec<T, N> transform_normal(
    const frame<T, N>& a, const vec<T, N>& b, bool non_rigid = false);

// Transforms points, vectors and directions by frames.
template <typename T, int N>
inline vec<T, N> transform_point_inverse(
    const frame<T, N>& a, const vec<T, N>& b);
template <typename T, int N>
inline vec<T, N> transform_vector_inverse(
    const frame<T, N>& a, const vec<T, N>& b);
template <typename T, int N>
inline vec<T, N> transform_direction_inverse(
    const frame<T, N>& a, const vec<T, N>& b);

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
inline frame<T, 3> rotation_frame(const quat<T, 4>& quat);
template <typename T>
inline frame<T, 3> rotation_frame(const mat<T, 3>& rot);

// Lookat frame. Z-axis can be inverted with inv_xz.
template <typename T>
inline frame<T, 3> lookat_frame(const vec<T, 3>& eye, const vec<T, 3>& center,
    const vec<T, 3>& up, bool inv_xz = false);

// OpenGL frustum, ortho and perspecgive matrices.
template <typename T>
inline mat<T, 4> frustum_mat(T l, T r, T b, T t, T n, T f);
template <typename T>
inline mat<T, 4> ortho_mat(T l, T r, T b, T t, T n, T f);
template <typename T>
inline mat<T, 4> ortho2d_mat(T left, T right, T bottom, T top);
template <typename T>
inline mat<T, 4> ortho_mat(T xmag, T ymag, T near, T far);
template <typename T>
inline mat<T, 4> perspective_mat(T fovy, T aspect, T near, T far);
template <typename T>
inline mat<T, 4> perspective_mat(T fovy, T aspect, T near);

// Rotation conversions.
template <typename T>
inline pair<vec<T, 3>, T> rotation_axisangle(const vec<T, 4>& quat);
template <typename T>
inline vec<T, 4> rotation_quat(const vec<T, 3>& axis, T angle);
template <typename T>
inline vec<T, 4> rotation_quat(const vec<T, 4>& axisangle);

}  // namespace yocto

// -----------------------------------------------------------------------------
// USER INTERFACE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Computes the image uv coordinates corresponding to the view parameters.
// Returns negative coordinates if out of the image.
template <typename T, typename I>
inline vec<I, 2> image_coords(const vec<T, 2>& mouse_pos,
    const vec<T, 2>& center, T scale, const vec<I, 2>& txt_size);

// Center image and autofit. Returns center and scale.
template <typename T, typename I>
inline pair<vec<T, 2>, T> camera_imview(const vec<T, 2>& center, T scale,
    const vec<I, 2>& imsize, const vec<I, 2>& winsize, bool zoom_to_fit);

// Turntable for UI navigation. Returns from and to.
template <typename T>
inline pair<vec<T, 3>, vec<T, 3>> camera_turntable(const vec<T, 3>& from,
    const vec<T, 3>& to, const vec<T, 3>& up, const vec<T, 2>& rotate, T dolly,
    const vec<T, 2>& pan);

// Turntable for UI navigation. Returns frame and focus.
template <typename T>
inline pair<frame<T, 3>, T> camera_turntable(const frame<T, 3>& frame, T focus,
    const vec<T, 2>& rotate, T dolly, const vec<T, 2>& pan);

// FPS camera for UI navigation for a frame parametrization. Returns frame.
template <typename T>
inline frame<T, 3> camera_fpscam(
    const frame<T, 3>& frame, const vec<T, 3>& transl, const vec<T, 2>& rotate);

// Computes the image uv coordinates corresponding to the view parameters.
// Returns negative coordinates if out of the image.
template <typename T, typename I>
[[deprecated]] inline vec<I, 2> get_image_coords(const vec<T, 2>& mouse_pos,
    const vec<T, 2>& center, T scale, const vec<I, 2>& txt_size);

// Center image and autofit.
template <typename T, typename I>
[[deprecated]] inline void update_imview(vec<T, 2>& center, T& scale,
    const vec<I, 2>& imsize, const vec<I, 2>& winsize, bool zoom_to_fit);

// Turntable for UI navigation.
template <typename T>
[[deprecated]] inline void update_turntable(vec<T, 3>& from, vec<T, 3>& to,
    const vec<T, 3>& up, const vec<T, 2>& rotate, T dolly,
    const vec<T, 2>& pan);

// Turntable for UI navigation.
template <typename T>
[[deprecated]] inline void update_turntable(frame<T, 3>& frame, T& focus,
    const vec<T, 2>& rotate, T dolly, const vec<T, 2>& pan);

// FPS camera for UI navigation for a frame parametrization.
template <typename T>
[[deprecated]] inline void update_fpscam(
    frame<T, 3>& frame, const vec<T, 3>& transl, const vec<T, 2>& rotate);

}  // namespace yocto

// -----------------------------------------------------------------------------
// PYTHON-LIKE ITERATORS
// -----------------------------------------------------------------------------
namespace yocto {

// Python range. Construct an object that iterates over an integer sequence.
template <typename T>
constexpr auto range(T max);
template <typename T>
constexpr auto range(T min, T max);
template <typename T>
constexpr auto range(T min, T max, T step);

// Python enumerate
template <typename Sequence, typename T = size_t>
constexpr auto enumerate(const Sequence& sequence, T start = 0);
template <typename Sequence, typename T = size_t>
constexpr auto enumerate(Sequence& sequence, T start = 0);

// Python zip
template <typename Sequence1, typename Sequence2>
constexpr auto zip(const Sequence1& sequence1, const Sequence2& sequence2);
template <typename Sequence1, typename Sequence2>
constexpr auto zip(Sequence1& sequence1, Sequence2& sequence2);
template <typename Sequence1, typename Sequence2>
constexpr auto zip(const Sequence1& sequence1, Sequence2& sequence2);
template <typename Sequence1, typename Sequence2>
constexpr auto zip(Sequence1& sequence1, const Sequence2& sequence2);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SIGNED-SIZE
// -----------------------------------------------------------------------------
namespace yocto {

template <typename T>
inline std::ptrdiff_t ssize(const T& container);

}

// -----------------------------------------------------------------------------
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// MATH CONSTANTS AND FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

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
  return a < 0 ? (T)-1 : (T)1;
}
template <typename T>
inline T sqr(T a) {
  return a * a;
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
inline bool isfinite(T a) {
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
// template <typename T>
// inline void swap(T& a, T& b) {
//   std::swap(a, b);
// }
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
  return u < a ? (T)0 : (T)1;
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
  return (a < (T)0.5) ? bias(a * 2, gain) / 2
                      : bias(a * 2 - 1, 1 - gain) / 2 + (T)0.5;
}
template <typename I>
inline I pow2(I a) {
  return 1 << a;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// VECTORS
// -----------------------------------------------------------------------------
namespace yocto {

// Vec2
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
inline T& vec<T, 3>::operator[](int i) {
  return (&x)[i];
}
template <typename T>
inline const T& vec<T, 3>::operator[](int i) const {
  return (&x)[i];
}

// Vec4
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
inline vec<T, 3> xyz(const vec<T, 4>& a) {
  return {a.x, a.y, a.z};
}

// Vector sequence operations.
template <typename T, int N>
inline int size(const vec<T, N>& a) {
  return N;
}
template <typename T, int N>
inline const T* begin(const vec<T, N>& a) {
  return &a.x;
}
template <typename T, int N>
inline const T* end(const vec<T, N>& a) {
  return &a.x + N;
}
template <typename T, int N>
inline T* begin(vec<T, N>& a) {
  return &a.x;
}
template <typename T, int N>
inline T* end(vec<T, N>& a) {
  return &a.x + N;
}
template <typename T, int N>
inline const T* data(const vec<T, N>& a) {
  return &a.x;
}
template <typename T, int N>
inline T* data(vec<T, N>& a) {
  return &a.x;
}

// Vector comparison operations.
template <typename T, int N>
inline bool operator==(const vec<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return a.x == b.x;
  } else if constexpr (N == 2) {
    return a.x == b.x && a.y == b.y;
  } else if constexpr (N == 3) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
  } else if constexpr (N == 4) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
  }
}
template <typename T, int N>
inline bool operator!=(const vec<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return a.x != b.x;
  } else if constexpr (N == 2) {
    return a.x != b.x || a.y != b.y;
  } else if constexpr (N == 3) {
    return a.x != b.x || a.y != b.y || a.z != b.z;
  } else if constexpr (N == 4) {
    return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w;
  }
}

// Vector operations.
template <typename T, int N>
inline vec<T, N> operator+(const vec<T, N>& a) {
  return a;
}
template <typename T, int N>
inline vec<T, N> operator-(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return {-a.x};
  } else if constexpr (N == 2) {
    return {-a.x, -a.y};
  } else if constexpr (N == 3) {
    return {-a.x, -a.y, -a.z};
  } else if constexpr (N == 4) {
    return {-a.x, -a.y, -a.z, -a.w};
  }
}

template <typename T, int N>
inline vec<T, N> operator+(const vec<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return {a.x + b.x};
  } else if constexpr (N == 2) {
    return {a.x + b.x, a.y + b.y};
  } else if constexpr (N == 3) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
  } else if constexpr (N == 4) {
    return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
  }
}
template <typename T, int N, typename T1>
inline vec<T, N> operator+(const vec<T, N>& a, T1 b) {
  if constexpr (N == 1) {
    return {a.x + b};
  } else if constexpr (N == 2) {
    return {a.x + b, a.y + b};
  } else if constexpr (N == 3) {
    return {a.x + b, a.y + b, a.z + b};
  } else if constexpr (N == 4) {
    return {a.x + b, a.y + b, a.z + b, a.w + b};
  }
}
template <typename T, int N, typename T1>
inline vec<T, N> operator+(T1 a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return {a + b.x};
  } else if constexpr (N == 2) {
    return {a + b.x, a + b.y};
  } else if constexpr (N == 3) {
    return {a + b.x, a + b.y, a + b.z};
  } else if constexpr (N == 4) {
    return {a + b.x, a + b.y, a + b.z, a + b.w};
  }
}
template <typename T, int N>
inline vec<T, N> operator-(const vec<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return {a.x - b.x};
  } else if constexpr (N == 2) {
    return {a.x - b.x, a.y - b.y};
  } else if constexpr (N == 3) {
    return {a.x - b.x, a.y - b.y, a.z - b.z};
  } else if constexpr (N == 4) {
    return {a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w};
  }
}
template <typename T, int N, typename T1>
inline vec<T, N> operator-(const vec<T, N>& a, T1 b) {
  if constexpr (N == 1) {
    return {a.x - b};
  } else if constexpr (N == 2) {
    return {a.x - b, a.y - b};
  } else if constexpr (N == 3) {
    return {a.x - b, a.y - b, a.z - b};
  } else if constexpr (N == 4) {
    return {a.x - b, a.y - b, a.z - b, a.w - b};
  }
}
template <typename T, int N, typename T1>
inline vec<T, N> operator-(T1 a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return {a - b.x};
  } else if constexpr (N == 2) {
    return {a - b.x, a - b.y};
  } else if constexpr (N == 3) {
    return {a - b.x, a - b.y, a - b.z};
  } else if constexpr (N == 4) {
    return {a - b.x, a - b.y, a - b.z, a - b.w};
  }
}
template <typename T, int N>
inline vec<T, N> operator*(const vec<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return {a.x * b.x};
  } else if constexpr (N == 2) {
    return {a.x * b.x, a.y * b.y};
  } else if constexpr (N == 3) {
    return {a.x * b.x, a.y * b.y, a.z * b.z};
  } else if constexpr (N == 4) {
    return {a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w};
  }
}
template <typename T, int N, typename T1>
inline vec<T, N> operator*(const vec<T, N>& a, T1 b) {
  if constexpr (N == 1) {
    return {a.x * b};
  } else if constexpr (N == 2) {
    return {a.x * b, a.y * b};
  } else if constexpr (N == 3) {
    return {a.x * b, a.y * b, a.z * b};
  } else if constexpr (N == 4) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
  }
}
template <typename T, int N, typename T1>
inline vec<T, N> operator*(T1 a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return {a * b.x};
  } else if constexpr (N == 2) {
    return {a * b.x, a * b.y};
  } else if constexpr (N == 3) {
    return {a * b.x, a * b.y, a * b.z};
  } else if constexpr (N == 4) {
    return {a * b.x, a * b.y, a * b.z, a * b.w};
  }
}
template <typename T, int N>
inline vec<T, N> operator/(const vec<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return {a.x / b.x};
  } else if constexpr (N == 2) {
    return {a.x / b.x, a.y / b.y};
  } else if constexpr (N == 3) {
    return {a.x / b.x, a.y / b.y, a.z / b.z};
  } else if constexpr (N == 4) {
    return {a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w};
  }
}
template <typename T, int N, typename T1>
inline vec<T, N> operator/(const vec<T, N>& a, T1 b) {
  if constexpr (N == 1) {
    return {a.x / b};
  } else if constexpr (N == 2) {
    return {a.x / b, a.y / b};
  } else if constexpr (N == 3) {
    return {a.x / b, a.y / b, a.z / b};
  } else if constexpr (N == 4) {
    return {a.x / b, a.y / b, a.z / b, a.w / b};
  }
}
template <typename T, int N, typename T1>
inline vec<T, N> operator/(T1 a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return {a / b.x};
  } else if constexpr (N == 2) {
    return {a / b.x, a / b.y};
  } else if constexpr (N == 3) {
    return {a / b.x, a / b.y, a / b.z};
  } else if constexpr (N == 4) {
    return {a / b.x, a / b.y, a / b.z, a / b.w};
  }
}

// Vector assignments
template <typename T, int N>
inline vec<T, N>& operator+=(vec<T, N>& a, const vec<T, N>& b) {
  return a = a + b;
}
template <typename T, int N, typename T1>
inline vec<T, N>& operator+=(vec<T, N>& a, T1 b) {
  return a = a + b;
}
template <typename T, int N>
inline vec<T, N>& operator-=(vec<T, N>& a, const vec<T, N>& b) {
  return a = a - b;
}
template <typename T, int N, typename T1>
inline vec<T, N>& operator-=(vec<T, N>& a, T1 b) {
  return a = a - b;
}
template <typename T, int N>
inline vec<T, N>& operator*=(vec<T, N>& a, const vec<T, N>& b) {
  return a = a * b;
}
template <typename T, int N, typename T1>
inline vec<T, N>& operator*=(vec<T, N>& a, T1 b) {
  return a = a * b;
}
template <typename T, int N>
inline vec<T, N>& operator/=(vec<T, N>& a, const vec<T, N>& b) {
  return a = a / b;
}
template <typename T, int N, typename T1>
inline vec<T, N>& operator/=(vec<T, N>& a, T1 b) {
  return a = a / b;
}

// Vector products and lengths.
template <typename T, int N>
inline T dot(const vec<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return a.x * b.x;
  } else if constexpr (N == 2) {
    return a.x * b.x + a.y * b.y;
  } else if constexpr (N == 3) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
  } else if constexpr (N == 4) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
  }
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
  return acos(clamp(dot(normalize(a), normalize(b)), (T)-1, (T)1));
}

// Orthogonal vectors.
template <typename T>
inline vec<T, 3> orthogonal(const vec<T, 3>& v) {
  // http://lolengine.net/blog/2013/09/21/picking-orthogonal-vector-combing-coconuts)
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

template <typename T, int N>
inline T length(const vec<T, N>& a) {
  return sqrt(dot(a, a));
}
template <typename T, int N>
inline T length_squared(const vec<T, N>& a) {
  return dot(a, a);
}
template <typename T, int N>
inline vec<T, N> normalize(const vec<T, N>& a) {
  auto l = length(a);
  return (l != 0) ? a / l : a;
}
template <typename T, int N>
inline T distance(const vec<T, N>& a, const vec<T, N>& b) {
  return length(a - b);
}
template <typename T, int N>
inline T distance_squared(const vec<T, N>& a, const vec<T, N>& b) {
  return dot(a - b, a - b);
}
template <typename T, int N>
inline T angle(const vec<T, N>& a, const vec<T, N>& b) {
  return acos(clamp(dot(normalize(a), normalize(b)), (T)-1, (T)1));
}

template <typename T>
inline vec<T, 4> slerp(const vec<T, 4>& a, const vec<T, 4>& b, T u) {
  // https://en.wikipedia.org/wiki/Slerp
  auto an = normalize(a), bn = normalize(b);
  auto d = dot(an, bn);
  if (d < 0) {
    bn = -bn;
    d  = -d;
  }
  if (d > (T)0.9995) return normalize(an + u * (bn - an));
  auto th = acos(clamp(d, (T)-1, (T)1));
  if (th == 0) return an;
  return an * (sin(th * (1 - u)) / sin(th)) + bn * (sin(th * u) / sin(th));
}

// Max element and clamp.
template <typename T, int N>
inline vec<T, N> max(const vec<T, N>& a, T b) {
  if constexpr (N == 1) {
    return {max(a.x, b)};
  } else if constexpr (N == 2) {
    return {max(a.x, b), max(a.y, b)};
  } else if constexpr (N == 3) {
    return {max(a.x, b), max(a.y, b), max(a.z, b)};
  } else if constexpr (N == 4) {
    return {max(a.x, b), max(a.y, b), max(a.z, b), max(a.w, b)};
  }
}
template <typename T, int N>
inline vec<T, N> min(const vec<T, N>& a, T b) {
  if constexpr (N == 1) {
    return {min(a.x, b)};
  } else if constexpr (N == 2) {
    return {min(a.x, b), min(a.y, b)};
  } else if constexpr (N == 3) {
    return {min(a.x, b), min(a.y, b), min(a.z, b)};
  } else if constexpr (N == 4) {
    return {min(a.x, b), min(a.y, b), min(a.z, b), min(a.w, b)};
  }
}
template <typename T, int N>
inline vec<T, N> max(const vec<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return {max(a.x, b.x)};
  } else if constexpr (N == 2) {
    return {max(a.x, b.x), max(a.y, b.y)};
  } else if constexpr (N == 3) {
    return {max(a.x, b.x), max(a.y, b.y), max(a.z, b.z)};
  } else if constexpr (N == 4) {
    return {max(a.x, b.x), max(a.y, b.y), max(a.z, b.z), max(a.w, b.w)};
  }
}
template <typename T, int N>
inline vec<T, N> min(const vec<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return {min(a.x, b.x)};
  } else if constexpr (N == 2) {
    return {min(a.x, b.x), min(a.y, b.y)};
  } else if constexpr (N == 3) {
    return {min(a.x, b.x), min(a.y, b.y), min(a.z, b.z)};
  } else if constexpr (N == 4) {
    return {min(a.x, b.x), min(a.y, b.y), min(a.z, b.z), min(a.w, b.w)};
  }
}
template <typename T, int N>
inline vec<T, N> clamp(const vec<T, N>& x, T min, T max) {
  if constexpr (N == 1) {
    return {clamp(x.x, min, max)};
  } else if constexpr (N == 2) {
    return {clamp(x.x, min, max), clamp(x.y, min, max)};
  } else if constexpr (N == 3) {
    return {clamp(x.x, min, max), clamp(x.y, min, max), clamp(x.z, min, max)};
  } else if constexpr (N == 4) {
    return {clamp(x.x, min, max), clamp(x.y, min, max), clamp(x.z, min, max),
        clamp(x.w, min, max)};
  }
}
template <typename T, int N>
inline vec<T, N> lerp(const vec<T, N>& a, const vec<T, N>& b, T u) {
  return a * (1 - u) + b * u;
}
template <typename T, int N>
inline vec<T, N> lerp(
    const vec<T, N>& a, const vec<T, N>& b, const vec<T, N>& u) {
  return a * (1 - u) + b * u;
}

template <typename T, int N>
inline T max(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return a.x;
  } else if constexpr (N == 2) {
    return max(a.x, a.y);
  } else if constexpr (N == 3) {
    return max(max(a.x, a.y), a.z);
  } else if constexpr (N == 4) {
    return max(max(max(a.x, a.y), a.z), a.w);
  }
}
template <typename T, int N>
inline T min(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return a.x;
  } else if constexpr (N == 2) {
    return min(a.x, a.y);
  } else if constexpr (N == 3) {
    return min(min(a.x, a.y), a.z);
  } else if constexpr (N == 4) {
    return min(min(min(a.x, a.y), a.z), a.w);
  }
}
template <typename T, int N>
inline T sum(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return a.x;
  } else if constexpr (N == 2) {
    return a.x + a.y;
  } else if constexpr (N == 3) {
    return a.x + a.y + a.z;
  } else if constexpr (N == 4) {
    return a.x + a.y + a.z + a.w;
  }
}
template <typename T, int N>
inline T mean(const vec<T, N>& a) {
  return sum(a) / N;
}

// Functions applied to vector elements
template <typename T, int N>
inline vec<T, N> abs(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return {abs(a.x)};
  } else if constexpr (N == 2) {
    return {abs(a.x), abs(a.y)};
  } else if constexpr (N == 3) {
    return {abs(a.x), abs(a.y), abs(a.z)};
  } else if constexpr (N == 4) {
    return {abs(a.x), abs(a.y), abs(a.z), abs(a.w)};
  }
}
template <typename T, int N>
inline vec<T, N> sqr(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return {sqr(a.x)};
  } else if constexpr (N == 2) {
    return {sqr(a.x), sqr(a.y)};
  } else if constexpr (N == 3) {
    return {sqr(a.x), sqr(a.y), sqr(a.z)};
  } else if constexpr (N == 4) {
    return {sqr(a.x), sqr(a.y), sqr(a.z), sqr(a.w)};
  }
}
template <typename T, int N>
inline vec<T, N> sqrt(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return {sqrt(a.x)};
  } else if constexpr (N == 2) {
    return {sqrt(a.x), sqrt(a.y)};
  } else if constexpr (N == 3) {
    return {sqrt(a.x), sqrt(a.y), sqrt(a.z)};
  } else if constexpr (N == 4) {
    return {sqrt(a.x), sqrt(a.y), sqrt(a.z), sqrt(a.w)};
  }
}
template <typename T, int N>
inline vec<T, N> exp(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return {exp(a.x)};
  } else if constexpr (N == 2) {
    return {exp(a.x), exp(a.y)};
  } else if constexpr (N == 3) {
    return {exp(a.x), exp(a.y), exp(a.z)};
  } else if constexpr (N == 4) {
    return {exp(a.x), exp(a.y), exp(a.z), exp(a.w)};
  }
}
template <typename T, int N>
inline vec<T, N> log(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return {log(a.x)};
  } else if constexpr (N == 2) {
    return {log(a.x), log(a.y)};
  } else if constexpr (N == 3) {
    return {log(a.x), log(a.y), log(a.z)};
  } else if constexpr (N == 4) {
    return {log(a.x), log(a.y), log(a.z), log(a.w)};
  }
}
template <typename T, int N>
inline vec<T, N> exp2(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return {exp2(a.x)};
  } else if constexpr (N == 2) {
    return {exp2(a.x), exp2(a.y)};
  } else if constexpr (N == 3) {
    return {exp2(a.x), exp2(a.y), exp2(a.z)};
  } else if constexpr (N == 4) {
    return {exp2(a.x), exp2(a.y), exp2(a.z), exp2(a.w)};
  }
}
template <typename T, int N>
inline vec<T, N> log2(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return {log2(a.x)};
  } else if constexpr (N == 2) {
    return {log2(a.x), log2(a.y)};
  } else if constexpr (N == 3) {
    return {log2(a.x), log2(a.y), log2(a.z)};
  } else if constexpr (N == 4) {
    return {log2(a.x), log2(a.y), log2(a.z), log2(a.w)};
  }
}
template <typename T, int N>
inline vec<T, N> pow(const vec<T, N>& a, T b) {
  if constexpr (N == 1) {
    return {pow(a.x, b)};
  } else if constexpr (N == 2) {
    return {pow(a.x, b), pow(a.y, b)};
  } else if constexpr (N == 3) {
    return {pow(a.x, b), pow(a.y, b), pow(a.z, b)};
  } else if constexpr (N == 4) {
    return {pow(a.x, b), pow(a.y, b), pow(a.z, b), pow(a.w, b)};
  }
}
template <typename T, int N>
inline vec<T, N> pow(const vec<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return {pow(a.x, b.x)};
  } else if constexpr (N == 2) {
    return {pow(a.x, b.x), pow(a.y, b.y)};
  } else if constexpr (N == 3) {
    return {pow(a.x, b.x), pow(a.y, b.y), pow(a.z, b.z)};
  } else if constexpr (N == 4) {
    return {pow(a.x, b.x), pow(a.y, b.y), pow(a.z, b.z), pow(a.w, b.w)};
  }
}
template <typename T, int N>
inline vec<T, N> gain(const vec<T, N>& a, T b) {
  if constexpr (N == 1) {
    return {gain(a.x, b)};
  } else if constexpr (N == 2) {
    return {gain(a.x, b), gain(a.y, b)};
  } else if constexpr (N == 3) {
    return {gain(a.x, b), gain(a.y, b), gain(a.z, b)};
  } else if constexpr (N == 4) {
    return {gain(a.x, b), gain(a.y, b), gain(a.z, b), gain(a.w, b)};
  }
}
// template <typename T, int N>
// inline void swap(vec<T, N>& a, vec<T, N>& b) {
//   std::swap(a, b);
// }
template <typename T, int N>
inline bool isfinite(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return isfinite(a.x);
  } else if constexpr (N == 2) {
    return isfinite(a.x) && isfinite(a.y);
  } else if constexpr (N == 3) {
    return isfinite(a.x) && isfinite(a.y) && isfinite(a.z);
  } else if constexpr (N == 4) {
    return isfinite(a.x) && isfinite(a.y) && isfinite(a.z) && isfinite(a.w);
  }
}

// Quaternion operatons represented as xi + yj + zk + w
// const auto identity_quat4f = vec4f{0, 0, 0, 1};
template <typename T>
inline vec<T, 4> quat_mul(const vec<T, 4>& a, T b) {
  return {a.x * b, a.y * b, a.z * b, a.w * b};
}
template <typename T>
inline vec<T, 4> quat_mul(const vec<T, 4>& a, const vec<T, 4>& b) {
  return {a.x * b.w + a.w * b.x + a.y * b.w - a.z * b.y,
      a.y * b.w + a.w * b.y + a.z * b.x - a.x * b.z,
      a.z * b.w + a.w * b.z + a.x * b.y - a.y * b.x,
      a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z};
}
template <typename T>
inline vec<T, 4> quat_conjugate(const vec<T, 4>& a) {
  return {-a.x, -a.y, -a.z, a.w};
}
template <typename T>
inline vec<T, 4> quat_inverse(const vec<T, 4>& a) {
  return quat_conjugate(a) / dot(a, a);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// MATRICES
// -----------------------------------------------------------------------------
namespace yocto {

// Small Fixed-size matrices stored in column major format.
template <typename T>
inline vec<T, 1>& mat<T, 1>::operator[](int i) {
  return (&x)[i];
}
template <typename T>
inline const vec<T, 1>& mat<T, 1>::operator[](int i) const {
  return (&x)[i];
}

// Small Fixed-size matrices stored in column major format.
template <typename T>
inline vec<T, 2>& mat<T, 2>::operator[](int i) {
  return (&x)[i];
}
template <typename T>
inline const vec<T, 2>& mat<T, 2>::operator[](int i) const {
  return (&x)[i];
}

// Small Fixed-size matrices stored in column major format.
template <typename T>
inline vec<T, 3>& mat<T, 3>::operator[](int i) {
  return (&x)[i];
}
template <typename T>
inline const vec<T, 3>& mat<T, 3>::operator[](int i) const {
  return (&x)[i];
}

// Small Fixed-size matrices stored in column major format.
template <typename T>
inline vec<T, 4>& mat<T, 4>::operator[](int i) {
  return (&x)[i];
}
template <typename T>
inline const vec<T, 4>& mat<T, 4>::operator[](int i) const {
  return (&x)[i];
}

// Matrix comparisons.
template <typename T, int N>
inline bool operator==(const mat<T, N>& a, const mat<T, N>& b) {
  if constexpr (N == 1) {
    return a.x == b.x;
  } else if constexpr (N == 2) {
    return a.x == b.x && a.y == b.y;
  } else if constexpr (N == 3) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
  } else if constexpr (N == 4) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
  }
}
template <typename T, int N>
inline bool operator!=(const mat<T, N>& a, const mat<T, N>& b) {
  return !(a == b);
}

// Matrix operations.
template <typename T, int N>
inline mat<T, N> operator+(const mat<T, N>& a, const mat<T, N>& b) {
  if constexpr (N == 1) {
    return {a.x + b.x};
  } else if constexpr (N == 2) {
    return {a.x + b.x, a.y + b.y};
  } else if constexpr (N == 3) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
  } else if constexpr (N == 4) {
    return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
  }
}
template <typename T, int N, typename T1>
inline mat<T, N> operator*(const mat<T, N>& a, T1 b) {
  if constexpr (N == 1) {
    return {a.x * b};
  } else if constexpr (N == 2) {
    return {a.x * b, a.y * b};
  } else if constexpr (N == 3) {
    return {a.x * b, a.y * b, a.z * b};
  } else if constexpr (N == 4) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
  }
}
template <typename T, int N>
inline vec<T, N> operator*(const mat<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return a.x * b.x;
  } else if constexpr (N == 2) {
    return a.x * b.x + a.y * b.y;
  } else if constexpr (N == 3) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
  } else if constexpr (N == 4) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
  }
}
template <typename T, int N>
inline vec<T, N> operator*(const vec<T, N>& a, const mat<T, N>& b) {
  if constexpr (N == 1) {
    return {dot(a, b.x)};
  } else if constexpr (N == 2) {
    return {dot(a, b.x), dot(a, b.y)};
  } else if constexpr (N == 3) {
    return {dot(a, b.x), dot(a, b.y), dot(a, b.z)};
  } else if constexpr (N == 4) {
    return {dot(a, b.x), dot(a, b.y), dot(a, b.z), dot(a, b.w)};
  }
}
template <typename T, int N>
inline mat<T, N> operator*(const mat<T, N>& a, const mat<T, N>& b) {
  if constexpr (N == 1) {
    return {a * b.x};
  } else if constexpr (N == 2) {
    return {a * b.x, a * b.y};
  } else if constexpr (N == 3) {
    return {a * b.x, a * b.y, a * b.z};
  } else if constexpr (N == 4) {
    return {a * b.x, a * b.y, a * b.z, a * b.w};
  }
}

// Matrix assignments.
template <typename T, int N>
inline mat<T, N>& operator+=(mat<T, N>& a, const mat<T, N>& b) {
  return a = a + b;
}
template <typename T, int N>
inline mat<T, N>& operator*=(mat<T, N>& a, const mat<T, N>& b) {
  return a = a * b;
}
template <typename T, int N, typename T1>
inline mat<T, N>& operator*=(mat<T, N>& a, T1 b) {
  return a = a * b;
}

// Matrix diagonals and transposes.
template <typename T, int N>
inline vec<T, N> diagonal(const mat<T, N>& a) {
  if constexpr (N == 1) {
    return {a.x.x};
  } else if constexpr (N == 2) {
    return {a.x.x, a.y.y};
  } else if constexpr (N == 3) {
    return {a.x.x, a.y.y, a.z.z};
  } else if constexpr (N == 4) {
    return {a.x.x, a.y.y, a.z.z, a.w.w};
  }
}
template <typename T, int N>
inline mat<T, N> transpose(const mat<T, N>& a) {
  if constexpr (N == 1) {
    return {{a.x.x}};
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
  }
}

// Matrix adjoints, determinants and inverses.
template <typename T, int N>
inline T determinant(const mat<T, N>& a) {
  if constexpr (N == 1) {
    return a.x;
  } else if constexpr (N == 2) {
    return cross(a.x, a.y);
  } else if constexpr (N == 3) {
    return dot(a.x, cross(a.y, a.z));
  } else if constexpr (N == 4) {
    return 0;  // TODO
  }
}
template <typename T, int N>
inline mat<T, N> adjoint(const mat<T, N>& a) {
  if constexpr (N == 1) {
    return {{a.x.x}};
  } else if constexpr (N == 2) {
    return {{a.y.y, -a.x.y}, {-a.y.x, a.x.x}};
  } else if constexpr (N == 3) {
    return transpose(
        mat<T, 3>{cross(a.y, a.z), cross(a.z, a.x), cross(a.x, a.y)});
  } else if constexpr (N == 4) {
    return {};  // TODO
  }
}
template <typename T, int N>
inline mat<T, N> inverse(const mat<T, N>& a) {
  return adjoint(a) * (1 / determinant(a));
}

// Constructs a basis from a direction
template <typename T>
inline mat<T, 3> basis_fromz(const vec<T, 3>& v) {
  // https://graphics.pixar.com/library/OrthonormalB/paper.pdf
  if constexpr (std::is_same_v<T, float>) {
    auto z    = normalize(v);
    auto sign = copysignf(1.0f, z.z);
    auto a    = -1.0f / (sign + z.z);
    auto b    = z.x * z.y * a;
    auto x    = vec<T, 3>{1.0f + sign * z.x * z.x * a, sign * b, -sign * z.x};
    auto y    = vec<T, 3>{b, sign + z.y * z.y * a, -z.y};
    return {x, y, z};
  } else if constexpr (std::is_same_v<T, float>) {
    // TODO: double
    return {};
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// RIGID BODY TRANSFORMS/FRAMES
// -----------------------------------------------------------------------------
namespace yocto {

// Rigid frames stored as a column-major affine transform matrix.
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
inline vec<T, 3>& frame<T, 3>::operator[](int i) {
  return (&x)[i];
}
template <typename T>
inline const vec<T, 3>& frame<T, 3>::operator[](int i) const {
  return (&x)[i];
}

// Frame properties
template <typename T, int N>
inline mat<T, N> rotation(const frame<T, N>& a) {
  if constexpr (N == 2) {
    return {a.x, a.y};
  } else if constexpr (N == 3) {
    return {a.x, a.y, a.z};
  }
}
template <typename T, int N>
inline vec<T, N> translation(const frame<T, N>& a) {
  if constexpr (N == 2) {
    return a.o;
  } else if constexpr (N == 3) {
    return a.o;
  }
}

// Frame construction
template <typename T, int N>
inline frame<T, N> make_frame(const mat<T, N>& m, const vec<T, N>& t) {
  if constexpr (N == 2) {
    return {m.x, m.y, t};
  } else if constexpr (N == 3) {
    return {m.x, m.y, m.z, t};
  }
}

// Frame/mat conversion
template <typename T, int N_>
inline frame<T, N_ - 1> mat_to_frame(const mat<T, N_>& m) {
  constexpr auto N = N_ - 1;
  if constexpr (N == 2) {
    return {{m.x.x, m.x.y}, {m.y.x, m.y.y}, {m.z.x, m.z.y}};
  } else if constexpr (N == 3) {
    return {{m.x.x, m.x.y, m.x.z}, {m.y.x, m.y.y, m.y.z}, {m.z.x, m.z.y, m.z.z},
        {m.w.x, m.w.y, m.w.z}};
  }
}
template <typename T, int N>
inline mat<T, N + 1> frame_to_mat(const frame<T, N>& f) {
  if constexpr (N == 2) {
    return {{f.x.x, f.x.y, 0}, {f.y.x, f.y.y, 0}, {f.o.x, f.o.y, 1}};
  } else if constexpr (N == 3) {
    return {{f.x.x, f.x.y, f.x.z, 0}, {f.y.x, f.y.y, f.y.z, 0},
        {f.z.x, f.z.y, f.z.z, 0}, {f.o.x, f.o.y, f.o.z, 1}};
  }
}

// Frame comparisons.
template <typename T, int N>
inline bool operator==(const frame<T, N>& a, const frame<T, N>& b) {
  if constexpr (N == 2) {
    return a.x == b.x && a.y == b.y && a.o == b.o;
  } else if constexpr (N == 3) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.o == b.o;
  }
}
template <typename T, int N>
inline bool operator!=(const frame<T, N>& a, const frame<T, N>& b) {
  return !(a == b);
}

// Frame composition, equivalent to affine matrix product.
template <typename T, int N>
inline frame<T, N> operator*(const frame<T, N>& a, const frame<T, N>& b) {
  return make_frame(rotation(a) * rotation(b), rotation(a) * b.o + a.o);
}
template <typename T, int N>
inline frame<T, N>& operator*=(frame<T, N>& a, const frame<T, N>& b) {
  return a = a * b;
}

// Frame inverse, equivalent to rigid affine inverse.
template <typename T, int N>
inline frame<T, N> inverse(const frame<T, N>& a, bool non_rigid) {
  if (non_rigid) {
    auto minv = inverse(rotation(a));
    return make_frame(minv, -(minv * a.o));
  } else {
    auto minv = transpose(rotation(a));
    return make_frame(minv, -(minv * a.o));
  }
}

// Frame construction from axis.
template <typename T>
inline frame<T, 3> frame_fromz(const vec<T, 3>& o, const vec<T, 3>& v) {
  // https://graphics.pixar.com/library/OrthonormalB/paper.pdf
  if constexpr (std::is_same_v<T, float>) {
    auto z    = normalize(v);
    auto sign = copysignf(1.0f, z.z);
    auto a    = -1.0f / (sign + z.z);
    auto b    = z.x * z.y * a;
    auto x    = vec<T, 3>{1.0f + sign * z.x * z.x * a, sign * b, -sign * z.x};
    auto y    = vec<T, 3>{b, sign + z.y * z.y * a, -z.y};
    return {x, y, z, o};
  } else if constexpr (std::is_same_v<T, double>) {
    // TODO: double
    return {};
  }
}
template <typename T>
inline frame<T, 3> frame_fromzx(
    const vec<T, 3>& o, const vec<T, 3>& z_, const vec<T, 3>& x_) {
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

// Quaternion operatons
template <typename T>
inline quat<T, 4> operator+(const quat<T, 4>& a, const quat<T, 4>& b) {
  return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}
template <typename T>
inline quat<T, 4> operator*(const quat<T, 4>& a, T b) {
  return {a.x * b, a.y * b, a.z * b, a.w * b};
}
template <typename T>
inline quat<T, 4> operator/(const quat<T, 4>& a, T b) {
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
inline T uangle(const quat<T, 4>& a, const quat<T, 4>& b) {
  auto d = dot(a, b);
  return d > 1 ? 0 : acos(d < -1 ? -1 : d);
}
template <typename T>
inline quat<T, 4> lerp(const quat<T, 4>& a, const quat<T, 4>& b, T t) {
  return a * (1 - t) + b * t;
}
template <typename T>
inline quat<T, 4> nlerp(const quat<T, 4>& a, const quat<T, 4>& b, T t) {
  return normalize(lerp(a, b, t));
}
template <typename T>
inline quat<T, 4> slerp(const quat<T, 4>& a, const quat<T, 4>& b, T t) {
  auto th = uangle(a, b);
  return th == 0
             ? a
             : a * (sin(th * (1 - t)) / sin(th)) + b * (sin(th * t) / sin(th));
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// TRANSFORMS
// -----------------------------------------------------------------------------
namespace yocto {

// Transforms points, vectors and directions by matrices.
template <typename T, int N>
inline vec<T, N> transform_point(const mat<T, N + 1>& a, const vec<T, N>& b) {
  if constexpr (N == 2) {
    auto tvb = a * vec<T, 3>{b.x, b.y, 1};
    return vec<T, 2>{tvb.x, tvb.y} / tvb.z;
  } else if constexpr (N == 3) {
    auto tvb = a * vec<T, 4>{b.x, b.y, b.z, 1};
    return vec<T, 3>{tvb.x, tvb.y, tvb.z} / tvb.w;
  }
}
template <typename T, int N>
inline vec<T, N> transform_vector(const mat<T, N + 1>& a, const vec<T, N>& b) {
  if constexpr (N == 2) {
    auto tvb = a * vec<T, 3>{b.x, b.y, 0};
    return vec<T, 2>{tvb.x, tvb.y} / tvb.z;
  } else if constexpr (N == 3) {
    auto tvb = a * vec<T, 4>{b.x, b.y, b.z, 0};
    return vec<T, 3>{tvb.x, tvb.y, tvb.z} / tvb.w;
  }
}
template <typename T, int N>
inline vec<T, N> transform_direction(
    const mat<T, N + 1>& a, const vec<T, N>& b) {
  return normalize(transform_vector(a, b));
}
template <typename T, int N>
inline vec<T, N> transform_normal(const mat<T, N + 1>& a, const vec<T, N>& b) {
  return normalize(transform_vector(transpose(inverse(a)), b));
}
template <typename T, int N>
inline vec<T, N> transform_vector(const mat<T, N>& a, const vec<T, N>& b) {
  return a * b;
}
template <typename T, int N>
inline vec<T, N> transform_direction(const mat<T, N>& a, const vec<T, N>& b) {
  return normalize(transform_vector(a, b));
}
template <typename T, int N>
inline vec<T, N> transform_normal(const mat<T, N>& a, const vec<T, N>& b) {
  return normalize(transform_vector(transpose(inverse(a)), b));
}

// Transforms points, vectors and directions by frames.
template <typename T, int N>
inline vec<T, N> transform_point(const frame<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 2) {
    return a.x * b.x + a.y * b.y + a.o;
  } else if constexpr (N == 3) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.o;
  }
}
template <typename T, int N>
inline vec<T, N> transform_vector(const frame<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 2) {
    return a.x * b.x + a.y * b.y;
  } else if constexpr (N == 3) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
  }
}
template <typename T, int N>
inline vec<T, N> transform_direction(const frame<T, N>& a, const vec<T, N>& b) {
  return normalize(transform_vector(a, b));
}
template <typename T, int N>
inline vec<T, N> transform_normal(
    const frame<T, N>& a, const vec<T, N>& b, bool non_rigid) {
  if (non_rigid) {
    return transform_normal(rotation(a), b);
  } else {
    return normalize(transform_vector(a, b));
  }
}

// Transforms points, vectors and directions by frames.
template <typename T, int N>
inline vec<T, N> transform_point_inverse(
    const frame<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 2) {
    return {dot(a.x, (b - a.o)), dot(a.y, (b - a.o))};
  } else if constexpr (N == 3) {
    return {dot(a.x, (b - a.o)), dot(a.y, (b - a.o)), dot(a.z, (b - a.o))};
  }
}
template <typename T, int N>
inline vec<T, N> transform_vector_inverse(
    const frame<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 2) {
    return {dot(a.x, b), dot(a.y, b)};
  } else if constexpr (N == 3) {
    return {dot(a.x, b), dot(a.y, b), dot(a.z, b)};
  }
}
template <typename T, int N>
inline vec<T, N> transform_direction_inverse(
    const frame<T, N>& a, const vec<T, N>& b) {
  return normalize(transform_vector_inverse(a, b));
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
inline frame<T, 3> rotation_frame(const mat<T, 3>& rot) {
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
inline mat<T, 4> frustum_mat(T l, T r, T b, T t, T n, T f) {
  return {{2 * n / (r - l), 0, 0, 0}, {0, 2 * n / (t - b), 0, 0},
      {(r + l) / (r - l), (t + b) / (t - b), -(f + n) / (f - n), -1},
      {0, 0, -2 * f * n / (f - n), 0}};
}
template <typename T>
inline mat<T, 4> ortho_mat(T l, T r, T b, T t, T n, T f) {
  return {{2 / (r - l), 0, 0, 0}, {0, 2 / (t - b), 0, 0},
      {0, 0, -2 / (f - n), 0},
      {-(r + l) / (r - l), -(t + b) / (t - b), -(f + n) / (f - n), 1}};
}
template <typename T>
inline mat<T, 4> ortho2d_mat(T left, T right, T bottom, T top) {
  return ortho_mat(left, right, bottom, top, -1, 1);
}
template <typename T>
inline mat<T, 4> ortho_mat(T xmag, T ymag, T near, T far) {
  return {{1 / xmag, 0, 0, 0}, {0, 1 / ymag, 0, 0}, {0, 0, 2 / (near - far), 0},
      {0, 0, (far + near) / (near - far), 1}};
}
template <typename T>
inline mat<T, 4> perspective_mat(T fovy, T aspect, T near, T far) {
  auto tg = tan(fovy / 2);
  return {{1 / (aspect * tg), 0, 0, 0}, {0, 1 / tg, 0, 0},
      {0, 0, (far + near) / (near - far), -1},
      {0, 0, 2 * far * near / (near - far), 0}};
}
template <typename T>
inline mat<T, 4> perspective_mat(T fovy, T aspect, T near) {
  auto tg = tan(fovy / 2);
  return {{1 / (aspect * tg), 0, 0, 0}, {0, 1 / tg, 0, 0}, {0, 0, -1, -1},
      {0, 0, 2 * near, 0}};
}

// Rotation conversions.
template <typename T>
inline pair<vec<T, 3>, T> rotation_axisangle(const vec<T, 4>& quat) {
  return {normalize(vec<T, 3>{quat.x, quat.y, quat.z}), 2 * acos(quat.w)};
}
template <typename T>
inline vec<T, 4> rotation_quat(const vec<T, 3>& axis, T angle) {
  auto len = length(axis);
  if (len == 0) return {0, 0, 0, 1};
  return vec<T, 4>{sin(angle / 2) * axis.x / len, sin(angle / 2) * axis.y / len,
      sin(angle / 2) * axis.z / len, cos(angle / 2)};
}
template <typename T>
inline vec<T, 4> rotation_quat(const vec<T, 4>& axisangle) {
  return rotation_quat(
      vec<T, 3>{axisangle.x, axisangle.y, axisangle.z}, axisangle.w);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// USER INTERFACE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Computes the image uv coordinates corresponding to the view parameters.
// Returns negative coordinates if out of the image.
template <typename T, typename I>
inline vec<I, 2> image_coords(const vec<T, 2>& mouse_pos,
    const vec<T, 2>& center, T scale, const vec<I, 2>& txt_size) {
  auto xyf = (mouse_pos - center) / scale;
  return vec<I, 2>{(int)round(xyf.x + txt_size.x / (T)2),
      (int)round(xyf.y + txt_size.y / (T)2)};
}

// Center image and autofit. Returns center and scale.
template <typename T, typename I>
inline pair<vec<T, 2>, T> camera_imview(const vec<T, 2>& center, T scale,
    const vec<I, 2>& imsize, const vec<I, 2>& winsize, bool zoom_to_fit) {
  if (zoom_to_fit) {
    return {{(T)winsize.x / 2, (T)winsize.y / 2},
        min(winsize.x / (T)imsize.x, winsize.y / (T)imsize.y)};
  } else {
    return {{(winsize.x >= imsize.x * scale) ? (T)winsize.x / 2 : center.x,
                (winsize.y >= imsize.y * scale) ? (T)winsize.y / 2 : center.y},
        scale};
  }
}

// Turntable for UI navigation. Returns from and to.
template <typename T>
inline pair<vec<T, 3>, vec<T, 3>> camera_turntable(const vec<T, 3>& from_,
    const vec<T, 3>& to_, const vec<T, 3>& up, const vec<T, 2>& rotate, T dolly,
    const vec<T, 2>& pan) {
  // copy values
  auto from = from_, to = to_;

  // rotate if necessary
  if (rotate != vec<T, 2>{0, 0}) {
    auto z     = normalize(to - from);
    auto lz    = length(to - from);
    auto phi   = atan2(z.z, z.x) + rotate.x;
    auto theta = acos(z.y) + rotate.y;
    theta      = clamp(theta, (T)0.001, (T)pi - (T)0.001);
    auto nz    = vec<T, 3>{sin(theta) * cos(phi) * lz, cos(theta) * lz,
        sin(theta) * sin(phi) * lz};
    from       = to - nz;
  }

  // dolly if necessary
  if (dolly != 0) {
    auto z  = normalize(to - from);
    auto lz = max((T)0.001, length(to - from) * (1 + dolly));
    z *= lz;
    from = to - z;
  }

  // pan if necessary
  if (pan != vec<T, 2>{0, 0}) {
    auto z = normalize(to - from);
    auto x = normalize(cross(up, z));
    auto y = normalize(cross(z, x));
    auto t = vec<T, 3>{pan.x * x.x + pan.y * y.x, pan.x * x.y + pan.y * y.y,
        pan.x * x.z + pan.y * y.z};
    from += t;
    to += t;
  }

  // done
  return {from, to};
}

// Turntable for UI navigation. Returns frame and focus.
template <typename T>
inline pair<frame<T, 3>, T> camera_turntable(const frame<T, 3>& frame_, T focus,
    const vec<T, 2>& rotate, T dolly, const vec<T, 2>& pan) {
  // copy values
  auto frame = frame_;

  // rotate if necessary
  if (rotate != vec<T, 2>{0, 0}) {
    auto phi   = atan2(frame.z.z, frame.z.x) + rotate.x;
    auto theta = acos(frame.z.y) + rotate.y;
    theta      = clamp(theta, (T)0.001, (T)pi - (T)0.001);
    auto new_z = vec<T, 3>{
        sin(theta) * cos(phi), cos(theta), sin(theta) * sin(phi)};
    auto new_center = frame.o - frame.z * focus;
    auto new_o      = new_center + new_z * focus;
    frame           = lookat_frame(new_o, new_center, {0, 1, 0});
    focus           = length(new_o - new_center);
  }

  // pan if necessary
  if (dolly != 0) {
    auto c  = frame.o - frame.z * focus;
    focus   = max(focus * (1 + dolly), (T)0.001);
    frame.o = c + frame.z * focus;
  }

  // pan if necessary
  if (pan != vec<T, 2>{0, 0}) {
    frame.o += frame.x * pan.x + frame.y * pan.y;
  }

  // done
  return {frame, focus};
}

// FPS camera for UI navigation for a frame parametrization. Returns frame.
template <typename T>
inline frame<T, 3> camera_fpscam(const frame<T, 3>& frame,
    const vec<T, 3>& transl, const vec<T, 2>& rotate) {
  // https://gamedev.stackexchange.com/questions/30644/how-to-keep-my-quaternion-using-fps-camera-from-tilting-and-messing-up
  auto y = vec<T, 3>{0, 1, 0};
  auto z = orthonormalize(frame.z, y);
  auto x = cross(y, z);

  auto rot = rotation_frame(vec<T, 3>{1, 0, 0}, rotate.y) *
             yocto::frame<T, 3>{frame.x, frame.y, frame.z, vec<T, 3>{0, 0, 0}} *
             rotation_frame(vec<T, 3>{0, 1, 0}, rotate.x);
  auto pos = frame.o + transl.x * x + transl.y * y + transl.z * z;

  return {rot.x, rot.y, rot.z, pos};
}

// Computes the image uv coordinates corresponding to the view parameters.
// Returns negative coordinates if out of the image.
template <typename T, typename I>
inline vec2i get_image_coords(const vec<T, 2>& mouse_pos,
    const vec<T, 2>& center, T scale, const vec<I, 2>& txt_size) {
  auto xyf = (mouse_pos - center) / scale;
  return vec2i{(int)round(xyf.x + txt_size.x / (T)2),
      (int)round(xyf.y + txt_size.y / (T)2)};
}

// Center image and autofit.
template <typename T, typename I>
inline void update_imview(vec<T, 2>& center, T& scale, const vec<I, 2>& imsize,
    const vec<I, 2>& winsize, bool zoom_to_fit) {
  if (zoom_to_fit) {
    scale  = min(winsize.x / (T)imsize.x, winsize.y / (T)imsize.y);
    center = {(T)winsize.x / 2, (T)winsize.y / 2};
  } else {
    if (winsize.x >= imsize.x * scale) center.x = (T)winsize.x / 2;
    if (winsize.y >= imsize.y * scale) center.y = (T)winsize.y / 2;
  }
}

// Turntable for UI navigation.
template <typename T>
inline void update_turntable(vec<T, 3>& from, vec<T, 3>& to, vec<T, 3>& up,
    const vec<T, 2>& rotate, T dolly, const vec<T, 2>& pan) {
  // rotate if necessary
  if (rotate != vec<T, 2>{0, 0}) {
    auto z     = normalize(to - from);
    auto lz    = length(to - from);
    auto phi   = atan2(z.z, z.x) + rotate.x;
    auto theta = acos(z.y) + rotate.y;
    theta      = clamp(theta, (T)0.001, (T)pi - (T)0.001);
    auto nz    = vec<T, 3>{sin(theta) * cos(phi) * lz, cos(theta) * lz,
        sin(theta) * sin(phi) * lz};
    from       = to - nz;
  }

  // dolly if necessary
  if (dolly != 0) {
    auto z  = normalize(to - from);
    auto lz = max(0.001f, length(to - from) * (1 + dolly));
    z *= lz;
    from = to - z;
  }

  // pan if necessary
  if (pan != vec<T, 2>{0, 0}) {
    auto z = normalize(to - from);
    auto x = normalize(cross(up, z));
    auto y = normalize(cross(z, x));
    auto t = vec<T, 3>{pan.x * x.x + pan.y * y.x, pan.x * x.y + pan.y * y.y,
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
  if (rotate != vec<T, 2>{0, 0}) {
    auto phi   = atan2(frame.z.z, frame.z.x) + rotate.x;
    auto theta = acos(frame.z.y) + rotate.y;
    theta      = clamp(theta, (T)0.001, (T)pi - (T)0.001);
    auto new_z = vec<T, 3>{
        sin(theta) * cos(phi), cos(theta), sin(theta) * sin(phi)};
    auto new_center = frame.o - frame.z * focus;
    auto new_o      = new_center + new_z * focus;
    frame           = lookat_frame(new_o, new_center, {0, 1, 0});
    focus           = length(new_o - new_center);
  }

  // pan if necessary
  if (dolly != 0) {
    auto c  = frame.o - frame.z * focus;
    focus   = max(focus * (1 + dolly), (T)0.001);
    frame.o = c + frame.z * focus;
  }

  // pan if necessary
  if (pan != vec<T, 2>{0, 0}) {
    frame.o += frame.x * pan.x + frame.y * pan.y;
  }
}

// FPS camera for UI navigation for a frame parametrization.
template <typename T>
inline void update_fpscam(
    frame<T, 3>& frame, const vec<T, 3>& transl, const vec<T, 2>& rotate) {
  // https://gamedev.stackexchange.com/questions/30644/how-to-keep-my-quaternion-using-fps-camera-from-tilting-and-messing-up
  auto y = vec<T, 3>{0, 1, 0};
  auto z = orthonormalize(frame.z, y);
  auto x = cross(y, z);

  auto rot = rotation_frame(vec<T, 3>{1, 0, 0}, rotate.y) *
             yocto::frame<T, 3>{frame.x, frame.y, frame.z, vec<T, 3>{0, 0, 0}} *
             rotation_frame(vec<T, 3>{0, 1, 0}, rotate.x);
  auto pos = frame.o + transl.x * x + transl.y * y + transl.z * z;

  frame = {rot.x, rot.y, rot.z, pos};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// PYTHON-LIKE ITERATORS
// -----------------------------------------------------------------------------
namespace yocto {

// Python `range()` equivalent. Construct an object to iterate over a sequence.
template <typename T>
constexpr auto range(T max) {
  return range((T)0, max);
}
template <typename T>
constexpr auto range(T min, T max) {
  struct range_iterator {
    T    index;
    void operator++() { ++index; }
    bool operator!=(const range_iterator& other) const {
      return index != other.index;
    }
    T operator*() const { return index; }
  };
  struct range_helper {
    T              begin_ = 0, end_ = 0;
    range_iterator begin() const { return {begin_}; }
    range_iterator end() const { return {end_}; }
  };
  return range_helper{min, max};
}
template <typename T>
constexpr auto range(T min, T max, T step) {
  struct range_iterator {
    T    index;
    T    step;
    void operator++() { index += step; }
    bool operator!=(const range_iterator& other) const {
      return index != other.index;
    }
    T operator*() const { return index; }
  };
  struct range_helper {
    T              begin_ = 0, end_ = 0, step_ = 0;
    range_iterator begin() const { return {begin_, step_}; }
    range_iterator end() const {
      return {begin_ + ((end_ - begin_) / step_) * step_, step_};
    }
  };
  return range_helper{min, max, step};
}

// Python enumerate
template <typename Sequence, typename T>
constexpr auto enumerate(const Sequence& sequence, T start) {
  using Iterator  = typename Sequence::const_iterator;
  using Reference = typename Sequence::const_reference;
  struct enumerate_iterator {
    T        index;
    Iterator iterator;
    bool     operator!=(const enumerate_iterator& other) const {
      return index != other.index;
    }
    void operator++() {
      ++index;
      ++iterator;
    }
    pair<const T&, Reference> operator*() const { return {index, *iterator}; }
  };
  struct enumerate_helper {
    const Sequence& sequence;
    T               begin_, end_;
    auto begin() { return enumerate_iterator{begin_, std::begin(sequence)}; }
    auto end() { return enumerate_iterator{end_, std::end(sequence)}; }
  };
  return enumerate_helper{sequence, 0, size(sequence)};
}

// Python enumerate
template <typename Sequence, typename T>
constexpr auto enumerate(Sequence& sequence, T start) {
  using Iterator  = typename Sequence::iterator;
  using Reference = typename Sequence::reference;
  struct enumerate_iterator {
    T        index;
    Iterator iterator;
    bool     operator!=(const enumerate_iterator& other) const {
      return index != other.index;
    }
    void operator++() {
      ++index;
      ++iterator;
    }
    pair<T&, Reference> operator*() const { return {index, *iterator}; }
  };
  struct enumerate_helper {
    Sequence& sequence;
    T         begin_, end_;
    auto begin() { return enumerate_iterator{begin_, std::begin(sequence)}; }
    auto end() { return enumerate_iterator{end_, std::end(sequence)}; }
  };
  return enumerate_helper{sequence, 0, size(sequence)};
}

// Python zip
template <typename Sequence1, typename Sequence2>
constexpr auto zip(const Sequence1& sequence1, const Sequence2& sequence2) {
  using Iterator1  = typename Sequence1::const_iterator;
  using Reference1 = typename Sequence1::const_reference;
  using Iterator2  = typename Sequence2::const_iterator;
  using Reference2 = typename Sequence2::const_reference;
  struct zip_iterator {
    Iterator1 iterator1;
    Iterator2 iterator2;
    bool      operator!=(const zip_iterator& other) const {
      return iterator1 != other.iterator1;
    }
    void operator++() {
      ++iterator1;
      ++iterator2;
    }
    pair<Reference1, Reference2> operator*() const {
      return {*iterator1, *iterator2};
    }
  };
  struct zip_helper {
    const Sequence1& sequence1;
    const Sequence2& sequence2;
    auto             begin() {
      return zip_iterator{std::begin(sequence1), std::begin(sequence2)};
    }
    auto end() {
      return zip_iterator{std::end(sequence1), std::end(sequence2)};
    }
  };
  return zip_helper{sequence1, sequence2};
}

// Python zip
template <typename Sequence1, typename Sequence2>
constexpr auto zip(Sequence1& sequence1, Sequence2& sequence2) {
  using Iterator1  = typename Sequence1::iterator;
  using Reference1 = typename Sequence1::reference;
  using Iterator2  = typename Sequence2::iterator;
  using Reference2 = typename Sequence2::reference;
  struct zip_iterator {
    Iterator1 iterator1;
    Iterator2 iterator2;
    bool      operator!=(const zip_iterator& other) const {
      return iterator1 != other.iterator1;
    }
    void operator++() {
      ++iterator1;
      ++iterator2;
    }
    pair<Reference1, Reference2> operator*() const {
      return {*iterator1, *iterator2};
    }
  };
  struct zip_helper {
    Sequence1& sequence1;
    Sequence2& sequence2;
    auto       begin() {
      return zip_iterator{std::begin(sequence1), std::begin(sequence2)};
    }
    auto end() {
      return zip_iterator{std::end(sequence1), std::end(sequence2)};
    }
  };
  return zip_helper{sequence1, sequence2};
}

// Python zip
template <typename Sequence1, typename Sequence2>
constexpr auto zip(const Sequence1& sequence1, Sequence2& sequence2) {
  using Iterator1  = typename Sequence1::const_iterator;
  using Reference1 = typename Sequence1::const_reference;
  using Iterator2  = typename Sequence2::iterator;
  using Reference2 = typename Sequence2::reference;
  struct zip_iterator {
    Iterator1 iterator1;
    Iterator2 iterator2;
    bool      operator!=(const zip_iterator& other) const {
      return iterator1 != other.iterator1;
    }
    void operator++() {
      ++iterator1;
      ++iterator2;
    }
    pair<Reference1, Reference2> operator*() const {
      return {*iterator1, *iterator2};
    }
  };
  struct zip_helper {
    const Sequence1& sequence1;
    Sequence2&       sequence2;
    auto             begin() {
      return zip_iterator{std::begin(sequence1), std::begin(sequence2)};
    }
    auto end() {
      return zip_iterator{std::end(sequence1), std::end(sequence2)};
    }
  };
  return zip_helper{sequence1, sequence2};
}

// Python zip
template <typename Sequence1, typename Sequence2>
constexpr auto zip(Sequence1& sequence1, const Sequence2& sequence2) {
  using Iterator1  = typename Sequence1::iterator;
  using Reference1 = typename Sequence1::reference;
  using Iterator2  = typename Sequence2::const_iterator;
  using Reference2 = typename Sequence2::const_reference;
  struct zip_iterator {
    Iterator1 iterator1;
    Iterator2 iterator2;
    bool      operator!=(const zip_iterator& other) const {
      return iterator1 != other.iterator1;
    }
    void operator++() {
      ++iterator1;
      ++iterator2;
    }
    pair<Reference1, Reference2> operator*() const {
      return {*iterator1, *iterator2};
    }
  };
  struct zip_helper {
    Sequence1&       sequence1;
    const Sequence2& sequence2;
    auto             begin() {
      return zip_iterator{std::begin(sequence1), std::begin(sequence2)};
    }
    auto end() {
      return zip_iterator{std::end(sequence1), std::end(sequence2)};
    }
  };
  return zip_helper{sequence1, sequence2};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SIGNED-SIZE
// -----------------------------------------------------------------------------
namespace yocto {

template <typename T>
inline std::ptrdiff_t ssize(const T& container) {
  return (std::ptrdiff_t)std::size(container);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// BACKWARD COMPATIBILITY
// -----------------------------------------------------------------------------
namespace yocto {

// Zero vector constants.
[[deprecated]] constexpr auto zero1f = vec1f{0};
[[deprecated]] constexpr auto zero2f = vec2f{0, 0};
[[deprecated]] constexpr auto zero3f = vec3f{0, 0, 0};
[[deprecated]] constexpr auto zero4f = vec4f{0, 0, 0, 0};
[[deprecated]] constexpr auto zero2i = vec2i{0, 0};
[[deprecated]] constexpr auto zero3i = vec3i{0, 0, 0};
[[deprecated]] constexpr auto zero4i = vec4i{0, 0, 0, 0};
[[deprecated]] constexpr auto zero4b = vec4b{0, 0, 0, 0};

// Indentity frames.
[[deprecated]] constexpr auto identity2x3f = frame2f{{1, 0}, {0, 1}, {0, 0}};
[[deprecated]] constexpr auto identity3x4f = frame3f{
    {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};

// Identity matrices constants.
[[deprecated]] constexpr auto identity2x2f = mat2f{{1, 0}, {0, 1}};
[[deprecated]] constexpr auto identity3x3f = mat3f{
    {1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
[[deprecated]] constexpr auto identity4x4f = mat4f{
    {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};

// Constants
[[deprecated]] constexpr auto identity_quat4f = quat4f{0, 0, 0, 1};

}  // namespace yocto

#endif
