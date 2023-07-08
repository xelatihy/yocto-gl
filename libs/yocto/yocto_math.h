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

#ifndef YOCTO_MATH_H_
#define YOCTO_MATH_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <array>
#include <cmath>
#include <cstdint>
#include <iosfwd>
#include <limits>
#include <type_traits>
#include <utility>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::pair;

}  // namespace yocto

// -----------------------------------------------------------------------------
// CUDA SUPPORT
// -----------------------------------------------------------------------------
#ifndef kernel
#ifdef __CUDACC__
#define kernel __device__
#else
#define kernel
#endif
#endif

// -----------------------------------------------------------------------------
// TYPE TRAITS
// -----------------------------------------------------------------------------
namespace yocto {

// Common type
template <typename... Ts>
using common_t = std::common_type_t<Ts...>;
#ifndef __CUDACC__
template <typename Func, typename... Ts>
using result_t = std::invoke_result_t<Func, Ts...>;
#else
template <typename Func, typename... Ts>
using result_t = std::result_of_t<Func(Ts...)>;
#endif
template <typename R>
using iterator_t = decltype(std::begin(std::declval<R&>()));
template <typename R>
using sentinel_t = decltype(std::end(std::declval<R&>()));
template <typename R>
using rsize_t = decltype(std::size(std::declval<R&>()));
#ifndef __CUDACC__
template <typename R>
using rdifference_t = std::iter_difference_t<iterator_t<R>>;
template <typename R>
using rvalue_t = std::iter_value_t<iterator_t<R>>;
template <typename R>
using rreference_t = std::iter_reference_t<iterator_t<R>>;
#endif

// Generic helpers
template <typename T>
struct type_constant {
  using type = T;
};

// Math specific traits
template <typename T>
struct is_vec : std::bool_constant<false> {};
template <typename T>
struct is_mat : std::bool_constant<false> {};
template <typename T>
struct is_frame : std::bool_constant<false> {};
template <typename T>
struct element_type : type_constant<void> {};
template <typename T>
struct element_dimension : std::integral_constant<size_t, 0> {};

// Helper variables
template <typename T>
constexpr auto vec_v = is_vec<std::remove_cvref_t<T>>::value;
template <typename T>
constexpr auto mat_v = is_mat<std::remove_cvref_t<T>>::value;
template <typename T>
constexpr auto frame_v = is_frame<std::remove_cvref_t<T>>::value;
#ifndef __CUDACC__
template <typename T>
using element_t = typename element_type<std::remove_cvref_t<T>>::type;
template <typename T>
constexpr auto element_d = element_dimension<std::remove_cvref_t<T>>::value;
#endif

}  // namespace yocto

// -----------------------------------------------------------------------------
// MATH CONSTANTS AND FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Basic types
using byte    = unsigned char;
using uint    = unsigned int;
using ushort  = unsigned short;
using index_t = ptrdiff_t;

// Numeric constants
template <typename T = double>
constexpr auto pi_t = (T)3.14159265358979323846;
constexpr auto pi   = 3.14159265358979323846;
constexpr auto pif  = (float)pi;

constexpr auto int_max  = std::numeric_limits<int>::max();
constexpr auto int_min  = std::numeric_limits<int>::lowest();
constexpr auto flt_max  = std::numeric_limits<float>::max();
constexpr auto flt_min  = std::numeric_limits<float>::lowest();
constexpr auto flt_eps  = std::numeric_limits<float>::epsilon();
constexpr auto uint_max = std::numeric_limits<unsigned int>::max();
constexpr auto uint_min = std::numeric_limits<unsigned int>::lowest();

template <typename T>
constexpr auto num_max = std::numeric_limits<T>::max();
template <typename T>
constexpr auto num_min = std::numeric_limits<T>::lowest();
template <typename T>
constexpr auto num_eps = std::numeric_limits<T>::epsilon();

// std imports
using std::data;
using std::size;
using std::swap;

template <typename T>
constexpr kernel T abs(T a) {
  return a < 0 ? -a : a;
}
template <typename T1, typename T2, typename T = common_t<T1, T2>>
constexpr kernel T min(T1 a, T2 b) {
  return (a < b) ? (T)a : (T)b;
}
template <typename T1, typename T2, typename T = common_t<T1, T2>>
constexpr kernel T max(T1 a, T2 b) {
  return (a > b) ? (T)a : (T)b;
}
template <typename T1, typename T2, typename T3,
    typename T = common_t<T1, T2, T3>>
constexpr kernel T clamp(T1 a, T2 min_, T3 max_) {
  return min(max(a, min_), max_);
}
template <typename T>
constexpr kernel T sign(T a) {
  return a < 0 ? (T)-1 : (T)1;
}
template <typename T>
constexpr kernel T sqr(T a) {
  return a * a;
}
template <typename T>
constexpr kernel T sqrt(T a) {
  return std::sqrt(a);
}
template <typename T>
constexpr kernel T sin(T a) {
  return std::sin(a);
}
template <typename T>
constexpr kernel T cos(T a) {
  return std::cos(a);
}
template <typename T>
constexpr kernel T tan(T a) {
  return std::tan(a);
}
template <typename T>
constexpr kernel T asin(T a) {
  return std::asin(a);
}
template <typename T>
constexpr kernel T acos(T a) {
  return std::acos(a);
}
template <typename T>
constexpr kernel T atan(T a) {
  return std::atan(a);
}
template <typename T>
constexpr kernel T log(T a) {
  return std::log(a);
}
template <typename T>
constexpr kernel T exp(T a) {
  return std::exp(a);
}
template <typename T>
constexpr kernel T log2(T a) {
  return std::log2(a);
}
template <typename T>
constexpr kernel T exp2(T a) {
  return std::exp2(a);
}
template <typename T1, typename T2, typename T = common_t<T1, T2>>
constexpr kernel T pow(T1 a, T2 b) {
  return std::pow((T)a, (T)b);
}
template <typename T>
constexpr kernel bool isfinite(T a) {
#ifndef __CUDACC__
  return std::isfinite(a);
#else
  return ::isfinite(a);
#endif
}
template <typename T1, typename T2, typename T = common_t<T1, T2>>
constexpr kernel T atan2(T1 a, T2 b) {
  return std::atan2((T)a, (T)b);
}
template <typename T>
constexpr kernel T round(T a) {
  return std::round(a);
}
template <typename T, typename T1>
constexpr kernel T fmod(T a, T1 b) {
  return std::fmod(a, (T)b);
}
template <typename T, typename T1>
constexpr kernel T mod(T a, T1 b) {
  if constexpr (std::is_floating_point_v<T> || std::is_floating_point_v<T1>) {
    auto m = fmod(a, b);
    return (m >= 0) ? m : m + b;
  } else {
    auto m = a % b;
    return (m >= 0) ? m : m + b;
  }
}
// template <typename T>
// inline void swap(T& a, T& b) {
//   std::swap(a, b);
// }
template <typename T>
constexpr kernel T radians(T a) {
  return a * (T)pi / 180;
}
template <typename T>
constexpr kernel T degrees(T a) {
  return a * 180 / (T)pi;
}
template <typename T1, typename T2, typename T3, typename T = common_t<T1, T2>>
constexpr kernel T lerp(T1 a, T2 b, T3 u) {
  return a * (1 - u) + b * u;
}
template <typename T1, typename T2, typename T = common_t<T1, T2>>
constexpr kernel T step(T1 a, T2 u) {
  return u < a ? (T)0 : (T)1;
}
template <typename T1, typename T2, typename T3,
    typename T = common_t<T1, T2, T3>>
constexpr kernel T smoothstep(T1 a, T2 b, T3 u) {
  auto t = clamp((u - a) / (b - a), (T)0, (T)1);
  return t * t * (3 - 2 * t);
}
template <typename T1, typename T2, typename T = common_t<T1, T2>>
constexpr kernel T bias(T1 a, T2 bias) {
  return a / ((1 / bias - 2) * (1 - a) + 1);
}
template <typename T1, typename T2, typename T = common_t<T1, T2>>
constexpr kernel T gain(T1 a, T2 gain) {
  return (a < (T)0.5) ? bias(a * 2, gain) / 2
                      : bias(a * 2 - 1, 1 - gain) / 2 + (T)0.5;
}
template <typename I>
constexpr kernel I pow2(I a) {
  return 1 << a;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// VECTORS
// -----------------------------------------------------------------------------
namespace yocto {

template <typename T, size_t N>
struct vec;

template <typename T>
struct vec<T, 1> {
  T x = 0;

  constexpr kernel vec() : x{0} {}
  constexpr kernel vec(T x_) : x{x_} {}
  template <typename T1>
  constexpr kernel vec(T1 x_) : x{(T)x_} {}

  constexpr vec(const vec& v)            = default;
  constexpr vec& operator=(const vec& v) = default;

  template <typename U>
  constexpr kernel explicit(!std::is_convertible_v<U, T>)
      vec(const vec<U, 1>& v) :
      x{(T)v.x} {}
  template <typename U>
  constexpr kernel explicit(!std::is_convertible_v<T, U>) operator vec<U, 1>() {
    return {(U)x};
  }

  constexpr kernel vec(const array<T, 1>& v) : x{v[0]} {}
  constexpr kernel operator array<T, 1>() { return {x}; }

  constexpr kernel T&       operator[](size_t i) { return (&x)[i]; }
  constexpr kernel const T& operator[](size_t i) const { return (&x)[i]; }
};

template <typename T>
struct vec<T, 2> {
  T x = 0;
  T y = 0;

  constexpr kernel vec() : x{0}, y{0} {}
  constexpr kernel explicit vec(T v) : x{v}, y{v} {}
  constexpr kernel vec(T x_, T y_) : x{x_}, y{y_} {}
  template <typename T1, typename T2>
  constexpr kernel vec(T1 x_, T2 y_) : x{(T)x_}, y{(T)y_} {}

  constexpr vec(const vec& v)            = default;
  constexpr vec& operator=(const vec& v) = default;

  template <typename U>
  constexpr kernel explicit(!std::is_convertible_v<U, T>)
      vec(const vec<U, 2>& v) :
      x{(T)v.x}, y{(T)v.y} {}
  template <typename U>
  constexpr kernel explicit(!std::is_convertible_v<T, U>) operator vec<U, 2>() {
    return {(U)x, (U)y};
  }

  constexpr kernel vec(const array<T, 2>& v) : x{v[0]}, y{v[1]} {}
  constexpr kernel operator array<T, 2>() { return {x, y}; }

  constexpr kernel T&       operator[](size_t i) { return (&x)[i]; }
  constexpr kernel const T& operator[](size_t i) const { return (&x)[i]; }
};

template <typename T>
struct vec<T, 3> {
  T x = 0;
  T y = 0;
  T z = 0;

  constexpr kernel vec() : x{0}, y{0}, z{0} {}
  constexpr kernel explicit vec(T v_) : x{v_}, y{v_}, z{v_} {}
  constexpr kernel vec(T x_, T y_, T z_) : x{x_}, y{y_}, z{z_} {}
  template <typename T1, typename T2, typename T3>
  constexpr kernel vec(T1 x_, T2 y_, T3 z_) : x{(T)x_}, y{(T)y_}, z{(T)z_} {}
  constexpr kernel vec(vec<T, 2> xy_, T z_) : x{xy_.x}, y{xy_.y}, z{z_} {}

  constexpr vec(const vec& v)            = default;
  constexpr vec& operator=(const vec& v) = default;

  template <typename U>
  constexpr kernel explicit(!std::is_convertible_v<U, T>)
      vec(const vec<U, 3>& v) :
      x{(T)v.x}, y{(T)v.y}, z{(T)v.z} {}
  template <typename U>
  constexpr kernel explicit(!std::is_convertible_v<T, U>) operator vec<U, 3>() {
    return {(U)x, (U)y, (U)z};
  }
  constexpr kernel vec(const array<T, 3>& v) : x{v[0]}, y{v[1]}, z{v[2]} {}
  constexpr kernel operator array<T, 3>() { return {x, y, z}; }

  constexpr kernel T&       operator[](size_t i) { return (&x)[i]; }
  constexpr kernel const T& operator[](size_t i) const { return (&x)[i]; }
};

template <typename T>
struct vec<T, 4> {
  T x = 0;
  T y = 0;
  T z = 0;
  T w = 0;

  constexpr kernel vec() : x{0}, y{0}, z{0}, w{0} {}
  constexpr kernel explicit vec(T v_) : x{v_}, y{v_}, z{v_}, w{v_} {}
  constexpr kernel vec(T x_, T y_, T z_, T w_) : x{x_}, y{y_}, z{z_}, w{w_} {}
  template <typename T1, typename T2, typename T3, typename T4>
  constexpr kernel vec(T1 x_, T2 y_, T3 z_, T4 w_) :
      x{(T)x_}, y{(T)y_}, z{(T)z_}, w{(T)w_} {}
  constexpr kernel vec(vec<T, 3> xyz_, T w_) :
      x{xyz_.x}, y{xyz_.y}, z{xyz_.z}, w{w_} {}

  constexpr vec(const vec& v)            = default;
  constexpr vec& operator=(const vec& v) = default;

  template <typename U>
  constexpr kernel explicit(!std::is_convertible_v<U, T>)
      vec(const vec<U, 4>& v) :
      x{(T)v.x}, y{(T)v.y}, z{(T)v.z}, w{(T)v.w} {}
  template <typename U>
  constexpr kernel explicit(!std::is_convertible_v<T, U>) operator vec<U, 4>() {
    return {(U)x, (U)y, (U)z, (U)w};
  }
  constexpr kernel vec(const array<T, 4>& v) :
      x{v[0]}, y{v[1]}, z{v[2]}, w{v[3]} {}
  constexpr kernel operator array<T, 4>() { return {x, y, z, w}; }

  constexpr kernel T&       operator[](size_t i) { return (&x)[i]; }
  constexpr kernel const T& operator[](size_t i) const { return (&x)[i]; }
};

// Deduction guides
#ifndef __CUDACC__
template <typename... Args>
vec(Args...) -> vec<common_t<Args...>, sizeof...(Args)>;
#endif

// Math specific traits
template <typename T, size_t N>
struct is_vec<vec<T, N>> : std::bool_constant<true> {};
template <typename T, size_t N>
struct element_type<vec<T, N>> : type_constant<T> {};
template <typename T, size_t N>
struct element_dimension<vec<T, N>> : std::integral_constant<size_t, N> {};

// Vector aliases
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
using vec1s = vec<size_t, 1>;
using vec2s = vec<size_t, 2>;
using vec3s = vec<size_t, 3>;
using vec4s = vec<size_t, 4>;
using vec1x = vec<index_t, 1>;
using vec2x = vec<index_t, 2>;
using vec3x = vec<index_t, 3>;
using vec4x = vec<index_t, 4>;

// Zero vector constants.
constexpr auto zero1f = vec1f{0};
constexpr auto zero2f = vec2f{0, 0};
constexpr auto zero3f = vec3f{0, 0, 0};
constexpr auto zero4f = vec4f{0, 0, 0, 0};
constexpr auto zero2i = vec2i{0, 0};
constexpr auto zero3i = vec3i{0, 0, 0};
constexpr auto zero4i = vec4i{0, 0, 0, 0};
constexpr auto zero4b = vec4b{0, 0, 0, 0};

// Generic constants
template <typename T, int N>
constexpr auto zero = vec<T, N>{0};
template <typename T, int N>
constexpr auto one = vec<T, N>{1};

// Element access
template <typename T>
constexpr kernel const vec<T, 2>& xy(const vec<T, 3>& a) {
  return (const vec<T, 2>&)a;
}
template <typename T>
constexpr kernel const vec<T, 2>& xy(const vec<T, 4>& a) {
  return (const vec<T, 2>&)a;
}
template <typename T>
constexpr kernel const vec<T, 3>& xyz(const vec<T, 4>& a) {
  return (const vec<T, 3>&)a;
}

// Vector sequence operations.
template <typename T, size_t N>
constexpr kernel bool empty(const vec<T, N>& a) {
  return false;
}
template <typename T, size_t N>
constexpr kernel size_t size(const vec<T, N>& a) {
  return N;
}
template <typename T, size_t N>
constexpr kernel ptrdiff_t ssize(const vec<T, N>& a) {
  return (ptrdiff_t)N;
}
template <typename T, size_t N>
constexpr kernel int isize(const vec<T, N>& a) {
  return (int)N;
}
template <typename T, size_t N>
constexpr kernel const T* begin(const vec<T, N>& a) {
  return &(a.x);
}
template <typename T, size_t N>
constexpr kernel const T* end(const vec<T, N>& a) {
  return &(a.x) + N;
}
template <typename T, size_t N>
constexpr kernel T* begin(vec<T, N>& a) {
  return &(a.x);
}
template <typename T, size_t N>
constexpr kernel T* end(vec<T, N>& a) {
  return &(a.x) + N;
}
template <typename T, size_t N>
constexpr kernel const T* data(const vec<T, N>& a) {
  return &(a.x);
}
template <typename T, size_t N>
constexpr kernel T* data(vec<T, N>& a) {
  return &(a.x);
}

// -----------------------------------------------------------------------------

// Vector comparison operations.
inline kernel bool operator==(const vec2f& a, const vec2f& b) {
  return a.x == b.x && a.y == b.y;
}
inline kernel bool operator==(const vec2f& a, float b) {
  return a.x == b && a.y == b;
}
inline kernel bool operator!=(const vec2f& a, const vec2f& b) {
  return a.x != b.x || a.y != b.y;
}
inline kernel bool operator!=(const vec2f& a, float b) {
  return a.x != b || a.y != b;
}
inline kernel bool operator<(const vec2f& a, const vec2f& b) {
  return a.x < b.x || (a.x == b.x && a.y < b.y);
}

// Vector operations.
inline kernel vec2f operator+(const vec2f& a) { return a; }
inline kernel vec2f operator-(const vec2f& a) { return {-a.x, -a.y}; }
inline kernel vec2f operator+(const vec2f& a, const vec2f& b) {
  return {a.x + b.x, a.y + b.y};
}
inline kernel vec2f operator+(const vec2f& a, float b) {
  return {a.x + b, a.y + b};
}
inline kernel vec2f operator+(float a, const vec2f& b) {
  return {a + b.x, a + b.y};
}
inline kernel vec2f operator-(const vec2f& a, const vec2f& b) {
  return {a.x - b.x, a.y - b.y};
}
inline kernel vec2f operator-(const vec2f& a, float b) {
  return {a.x - b, a.y - b};
}
inline kernel vec2f operator-(float a, const vec2f& b) {
  return {a - b.x, a - b.y};
}
inline kernel vec2f operator*(const vec2f& a, const vec2f& b) {
  return {a.x * b.x, a.y * b.y};
}
inline kernel vec2f operator*(const vec2f& a, float b) {
  return {a.x * b, a.y * b};
}
inline kernel vec2f operator*(float a, const vec2f& b) {
  return {a * b.x, a * b.y};
}
inline kernel vec2f operator/(const vec2f& a, const vec2f& b) {
  return {a.x / b.x, a.y / b.y};
}
inline kernel vec2f operator/(const vec2f& a, float b) {
  return {a.x / b, a.y / b};
}
inline kernel vec2f operator/(float a, const vec2f& b) {
  return {a / b.x, a / b.y};
}

// Vector assignments
inline kernel vec2f& operator+=(vec2f& a, const vec2f& b) { return a = a + b; }
inline kernel vec2f& operator+=(vec2f& a, float b) { return a = a + b; }
inline kernel vec2f& operator-=(vec2f& a, const vec2f& b) { return a = a - b; }
inline kernel vec2f& operator-=(vec2f& a, float b) { return a = a - b; }
inline kernel vec2f& operator*=(vec2f& a, const vec2f& b) { return a = a * b; }
inline kernel vec2f& operator*=(vec2f& a, float b) { return a = a * b; }
inline kernel vec2f& operator/=(vec2f& a, const vec2f& b) { return a = a / b; }
inline kernel vec2f& operator/=(vec2f& a, float b) { return a = a / b; }

// Vector products and lengths.
inline kernel float dot(const vec2f& a, const vec2f& b) {
  return a.x * b.x + a.y * b.y;
}
inline kernel float cross(const vec2f& a, const vec2f& b) {
  return a.x * b.y - a.y * b.x;
}

inline kernel float length(const vec2f& a) { return sqrt(dot(a, a)); }
inline kernel float length2(const vec2f& a) { return dot(a, a); }
[[deprecated]] inline kernel float length_squared(const vec2f& a) {
  return dot(a, a);
}
inline kernel vec2f normalize(const vec2f& a) {
  auto l = length(a);
  return (l != 0) ? a / l : a;
}
inline kernel float distance(const vec2f& a, const vec2f& b) {
  return length(a - b);
}
inline kernel float distance2(const vec2f& a, const vec2f& b) {
  return dot(a - b, a - b);
}
[[deprecated]] inline kernel float distance_squared(
    const vec2f& a, const vec2f& b) {
  return dot(a - b, a - b);
}
inline kernel float angle(const vec2f& a, const vec2f& b) {
  return acos(clamp(dot(normalize(a), normalize(b)), -1.0f, 1.0f));
}

// Max element and clamp.
inline kernel vec2f max(const vec2f& a, const vec2f& b) {
  return {max(a.x, b.x), max(a.y, b.y)};
}
inline kernel vec2f max(const vec2f& a, float b) {
  return {max(a.x, b), max(a.y, b)};
}
inline kernel vec2f min(const vec2f& a, const vec2f& b) {
  return {min(a.x, b.x), min(a.y, b.y)};
}
inline kernel vec2f min(const vec2f& a, float b) {
  return {min(a.x, b), min(a.y, b)};
}
constexpr kernel vec2f clamp(
    const vec2f& x, const vec2f& min, const vec2f& max) {
  return {clamp(x.x, min.x, max.x), clamp(x.y, min.y, max.y)};
}
constexpr kernel vec2f clamp(const vec2f& x, float min, float max) {
  return {clamp(x.x, min, max), clamp(x.y, min, max)};
}
constexpr kernel vec2f clamp(const vec2f& x, const vec2f& min, float max) {
  return {clamp(x.x, min.x, max), clamp(x.y, min.y, max)};
}
constexpr kernel vec2f clamp(const vec2f& x, float min, const vec2f& max) {
  return {clamp(x.x, min, max.x), clamp(x.y, min, max.y)};
}
constexpr kernel vec2f lerp(const vec2f& a, const vec2f& b, float u) {
  return a * (1 - u) + b * u;
}
constexpr kernel vec2f lerp(const vec2f& a, const vec2f& b, const vec2f& u) {
  return a * (1 - u) + b * u;
}

inline kernel float max(const vec2f& a) { return max(a.x, a.y); }
inline kernel float min(const vec2f& a) { return min(a.x, a.y); }
inline kernel int   argmax(const vec2f& a) { return a.x >= a.y ? 0 : 1; }
inline kernel int   argmin(const vec2f& a) { return a.x <= a.y ? 0 : 1; }
inline kernel int   find_index(const vec2f& a, float b) {
  return a.x == b ? 0 : (a.y == b ? 1 : -1);
}

inline kernel float sum(const vec2f& a) { return a.x + a.y; }
inline kernel float prod(const vec2f& a) { return a.x * a.y; }
inline kernel float mean(const vec2f& a) { return sum(a) / 2; }

// Functions applied to vector elements
inline kernel vec2f abs(const vec2f& a) { return {abs(a.x), abs(a.y)}; }
inline kernel vec2f sqr(const vec2f& a) { return {sqr(a.x), sqr(a.y)}; }
inline kernel vec2f sqrt(const vec2f& a) { return {sqrt(a.x), sqrt(a.y)}; }
inline kernel vec2f exp(const vec2f& a) { return {exp(a.x), exp(a.y)}; }
inline kernel vec2f log(const vec2f& a) { return {log(a.x), log(a.y)}; }
inline kernel vec2f exp2(const vec2f& a) { return {exp2(a.x), exp2(a.y)}; }
inline kernel vec2f log2(const vec2f& a) { return {log2(a.x), log2(a.y)}; }
inline kernel vec2f pow(const vec2f& a, const vec2f& b) {
  return {pow(a.x, b.x), pow(a.y, b.y)};
}
inline kernel vec2f pow(const vec2f& a, float b) {
  return {pow(a.x, b), pow(a.y, b)};
}
inline kernel vec2f round(const vec2f& a) { return {round(a.x), round(a.y)}; }
inline kernel vec2f fmod(const vec2f& a, const vec2f& b) {
  return {fmod(a.x, b.x), fmod(a.y, b.y)};
}
inline kernel vec2f fmod(const vec2f& a, float b) {
  return {fmod(a.x, b), fmod(a.y, b)};
}
inline kernel vec2f mod(const vec2f& a, const vec2f& b) {
  return {mod(a.x, b.x), mod(a.y, b.y)};
}
inline kernel vec2f mod(const vec2f& a, float b) {
  return {mod(a.x, b), mod(a.y, b)};
}
inline kernel vec2f bias(const vec2f& a, float b) {
  return {bias(a.x, b), bias(a.y, b)};
}
inline kernel vec2f gain(const vec2f& a, float b) {
  return {gain(a.x, b), gain(a.y, b)};
}
inline kernel bool isfinite(const vec2f& a) {
  return isfinite(a.x) && isfinite(a.y);
}

// Conversion between ranges
inline kernel vec2f unit_to_uv(const vec2f& a) { return (a + 1) / 2; }
inline kernel vec2f uv_to_unit(const vec2f& a) { return a * 2 - 1; }

// -----------------------------------------------------------------------------

// Vector comparison operations.
inline kernel bool operator==(const vec3f& a, const vec3f& b) {
  return a.x == b.x && a.y == b.y && a.z == b.z;
}
inline kernel bool operator==(const vec3f& a, float b) {
  return a.x == b && a.y == b && a.z == b;
}
inline kernel bool operator!=(const vec3f& a, const vec3f& b) {
  return a.x != b.x || a.y != b.y || a.z != b.z;
}
inline kernel bool operator!=(const vec3f& a, float b) {
  return a.x != b || a.y != b || a.z != b;
}
inline kernel bool operator<(const vec3f& a, const vec3f& b) {
  return a.x < b.x || (a.x == b.x && (a.y < b.y || (a.y == b.y && a.z < b.z)));
}

// Vector operations.
inline kernel vec3f operator+(const vec3f& a) { return a; }
inline kernel vec3f operator-(const vec3f& a) { return {-a.x, -a.y, -a.z}; }
inline kernel vec3f operator+(const vec3f& a, const vec3f& b) {
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}
inline kernel vec3f operator+(const vec3f& a, float b) {
  return {a.x + b, a.y + b, a.z + b};
}
inline kernel vec3f operator+(float a, const vec3f& b) {
  return {a + b.x, a + b.y, a + b.z};
}
inline kernel vec3f operator-(const vec3f& a, const vec3f& b) {
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}
inline kernel vec3f operator-(const vec3f& a, float b) {
  return {a.x - b, a.y - b, a.z - b};
}
inline kernel vec3f operator-(float a, const vec3f& b) {
  return {a - b.x, a - b.y, a - b.z};
}
inline kernel vec3f operator*(const vec3f& a, const vec3f& b) {
  return {a.x * b.x, a.y * b.y, a.z * b.z};
}
inline kernel vec3f operator*(const vec3f& a, float b) {
  return {a.x * b, a.y * b, a.z * b};
}
inline kernel vec3f operator*(float a, const vec3f& b) {
  return {a * b.x, a * b.y, a * b.z};
}
inline kernel vec3f operator/(const vec3f& a, const vec3f& b) {
  return {a.x / b.x, a.y / b.y, a.z / b.z};
}
inline kernel vec3f operator/(const vec3f& a, float b) {
  return {a.x / b, a.y / b, a.z / b};
}
inline kernel vec3f operator/(float a, const vec3f& b) {
  return {a / b.x, a / b.y, a / b.z};
}

// Vector assignments
inline kernel vec3f& operator+=(vec3f& a, const vec3f& b) { return a = a + b; }
inline kernel vec3f& operator+=(vec3f& a, float b) { return a = a + b; }
inline kernel vec3f& operator-=(vec3f& a, const vec3f& b) { return a = a - b; }
inline kernel vec3f& operator-=(vec3f& a, float b) { return a = a - b; }
inline kernel vec3f& operator*=(vec3f& a, const vec3f& b) { return a = a * b; }
inline kernel vec3f& operator*=(vec3f& a, float b) { return a = a * b; }
inline kernel vec3f& operator/=(vec3f& a, const vec3f& b) { return a = a / b; }
inline kernel vec3f& operator/=(vec3f& a, float b) { return a = a / b; }

// Vector products and lengths.
inline kernel float dot(const vec3f& a, const vec3f& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline kernel vec3f cross(const vec3f& a, const vec3f& b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline kernel float length(const vec3f& a) { return sqrt(dot(a, a)); }
inline kernel float length2(const vec3f& a) { return dot(a, a); }
[[deprecated]] inline kernel float length_squared(const vec3f& a) {
  return dot(a, a);
}
inline kernel vec3f normalize(const vec3f& a) {
  auto l = length(a);
  return (l != 0) ? a / l : a;
}
inline kernel float distance(const vec3f& a, const vec3f& b) {
  return length(a - b);
}
inline kernel float distance2(const vec3f& a, const vec3f& b) {
  return dot(a - b, a - b);
}
[[deprecated]] inline kernel float distance_squared(
    const vec3f& a, const vec3f& b) {
  return dot(a - b, a - b);
}
inline kernel float angle(const vec3f& a, const vec3f& b) {
  return acos(clamp(dot(normalize(a), normalize(b)), -1.0f, 1.0f));
}
inline kernel vec3f slerp(const vec3f& a, const vec3f& b, float u) {
  // https://en.wikipedia.org/wiki/Slerp
  auto an = normalize(a), bn = normalize(b);
  auto d = dot(an, bn);
  if (d < 0) {
    bn = -bn;
    d  = -d;
  }
  if (d > (float)0.9995) return normalize(an + u * (bn - an));
  auto theta = acos(clamp(d, (float)-1, (float)1));
  if (theta == 0) return an;
  return an * (sin(theta * (1 - u)) / sin(theta)) +
         bn * (sin(theta * u) / sin(theta));
}

// Orthogonal vectors.
inline kernel vec3f orthogonal(const vec3f& v) {
  // http://lolengine.net/blog/2013/09/21/picking-orthogonal-vector-combing-coconuts)
  return abs(v.x) > abs(v.z) ? vec3f{-v.y, v.x, 0} : vec3f{0, -v.z, v.y};
}
inline kernel vec3f orthonormalize(const vec3f& a, const vec3f& b) {
  return normalize(a - b * dot(a, b));
}

// Reflected and refracted vector.
inline kernel vec3f reflect(const vec3f& w, const vec3f& n) {
  return -w + 2 * dot(n, w) * n;
}
inline kernel vec3f refract(const vec3f& w, const vec3f& n, float inv_eta) {
  auto cosine = dot(n, w);
  auto k      = 1 + inv_eta * inv_eta * (cosine * cosine - 1);
  if (k < 0) return {0, 0, 0};  // tir
  return -w * inv_eta + (inv_eta * cosine - sqrt(k)) * n;
}

// Max element and clamp.
inline kernel vec3f max(const vec3f& a, const vec3f& b) {
  return {max(a.x, b.x), max(a.y, b.y), max(a.z, b.z)};
}
inline kernel vec3f max(const vec3f& a, float b) {
  return {max(a.x, b), max(a.y, b), max(a.z, b)};
}
inline kernel vec3f min(const vec3f& a, const vec3f& b) {
  return {min(a.x, b.x), min(a.y, b.y), min(a.z, b.z)};
}
inline kernel vec3f min(const vec3f& a, float b) {
  return {min(a.x, b), min(a.y, b), min(a.z, b)};
}
inline kernel vec3f clamp(const vec3f& x, const vec3f& min, const vec3f& max) {
  return {clamp(x.x, min.x, max.x), clamp(x.y, min.y, max.y),
      clamp(x.z, min.z, max.z)};
}
inline kernel vec3f clamp(const vec3f& x, float min, float max) {
  return {clamp(x.x, min, max), clamp(x.y, min, max), clamp(x.z, min, max)};
}
inline kernel vec3f clamp(const vec3f& x, const vec3f& min, float max) {
  return {
      clamp(x.x, min.x, max), clamp(x.y, min.y, max), clamp(x.z, min.z, max)};
}
inline kernel vec3f clamp(const vec3f& x, float min, const vec3f& max) {
  return {
      clamp(x.x, min, max.x), clamp(x.y, min, max.y), clamp(x.z, min, max.z)};
}
inline kernel vec3f lerp(const vec3f& a, const vec3f& b, float u) {
  return a * (1 - u) + b * u;
}
inline kernel vec3f lerp(const vec3f& a, const vec3f& b, const vec3f& u) {
  return a * (1 - u) + b * u;
}

inline kernel float max(const vec3f& a) { return max(max(a.x, a.y), a.z); }
inline kernel float min(const vec3f& a) { return min(min(a.x, a.y), a.z); }
inline kernel int   argmax(const vec3f& a) {
  return a.x >= a.y ? (a.x >= a.z ? 0 : 2) : (a.y >= a.z ? 1 : 2);
}
inline kernel int argmin(const vec3f& a) {
  return a.x <= a.y ? (a.x <= a.z ? 0 : 2) : (a.y <= a.z ? 1 : 2);
}
inline kernel int find_index(const vec3f& a, float b) {
  return a.x == b ? 0 : (a.y == b ? 1 : (a.z == b ? 2 : -1));
}

inline kernel float sum(const vec3f& a) { return a.x + a.y + a.z; }
inline kernel float prod(const vec3f& a) { return a.x * a.y * a.z; }
inline kernel float mean(const vec3f& a) { return sum(a) / 3; }

// Functions applied to vector elements
inline kernel vec3f abs(const vec3f& a) {
  return {abs(a.x), abs(a.y), abs(a.z)};
}
inline kernel vec3f sqr(const vec3f& a) {
  return {sqr(a.x), sqr(a.y), sqr(a.z)};
}
inline kernel vec3f sqrt(const vec3f& a) {
  return {sqrt(a.x), sqrt(a.y), sqrt(a.z)};
}
inline kernel vec3f exp(const vec3f& a) {
  return {exp(a.x), exp(a.y), exp(a.z)};
}
inline kernel vec3f log(const vec3f& a) {
  return {log(a.x), log(a.y), log(a.z)};
}
inline kernel vec3f exp2(const vec3f& a) {
  return {exp2(a.x), exp2(a.y), exp2(a.z)};
}
inline kernel vec3f log2(const vec3f& a) {
  return {log2(a.x), log2(a.y), log2(a.z)};
}
inline kernel vec3f pow(const vec3f& a, const vec3f& b) {
  return {pow(a.x, b.x), pow(a.y, b.y), pow(a.z, b.z)};
}
inline kernel vec3f pow(const vec3f& a, float b) {
  return {pow(a.x, b), pow(a.y, b), pow(a.z, b)};
}
inline kernel vec3f round(const vec3f& a) {
  return {round(a.x), round(a.y), round(a.z)};
}
inline kernel vec3f fmod(const vec3f& a, const vec3f& b) {
  return {fmod(a.x, b.x), fmod(a.y, b.y), fmod(a.z, b.z)};
}
inline kernel vec3f fmod(const vec3f& a, float b) {
  return {fmod(a.x, b), fmod(a.y, b), fmod(a.z, b)};
}
inline kernel vec3f mod(const vec3f& a, const vec3f& b) {
  return {mod(a.x, b.x), mod(a.y, b.y), mod(a.z, b.z)};
}
inline kernel vec3f mod(const vec3f& a, float b) {
  return {mod(a.x, b), mod(a.y, b), mod(a.z, b)};
}
inline kernel vec3f bias(const vec3f& a, float b) {
  return {bias(a.x, b), bias(a.y, b), bias(a.z, b)};
}
inline kernel vec3f gain(const vec3f& a, float b) {
  return {gain(a.x, b), gain(a.y, b), gain(a.z, b)};
}
inline kernel bool isfinite(const vec3f& a) {
  return isfinite(a.x) && isfinite(a.y) && isfinite(a.z);
}

// Conversion between coordinates
inline kernel vec2f cartesian_to_sphericaluv(const vec3f& w) {
  auto phi = atan2(w.y, w.x), theta = acos(clamp(w.z, -1, 1));
  return {mod(phi / (2 * (float)pi), 1), theta / (float)pi};
}
inline kernel vec2f cartesiany_to_sphericaluv(const vec3f& w) {
  auto phi = atan2(w.z, w.x), theta = acos(clamp(w.y, -1, 1));
  return {mod(phi / (2 * (float)pi), 1), theta / (float)pi};
}
inline kernel vec3f sphericaluv_to_cartesian(const vec2f& uv) {
  auto phi = uv.x * 2 * (float)pi, theta = uv.y * (float)pi;
  return {cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)};
}
inline kernel vec3f sphericaluv_to_cartesiany(const vec2f& uv) {
  auto phi = uv.x * 2 * (float)pi, theta = uv.y * (float)pi;
  return {cos(phi) * sin(theta), cos(theta), sin(phi) * sin(theta)};
}

// ------------------------------------------------------------------

// Vector comparison operations.
inline kernel bool operator==(const vec4f& a, const vec4f& b) {
  return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}
inline kernel bool operator==(const vec4f& a, float b) {
  return a.x == b && a.y == b && a.z == b && a.w == b;
}
inline kernel bool operator!=(const vec4f& a, const vec4f& b) {
  return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w;
}
inline kernel bool operator!=(const vec4f& a, float b) {
  return a.x != b || a.y != b || a.z != b || a.w != b;
}
inline kernel bool operator<(const vec4f& a, const vec4f& b) {
  return a.x < b.x ||
         (a.x == b.x &&
             (a.y < b.y ||
                 (a.y == b.y && (a.z < b.z || (a.z == b.z && a.w < b.w)))));
}

// Vector operations.
inline kernel vec4f operator+(const vec4f& a) { return a; }
inline kernel vec4f operator-(const vec4f& a) {
  return {-a.x, -a.y, -a.z, -a.w};
}
inline kernel vec4f operator+(const vec4f& a, const vec4f& b) {
  return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}
inline kernel vec4f operator+(const vec4f& a, float b) {
  return {a.x + b, a.y + b, a.z + b, a.w + b};
}
inline kernel vec4f operator+(float a, const vec4f& b) {
  return {a + b.x, a + b.y, a + b.z, a + b.w};
}
inline kernel vec4f operator-(const vec4f& a, const vec4f& b) {
  return {a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w};
}
inline kernel vec4f operator-(const vec4f& a, float b) {
  return {a.x - b, a.y - b, a.z - b, a.w - b};
}
inline kernel vec4f operator-(float a, const vec4f& b) {
  return {a - b.x, a - b.y, a - b.z, a - b.w};
}
inline kernel vec4f operator*(const vec4f& a, const vec4f& b) {
  return {a.x * b.x, a.y * b.y, a.z * b.z, a.w * b.w};
}
inline kernel vec4f operator*(const vec4f& a, float b) {
  return {a.x * b, a.y * b, a.z * b, a.w * b};
}
inline kernel vec4f operator*(float a, const vec4f& b) {
  return {a * b.x, a * b.y, a * b.z, a * b.w};
}
inline kernel vec4f operator/(const vec4f& a, const vec4f& b) {
  return {a.x / b.x, a.y / b.y, a.z / b.z, a.w / b.w};
}
inline kernel vec4f operator/(const vec4f& a, float b) {
  return {a.x / b, a.y / b, a.z / b, a.w / b};
}
inline kernel vec4f operator/(float a, const vec4f& b) {
  return {a / b.x, a / b.y, a / b.z, a / b.w};
}

// Vector assignments
inline kernel vec4f& operator+=(vec4f& a, const vec4f& b) { return a = a + b; }
inline kernel vec4f& operator+=(vec4f& a, float b) { return a = a + b; }
inline kernel vec4f& operator-=(vec4f& a, const vec4f& b) { return a = a - b; }
inline kernel vec4f& operator-=(vec4f& a, float b) { return a = a - b; }
inline kernel vec4f& operator*=(vec4f& a, const vec4f& b) { return a = a * b; }
inline kernel vec4f& operator*=(vec4f& a, float b) { return a = a * b; }
inline kernel vec4f& operator/=(vec4f& a, const vec4f& b) { return a = a / b; }
inline kernel vec4f& operator/=(vec4f& a, float b) { return a = a / b; }

// Vector products and lengths.
inline kernel float dot(const vec4f& a, const vec4f& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}
inline kernel float length(const vec4f& a) { return sqrt(dot(a, a)); }
inline kernel float length2(const vec4f& a) { return dot(a, a); }
[[deprecated]] inline kernel float length_squared(const vec4f& a) {
  return dot(a, a);
}
inline kernel vec4f normalize(const vec4f& a) {
  auto l = length(a);
  return (l != 0) ? a / l : a;
}
inline kernel float distance(const vec4f& a, const vec4f& b) {
  return length(a - b);
}
inline kernel float distance2(const vec4f& a, const vec4f& b) {
  return dot(a - b, a - b);
}
[[deprecated]] inline kernel float distance_squared(
    const vec4f& a, const vec4f& b) {
  return dot(a - b, a - b);
}
inline kernel float angle(const vec4f& a, const vec4f& b) {
  return acos(clamp(dot(normalize(a), normalize(b)), (float)-1, (float)1));
}
inline kernel vec4f slerp(const vec4f& a, const vec4f& b, float u) {
  // https://en.wikipedia.org/wiki/Slerp
  auto an = normalize(a), bn = normalize(b);
  auto d = dot(an, bn);
  if (d < 0) {
    bn = -bn;
    d  = -d;
  }
  if (d > (float)0.9995) return normalize(an + u * (bn - an));
  auto theta = acos(clamp(d, (float)-1, (float)1));
  if (theta == 0) return an;
  return an * (sin(theta * (1 - u)) / sin(theta)) +
         bn * (sin(theta * u) / sin(theta));
}

// Max element and clamp.
inline kernel vec4f max(const vec4f& a, const vec4f& b) {
  return {max(a.x, b.x), max(a.y, b.y), max(a.z, b.z), max(a.w, b.w)};
}
inline kernel vec4f max(const vec4f& a, float b) {
  return {max(a.x, b), max(a.y, b), max(a.z, b), max(a.w, b)};
}
inline kernel vec4f min(const vec4f& a, const vec4f& b) {
  return {min(a.x, b.x), min(a.y, b.y), min(a.z, b.z), min(a.w, b.w)};
}
inline kernel vec4f min(const vec4f& a, float b) {
  return {min(a.x, b), min(a.y, b), min(a.z, b), min(a.w, b)};
}
inline kernel vec4f clamp(const vec4f& x, const vec4f& min, const vec4f& max) {
  return {clamp(x.x, min.x, max.x), clamp(x.y, min.y, max.y),
      clamp(x.z, min.z, max.z), clamp(x.w, min.w, max.w)};
}
inline kernel vec4f clamp(const vec4f& x, float min, float max) {
  return {clamp(x.x, min, max), clamp(x.y, min, max), clamp(x.z, min, max),
      clamp(x.w, min, max)};
}
inline kernel vec4f clamp(const vec4f& x, const vec4f& min, float max) {
  return {clamp(x.x, min.x, max), clamp(x.y, min.y, max),
      clamp(x.z, min.z, max), clamp(x.w, min.w, max)};
}
inline kernel vec4f clamp(const vec4f& x, float min, const vec4f& max) {
  return {clamp(x.x, min, max.x), clamp(x.y, min, max.y),
      clamp(x.z, min, max.z), clamp(x.w, min, max.w)};
}
inline kernel vec4f lerp(const vec4f& a, const vec4f& b, float u) {
  return a * (1 - u) + b * u;
}
inline kernel vec4f lerp(const vec4f& a, const vec4f& b, const vec4f& u) {
  return a * (1 - u) + b * u;
}

inline kernel float max(const vec4f& a) {
  return max(max(max(a.x, a.y), a.z), a.w);
}
inline kernel float min(const vec4f& a) {
  return min(min(min(a.x, a.y), a.z), a.w);
}
inline kernel int argmax(const vec4f& a) {
  if (a.w >= a.x && a.w >= a.y && a.w >= a.z)
    return 3;
  else
    return argmax(xyz(a));
}
inline kernel int argmin(const vec4f& a) {
  if (a.w <= a.x && a.w <= a.y && a.w <= a.z)
    return 3;
  else
    return argmin(xyz(a));
}
inline kernel int find_index(const vec4f& a, float b) {
  return a.x == b ? 0 : (a.y == b ? 1 : (a.z == b ? 2 : (a.w == b ? 3 : -1)));
}

inline kernel float sum(const vec4f& a) { return a.x + a.y + a.z + a.w; }
inline kernel float prod(const vec4f& a) { return a.x * a.y * a.z * a.w; }
inline kernel float mean(const vec4f& a) { return sum(a) / 4; }

// Functions applied to vector elements
inline kernel vec4f abs(const vec4f& a) {
  return {abs(a.x), abs(a.y), abs(a.z), abs(a.w)};
}
inline kernel vec4f sqr(const vec4f& a) {
  return {sqr(a.x), sqr(a.y), sqr(a.z), sqr(a.w)};
}
inline kernel vec4f sqrt(const vec4f& a) {
  return {sqrt(a.x), sqrt(a.y), sqrt(a.z), sqrt(a.w)};
}
inline kernel vec4f exp(const vec4f& a) {
  return {exp(a.x), exp(a.y), exp(a.z), exp(a.w)};
}
inline kernel vec4f log(const vec4f& a) {
  return {log(a.x), log(a.y), log(a.z), log(a.w)};
}
inline kernel vec4f exp2(const vec4f& a) {
  return {exp2(a.x), exp2(a.y), exp2(a.z), exp2(a.w)};
}
inline kernel vec4f log2(const vec4f& a) {
  return {log2(a.x), log2(a.y), log2(a.z), log2(a.w)};
}
inline kernel vec4f pow(const vec4f& a, const vec4f& b) {
  return {pow(a.x, b.x), pow(a.y, b.y), pow(a.z, b.z), pow(a.w, b.w)};
}
inline kernel vec4f pow(const vec4f& a, float b) {
  return {pow(a.x, b), pow(a.y, b), pow(a.z, b), pow(a.w, b)};
}
inline kernel vec4f round(const vec4f& a) {
  return {round(a.x), round(a.y), round(a.z), round(a.w)};
}
inline kernel vec4f fmod(const vec4f& a, const vec4f& b) {
  return {fmod(a.x, b.x), fmod(a.y, b.y), fmod(a.z, b.z), fmod(a.w, b.w)};
}
inline kernel vec4f fmod(const vec4f& a, float b) {
  return {fmod(a.x, b), fmod(a.y, b), fmod(a.z, b), fmod(a.w, b)};
}
inline kernel vec4f mod(const vec4f& a, const vec4f& b) {
  return {mod(a.x, b.x), mod(a.y, b.y), mod(a.z, b.z), mod(a.w, b.w)};
}
inline kernel vec4f mod(const vec4f& a, float b) {
  return {mod(a.x, b), mod(a.y, b), mod(a.z, b), mod(a.w, b)};
}
inline kernel vec4f bias(const vec4f& a, float b) {
  return {bias(a.x, b), bias(a.y, b), bias(a.z, b), bias(a.w, b)};
}
inline kernel vec4f gain(const vec4f& a, float b) {
  return {gain(a.x, b), gain(a.y, b), gain(a.z, b), gain(a.w, b)};
}
inline kernel bool isfinite(const vec4f& a) {
  return isfinite(a.x) && isfinite(a.y) && isfinite(a.z) && isfinite(a.w);
}

// Quaternion operations represented as xi + yj + zk + w
// const auto identity_quat4f = vec4f{0, 0, 0, 1};
inline kernel vec4f quat_mul(const vec4f& a, float b) {
  return {a.x * b, a.y * b, a.z * b, a.w * b};
}
inline kernel vec4f quat_mul(const vec4f& a, const vec4f& b) {
  return {a.x * b.w + a.w * b.x + a.y * b.w - a.z * b.y,
      a.y * b.w + a.w * b.y + a.z * b.x - a.x * b.z,
      a.z * b.w + a.w * b.z + a.x * b.y - a.y * b.x,
      a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z};
}
inline kernel vec4f quat_conjugate(const vec4f& a) {
  return {-a.x, -a.y, -a.z, a.w};
}
inline kernel vec4f quat_inverse(const vec4f& a) {
  return quat_conjugate(a) / dot(a, a);
}

// ------------------------------------------------------------------

// Vector comparison operations.
template <typename T1, typename T2, size_t N>
constexpr kernel bool operator==(const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return a.x == b.x;
  } else if constexpr (N == 2) {
    return a.x == b.x && a.y == b.y;
  } else if constexpr (N == 3) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
  } else if constexpr (N == 4) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
  } else {
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel bool operator==(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return a.x == b;
  } else if constexpr (N == 2) {
    return a.x == b && a.y == b;
  } else if constexpr (N == 3) {
    return a.x == b && a.y == b && a.z == b;
  } else if constexpr (N == 4) {
    return a.x == b && a.y == b && a.z == b && a.w == b;
  } else {
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel bool operator!=(const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N != 1) {
    return a.x != b.x;
  } else if constexpr (N != 2) {
    return a.x != b.x || a.y != b.y;
  } else if constexpr (N != 3) {
    return a.x != b.x || a.y != b.y || a.z != b.z;
  } else if constexpr (N != 4) {
    return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w;
  } else {
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel bool operator!=(const vec<T1, N>& a, T2 b) {
  if constexpr (N != 1) {
    return a.x != b;
  } else if constexpr (N != 2) {
    return a.x != b || a.y != b;
  } else if constexpr (N != 3) {
    return a.x != b || a.y != b || a.z != b;
  } else if constexpr (N != 4) {
    return a.x != b || a.y != b || a.z != b || a.w != b;
  } else {
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel bool operator<(const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return a.x < b.x;
  } else if constexpr (N == 2) {
    return a.x < b.x || (a.x == b.x && a.y < b.y);
  } else if constexpr (N == 3) {
    return a.x < b.x ||
           (a.x == b.x && (a.y < b.y || (a.y == b.y && a.z < b.z)));
  } else if constexpr (N == 4) {
    return a.x < b.x ||
           (a.x == b.x &&
               (a.y < b.y ||
                   (a.y == b.y && (a.z < b.z || (a.z == b.z && a.w < b.w)))));
  } else {
  }
}

// Vector operations.
template <typename T, size_t N>
constexpr kernel vec<T, N> operator+(const vec<T, N>& a) {
  return a;
}
template <typename T, size_t N>
constexpr kernel vec<T, N> operator-(const vec<T, N>& a) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator+(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator+(const vec<T1, N>& a, T2 b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator+(T1 a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator-(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator-(const vec<T1, N>& a, T2 b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator-(T1 a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator*(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator*(const vec<T1, N>& a, T2 b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator*(T1 a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator/(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator/(const vec<T1, N>& a, T2 b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator/(T1 a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator%(const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x % b.x};
  } else if constexpr (N == 2) {
    return {a.x % b.x, a.y % b.y};
  } else if constexpr (N == 3) {
    return {a.x % b.x, a.y % b.y, a.z % b.z};
  } else if constexpr (N == 4) {
    return {a.x % b.x, a.y % b.y, a.z % b.z, a.w % b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator%(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return {a.x % b};
  } else if constexpr (N == 2) {
    return {a.x % b, a.y % b};
  } else if constexpr (N == 3) {
    return {a.x % b, a.y % b, a.z % b};
  } else if constexpr (N == 4) {
    return {a.x % b, a.y % b, a.z % b, a.w % b};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator%(T1 a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x % b.x};
  } else if constexpr (N == 2) {
    return {a.x % b.x, a.y % b.y};
  } else if constexpr (N == 3) {
    return {a.x % b.x, a.y % b.y, a.z % b.z};
  } else if constexpr (N == 4) {
    return {a.x % b.x, a.y % b.y, a.z % b.z, a.w % b.w};
  }
}
template <typename T, size_t N>
constexpr kernel vec<T, N> operator~(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return {~a.x};
  } else if constexpr (N == 2) {
    return {~a.x, ~a.y};
  } else if constexpr (N == 3) {
    return {~a.x, ~a.y, ~a.z};
  } else if constexpr (N == 4) {
    return {~a.x, ~a.y, ~a.z, ~a.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator&(const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x & b.x};
  } else if constexpr (N == 2) {
    return {a.x & b.x, a.y & b.y};
  } else if constexpr (N == 3) {
    return {a.x & b.x, a.y & b.y, a.z & b.z};
  } else if constexpr (N == 4) {
    return {a.x & b.x, a.y & b.y, a.z & b.z, a.w & b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator&(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return {a.x & b};
  } else if constexpr (N == 2) {
    return {a.x & b, a.y & b};
  } else if constexpr (N == 3) {
    return {a.x & b, a.y & b, a.z & b};
  } else if constexpr (N == 4) {
    return {a.x & b, a.y & b, a.z & b, a.w & b};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator&(T1 a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x & b.x};
  } else if constexpr (N == 2) {
    return {a.x & b.x, a.y & b.y};
  } else if constexpr (N == 3) {
    return {a.x & b.x, a.y & b.y, a.z & b.z};
  } else if constexpr (N == 4) {
    return {a.x & b.x, a.y & b.y, a.z & b.z, a.w & b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator|(const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x | b.x};
  } else if constexpr (N == 2) {
    return {a.x | b.x, a.y | b.y};
  } else if constexpr (N == 3) {
    return {a.x | b.x, a.y | b.y, a.z | b.z};
  } else if constexpr (N == 4) {
    return {a.x | b.x, a.y | b.y, a.z | b.z, a.w | b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator|(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return {a.x | b};
  } else if constexpr (N == 2) {
    return {a.x | b, a.y | b};
  } else if constexpr (N == 3) {
    return {a.x | b, a.y | b, a.z | b};
  } else if constexpr (N == 4) {
    return {a.x | b, a.y | b, a.z | b, a.w | b};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator|(T1 a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x | b.x};
  } else if constexpr (N == 2) {
    return {a.x | b.x, a.y | b.y};
  } else if constexpr (N == 3) {
    return {a.x | b.x, a.y | b.y, a.z | b.z};
  } else if constexpr (N == 4) {
    return {a.x | b.x, a.y | b.y, a.z | b.z, a.w | b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator^(const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x ^ b.x};
  } else if constexpr (N == 2) {
    return {a.x ^ b.x, a.y ^ b.y};
  } else if constexpr (N == 3) {
    return {a.x ^ b.x, a.y ^ b.y, a.z ^ b.z};
  } else if constexpr (N == 4) {
    return {a.x ^ b.x, a.y ^ b.y, a.z ^ b.z, a.w ^ b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator^(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return {a.x ^ b};
  } else if constexpr (N == 2) {
    return {a.x ^ b, a.y ^ b};
  } else if constexpr (N == 3) {
    return {a.x ^ b, a.y ^ b, a.z ^ b};
  } else if constexpr (N == 4) {
    return {a.x ^ b, a.y ^ b, a.z ^ b, a.w ^ b};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator^(T1 a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x ^ b.x};
  } else if constexpr (N == 2) {
    return {a.x ^ b.x, a.y ^ b.y};
  } else if constexpr (N == 3) {
    return {a.x ^ b.x, a.y ^ b.y, a.z ^ b.z};
  } else if constexpr (N == 4) {
    return {a.x ^ b.x, a.y ^ b.y, a.z ^ b.z, a.w ^ b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator>>(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x >> b.x};
  } else if constexpr (N == 2) {
    return {a.x >> b.x, a.y >> b.y};
  } else if constexpr (N == 3) {
    return {a.x >> b.x, a.y >> b.y, a.z >> b.z};
  } else if constexpr (N == 4) {
    return {a.x >> b.x, a.y >> b.y, a.z >> b.z, a.w >> b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator>>(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return {a.x >> b};
  } else if constexpr (N == 2) {
    return {a.x >> b, a.y >> b};
  } else if constexpr (N == 3) {
    return {a.x >> b, a.y >> b, a.z >> b};
  } else if constexpr (N == 4) {
    return {a.x >> b, a.y >> b, a.z >> b, a.w >> b};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator>>(T1 a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x >> b.x};
  } else if constexpr (N == 2) {
    return {a.x >> b.x, a.y >> b.y};
  } else if constexpr (N == 3) {
    return {a.x >> b.x, a.y >> b.y, a.z >> b.z};
  } else if constexpr (N == 4) {
    return {a.x >> b.x, a.y >> b.y, a.z >> b.z, a.w >> b.w};
  }
}
template <typename T, size_t N>
constexpr kernel vec<T, N> operator<<(const vec<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return {a.x << b.x};
  } else if constexpr (N == 2) {
    return {a.x << b.x, a.y << b.y};
  } else if constexpr (N == 3) {
    return {a.x << b.x, a.y << b.y, a.z << b.z};
  } else if constexpr (N == 4) {
    return {a.x << b.x, a.y << b.y, a.z << b.z, a.w << b.w};
  }
}
template <typename T, size_t N>
constexpr kernel vec<T, N> operator<<(const vec<T, N>& a, T b) {
  if constexpr (N == 1) {
    return {a.x << b};
  } else if constexpr (N == 2) {
    return {a.x << b, a.y << b};
  } else if constexpr (N == 3) {
    return {a.x << b, a.y << b, a.z << b};
  } else if constexpr (N == 4) {
    return {a.x << b, a.y << b, a.z << b, a.w << b};
  }
}
template <typename T, size_t N>
constexpr kernel vec<T, N> operator<<(T a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return {a.x << b.x};
  } else if constexpr (N == 2) {
    return {a.x << b.x, a.y << b.y};
  } else if constexpr (N == 3) {
    return {a.x << b.x, a.y << b.y, a.z << b.z};
  } else if constexpr (N == 4) {
    return {a.x << b.x, a.y << b.y, a.z << b.z, a.w << b.w};
  }
}
template <size_t N>
constexpr kernel vec<bool, N> operator!(const vec<bool, N>& a) {
  if constexpr (N == 1) {
    return {!a.x};
  } else if constexpr (N == 2) {
    return {!a.x, !a.y};
  } else if constexpr (N == 3) {
    return {!a.x, !a.y, !a.z};
  } else if constexpr (N == 4) {
    return {!a.x, !a.y, !a.z, !a.w};
  }
}
template <size_t N>
constexpr kernel vec<bool, N> operator&&(
    const vec<bool, N>& a, const vec<bool, N>& b) {
  if constexpr (N == 1) {
    return {a.x && b.x};
  } else if constexpr (N == 2) {
    return {a.x && b.x, a.y && b.y};
  } else if constexpr (N == 3) {
    return {a.x && b.x, a.y && b.y, a.z && b.z};
  } else if constexpr (N == 4) {
    return {a.x && b.x, a.y && b.y, a.z && b.z, a.w && b.w};
  }
}
template <size_t N>
constexpr kernel vec<bool, N> operator&&(const vec<bool, N>& a, bool b) {
  if constexpr (N == 1) {
    return {a.x && b};
  } else if constexpr (N == 2) {
    return {a.x && b, a.y && b};
  } else if constexpr (N == 3) {
    return {a.x && b, a.y && b, a.z && b};
  } else if constexpr (N == 4) {
    return {a.x && b, a.y && b, a.z && b, a.w && b};
  }
}
template <size_t N>
constexpr kernel vec<bool, N> operator&&(bool a, const vec<bool, N>& b) {
  if constexpr (N == 1) {
    return {a && b.x};
  } else if constexpr (N == 2) {
    return {a && b.x, a && b.y};
  } else if constexpr (N == 3) {
    return {a && b.x, a && b.y, a && b.z};
  } else if constexpr (N == 4) {
    return {a && b.x, a && b.y, a && b.z, a && b.w};
  }
}
template <size_t N>
constexpr kernel vec<bool, N> operator||(
    const vec<bool, N>& a, const vec<bool, N>& b) {
  if constexpr (N == 1) {
    return {a.x || b.x};
  } else if constexpr (N == 2) {
    return {a.x || b.x, a.y || b.y};
  } else if constexpr (N == 3) {
    return {a.x || b.x, a.y || b.y, a.z || b.z};
  } else if constexpr (N == 4) {
    return {a.x || b.x, a.y || b.y, a.z || b.z, a.w || b.w};
  }
}
template <size_t N>
constexpr kernel vec<bool, N> operator||(const vec<bool, N>& a, bool b) {
  if constexpr (N == 1) {
    return {a.x || b};
  } else if constexpr (N == 2) {
    return {a.x || b, a.y || b};
  } else if constexpr (N == 3) {
    return {a.x || b, a.y || b, a.z || b};
  } else if constexpr (N == 4) {
    return {a.x || b, a.y || b, a.z || b, a.w || b};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator||(T1 a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a || b.x};
  } else if constexpr (N == 2) {
    return {a || b.x, a || b.y};
  } else if constexpr (N == 3) {
    return {a || b.x, a || b.y, a || b.z};
  } else if constexpr (N == 4) {
    return {a || b.x, a || b.y, a || b.z, a || b.w};
  }
}

// Vector assignments
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator+=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a + b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator+=(vec<T, N>& a, T1 b) {
  return a = a + b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator-=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a - b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator-=(vec<T, N>& a, T1 b) {
  return a = a - b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator*=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a * b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator*=(vec<T, N>& a, T1 b) {
  return a = a * b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator/=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a / b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator/=(vec<T, N>& a, T1 b) {
  return a = a / b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator%=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a % b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator%=(vec<T, N>& a, T1 b) {
  return a = a % b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator&=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a & b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator&=(vec<T, N>& a, T1 b) {
  return a = a & b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator|=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a | b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator|=(vec<T, N>& a, T1 b) {
  return a = a | b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator^=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a ^ b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator^=(vec<T, N>& a, T1 b) {
  return a = a ^ b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator>>=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a >> b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator>>=(vec<T, N>& a, T1 b) {
  return a = a >> b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator<<=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a << b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator<<=(vec<T, N>& a, T1 b) {
  return a = a << b;
}

// Vector products and lengths.
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel T dot(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, typename T = common_t<T1, T2>>
constexpr kernel T cross(const vec<T1, 2>& a, const vec<T2, 2>& b) {
  return a.x * b.y - a.y * b.x;
}
template <typename T1, typename T2, typename T = common_t<T1, T2>>
constexpr kernel vec<T, 3> cross(const vec<T1, 3>& a, const vec<T2, 3>& b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
template <typename T1, typename T2, typename T = common_t<T1, T2>>
constexpr kernel T angle(const vec<T1, 3>& a, const vec<T2, 3>& b) {
  return acos(clamp(dot(normalize(a), normalize(b)), (T)-1, (T)1));
}

// Orthogonal vectors.
template <typename T>
constexpr kernel vec<T, 3> orthogonal(const vec<T, 3>& v) {
  // http://lolengine.net/blog/2013/09/21/picking-orthogonal-vector-combing-coconuts)
  return abs(v.x) > abs(v.z) ? vec<T, 3>{-v.y, v.x, 0}
                             : vec<T, 3>{0, -v.z, v.y};
}
template <typename T1, typename T2, typename T = common_t<T1, T2>>
constexpr kernel vec<T, 3> orthonormalize(
    const vec<T1, 3>& a, const vec<T2, 3>& b) {
  return normalize(a - b * dot(a, b));
}

// Reflected and refracted vector.
template <typename T1, typename T2, typename T = common_t<T1, T2>>
constexpr kernel vec<T, 3> reflect(const vec<T1, 3>& w, const vec<T2, 3>& n) {
  return -w + 2 * dot(n, w) * n;
}
template <typename T1, typename T2, typename T3,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, 3> refract(
    const vec<T1, 3>& w, const vec<T2, 3>& n, T3 inv_eta) {
  auto cosine = dot(n, w);
  auto k      = 1 + inv_eta * inv_eta * (cosine * cosine - 1);
  if (k < 0) return {0, 0, 0};  // tir
  return -w * inv_eta + (inv_eta * cosine - sqrt(k)) * n;
}

template <typename T, size_t N>
constexpr kernel T length(const vec<T, N>& a) {
  return sqrt(dot(a, a));
}
template <typename T, size_t N>
constexpr kernel T length2(const vec<T, N>& a) {
  return dot(a, a);
}
template <typename T, size_t N>
[[deprecated]] constexpr kernel T length_squared(const vec<T, N>& a) {
  return dot(a, a);
}
template <typename T, size_t N>
constexpr kernel vec<T, N> normalize(const vec<T, N>& a) {
  auto l = length(a);
  return (l != 0) ? a / l : a;
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel T distance(const vec<T1, N>& a, const vec<T2, N>& b) {
  return length(a - b);
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel T distance2(const vec<T1, N>& a, const vec<T2, N>& b) {
  return dot(a - b, a - b);
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
[[deprecated]] constexpr kernel T distance_squared(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  return dot(a - b, a - b);
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel T angle(const vec<T1, N>& a, const vec<T2, N>& b) {
  return acos(clamp(dot(normalize(a), normalize(b)), (T)-1, (T)1));
}
template <typename T1, typename T2, typename T3, size_t N,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, N> slerp(
    const vec<T1, N>& a, const vec<T2, N>& b, T3 u) {
  static_assert(N != 1);
  // https://en.wikipedia.org/wiki/Slerp
  auto an = normalize(a), bn = normalize(b);
  auto d = dot(an, bn);
  if (d < 0) {
    bn = -bn;
    d  = -d;
  }
  if (d > (T)0.9995) return normalize(an + u * (bn - an));
  auto theta = acos(clamp(d, (T)-1, (T)1));
  if (theta == 0) return an;
  return an * (sin(theta * (1 - u)) / sin(theta)) +
         bn * (sin(theta * u) / sin(theta));
}

// Max element and clamp.
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> max(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> max(const vec<T1, N>& a, T2 b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> min(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> min(const vec<T1, N>& a, T2 b) {
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
template <typename T1, typename T2, typename T3, size_t N,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, N> clamp(
    const vec<T1, N>& x, const vec<T2, N>& min, const vec<T3, N>& max) {
  if constexpr (N == 1) {
    return {clamp(x.x, min.x, max.x)};
  } else if constexpr (N == 2) {
    return {clamp(x.x, min.x, max.x), clamp(x.y, min.y, max.y)};
  } else if constexpr (N == 3) {
    return {clamp(x.x, min.x, max.x), clamp(x.y, min.y, max.y),
        clamp(x.z, min.z, max.z)};
  } else if constexpr (N == 4) {
    return {clamp(x.x, min.x, max.x), clamp(x.y, min.y, max.y),
        clamp(x.z, min.z, max.z), clamp(x.w, min.w, max.w)};
  }
}
template <typename T1, typename T2, typename T3, size_t N,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, N> clamp(const vec<T1, N>& x, T2 min, T3 max) {
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
template <typename T1, typename T2, typename T3, size_t N,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, N> clamp(
    const vec<T1, N>& x, const vec<T2, N>& min, T3 max) {
  if constexpr (N == 1) {
    return {clamp(x.x, min.x, max)};
  } else if constexpr (N == 2) {
    return {clamp(x.x, min.x, max), clamp(x.y, min.y, max)};
  } else if constexpr (N == 3) {
    return {
        clamp(x.x, min.x, max), clamp(x.y, min.y, max), clamp(x.z, min.z, max)};
  } else if constexpr (N == 4) {
    return {clamp(x.x, min.x, max), clamp(x.y, min.y, max),
        clamp(x.z, min.z, max), clamp(x.w, min.w, max)};
  }
}
template <typename T1, typename T2, typename T3, size_t N,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, N> clamp(
    const vec<T1, N>& x, T2 min, const vec<T3, N>& max) {
  if constexpr (N == 1) {
    return {clamp(x.x, min, max.x)};
  } else if constexpr (N == 2) {
    return {clamp(x.x, min, max.x), clamp(x.y, min, max.y)};
  } else if constexpr (N == 3) {
    return {clamp(x.x, min, max.x), clamp(x.y, min, max.y),
        clamp(x.z, min.z, max.z)};
  } else if constexpr (N == 4) {
    return {clamp(x.x, min, max.x), clamp(x.y, min, max.y),
        clamp(x.z, min, max.z), clamp(x.w, min, max.w)};
  }
}
template <typename T1, typename T2, typename T3, size_t N,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, N> lerp(
    const vec<T1, N>& a, const vec<T2, N>& b, T3 u) {
  return a * (1 - u) + b * u;
}
template <typename T1, typename T2, typename T3, size_t N,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, N> lerp(
    const vec<T1, N>& a, const vec<T2, N>& b, const vec<T3, N>& u) {
  return a * (1 - u) + b * u;
}

template <typename T, size_t N>
constexpr kernel T max(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel T min(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel size_t argmax(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return 0;
  } else if constexpr (N == 2) {
    return a.x >= a.y ? 0 : 1;
  } else if constexpr (N == 3) {
    return a.x >= a.y ? (a.x >= a.z ? 0 : 2) : (a.y >= a.z ? 1 : 2);
  } else if constexpr (N == 4) {
    if (a.w >= a.x && a.w >= a.y && a.w >= a.z)
      return 3;
    else
      return argmax(xyz(a));
  }
}
template <typename T, size_t N>
constexpr kernel size_t argmin(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return 0;
  } else if constexpr (N == 2) {
    return a.x <= a.y ? 0 : 1;
  } else if constexpr (N == 3) {
    return a.x <= a.y ? (a.x <= a.z ? 0 : 2) : (a.y <= a.z ? 1 : 2);
  } else if constexpr (N == 4) {
    if (a.w <= a.x && a.w <= a.y && a.w <= a.z)
      return 3;
    else
      return argmin(xyz(a));
  }
}
template <typename T, size_t N, typename I = int>
constexpr kernel I find_index(const vec<T, N>& a, T b) {
  if constexpr (N == 1) {
    return a.x == b ? 0 : -1;
  } else if constexpr (N == 2) {
    return a.x == b ? 0 : (a.y == b ? 1 : -1);
  } else if constexpr (N == 3) {
    return a.x == b ? 0 : (a.y == b ? 1 : (a.z == b ? 2 : -1));
  } else if constexpr (N == 4) {
    return a.x == b ? 0 : (a.y == b ? 1 : (a.z == b ? 2 : (a.w == b ? 3 : -1)));
  }
}

template <typename T, size_t N>
constexpr kernel T sum(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel T prod(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return a.x;
  } else if constexpr (N == 2) {
    return a.x * a.y;
  } else if constexpr (N == 3) {
    return a.x * a.y * a.z;
  } else if constexpr (N == 4) {
    return a.x * a.y * a.z * a.w;
  }
}
template <typename T, size_t N>
constexpr kernel T mean(const vec<T, N>& a) {
  return sum(a) / N;
}
template <size_t N>
constexpr kernel bool all(const vec<bool, N>& a) {
  if constexpr (N == 1) {
    return a.x;
  } else if constexpr (N == 2) {
    return a.x && a.y;
  } else if constexpr (N == 3) {
    return a.x && a.y && a.z;
  } else if constexpr (N == 4) {
    return a.x && a.y && a.z && a.w;
  }
}
template <typename T, size_t N>
constexpr kernel T any(const vec<bool, N>& a) {
  if constexpr (N == 1) {
    return a.x;
  } else if constexpr (N == 2) {
    return a.x || a.y;
  } else if constexpr (N == 3) {
    return a.x || a.y || a.z;
  } else if constexpr (N == 4) {
    return a.x || a.y || a.z || a.w;
  }
}

// Functions applied to vector elements
template <typename T, size_t N>
constexpr kernel vec<T, N> abs(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel vec<T, N> sqr(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel vec<T, N> sqrt(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel vec<T, N> exp(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel vec<T, N> log(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel vec<T, N> exp2(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel vec<T, N> log2(const vec<T, N>& a) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> pow(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> pow(const vec<T1, N>& a, T2 b) {
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
template <typename T, size_t N>
constexpr kernel vec<T, N> round(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return {round(a.x)};
  } else if constexpr (N == 2) {
    return {round(a.x), round(a.y)};
  } else if constexpr (N == 3) {
    return {round(a.x), round(a.y), round(a.z)};
  } else if constexpr (N == 4) {
    return {round(a.x), round(a.y), round(a.z), round(a.w)};
  }
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N> fmod(const vec<T, N>& a, const vec<T1, N>& b) {
  if constexpr (N == 1) {
    return {fmod(a.x, b.x)};
  } else if constexpr (N == 2) {
    return {fmod(a.x, b.x), fmod(a.y, b.y)};
  } else if constexpr (N == 3) {
    return {fmod(a.x, b.x), fmod(a.y, b.y), fmod(a.z, b.z)};
  } else if constexpr (N == 4) {
    return {fmod(a.x, b.x), fmod(a.y, b.y), fmod(a.z, b.z), fmod(a.w, b.w)};
  }
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N> fmod(const vec<T, N>& a, T1 b) {
  if constexpr (N == 1) {
    return {fmod(a.x, b)};
  } else if constexpr (N == 2) {
    return {fmod(a.x, b), fmod(a.y, b)};
  } else if constexpr (N == 3) {
    return {fmod(a.x, b), fmod(a.y, b), fmod(a.z, b)};
  } else if constexpr (N == 4) {
    return {fmod(a.x, b), fmod(a.y, b), fmod(a.z, b), fmod(a.w, b)};
  }
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N> mod(const vec<T, N>& a, const vec<T1, N>& b) {
  if constexpr (N == 1) {
    return {mod(a.x, b.x)};
  } else if constexpr (N == 2) {
    return {mod(a.x, b.x), mod(a.y, b.y)};
  } else if constexpr (N == 3) {
    return {mod(a.x, b.x), mod(a.y, b.y), mod(a.z, b.z)};
  } else if constexpr (N == 4) {
    return {mod(a.x, b.x), mod(a.y, b.y), mod(a.z, b.z), mod(a.w, b.w)};
  }
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N> mod(const vec<T, N>& a, T1 b) {
  if constexpr (N == 1) {
    return {mod(a.x, b)};
  } else if constexpr (N == 2) {
    return {mod(a.x, b), mod(a.y, b)};
  } else if constexpr (N == 3) {
    return {mod(a.x, b), mod(a.y, b), mod(a.z, b)};
  } else if constexpr (N == 4) {
    return {mod(a.x, b), mod(a.y, b), mod(a.z, b), mod(a.w, b)};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> bias(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return {bias(a.x, b)};
  } else if constexpr (N == 2) {
    return {bias(a.x, b), bias(a.y, b)};
  } else if constexpr (N == 3) {
    return {bias(a.x, b), bias(a.y, b), bias(a.z, b)};
  } else if constexpr (N == 4) {
    return {bias(a.x, b), bias(a.y, b), bias(a.z, b), bias(a.w, b)};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> gain(const vec<T1, N>& a, T2 b) {
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
// template <typename T, size_t N>
// inline void swap(vec<T, N>& a, vec<T, N>& b) {
//   std::swap(a, b);
// }
template <typename T, size_t N>
constexpr kernel bool isfinite(const vec<T, N>& a) {
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

// Conversion between ranges
template <typename T>
constexpr kernel T unit_to_uv(T a) {
  return (a + 1) / 2;
}
template <typename T, size_t N>
constexpr kernel vec<T, N> unit_to_uv(const vec<T, N>& a) {
  return (a + 1) / 2;
}
template <typename T>
constexpr kernel T uv_to_unit(T a) {
  return a * 2 - 1;
}
template <typename T, size_t N>
constexpr kernel vec<T, N> uv_to_unit(const vec<T, N>& a) {
  return a * 2 - 1;
}

// Conversion between coordinates
template <typename T>
constexpr kernel vec<T, 2> cartesian_to_sphericaluv(const vec<T, 3>& w) {
  auto phi = atan2(w.y, w.x), theta = acos(clamp(w.z, -1, 1));
  return {mod(phi / (2 * (T)pi), 1), theta / (T)pi};
}
template <typename T>
constexpr kernel vec<T, 2> cartesiany_to_sphericaluv(const vec<T, 3>& w) {
  auto phi = atan2(w.z, w.x), theta = acos(clamp(w.y, -1, 1));
  return {mod(phi / (2 * (T)pi), 1), theta / (T)pi};
}
template <typename T>
constexpr kernel vec<T, 3> sphericaluv_to_cartesian(const vec<T, 2>& uv) {
  auto phi = uv.x * 2 * (T)pi, theta = uv.y * (T)pi;
  return {cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)};
}
template <typename T>
constexpr kernel vec<T, 3> sphericaluv_to_cartesiany(const vec<T, 2>& uv) {
  auto phi = uv.x * 2 * (T)pi, theta = uv.y * (T)pi;
  return {cos(phi) * sin(theta), cos(theta), sin(phi) * sin(theta)};
}

// Quaternion operations represented as xi + yj + zk + w
// const auto identity_quat4f = vec4f{0, 0, 0, 1};
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, 4> quat_mul(const vec<T1, 4>& a, T2 b) {
  return {a.x * b, a.y * b, a.z * b, a.w * b};
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, 4> quat_mul(const vec<T1, 4>& a, const vec<T2, 4>& b) {
  return {a.x * b.w + a.w * b.x + a.y * b.w - a.z * b.y,
      a.y * b.w + a.w * b.y + a.z * b.x - a.x * b.z,
      a.z * b.w + a.w * b.z + a.x * b.y - a.y * b.x,
      a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z};
}
template <typename T>
constexpr kernel vec<T, 4> quat_conjugate(const vec<T, 4>& a) {
  return {-a.x, -a.y, -a.z, a.w};
}
template <typename T>
constexpr kernel vec<T, 4> quat_inverse(const vec<T, 4>& a) {
  return quat_conjugate(a) / dot(a, a);
}

// Component-wise comparison operations.
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_equal(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x == b.x};
  } else if constexpr (N == 2) {
    return {a.x == b.x, a.y == b.y};
  } else if constexpr (N == 3) {
    return {a.x == b.x, a.y == b.y, a.z == b.z};
  } else if constexpr (N == 4) {
    return {a.x == b.x, a.y == b.y, a.z == b.z, a.w == b.w};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_equal(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return {a.x == b};
  } else if constexpr (N == 2) {
    return {a.x == b, a.y == b};
  } else if constexpr (N == 3) {
    return {a.x == b, a.y == b, a.z == b};
  } else if constexpr (N == 4) {
    return {a.x == b, a.y == b, a.z == b, a.w == b};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_not_equal(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x != b.x};
  } else if constexpr (N == 2) {
    return {a.x != b.x, a.y != b.y};
  } else if constexpr (N == 3) {
    return {a.x != b.x, a.y != b.y, a.z != b.z};
  } else if constexpr (N == 4) {
    return {a.x != b.x, a.y != b.y, a.z != b.z, a.w != b.w};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_not_equal(
    const vec<T1, N>& a, const T2 b) {
  if constexpr (N == 1) {
    return {a.x != b};
  } else if constexpr (N == 2) {
    return {a.x != b, a.y != b};
  } else if constexpr (N == 3) {
    return {a.x != b, a.y != b, a.z != b};
  } else if constexpr (N == 4) {
    return {a.x != b, a.y != b, a.z != b, a.w != b};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_less(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x < b.x};
  } else if constexpr (N == 2) {
    return {a.x < b.x, a.y < b.y};
  } else if constexpr (N == 3) {
    return {a.x < b.x, a.y < b.y, a.z < b.z};
  } else if constexpr (N == 4) {
    return {a.x < b.x, a.y < b.y, a.z < b.z, a.w < b.w};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_less(const vec<T1, N>& a, const T2 b) {
  if constexpr (N == 1) {
    return {a.x < b};
  } else if constexpr (N == 2) {
    return {a.x < b, a.y < b};
  } else if constexpr (N == 3) {
    return {a.x < b, a.y < b, a.z < b};
  } else if constexpr (N == 4) {
    return {a.x < b, a.y < b, a.z < b, a.w < b};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_greater(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x > b.x};
  } else if constexpr (N == 2) {
    return {a.x > b.x, a.y > b.y};
  } else if constexpr (N == 3) {
    return {a.x > b.x, a.y > b.y, a.z > b.z};
  } else if constexpr (N == 4) {
    return {a.x > b.x, a.y > b.y, a.z > b.z, a.w > b.w};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_greater(
    const vec<T1, N>& a, const T2 b) {
  if constexpr (N == 1) {
    return {a.x > b};
  } else if constexpr (N == 2) {
    return {a.x > b, a.y > b};
  } else if constexpr (N == 3) {
    return {a.x > b, a.y > b, a.z > b};
  } else if constexpr (N == 4) {
    return {a.x > b, a.y > b, a.z > b, a.w > b};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_less_equal(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x <= b.x};
  } else if constexpr (N == 2) {
    return {a.x <= b.x, a.y <= b.y};
  } else if constexpr (N == 3) {
    return {a.x <= b.x, a.y <= b.y, a.z <= b.z};
  } else if constexpr (N == 4) {
    return {a.x <= b.x, a.y <= b.y, a.z <= b.z, a.w <= b.w};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_less_equal(
    const vec<T1, N>& a, const T2 b) {
  if constexpr (N == 1) {
    return {a.x <= b};
  } else if constexpr (N == 2) {
    return {a.x <= b, a.y <= b};
  } else if constexpr (N == 3) {
    return {a.x <= b, a.y <= b, a.z <= b};
  } else if constexpr (N == 4) {
    return {a.x <= b, a.y <= b, a.z <= b, a.w <= b};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_greater_equal(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x >= b.x};
  } else if constexpr (N == 2) {
    return {a.x >= b.x, a.y >= b.y};
  } else if constexpr (N == 3) {
    return {a.x >= b.x, a.y >= b.y, a.z >= b.z};
  } else if constexpr (N == 4) {
    return {a.x >= b.x, a.y >= b.y, a.z >= b.z, a.w >= b.w};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_greater_equal(
    const vec<T1, N>& a, const T2 b) {
  if constexpr (N == 1) {
    return {a.x >= b};
  } else if constexpr (N == 2) {
    return {a.x >= b, a.y >= b};
  } else if constexpr (N == 3) {
    return {a.x >= b, a.y >= b, a.z >= b};
  } else if constexpr (N == 4) {
    return {a.x >= b, a.y >= b, a.z >= b, a.w >= b};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> select(
    const vec<bool, N>& a, const vec<T1, N>& b, const vec<T2, N>& c) {
  if constexpr (N == 1) {
    return {a.x ? b.x : c.x};
  } else if constexpr (N == 2) {
    return {a.x ? b.x : c.x, a.y ? b.y : c.y};
  } else if constexpr (N == 3) {
    return {a.x ? b.x : c.x, a.y ? b.y : c.y, a.z ? b.z : c.z};
  } else if constexpr (N == 4) {
    return {a.x ? b.x : c.x, a.y ? b.y : c.y, a.z ? b.z : c.z, a.w ? b.w : c.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> select(
    const vec<bool, N>& a, T1 b, const vec<T2, N>& c) {
  if constexpr (N == 1) {
    return {a.x ? b : c.x};
  } else if constexpr (N == 2) {
    return {a.x ? b : c.x, a.y ? b : c.y};
  } else if constexpr (N == 3) {
    return {a.x ? b : c.x, a.y ? b : c.y, a.z ? b : c.z};
  } else if constexpr (N == 4) {
    return {a.x ? b : c.x, a.y ? b : c.y, a.z ? b : c.z, a.w ? b : c.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> select(
    const vec<bool, N>& a, const vec<T1, N>& b, T2 c) {
  if constexpr (N == 1) {
    return {a.x ? b.x : c};
  } else if constexpr (N == 2) {
    return {a.x ? b.x : c, a.y ? b.y : c};
  } else if constexpr (N == 3) {
    return {a.x ? b.x : c, a.y ? b.y : c, a.z ? b.z : c};
  } else if constexpr (N == 4) {
    return {a.x ? b.x : c, a.y ? b.y : c, a.z ? b.z : c, a.w ? b.w : c};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> select(const vec<bool, N>& a, T1 b, T2 c) {
  if constexpr (N == 1) {
    return {a.x ? b : c};
  } else if constexpr (N == 2) {
    return {a.x ? b : c, a.y ? b : c};
  } else if constexpr (N == 3) {
    return {a.x ? b : c, a.y ? b : c, a.z ? b : c};
  } else if constexpr (N == 4) {
    return {a.x ? b : c, a.y ? b : c, a.z ? b : c, a.w ? b : c};
  }
}

#if 0

// ------------------------------------------------------------------

// Vector comparison operations.
template <typename T1, typename T2, size_t N>
constexpr kernel bool operator==(const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return a.x == b.x;
  } else if constexpr (N == 2) {
    return a.x == b.x && a.y == b.y;
  } else if constexpr (N == 3) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
  } else if constexpr (N == 4) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
  } else {
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel bool operator==(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return a.x == b;
  } else if constexpr (N == 2) {
    return a.x == b && a.y == b;
  } else if constexpr (N == 3) {
    return a.x == b && a.y == b && a.z == b;
  } else if constexpr (N == 4) {
    return a.x == b && a.y == b && a.z == b && a.w == b;
  } else {
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel bool operator!=(const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N != 1) {
    return a.x != b.x;
  } else if constexpr (N != 2) {
    return a.x != b.x || a.y != b.y;
  } else if constexpr (N != 3) {
    return a.x != b.x || a.y != b.y || a.z != b.z;
  } else if constexpr (N != 4) {
    return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w;
  } else {
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel bool operator!=(const vec<T1, N>& a, T2 b) {
  if constexpr (N != 1) {
    return a.x != b;
  } else if constexpr (N != 2) {
    return a.x != b || a.y != b;
  } else if constexpr (N != 3) {
    return a.x != b || a.y != b || a.z != b;
  } else if constexpr (N != 4) {
    return a.x != b || a.y != b || a.z != b || a.w != b;
  } else {
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel bool operator<(const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return a.x < b.x;
  } else if constexpr (N == 2) {
    return a.x < b.x || (a.x == b.x && a.y < b.y);
  } else if constexpr (N == 3) {
    return a.x < b.x ||
           (a.x == b.x && (a.y < b.y || (a.y == b.y && a.z < b.z)));
  } else if constexpr (N == 4) {
    return a.x < b.x ||
           (a.x == b.x &&
               (a.y < b.y ||
                   (a.y == b.y && (a.z < b.z || (a.z == b.z && a.w < b.w)))));
  } else {
  }
}

// Vector operations.
template <typename T, size_t N>
constexpr kernel vec<T, N> operator+(const vec<T, N>& a) {
  return a;
}
template <typename T, size_t N>
constexpr kernel vec<T, N> operator-(const vec<T, N>& a) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator+(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator+(const vec<T1, N>& a, T2 b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator+(T1 a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator-(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator-(const vec<T1, N>& a, T2 b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator-(T1 a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator*(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator*(const vec<T1, N>& a, T2 b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator*(T1 a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator/(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator/(const vec<T1, N>& a, T2 b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator/(T1 a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator%(const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x % b.x};
  } else if constexpr (N == 2) {
    return {a.x % b.x, a.y % b.y};
  } else if constexpr (N == 3) {
    return {a.x % b.x, a.y % b.y, a.z % b.z};
  } else if constexpr (N == 4) {
    return {a.x % b.x, a.y % b.y, a.z % b.z, a.w % b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator%(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return {a.x % b};
  } else if constexpr (N == 2) {
    return {a.x % b, a.y % b};
  } else if constexpr (N == 3) {
    return {a.x % b, a.y % b, a.z % b};
  } else if constexpr (N == 4) {
    return {a.x % b, a.y % b, a.z % b, a.w % b};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator%(T1 a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x % b.x};
  } else if constexpr (N == 2) {
    return {a.x % b.x, a.y % b.y};
  } else if constexpr (N == 3) {
    return {a.x % b.x, a.y % b.y, a.z % b.z};
  } else if constexpr (N == 4) {
    return {a.x % b.x, a.y % b.y, a.z % b.z, a.w % b.w};
  }
}
template <typename T, size_t N>
constexpr kernel vec<T, N> operator~(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return {~a.x};
  } else if constexpr (N == 2) {
    return {~a.x, ~a.y};
  } else if constexpr (N == 3) {
    return {~a.x, ~a.y, ~a.z};
  } else if constexpr (N == 4) {
    return {~a.x, ~a.y, ~a.z, ~a.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator&(const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x & b.x};
  } else if constexpr (N == 2) {
    return {a.x & b.x, a.y & b.y};
  } else if constexpr (N == 3) {
    return {a.x & b.x, a.y & b.y, a.z & b.z};
  } else if constexpr (N == 4) {
    return {a.x & b.x, a.y & b.y, a.z & b.z, a.w & b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator&(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return {a.x & b};
  } else if constexpr (N == 2) {
    return {a.x & b, a.y & b};
  } else if constexpr (N == 3) {
    return {a.x & b, a.y & b, a.z & b};
  } else if constexpr (N == 4) {
    return {a.x & b, a.y & b, a.z & b, a.w & b};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator&(T1 a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x & b.x};
  } else if constexpr (N == 2) {
    return {a.x & b.x, a.y & b.y};
  } else if constexpr (N == 3) {
    return {a.x & b.x, a.y & b.y, a.z & b.z};
  } else if constexpr (N == 4) {
    return {a.x & b.x, a.y & b.y, a.z & b.z, a.w & b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator|(const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x | b.x};
  } else if constexpr (N == 2) {
    return {a.x | b.x, a.y | b.y};
  } else if constexpr (N == 3) {
    return {a.x | b.x, a.y | b.y, a.z | b.z};
  } else if constexpr (N == 4) {
    return {a.x | b.x, a.y | b.y, a.z | b.z, a.w | b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator|(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return {a.x | b};
  } else if constexpr (N == 2) {
    return {a.x | b, a.y | b};
  } else if constexpr (N == 3) {
    return {a.x | b, a.y | b, a.z | b};
  } else if constexpr (N == 4) {
    return {a.x | b, a.y | b, a.z | b, a.w | b};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator|(T1 a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x | b.x};
  } else if constexpr (N == 2) {
    return {a.x | b.x, a.y | b.y};
  } else if constexpr (N == 3) {
    return {a.x | b.x, a.y | b.y, a.z | b.z};
  } else if constexpr (N == 4) {
    return {a.x | b.x, a.y | b.y, a.z | b.z, a.w | b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator^(const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x ^ b.x};
  } else if constexpr (N == 2) {
    return {a.x ^ b.x, a.y ^ b.y};
  } else if constexpr (N == 3) {
    return {a.x ^ b.x, a.y ^ b.y, a.z ^ b.z};
  } else if constexpr (N == 4) {
    return {a.x ^ b.x, a.y ^ b.y, a.z ^ b.z, a.w ^ b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator^(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return {a.x ^ b};
  } else if constexpr (N == 2) {
    return {a.x ^ b, a.y ^ b};
  } else if constexpr (N == 3) {
    return {a.x ^ b, a.y ^ b, a.z ^ b};
  } else if constexpr (N == 4) {
    return {a.x ^ b, a.y ^ b, a.z ^ b, a.w ^ b};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator^(T1 a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x ^ b.x};
  } else if constexpr (N == 2) {
    return {a.x ^ b.x, a.y ^ b.y};
  } else if constexpr (N == 3) {
    return {a.x ^ b.x, a.y ^ b.y, a.z ^ b.z};
  } else if constexpr (N == 4) {
    return {a.x ^ b.x, a.y ^ b.y, a.z ^ b.z, a.w ^ b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator>>(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x >> b.x};
  } else if constexpr (N == 2) {
    return {a.x >> b.x, a.y >> b.y};
  } else if constexpr (N == 3) {
    return {a.x >> b.x, a.y >> b.y, a.z >> b.z};
  } else if constexpr (N == 4) {
    return {a.x >> b.x, a.y >> b.y, a.z >> b.z, a.w >> b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator>>(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return {a.x >> b};
  } else if constexpr (N == 2) {
    return {a.x >> b, a.y >> b};
  } else if constexpr (N == 3) {
    return {a.x >> b, a.y >> b, a.z >> b};
  } else if constexpr (N == 4) {
    return {a.x >> b, a.y >> b, a.z >> b, a.w >> b};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator>>(T1 a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x >> b.x};
  } else if constexpr (N == 2) {
    return {a.x >> b.x, a.y >> b.y};
  } else if constexpr (N == 3) {
    return {a.x >> b.x, a.y >> b.y, a.z >> b.z};
  } else if constexpr (N == 4) {
    return {a.x >> b.x, a.y >> b.y, a.z >> b.z, a.w >> b.w};
  }
}
template <typename T, size_t N>
constexpr kernel vec<T, N> operator<<(const vec<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return {a.x << b.x};
  } else if constexpr (N == 2) {
    return {a.x << b.x, a.y << b.y};
  } else if constexpr (N == 3) {
    return {a.x << b.x, a.y << b.y, a.z << b.z};
  } else if constexpr (N == 4) {
    return {a.x << b.x, a.y << b.y, a.z << b.z, a.w << b.w};
  }
}
template <typename T, size_t N>
constexpr kernel vec<T, N> operator<<(const vec<T, N>& a, T b) {
  if constexpr (N == 1) {
    return {a.x << b};
  } else if constexpr (N == 2) {
    return {a.x << b, a.y << b};
  } else if constexpr (N == 3) {
    return {a.x << b, a.y << b, a.z << b};
  } else if constexpr (N == 4) {
    return {a.x << b, a.y << b, a.z << b, a.w << b};
  }
}
template <typename T, size_t N>
constexpr kernel vec<T, N> operator<<(T a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return {a.x << b.x};
  } else if constexpr (N == 2) {
    return {a.x << b.x, a.y << b.y};
  } else if constexpr (N == 3) {
    return {a.x << b.x, a.y << b.y, a.z << b.z};
  } else if constexpr (N == 4) {
    return {a.x << b.x, a.y << b.y, a.z << b.z, a.w << b.w};
  }
}
template <size_t N>
constexpr kernel vec<bool, N> operator!(const vec<bool, N>& a) {
  if constexpr (N == 1) {
    return {!a.x};
  } else if constexpr (N == 2) {
    return {!a.x, !a.y};
  } else if constexpr (N == 3) {
    return {!a.x, !a.y, !a.z};
  } else if constexpr (N == 4) {
    return {!a.x, !a.y, !a.z, !a.w};
  }
}
template <size_t N>
constexpr kernel vec<bool, N> operator&&(
    const vec<bool, N>& a, const vec<bool, N>& b) {
  if constexpr (N == 1) {
    return {a.x && b.x};
  } else if constexpr (N == 2) {
    return {a.x && b.x, a.y && b.y};
  } else if constexpr (N == 3) {
    return {a.x && b.x, a.y && b.y, a.z && b.z};
  } else if constexpr (N == 4) {
    return {a.x && b.x, a.y && b.y, a.z && b.z, a.w && b.w};
  }
}
template <size_t N>
constexpr kernel vec<bool, N> operator&&(const vec<bool, N>& a, bool b) {
  if constexpr (N == 1) {
    return {a.x && b};
  } else if constexpr (N == 2) {
    return {a.x && b, a.y && b};
  } else if constexpr (N == 3) {
    return {a.x && b, a.y && b, a.z && b};
  } else if constexpr (N == 4) {
    return {a.x && b, a.y && b, a.z && b, a.w && b};
  }
}
template <size_t N>
constexpr kernel vec<bool, N> operator&&(bool a, const vec<bool, N>& b) {
  if constexpr (N == 1) {
    return {a && b.x};
  } else if constexpr (N == 2) {
    return {a && b.x, a && b.y};
  } else if constexpr (N == 3) {
    return {a && b.x, a && b.y, a && b.z};
  } else if constexpr (N == 4) {
    return {a && b.x, a && b.y, a && b.z, a && b.w};
  }
}
template <size_t N>
constexpr kernel vec<bool, N> operator||(
    const vec<bool, N>& a, const vec<bool, N>& b) {
  if constexpr (N == 1) {
    return {a.x || b.x};
  } else if constexpr (N == 2) {
    return {a.x || b.x, a.y || b.y};
  } else if constexpr (N == 3) {
    return {a.x || b.x, a.y || b.y, a.z || b.z};
  } else if constexpr (N == 4) {
    return {a.x || b.x, a.y || b.y, a.z || b.z, a.w || b.w};
  }
}
template <size_t N>
constexpr kernel vec<bool, N> operator||(const vec<bool, N>& a, bool b) {
  if constexpr (N == 1) {
    return {a.x || b};
  } else if constexpr (N == 2) {
    return {a.x || b, a.y || b};
  } else if constexpr (N == 3) {
    return {a.x || b, a.y || b, a.z || b};
  } else if constexpr (N == 4) {
    return {a.x || b, a.y || b, a.z || b, a.w || b};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator||(T1 a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a || b.x};
  } else if constexpr (N == 2) {
    return {a || b.x, a || b.y};
  } else if constexpr (N == 3) {
    return {a || b.x, a || b.y, a || b.z};
  } else if constexpr (N == 4) {
    return {a || b.x, a || b.y, a || b.z, a || b.w};
  }
}

// Vector assignments
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator+=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a + b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator+=(vec<T, N>& a, T1 b) {
  return a = a + b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator-=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a - b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator-=(vec<T, N>& a, T1 b) {
  return a = a - b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator*=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a * b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator*=(vec<T, N>& a, T1 b) {
  return a = a * b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator/=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a / b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator/=(vec<T, N>& a, T1 b) {
  return a = a / b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator%=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a % b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator%=(vec<T, N>& a, T1 b) {
  return a = a % b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator&=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a & b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator&=(vec<T, N>& a, T1 b) {
  return a = a & b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator|=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a | b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator|=(vec<T, N>& a, T1 b) {
  return a = a | b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator^=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a ^ b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator^=(vec<T, N>& a, T1 b) {
  return a = a ^ b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator>>=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a >> b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator>>=(vec<T, N>& a, T1 b) {
  return a = a >> b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator<<=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a << b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator<<=(vec<T, N>& a, T1 b) {
  return a = a << b;
}

// Vector products and lengths.
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel T dot(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, typename T = common_t<T1, T2>>
constexpr kernel T cross(const vec<T1, 2>& a, const vec<T2, 2>& b) {
  return a.x * b.y - a.y * b.x;
}
template <typename T1, typename T2, typename T = common_t<T1, T2>>
constexpr kernel vec<T, 3> cross(const vec<T1, 3>& a, const vec<T2, 3>& b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
template <typename T1, typename T2, typename T = common_t<T1, T2>>
constexpr kernel T angle(const vec<T1, 3>& a, const vec<T2, 3>& b) {
  return acos(clamp(dot(normalize(a), normalize(b)), (T)-1, (T)1));
}

// Orthogonal vectors.
template <typename T>
constexpr kernel vec<T, 3> orthogonal(const vec<T, 3>& v) {
  // http://lolengine.net/blog/2013/09/21/picking-orthogonal-vector-combing-coconuts)
  return abs(v.x) > abs(v.z) ? vec<T, 3>{-v.y, v.x, 0}
                             : vec<T, 3>{0, -v.z, v.y};
}
template <typename T1, typename T2, typename T = common_t<T1, T2>>
constexpr kernel vec<T, 3> orthonormalize(
    const vec<T1, 3>& a, const vec<T2, 3>& b) {
  return normalize(a - b * dot(a, b));
}

// Reflected and refracted vector.
template <typename T1, typename T2, typename T = common_t<T1, T2>>
constexpr kernel vec<T, 3> reflect(const vec<T1, 3>& w, const vec<T2, 3>& n) {
  return -w + 2 * dot(n, w) * n;
}
template <typename T1, typename T2, typename T3,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, 3> refract(
    const vec<T1, 3>& w, const vec<T2, 3>& n, T3 inv_eta) {
  auto cosine = dot(n, w);
  auto k      = 1 + inv_eta * inv_eta * (cosine * cosine - 1);
  if (k < 0) return {0, 0, 0};  // tir
  return -w * inv_eta + (inv_eta * cosine - sqrt(k)) * n;
}

template <typename T, size_t N>
constexpr kernel T length(const vec<T, N>& a) {
  return sqrt(dot(a, a));
}
template <typename T, size_t N>
constexpr kernel T length2(const vec<T, N>& a) {
  return dot(a, a);
}
template <typename T, size_t N>
[[deprecated]] constexpr kernel T length_squared(const vec<T, N>& a) {
  return dot(a, a);
}
template <typename T, size_t N>
constexpr kernel vec<T, N> normalize(const vec<T, N>& a) {
  auto l = length(a);
  return (l != 0) ? a / l : a;
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel T distance(const vec<T1, N>& a, const vec<T2, N>& b) {
  return length(a - b);
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel T distance2(const vec<T1, N>& a, const vec<T2, N>& b) {
  return dot(a - b, a - b);
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
[[deprecated]] constexpr kernel T distance_squared(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  return dot(a - b, a - b);
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel T angle(const vec<T1, N>& a, const vec<T2, N>& b) {
  return acos(clamp(dot(normalize(a), normalize(b)), (T)-1, (T)1));
}
template <typename T1, typename T2, typename T3, size_t N,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, N> slerp(
    const vec<T1, N>& a, const vec<T2, N>& b, T3 u) {
  static_assert(N != 1);
  // https://en.wikipedia.org/wiki/Slerp
  auto an = normalize(a), bn = normalize(b);
  auto d = dot(an, bn);
  if (d < 0) {
    bn = -bn;
    d  = -d;
  }
  if (d > (T)0.9995) return normalize(an + u * (bn - an));
  auto theta = acos(clamp(d, (T)-1, (T)1));
  if (theta == 0) return an;
  return an * (sin(theta * (1 - u)) / sin(theta)) +
         bn * (sin(theta * u) / sin(theta));
}

// Max element and clamp.
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> max(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> max(const vec<T1, N>& a, T2 b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> min(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> min(const vec<T1, N>& a, T2 b) {
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
template <typename T1, typename T2, typename T3, size_t N,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, N> clamp(
    const vec<T1, N>& x, const vec<T2, N>& min, const vec<T3, N>& max) {
  if constexpr (N == 1) {
    return {clamp(x.x, min.x, max.x)};
  } else if constexpr (N == 2) {
    return {clamp(x.x, min.x, max.x), clamp(x.y, min.y, max.y)};
  } else if constexpr (N == 3) {
    return {clamp(x.x, min.x, max.x), clamp(x.y, min.y, max.y),
        clamp(x.z, min.z, max.z)};
  } else if constexpr (N == 4) {
    return {clamp(x.x, min.x, max.x), clamp(x.y, min.y, max.y),
        clamp(x.z, min.z, max.z), clamp(x.w, min.w, max.w)};
  }
}
template <typename T1, typename T2, typename T3, size_t N,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, N> clamp(const vec<T1, N>& x, T2 min, T3 max) {
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
template <typename T1, typename T2, typename T3, size_t N,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, N> clamp(
    const vec<T1, N>& x, const vec<T2, N>& min, T3 max) {
  if constexpr (N == 1) {
    return {clamp(x.x, min.x, max)};
  } else if constexpr (N == 2) {
    return {clamp(x.x, min.x, max), clamp(x.y, min.y, max)};
  } else if constexpr (N == 3) {
    return {
        clamp(x.x, min.x, max), clamp(x.y, min.y, max), clamp(x.z, min.z, max)};
  } else if constexpr (N == 4) {
    return {clamp(x.x, min.x, max), clamp(x.y, min.y, max),
        clamp(x.z, min.z, max), clamp(x.w, min.w, max)};
  }
}
template <typename T1, typename T2, typename T3, size_t N,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, N> clamp(
    const vec<T1, N>& x, T2 min, const vec<T3, N>& max) {
  if constexpr (N == 1) {
    return {clamp(x.x, min, max.x)};
  } else if constexpr (N == 2) {
    return {clamp(x.x, min, max.x), clamp(x.y, min, max.y)};
  } else if constexpr (N == 3) {
    return {clamp(x.x, min, max.x), clamp(x.y, min, max.y),
        clamp(x.z, min.z, max.z)};
  } else if constexpr (N == 4) {
    return {clamp(x.x, min, max.x), clamp(x.y, min, max.y),
        clamp(x.z, min, max.z), clamp(x.w, min, max.w)};
  }
}
template <typename T1, typename T2, typename T3, size_t N,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, N> lerp(
    const vec<T1, N>& a, const vec<T2, N>& b, T3 u) {
  return a * (1 - u) + b * u;
}
template <typename T1, typename T2, typename T3, size_t N,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, N> lerp(
    const vec<T1, N>& a, const vec<T2, N>& b, const vec<T3, N>& u) {
  return a * (1 - u) + b * u;
}

template <typename T, size_t N>
constexpr kernel T max(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel T min(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel size_t argmax(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return 0;
  } else if constexpr (N == 2) {
    return a.x >= a.y ? 0 : 1;
  } else if constexpr (N == 3) {
    return a.x >= a.y ? (a.x >= a.z ? 0 : 2) : (a.y >= a.z ? 1 : 2);
  } else if constexpr (N == 4) {
    if (a.w >= a.x && a.w >= a.y && a.w >= a.z)
      return 3;
    else
      return argmax(xyz(a));
  }
}
template <typename T, size_t N>
constexpr kernel size_t argmin(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return 0;
  } else if constexpr (N == 2) {
    return a.x <= a.y ? 0 : 1;
  } else if constexpr (N == 3) {
    return a.x <= a.y ? (a.x <= a.z ? 0 : 2) : (a.y <= a.z ? 1 : 2);
  } else if constexpr (N == 4) {
    if (a.w <= a.x && a.w <= a.y && a.w <= a.z)
      return 3;
    else
      return argmin(xyz(a));
  }
}
template <typename T, size_t N, typename I = int>
constexpr kernel I find_index(const vec<T, N>& a, T b) {
  if constexpr (N == 1) {
    return a.x == b ? 0 : -1;
  } else if constexpr (N == 2) {
    return a.x == b ? 0 : (a.y == b ? 1 : -1);
  } else if constexpr (N == 3) {
    return a.x == b ? 0 : (a.y == b ? 1 : (a.z == b ? 2 : -1));
  } else if constexpr (N == 4) {
    return a.x == b ? 0 : (a.y == b ? 1 : (a.z == b ? 2 : (a.w == b ? 3 : -1)));
  }
}

template <typename T, size_t N>
constexpr kernel T sum(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel T prod(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return a.x;
  } else if constexpr (N == 2) {
    return a.x * a.y;
  } else if constexpr (N == 3) {
    return a.x * a.y * a.z;
  } else if constexpr (N == 4) {
    return a.x * a.y * a.z * a.w;
  }
}
template <typename T, size_t N>
constexpr kernel T mean(const vec<T, N>& a) {
  return sum(a) / N;
}
template <size_t N>
constexpr kernel bool all(const vec<bool, N>& a) {
  if constexpr (N == 1) {
    return a.x;
  } else if constexpr (N == 2) {
    return a.x && a.y;
  } else if constexpr (N == 3) {
    return a.x && a.y && a.z;
  } else if constexpr (N == 4) {
    return a.x && a.y && a.z && a.w;
  }
}
template <typename T, size_t N>
constexpr kernel T any(const vec<bool, N>& a) {
  if constexpr (N == 1) {
    return a.x;
  } else if constexpr (N == 2) {
    return a.x || a.y;
  } else if constexpr (N == 3) {
    return a.x || a.y || a.z;
  } else if constexpr (N == 4) {
    return a.x || a.y || a.z || a.w;
  }
}

// Functions applied to vector elements
template <typename T, size_t N>
constexpr kernel vec<T, N> abs(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel vec<T, N> sqr(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel vec<T, N> sqrt(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel vec<T, N> exp(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel vec<T, N> log(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel vec<T, N> exp2(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel vec<T, N> log2(const vec<T, N>& a) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> pow(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> pow(const vec<T1, N>& a, T2 b) {
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
template <typename T, size_t N>
constexpr kernel vec<T, N> round(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return {round(a.x)};
  } else if constexpr (N == 2) {
    return {round(a.x), round(a.y)};
  } else if constexpr (N == 3) {
    return {round(a.x), round(a.y), round(a.z)};
  } else if constexpr (N == 4) {
    return {round(a.x), round(a.y), round(a.z), round(a.w)};
  }
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N> fmod(const vec<T, N>& a, const vec<T1, N>& b) {
  if constexpr (N == 1) {
    return {fmod(a.x, b.x)};
  } else if constexpr (N == 2) {
    return {fmod(a.x, b.x), fmod(a.y, b.y)};
  } else if constexpr (N == 3) {
    return {fmod(a.x, b.x), fmod(a.y, b.y), fmod(a.z, b.z)};
  } else if constexpr (N == 4) {
    return {fmod(a.x, b.x), fmod(a.y, b.y), fmod(a.z, b.z), fmod(a.w, b.w)};
  }
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N> fmod(const vec<T, N>& a, T1 b) {
  if constexpr (N == 1) {
    return {fmod(a.x, b)};
  } else if constexpr (N == 2) {
    return {fmod(a.x, b), fmod(a.y, b)};
  } else if constexpr (N == 3) {
    return {fmod(a.x, b), fmod(a.y, b), fmod(a.z, b)};
  } else if constexpr (N == 4) {
    return {fmod(a.x, b), fmod(a.y, b), fmod(a.z, b), fmod(a.w, b)};
  }
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N> mod(const vec<T, N>& a, const vec<T1, N>& b) {
  if constexpr (N == 1) {
    return {mod(a.x, b.x)};
  } else if constexpr (N == 2) {
    return {mod(a.x, b.x), mod(a.y, b.y)};
  } else if constexpr (N == 3) {
    return {mod(a.x, b.x), mod(a.y, b.y), mod(a.z, b.z)};
  } else if constexpr (N == 4) {
    return {mod(a.x, b.x), mod(a.y, b.y), mod(a.z, b.z), mod(a.w, b.w)};
  }
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N> mod(const vec<T, N>& a, T1 b) {
  if constexpr (N == 1) {
    return {mod(a.x, b)};
  } else if constexpr (N == 2) {
    return {mod(a.x, b), mod(a.y, b)};
  } else if constexpr (N == 3) {
    return {mod(a.x, b), mod(a.y, b), mod(a.z, b)};
  } else if constexpr (N == 4) {
    return {mod(a.x, b), mod(a.y, b), mod(a.z, b), mod(a.w, b)};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> bias(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return {bias(a.x, b)};
  } else if constexpr (N == 2) {
    return {bias(a.x, b), bias(a.y, b)};
  } else if constexpr (N == 3) {
    return {bias(a.x, b), bias(a.y, b), bias(a.z, b)};
  } else if constexpr (N == 4) {
    return {bias(a.x, b), bias(a.y, b), bias(a.z, b), bias(a.w, b)};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> gain(const vec<T1, N>& a, T2 b) {
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
// template <typename T, size_t N>
// inline void swap(vec<T, N>& a, vec<T, N>& b) {
//   std::swap(a, b);
// }
template <typename T, size_t N>
constexpr kernel bool isfinite(const vec<T, N>& a) {
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

// Conversion between ranges
template <typename T>
constexpr kernel T unit_to_uv(T a) {
  return (a + 1) / 2;
}
template <typename T, size_t N>
constexpr kernel vec<T, N> unit_to_uv(const vec<T, N>& a) {
  return (a + 1) / 2;
}
template <typename T>
constexpr kernel T uv_to_unit(T a) {
  return a * 2 - 1;
}
template <typename T, size_t N>
constexpr kernel vec<T, N> uv_to_unit(const vec<T, N>& a) {
  return a * 2 - 1;
}

// Conversion between coordinates
template <typename T>
constexpr kernel vec<T, 2> cartesian_to_sphericaluv(const vec<T, 3>& w) {
  auto phi = atan2(w.y, w.x), theta = acos(clamp(w.z, -1, 1));
  return {mod(phi / (2 * (T)pi), 1), theta / (T)pi};
}
template <typename T>
constexpr kernel vec<T, 2> cartesiany_to_sphericaluv(const vec<T, 3>& w) {
  auto phi = atan2(w.z, w.x), theta = acos(clamp(w.y, -1, 1));
  return {mod(phi / (2 * (T)pi), 1), theta / (T)pi};
}
template <typename T>
constexpr kernel vec<T, 3> sphericaluv_to_cartesian(const vec<T, 2>& uv) {
  auto phi = uv.x * 2 * (T)pi, theta = uv.y * (T)pi;
  return {cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)};
}
template <typename T>
constexpr kernel vec<T, 3> sphericaluv_to_cartesiany(const vec<T, 2>& uv) {
  auto phi = uv.x * 2 * (T)pi, theta = uv.y * (T)pi;
  return {cos(phi) * sin(theta), cos(theta), sin(phi) * sin(theta)};
}

// Quaternion operations represented as xi + yj + zk + w
// const auto identity_quat4f = vec4f{0, 0, 0, 1};
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, 4> quat_mul(const vec<T1, 4>& a, T2 b) {
  return {a.x * b, a.y * b, a.z * b, a.w * b};
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, 4> quat_mul(const vec<T1, 4>& a, const vec<T2, 4>& b) {
  return {a.x * b.w + a.w * b.x + a.y * b.w - a.z * b.y,
      a.y * b.w + a.w * b.y + a.z * b.x - a.x * b.z,
      a.z * b.w + a.w * b.z + a.x * b.y - a.y * b.x,
      a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z};
}
template <typename T>
constexpr kernel vec<T, 4> quat_conjugate(const vec<T, 4>& a) {
  return {-a.x, -a.y, -a.z, a.w};
}
template <typename T>
constexpr kernel vec<T, 4> quat_inverse(const vec<T, 4>& a) {
  return quat_conjugate(a) / dot(a, a);
}

// Component-wise comparison operations.
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_equal(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x == b.x};
  } else if constexpr (N == 2) {
    return {a.x == b.x, a.y == b.y};
  } else if constexpr (N == 3) {
    return {a.x == b.x, a.y == b.y, a.z == b.z};
  } else if constexpr (N == 4) {
    return {a.x == b.x, a.y == b.y, a.z == b.z, a.w == b.w};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_equal(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return {a.x == b};
  } else if constexpr (N == 2) {
    return {a.x == b, a.y == b};
  } else if constexpr (N == 3) {
    return {a.x == b, a.y == b, a.z == b};
  } else if constexpr (N == 4) {
    return {a.x == b, a.y == b, a.z == b, a.w == b};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_not_equal(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x != b.x};
  } else if constexpr (N == 2) {
    return {a.x != b.x, a.y != b.y};
  } else if constexpr (N == 3) {
    return {a.x != b.x, a.y != b.y, a.z != b.z};
  } else if constexpr (N == 4) {
    return {a.x != b.x, a.y != b.y, a.z != b.z, a.w != b.w};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_not_equal(
    const vec<T1, N>& a, const T2 b) {
  if constexpr (N == 1) {
    return {a.x != b};
  } else if constexpr (N == 2) {
    return {a.x != b, a.y != b};
  } else if constexpr (N == 3) {
    return {a.x != b, a.y != b, a.z != b};
  } else if constexpr (N == 4) {
    return {a.x != b, a.y != b, a.z != b, a.w != b};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_less(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x < b.x};
  } else if constexpr (N == 2) {
    return {a.x < b.x, a.y < b.y};
  } else if constexpr (N == 3) {
    return {a.x < b.x, a.y < b.y, a.z < b.z};
  } else if constexpr (N == 4) {
    return {a.x < b.x, a.y < b.y, a.z < b.z, a.w < b.w};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_less(const vec<T1, N>& a, const T2 b) {
  if constexpr (N == 1) {
    return {a.x < b};
  } else if constexpr (N == 2) {
    return {a.x < b, a.y < b};
  } else if constexpr (N == 3) {
    return {a.x < b, a.y < b, a.z < b};
  } else if constexpr (N == 4) {
    return {a.x < b, a.y < b, a.z < b, a.w < b};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_greater(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x > b.x};
  } else if constexpr (N == 2) {
    return {a.x > b.x, a.y > b.y};
  } else if constexpr (N == 3) {
    return {a.x > b.x, a.y > b.y, a.z > b.z};
  } else if constexpr (N == 4) {
    return {a.x > b.x, a.y > b.y, a.z > b.z, a.w > b.w};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_greater(
    const vec<T1, N>& a, const T2 b) {
  if constexpr (N == 1) {
    return {a.x > b};
  } else if constexpr (N == 2) {
    return {a.x > b, a.y > b};
  } else if constexpr (N == 3) {
    return {a.x > b, a.y > b, a.z > b};
  } else if constexpr (N == 4) {
    return {a.x > b, a.y > b, a.z > b, a.w > b};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_less_equal(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x <= b.x};
  } else if constexpr (N == 2) {
    return {a.x <= b.x, a.y <= b.y};
  } else if constexpr (N == 3) {
    return {a.x <= b.x, a.y <= b.y, a.z <= b.z};
  } else if constexpr (N == 4) {
    return {a.x <= b.x, a.y <= b.y, a.z <= b.z, a.w <= b.w};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_less_equal(
    const vec<T1, N>& a, const T2 b) {
  if constexpr (N == 1) {
    return {a.x <= b};
  } else if constexpr (N == 2) {
    return {a.x <= b, a.y <= b};
  } else if constexpr (N == 3) {
    return {a.x <= b, a.y <= b, a.z <= b};
  } else if constexpr (N == 4) {
    return {a.x <= b, a.y <= b, a.z <= b, a.w <= b};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_greater_equal(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x >= b.x};
  } else if constexpr (N == 2) {
    return {a.x >= b.x, a.y >= b.y};
  } else if constexpr (N == 3) {
    return {a.x >= b.x, a.y >= b.y, a.z >= b.z};
  } else if constexpr (N == 4) {
    return {a.x >= b.x, a.y >= b.y, a.z >= b.z, a.w >= b.w};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_greater_equal(
    const vec<T1, N>& a, const T2 b) {
  if constexpr (N == 1) {
    return {a.x >= b};
  } else if constexpr (N == 2) {
    return {a.x >= b, a.y >= b};
  } else if constexpr (N == 3) {
    return {a.x >= b, a.y >= b, a.z >= b};
  } else if constexpr (N == 4) {
    return {a.x >= b, a.y >= b, a.z >= b, a.w >= b};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> select(
    const vec<bool, N>& a, const vec<T1, N>& b, const vec<T2, N>& c) {
  if constexpr (N == 1) {
    return {a.x ? b.x : c.x};
  } else if constexpr (N == 2) {
    return {a.x ? b.x : c.x, a.y ? b.y : c.y};
  } else if constexpr (N == 3) {
    return {a.x ? b.x : c.x, a.y ? b.y : c.y, a.z ? b.z : c.z};
  } else if constexpr (N == 4) {
    return {a.x ? b.x : c.x, a.y ? b.y : c.y, a.z ? b.z : c.z, a.w ? b.w : c.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> select(
    const vec<bool, N>& a, T1 b, const vec<T2, N>& c) {
  if constexpr (N == 1) {
    return {a.x ? b : c.x};
  } else if constexpr (N == 2) {
    return {a.x ? b : c.x, a.y ? b : c.y};
  } else if constexpr (N == 3) {
    return {a.x ? b : c.x, a.y ? b : c.y, a.z ? b : c.z};
  } else if constexpr (N == 4) {
    return {a.x ? b : c.x, a.y ? b : c.y, a.z ? b : c.z, a.w ? b : c.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> select(
    const vec<bool, N>& a, const vec<T1, N>& b, T2 c) {
  if constexpr (N == 1) {
    return {a.x ? b.x : c};
  } else if constexpr (N == 2) {
    return {a.x ? b.x : c, a.y ? b.y : c};
  } else if constexpr (N == 3) {
    return {a.x ? b.x : c, a.y ? b.y : c, a.z ? b.z : c};
  } else if constexpr (N == 4) {
    return {a.x ? b.x : c, a.y ? b.y : c, a.z ? b.z : c, a.w ? b.w : c};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> select(const vec<bool, N>& a, T1 b, T2 c) {
  if constexpr (N == 1) {
    return {a.x ? b : c};
  } else if constexpr (N == 2) {
    return {a.x ? b : c, a.y ? b : c};
  } else if constexpr (N == 3) {
    return {a.x ? b : c, a.y ? b : c, a.z ? b : c};
  } else if constexpr (N == 4) {
    return {a.x ? b : c, a.y ? b : c, a.z ? b : c, a.w ? b : c};
  }
}




// ------------------------------------------------------------------

// Vector comparison operations.
template <typename T1, typename T2, size_t N>
constexpr kernel bool operator==(const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return a.x == b.x;
  } else if constexpr (N == 2) {
    return a.x == b.x && a.y == b.y;
  } else if constexpr (N == 3) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
  } else if constexpr (N == 4) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
  } else {
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel bool operator==(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return a.x == b;
  } else if constexpr (N == 2) {
    return a.x == b && a.y == b;
  } else if constexpr (N == 3) {
    return a.x == b && a.y == b && a.z == b;
  } else if constexpr (N == 4) {
    return a.x == b && a.y == b && a.z == b && a.w == b;
  } else {
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel bool operator!=(const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N != 1) {
    return a.x != b.x;
  } else if constexpr (N != 2) {
    return a.x != b.x || a.y != b.y;
  } else if constexpr (N != 3) {
    return a.x != b.x || a.y != b.y || a.z != b.z;
  } else if constexpr (N != 4) {
    return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w;
  } else {
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel bool operator!=(const vec<T1, N>& a, T2 b) {
  if constexpr (N != 1) {
    return a.x != b;
  } else if constexpr (N != 2) {
    return a.x != b || a.y != b;
  } else if constexpr (N != 3) {
    return a.x != b || a.y != b || a.z != b;
  } else if constexpr (N != 4) {
    return a.x != b || a.y != b || a.z != b || a.w != b;
  } else {
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel bool operator<(const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return a.x < b.x;
  } else if constexpr (N == 2) {
    return a.x < b.x || (a.x == b.x && a.y < b.y);
  } else if constexpr (N == 3) {
    return a.x < b.x ||
           (a.x == b.x && (a.y < b.y || (a.y == b.y && a.z < b.z)));
  } else if constexpr (N == 4) {
    return a.x < b.x ||
           (a.x == b.x &&
               (a.y < b.y ||
                   (a.y == b.y && (a.z < b.z || (a.z == b.z && a.w < b.w)))));
  } else {
  }
}

// Vector operations.
template <typename T, size_t N>
constexpr kernel vec<T, N> operator+(const vec<T, N>& a) {
  return a;
}
template <typename T, size_t N>
constexpr kernel vec<T, N> operator-(const vec<T, N>& a) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator+(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator+(const vec<T1, N>& a, T2 b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator+(T1 a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator-(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator-(const vec<T1, N>& a, T2 b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator-(T1 a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator*(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator*(const vec<T1, N>& a, T2 b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator*(T1 a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator/(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator/(const vec<T1, N>& a, T2 b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator/(T1 a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator%(const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x % b.x};
  } else if constexpr (N == 2) {
    return {a.x % b.x, a.y % b.y};
  } else if constexpr (N == 3) {
    return {a.x % b.x, a.y % b.y, a.z % b.z};
  } else if constexpr (N == 4) {
    return {a.x % b.x, a.y % b.y, a.z % b.z, a.w % b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator%(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return {a.x % b};
  } else if constexpr (N == 2) {
    return {a.x % b, a.y % b};
  } else if constexpr (N == 3) {
    return {a.x % b, a.y % b, a.z % b};
  } else if constexpr (N == 4) {
    return {a.x % b, a.y % b, a.z % b, a.w % b};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator%(T1 a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x % b.x};
  } else if constexpr (N == 2) {
    return {a.x % b.x, a.y % b.y};
  } else if constexpr (N == 3) {
    return {a.x % b.x, a.y % b.y, a.z % b.z};
  } else if constexpr (N == 4) {
    return {a.x % b.x, a.y % b.y, a.z % b.z, a.w % b.w};
  }
}
template <typename T, size_t N>
constexpr kernel vec<T, N> operator~(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return {~a.x};
  } else if constexpr (N == 2) {
    return {~a.x, ~a.y};
  } else if constexpr (N == 3) {
    return {~a.x, ~a.y, ~a.z};
  } else if constexpr (N == 4) {
    return {~a.x, ~a.y, ~a.z, ~a.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator&(const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x & b.x};
  } else if constexpr (N == 2) {
    return {a.x & b.x, a.y & b.y};
  } else if constexpr (N == 3) {
    return {a.x & b.x, a.y & b.y, a.z & b.z};
  } else if constexpr (N == 4) {
    return {a.x & b.x, a.y & b.y, a.z & b.z, a.w & b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator&(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return {a.x & b};
  } else if constexpr (N == 2) {
    return {a.x & b, a.y & b};
  } else if constexpr (N == 3) {
    return {a.x & b, a.y & b, a.z & b};
  } else if constexpr (N == 4) {
    return {a.x & b, a.y & b, a.z & b, a.w & b};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator&(T1 a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x & b.x};
  } else if constexpr (N == 2) {
    return {a.x & b.x, a.y & b.y};
  } else if constexpr (N == 3) {
    return {a.x & b.x, a.y & b.y, a.z & b.z};
  } else if constexpr (N == 4) {
    return {a.x & b.x, a.y & b.y, a.z & b.z, a.w & b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator|(const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x | b.x};
  } else if constexpr (N == 2) {
    return {a.x | b.x, a.y | b.y};
  } else if constexpr (N == 3) {
    return {a.x | b.x, a.y | b.y, a.z | b.z};
  } else if constexpr (N == 4) {
    return {a.x | b.x, a.y | b.y, a.z | b.z, a.w | b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator|(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return {a.x | b};
  } else if constexpr (N == 2) {
    return {a.x | b, a.y | b};
  } else if constexpr (N == 3) {
    return {a.x | b, a.y | b, a.z | b};
  } else if constexpr (N == 4) {
    return {a.x | b, a.y | b, a.z | b, a.w | b};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator|(T1 a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x | b.x};
  } else if constexpr (N == 2) {
    return {a.x | b.x, a.y | b.y};
  } else if constexpr (N == 3) {
    return {a.x | b.x, a.y | b.y, a.z | b.z};
  } else if constexpr (N == 4) {
    return {a.x | b.x, a.y | b.y, a.z | b.z, a.w | b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator^(const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x ^ b.x};
  } else if constexpr (N == 2) {
    return {a.x ^ b.x, a.y ^ b.y};
  } else if constexpr (N == 3) {
    return {a.x ^ b.x, a.y ^ b.y, a.z ^ b.z};
  } else if constexpr (N == 4) {
    return {a.x ^ b.x, a.y ^ b.y, a.z ^ b.z, a.w ^ b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator^(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return {a.x ^ b};
  } else if constexpr (N == 2) {
    return {a.x ^ b, a.y ^ b};
  } else if constexpr (N == 3) {
    return {a.x ^ b, a.y ^ b, a.z ^ b};
  } else if constexpr (N == 4) {
    return {a.x ^ b, a.y ^ b, a.z ^ b, a.w ^ b};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator^(T1 a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x ^ b.x};
  } else if constexpr (N == 2) {
    return {a.x ^ b.x, a.y ^ b.y};
  } else if constexpr (N == 3) {
    return {a.x ^ b.x, a.y ^ b.y, a.z ^ b.z};
  } else if constexpr (N == 4) {
    return {a.x ^ b.x, a.y ^ b.y, a.z ^ b.z, a.w ^ b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator>>(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x >> b.x};
  } else if constexpr (N == 2) {
    return {a.x >> b.x, a.y >> b.y};
  } else if constexpr (N == 3) {
    return {a.x >> b.x, a.y >> b.y, a.z >> b.z};
  } else if constexpr (N == 4) {
    return {a.x >> b.x, a.y >> b.y, a.z >> b.z, a.w >> b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator>>(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return {a.x >> b};
  } else if constexpr (N == 2) {
    return {a.x >> b, a.y >> b};
  } else if constexpr (N == 3) {
    return {a.x >> b, a.y >> b, a.z >> b};
  } else if constexpr (N == 4) {
    return {a.x >> b, a.y >> b, a.z >> b, a.w >> b};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator>>(T1 a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x >> b.x};
  } else if constexpr (N == 2) {
    return {a.x >> b.x, a.y >> b.y};
  } else if constexpr (N == 3) {
    return {a.x >> b.x, a.y >> b.y, a.z >> b.z};
  } else if constexpr (N == 4) {
    return {a.x >> b.x, a.y >> b.y, a.z >> b.z, a.w >> b.w};
  }
}
template <typename T, size_t N>
constexpr kernel vec<T, N> operator<<(const vec<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return {a.x << b.x};
  } else if constexpr (N == 2) {
    return {a.x << b.x, a.y << b.y};
  } else if constexpr (N == 3) {
    return {a.x << b.x, a.y << b.y, a.z << b.z};
  } else if constexpr (N == 4) {
    return {a.x << b.x, a.y << b.y, a.z << b.z, a.w << b.w};
  }
}
template <typename T, size_t N>
constexpr kernel vec<T, N> operator<<(const vec<T, N>& a, T b) {
  if constexpr (N == 1) {
    return {a.x << b};
  } else if constexpr (N == 2) {
    return {a.x << b, a.y << b};
  } else if constexpr (N == 3) {
    return {a.x << b, a.y << b, a.z << b};
  } else if constexpr (N == 4) {
    return {a.x << b, a.y << b, a.z << b, a.w << b};
  }
}
template <typename T, size_t N>
constexpr kernel vec<T, N> operator<<(T a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return {a.x << b.x};
  } else if constexpr (N == 2) {
    return {a.x << b.x, a.y << b.y};
  } else if constexpr (N == 3) {
    return {a.x << b.x, a.y << b.y, a.z << b.z};
  } else if constexpr (N == 4) {
    return {a.x << b.x, a.y << b.y, a.z << b.z, a.w << b.w};
  }
}
template <size_t N>
constexpr kernel vec<bool, N> operator!(const vec<bool, N>& a) {
  if constexpr (N == 1) {
    return {!a.x};
  } else if constexpr (N == 2) {
    return {!a.x, !a.y};
  } else if constexpr (N == 3) {
    return {!a.x, !a.y, !a.z};
  } else if constexpr (N == 4) {
    return {!a.x, !a.y, !a.z, !a.w};
  }
}
template <size_t N>
constexpr kernel vec<bool, N> operator&&(
    const vec<bool, N>& a, const vec<bool, N>& b) {
  if constexpr (N == 1) {
    return {a.x && b.x};
  } else if constexpr (N == 2) {
    return {a.x && b.x, a.y && b.y};
  } else if constexpr (N == 3) {
    return {a.x && b.x, a.y && b.y, a.z && b.z};
  } else if constexpr (N == 4) {
    return {a.x && b.x, a.y && b.y, a.z && b.z, a.w && b.w};
  }
}
template <size_t N>
constexpr kernel vec<bool, N> operator&&(const vec<bool, N>& a, bool b) {
  if constexpr (N == 1) {
    return {a.x && b};
  } else if constexpr (N == 2) {
    return {a.x && b, a.y && b};
  } else if constexpr (N == 3) {
    return {a.x && b, a.y && b, a.z && b};
  } else if constexpr (N == 4) {
    return {a.x && b, a.y && b, a.z && b, a.w && b};
  }
}
template <size_t N>
constexpr kernel vec<bool, N> operator&&(bool a, const vec<bool, N>& b) {
  if constexpr (N == 1) {
    return {a && b.x};
  } else if constexpr (N == 2) {
    return {a && b.x, a && b.y};
  } else if constexpr (N == 3) {
    return {a && b.x, a && b.y, a && b.z};
  } else if constexpr (N == 4) {
    return {a && b.x, a && b.y, a && b.z, a && b.w};
  }
}
template <size_t N>
constexpr kernel vec<bool, N> operator||(
    const vec<bool, N>& a, const vec<bool, N>& b) {
  if constexpr (N == 1) {
    return {a.x || b.x};
  } else if constexpr (N == 2) {
    return {a.x || b.x, a.y || b.y};
  } else if constexpr (N == 3) {
    return {a.x || b.x, a.y || b.y, a.z || b.z};
  } else if constexpr (N == 4) {
    return {a.x || b.x, a.y || b.y, a.z || b.z, a.w || b.w};
  }
}
template <size_t N>
constexpr kernel vec<bool, N> operator||(const vec<bool, N>& a, bool b) {
  if constexpr (N == 1) {
    return {a.x || b};
  } else if constexpr (N == 2) {
    return {a.x || b, a.y || b};
  } else if constexpr (N == 3) {
    return {a.x || b, a.y || b, a.z || b};
  } else if constexpr (N == 4) {
    return {a.x || b, a.y || b, a.z || b, a.w || b};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator||(T1 a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a || b.x};
  } else if constexpr (N == 2) {
    return {a || b.x, a || b.y};
  } else if constexpr (N == 3) {
    return {a || b.x, a || b.y, a || b.z};
  } else if constexpr (N == 4) {
    return {a || b.x, a || b.y, a || b.z, a || b.w};
  }
}

// Vector assignments
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator+=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a + b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator+=(vec<T, N>& a, T1 b) {
  return a = a + b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator-=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a - b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator-=(vec<T, N>& a, T1 b) {
  return a = a - b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator*=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a * b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator*=(vec<T, N>& a, T1 b) {
  return a = a * b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator/=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a / b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator/=(vec<T, N>& a, T1 b) {
  return a = a / b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator%=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a % b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator%=(vec<T, N>& a, T1 b) {
  return a = a % b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator&=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a & b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator&=(vec<T, N>& a, T1 b) {
  return a = a & b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator|=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a | b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator|=(vec<T, N>& a, T1 b) {
  return a = a | b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator^=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a ^ b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator^=(vec<T, N>& a, T1 b) {
  return a = a ^ b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator>>=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a >> b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator>>=(vec<T, N>& a, T1 b) {
  return a = a >> b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator<<=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a << b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator<<=(vec<T, N>& a, T1 b) {
  return a = a << b;
}

// Vector products and lengths.
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel T dot(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, typename T = common_t<T1, T2>>
constexpr kernel T cross(const vec<T1, 2>& a, const vec<T2, 2>& b) {
  return a.x * b.y - a.y * b.x;
}
template <typename T1, typename T2, typename T = common_t<T1, T2>>
constexpr kernel vec<T, 3> cross(const vec<T1, 3>& a, const vec<T2, 3>& b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
template <typename T1, typename T2, typename T = common_t<T1, T2>>
constexpr kernel T angle(const vec<T1, 3>& a, const vec<T2, 3>& b) {
  return acos(clamp(dot(normalize(a), normalize(b)), (T)-1, (T)1));
}

// Orthogonal vectors.
template <typename T>
constexpr kernel vec<T, 3> orthogonal(const vec<T, 3>& v) {
  // http://lolengine.net/blog/2013/09/21/picking-orthogonal-vector-combing-coconuts)
  return abs(v.x) > abs(v.z) ? vec<T, 3>{-v.y, v.x, 0}
                             : vec<T, 3>{0, -v.z, v.y};
}
template <typename T1, typename T2, typename T = common_t<T1, T2>>
constexpr kernel vec<T, 3> orthonormalize(
    const vec<T1, 3>& a, const vec<T2, 3>& b) {
  return normalize(a - b * dot(a, b));
}

// Reflected and refracted vector.
template <typename T1, typename T2, typename T = common_t<T1, T2>>
constexpr kernel vec<T, 3> reflect(const vec<T1, 3>& w, const vec<T2, 3>& n) {
  return -w + 2 * dot(n, w) * n;
}
template <typename T1, typename T2, typename T3,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, 3> refract(
    const vec<T1, 3>& w, const vec<T2, 3>& n, T3 inv_eta) {
  auto cosine = dot(n, w);
  auto k      = 1 + inv_eta * inv_eta * (cosine * cosine - 1);
  if (k < 0) return {0, 0, 0};  // tir
  return -w * inv_eta + (inv_eta * cosine - sqrt(k)) * n;
}

template <typename T, size_t N>
constexpr kernel T length(const vec<T, N>& a) {
  return sqrt(dot(a, a));
}
template <typename T, size_t N>
constexpr kernel T length2(const vec<T, N>& a) {
  return dot(a, a);
}
template <typename T, size_t N>
[[deprecated]] constexpr kernel T length_squared(const vec<T, N>& a) {
  return dot(a, a);
}
template <typename T, size_t N>
constexpr kernel vec<T, N> normalize(const vec<T, N>& a) {
  auto l = length(a);
  return (l != 0) ? a / l : a;
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel T distance(const vec<T1, N>& a, const vec<T2, N>& b) {
  return length(a - b);
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel T distance2(const vec<T1, N>& a, const vec<T2, N>& b) {
  return dot(a - b, a - b);
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
[[deprecated]] constexpr kernel T distance_squared(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  return dot(a - b, a - b);
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel T angle(const vec<T1, N>& a, const vec<T2, N>& b) {
  return acos(clamp(dot(normalize(a), normalize(b)), (T)-1, (T)1));
}
template <typename T1, typename T2, typename T3, size_t N,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, N> slerp(
    const vec<T1, N>& a, const vec<T2, N>& b, T3 u) {
  static_assert(N != 1);
  // https://en.wikipedia.org/wiki/Slerp
  auto an = normalize(a), bn = normalize(b);
  auto d = dot(an, bn);
  if (d < 0) {
    bn = -bn;
    d  = -d;
  }
  if (d > (T)0.9995) return normalize(an + u * (bn - an));
  auto theta = acos(clamp(d, (T)-1, (T)1));
  if (theta == 0) return an;
  return an * (sin(theta * (1 - u)) / sin(theta)) +
         bn * (sin(theta * u) / sin(theta));
}

// Max element and clamp.
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> max(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> max(const vec<T1, N>& a, T2 b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> min(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> min(const vec<T1, N>& a, T2 b) {
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
template <typename T1, typename T2, typename T3, size_t N,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, N> clamp(
    const vec<T1, N>& x, const vec<T2, N>& min, const vec<T3, N>& max) {
  if constexpr (N == 1) {
    return {clamp(x.x, min.x, max.x)};
  } else if constexpr (N == 2) {
    return {clamp(x.x, min.x, max.x), clamp(x.y, min.y, max.y)};
  } else if constexpr (N == 3) {
    return {clamp(x.x, min.x, max.x), clamp(x.y, min.y, max.y),
        clamp(x.z, min.z, max.z)};
  } else if constexpr (N == 4) {
    return {clamp(x.x, min.x, max.x), clamp(x.y, min.y, max.y),
        clamp(x.z, min.z, max.z), clamp(x.w, min.w, max.w)};
  }
}
template <typename T1, typename T2, typename T3, size_t N,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, N> clamp(const vec<T1, N>& x, T2 min, T3 max) {
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
template <typename T1, typename T2, typename T3, size_t N,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, N> clamp(
    const vec<T1, N>& x, const vec<T2, N>& min, T3 max) {
  if constexpr (N == 1) {
    return {clamp(x.x, min.x, max)};
  } else if constexpr (N == 2) {
    return {clamp(x.x, min.x, max), clamp(x.y, min.y, max)};
  } else if constexpr (N == 3) {
    return {
        clamp(x.x, min.x, max), clamp(x.y, min.y, max), clamp(x.z, min.z, max)};
  } else if constexpr (N == 4) {
    return {clamp(x.x, min.x, max), clamp(x.y, min.y, max),
        clamp(x.z, min.z, max), clamp(x.w, min.w, max)};
  }
}
template <typename T1, typename T2, typename T3, size_t N,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, N> clamp(
    const vec<T1, N>& x, T2 min, const vec<T3, N>& max) {
  if constexpr (N == 1) {
    return {clamp(x.x, min, max.x)};
  } else if constexpr (N == 2) {
    return {clamp(x.x, min, max.x), clamp(x.y, min, max.y)};
  } else if constexpr (N == 3) {
    return {clamp(x.x, min, max.x), clamp(x.y, min, max.y),
        clamp(x.z, min.z, max.z)};
  } else if constexpr (N == 4) {
    return {clamp(x.x, min, max.x), clamp(x.y, min, max.y),
        clamp(x.z, min, max.z), clamp(x.w, min, max.w)};
  }
}
template <typename T1, typename T2, typename T3, size_t N,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, N> lerp(
    const vec<T1, N>& a, const vec<T2, N>& b, T3 u) {
  return a * (1 - u) + b * u;
}
template <typename T1, typename T2, typename T3, size_t N,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, N> lerp(
    const vec<T1, N>& a, const vec<T2, N>& b, const vec<T3, N>& u) {
  return a * (1 - u) + b * u;
}

template <typename T, size_t N>
constexpr kernel T max(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel T min(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel size_t argmax(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return 0;
  } else if constexpr (N == 2) {
    return a.x >= a.y ? 0 : 1;
  } else if constexpr (N == 3) {
    return a.x >= a.y ? (a.x >= a.z ? 0 : 2) : (a.y >= a.z ? 1 : 2);
  } else if constexpr (N == 4) {
    if (a.w >= a.x && a.w >= a.y && a.w >= a.z)
      return 3;
    else
      return argmax(xyz(a));
  }
}
template <typename T, size_t N>
constexpr kernel size_t argmin(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return 0;
  } else if constexpr (N == 2) {
    return a.x <= a.y ? 0 : 1;
  } else if constexpr (N == 3) {
    return a.x <= a.y ? (a.x <= a.z ? 0 : 2) : (a.y <= a.z ? 1 : 2);
  } else if constexpr (N == 4) {
    if (a.w <= a.x && a.w <= a.y && a.w <= a.z)
      return 3;
    else
      return argmin(xyz(a));
  }
}
template <typename T, size_t N, typename I = int>
constexpr kernel I find_index(const vec<T, N>& a, T b) {
  if constexpr (N == 1) {
    return a.x == b ? 0 : -1;
  } else if constexpr (N == 2) {
    return a.x == b ? 0 : (a.y == b ? 1 : -1);
  } else if constexpr (N == 3) {
    return a.x == b ? 0 : (a.y == b ? 1 : (a.z == b ? 2 : -1));
  } else if constexpr (N == 4) {
    return a.x == b ? 0 : (a.y == b ? 1 : (a.z == b ? 2 : (a.w == b ? 3 : -1)));
  }
}

template <typename T, size_t N>
constexpr kernel T sum(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel T prod(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return a.x;
  } else if constexpr (N == 2) {
    return a.x * a.y;
  } else if constexpr (N == 3) {
    return a.x * a.y * a.z;
  } else if constexpr (N == 4) {
    return a.x * a.y * a.z * a.w;
  }
}
template <typename T, size_t N>
constexpr kernel T mean(const vec<T, N>& a) {
  return sum(a) / N;
}
template <size_t N>
constexpr kernel bool all(const vec<bool, N>& a) {
  if constexpr (N == 1) {
    return a.x;
  } else if constexpr (N == 2) {
    return a.x && a.y;
  } else if constexpr (N == 3) {
    return a.x && a.y && a.z;
  } else if constexpr (N == 4) {
    return a.x && a.y && a.z && a.w;
  }
}
template <typename T, size_t N>
constexpr kernel T any(const vec<bool, N>& a) {
  if constexpr (N == 1) {
    return a.x;
  } else if constexpr (N == 2) {
    return a.x || a.y;
  } else if constexpr (N == 3) {
    return a.x || a.y || a.z;
  } else if constexpr (N == 4) {
    return a.x || a.y || a.z || a.w;
  }
}

// Functions applied to vector elements
template <typename T, size_t N>
constexpr kernel vec<T, N> abs(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel vec<T, N> sqr(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel vec<T, N> sqrt(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel vec<T, N> exp(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel vec<T, N> log(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel vec<T, N> exp2(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel vec<T, N> log2(const vec<T, N>& a) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> pow(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> pow(const vec<T1, N>& a, T2 b) {
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
template <typename T, size_t N>
constexpr kernel vec<T, N> round(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return {round(a.x)};
  } else if constexpr (N == 2) {
    return {round(a.x), round(a.y)};
  } else if constexpr (N == 3) {
    return {round(a.x), round(a.y), round(a.z)};
  } else if constexpr (N == 4) {
    return {round(a.x), round(a.y), round(a.z), round(a.w)};
  }
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N> fmod(const vec<T, N>& a, const vec<T1, N>& b) {
  if constexpr (N == 1) {
    return {fmod(a.x, b.x)};
  } else if constexpr (N == 2) {
    return {fmod(a.x, b.x), fmod(a.y, b.y)};
  } else if constexpr (N == 3) {
    return {fmod(a.x, b.x), fmod(a.y, b.y), fmod(a.z, b.z)};
  } else if constexpr (N == 4) {
    return {fmod(a.x, b.x), fmod(a.y, b.y), fmod(a.z, b.z), fmod(a.w, b.w)};
  }
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N> fmod(const vec<T, N>& a, T1 b) {
  if constexpr (N == 1) {
    return {fmod(a.x, b)};
  } else if constexpr (N == 2) {
    return {fmod(a.x, b), fmod(a.y, b)};
  } else if constexpr (N == 3) {
    return {fmod(a.x, b), fmod(a.y, b), fmod(a.z, b)};
  } else if constexpr (N == 4) {
    return {fmod(a.x, b), fmod(a.y, b), fmod(a.z, b), fmod(a.w, b)};
  }
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N> mod(const vec<T, N>& a, const vec<T1, N>& b) {
  if constexpr (N == 1) {
    return {mod(a.x, b.x)};
  } else if constexpr (N == 2) {
    return {mod(a.x, b.x), mod(a.y, b.y)};
  } else if constexpr (N == 3) {
    return {mod(a.x, b.x), mod(a.y, b.y), mod(a.z, b.z)};
  } else if constexpr (N == 4) {
    return {mod(a.x, b.x), mod(a.y, b.y), mod(a.z, b.z), mod(a.w, b.w)};
  }
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N> mod(const vec<T, N>& a, T1 b) {
  if constexpr (N == 1) {
    return {mod(a.x, b)};
  } else if constexpr (N == 2) {
    return {mod(a.x, b), mod(a.y, b)};
  } else if constexpr (N == 3) {
    return {mod(a.x, b), mod(a.y, b), mod(a.z, b)};
  } else if constexpr (N == 4) {
    return {mod(a.x, b), mod(a.y, b), mod(a.z, b), mod(a.w, b)};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> bias(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return {bias(a.x, b)};
  } else if constexpr (N == 2) {
    return {bias(a.x, b), bias(a.y, b)};
  } else if constexpr (N == 3) {
    return {bias(a.x, b), bias(a.y, b), bias(a.z, b)};
  } else if constexpr (N == 4) {
    return {bias(a.x, b), bias(a.y, b), bias(a.z, b), bias(a.w, b)};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> gain(const vec<T1, N>& a, T2 b) {
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
// template <typename T, size_t N>
// inline void swap(vec<T, N>& a, vec<T, N>& b) {
//   std::swap(a, b);
// }
template <typename T, size_t N>
constexpr kernel bool isfinite(const vec<T, N>& a) {
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

// Conversion between ranges
template <typename T>
constexpr kernel T unit_to_uv(T a) {
  return (a + 1) / 2;
}
template <typename T, size_t N>
constexpr kernel vec<T, N> unit_to_uv(const vec<T, N>& a) {
  return (a + 1) / 2;
}
template <typename T>
constexpr kernel T uv_to_unit(T a) {
  return a * 2 - 1;
}
template <typename T, size_t N>
constexpr kernel vec<T, N> uv_to_unit(const vec<T, N>& a) {
  return a * 2 - 1;
}

// Conversion between coordinates
template <typename T>
constexpr kernel vec<T, 2> cartesian_to_sphericaluv(const vec<T, 3>& w) {
  auto phi = atan2(w.y, w.x), theta = acos(clamp(w.z, -1, 1));
  return {mod(phi / (2 * (T)pi), 1), theta / (T)pi};
}
template <typename T>
constexpr kernel vec<T, 2> cartesiany_to_sphericaluv(const vec<T, 3>& w) {
  auto phi = atan2(w.z, w.x), theta = acos(clamp(w.y, -1, 1));
  return {mod(phi / (2 * (T)pi), 1), theta / (T)pi};
}
template <typename T>
constexpr kernel vec<T, 3> sphericaluv_to_cartesian(const vec<T, 2>& uv) {
  auto phi = uv.x * 2 * (T)pi, theta = uv.y * (T)pi;
  return {cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)};
}
template <typename T>
constexpr kernel vec<T, 3> sphericaluv_to_cartesiany(const vec<T, 2>& uv) {
  auto phi = uv.x * 2 * (T)pi, theta = uv.y * (T)pi;
  return {cos(phi) * sin(theta), cos(theta), sin(phi) * sin(theta)};
}

// Quaternion operations represented as xi + yj + zk + w
// const auto identity_quat4f = vec4f{0, 0, 0, 1};
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, 4> quat_mul(const vec<T1, 4>& a, T2 b) {
  return {a.x * b, a.y * b, a.z * b, a.w * b};
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, 4> quat_mul(const vec<T1, 4>& a, const vec<T2, 4>& b) {
  return {a.x * b.w + a.w * b.x + a.y * b.w - a.z * b.y,
      a.y * b.w + a.w * b.y + a.z * b.x - a.x * b.z,
      a.z * b.w + a.w * b.z + a.x * b.y - a.y * b.x,
      a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z};
}
template <typename T>
constexpr kernel vec<T, 4> quat_conjugate(const vec<T, 4>& a) {
  return {-a.x, -a.y, -a.z, a.w};
}
template <typename T>
constexpr kernel vec<T, 4> quat_inverse(const vec<T, 4>& a) {
  return quat_conjugate(a) / dot(a, a);
}

// Component-wise comparison operations.
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_equal(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x == b.x};
  } else if constexpr (N == 2) {
    return {a.x == b.x, a.y == b.y};
  } else if constexpr (N == 3) {
    return {a.x == b.x, a.y == b.y, a.z == b.z};
  } else if constexpr (N == 4) {
    return {a.x == b.x, a.y == b.y, a.z == b.z, a.w == b.w};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_equal(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return {a.x == b};
  } else if constexpr (N == 2) {
    return {a.x == b, a.y == b};
  } else if constexpr (N == 3) {
    return {a.x == b, a.y == b, a.z == b};
  } else if constexpr (N == 4) {
    return {a.x == b, a.y == b, a.z == b, a.w == b};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_not_equal(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x != b.x};
  } else if constexpr (N == 2) {
    return {a.x != b.x, a.y != b.y};
  } else if constexpr (N == 3) {
    return {a.x != b.x, a.y != b.y, a.z != b.z};
  } else if constexpr (N == 4) {
    return {a.x != b.x, a.y != b.y, a.z != b.z, a.w != b.w};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_not_equal(
    const vec<T1, N>& a, const T2 b) {
  if constexpr (N == 1) {
    return {a.x != b};
  } else if constexpr (N == 2) {
    return {a.x != b, a.y != b};
  } else if constexpr (N == 3) {
    return {a.x != b, a.y != b, a.z != b};
  } else if constexpr (N == 4) {
    return {a.x != b, a.y != b, a.z != b, a.w != b};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_less(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x < b.x};
  } else if constexpr (N == 2) {
    return {a.x < b.x, a.y < b.y};
  } else if constexpr (N == 3) {
    return {a.x < b.x, a.y < b.y, a.z < b.z};
  } else if constexpr (N == 4) {
    return {a.x < b.x, a.y < b.y, a.z < b.z, a.w < b.w};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_less(const vec<T1, N>& a, const T2 b) {
  if constexpr (N == 1) {
    return {a.x < b};
  } else if constexpr (N == 2) {
    return {a.x < b, a.y < b};
  } else if constexpr (N == 3) {
    return {a.x < b, a.y < b, a.z < b};
  } else if constexpr (N == 4) {
    return {a.x < b, a.y < b, a.z < b, a.w < b};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_greater(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x > b.x};
  } else if constexpr (N == 2) {
    return {a.x > b.x, a.y > b.y};
  } else if constexpr (N == 3) {
    return {a.x > b.x, a.y > b.y, a.z > b.z};
  } else if constexpr (N == 4) {
    return {a.x > b.x, a.y > b.y, a.z > b.z, a.w > b.w};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_greater(
    const vec<T1, N>& a, const T2 b) {
  if constexpr (N == 1) {
    return {a.x > b};
  } else if constexpr (N == 2) {
    return {a.x > b, a.y > b};
  } else if constexpr (N == 3) {
    return {a.x > b, a.y > b, a.z > b};
  } else if constexpr (N == 4) {
    return {a.x > b, a.y > b, a.z > b, a.w > b};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_less_equal(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x <= b.x};
  } else if constexpr (N == 2) {
    return {a.x <= b.x, a.y <= b.y};
  } else if constexpr (N == 3) {
    return {a.x <= b.x, a.y <= b.y, a.z <= b.z};
  } else if constexpr (N == 4) {
    return {a.x <= b.x, a.y <= b.y, a.z <= b.z, a.w <= b.w};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_less_equal(
    const vec<T1, N>& a, const T2 b) {
  if constexpr (N == 1) {
    return {a.x <= b};
  } else if constexpr (N == 2) {
    return {a.x <= b, a.y <= b};
  } else if constexpr (N == 3) {
    return {a.x <= b, a.y <= b, a.z <= b};
  } else if constexpr (N == 4) {
    return {a.x <= b, a.y <= b, a.z <= b, a.w <= b};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_greater_equal(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x >= b.x};
  } else if constexpr (N == 2) {
    return {a.x >= b.x, a.y >= b.y};
  } else if constexpr (N == 3) {
    return {a.x >= b.x, a.y >= b.y, a.z >= b.z};
  } else if constexpr (N == 4) {
    return {a.x >= b.x, a.y >= b.y, a.z >= b.z, a.w >= b.w};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_greater_equal(
    const vec<T1, N>& a, const T2 b) {
  if constexpr (N == 1) {
    return {a.x >= b};
  } else if constexpr (N == 2) {
    return {a.x >= b, a.y >= b};
  } else if constexpr (N == 3) {
    return {a.x >= b, a.y >= b, a.z >= b};
  } else if constexpr (N == 4) {
    return {a.x >= b, a.y >= b, a.z >= b, a.w >= b};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> select(
    const vec<bool, N>& a, const vec<T1, N>& b, const vec<T2, N>& c) {
  if constexpr (N == 1) {
    return {a.x ? b.x : c.x};
  } else if constexpr (N == 2) {
    return {a.x ? b.x : c.x, a.y ? b.y : c.y};
  } else if constexpr (N == 3) {
    return {a.x ? b.x : c.x, a.y ? b.y : c.y, a.z ? b.z : c.z};
  } else if constexpr (N == 4) {
    return {a.x ? b.x : c.x, a.y ? b.y : c.y, a.z ? b.z : c.z, a.w ? b.w : c.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> select(
    const vec<bool, N>& a, T1 b, const vec<T2, N>& c) {
  if constexpr (N == 1) {
    return {a.x ? b : c.x};
  } else if constexpr (N == 2) {
    return {a.x ? b : c.x, a.y ? b : c.y};
  } else if constexpr (N == 3) {
    return {a.x ? b : c.x, a.y ? b : c.y, a.z ? b : c.z};
  } else if constexpr (N == 4) {
    return {a.x ? b : c.x, a.y ? b : c.y, a.z ? b : c.z, a.w ? b : c.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> select(
    const vec<bool, N>& a, const vec<T1, N>& b, T2 c) {
  if constexpr (N == 1) {
    return {a.x ? b.x : c};
  } else if constexpr (N == 2) {
    return {a.x ? b.x : c, a.y ? b.y : c};
  } else if constexpr (N == 3) {
    return {a.x ? b.x : c, a.y ? b.y : c, a.z ? b.z : c};
  } else if constexpr (N == 4) {
    return {a.x ? b.x : c, a.y ? b.y : c, a.z ? b.z : c, a.w ? b.w : c};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> select(const vec<bool, N>& a, T1 b, T2 c) {
  if constexpr (N == 1) {
    return {a.x ? b : c};
  } else if constexpr (N == 2) {
    return {a.x ? b : c, a.y ? b : c};
  } else if constexpr (N == 3) {
    return {a.x ? b : c, a.y ? b : c, a.z ? b : c};
  } else if constexpr (N == 4) {
    return {a.x ? b : c, a.y ? b : c, a.z ? b : c, a.w ? b : c};
  }
}


// ------------------------------------------------------------------

// Vector comparison operations.
template <typename T1, typename T2, size_t N>
constexpr kernel bool operator==(const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return a.x == b.x;
  } else if constexpr (N == 2) {
    return a.x == b.x && a.y == b.y;
  } else if constexpr (N == 3) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
  } else if constexpr (N == 4) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
  } else {
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel bool operator==(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return a.x == b;
  } else if constexpr (N == 2) {
    return a.x == b && a.y == b;
  } else if constexpr (N == 3) {
    return a.x == b && a.y == b && a.z == b;
  } else if constexpr (N == 4) {
    return a.x == b && a.y == b && a.z == b && a.w == b;
  } else {
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel bool operator!=(const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N != 1) {
    return a.x != b.x;
  } else if constexpr (N != 2) {
    return a.x != b.x || a.y != b.y;
  } else if constexpr (N != 3) {
    return a.x != b.x || a.y != b.y || a.z != b.z;
  } else if constexpr (N != 4) {
    return a.x != b.x || a.y != b.y || a.z != b.z || a.w != b.w;
  } else {
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel bool operator!=(const vec<T1, N>& a, T2 b) {
  if constexpr (N != 1) {
    return a.x != b;
  } else if constexpr (N != 2) {
    return a.x != b || a.y != b;
  } else if constexpr (N != 3) {
    return a.x != b || a.y != b || a.z != b;
  } else if constexpr (N != 4) {
    return a.x != b || a.y != b || a.z != b || a.w != b;
  } else {
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel bool operator<(const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return a.x < b.x;
  } else if constexpr (N == 2) {
    return a.x < b.x || (a.x == b.x && a.y < b.y);
  } else if constexpr (N == 3) {
    return a.x < b.x ||
           (a.x == b.x && (a.y < b.y || (a.y == b.y && a.z < b.z)));
  } else if constexpr (N == 4) {
    return a.x < b.x ||
           (a.x == b.x &&
               (a.y < b.y ||
                   (a.y == b.y && (a.z < b.z || (a.z == b.z && a.w < b.w)))));
  } else {
  }
}

// Vector operations.
template <typename T, size_t N>
constexpr kernel vec<T, N> operator+(const vec<T, N>& a) {
  return a;
}
template <typename T, size_t N>
constexpr kernel vec<T, N> operator-(const vec<T, N>& a) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator+(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator+(const vec<T1, N>& a, T2 b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator+(T1 a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator-(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator-(const vec<T1, N>& a, T2 b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator-(T1 a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator*(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator*(const vec<T1, N>& a, T2 b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator*(T1 a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator/(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator/(const vec<T1, N>& a, T2 b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator/(T1 a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator%(const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x % b.x};
  } else if constexpr (N == 2) {
    return {a.x % b.x, a.y % b.y};
  } else if constexpr (N == 3) {
    return {a.x % b.x, a.y % b.y, a.z % b.z};
  } else if constexpr (N == 4) {
    return {a.x % b.x, a.y % b.y, a.z % b.z, a.w % b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator%(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return {a.x % b};
  } else if constexpr (N == 2) {
    return {a.x % b, a.y % b};
  } else if constexpr (N == 3) {
    return {a.x % b, a.y % b, a.z % b};
  } else if constexpr (N == 4) {
    return {a.x % b, a.y % b, a.z % b, a.w % b};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator%(T1 a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x % b.x};
  } else if constexpr (N == 2) {
    return {a.x % b.x, a.y % b.y};
  } else if constexpr (N == 3) {
    return {a.x % b.x, a.y % b.y, a.z % b.z};
  } else if constexpr (N == 4) {
    return {a.x % b.x, a.y % b.y, a.z % b.z, a.w % b.w};
  }
}
template <typename T, size_t N>
constexpr kernel vec<T, N> operator~(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return {~a.x};
  } else if constexpr (N == 2) {
    return {~a.x, ~a.y};
  } else if constexpr (N == 3) {
    return {~a.x, ~a.y, ~a.z};
  } else if constexpr (N == 4) {
    return {~a.x, ~a.y, ~a.z, ~a.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator&(const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x & b.x};
  } else if constexpr (N == 2) {
    return {a.x & b.x, a.y & b.y};
  } else if constexpr (N == 3) {
    return {a.x & b.x, a.y & b.y, a.z & b.z};
  } else if constexpr (N == 4) {
    return {a.x & b.x, a.y & b.y, a.z & b.z, a.w & b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator&(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return {a.x & b};
  } else if constexpr (N == 2) {
    return {a.x & b, a.y & b};
  } else if constexpr (N == 3) {
    return {a.x & b, a.y & b, a.z & b};
  } else if constexpr (N == 4) {
    return {a.x & b, a.y & b, a.z & b, a.w & b};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator&(T1 a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x & b.x};
  } else if constexpr (N == 2) {
    return {a.x & b.x, a.y & b.y};
  } else if constexpr (N == 3) {
    return {a.x & b.x, a.y & b.y, a.z & b.z};
  } else if constexpr (N == 4) {
    return {a.x & b.x, a.y & b.y, a.z & b.z, a.w & b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator|(const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x | b.x};
  } else if constexpr (N == 2) {
    return {a.x | b.x, a.y | b.y};
  } else if constexpr (N == 3) {
    return {a.x | b.x, a.y | b.y, a.z | b.z};
  } else if constexpr (N == 4) {
    return {a.x | b.x, a.y | b.y, a.z | b.z, a.w | b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator|(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return {a.x | b};
  } else if constexpr (N == 2) {
    return {a.x | b, a.y | b};
  } else if constexpr (N == 3) {
    return {a.x | b, a.y | b, a.z | b};
  } else if constexpr (N == 4) {
    return {a.x | b, a.y | b, a.z | b, a.w | b};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator|(T1 a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x | b.x};
  } else if constexpr (N == 2) {
    return {a.x | b.x, a.y | b.y};
  } else if constexpr (N == 3) {
    return {a.x | b.x, a.y | b.y, a.z | b.z};
  } else if constexpr (N == 4) {
    return {a.x | b.x, a.y | b.y, a.z | b.z, a.w | b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator^(const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x ^ b.x};
  } else if constexpr (N == 2) {
    return {a.x ^ b.x, a.y ^ b.y};
  } else if constexpr (N == 3) {
    return {a.x ^ b.x, a.y ^ b.y, a.z ^ b.z};
  } else if constexpr (N == 4) {
    return {a.x ^ b.x, a.y ^ b.y, a.z ^ b.z, a.w ^ b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator^(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return {a.x ^ b};
  } else if constexpr (N == 2) {
    return {a.x ^ b, a.y ^ b};
  } else if constexpr (N == 3) {
    return {a.x ^ b, a.y ^ b, a.z ^ b};
  } else if constexpr (N == 4) {
    return {a.x ^ b, a.y ^ b, a.z ^ b, a.w ^ b};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator^(T1 a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x ^ b.x};
  } else if constexpr (N == 2) {
    return {a.x ^ b.x, a.y ^ b.y};
  } else if constexpr (N == 3) {
    return {a.x ^ b.x, a.y ^ b.y, a.z ^ b.z};
  } else if constexpr (N == 4) {
    return {a.x ^ b.x, a.y ^ b.y, a.z ^ b.z, a.w ^ b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator>>(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x >> b.x};
  } else if constexpr (N == 2) {
    return {a.x >> b.x, a.y >> b.y};
  } else if constexpr (N == 3) {
    return {a.x >> b.x, a.y >> b.y, a.z >> b.z};
  } else if constexpr (N == 4) {
    return {a.x >> b.x, a.y >> b.y, a.z >> b.z, a.w >> b.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator>>(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return {a.x >> b};
  } else if constexpr (N == 2) {
    return {a.x >> b, a.y >> b};
  } else if constexpr (N == 3) {
    return {a.x >> b, a.y >> b, a.z >> b};
  } else if constexpr (N == 4) {
    return {a.x >> b, a.y >> b, a.z >> b, a.w >> b};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator>>(T1 a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x >> b.x};
  } else if constexpr (N == 2) {
    return {a.x >> b.x, a.y >> b.y};
  } else if constexpr (N == 3) {
    return {a.x >> b.x, a.y >> b.y, a.z >> b.z};
  } else if constexpr (N == 4) {
    return {a.x >> b.x, a.y >> b.y, a.z >> b.z, a.w >> b.w};
  }
}
template <typename T, size_t N>
constexpr kernel vec<T, N> operator<<(const vec<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return {a.x << b.x};
  } else if constexpr (N == 2) {
    return {a.x << b.x, a.y << b.y};
  } else if constexpr (N == 3) {
    return {a.x << b.x, a.y << b.y, a.z << b.z};
  } else if constexpr (N == 4) {
    return {a.x << b.x, a.y << b.y, a.z << b.z, a.w << b.w};
  }
}
template <typename T, size_t N>
constexpr kernel vec<T, N> operator<<(const vec<T, N>& a, T b) {
  if constexpr (N == 1) {
    return {a.x << b};
  } else if constexpr (N == 2) {
    return {a.x << b, a.y << b};
  } else if constexpr (N == 3) {
    return {a.x << b, a.y << b, a.z << b};
  } else if constexpr (N == 4) {
    return {a.x << b, a.y << b, a.z << b, a.w << b};
  }
}
template <typename T, size_t N>
constexpr kernel vec<T, N> operator<<(T a, const vec<T, N>& b) {
  if constexpr (N == 1) {
    return {a.x << b.x};
  } else if constexpr (N == 2) {
    return {a.x << b.x, a.y << b.y};
  } else if constexpr (N == 3) {
    return {a.x << b.x, a.y << b.y, a.z << b.z};
  } else if constexpr (N == 4) {
    return {a.x << b.x, a.y << b.y, a.z << b.z, a.w << b.w};
  }
}
template <size_t N>
constexpr kernel vec<bool, N> operator!(const vec<bool, N>& a) {
  if constexpr (N == 1) {
    return {!a.x};
  } else if constexpr (N == 2) {
    return {!a.x, !a.y};
  } else if constexpr (N == 3) {
    return {!a.x, !a.y, !a.z};
  } else if constexpr (N == 4) {
    return {!a.x, !a.y, !a.z, !a.w};
  }
}
template <size_t N>
constexpr kernel vec<bool, N> operator&&(
    const vec<bool, N>& a, const vec<bool, N>& b) {
  if constexpr (N == 1) {
    return {a.x && b.x};
  } else if constexpr (N == 2) {
    return {a.x && b.x, a.y && b.y};
  } else if constexpr (N == 3) {
    return {a.x && b.x, a.y && b.y, a.z && b.z};
  } else if constexpr (N == 4) {
    return {a.x && b.x, a.y && b.y, a.z && b.z, a.w && b.w};
  }
}
template <size_t N>
constexpr kernel vec<bool, N> operator&&(const vec<bool, N>& a, bool b) {
  if constexpr (N == 1) {
    return {a.x && b};
  } else if constexpr (N == 2) {
    return {a.x && b, a.y && b};
  } else if constexpr (N == 3) {
    return {a.x && b, a.y && b, a.z && b};
  } else if constexpr (N == 4) {
    return {a.x && b, a.y && b, a.z && b, a.w && b};
  }
}
template <size_t N>
constexpr kernel vec<bool, N> operator&&(bool a, const vec<bool, N>& b) {
  if constexpr (N == 1) {
    return {a && b.x};
  } else if constexpr (N == 2) {
    return {a && b.x, a && b.y};
  } else if constexpr (N == 3) {
    return {a && b.x, a && b.y, a && b.z};
  } else if constexpr (N == 4) {
    return {a && b.x, a && b.y, a && b.z, a && b.w};
  }
}
template <size_t N>
constexpr kernel vec<bool, N> operator||(
    const vec<bool, N>& a, const vec<bool, N>& b) {
  if constexpr (N == 1) {
    return {a.x || b.x};
  } else if constexpr (N == 2) {
    return {a.x || b.x, a.y || b.y};
  } else if constexpr (N == 3) {
    return {a.x || b.x, a.y || b.y, a.z || b.z};
  } else if constexpr (N == 4) {
    return {a.x || b.x, a.y || b.y, a.z || b.z, a.w || b.w};
  }
}
template <size_t N>
constexpr kernel vec<bool, N> operator||(const vec<bool, N>& a, bool b) {
  if constexpr (N == 1) {
    return {a.x || b};
  } else if constexpr (N == 2) {
    return {a.x || b, a.y || b};
  } else if constexpr (N == 3) {
    return {a.x || b, a.y || b, a.z || b};
  } else if constexpr (N == 4) {
    return {a.x || b, a.y || b, a.z || b, a.w || b};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator||(T1 a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a || b.x};
  } else if constexpr (N == 2) {
    return {a || b.x, a || b.y};
  } else if constexpr (N == 3) {
    return {a || b.x, a || b.y, a || b.z};
  } else if constexpr (N == 4) {
    return {a || b.x, a || b.y, a || b.z, a || b.w};
  }
}

// Vector assignments
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator+=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a + b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator+=(vec<T, N>& a, T1 b) {
  return a = a + b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator-=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a - b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator-=(vec<T, N>& a, T1 b) {
  return a = a - b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator*=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a * b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator*=(vec<T, N>& a, T1 b) {
  return a = a * b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator/=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a / b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator/=(vec<T, N>& a, T1 b) {
  return a = a / b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator%=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a % b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator%=(vec<T, N>& a, T1 b) {
  return a = a % b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator&=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a & b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator&=(vec<T, N>& a, T1 b) {
  return a = a & b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator|=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a | b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator|=(vec<T, N>& a, T1 b) {
  return a = a | b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator^=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a ^ b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator^=(vec<T, N>& a, T1 b) {
  return a = a ^ b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator>>=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a >> b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator>>=(vec<T, N>& a, T1 b) {
  return a = a >> b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator<<=(vec<T, N>& a, const vec<T1, N>& b) {
  return a = a << b;
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N>& operator<<=(vec<T, N>& a, T1 b) {
  return a = a << b;
}

// Vector products and lengths.
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel T dot(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, typename T = common_t<T1, T2>>
constexpr kernel T cross(const vec<T1, 2>& a, const vec<T2, 2>& b) {
  return a.x * b.y - a.y * b.x;
}
template <typename T1, typename T2, typename T = common_t<T1, T2>>
constexpr kernel vec<T, 3> cross(const vec<T1, 3>& a, const vec<T2, 3>& b) {
  return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
template <typename T1, typename T2, typename T = common_t<T1, T2>>
constexpr kernel T angle(const vec<T1, 3>& a, const vec<T2, 3>& b) {
  return acos(clamp(dot(normalize(a), normalize(b)), (T)-1, (T)1));
}

// Orthogonal vectors.
template <typename T>
constexpr kernel vec<T, 3> orthogonal(const vec<T, 3>& v) {
  // http://lolengine.net/blog/2013/09/21/picking-orthogonal-vector-combing-coconuts)
  return abs(v.x) > abs(v.z) ? vec<T, 3>{-v.y, v.x, 0}
                             : vec<T, 3>{0, -v.z, v.y};
}
template <typename T1, typename T2, typename T = common_t<T1, T2>>
constexpr kernel vec<T, 3> orthonormalize(
    const vec<T1, 3>& a, const vec<T2, 3>& b) {
  return normalize(a - b * dot(a, b));
}

// Reflected and refracted vector.
template <typename T1, typename T2, typename T = common_t<T1, T2>>
constexpr kernel vec<T, 3> reflect(const vec<T1, 3>& w, const vec<T2, 3>& n) {
  return -w + 2 * dot(n, w) * n;
}
template <typename T1, typename T2, typename T3,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, 3> refract(
    const vec<T1, 3>& w, const vec<T2, 3>& n, T3 inv_eta) {
  auto cosine = dot(n, w);
  auto k      = 1 + inv_eta * inv_eta * (cosine * cosine - 1);
  if (k < 0) return {0, 0, 0};  // tir
  return -w * inv_eta + (inv_eta * cosine - sqrt(k)) * n;
}

template <typename T, size_t N>
constexpr kernel T length(const vec<T, N>& a) {
  return sqrt(dot(a, a));
}
template <typename T, size_t N>
constexpr kernel T length2(const vec<T, N>& a) {
  return dot(a, a);
}
template <typename T, size_t N>
[[deprecated]] constexpr kernel T length_squared(const vec<T, N>& a) {
  return dot(a, a);
}
template <typename T, size_t N>
constexpr kernel vec<T, N> normalize(const vec<T, N>& a) {
  auto l = length(a);
  return (l != 0) ? a / l : a;
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel T distance(const vec<T1, N>& a, const vec<T2, N>& b) {
  return length(a - b);
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel T distance2(const vec<T1, N>& a, const vec<T2, N>& b) {
  return dot(a - b, a - b);
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
[[deprecated]] constexpr kernel T distance_squared(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  return dot(a - b, a - b);
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel T angle(const vec<T1, N>& a, const vec<T2, N>& b) {
  return acos(clamp(dot(normalize(a), normalize(b)), (T)-1, (T)1));
}
template <typename T1, typename T2, typename T3, size_t N,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, N> slerp(
    const vec<T1, N>& a, const vec<T2, N>& b, T3 u) {
  static_assert(N != 1);
  // https://en.wikipedia.org/wiki/Slerp
  auto an = normalize(a), bn = normalize(b);
  auto d = dot(an, bn);
  if (d < 0) {
    bn = -bn;
    d  = -d;
  }
  if (d > (T)0.9995) return normalize(an + u * (bn - an));
  auto theta = acos(clamp(d, (T)-1, (T)1));
  if (theta == 0) return an;
  return an * (sin(theta * (1 - u)) / sin(theta)) +
         bn * (sin(theta * u) / sin(theta));
}

// Max element and clamp.
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> max(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> max(const vec<T1, N>& a, T2 b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> min(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> min(const vec<T1, N>& a, T2 b) {
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
template <typename T1, typename T2, typename T3, size_t N,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, N> clamp(
    const vec<T1, N>& x, const vec<T2, N>& min, const vec<T3, N>& max) {
  if constexpr (N == 1) {
    return {clamp(x.x, min.x, max.x)};
  } else if constexpr (N == 2) {
    return {clamp(x.x, min.x, max.x), clamp(x.y, min.y, max.y)};
  } else if constexpr (N == 3) {
    return {clamp(x.x, min.x, max.x), clamp(x.y, min.y, max.y),
        clamp(x.z, min.z, max.z)};
  } else if constexpr (N == 4) {
    return {clamp(x.x, min.x, max.x), clamp(x.y, min.y, max.y),
        clamp(x.z, min.z, max.z), clamp(x.w, min.w, max.w)};
  }
}
template <typename T1, typename T2, typename T3, size_t N,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, N> clamp(const vec<T1, N>& x, T2 min, T3 max) {
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
template <typename T1, typename T2, typename T3, size_t N,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, N> clamp(
    const vec<T1, N>& x, const vec<T2, N>& min, T3 max) {
  if constexpr (N == 1) {
    return {clamp(x.x, min.x, max)};
  } else if constexpr (N == 2) {
    return {clamp(x.x, min.x, max), clamp(x.y, min.y, max)};
  } else if constexpr (N == 3) {
    return {
        clamp(x.x, min.x, max), clamp(x.y, min.y, max), clamp(x.z, min.z, max)};
  } else if constexpr (N == 4) {
    return {clamp(x.x, min.x, max), clamp(x.y, min.y, max),
        clamp(x.z, min.z, max), clamp(x.w, min.w, max)};
  }
}
template <typename T1, typename T2, typename T3, size_t N,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, N> clamp(
    const vec<T1, N>& x, T2 min, const vec<T3, N>& max) {
  if constexpr (N == 1) {
    return {clamp(x.x, min, max.x)};
  } else if constexpr (N == 2) {
    return {clamp(x.x, min, max.x), clamp(x.y, min, max.y)};
  } else if constexpr (N == 3) {
    return {clamp(x.x, min, max.x), clamp(x.y, min, max.y),
        clamp(x.z, min.z, max.z)};
  } else if constexpr (N == 4) {
    return {clamp(x.x, min, max.x), clamp(x.y, min, max.y),
        clamp(x.z, min, max.z), clamp(x.w, min, max.w)};
  }
}
template <typename T1, typename T2, typename T3, size_t N,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, N> lerp(
    const vec<T1, N>& a, const vec<T2, N>& b, T3 u) {
  return a * (1 - u) + b * u;
}
template <typename T1, typename T2, typename T3, size_t N,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, N> lerp(
    const vec<T1, N>& a, const vec<T2, N>& b, const vec<T3, N>& u) {
  return a * (1 - u) + b * u;
}

template <typename T, size_t N>
constexpr kernel T max(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel T min(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel size_t argmax(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return 0;
  } else if constexpr (N == 2) {
    return a.x >= a.y ? 0 : 1;
  } else if constexpr (N == 3) {
    return a.x >= a.y ? (a.x >= a.z ? 0 : 2) : (a.y >= a.z ? 1 : 2);
  } else if constexpr (N == 4) {
    if (a.w >= a.x && a.w >= a.y && a.w >= a.z)
      return 3;
    else
      return argmax(xyz(a));
  }
}
template <typename T, size_t N>
constexpr kernel size_t argmin(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return 0;
  } else if constexpr (N == 2) {
    return a.x <= a.y ? 0 : 1;
  } else if constexpr (N == 3) {
    return a.x <= a.y ? (a.x <= a.z ? 0 : 2) : (a.y <= a.z ? 1 : 2);
  } else if constexpr (N == 4) {
    if (a.w <= a.x && a.w <= a.y && a.w <= a.z)
      return 3;
    else
      return argmin(xyz(a));
  }
}
template <typename T, size_t N, typename I = int>
constexpr kernel I find_index(const vec<T, N>& a, T b) {
  if constexpr (N == 1) {
    return a.x == b ? 0 : -1;
  } else if constexpr (N == 2) {
    return a.x == b ? 0 : (a.y == b ? 1 : -1);
  } else if constexpr (N == 3) {
    return a.x == b ? 0 : (a.y == b ? 1 : (a.z == b ? 2 : -1));
  } else if constexpr (N == 4) {
    return a.x == b ? 0 : (a.y == b ? 1 : (a.z == b ? 2 : (a.w == b ? 3 : -1)));
  }
}

template <typename T, size_t N>
constexpr kernel T sum(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel T prod(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return a.x;
  } else if constexpr (N == 2) {
    return a.x * a.y;
  } else if constexpr (N == 3) {
    return a.x * a.y * a.z;
  } else if constexpr (N == 4) {
    return a.x * a.y * a.z * a.w;
  }
}
template <typename T, size_t N>
constexpr kernel T mean(const vec<T, N>& a) {
  return sum(a) / N;
}
template <size_t N>
constexpr kernel bool all(const vec<bool, N>& a) {
  if constexpr (N == 1) {
    return a.x;
  } else if constexpr (N == 2) {
    return a.x && a.y;
  } else if constexpr (N == 3) {
    return a.x && a.y && a.z;
  } else if constexpr (N == 4) {
    return a.x && a.y && a.z && a.w;
  }
}
template <typename T, size_t N>
constexpr kernel T any(const vec<bool, N>& a) {
  if constexpr (N == 1) {
    return a.x;
  } else if constexpr (N == 2) {
    return a.x || a.y;
  } else if constexpr (N == 3) {
    return a.x || a.y || a.z;
  } else if constexpr (N == 4) {
    return a.x || a.y || a.z || a.w;
  }
}

// Functions applied to vector elements
template <typename T, size_t N>
constexpr kernel vec<T, N> abs(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel vec<T, N> sqr(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel vec<T, N> sqrt(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel vec<T, N> exp(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel vec<T, N> log(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel vec<T, N> exp2(const vec<T, N>& a) {
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
template <typename T, size_t N>
constexpr kernel vec<T, N> log2(const vec<T, N>& a) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> pow(const vec<T1, N>& a, const vec<T2, N>& b) {
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
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> pow(const vec<T1, N>& a, T2 b) {
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
template <typename T, size_t N>
constexpr kernel vec<T, N> round(const vec<T, N>& a) {
  if constexpr (N == 1) {
    return {round(a.x)};
  } else if constexpr (N == 2) {
    return {round(a.x), round(a.y)};
  } else if constexpr (N == 3) {
    return {round(a.x), round(a.y), round(a.z)};
  } else if constexpr (N == 4) {
    return {round(a.x), round(a.y), round(a.z), round(a.w)};
  }
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N> fmod(const vec<T, N>& a, const vec<T1, N>& b) {
  if constexpr (N == 1) {
    return {fmod(a.x, b.x)};
  } else if constexpr (N == 2) {
    return {fmod(a.x, b.x), fmod(a.y, b.y)};
  } else if constexpr (N == 3) {
    return {fmod(a.x, b.x), fmod(a.y, b.y), fmod(a.z, b.z)};
  } else if constexpr (N == 4) {
    return {fmod(a.x, b.x), fmod(a.y, b.y), fmod(a.z, b.z), fmod(a.w, b.w)};
  }
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N> fmod(const vec<T, N>& a, T1 b) {
  if constexpr (N == 1) {
    return {fmod(a.x, b)};
  } else if constexpr (N == 2) {
    return {fmod(a.x, b), fmod(a.y, b)};
  } else if constexpr (N == 3) {
    return {fmod(a.x, b), fmod(a.y, b), fmod(a.z, b)};
  } else if constexpr (N == 4) {
    return {fmod(a.x, b), fmod(a.y, b), fmod(a.z, b), fmod(a.w, b)};
  }
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N> mod(const vec<T, N>& a, const vec<T1, N>& b) {
  if constexpr (N == 1) {
    return {mod(a.x, b.x)};
  } else if constexpr (N == 2) {
    return {mod(a.x, b.x), mod(a.y, b.y)};
  } else if constexpr (N == 3) {
    return {mod(a.x, b.x), mod(a.y, b.y), mod(a.z, b.z)};
  } else if constexpr (N == 4) {
    return {mod(a.x, b.x), mod(a.y, b.y), mod(a.z, b.z), mod(a.w, b.w)};
  }
}
template <typename T, typename T1, size_t N>
constexpr kernel vec<T, N> mod(const vec<T, N>& a, T1 b) {
  if constexpr (N == 1) {
    return {mod(a.x, b)};
  } else if constexpr (N == 2) {
    return {mod(a.x, b), mod(a.y, b)};
  } else if constexpr (N == 3) {
    return {mod(a.x, b), mod(a.y, b), mod(a.z, b)};
  } else if constexpr (N == 4) {
    return {mod(a.x, b), mod(a.y, b), mod(a.z, b), mod(a.w, b)};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> bias(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return {bias(a.x, b)};
  } else if constexpr (N == 2) {
    return {bias(a.x, b), bias(a.y, b)};
  } else if constexpr (N == 3) {
    return {bias(a.x, b), bias(a.y, b), bias(a.z, b)};
  } else if constexpr (N == 4) {
    return {bias(a.x, b), bias(a.y, b), bias(a.z, b), bias(a.w, b)};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> gain(const vec<T1, N>& a, T2 b) {
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
// template <typename T, size_t N>
// inline void swap(vec<T, N>& a, vec<T, N>& b) {
//   std::swap(a, b);
// }
template <typename T, size_t N>
constexpr kernel bool isfinite(const vec<T, N>& a) {
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

// Conversion between ranges
template <typename T>
constexpr kernel T unit_to_uv(T a) {
  return (a + 1) / 2;
}
template <typename T, size_t N>
constexpr kernel vec<T, N> unit_to_uv(const vec<T, N>& a) {
  return (a + 1) / 2;
}
template <typename T>
constexpr kernel T uv_to_unit(T a) {
  return a * 2 - 1;
}
template <typename T, size_t N>
constexpr kernel vec<T, N> uv_to_unit(const vec<T, N>& a) {
  return a * 2 - 1;
}

// Conversion between coordinates
template <typename T>
constexpr kernel vec<T, 2> cartesian_to_sphericaluv(const vec<T, 3>& w) {
  auto phi = atan2(w.y, w.x), theta = acos(clamp(w.z, -1, 1));
  return {mod(phi / (2 * (T)pi), 1), theta / (T)pi};
}
template <typename T>
constexpr kernel vec<T, 2> cartesiany_to_sphericaluv(const vec<T, 3>& w) {
  auto phi = atan2(w.z, w.x), theta = acos(clamp(w.y, -1, 1));
  return {mod(phi / (2 * (T)pi), 1), theta / (T)pi};
}
template <typename T>
constexpr kernel vec<T, 3> sphericaluv_to_cartesian(const vec<T, 2>& uv) {
  auto phi = uv.x * 2 * (T)pi, theta = uv.y * (T)pi;
  return {cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)};
}
template <typename T>
constexpr kernel vec<T, 3> sphericaluv_to_cartesiany(const vec<T, 2>& uv) {
  auto phi = uv.x * 2 * (T)pi, theta = uv.y * (T)pi;
  return {cos(phi) * sin(theta), cos(theta), sin(phi) * sin(theta)};
}

// Quaternion operations represented as xi + yj + zk + w
// const auto identity_quat4f = vec4f{0, 0, 0, 1};
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, 4> quat_mul(const vec<T1, 4>& a, T2 b) {
  return {a.x * b, a.y * b, a.z * b, a.w * b};
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, 4> quat_mul(const vec<T1, 4>& a, const vec<T2, 4>& b) {
  return {a.x * b.w + a.w * b.x + a.y * b.w - a.z * b.y,
      a.y * b.w + a.w * b.y + a.z * b.x - a.x * b.z,
      a.z * b.w + a.w * b.z + a.x * b.y - a.y * b.x,
      a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z};
}
template <typename T>
constexpr kernel vec<T, 4> quat_conjugate(const vec<T, 4>& a) {
  return {-a.x, -a.y, -a.z, a.w};
}
template <typename T>
constexpr kernel vec<T, 4> quat_inverse(const vec<T, 4>& a) {
  return quat_conjugate(a) / dot(a, a);
}

// Component-wise comparison operations.
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_equal(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x == b.x};
  } else if constexpr (N == 2) {
    return {a.x == b.x, a.y == b.y};
  } else if constexpr (N == 3) {
    return {a.x == b.x, a.y == b.y, a.z == b.z};
  } else if constexpr (N == 4) {
    return {a.x == b.x, a.y == b.y, a.z == b.z, a.w == b.w};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_equal(const vec<T1, N>& a, T2 b) {
  if constexpr (N == 1) {
    return {a.x == b};
  } else if constexpr (N == 2) {
    return {a.x == b, a.y == b};
  } else if constexpr (N == 3) {
    return {a.x == b, a.y == b, a.z == b};
  } else if constexpr (N == 4) {
    return {a.x == b, a.y == b, a.z == b, a.w == b};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_not_equal(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x != b.x};
  } else if constexpr (N == 2) {
    return {a.x != b.x, a.y != b.y};
  } else if constexpr (N == 3) {
    return {a.x != b.x, a.y != b.y, a.z != b.z};
  } else if constexpr (N == 4) {
    return {a.x != b.x, a.y != b.y, a.z != b.z, a.w != b.w};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_not_equal(
    const vec<T1, N>& a, const T2 b) {
  if constexpr (N == 1) {
    return {a.x != b};
  } else if constexpr (N == 2) {
    return {a.x != b, a.y != b};
  } else if constexpr (N == 3) {
    return {a.x != b, a.y != b, a.z != b};
  } else if constexpr (N == 4) {
    return {a.x != b, a.y != b, a.z != b, a.w != b};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_less(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x < b.x};
  } else if constexpr (N == 2) {
    return {a.x < b.x, a.y < b.y};
  } else if constexpr (N == 3) {
    return {a.x < b.x, a.y < b.y, a.z < b.z};
  } else if constexpr (N == 4) {
    return {a.x < b.x, a.y < b.y, a.z < b.z, a.w < b.w};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_less(const vec<T1, N>& a, const T2 b) {
  if constexpr (N == 1) {
    return {a.x < b};
  } else if constexpr (N == 2) {
    return {a.x < b, a.y < b};
  } else if constexpr (N == 3) {
    return {a.x < b, a.y < b, a.z < b};
  } else if constexpr (N == 4) {
    return {a.x < b, a.y < b, a.z < b, a.w < b};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_greater(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x > b.x};
  } else if constexpr (N == 2) {
    return {a.x > b.x, a.y > b.y};
  } else if constexpr (N == 3) {
    return {a.x > b.x, a.y > b.y, a.z > b.z};
  } else if constexpr (N == 4) {
    return {a.x > b.x, a.y > b.y, a.z > b.z, a.w > b.w};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_greater(
    const vec<T1, N>& a, const T2 b) {
  if constexpr (N == 1) {
    return {a.x > b};
  } else if constexpr (N == 2) {
    return {a.x > b, a.y > b};
  } else if constexpr (N == 3) {
    return {a.x > b, a.y > b, a.z > b};
  } else if constexpr (N == 4) {
    return {a.x > b, a.y > b, a.z > b, a.w > b};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_less_equal(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x <= b.x};
  } else if constexpr (N == 2) {
    return {a.x <= b.x, a.y <= b.y};
  } else if constexpr (N == 3) {
    return {a.x <= b.x, a.y <= b.y, a.z <= b.z};
  } else if constexpr (N == 4) {
    return {a.x <= b.x, a.y <= b.y, a.z <= b.z, a.w <= b.w};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_less_equal(
    const vec<T1, N>& a, const T2 b) {
  if constexpr (N == 1) {
    return {a.x <= b};
  } else if constexpr (N == 2) {
    return {a.x <= b, a.y <= b};
  } else if constexpr (N == 3) {
    return {a.x <= b, a.y <= b, a.z <= b};
  } else if constexpr (N == 4) {
    return {a.x <= b, a.y <= b, a.z <= b, a.w <= b};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_greater_equal(
    const vec<T1, N>& a, const vec<T2, N>& b) {
  if constexpr (N == 1) {
    return {a.x >= b.x};
  } else if constexpr (N == 2) {
    return {a.x >= b.x, a.y >= b.y};
  } else if constexpr (N == 3) {
    return {a.x >= b.x, a.y >= b.y, a.z >= b.z};
  } else if constexpr (N == 4) {
    return {a.x >= b.x, a.y >= b.y, a.z >= b.z, a.w >= b.w};
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel vec<bool, N> component_greater_equal(
    const vec<T1, N>& a, const T2 b) {
  if constexpr (N == 1) {
    return {a.x >= b};
  } else if constexpr (N == 2) {
    return {a.x >= b, a.y >= b};
  } else if constexpr (N == 3) {
    return {a.x >= b, a.y >= b, a.z >= b};
  } else if constexpr (N == 4) {
    return {a.x >= b, a.y >= b, a.z >= b, a.w >= b};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> select(
    const vec<bool, N>& a, const vec<T1, N>& b, const vec<T2, N>& c) {
  if constexpr (N == 1) {
    return {a.x ? b.x : c.x};
  } else if constexpr (N == 2) {
    return {a.x ? b.x : c.x, a.y ? b.y : c.y};
  } else if constexpr (N == 3) {
    return {a.x ? b.x : c.x, a.y ? b.y : c.y, a.z ? b.z : c.z};
  } else if constexpr (N == 4) {
    return {a.x ? b.x : c.x, a.y ? b.y : c.y, a.z ? b.z : c.z, a.w ? b.w : c.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> select(
    const vec<bool, N>& a, T1 b, const vec<T2, N>& c) {
  if constexpr (N == 1) {
    return {a.x ? b : c.x};
  } else if constexpr (N == 2) {
    return {a.x ? b : c.x, a.y ? b : c.y};
  } else if constexpr (N == 3) {
    return {a.x ? b : c.x, a.y ? b : c.y, a.z ? b : c.z};
  } else if constexpr (N == 4) {
    return {a.x ? b : c.x, a.y ? b : c.y, a.z ? b : c.z, a.w ? b : c.w};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> select(
    const vec<bool, N>& a, const vec<T1, N>& b, T2 c) {
  if constexpr (N == 1) {
    return {a.x ? b.x : c};
  } else if constexpr (N == 2) {
    return {a.x ? b.x : c, a.y ? b.y : c};
  } else if constexpr (N == 3) {
    return {a.x ? b.x : c, a.y ? b.y : c, a.z ? b.z : c};
  } else if constexpr (N == 4) {
    return {a.x ? b.x : c, a.y ? b.y : c, a.z ? b.z : c, a.w ? b.w : c};
  }
}
template <typename T1, typename T2, size_t N, typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> select(const vec<bool, N>& a, T1 b, T2 c) {
  if constexpr (N == 1) {
    return {a.x ? b : c};
  } else if constexpr (N == 2) {
    return {a.x ? b : c, a.y ? b : c};
  } else if constexpr (N == 3) {
    return {a.x ? b : c, a.y ? b : c, a.z ? b : c};
  } else if constexpr (N == 4) {
    return {a.x ? b : c, a.y ? b : c, a.z ? b : c, a.w ? b : c};
  }
}

#endif

}  // namespace yocto

// -----------------------------------------------------------------------------
// MATRICES
// -----------------------------------------------------------------------------
namespace yocto {

// Small Fixed-size matrices stored in column major format.
struct mat2x2f {
  vec2f x = {0, 0};
  vec2f y = {0, 0};

  constexpr kernel mat2x2f() : x{0, 0}, y{0, 0} {}
  constexpr kernel mat2x2f(const vec2f& x_, const vec2f& y_) : x{x_}, y{y_} {}

  constexpr kernel vec2f&       operator[](size_t i) { return (&x)[i]; }
  constexpr kernel const vec2f& operator[](size_t i) const { return (&x)[i]; }
};

// Small Fixed-size matrices stored in column major format.
struct mat3x3f {
  vec3f x = {0, 0, 0};
  vec3f y = {0, 0, 0};
  vec3f z = {0, 0, 0};

  constexpr kernel mat3x3f() : x{0, 0, 0}, y{0, 0, 0}, z{0, 0, 0} {}
  constexpr kernel mat3x3f(const vec3f& x_, const vec3f& y_, const vec3f& z_) :
      x{x_}, y{y_}, z{z_} {}

  constexpr kernel vec3f&       operator[](size_t i) { return (&x)[i]; }
  constexpr kernel const vec3f& operator[](size_t i) const { return (&x)[i]; }
};

// Small Fixed-size matrices stored in column major format.
struct mat4x4f {
  vec4f x = {0, 0, 0, 0};
  vec4f y = {0, 0, 0, 0};
  vec4f z = {0, 0, 0, 0};
  vec4f w = {0, 0, 0, 0};

  constexpr kernel mat4x4f() :
      x{0, 0, 0, 0}, y{0, 0, 0, 0}, z{0, 0, 0, 0}, w{0, 0, 0, 0} {}
  constexpr kernel mat4x4f(
      const vec4f& x_, const vec4f& y_, const vec4f& z_, const vec4f& w_) :
      x{x_}, y{y_}, z{z_}, w{w_} {}

  constexpr kernel vec4f&       operator[](size_t i) { return (&x)[i]; }
  constexpr kernel const vec4f& operator[](size_t i) const { return (&x)[i]; }
};

// Matrix aliases
using mat2f = mat2x2f;
using mat3f = mat3x3f;
using mat4f = mat4x4f;

// Identity matrices constants.
inline auto identity2x2f = mat2x2f{{1, 0}, {0, 1}};
inline auto identity3x3f = mat3x3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
inline auto identity4x4f = mat4x4f{
    {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};

// Rows and columns
inline kernel vec2f row(const mat2x2f& a, int i) { return {a.x[i], a.y[i]}; }
inline kernel const vec2f& col(const mat2x2f& a, int j) { return (&a.x)[j]; }

// Data access
inline kernel float*       data(mat2x2f& a) { return &a.x.x; }
inline kernel const float* data(const mat2x2f& a) { return &a.x.x; }

// Matrix comparisons.
inline kernel bool operator==(const mat2x2f& a, const mat2x2f& b) {
  return a.x == b.x && a.y == b.y;
}
inline kernel bool operator!=(const mat2x2f& a, const mat2x2f& b) {
  return !(a == b);
}

// Matrix operations.
inline kernel mat2x2f operator+(const mat2x2f& a, const mat2x2f& b) {
  return {a.x + b.x, a.y + b.y};
}
inline kernel mat2x2f operator-(const mat2x2f& a, const mat2x2f& b) {
  return {a.x - b.x, a.y - b.y};
}
inline kernel mat2x2f operator*(const mat2x2f& a, float b) {
  return {a.x * b, a.y * b};
}
inline kernel mat2x2f operator*(float a, const mat2x2f& b) {
  return {a * b.x, a * b.y};
}
inline kernel mat2x2f operator/(const mat2x2f& a, float b) {
  return {a.x / b, a.y / b};
}
inline kernel vec2f operator*(const mat2x2f& a, const vec2f& b) {
  return a.x * b.x + a.y * b.y;
}
inline kernel vec2f operator*(const vec2f& a, const mat2x2f& b) {
  return {dot(a, b.x), dot(a, b.y)};
}
inline kernel mat2x2f operator*(const mat2x2f& a, const mat2x2f& b) {
  return {a * b.x, a * b.y};
}

// Matrix assignments.
template <typename T, typename T1, size_t N, size_t M>
inline kernel mat2x2f& operator+=(mat2x2f& a, const mat2x2f& b) {
  return a = a + b;
}
template <typename T, typename T1, size_t N, size_t M>
inline kernel mat2x2f& operator-=(mat2x2f& a, const mat2x2f& b) {
  return a = a - b;
}
inline kernel mat2x2f& operator*=(mat2x2f& a, const mat2x2f& b) {
  return a = a * b;
}
inline kernel mat2x2f& operator*=(mat2x2f& a, float b) { return a = a * b; }
inline kernel mat2x2f& operator/=(mat2x2f& a, float b) { return a = a / b; }

// Matrix diagonals and transposes.
inline kernel vec2f   diagonal(const mat2x2f& a) { return {a.x.x, a.y.y}; }
inline kernel mat2x2f transpose(const mat2x2f& a) {
  return {{a.x.x, a.y.x}, {a.x.y, a.y.y}};
}

// Matrix adjoints, determinants and inverses.
inline kernel float   determinant(const mat2x2f& a) { return cross(a.x, a.y); }
inline kernel mat2x2f adjoint(const mat2x2f& a) {
  return {{a.y.y, -a.x.y}, {-a.y.x, a.x.x}};
}
inline kernel mat2x2f inverse(const mat2x2f& a) {
  return adjoint(a) * (1 / determinant(a));
}

// Rows and columns
inline kernel vec3f row(const mat3x3f& a, int i) {
  return {a.x[i], a.y[i], a.z[i]};
}
inline kernel const vec3f& col(const mat3x3f& a, int j) { return (&a.x)[j]; }

// Data access
inline kernel float*       data(mat3x3f& a) { return &a.x.x; }
inline kernel const float* data(const mat3x3f& a) { return &a.x.x; }

// Matrix comparisons.
inline kernel bool operator==(const mat3x3f& a, const mat3x3f& b) {
  return a.x == b.x && a.y == b.y && a.z == b.z;
}
inline kernel bool operator!=(const mat3x3f& a, const mat3x3f& b) {
  return !(a == b);
}

// Matrix operations.
inline kernel mat3x3f operator+(const mat3x3f& a, const mat3x3f& b) {
  return {a.x + b.x, a.y + b.y, a.z + b.z};
}
inline kernel mat3x3f operator-(const mat3x3f& a, const mat3x3f& b) {
  return {a.x - b.x, a.y - b.y, a.z - b.z};
}
inline kernel mat3x3f operator*(const mat3x3f& a, float b) {
  return {a.x * b, a.y * b, a.z * b};
}
inline kernel mat3x3f operator*(float a, const mat3x3f& b) {
  return {a * b.x, a * b.y, a * b.z};
}
inline kernel mat3x3f operator/(const mat3x3f& a, float b) {
  return {a.x / b, a.y / b, a.z / b};
}
inline kernel vec3f operator*(const mat3x3f& a, const vec3f& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline kernel vec3f operator*(const vec3f& a, const mat3x3f& b) {
  return {dot(a, b.x), dot(a, b.y), dot(a, b.z)};
}
inline kernel mat3x3f operator*(const mat3x3f& a, const mat3x3f& b) {
  return {a * b.x, a * b.y, a * b.z};
}

// Matrix assignments.
inline kernel mat3x3f& operator+=(mat3x3f& a, const mat3x3f& b) {
  return a = a + b;
}
inline kernel mat3x3f& operator-=(mat3x3f& a, const mat3x3f& b) {
  return a = a - b;
}
inline kernel mat3x3f& operator*=(mat3x3f& a, const mat3x3f& b) {
  return a = a * b;
}
inline kernel mat3x3f& operator*=(mat3x3f& a, float b) { return a = a * b; }
inline kernel mat3x3f& operator/=(mat3x3f& a, float b) { return a = a / b; }

// Matrix diagonals and transposes.
inline kernel vec3f diagonal(const mat3x3f& a) { return {a.x.x, a.y.y, a.z.z}; }
inline kernel mat3x3f transpose(const mat3x3f& a) {
  return {
      {a.x.x, a.y.x, a.z.x},
      {a.x.y, a.y.y, a.z.y},
      {a.x.z, a.y.z, a.z.z},
  };
}

// Matrix adjoints, determinants and inverses.
inline kernel float determinant(const mat3x3f& a) {
  return dot(a.x, cross(a.y, a.z));
}
inline kernel mat3x3f adjoint(const mat3x3f& a) {
  return transpose(mat3x3f{cross(a.y, a.z), cross(a.z, a.x), cross(a.x, a.y)});
}
inline kernel mat3x3f inverse(const mat3x3f& a) {
  return adjoint(a) * (1 / determinant(a));
}

// Constructs a basis from a direction
inline kernel mat3x3f basis_fromz(const vec3f& v) {
  // https://graphics.pixar.com/library/OrthonormalB/paper.pdf
  auto z    = normalize(v);
  auto sign = copysignf(1.0f, z.z);
  auto a    = -1.0f / (sign + z.z);
  auto b    = z.x * z.y * a;
  auto x    = vec3f{1.0f + sign * z.x * z.x * a, sign * b, -sign * z.x};
  auto y    = vec3f{b, sign + z.y * z.y * a, -z.y};
  return {x, y, z};
}

// Rows and columns
inline kernel vec4f row(const mat4x4f& a, int i) {
  return {a.x[i], a.y[i], a.z[i], a.w[i]};
}
inline kernel const vec4f& col(const mat4x4f& a, int j) { return (&a.x)[j]; }

// Data access
inline kernel float*       data(mat4x4f& a) { return &a.x.x; }
inline kernel const float* data(const mat4x4f& a) { return &a.x.x; }

// Matrix comparisons.
inline kernel bool operator==(const mat4x4f& a, const mat4x4f& b) {
  return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
}
inline kernel bool operator!=(const mat4x4f& a, const mat4x4f& b) {
  return !(a == b);
}

// Matrix operations.
inline kernel mat4x4f operator+(const mat4x4f& a, const mat4x4f& b) {
  return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}
inline kernel mat4x4f operator-(const mat4x4f& a, const mat4x4f& b) {
  return {a.x - b.x, a.y - b.y, a.z - b.z, a.w - b.w};
}
inline kernel mat4x4f operator*(const mat4x4f& a, float b) {
  return {a.x * b, a.y * b, a.z * b, a.w * b};
}
inline kernel mat4x4f operator*(float a, const mat4x4f& b) {
  return {a * b.x, a * b.y, a * b.z, a * b.w};
}
inline kernel mat4x4f operator/(const mat4x4f& a, float b) {
  return {a.x / b, a.y / b, a.z / b, a.w / b};
}
inline kernel vec4f operator*(const mat4x4f& a, const vec4f& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}
inline kernel vec4f operator*(const vec4f& a, const mat4x4f& b) {
  return {dot(a, b.x), dot(a, b.y), dot(a, b.z), dot(a, b.w)};
}
inline kernel mat4x4f operator*(const mat4x4f& a, const mat4x4f& b) {
  return {a * b.x, a * b.y, a * b.z, a * b.w};
}

// Matrix assignments.
inline kernel mat4x4f& operator+=(mat4x4f& a, const mat4x4f& b) {
  return a = a + b;
}
inline kernel mat4x4f& operator-=(mat4x4f& a, const mat4x4f& b) {
  return a = a - b;
}
inline kernel mat4x4f& operator*=(mat4x4f& a, const mat4x4f& b) {
  return a = a * b;
}
inline kernel mat4x4f& operator*=(mat4x4f& a, float b) { return a = a * b; }
inline kernel mat4x4f& operator/=(mat4x4f& a, float b) { return a = a / b; }

// Matrix diagonals and transposes.
inline kernel vec4f diagonal(const mat4x4f& a) {
  return {a.x.x, a.y.y, a.z.z, a.w.w};
}
inline kernel mat4x4f transpose(const mat4x4f& a) {
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

  constexpr kernel frame2f() : x{1, 0}, y{0, 1}, o{0, 0} {}
  constexpr kernel frame2f(const vec2f& x_, const vec2f& y_, const vec2f& o_) :
      x{x_}, y{y_}, o{o_} {}
  constexpr kernel frame2f(const mat2x2f& xy_, const vec2f& o_) :
      x{xy_.x}, y{xy_.y}, o{o_} {}

  explicit constexpr kernel frame2f(const mat3x3f& m) :
      x{xy(m.x)}, y{xy(m.y)}, o{xy(m.z)} {}
  explicit constexpr kernel operator mat3x3f() const {
    return {{x, 0}, {y, 0}, {o, 1}};
  }

  constexpr kernel vec2f&       operator[](size_t i) { return (&x)[i]; }
  constexpr kernel const vec2f& operator[](size_t i) const { return (&x)[i]; }
};

// Rigid frames stored as a column-major affine transform matrix.
struct frame3f {
  vec3f x = {1, 0, 0};
  vec3f y = {0, 0, 1};
  vec3f z = {0, 0, 1};
  vec3f o = {0, 0, 0};

  constexpr kernel frame3f() : x{1, 0, 0}, y{0, 1, 0}, z{0, 0, 1}, o{0, 0, 0} {}
  constexpr kernel frame3f(
      const vec3f& x_, const vec3f& y_, const vec3f& z_, const vec3f& o_) :
      x{x_}, y{y_}, z{z_}, o{o_} {}
  constexpr kernel frame3f(const mat3x3f& xyz_, const vec3f& o_) :
      x{xyz_.x}, y{xyz_.y}, z{xyz_.z}, o{o_} {}

  explicit constexpr kernel frame3f(const mat4x4f& m) :
      x{xyz(m.x)}, y{xyz(m.y)}, z{xyz(m.z)}, o{xyz(m.w)} {}
  explicit constexpr kernel operator mat4x4f() const {
    return {{x, 0}, {y, 0}, {z, 0}, {o, 1}};
  }

  constexpr kernel vec3f&       operator[](size_t i) { return (&x)[i]; }
  constexpr kernel const vec3f& operator[](size_t i) const { return (&x)[i]; }
};

// Identity frames.
constexpr auto identity2x3f = frame2f{{1, 0}, {0, 1}, {0, 0}};
constexpr auto identity3x4f = frame3f{
    {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};

// Frame properties
inline kernel const mat2f& rotation(const frame2f& a) {
  return (const mat2f&)a;
}
inline kernel const vec2f& translation(const frame2f& a) { return a.o; }

// Frame/mat conversion
inline kernel frame2f to_frame(const mat3f& m) {
  return {xy(m.x), xy(m.y), xy(m.z)};
}
inline kernel mat3x3f to_mat(const frame2f& f) {
  return {{f.x, 0}, {f.y, 0}, {f.o, 1}};
}

// Frame comparisons.
inline kernel bool operator==(const frame2f& a, const frame2f& b) {
  return a.x == b.x && a.y == b.y && a.o == b.o;
}
inline kernel bool operator!=(const frame2f& a, const frame2f& b) {
  return !(a == b);
}

// Frame composition, equivalent to affine matrix product.
inline kernel frame2f operator*(const frame2f& a, const frame2f& b) {
  auto &ma = rotation(a), &mb = rotation(b);
  return {ma * mb, ma * b.o + a.o};
}
inline kernel frame2f& operator*=(frame2f& a, const frame2f& b) {
  return a = a * b;
}

// Frame inverse, equivalent to rigid affine inverse.
inline kernel frame2f inverse(const frame2f& a, bool non_rigid = false) {
  if (non_rigid) {
    auto minv = inverse(rotation(a));
    return {minv, -(minv * a.o)};
  } else {
    auto minv = transpose(rotation(a));
    return {minv, -(minv * a.o)};
  }
}

// Frame properties
inline kernel const mat3x3f& rotation(const frame3f& a) {
  return (const mat3x3f&)a;
}
inline kernel const vec3f& translation(const frame3f& a) { return a.o; }

// Frame/mat conversion
inline kernel frame3f to_frame(const mat4f& m) {
  return {xyz(m.x), xyz(m.y), xyz(m.z), xyz(m.w)};
}
inline kernel mat4x4f to_mat(const frame3f& f) {
  return {{f.x, 0}, {f.y, 0}, {f.z, 0}, {f.o, 1}};
}

// Frame comparisons.
inline kernel bool operator==(const frame3f& a, const frame3f& b) {
  return a.x == b.x && a.y == b.y && a.z == b.z && a.o == b.o;
}
inline kernel bool operator!=(const frame3f& a, const frame3f& b) {
  return !(a == b);
}

// Frame composition, equivalent to affine matrix product.
inline kernel frame3f operator*(const frame3f& a, const frame3f& b) {
  auto &ma = rotation(a), &mb = rotation(b);
  return {ma * mb, ma * b.o + a.o};
}
inline kernel frame3f& operator*=(frame3f& a, const frame3f& b) {
  return a = a * b;
}

// Frame inverse, equivalent to rigid affine inverse.
inline kernel frame3f inverse(const frame3f& a, bool non_rigid = false) {
  if (non_rigid) {
    auto minv = inverse(rotation(a));
    return {minv, -(minv * a.o)};
  } else {
    auto minv = transpose(rotation(a));
    return {minv, -(minv * a.o)};
  }
}

// Frame construction from axis.
inline kernel frame3f frame_fromz(const vec3f& o, const vec3f& v) {
  // https://graphics.pixar.com/library/OrthonormalB/paper.pdf
  auto z    = normalize(v);
  auto sign = copysignf(1.0f, z.z);
  auto a    = -1 / (sign + z.z);
  auto b    = z.x * z.y * a;
  auto x    = vec3f{1 + sign * z.x * z.x * a, sign * b, -sign * z.x};
  auto y    = vec3f{b, sign + z.y * z.y * a, -z.y};
  return {x, y, z, o};
}
inline kernel frame3f frame_fromzx(
    const vec3f& o, const vec3f& z_, const vec3f& x_) {
  auto z = normalize(z_);
  auto x = orthonormalize(x_, z);
  auto y = normalize(cross(z, x));
  return {x, y, z, o};
}
inline kernel frame3f orthonormalize(const frame3f& frame_) {
  auto z = normalize(frame_.z);
  auto x = orthonormalize(frame_.x, z);
  auto y = normalize(cross(z, x));
  auto o = frame_.o;
  return frame3f{x, y, z, o};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// QUATERNIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Quaternions to represent rotations
struct quat4f {
  float x = 0, y = 0, z = 0, w = 0;

  constexpr kernel quat4f() : x{0}, y{0}, z{0}, w{1} {}
  constexpr kernel quat4f(float x_, float y_, float z_, float w_) :
      x{x_}, y{y_}, z{z_}, w{w_} {}

  constexpr kernel float&       operator[](size_t i) { return (&x)[i]; }
  constexpr kernel const float& operator[](size_t i) const { return (&x)[i]; }
};

constexpr auto identityq4f = quat4f{0, 0, 0, 1};

// Quaternion operations
inline kernel quat4f operator+(const quat4f& a, const quat4f& b) {
  return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}
inline kernel quat4f operator*(const quat4f& a, float b) {
  return {a.x * b, a.y * b, a.z * b, a.w * b};
}
inline kernel quat4f operator/(const quat4f& a, float b) {
  return {a.x / b, a.y / b, a.z / b, a.w / b};
}
inline kernel quat4f operator*(const quat4f& a, const quat4f& b) {
  return {a.x * b.w + a.w * b.x + a.y * b.w - a.z * b.y,
      a.y * b.w + a.w * b.y + a.z * b.x - a.x * b.z,
      a.z * b.w + a.w * b.z + a.x * b.y - a.y * b.x,
      a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z};
}

// Quaternion operations
inline kernel float dot(const quat4f& a, const quat4f& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}
inline kernel float  length(const quat4f& a) { return sqrt(dot(a, a)); }
inline kernel quat4f normalize(const quat4f& a) {
  auto l = length(a);
  return (l != 0) ? a / l : a;
}
inline kernel quat4f conjugate(const quat4f& a) {
  return {-a.x, -a.y, -a.z, a.w};
}
inline kernel quat4f inverse(const quat4f& a) {
  return conjugate(a) / dot(a, a);
}
inline kernel float uangle(const quat4f& a, const quat4f& b) {
  auto d = dot(a, b);
  return d > 1 ? 0 : acos(d < -1 ? -1 : d);
}
inline kernel quat4f lerp(const quat4f& a, const quat4f& b, float t) {
  return a * (1 - t) + b * t;
}
inline kernel quat4f nlerp(const quat4f& a, const quat4f& b, float t) {
  return normalize(lerp(a, b, t));
}
inline kernel quat4f slerp(const quat4f& a, const quat4f& b, float t) {
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
inline kernel vec2f transform_point(const mat3x3f& a, const vec2f& b) {
  auto tvb = a * vec3f{b.x, b.y, 1};
  return vec2f{tvb.x, tvb.y} / tvb.z;
}
inline kernel vec2f transform_vector(const mat3x3f& a, const vec2f& b) {
  auto tvb = a * vec3f{b, 0};
  return vec2f{tvb.x, tvb.y} / tvb.z;
}
inline kernel vec2f transform_direction(const mat3x3f& a, const vec2f& b) {
  return normalize(transform_vector(a, b));
}
inline kernel vec2f transform_normal(const mat3x3f& a, const vec2f& b) {
  return normalize(transform_vector(transpose(inverse(a)), b));
}
inline kernel vec2f transform_vector(const mat2x2f& a, const vec2f& b) {
  return a * b;
}
inline kernel vec2f transform_direction(const mat2x2f& a, const vec2f& b) {
  return normalize(transform_vector(a, b));
}
inline kernel vec2f transform_normal(const mat2x2f& a, const vec2f& b) {
  return normalize(transform_vector(transpose(inverse(a)), b));
}

// Transforms points, vectors and directions by matrices.
inline kernel vec3f transform_point(const mat4x4f& a, const vec3f& b) {
  auto tvb = a * vec4f{b, 1};
  return vec3f{tvb.x, tvb.y, tvb.z} / tvb.w;
}
inline kernel vec3f transform_vector(const mat4x4f& a, const vec3f& b) {
  auto tvb = a * vec4f{b, 0};
  return vec3f{tvb.x, tvb.y, tvb.z} / tvb.w;
}
inline kernel vec3f transform_direction(const mat4x4f& a, const vec3f& b) {
  return normalize(transform_vector(a, b));
}
inline kernel vec3f transform_vector(const mat3x3f& a, const vec3f& b) {
  return a * b;
}
inline kernel vec3f transform_direction(const mat3x3f& a, const vec3f& b) {
  return normalize(transform_vector(a, b));
}
inline kernel vec3f transform_normal(const mat3x3f& a, const vec3f& b) {
  return normalize(transform_vector(transpose(inverse(a)), b));
}

// Transforms points, vectors and directions by frames.
inline kernel vec2f transform_point(const frame2f& a, const vec2f& b) {
  return a.x * b.x + a.y * b.y + a.o;
}
inline kernel vec2f transform_vector(const frame2f& a, const vec2f& b) {
  return a.x * b.x + a.y * b.y;
}
inline kernel vec2f transform_direction(const frame2f& a, const vec2f& b) {
  return normalize(transform_vector(a, b));
}
inline kernel vec2f transform_normal(
    const frame2f& a, const vec2f& b, bool non_rigid = false) {
  if (non_rigid) {
    return transform_normal(rotation(a), b);
  } else {
    return normalize(transform_vector(a, b));
  }
}

// Transforms points, vectors and directions by frames.
inline kernel vec3f transform_point(const frame3f& a, const vec3f& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z + a.o;
}
inline kernel vec3f transform_vector(const frame3f& a, const vec3f& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline kernel vec3f transform_direction(const frame3f& a, const vec3f& b) {
  return normalize(transform_vector(a, b));
}
inline kernel vec3f transform_normal(
    const frame3f& a, const vec3f& b, bool non_rigid = false) {
  if (non_rigid) {
    return transform_normal(rotation(a), b);
  } else {
    return normalize(transform_vector(a, b));
  }
}

// Transforms points, vectors and directions by frames.
inline kernel vec2f transform_point_inverse(const frame2f& a, const vec2f& b) {
  return {dot(a.x, b - a.o), dot(a.y, b - a.o)};
}
inline kernel vec2f transform_vector_inverse(const frame2f& a, const vec2f& b) {
  return {dot(a.x, b), dot(a.y, b)};
}
inline kernel vec2f transform_direction_inverse(
    const frame2f& a, const vec2f& b) {
  return normalize(transform_vector_inverse(a, b));
}

// Transforms points, vectors and directions by frames.
inline kernel vec3f transform_point_inverse(const frame3f& a, const vec3f& b) {
  return {dot(a.x, b - a.o), dot(a.y, b - a.o), dot(a.z, b - a.o)};
}
inline kernel vec3f transform_vector_inverse(const frame3f& a, const vec3f& b) {
  return {dot(a.x, b), dot(a.y, b), dot(a.z, b)};
}
inline kernel vec3f transform_direction_inverse(
    const frame3f& a, const vec3f& b) {
  return normalize(transform_vector_inverse(a, b));
}

// Translation, scaling and rotations transforms.
inline kernel frame3f translation_frame(const vec3f& a) {
  return {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, a};
}
inline kernel frame3f scaling_frame(const vec3f& a) {
  return {{a.x, 0, 0}, {0, a.y, 0}, {0, 0, a.z}, {0, 0, 0}};
}
inline kernel frame3f scaling_frame(float a) {
  return scaling_frame(vec3f{a, a, a});
}
inline kernel frame3f rotation_frame(const vec3f& axis, float angle) {
  auto s = sin(angle), c = cos(angle);
  auto vv = normalize(axis);
  return {
      {c + (1 - c) * vv.x * vv.x, (1 - c) * vv.x * vv.y + s * vv.z,
          (1 - c) * vv.x * vv.z - s * vv.y},
      {(1 - c) * vv.x * vv.y - s * vv.z, c + (1 - c) * vv.y * vv.y,
          (1 - c) * vv.y * vv.z + s * vv.x},
      {(1 - c) * vv.x * vv.z + s * vv.y, (1 - c) * vv.y * vv.z - s * vv.x,
          c + (1 - c) * vv.z * vv.z},
      {0, 0, 0},
  };
}
inline kernel frame3f rotation_frame(const vec4f& quat) {
  auto v = quat;
  return {
      {v.w * v.w + v.x * v.x - v.y * v.y - v.z * v.z,
          (v.x * v.y + v.z * v.w) * 2, (v.z * v.x - v.y * v.w) * 2},
      {(v.x * v.y - v.z * v.w) * 2,
          v.w * v.w - v.x * v.x + v.y * v.y - v.z * v.z,
          (v.y * v.z + v.x * v.w) * 2},
      {(v.z * v.x + v.y * v.w) * 2, (v.y * v.z - v.x * v.w) * 2,
          v.w * v.w - v.x * v.x - v.y * v.y + v.z * v.z},
      {0, 0, 0},
  };
}
inline kernel frame3f rotation_frame(const quat4f& quat) {
  auto v = quat;
  return {
      {v.w * v.w + v.x * v.x - v.y * v.y - v.z * v.z,
          (v.x * v.y + v.z * v.w) * 2, (v.z * v.x - v.y * v.w) * 2},
      {(v.x * v.y - v.z * v.w) * 2,
          v.w * v.w - v.x * v.x + v.y * v.y - v.z * v.z,
          (v.y * v.z + v.x * v.w) * 2},
      {(v.z * v.x + v.y * v.w) * 2, (v.y * v.z - v.x * v.w) * 2,
          v.w * v.w - v.x * v.x - v.y * v.y + v.z * v.z},
      {0, 0, 0},
  };
}
inline kernel frame3f rotation_frame(const mat3x3f& rot) {
  return {rot, {0, 0, 0}};
}

// Lookat frame. Z-axis can be inverted with inv_xz.
inline kernel frame3f lookat_frame(const vec3f& eye, const vec3f& center,
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
inline kernel mat4x4f frustum_mat(
    float l, float r, float b, float t, float n, float f) {
  return {{2 * n / (r - l), 0, 0, 0}, {0, 2 * n / (t - b), 0, 0},
      {(r + l) / (r - l), (t + b) / (t - b), -(f + n) / (f - n), -1},
      {0, 0, -2 * f * n / (f - n), 0}};
}
inline kernel mat4x4f ortho_mat(
    float l, float r, float b, float t, float n, float f) {
  return {{2 / (r - l), 0, 0, 0}, {0, 2 / (t - b), 0, 0},
      {0, 0, -2 / (f - n), 0},
      {-(r + l) / (r - l), -(t + b) / (t - b), -(f + n) / (f - n), 1}};
}
inline kernel mat4x4f ortho2d_mat(
    float left, float right, float bottom, float top) {
  return ortho_mat(left, right, bottom, top, -1, 1);
}
inline kernel mat4x4f ortho_mat(float xmag, float ymag, float near, float far) {
  return {{1 / xmag, 0, 0, 0}, {0, 1 / ymag, 0, 0}, {0, 0, 2 / (near - far), 0},
      {0, 0, (far + near) / (near - far), 1}};
}
inline kernel mat4x4f perspective_mat(
    float fovy, float aspect, float near, float far) {
  auto tg = tan(fovy / 2);
  return {{1 / (aspect * tg), 0, 0, 0}, {0, 1 / tg, 0, 0},
      {0, 0, (far + near) / (near - far), -1},
      {0, 0, 2 * far * near / (near - far), 0}};
}
inline kernel mat4x4f perspective_mat(float fovy, float aspect, float near) {
  auto tg = tan(fovy / 2);
  return {{1 / (aspect * tg), 0, 0, 0}, {0, 1 / tg, 0, 0}, {0, 0, -1, -1},
      {0, 0, 2 * near, 0}};
}

// Rotation conversions.
inline kernel pair<vec3f, float> rotation_axisangle(const vec4f& quat) {
  return {normalize(vec3f{quat.x, quat.y, quat.z}), 2 * acos(quat.w)};
}
inline kernel vec4f rotation_quat(const vec3f& axis, float angle) {
  auto len = length(axis);
  if (len == 0) return {0, 0, 0, 1};
  return vec4f{sin(angle / 2) * axis.x / len, sin(angle / 2) * axis.y / len,
      sin(angle / 2) * axis.z / len, cos(angle / 2)};
}
inline kernel vec4f rotation_quat(const vec4f& axisangle) {
  return rotation_quat(
      vec3f{axisangle.x, axisangle.y, axisangle.z}, axisangle.w);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// USER INTERFACE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Computes the aspect ratio.
inline kernel float aspect_ratio(const vec2i& size) {
  return (float)size.x / (float)size.y;
}

// Flip u from [0,1] to [1,0]
inline kernel vec2f flip_u(const vec2f& uv) { return {1 - uv.x, uv.y}; }
// Flip v from [0,1] to [1,0]
inline kernel vec2f flip_v(const vec2f& uv) { return {uv.x, 1 - uv.y}; }

// Computes the image uv coordinates corresponding to the view parameters.
// Returns negative coordinates if out of the image.
inline kernel vec2i image_coords(const vec2f& mouse_pos, const vec2f& center,
    float scale, const vec2i& size, bool clamped = true) {
  auto xy = (mouse_pos - center) / scale;
  auto ij = (vec2i)round(xy + size / 2.0f);
  return clamped ? clamp(ij, 0, size) : ij;
}

// Center image and autofit. Returns center and scale.
inline kernel pair<vec2f, float> camera_imview(const vec2f& center, float scale,
    const vec2i& imsize, const vec2i& winsize, bool zoom_to_fit) {
  if (zoom_to_fit) {
    return {(vec2f)winsize / 2, min(winsize / (vec2f)imsize)};
  } else {
    return {select(component_greater_equal(winsize, imsize * scale),
                (vec2f)winsize / 2, center),
        scale};
  }
}

// Turntable for UI navigation. Returns from and to.
inline kernel pair<vec3f, vec3f> camera_turntable(const vec3f& from_,
    const vec3f& to_, const vec3f& up, const vec2f& rotate, float dolly,
    const vec2f& pan) {
  // copy values
  auto from = from_, to = to_;

  // rotate if necessary
  if (rotate != vec2f{0, 0}) {
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
  if (dolly != 0) {
    auto z  = normalize(to - from);
    auto lz = max(0.001f, length(to - from) * (1 + dolly));
    z *= lz;
    from = to - z;
  }

  // pan if necessary
  if (pan != vec2f{0, 0}) {
    auto z = normalize(to - from);
    auto x = normalize(cross(up, z));
    auto y = normalize(cross(z, x));
    auto t = vec3f{pan.x * x.x + pan.y * y.x, pan.x * x.y + pan.y * y.y,
        pan.x * x.z + pan.y * y.z};
    from += t;
    to += t;
  }

  // done
  return {from, to};
}

// Turntable for UI navigation. Returns frame and focus.
inline kernel pair<frame3f, float> camera_turntable(const frame3f& frame_,
    float focus, const vec2f& rotate, float dolly, const vec2f& pan) {
  // copy values
  auto frame = frame_;

  // rotate if necessary
  if (rotate != vec2f{0, 0}) {
    auto phi   = atan2(frame.z.z, frame.z.x) + rotate.x;
    auto theta = acos(frame.z.y) + rotate.y;
    theta      = clamp(theta, 0.001f, pif - 0.001f);
    auto new_z = vec3f{
        sin(theta) * cos(phi), cos(theta), sin(theta) * sin(phi)};
    auto new_center = frame.o - frame.z * focus;
    auto new_o      = new_center + new_z * focus;
    frame           = lookat_frame(new_o, new_center, vec3f{0, 1, 0});
    focus           = length(new_o - new_center);
  }

  // pan if necessary
  if (dolly != 0) {
    auto c  = frame.o - frame.z * focus;
    focus   = max(focus * (1 + dolly), 0.001f);
    frame.o = c + frame.z * focus;
  }

  // pan if necessary
  if (pan != vec2f{0, 0}) {
    frame.o += frame.x * pan.x + frame.y * pan.y;
  }

  // done
  return {frame, focus};
}

// FPS camera for UI navigation for a frame parametrization. Returns frame.
inline kernel frame3f camera_fpscam(
    const frame3f& frame, const vec3f& transl, const vec2f& rotate) {
  // https://gamedev.stackexchange.com/questions/30644/how-to-keep-my-quaternion-using-fps-camera-from-tilting-and-messing-up
  auto y = vec3f{0, 1, 0};
  auto z = orthonormalize(frame.z, y);
  auto x = cross(y, z);

  auto rot = rotation_frame(vec3f{1, 0, 0}, rotate.y) *
             yocto::frame3f{frame.x, frame.y, frame.z, vec3f{0, 0, 0}} *
             rotation_frame(vec3f{0, 1, 0}, rotate.x);
  auto pos = frame.o + transl.x * x + transl.y * y + transl.z * z;

  return {rot.x, rot.y, rot.z, pos};
}

// Computes the image uv coordinates corresponding to the view parameters.
// Returns negative coordinates if out of the image.
inline kernel vec2i get_image_coords(const vec2f& mouse_pos,
    const vec2f& center, float scale, const vec2i& txt_size) {
  auto xy = (mouse_pos - center) / scale;
  return (vec2i)round(xy + txt_size / 2.0f);
}

// Center image and autofit.
inline kernel void update_imview(vec2f& center, float& scale,
    const vec2i& imsize, const vec2i& winsize, bool zoom_to_fit) {
  if (zoom_to_fit) {
    scale  = min((vec2f)winsize / imsize);
    center = (vec2f)winsize / 2;
  } else {
    if (winsize.x >= imsize.x * scale) center.x = winsize.x / 2.0f;
    if (winsize.y >= imsize.y * scale) center.y = winsize.y / 2.0f;
  }
}

// Turntable for UI navigation.
inline kernel void update_turntable(vec3f& from, vec3f& to, vec3f& up,
    const vec2f& rotate, float dolly, const vec2f& pan) {
  // rotate if necessary
  if (rotate != vec2f{0, 0}) {
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
  if (dolly != 0) {
    auto z  = normalize(to - from);
    auto lz = max(0.001f, length(to - from) * (1 + dolly));
    z *= lz;
    from = to - z;
  }

  // pan if necessary
  if (pan != vec2f{0, 0}) {
    auto z = normalize(to - from);
    auto x = normalize(cross(up, z));
    auto y = normalize(cross(z, x));
    auto t = pan.x * x + pan.y * y;
    from += t;
    to += t;
  }
}

// Turntable for UI navigation.
inline kernel void update_turntable(frame3f& frame, float& focus,
    const vec2f& rotate, float dolly, const vec2f& pan) {
  // rotate if necessary
  if (rotate != vec2f{0, 0}) {
    auto phi   = atan2(frame.z.z, frame.z.x) + rotate.x;
    auto theta = acos(frame.z.y) + rotate.y;
    theta      = clamp(theta, 0.001f, pif - 0.001f);
    auto new_z = vec3f{
        sin(theta) * cos(phi), cos(theta), sin(theta) * sin(phi)};
    auto new_center = frame.o - frame.z * focus;
    auto new_o      = new_center + new_z * focus;
    frame           = lookat_frame(new_o, new_center, {0, 1, 0});
    focus           = length(new_o - new_center);
  }

  // pan if necessary
  if (dolly != 0) {
    auto c  = frame.o - frame.z * focus;
    focus   = max(focus * (1 + dolly), 0.001f);
    frame.o = c + frame.z * focus;
  }

  // pan if necessary
  if (pan != vec2f{0, 0}) {
    frame.o += frame.x * pan.x + frame.y * pan.y;
  }
}

// FPS camera for UI navigation for a frame parametrization.
inline kernel void update_fpscam(
    frame3f& frame, const vec3f& transl, const vec2f& rotate) {
  // https://gamedev.stackexchange.com/questions/30644/how-to-keep-my-quaternion-using-fps-camera-from-tilting-and-messing-up
  auto y = vec3f{0, 1, 0};
  auto z = orthonormalize(frame.z, y);
  auto x = cross(y, z);

  auto rot = rotation_frame(vec3f{1, 0, 0}, rotate.y) *
             yocto::frame3f{frame.x, frame.y, frame.z, vec3f{0, 0, 0}} *
             rotation_frame(vec3f{0, 1, 0}, rotate.x);
  auto pos = frame.o + transl.x * x + transl.y * y + transl.z * z;

  frame = {rot.x, rot.y, rot.z, pos};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// SIGNED-SIZE
// -----------------------------------------------------------------------------
namespace yocto {

template <typename Sequence>
constexpr kernel std::ptrdiff_t ssize(const Sequence& sequence) {
  return (std::ptrdiff_t)std::size(sequence);
}
template <typename Sequence>
constexpr kernel int isize(const Sequence& sequence) {
  return (int)std::size(sequence);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IOSTREAM SUPPORT
// -----------------------------------------------------------------------------
namespace yocto {

template <typename C>
inline std::basic_istream<C>& operator>>(
    std::basic_istream<C>& s, const char a) {
  constexpr auto max_size = std::numeric_limits<ptrdiff_t>::max();
  return s.ignore(max_size, a);
}
template <typename C>
inline std::basic_istream<C>& operator>>(
    std::basic_istream<C>& s, const char* a) {
  constexpr auto max_size = std::numeric_limits<ptrdiff_t>::max();
  return s.ignore(max_size, a[0]);
}

template <typename C, typename T, size_t N>
inline std::basic_ostream<C>& operator<<(
    std::basic_ostream<C>& s, const vec<T, N>& a) {
  if constexpr (N == 1) {
    return s << "[" << a.x << "]";
  } else if constexpr (N == 2) {
    return s << "[" << a.x << "," << a.y << "]";
  } else if constexpr (N == 3) {
    return s << "[" << a.x << "," << a.y << "," << a.z << "]";
  } else if constexpr (N == 4) {
    return s << "[" << a.x << "," << a.y << "," << a.z << "," << a.w << "]";
  } else {
  }
}

template <typename C, typename T, size_t N>
inline std::basic_istream<C>& operator>>(
    std::basic_istream<C>& s, vec<T, N>& a) {
  if constexpr (N == 1) {
    return s >> "[" >> a.x >> "]";
  } else if constexpr (N == 2) {
    return s >> "[" >> a.x >> "," >> a.y >> "]";
  } else if constexpr (N == 3) {
    return s >> "[" >> a.x >> "," >> a.y >> "," >> a.z >> "]";
  } else if constexpr (N == 4) {
    return s >> "[" >> a.x >> "," >> a.y >> "," >> a.z >> "," >> a.w >> "]";
  } else {
  }
}

template <typename C>
inline std::basic_ostream<C>& operator<<(
    std::basic_ostream<C>& out, const mat2x2f& a) {
  return out << "[" << a.x << "," << a.y << "]";
}

template <typename C>
inline std::basic_istream<C>& operator>>(
    std::basic_istream<C>& out, mat2x2f& a) {
  return out >> "[" >> a.x >> "," >> a.y >> "]";
}

template <typename C>
inline std::basic_ostream<C>& operator<<(
    std::basic_ostream<C>& out, const mat3x3f& a) {
  return out << "[" << a.x << "," << a.y << "," << a.z << "]";
}

template <typename C>
inline std::basic_istream<C>& operator>>(
    std::basic_istream<C>& out, mat3x3f& a) {
  return out >> "[" >> a.x >> "," >> a.y >> "," >> a.z >> "]";
}

template <typename C>
inline std::basic_ostream<C>& operator<<(
    std::basic_ostream<C>& out, const mat4x4f& a) {
  return out << "[" << a.x << "," << a.y << "," << a.z << "," << a.w << "]";
}

template <typename C>
inline std::basic_istream<C>& operator>>(
    std::basic_istream<C>& out, mat4x4f& a) {
  return out >> "[" >> a.x >> "," >> a.y >> "," >> a.z >> "," >> a.w >> "]";
}

template <typename C>
inline std::basic_ostream<C>& operator<<(
    std::basic_ostream<C>& out, const frame2f& a) {
  return out << "[" << a.x << "," << a.y << "," << a.o << "]";
}
template <typename C>
inline std::basic_istream<C>& operator>>(
    std::basic_istream<C>& out, frame2f& a) {
  return out >> "[" >> a.x >> "," >> a.y >> "," >> a.o >> "]";
}

template <typename C>
inline std::basic_ostream<C>& operator<<(
    std::basic_ostream<C>& out, const frame3f& a) {
  return out << "[" << a.x << "," << a.y << "," << a.z << "," << a.o << "]";
}
template <typename C>
inline std::basic_istream<C>& operator>>(
    std::basic_istream<C>& out, frame3f& a) {
  return out >> "[" >> a.x >> "," >> a.y >> "," >> a.z >> "," >> a.o >> "]";
}

template <typename C>
inline std::basic_ostream<C>& operator<<(
    std::basic_ostream<C>& out, const quat4f& a) {
  return out << "[" << a.x << "," << a.y << "," << a.z << "," << a.w << "]";
}
template <typename C>
inline std::basic_istream<C>& operator>>(
    std::basic_istream<C>& out, quat4f& a) {
  return out >> "[" >> a.x >> "," >> a.y >> "," >> a.z >> "," >> a.w >> "]";
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// CUDA SUPPORT
// -----------------------------------------------------------------------------
#ifdef kernel
#undef kernel
#endif

#endif
