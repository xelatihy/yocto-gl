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
// MATH CONSTANTS AND FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

using byte   = unsigned char;
using uint   = unsigned int;
using ushort = unsigned short;

// Common type
template <typename... Ts>
using common_t = std::common_type_t<Ts...>;

template <typename T = double>
constexpr auto pi_t = (T)3.14159265358979323846;
constexpr auto pi   = 3.14159265358979323846;
constexpr auto pif  = (float)pi;

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
  return std::fmod(a, b);
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
  union {  // clang-format off
    array<T, 1> d;
    T x;
  };  // clang-format on

  constexpr kernel vec() : x{0} {}
  constexpr kernel vec(T x_) : x{x_} {}

  template <typename U>
  constexpr kernel explicit(!std::is_convertible_v<U, T>)
      vec(const vec<U, 1>& v)
      : x{(T)v.x} {}
  template <typename U>
  constexpr kernel explicit(!std::is_convertible_v<T, U>) operator vec<U, 1>() {
    return {(U)x};
  }

  constexpr kernel vec(const array<T, 1>& v) : x{v[0]} {}
  constexpr kernel operator array<T, 1>() { return {x}; }

  template <typename U>
  constexpr kernel explicit(!std::is_convertible_v<U, T>)
      vec(const array<U, 1>& v)
      : x{(T)v[0]} {}
  template <typename U>
  constexpr kernel explicit(!std::is_convertible_v<T, U>)
  operator array<U, 1>() {
    return {(U)x};
  }

  constexpr kernel T&       operator[](size_t i) { return d[i]; }
  constexpr kernel const T& operator[](size_t i) const { return d[i]; }
};

template <typename T>
struct vec<T, 2> {
  union {  // clang-format off
    array<T, 2> d;
    struct { T x, y; };
  };  // clang-format on

  constexpr kernel vec() : x{0}, y{0} {}
  constexpr kernel explicit vec(T v_) : x{v_}, y{v_} {}
  constexpr kernel vec(T x_, T y_) : x{x_}, y{y_} {}

  template <typename U>
  constexpr kernel explicit(!std::is_convertible_v<U, T>)
      vec(const vec<U, 2>& v)
      : x{(T)v.x}, y{(T)v.y} {}
  template <typename U>
  constexpr kernel explicit(!std::is_convertible_v<T, U>) operator vec<U, 2>() {
    return {(U)x, (U)y};
  }

  constexpr kernel vec(const array<T, 2>& v) : x{v[0]}, y{v[1]} {}
  constexpr kernel operator array<T, 2>() { return {x, y}; }

  template <typename U>
  constexpr kernel explicit(!std::is_convertible_v<U, T>)
      vec(const array<U, 2>& v)
      : x{(T)v[0]}, y{(T)v[1]} {}
  template <typename U>
  constexpr kernel explicit(!std::is_convertible_v<T, U>)
  operator array<U, 2>() {
    return {(U)x, (U)y};
  }

  constexpr kernel T&       operator[](size_t i) { return d[i]; }
  constexpr kernel const T& operator[](size_t i) const { return d[i]; }
};

template <typename T>
struct vec<T, 3> {
  union {  // clang-format off
    array<T, 3> d;
    struct { T x, y, z; };
  };  // clang-format on

  constexpr kernel vec() : x{0}, y{0}, z{0} {}
  constexpr kernel explicit vec(T v_) : x{v_}, y{v_}, z{v_} {}
  constexpr kernel vec(T x_, T y_, T z_) : x{x_}, y{y_}, z{z_} {}
  constexpr kernel vec(vec<T, 2> xy_, T z_) : x{xy_.x}, y{xy_.y}, z{z_} {}

  template <typename U>
  constexpr kernel explicit(!std::is_convertible_v<U, T>)
      vec(const vec<U, 3>& v)
      : x{(T)v.x}, y{(T)v.y}, z{(T)v.z} {}
  template <typename U>
  constexpr kernel explicit(!std::is_convertible_v<T, U>) operator vec<U, 3>() {
    return {(U)x, (U)y, (U)z};
  }
  constexpr kernel vec(const array<T, 3>& v) : x{v[0]}, y{v[1]}, z{v[2]} {}
  constexpr kernel operator array<T, 3>() { return {x, y, z}; }

  template <typename U>
  constexpr kernel explicit(!std::is_convertible_v<U, T>)
      vec(const array<U, 3>& v)
      : x{(T)v[0]}, y{(T)v[1]}, z{(T)v[2]} {}
  template <typename U>
  constexpr kernel explicit(!std::is_convertible_v<T, U>)
  operator array<U, 3>() {
    return {(U)x, (U)y, (U)z};
  }

  constexpr kernel T&       operator[](size_t i) { return d[i]; }
  constexpr kernel const T& operator[](size_t i) const { return d[i]; }
};

template <typename T>
struct vec<T, 4> {
  union {  // clang-format off
    array<T, 4> d;
    struct { T x, y, z, w; };
  };  // clang-format on

  constexpr kernel vec() : x{0}, y{0}, z{0}, w{0} {}
  constexpr kernel explicit vec(T v_) : x{v_}, y{v_}, z{v_}, w{v_} {}
  constexpr kernel vec(T x_, T y_, T z_, T w_) : x{x_}, y{y_}, z{z_}, w{w_} {}
  constexpr kernel vec(vec<T, 3> xyz_, T w_)
      : x{xyz_.x}, y{xyz_.y}, z{xyz_.z}, w{w_} {}

  template <typename U>
  constexpr kernel explicit(!std::is_convertible_v<U, T>)
      vec(const vec<U, 4>& v)
      : x{(T)v.x}, y{(T)v.y}, z{(T)v.z}, w{(T)v.w} {}
  template <typename U>
  constexpr kernel explicit(!std::is_convertible_v<T, U>) operator vec<U, 4>() {
    return {(U)x, (U)y, (U)z, (U)w};
  }
  constexpr kernel vec(const array<T, 4>& v)
      : x{v[0]}, y{v[1]}, z{v[2]}, w{v[3]} {}
  constexpr kernel operator array<T, 4>() { return {x, y, z, w}; }

  template <typename U>
  constexpr kernel explicit(!std::is_convertible_v<U, T>)
      vec(const array<U, 4>& v)
      : x{(T)v[0]}, y{(T)v[1]}, z{(T)v[2]}, w{(T)v[3]} {}
  template <typename U>
  constexpr kernel explicit(!std::is_convertible_v<T, U>)
  operator array<U, 4>() {
    return {(U)x, (U)y, (U)z, (U)w};
  }

  constexpr kernel T&       operator[](size_t i) { return d[i]; }
  constexpr kernel const T& operator[](size_t i) const { return d[i]; }
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
using vec1b = vec<byte, 3>;
using vec2b = vec<byte, 4>;
using vec3b = vec<byte, 3>;
using vec4b = vec<byte, 4>;
using vec1s = vec<size_t, 1>;
using vec2s = vec<size_t, 2>;
using vec3s = vec<size_t, 3>;
using vec4s = vec<size_t, 4>;

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
constexpr kernel vec<T, 3> xyz(const vec<T, 4>& a) {
  return {a.x, a.y, a.z};
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
constexpr kernel const T* begin(const vec<T, N>& a) {
  return &a.x;
}
template <typename T, size_t N>
constexpr kernel const T* end(const vec<T, N>& a) {
  return &a.x + N;
}
template <typename T, size_t N>
constexpr kernel T* begin(vec<T, N>& a) {
  return &a.x;
}
template <typename T, size_t N>
constexpr kernel T* end(vec<T, N>& a) {
  return &a.x + N;
}
template <typename T, size_t N>
constexpr kernel const T* data(const vec<T, N>& a) {
  return &a.x;
}
template <typename T, size_t N>
constexpr kernel T* data(vec<T, N>& a) {
  return &a.x;
}
template <size_t I, typename T, size_t N>
constexpr kernel T& get(vec<T, N>& a) noexcept {
  return a.d[I];
}
template <size_t I, typename T, size_t N>
constexpr kernel T&& get(vec<T, N>&& a) noexcept {
  return (T &&)(a.d[I]);
}
template <size_t I, typename T, size_t N>
constexpr kernel const T& get(const vec<T, N>& a) noexcept {
  return a.d[I];
}
template <size_t I, typename T, size_t N>
constexpr kernel const T&& get(const vec<T, N>&& a) noexcept {
  return (const T&&)(a.d[I]);
}

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
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel bool operator!=(const vec<T1, N>& a, const vec<T2, N>& b) {
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

template <typename T1, typename T2, typename T3,
    typename T = common_t<T1, T2, T3>>
constexpr kernel vec<T, 4> slerp(
    const vec<T1, 4>& a, const vec<T2, 4>& b, T3 u) {
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
        clamp(x.x, min.x, max), clamp(x.y, min.y, max), clamp(x.z, min, max.z)};
  } else if constexpr (N == 4) {
    return {clamp(x.x, min.x, max), clamp(x.y, min.y, max),
        clamp(x.z, min, max.z), clamp(x.w, min, max.w)};
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
    return {
        clamp(x.x, min, max.x), clamp(x.y, min, max.y), clamp(x.z, min, max.z)};
  } else if constexpr (N == 4) {
    return {clamp(x.x, min, max.x), clamp(x.y, min, max.y),
        clamp(x.z, min, max.z), clamp(x.w, min, max.w)};
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
      return 0;
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
      return 0;
    else
      return argmin(xyz(a));
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
constexpr kernel T mean(const vec<T, N>& a) {
  return sum(a) / N;
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

}  // namespace yocto

// -----------------------------------------------------------------------------
// CUSTOMIZATION POINTS FOR STD
// -----------------------------------------------------------------------------
namespace std {

// Structured binding support
template <typename T, size_t N>
struct tuple_size<yocto::vec<T, N>> {
  static const size_t value = N;
};
template <size_t I, typename T, size_t N>
struct tuple_element<I, yocto::vec<T, N>> {
  using type = T;
};

}  // namespace std

// -----------------------------------------------------------------------------
// MATRICES
// -----------------------------------------------------------------------------
namespace yocto {

// Small Fixed-size matrices stored in column major format.
template <typename T, size_t N, size_t M>
struct mat;

// Small Fixed-size matrices stored in column major format.
template <typename T, size_t N>
struct mat<T, N, 1> {
  union {  // clang-format off
    array<T, N * 1> d;
    array<vec<T, N>, 1> cols;
    vec<T, N> x;
  };  // clang-format on

  constexpr kernel mat() : x{1} {}
  constexpr kernel mat(const vec<T, N>& x_) : x{x_} {}

  constexpr kernel vec<T, N>& operator[](size_t i) { return cols[i]; }
  constexpr kernel const vec<T, N>& operator[](size_t i) const {
    return cols[i];
  }
};

// Small Fixed-size matrices stored in column major format.
template <typename T, size_t N>
struct mat<T, N, 2> {
  union {  // clang-format off
    array<T, N * 2> d;
    array<vec<T, N>, 2> cols;
    struct { vec<T, N> x, y; };
  };  // clang-format on

  constexpr kernel mat() : x{1, 0}, y{0, 1} {}
  constexpr kernel mat(const vec<T, N>& x_, const vec<T, N>& y_)
      : x{x_}, y{y_} {}

  constexpr kernel vec<T, N>& operator[](size_t i) { return cols[i]; }
  constexpr kernel const vec<T, N>& operator[](size_t i) const {
    return cols[i];
  }
};

// Small Fixed-size matrices stored in column major format.
template <typename T, size_t N>
struct mat<T, N, 3> {
  union {  // clang-format off
    array<T, N * 3> d;
    array<vec<T, N>, 3> cols;
    struct { vec<T, N> x, y, z; };
  };  // clang-format on

  constexpr kernel mat() : x{1, 0, 0}, y{0, 1, 0}, z{0, 0, 1} {}
  constexpr kernel mat(
      const vec<T, N>& x_, const vec<T, N>& y_, const vec<T, N>& z_)
      : x{x_}, y{y_}, z{z_} {}

  constexpr kernel vec<T, N>& operator[](size_t i) { return cols[i]; }
  constexpr kernel const vec<T, N>& operator[](size_t i) const {
    return cols[i];
  }
};

// Small Fixed-size matrices stored in column major format.
template <typename T, size_t N>
struct mat<T, N, 4> {
  union {  // clang-format off
    array<T, N * 4> d;
    array<vec<T, N>, 4> cols;
    struct { vec<T, N> x, y, z, w; };
  };  // clang-format on

  constexpr kernel mat()
      : x{1, 0, 0, 0}, y{0, 1, 0, 0}, z{0, 0, 1, 0}, w{0, 0, 0, 1} {}
  constexpr kernel mat(const vec<T, N>& x_, const vec<T, N>& y_,
      const vec<T, N>& z_, const vec<T, N>& w_)
      : x{x_}, y{y_}, z{z_}, w{w_} {}

  constexpr kernel vec<T, N>& operator[](size_t i) { return cols[i]; }
  constexpr kernel const vec<T, N>& operator[](size_t i) const {
    return cols[i];
  }
};

// Matrix aliases
using mat1f   = mat<float, 1, 1>;
using mat2f   = mat<float, 2, 2>;
using mat3f   = mat<float, 3, 3>;
using mat4f   = mat<float, 4, 4>;
using mat1x1f = mat<float, 1, 1>;
using mat2x2f = mat<float, 2, 2>;
using mat3x3f = mat<float, 3, 3>;
using mat4x4f = mat<float, 4, 4>;

// Identity matrices constants.
constexpr auto identity2x2f = mat2x2f{{1, 0}, {0, 1}};
constexpr auto identity3x3f = mat3x3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
constexpr auto identity4x4f = mat4x4f{
    {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};

// Matrix comparisons.
template <typename T1, typename T2, size_t N, size_t M>
constexpr kernel bool operator==(
    const mat<T1, N, M>& a, const mat<T2, N, M>& b) {
  if constexpr (M == 1) {
    return a.x == b.x;
  } else if constexpr (M == 2) {
    return a.x == b.x && a.y == b.y;
  } else if constexpr (M == 3) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
  } else if constexpr (M == 4) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.w == b.w;
  }
}
template <typename T1, typename T2, size_t N, size_t M>
constexpr kernel bool operator!=(
    const mat<T1, N, M>& a, const mat<T2, N, M>& b) {
  return !(a == b);
}

// Matrix operations.
template <typename T1, typename T2, size_t N, size_t M,
    typename T = common_t<T1, T2>>
constexpr kernel mat<T, N, M> operator+(
    const mat<T1, N, M>& a, const mat<T2, N, M>& b) {
  if constexpr (M == 1) {
    return {a.x + b.x};
  } else if constexpr (M == 2) {
    return {a.x + b.x, a.y + b.y};
  } else if constexpr (M == 3) {
    return {a.x + b.x, a.y + b.y, a.z + b.z};
  } else if constexpr (M == 4) {
    return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
  }
}
template <typename T1, typename T2, size_t N, size_t M,
    typename T = common_t<T1, T2>>
constexpr kernel mat<T, N, M> operator*(const mat<T1, N, M>& a, T2 b) {
  if constexpr (M == 1) {
    return {a.x * b};
  } else if constexpr (M == 2) {
    return {a.x * b, a.y * b};
  } else if constexpr (M == 3) {
    return {a.x * b, a.y * b, a.z * b};
  } else if constexpr (M == 4) {
    return {a.x * b, a.y * b, a.z * b, a.w * b};
  }
}
template <typename T1, typename T2, size_t N, size_t M,
    typename T = common_t<T1, T2>>
constexpr kernel vec<T, N> operator*(
    const mat<T1, N, M>& a, const vec<T2, M>& b) {
  if constexpr (M == 1) {
    return a.x * b.x;
  } else if constexpr (M == 2) {
    return a.x * b.x + a.y * b.y;
  } else if constexpr (M == 3) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
  } else if constexpr (M == 4) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
  }
}
template <typename T1, typename T2, size_t N, size_t M,
    typename T = common_t<T1, T2>>
constexpr kernel vec<T, M> operator*(
    const vec<T1, N>& a, const mat<T2, N, M>& b) {
  if constexpr (M == 1) {
    return {dot(a, b.x)};
  } else if constexpr (M == 2) {
    return {dot(a, b.x), dot(a, b.y)};
  } else if constexpr (M == 3) {
    return {dot(a, b.x), dot(a, b.y), dot(a, b.z)};
  } else if constexpr (M == 4) {
    return {dot(a, b.x), dot(a, b.y), dot(a, b.z), dot(a, b.w)};
  }
}
template <typename T1, typename T2, size_t N, size_t M, size_t K,
    typename T = common_t<T1, T2>>
constexpr kernel mat<T, N, M> operator*(
    const mat<T1, N, K>& a, const mat<T2, K, M>& b) {
  if constexpr (M == 1) {
    return {a * b.x};
  } else if constexpr (M == 2) {
    return {a * b.x, a * b.y};
  } else if constexpr (M == 3) {
    return {a * b.x, a * b.y, a * b.z};
  } else if constexpr (M == 4) {
    return {a * b.x, a * b.y, a * b.z, a * b.w};
  }
}

// Matrix assignments.
template <typename T, typename T1, size_t N, size_t M>
constexpr kernel mat<T, N, M>& operator+=(
    mat<T, N, M>& a, const mat<T1, N, M>& b) {
  return a = a + b;
}
template <typename T, typename T1, size_t N>
constexpr kernel mat<T, N, N>& operator*=(
    mat<T, N, N>& a, const mat<T1, N, N>& b) {
  return a = a * b;
}
template <typename T, typename T1, size_t N, size_t M>
constexpr kernel mat<T, N, M>& operator*=(mat<T, N, M>& a, T1 b) {
  return a = a * b;
}

// Matrix diagonals and transposes.
template <typename T, size_t N>
constexpr kernel vec<T, N> diagonal(const mat<T, N, N>& a) {
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
template <typename T, size_t N>
constexpr kernel mat<T, N, N> transpose(const mat<T, N, N>& a) {
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
template <typename T, size_t N>
constexpr kernel T determinant(const mat<T, N, N>& a) {
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
template <typename T, size_t N>
constexpr kernel mat<T, N, N> adjoint(const mat<T, N, N>& a) {
  if constexpr (N == 1) {
    return {{a.x.x}};
  } else if constexpr (N == 2) {
    return {{a.y.y, -a.x.y}, {-a.y.x, a.x.x}};
  } else if constexpr (N == 3) {
    return transpose(
        mat<T, 3, 3>{cross(a.y, a.z), cross(a.z, a.x), cross(a.x, a.y)});
  } else if constexpr (N == 4) {
    return {};  // TODO
  }
}
template <typename T, size_t N>
constexpr kernel mat<T, N, N> inverse(const mat<T, N, N>& a) {
  return adjoint(a) * (1 / determinant(a));
}

// Constructs a basis from a direction
template <typename T>
constexpr kernel mat<T, 3, 3> basis_fromz(const vec<T, 3>& v) {
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
template <typename T, size_t N>
struct frame;

// Rigid frames stored as a column-major affine transform matrix.
template <typename T>
struct frame<T, 2> {
  union {  // clang-format off
    array<vec<T, 2>, 3> cols;
    struct { vec<T, 2> x, y, o; };
    struct { mat<T, 2, 2> m; vec<T, 2> t; };
  };  // clang-format on

  constexpr kernel frame() : x{1, 0}, y{0, 1}, o{0, 0} {}
  constexpr kernel frame(
      const vec<T, 2>& x_, const vec<T, 2>& y_, const vec<T, 2>& o_)
      : x{x_}, y{y_}, o{o_} {}
  constexpr kernel frame(const mat<T, 2, 2>& xy_, const vec<T, 2>& o_)
      : x{xy_.x}, y{xy_.y}, o{o_} {}

  explicit constexpr kernel frame(const mat<T, 3, 3>& m)
      : x{m.x.x, m.x.y}, y{m.y.x, m.y.y}, o{m.z.x, m.z.y} {}
  explicit constexpr kernel operator mat<T, 3, 3>() const {
    return {{x.x, x.y, 0}, {y.x, y.y, 0}, {o.x, o.y, 1}};
  }

  constexpr kernel vec<T, 2>& operator[](size_t i) { return cols[i]; }
  constexpr kernel const vec<T, 2>& operator[](size_t i) const {
    return cols[i];
  }
};

// Rigid frames stored as a column-major affine transform matrix.
template <typename T>
struct frame<T, 3> {
  union {  // clang-format off
    array<vec<T, 3>, 4> cols;
    struct { vec<T, 3> x, y, z, o; };
    struct { mat<T, 3, 3> m; vec<T, 3> t; };
  };      // clang-format off

  constexpr  kernel frame() : x{1, 0, 0}, y{0, 1, 0}, z{0, 0, 1}, o{0, 0, 0} {}
  constexpr  kernel frame(const vec<T, 3>& x_, const vec<T, 3>& y_, const vec<T, 3>& z_,
      const vec<T, 3>& o_)
      : x{x_}, y{y_}, z{z_}, o{o_} {}
  constexpr  kernel frame(const mat<T, 3, 3>& xyz_, const vec<T, 3>& o_)
      : x{xyz_.x}, y{xyz_.y}, z{xyz_.z}, o{o_} {}

  explicit constexpr  kernel frame(const mat<T, 4, 4>& m)
      : x{m.x.x, m.x.y, m.x.z}
      , y{m.y.x, m.y.y, m.y.z}
      , z{m.z.x, m.z.y, m.z.z}
      , o{m.w.x, m.w.y, m.w.z} {}
  explicit constexpr  kernel operator mat<T, 4, 4>() const {
    return {{x.x, x.y, x.z, 0}, {y.x, y.y, y.z, 0}, {z.x, z.y, z.z, 0},
        {o.x, o.y, o.z, 1}};
  }

  constexpr  kernel vec<T, 3>&       operator[](size_t i) { return cols[i]; }
  constexpr  kernel const vec<T, 3>& operator[](size_t i) const { return cols[i]; }
};

// Frame aliases
using frame2f = frame<float, 2>;
using frame3f = frame<float, 3>;

// Identity frames.
constexpr auto identity2x3f = frame2f{{1, 0}, {0, 1}, {0, 0}};
constexpr auto identity3x4f = frame3f{
    {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};

// Frame properties
template <typename T, size_t N>
constexpr kernel mat<T, N, N> rotation(const frame<T, N>& a) {
  if constexpr (N == 2) {
    return {a.x, a.y};
  } else if constexpr (N == 3) {
    return {a.x, a.y, a.z};
  }
}
template <typename T, size_t N>
constexpr kernel vec<T, N> translation(const frame<T, N>& a) {
  if constexpr (N == 2) {
    return a.o;
  } else if constexpr (N == 3) {
    return a.o;
  }
}

// Frame/mat conversion
template <typename T, size_t N>
constexpr kernel frame<T, N - 1> mat_to_frame(const mat<T, N, N>& m) {
  return (frame<T, N - 1>)m;
}
template <typename T, size_t N>
constexpr kernel mat<T, N + 1, N + 1> frame_to_mat(const frame<T, N>& f) {
  return (mat<T, N + 1, N + 1>)f;
}

// Frame comparisons.
template <typename T1, typename T2, size_t N>
constexpr kernel bool operator==(const frame<T1, N>& a, const frame<T2, N>& b) {
  if constexpr (N == 2) {
    return a.x == b.x && a.y == b.y && a.o == b.o;
  } else if constexpr (N == 3) {
    return a.x == b.x && a.y == b.y && a.z == b.z && a.o == b.o;
  }
}
template <typename T1, typename T2, size_t N>
constexpr kernel bool operator!=(const frame<T1, N>& a, const frame<T2, N>& b) {
  return !(a == b);
}

// Frame composition, equivalent to affine matrix product.
template <typename T, size_t N>
constexpr kernel frame<T, N> operator*(const frame<T, N>& a, const frame<T, N>& b) {
  return {a.m * b.m, a.m * b.t + a.t};
}
template <typename T, size_t N>
constexpr kernel frame<T, N>& operator*=(frame<T, N>& a, const frame<T, N>& b) {
  return a = a * b;
}

// Frame inverse, equivalent to rigid affine inverse.
template <typename T, size_t N>
constexpr kernel frame<T, N> inverse(const frame<T, N>& a, bool non_rigid = false) {
  if (non_rigid) {
    auto minv = inverse(a.m);
    return {minv, -(minv * a.t)};
  } else {
    auto minv = transpose(a.m);
    return {minv, -(minv * a.t)};
  }
}

// Frame construction from axis.
template <typename T>
constexpr kernel frame<T, 3> frame_fromz(const vec<T, 3>& o, const vec<T, 3>& v) {
  // https://graphics.pixar.com/library/OrthonormalB/paper.pdf
  if constexpr (std::is_same_v<T, float>) {
    auto z    = normalize(v);
    auto sign = copysignf((T)1, z.z);
    auto a    = -1 / (sign + z.z);
    auto b    = z.x * z.y * a;
    auto x    = vec<T, 3>{1 + sign * z.x * z.x * a, sign * b, -sign * z.x};
    auto y    = vec<T, 3>{b, sign + z.y * z.y * a, -z.y};
    return {x, y, z, o};
  } else if constexpr (std::is_same_v<T, double>) {
    // TODO: double
    return {};
  }
}
template <typename T>
constexpr kernel frame<T, 3> frame_fromzx(
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

// Quaternions to represent rotations
template <typename T, size_t N>
struct quat;

// Quaternions to represent rotations
template <typename T>
struct quat<T, 4> {
  union {  // clang-format off
    array<T, 4> d;
    struct { T x, y, z, w; };
  };  // clang-format on

  constexpr kernel quat() : x{0}, y{0}, z{0}, w{1} {}
  constexpr kernel quat(T x_, T y_, T z_, T w_) : x{x_}, y{y_}, z{z_}, w{w_} {}

  constexpr kernel T&       operator[](size_t i) { return d[i]; }
  constexpr kernel const T& operator[](size_t i) const { return d[i]; }
};

// Quaternion aliases
using quat4f = quat<float, 4>;

[[deprecated]] constexpr auto identity_quat4f = quat4f{0, 0, 0, 1};

// Quaternion operations
template <typename T>
constexpr kernel quat<T, 4> operator+(
    const quat<T, 4>& a, const quat<T, 4>& b) {
  return {a.x + b.x, a.y + b.y, a.z + b.z, a.w + b.w};
}
template <typename T>
constexpr kernel quat<T, 4> operator*(const quat<T, 4>& a, T b) {
  return {a.x * b, a.y * b, a.z * b, a.w * b};
}
template <typename T>
constexpr kernel quat<T, 4> operator/(const quat<T, 4>& a, T b) {
  return {a.x / b, a.y / b, a.z / b, a.w / b};
}
template <typename T>
constexpr kernel quat<T, 4> operator*(
    const quat<T, 4>& a, const quat<T, 4>& b) {
  return {a.x * b.w + a.w * b.x + a.y * b.w - a.z * b.y,
      a.y * b.w + a.w * b.y + a.z * b.x - a.x * b.z,
      a.z * b.w + a.w * b.z + a.x * b.y - a.y * b.x,
      a.w * b.w - a.x * b.x - a.y * b.y - a.z * b.z};
}

// Quaternion operations
template <typename T>
constexpr kernel T dot(const quat<T, 4>& a, const quat<T, 4>& b) {
  return a.x * b.x + a.y * b.y + a.z * b.z + a.w * b.w;
}
template <typename T>
constexpr kernel T length(const quat<T, 4>& a) {
  return sqrt(dot(a, a));
}
template <typename T>
constexpr kernel quat<T, 4> normalize(const quat<T, 4>& a) {
  auto l = length(a);
  return (l != 0) ? a / l : a;
}
template <typename T>
constexpr kernel quat<T, 4> conjugate(const quat<T, 4>& a) {
  return {-a.x, -a.y, -a.z, a.w};
}
template <typename T>
constexpr kernel quat<T, 4> inverse(const quat<T, 4>& a) {
  return conjugate(a) / dot(a, a);
}
template <typename T>
constexpr kernel T uangle(const quat<T, 4>& a, const quat<T, 4>& b) {
  auto d = dot(a, b);
  return d > 1 ? 0 : acos(d < -1 ? -1 : d);
}
template <typename T>
constexpr kernel quat<T, 4> lerp(
    const quat<T, 4>& a, const quat<T, 4>& b, T t) {
  return a * (1 - t) + b * t;
}
template <typename T>
constexpr kernel quat<T, 4> nlerp(
    const quat<T, 4>& a, const quat<T, 4>& b, T t) {
  return normalize(lerp(a, b, t));
}
template <typename T>
constexpr kernel quat<T, 4> slerp(
    const quat<T, 4>& a, const quat<T, 4>& b, T t) {
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
template <typename T, size_t N>
constexpr kernel vec<T, N> transform_point(
    const mat<T, N + 1, N + 1>& a, const vec<T, N>& b) {
  if constexpr (N == 2) {
    auto tvb = a * vec<T, 3>{b.x, b.y, 1};
    return vec<T, 2>{tvb.x, tvb.y} / tvb.z;
  } else if constexpr (N == 3) {
    auto tvb = a * vec<T, 4>{b.x, b.y, b.z, 1};
    return vec<T, 3>{tvb.x, tvb.y, tvb.z} / tvb.w;
  }
}
template <typename T, size_t N>
constexpr kernel vec<T, N> transform_vector(
    const mat<T, N + 1, N + 1>& a, const vec<T, N>& b) {
  if constexpr (N == 2) {
    auto tvb = a * vec<T, 3>{b.x, b.y, 0};
    return vec<T, 2>{tvb.x, tvb.y} / tvb.z;
  } else if constexpr (N == 3) {
    auto tvb = a * vec<T, 4>{b.x, b.y, b.z, 0};
    return vec<T, 3>{tvb.x, tvb.y, tvb.z} / tvb.w;
  }
}
template <typename T, size_t N>
constexpr kernel vec<T, N> transform_direction(
    const mat<T, N + 1, N + 1>& a, const vec<T, N>& b) {
  return normalize(transform_vector(a, b));
}
template <typename T, size_t N>
constexpr kernel vec<T, N> transform_normal(
    const mat<T, N + 1, N + 1>& a, const vec<T, N>& b) {
  return normalize(transform_vector(transpose(inverse(a)), b));
}
template <typename T, size_t N>
constexpr kernel vec<T, N> transform_vector(
    const mat<T, N, N>& a, const vec<T, N>& b) {
  return a * b;
}
template <typename T, size_t N>
constexpr kernel vec<T, N> transform_direction(
    const mat<T, N, N>& a, const vec<T, N>& b) {
  return normalize(transform_vector(a, b));
}
template <typename T, size_t N>
constexpr kernel vec<T, N> transform_normal(
    const mat<T, N, N>& a, const vec<T, N>& b) {
  return normalize(transform_vector(transpose(inverse(a)), b));
}

// Transforms points, vectors and directions by frames.
template <typename T, size_t N>
constexpr kernel vec<T, N> transform_point(
    const frame<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 2) {
    return a.x * b.x + a.y * b.y + a.o;
  } else if constexpr (N == 3) {
    return a.x * b.x + a.y * b.y + a.z * b.z + a.o;
  }
}
template <typename T, size_t N>
constexpr kernel vec<T, N> transform_vector(
    const frame<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 2) {
    return a.x * b.x + a.y * b.y;
  } else if constexpr (N == 3) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
  }
}
template <typename T, size_t N>
constexpr kernel vec<T, N> transform_direction(
    const frame<T, N>& a, const vec<T, N>& b) {
  return normalize(transform_vector(a, b));
}
template <typename T, size_t N>
constexpr kernel vec<T, N> transform_normal(
    const frame<T, N>& a, const vec<T, N>& b, bool non_rigid = false) {
  if (non_rigid) {
    return transform_normal(rotation(a), b);
  } else {
    return normalize(transform_vector(a, b));
  }
}

// Transforms points, vectors and directions by frames.
template <typename T, size_t N>
constexpr kernel vec<T, N> transform_point_inverse(
    const frame<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 2) {
    return {dot(a.x, (b - a.o)), dot(a.y, (b - a.o))};
  } else if constexpr (N == 3) {
    return {dot(a.x, (b - a.o)), dot(a.y, (b - a.o)), dot(a.z, (b - a.o))};
  }
}
template <typename T, size_t N>
constexpr kernel vec<T, N> transform_vector_inverse(
    const frame<T, N>& a, const vec<T, N>& b) {
  if constexpr (N == 2) {
    return {dot(a.x, b), dot(a.y, b)};
  } else if constexpr (N == 3) {
    return {dot(a.x, b), dot(a.y, b), dot(a.z, b)};
  }
}
template <typename T, size_t N>
constexpr kernel vec<T, N> transform_direction_inverse(
    const frame<T, N>& a, const vec<T, N>& b) {
  return normalize(transform_vector_inverse(a, b));
}

// Translation, scaling and rotations transforms.
template <typename T>
constexpr kernel frame<T, 3> translation_frame(const vec<T, 3>& a) {
  return {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, a};
}
template <typename T>
constexpr kernel frame<T, 3> scaling_frame(const vec<T, 3>& a) {
  return {{a.x, 0, 0}, {0, a.y, 0}, {0, 0, a.z}, {0, 0, 0}};
}
template <typename T>
constexpr kernel frame<T, 3> rotation_frame(const vec<T, 3>& axis, T angle) {
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
constexpr kernel frame<T, 3> rotation_frame(const vec<T, 4>& quat) {
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
constexpr kernel frame<T, 3> rotation_frame(const quat<T, 4>& quat) {
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
constexpr kernel frame<T, 3> rotation_frame(const mat<T, 3, 3>& rot) {
  return {rot.x, rot.y, rot.z, {0, 0, 0}};
}

// Lookat frame. Z-axis can be inverted with inv_xz.
template <typename T>
constexpr kernel frame<T, 3> lookat_frame(const vec<T, 3>& eye,
    const vec<T, 3>& center, const vec<T, 3>& up, bool inv_xz = false) {
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
constexpr kernel mat<T, 4, 4> frustum_mat(T l, T r, T b, T t, T n, T f) {
  return {{2 * n / (r - l), 0, 0, 0}, {0, 2 * n / (t - b), 0, 0},
      {(r + l) / (r - l), (t + b) / (t - b), -(f + n) / (f - n), -1},
      {0, 0, -2 * f * n / (f - n), 0}};
}
template <typename T>
constexpr kernel mat<T, 4, 4> ortho_mat(T l, T r, T b, T t, T n, T f) {
  return {{2 / (r - l), 0, 0, 0}, {0, 2 / (t - b), 0, 0},
      {0, 0, -2 / (f - n), 0},
      {-(r + l) / (r - l), -(t + b) / (t - b), -(f + n) / (f - n), 1}};
}
template <typename T>
constexpr kernel mat<T, 4, 4> ortho2d_mat(T left, T right, T bottom, T top) {
  return ortho_mat(left, right, bottom, top, -1, 1);
}
template <typename T>
constexpr kernel mat<T, 4, 4> ortho_mat(T xmag, T ymag, T near, T far) {
  return {{1 / xmag, 0, 0, 0}, {0, 1 / ymag, 0, 0}, {0, 0, 2 / (near - far), 0},
      {0, 0, (far + near) / (near - far), 1}};
}
template <typename T>
constexpr kernel mat<T, 4, 4> perspective_mat(T fovy, T aspect, T near, T far) {
  auto tg = tan(fovy / 2);
  return {{1 / (aspect * tg), 0, 0, 0}, {0, 1 / tg, 0, 0},
      {0, 0, (far + near) / (near - far), -1},
      {0, 0, 2 * far * near / (near - far), 0}};
}
template <typename T>
constexpr kernel mat<T, 4, 4> perspective_mat(T fovy, T aspect, T near) {
  auto tg = tan(fovy / 2);
  return {{1 / (aspect * tg), 0, 0, 0}, {0, 1 / tg, 0, 0}, {0, 0, -1, -1},
      {0, 0, 2 * near, 0}};
}

// Rotation conversions.
template <typename T>
constexpr kernel pair<vec<T, 3>, T> rotation_axisangle(const vec<T, 4>& quat) {
  return {normalize(vec<T, 3>{quat.x, quat.y, quat.z}), 2 * acos(quat.w)};
}
template <typename T>
constexpr kernel vec<T, 4> rotation_quat(const vec<T, 3>& axis, T angle) {
  auto len = length(axis);
  if (len == 0) return {0, 0, 0, 1};
  return vec<T, 4>{sin(angle / 2) * axis.x / len, sin(angle / 2) * axis.y / len,
      sin(angle / 2) * axis.z / len, cos(angle / 2)};
}
template <typename T>
constexpr kernel vec<T, 4> rotation_quat(const vec<T, 4>& axisangle) {
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
constexpr kernel vec<I, 2> image_coords(const vec<T, 2>& mouse_pos,
    const vec<T, 2>& center, T scale, const vec<I, 2>& size,
    bool clamped = true) {
  auto xy = (mouse_pos - center) / scale;
  auto ij = (vec<I, 2>)round(xy + size / (T)2);
  return clamped ? clamp(ij, vec<I, 2>{0, 0}, size) : ij;
}

// Center image and autofit. Returns center and scale.
template <typename T, typename I>
constexpr kernel pair<vec<T, 2>, T> camera_imview(const vec<T, 2>& center,
    T scale, const vec<I, 2>& imsize, const vec<I, 2>& winsize,
    bool zoom_to_fit) {
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
constexpr kernel pair<vec<T, 3>, vec<T, 3>> camera_turntable(
    const vec<T, 3>& from_, const vec<T, 3>& to_, const vec<T, 3>& up,
    const vec<T, 2>& rotate, T dolly, const vec<T, 2>& pan) {
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
constexpr kernel pair<frame<T, 3>, T> camera_turntable(
    const frame<T, 3>& frame_, T focus, const vec<T, 2>& rotate, T dolly,
    const vec<T, 2>& pan) {
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
    frame           = lookat_frame(new_o, new_center, vec<T, 3>{0, 1, 0});
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
constexpr kernel frame<T, 3> camera_fpscam(const frame<T, 3>& frame,
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
constexpr kernel vec2i get_image_coords(const vec<T, 2>& mouse_pos,
    const vec<T, 2>& center, T scale, const vec<I, 2>& txt_size) {
  auto xyf = (mouse_pos - center) / scale;
  return vec2i{(int)round(xyf.x + txt_size.x / (T)2),
      (int)round(xyf.y + txt_size.y / (T)2)};
}

// Center image and autofit.
template <typename T, typename I>
constexpr kernel void update_imview(vec<T, 2>& center, T& scale,
    const vec<I, 2>& imsize, const vec<I, 2>& winsize, bool zoom_to_fit) {
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
constexpr kernel void update_turntable(vec<T, 3>& from, vec<T, 3>& to,
    vec<T, 3>& up, const vec<T, 2>& rotate, T dolly, const vec<T, 2>& pan) {
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
constexpr kernel void update_turntable(frame<T, 3>& frame, T& focus,
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
constexpr kernel void update_fpscam(
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

// Python range. Construct an object that iterates over an integer sequence.
template <typename T>
constexpr kernel auto range(T max);
template <typename T>
constexpr kernel auto range(T min, T max);
template <typename T>
constexpr kernel auto range(T min, T max, T step);

// Python range in 2d. Construct an object that iterates over a 2d integer
// sequence.
template <typename T>
constexpr kernel auto range(vec<T, 2> max);
template <typename T>
constexpr kernel auto range(array<T, 2> max);

// Python enumerate
template <typename Sequence, typename T = size_t>
constexpr kernel auto enumerate(const Sequence& sequence, T start = 0);
template <typename Sequence, typename T = size_t>
constexpr kernel auto enumerate(Sequence& sequence, T start = 0);

// Python zip
template <typename Sequence1, typename Sequence2>
constexpr kernel auto zip(
    const Sequence1& sequence1, const Sequence2& sequence2);
template <typename Sequence1, typename Sequence2>
constexpr kernel auto zip(Sequence1& sequence1, Sequence2& sequence2);
template <typename Sequence1, typename Sequence2>
constexpr kernel auto zip(const Sequence1& sequence1, Sequence2& sequence2);
template <typename Sequence1, typename Sequence2>
constexpr kernel auto zip(Sequence1& sequence1, const Sequence2& sequence2);

// Implementation of Python range.
template <typename T>
constexpr kernel auto range(T max) {
  return range((T)0, max);
}
template <typename T>
constexpr kernel auto range(T min, T max) {
  struct range_iterator {
    T                     index;
    constexpr kernel void operator++() { ++index; }
    constexpr kernel bool operator!=(const range_iterator& other) const {
      return index != other.index;
    }
    constexpr kernel T operator*() const { return index; }
  };
  struct range_helper {
    T                               begin_ = 0, end_ = 0;
    constexpr kernel range_iterator begin() const { return {begin_}; }
    constexpr kernel range_iterator end() const { return {end_}; }
  };
  return range_helper{min, max};
}
template <typename T>
constexpr kernel auto range(T min, T max, T step) {
  struct range_iterator {
    T                     index;
    T                     step;
    constexpr kernel void operator++() { index += step; }
    constexpr kernel bool operator!=(const range_iterator& other) const {
      return index != other.index;
    }
    constexpr kernel T operator*() const { return index; }
  };
  struct range_helper {
    T                               begin_ = 0, end_ = 0, step_ = 0;
    constexpr kernel range_iterator begin() const { return {begin_, step_}; }
    constexpr kernel range_iterator end() const {
      return {begin_ + ((end_ - begin_) / step_) * step_, step_};
    }
  };
  return range_helper{min, max, step};
}

// Python range in 2d. Construct an object that iterates over a 2d integer
// sequence.
template <typename T>
constexpr kernel auto range(vec<T, 2> max) {
  struct range_sentinel {};
  struct range_iterator {
    vec<T, 2>             index, end;
    constexpr kernel void operator++() {
      ++index[0];
      if (index[0] >= end[0]) {
        index[0] = 0;
        index[1]++;
      }
    }
    constexpr kernel bool operator!=(const range_sentinel&) const {
      return index[1] != end[1];
    }
    constexpr kernel vec<T, 2> operator*() const { return index; }
  };
  struct range_sequence {
    vec<T, 2>                       end_ = {0, 0};
    constexpr kernel range_iterator begin() const { return {{0, 0}, end_}; }
    constexpr kernel range_sentinel end() const { return {}; }
  };
  return range_sequence{max};
}
template <typename T>
constexpr kernel auto range(array<T, 2> max) {
  return range(vec<T, 2>{max});
}

// Implementation of Python enumerate.
template <typename Sequence, typename T>
constexpr kernel auto enumerate(const Sequence& sequence, T start) {
  using Iterator  = typename Sequence::const_iterator;
  using Reference = typename Sequence::const_reference;
  struct enumerate_iterator {
    T                     index;
    Iterator              iterator;
    constexpr kernel bool operator!=(const enumerate_iterator& other) const {
      return index != other.index;
    }
    constexpr kernel void operator++() {
      ++index;
      ++iterator;
    }
    constexpr kernel pair<const T&, Reference> operator*() const {
      return {index, *iterator};
    }
  };
  struct enumerate_helper {
    const Sequence&       sequence;
    T                     begin_, end_;
    constexpr kernel auto begin() {
      return enumerate_iterator{begin_, std::begin(sequence)};
    }
    constexpr kernel auto end() {
      return enumerate_iterator{end_, std::end(sequence)};
    }
  };
  return enumerate_helper{sequence, 0, size(sequence)};
}

// Python enumerate
template <typename Sequence, typename T>
constexpr kernel auto enumerate(Sequence& sequence, T start) {
  using Iterator  = typename Sequence::iterator;
  using Reference = typename Sequence::reference;
  struct enumerate_iterator {
    T                     index;
    Iterator              iterator;
    constexpr kernel bool operator!=(const enumerate_iterator& other) const {
      return index != other.index;
    }
    constexpr kernel void operator++() {
      ++index;
      ++iterator;
    }
    constexpr kernel pair<T&, Reference> operator*() const {
      return {index, *iterator};
    }
  };
  struct enumerate_helper {
    Sequence&             sequence;
    T                     begin_, end_;
    constexpr kernel auto begin() {
      return enumerate_iterator{begin_, std::begin(sequence)};
    }
    constexpr kernel auto end() {
      return enumerate_iterator{end_, std::end(sequence)};
    }
  };
  return enumerate_helper{sequence, 0, size(sequence)};
}

// Python zip
template <typename Sequence1, typename Sequence2>
constexpr kernel auto zip(
    const Sequence1& sequence1, const Sequence2& sequence2) {
  using Iterator1  = typename Sequence1::const_iterator;
  using Reference1 = typename Sequence1::const_reference;
  using Iterator2  = typename Sequence2::const_iterator;
  using Reference2 = typename Sequence2::const_reference;
  struct zip_iterator {
    Iterator1             iterator1;
    Iterator2             iterator2;
    constexpr kernel bool operator!=(const zip_iterator& other) const {
      return iterator1 != other.iterator1;
    }
    constexpr kernel void operator++() {
      ++iterator1;
      ++iterator2;
    }
    constexpr kernel pair<Reference1, Reference2> operator*() const {
      return {*iterator1, *iterator2};
    }
  };
  struct zip_helper {
    const Sequence1&      sequence1;
    const Sequence2&      sequence2;
    constexpr kernel auto begin() {
      return zip_iterator{std::begin(sequence1), std::begin(sequence2)};
    }
    constexpr kernel auto end() {
      return zip_iterator{std::end(sequence1), std::end(sequence2)};
    }
  };
  return zip_helper{sequence1, sequence2};
}

// Implementation of Python zip
template <typename Sequence1, typename Sequence2>
constexpr kernel auto zip(Sequence1& sequence1, Sequence2& sequence2) {
  using Iterator1  = typename Sequence1::iterator;
  using Reference1 = typename Sequence1::reference;
  using Iterator2  = typename Sequence2::iterator;
  using Reference2 = typename Sequence2::reference;
  struct zip_iterator {
    Iterator1             iterator1;
    Iterator2             iterator2;
    constexpr kernel bool operator!=(const zip_iterator& other) const {
      return iterator1 != other.iterator1;
    }
    constexpr kernel void operator++() {
      ++iterator1;
      ++iterator2;
    }
    constexpr kernel pair<Reference1, Reference2> operator*() const {
      return {*iterator1, *iterator2};
    }
  };
  struct zip_helper {
    Sequence1&            sequence1;
    Sequence2&            sequence2;
    constexpr kernel auto begin() {
      return zip_iterator{std::begin(sequence1), std::begin(sequence2)};
    }
    constexpr kernel auto end() {
      return zip_iterator{std::end(sequence1), std::end(sequence2)};
    }
  };
  return zip_helper{sequence1, sequence2};
}

// Implementation of Python zip
template <typename Sequence1, typename Sequence2>
constexpr kernel auto zip(const Sequence1& sequence1, Sequence2& sequence2) {
  using Iterator1  = typename Sequence1::const_iterator;
  using Reference1 = typename Sequence1::const_reference;
  using Iterator2  = typename Sequence2::iterator;
  using Reference2 = typename Sequence2::reference;
  struct zip_iterator {
    Iterator1             iterator1;
    Iterator2             iterator2;
    constexpr kernel bool operator!=(const zip_iterator& other) const {
      return iterator1 != other.iterator1;
    }
    constexpr kernel void operator++() {
      ++iterator1;
      ++iterator2;
    }
    constexpr kernel pair<Reference1, Reference2> operator*() const {
      return {*iterator1, *iterator2};
    }
  };
  struct zip_helper {
    const Sequence1&      sequence1;
    Sequence2&            sequence2;
    constexpr kernel auto begin() {
      return zip_iterator{std::begin(sequence1), std::begin(sequence2)};
    }
    constexpr kernel auto end() {
      return zip_iterator{std::end(sequence1), std::end(sequence2)};
    }
  };
  return zip_helper{sequence1, sequence2};
}

// Implementation of Python zip
template <typename Sequence1, typename Sequence2>
constexpr kernel auto zip(Sequence1& sequence1, const Sequence2& sequence2) {
  using Iterator1  = typename Sequence1::iterator;
  using Reference1 = typename Sequence1::reference;
  using Iterator2  = typename Sequence2::const_iterator;
  using Reference2 = typename Sequence2::const_reference;
  struct zip_iterator {
    Iterator1             iterator1;
    Iterator2             iterator2;
    constexpr kernel bool operator!=(const zip_iterator& other) const {
      return iterator1 != other.iterator1;
    }
    constexpr kernel void operator++() {
      ++iterator1;
      ++iterator2;
    }
    constexpr kernel pair<Reference1, Reference2> operator*() const {
      return {*iterator1, *iterator2};
    }
  };
  struct zip_helper {
    Sequence1&            sequence1;
    const Sequence2&      sequence2;
    constexpr kernel auto begin() {
      return zip_iterator{std::begin(sequence1), std::begin(sequence2)};
    }
    constexpr kernel auto end() {
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
constexpr kernel std::ptrdiff_t ssize(const T& container) {
  return (std::ptrdiff_t)std::size(container);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// CUDA SUPPORT
// -----------------------------------------------------------------------------
#ifdef kernel
#undef kernel
#endif

#endif
