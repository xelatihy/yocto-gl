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
#include <tuple>
#include <type_traits>
#include <utility>

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
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
// MATH CONCEPTS AND TYPE TRAITS
// -----------------------------------------------------------------------------
namespace yocto {

// Index type
using index_t = std::size_t;

// Vec type traits
template <typename T>
struct is_vec_s : std::bool_constant<false> {};
template <typename T>
constexpr auto is_vec = is_vec_s<std::remove_cvref_t<T>>::value;
template <typename T>
struct vec_size_s : std::integral_constant<index_t, 0> {};
template <typename T>
constexpr auto vec_size = vec_size_s<std::remove_cvref_t<T>>::value;
template <typename T>
struct vec_etype_s {
  using type = void;
};
template <typename T>
using vec_etype = typename vec_etype_s<std::remove_cvref_t<T>>::type;

// Mat type traits
template <typename T>
struct is_mat_s : std::bool_constant<false> {};
template <typename T>
constexpr auto is_mat = is_mat_s<std::remove_cvref_t<T>>::value;
template <typename T>
struct mat_cols_s : std::integral_constant<index_t, 0> {};
template <typename T>
constexpr auto mat_cols = mat_cols_s<std::remove_cvref_t<T>>::value;
template <typename T>
struct mat_rows_s : std::integral_constant<index_t, 0> {};
template <typename T>
constexpr auto mat_rows = mat_rows_s<std::remove_cvref_t<T>>::value;
template <typename T>
struct mat_etype_s {
  using type = void;
};
template <typename T>
using mat_etype = typename mat_etype_s<std::remove_cvref_t<T>>::type;

// Frame type traits
template <typename T>
struct is_frame_s : std::bool_constant<false> {};
template <typename T>
constexpr auto is_frame = is_frame_s<std::remove_cvref_t<T>>::value;
template <typename T>
struct frame_cols_s : std::integral_constant<index_t, 0> {};
template <typename T>
constexpr auto frame_cols = frame_cols_s<std::remove_cvref_t<T>>::value;
template <typename T>
struct frame_rows_s : std::integral_constant<index_t, 0> {};
template <typename T>
constexpr auto frame_rows = frame_rows_s<std::remove_cvref_t<T>>::value;
template <typename T>
struct frame_etype_s {
  using type = void;
};
template <typename T>
using frame_etype = typename frame_etype_s<std::remove_cvref_t<T>>::type;

// Quat type traits
template <typename T>
struct is_quat_s : std::bool_constant<false> {};
template <typename T>
constexpr auto is_quat = is_quat_s<std::remove_cvref_t<T>>::value;
template <typename T>
struct quat_size_s : std::integral_constant<index_t, 0> {};
template <typename T>
constexpr auto quat_size = quat_size_s<std::remove_cvref_t<T>>::value;
template <typename T>
struct quat_etype_s {
  using type = void;
};
template <typename T>
using quat_etype = typename quat_etype_s<std::remove_cvref_t<T>>::type;

// Number type traits
template <typename T>
constexpr auto is_integral = std::is_integral_v<T>;
template <typename T>
constexpr auto is_real = std::is_floating_point_v<T>;
template <typename T>
constexpr auto is_number = std::is_integral_v<T> || std::is_floating_point_v<T>;
template <typename T>
constexpr auto is_bool = std::is_same_v<T, bool>;

// Vector type traits
template <typename T>
constexpr auto is_any_vec = is_vec<T>;
template <typename T>
constexpr auto is_int_vec = is_vec<T> && is_integral<vec_etype<T>>;
template <typename T>
constexpr auto is_real_vec = is_vec<T> && is_real<vec_etype<T>>;
template <typename T>
constexpr auto is_num_vec = is_vec<T> && is_number<vec_etype<T>>;
template <typename T>
constexpr auto is_bool_vec = is_vec<T> && is_bool<vec_etype<T>>;

// Matrix type traits
template <typename T>
constexpr auto is_real_mat = is_mat<T> && is_real<mat_etype<T>>;
template <typename T>
constexpr auto is_num_mat = is_mat<T> && is_number<mat_etype<T>>;

// Frame type traits
template <typename T>
constexpr auto is_real_frame = is_frame<T> && is_real<frame_etype<T>>;
template <typename T>
constexpr auto is_num_frame = is_frame<T> && is_number<frame_etype<T>>;

// Vector type traits
template <typename T>
constexpr auto is_real_quat = is_quat<T> && is_real<quat_etype<T>>;
template <typename T>
constexpr auto is_num_quat = is_quat<T> && is_number<quat_etype<T>>;

// Concepts
template <typename T>
concept integral = is_integral<T>;
template <typename T>
concept real = is_real<T>;
template <typename T>
concept num = is_number<T>;

template <typename T>
concept anyN = is_any_vec<T>;
template <typename T>
concept inteN = is_int_vec<T>;
template <typename T>
concept realN = is_real_vec<T>;
template <typename T>
concept numN = is_num_vec<T>;
template <typename T>
concept boolN = is_bool_vec<T>;

template <typename T>
concept num2 = is_num_vec<T> && vec_size<T> == 2;
template <typename T>
concept num3 = is_num_vec<T> && vec_size<T> == 3;
template <typename T>
concept num4 = is_num_vec<T> && vec_size<T> == 4;

template <typename T>
concept inte2 = is_int_vec<T> && vec_size<T> == 2;
template <typename T>
concept inte3 = is_int_vec<T> && vec_size<T> == 3;
template <typename T>
concept inte4 = is_int_vec<T> && vec_size<T> == 4;

template <typename T>
concept real2 = is_real_vec<T> && vec_size<T> == 2;
template <typename T>
concept real3 = is_real_vec<T> && vec_size<T> == 3;
template <typename T>
concept real4 = is_real_vec<T> && vec_size<T> == 4;

template <typename T>
concept realNxM = is_real_mat<T>;
template <typename T>
concept numNxM = is_num_mat<T>;

template <typename T>
concept realNxF = is_real_frame<T>;
template <typename T>
concept numNxF = is_num_frame<T>;

template <typename T>
concept num2xF = is_num_frame<T> && frame_rows<T> == 2;
template <typename T>
concept num3xF = is_num_frame<T> && frame_rows<T> == 3;

template <typename T>
concept numQ = is_num_quat<T> && quat_size<T> == 4;
template <typename T>
concept realQ = is_real_quat<T> && quat_size<T> == 4;

template <index_t N, typename... Ts>
concept same_size = (sizeof...(Ts) == N);

// Helpers
template <index_t... Is>
using index_seq = std::integer_sequence<index_t, Is...>;
template <index_t N>
using indices = std::make_integer_sequence<index_t, N>;

}  // namespace yocto

// -----------------------------------------------------------------------------
// MATH CONSTANTS AND FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

using byte   = unsigned char;
using uint   = unsigned int;
using ushort = unsigned short;
using index  = long;

// Common type
template <typename... Ts>
using common_type = std::common_type_t<Ts...>;
template <size_t... Is>
struct common_size_s;
template <size_t I>
struct common_size_s<I> : std::integral_constant<size_t, I> {};
template <size_t I, size_t... Is>
struct common_size_s<I, Is...> : std::integral_constant<size_t, I> {};
template <size_t... Is>
constexpr auto common_size = common_size_s<Is...>::value;
#ifndef __CUDACC__
template <typename Func, typename... Ts>
using result_t = std::invoke_result_t<Func, Ts...>;
#else
template <typename Func, typename... Ts>
using result_t = std::result_of_t<Func(Ts...)>;
#endif

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

constexpr kernel auto abs(num auto a) { return a < 0 ? -a : a; }
constexpr kernel auto min(num auto a, num auto b) {
  using T = common_type<decltype(a), decltype(b)>;
  return (a < b) ? (T)a : (T)b;
}
constexpr kernel auto max(num auto a, num auto b) {
  using T = common_type<decltype(a), decltype(b)>;
  return (a > b) ? (T)a : (T)b;
}
constexpr kernel auto clamp(num auto a, num auto min_, num auto max_) {
  using T = decltype(a);
  return min(max(a, (T)min_), (T)max_);
}
constexpr kernel auto sign(num auto a) {
  using T = decltype(a);
  return a < 0 ? (T)-1 : (T)1;
}
constexpr kernel auto sqr(num auto a) { return a * a; }
constexpr kernel auto sqrt(real auto a) { return std::sqrt(a); }
constexpr kernel auto sin(real auto a) { return std::sin(a); }
constexpr kernel auto cos(real auto a) { return std::cos(a); }
constexpr kernel auto tan(real auto a) { return std::tan(a); }
constexpr kernel auto asin(real auto a) { return std::asin(a); }
constexpr kernel auto acos(real auto a) { return std::acos(a); }
constexpr kernel auto atan(real auto a) { return std::atan(a); }
constexpr kernel auto log(real auto a) { return std::log(a); }
constexpr kernel auto exp(real auto a) { return std::exp(a); }
constexpr kernel auto log2(real auto a) { return std::log2(a); }
constexpr kernel auto exp2(real auto a) { return std::exp2(a); }
constexpr kernel auto pow(num auto a, num auto b) {
  using T = common_type<decltype(a), decltype(b)>;
  return std::pow((T)a, (T)b);
}
constexpr kernel bool isfinite(real auto a) { return std::isfinite(a); }
constexpr kernel auto atan2(num auto a, num auto b) {
  using T = common_type<decltype(a), decltype(b)>;
  return std::atan2((T)a, (T)b);
}
constexpr kernel auto round(real auto a) { return std::round(a); }
constexpr kernel auto fmod(real auto a, num auto b) {
  return std::fmod(a, (decltype(a))b);
}
constexpr kernel auto mod(real auto a, num auto b) {
  auto m = fmod(a, b);
  return (m >= 0) ? m : m + b;
}
constexpr kernel auto mod(integral auto a, integral auto b) {
  auto m = a % b;
  return (m >= 0) ? m : m + b;
}
// template <typename T>
// inline void swap(T& a, T& b) {
//   std::swap(a, b);
// }
constexpr kernel auto radians(real auto a) { return a * (decltype(a))pi / 180; }
constexpr kernel auto degrees(real auto a) { return a * 180 / (decltype(a))pi; }
constexpr kernel auto lerp(num auto a, num auto b, num auto u) {
  return a * (1 - u) + b * u;
}
constexpr kernel auto step(num auto a, num auto u) {
  using T = common_type<decltype(a), decltype(u)>;
  return u < a ? (T)0 : (T)1;
}
constexpr kernel auto smoothstep(num auto a, num auto b, num auto u) {
  using T = common_type<decltype(a), decltype(b), decltype(u)>;
  auto t  = clamp((u - a) / (b - a), (T)0, (T)1);
  return t * t * (3 - 2 * t);
}
constexpr kernel auto bias(num auto a, num auto bias) {
  return a / ((1 / bias - 2) * (1 - a) + 1);
}
constexpr kernel auto gain(num auto a, num auto gain) {
  using T = common_type<decltype(a), decltype(gain)>;
  return (a < (T)0.5) ? bias(a * 2, gain) / 2
                      : bias(a * 2 - 1, 1 - gain) / 2 + (T)0.5;
}
constexpr kernel auto pow2(integral auto a) { return 1 << a; }

}  // namespace yocto

// -----------------------------------------------------------------------------
// VECTORS
// -----------------------------------------------------------------------------
namespace yocto {

// Static vectors
template <typename T, index_t N>
struct vec {
  T d[N];

  constexpr kernel vec() : d{0} {}

  constexpr kernel explicit vec(T v) {
    for (auto& e : d) e = v;
  }
  template <num... Ts>
    requires same_size<N, Ts...>
  constexpr kernel vec(Ts... v) : d{(T)v...} {}

  constexpr kernel vec(const vec<T, N - 1>& first, T last) {
    for (auto idx = 0; idx < N - 1; idx++) d[idx] = (T)first.d[idx];
    d[N - 1] = last;
  }

  constexpr kernel      vec(const vec& v)       = default;
  constexpr kernel vec& operator=(const vec& v) = default;

  template <typename T1>
  constexpr kernel explicit(!std::is_convertible_v<T1, T>)
      vec(const vec<T1, N>& v)
      : vec(indices<N>(), v) {}
  template <typename T1>
  constexpr kernel explicit(!std::is_convertible_v<T, T1>)
  operator vec<T1, N>() {
    return to_vec<T1>(indices<N>());
  }

  constexpr kernel bool   empty() const { return false; }
  constexpr kernel size_t size() const { return N; }

  constexpr kernel T&       operator[](index_t i) { return d[i]; }
  constexpr kernel const T& operator[](index_t i) const { return d[i]; }

  constexpr kernel const T* begin() const { return d; }
  constexpr kernel const T* end() const { return d + N; }
  constexpr kernel T*       begin() { return d; }
  constexpr kernel T*       end() { return d + N; }

  constexpr kernel const T* data() const { return d; }
  constexpr kernel T*       data() { return d; }

  constexpr const T& front() const { return d[0]; }
  constexpr T&       front() { return d[0]; }
  constexpr const T& back() const { return d[N - 1]; }
  constexpr T&       back() { return d[N - 1]; }
  template <index_t M = N - 1>
  constexpr const vec<T, M>& first() const {
    static_assert(M < N);
    return *(vec<T, M>*)this;
  }
  template <index_t M = N - 1>
  constexpr vec<T, M>& first() {
    static_assert(M < N);
    return *(vec<T, M>*)this;
  }

  constexpr kernel T&       x() { return d[0]; }
  constexpr kernel const T& x() const { return d[0]; }
  constexpr kernel T&       y() { return d[1]; }
  constexpr kernel const T& y() const { return d[1]; }
  constexpr kernel T&       z() { return d[2]; }
  constexpr kernel const T& z() const { return d[2]; }
  constexpr kernel T&       w() { return d[3]; }
  constexpr kernel const T& w() const { return d[3]; }
  constexpr kernel vec<T, 2>& xy() { return *(vec<T, 2>)this; }
  constexpr kernel const vec<T, 2>& xy() const { return *(vec<T, 2>)this; }
  constexpr kernel vec<T, 3>& xyz() { return *(vec<T, 3>)this; }
  constexpr kernel const vec<T, 3>& xyz() const { return *(vec<T, 3>)this; }

  constexpr kernel T&       r() { return d[0]; }
  constexpr kernel const T& r() const { return d[0]; }
  constexpr kernel T&       g() { return d[1]; }
  constexpr kernel const T& g() const { return d[1]; }
  constexpr kernel T&       b() { return d[2]; }
  constexpr kernel const T& b() const { return d[2]; }
  constexpr kernel T&       a() { return d[3]; }
  constexpr kernel const T& a() const { return d[3]; }
  constexpr kernel vec<T, 3>& rgb() { return *(vec<T, 3>)this; }
  constexpr kernel const vec<T, 3>& rgb() const { return *(vec<T, 3>)this; }

 private:
  template <typename T1, index_t... Is>
  constexpr kernel vec(index_seq<Is...>, const vec<T1, N>& v)
      : d{(T)v[Is]...} {}
  template <typename T1, index_t... Is>
  constexpr kernel vec<T1, N> to_vec(index_seq<Is...>) {
    return {(T1)d[Is]...};
  }
};

template <typename T>
struct vec<T, 1> {
  T d[1];

  constexpr kernel vec() : d{0} {}
  constexpr kernel vec(T x_) : d{x_} {}

  template <typename U>
  constexpr kernel explicit(!std::is_convertible_v<U, T>)
      vec(const vec<U, 1>& v)
      : d{(T)v.d[0]} {}
  template <typename U>
  constexpr kernel explicit(!std::is_convertible_v<T, U>) operator vec<U, 1>() {
    return {(U)d[0]};
  }

  constexpr kernel vec(const std::array<T, 1>& v) : d{v[0]} {}
  constexpr kernel operator std::array<T, 1>() { return {d[0]}; }

  constexpr kernel bool   empty() const { return false; }
  constexpr kernel size_t size() const { return 1; }

  constexpr kernel T&       operator[](index_t i) { return d[i]; }
  constexpr kernel const T& operator[](index_t i) const { return d[i]; }

  constexpr kernel const T* begin() const { return d; }
  constexpr kernel const T* end() const { return d + 1; }
  constexpr kernel T*       begin() { return d; }
  constexpr kernel T*       end() { return d + 1; }

  constexpr kernel const T* data() const { return d; }
  constexpr kernel T*       data() { return d; }

  constexpr kernel T&       x() { return d[0]; }
  constexpr kernel const T& x() const { return d[0]; }
};

// Deduction guides
template <num... Args>
vec(Args...) -> vec<common_type<Args...>, sizeof...(Args)>;
template <numN V, num N>
vec(V, N) -> vec<common_type<vec_etype<V>, N>, vec_size<V> + 1>;

// Vector type traits
template <typename T, index_t N>
struct is_vec_s<vec<T, N>> : std::bool_constant<true> {};
template <typename T, index_t N>
struct vec_size_s<vec<T, N>> : std::integral_constant<index_t, N> {};
template <typename T, index_t N>
struct vec_etype_s<vec<T, N>> {
  using type = T;
};

// Vector aliases
template <typename T>
using vec2 = vec<T, 2>;
template <typename T>
using vec3 = vec<T, 3>;
template <typename T>
using vec4 = vec<T, 4>;

// Generic constants
template <typename T, int N>
constexpr auto zero_vec = vec<T, N>{0};
template <typename T, int N>
constexpr auto one_vec = vec<T, N>{1};

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

// Element access
template <typename T>
constexpr kernel vec<T, 2> xy(const vec<T, 3>& a) {
  return {a[0], a[1]};
}
template <typename T>
constexpr kernel vec<T, 2> xy(const vec<T, 4>& a) {
  return {a[0], a[1]};
}
template <typename T>
constexpr kernel vec<T, 3> xyz(const vec<T, 4>& a) {
  return {a[0], a[1], a[2]};
}

// Vector sequence operations.
constexpr kernel bool    empty(anyN auto const& a) { return false; }
constexpr kernel index_t size(anyN auto const& a) { return a.size(); }
constexpr kernel auto    begin(anyN auto const& a) { return a.begin(); }
constexpr kernel auto    end(anyN auto const& a) { return a.end(); }
constexpr kernel auto    begin(anyN auto& a) { return a.begin(); }
constexpr kernel auto    end(anyN auto& a) { return a.end(); }
constexpr kernel auto    data(anyN auto const& a) { return a.data(); }
constexpr kernel auto    data(anyN auto& a) { return a.data(); }
template <std::size_t I>
constexpr kernel auto const& get(anyN auto const& a) noexcept {
  return a.d[I];
}
template <std::size_t I>
constexpr kernel auto& get(anyN auto& a) noexcept {
  return a.d[I];
}
template <std::size_t I>
constexpr kernel auto const&& get(anyN auto const&& a) noexcept {
  return a.d[I];
}
template <std::size_t I>
constexpr kernel auto&& get(anyN auto&& a) noexcept {
  return a.d[I];
}

// Implementation of operations using map and fold
template <typename T1, index_t N, typename Func, index_t... Is,
    typename T = result_t<Func, T1>>
constexpr kernel vec<T, N> map(
    index_seq<Is...>, const vec<T1, N>& a, Func&& func) {
  return vec{func(a[Is])...};
}
template <typename T1, index_t N, typename Func,
    typename T = result_t<Func, T1>>
constexpr kernel vec<T, N> map(const vec<T1, N>& a, Func&& func) {
  return map(indices<N>(), a, std::forward<Func>(func));
}
template <typename T1, typename T2, index_t N, typename Func, index_t... Is,
    typename T = result_t<Func, T1, T2>>
constexpr kernel vec<T, N> map(
    index_seq<Is...>, const vec<T1, N>& a, const vec<T2, N>& b, Func&& func) {
  return {func(a[Is], b[Is])...};
}
template <typename T1, typename T2, index_t N, typename Func,
    typename T = result_t<Func, T1, T2>>
constexpr kernel vec<T, N> map(
    const vec<T1, N>& a, const vec<T2, N>& b, Func&& func) {
  return map(indices<N>(), a, b, std::forward<Func>(func));
}
template <typename T1, typename T2, index_t N, typename Func, index_t... Is,
    typename T = result_t<Func, T1, T2>>
constexpr kernel vec<T, N> map(
    index_seq<Is...>, const vec<T1, N>& a, T2 b, Func&& func) {
  return {func(a[Is], b)...};
}
template <typename T1, typename T2, index_t N, typename Func,
    typename T = result_t<Func, T1, T2>>
constexpr kernel vec<T, N> map(const vec<T1, N>& a, T2 b, Func&& func) {
  return map(indices<N>(), a, b, std::forward<Func>(func));
}
template <typename T1, typename T2, index_t N, typename Func, index_t... Is,
    typename T = result_t<Func, T1, T2>>
constexpr kernel vec<T, N> map(
    index_seq<Is...>, T1 a, const vec<T2, N>& b, Func&& func) {
  return {func(a, b[Is])...};
}
template <typename T1, typename T2, index_t N, typename Func,
    typename T = result_t<Func, T1, T2>>
constexpr kernel vec<T, N> map(T1 a, const vec<T2, N>& b, Func&& func) {
  return map(indices<N>(), a, b, std::forward<Func>(func));
}
template <typename T1, typename T2, typename T3, index_t N, typename Func,
    index_t... Is, typename T = result_t<Func, T1, T2, T3>>
constexpr kernel vec<T, N> map(index_seq<Is...>, const vec<T1, N>& a,
    const vec<T2, N>& b, const vec<T3, N>& c, Func&& func) {
  return {func(a[Is], b[Is], c[Is])...};
}
template <typename T1, typename T2, typename T3, index_t N, typename Func,
    typename T = result_t<Func, T1, T2, T3>>
constexpr kernel vec<T, N> map(const vec<T1, N>& a, const vec<T2, N>& b,
    const vec<T3, N>& c, Func&& func) {
  return map(indices<N>(), a, b, c, std::forward<Func>(func));
}
template <typename T1, typename T2, typename T3, index_t N, typename Func,
    index_t... Is, typename T = result_t<Func, T1, T2, T3>>
constexpr kernel vec<T, N> map(
    index_seq<Is...>, const vec<T1, N>& a, T2 b, T3 c, Func&& func) {
  return {func(a[Is], b, c)...};
}
template <typename T1, typename T2, typename T3, index_t N, typename Func,
    typename T = result_t<Func, T1, T2, T3>>
constexpr kernel vec<T, N> map(const vec<T1, N>& a, T2 b, T3 c, Func&& func) {
  return map(indices<N>(), a, b, c, std::forward<Func>(func));
}
template <typename T1, typename T2, typename T3, index_t N, typename Func,
    index_t... Is, typename T = result_t<Func, T1, T2, T3>>
constexpr kernel vec<T, N> map(index_seq<Is...>, const vec<T1, N>& a,
    const vec<T2, N>& b, T3 c, Func&& func) {
  return {func(a[Is], b[Is], c)...};
}
template <typename T1, typename T2, typename T3, index_t N, typename Func,
    typename T = result_t<Func, T1, T2, T3>>
constexpr kernel vec<T, N> map(
    const vec<T1, N>& a, const vec<T2, N>& b, T3 c, Func&& func) {
  return map(indices<N>(), a, b, c, std::forward<Func>(func));
}
template <typename T1, typename T2, typename T3, index_t N, typename Func,
    index_t... Is, typename T = result_t<Func, T1, T2, T3>>
constexpr kernel vec<T, N> map(index_seq<Is...>, const vec<T1, N>& a, T2 b,
    const vec<T3, N>& c, Func&& func) {
  return {func(a[Is], b, c[Is])...};
}
template <typename T1, typename T2, typename T3, index_t N, typename Func,
    typename T = result_t<Func, T1, T2, T3>>
constexpr kernel vec<T, N> map(
    const vec<T1, N>& a, T2 b, const vec<T3, N>& c, Func&& func) {
  return map(indices<N>(), a, b, c, std::forward<Func>(func));
}
template <typename T, index_t N, typename Func>
constexpr kernel T fold(const vec<T, N>& a, Func&& func) {
  if constexpr (N == 1) {
    return a[0];
  } else if constexpr (N == 2) {
    return func(a[0], a[1]);
  } else if constexpr (N == 3) {
    return func(func(a[0], a[1]), a[2]);
  } else if constexpr (N == 4) {
    return func(func(func(a[0], a[1]), a[2]), a[3]);
  } else {
    return func(
        fold((const vec<T, N - 1>&)a, std::forward<Func>(func)), a[N - 1]);
  }
}
template <typename T, index_t N, index_t... Is>
constexpr kernel T fold_sum(index_seq<Is...>, const vec<T, N>& a) {
  return (a[Is] + ...);
}
template <typename T, index_t N>
constexpr kernel T fold_sum(const vec<T, N>& a) {
  return fold_sum(indices<N>(), a);
}
template <typename T, index_t N, index_t... Is>
constexpr kernel T fold_prod(index_seq<Is...>, const vec<T, N>& a) {
  return (a[Is] * ...);
}
template <typename T, index_t N>
constexpr kernel T fold_prod(const vec<T, N>& a) {
  return fold_prod(indices<N>(), a);
}
template <index_t N, index_t... Is>
constexpr kernel bool fold_and(index_seq<Is...>, const vec<bool, N>& a) {
  return (... && a[Is]);
}
template <index_t N>
constexpr kernel bool fold_and(const vec<bool, N>& a) {
  return fold_and(indices<N>(), a);
}
template <index_t N, index_t... Is>
constexpr kernel bool fold_or(index_seq<Is...>, const vec<bool, N>& a) {
  return (... || a[Is]);
}
template <index_t N>
constexpr kernel bool fold_or(const vec<bool, N>& a) {
  return fold_and(indices<N>(), a);
}

// Vector comparison operations.
constexpr kernel bool operator==(numN auto const& a, numN auto const& b) {
  return fold_and(map(a, b, [](num auto a, num auto b) { return a == b; }));
}
constexpr kernel bool operator==(numN auto const& a, num auto b) {
  return fold_and(map(a, b, [](num auto a, num auto b) { return a == b; }));
}
constexpr kernel bool operator!=(numN auto const& a, numN auto const& b) {
  return fold_or(map(a, b, [](num auto a, num auto b) { return a != b; }));
}
constexpr kernel bool operator!=(numN auto const& a, num auto b) {
  return fold_or(map(a, b, [](num auto a, num auto b) { return a != b; }));
}

// Vector operations.
constexpr kernel auto operator+(numN auto const& a) { return a; }
constexpr kernel auto operator-(numN auto const& a) {
  return map(a, [](auto a) { return -a; });
}
constexpr kernel auto operator+(numN auto const& a, numN auto const& b) {
  return map(a, b, [](auto a, auto b) { return a + b; });
}
constexpr kernel auto operator+(numN auto const& a, num auto b) {
  return map(a, b, [](auto a, auto b) { return a + b; });
}
constexpr kernel auto operator+(num auto a, numN auto const& b) {
  return map(a, b, [](auto a, auto b) { return a + b; });
}
constexpr kernel auto operator-(numN auto const& a, numN auto const& b) {
  return map(a, b, [](auto a, auto b) { return a - b; });
}
constexpr kernel auto operator-(numN auto const& a, num auto b) {
  return map(a, b, [](auto a, auto b) { return a - b; });
}
constexpr kernel auto operator-(num auto a, numN auto const& b) {
  return map(a, b, [](auto a, auto b) { return a - b; });
}
constexpr kernel auto operator*(numN auto const& a, numN auto const& b) {
  return map(a, b, [](auto a, auto b) { return a * b; });
}
constexpr kernel auto operator*(numN auto const& a, num auto b) {
  return map(a, b, [](auto a, auto b) { return a * b; });
}
constexpr kernel auto operator*(num auto a, numN auto const& b) {
  return map(a, b, [](auto a, auto b) { return a * b; });
}
constexpr kernel auto operator/(numN auto const& a, numN auto const& b) {
  return map(a, b, [](auto a, auto b) { return a / b; });
}
constexpr kernel auto operator/(numN auto const& a, num auto b) {
  return map(a, b, [](auto a, auto b) { return a / b; });
}
constexpr kernel auto operator/(num auto a, numN auto const& b) {
  return map(a, b, [](auto a, auto b) { return a / b; });
}
constexpr kernel auto operator%(inteN auto const& a, inteN auto const& b) {
  return map(a, b, [](auto a, auto b) { return a % b; });
}
constexpr kernel auto operator%(inteN auto const& a, integral auto b) {
  return map(a, b, [](auto a, auto b) { return a % b; });
}
constexpr kernel auto operator%(num auto a, inteN auto const& b) {
  return map(a, b, [](auto a, auto b) { return a % b; });
}
constexpr kernel auto operator^(inteN auto const& a, inteN auto const& b) {
  return map(a, b, [](auto a, auto b) { return a ^ b; });
}
constexpr kernel auto operator^(inteN auto const& a, integral auto b) {
  return map(a, b, [](auto a, auto b) { return a ^ b; });
}
constexpr kernel auto operator^(integral auto a, inteN auto const& b) {
  return map(a, b, [](auto a, auto b) { return a ^ b; });
}
constexpr kernel auto operator>>(inteN auto const& a, inteN auto const& b) {
  return map(a, b, [](auto a, auto b) { return a >> b; });
}
constexpr kernel auto operator>>(inteN auto const& a, integral auto b) {
  return map(a, b, [](auto a, auto b) { return a >> b; });
}
constexpr kernel auto operator>>(integral auto a, inteN auto const& b) {
  return map(a, b, [](auto a, auto b) { return a >> b; });
}
constexpr kernel auto operator<<(inteN auto const& a, inteN auto const& b) {
  return map(a, b, [](auto a, auto b) { return a << b; });
}
constexpr kernel auto operator<<(inteN auto const& a, integral auto b) {
  return map(a, b, [](auto a, auto b) { return a << b; });
}
constexpr kernel auto operator<<(integral auto a, inteN auto const& b) {
  return map(a, b, [](auto a, auto b) { return a << b; });
}

// Vector assignments
constexpr kernel auto& operator+=(numN auto& a, numN auto const& b) {
  return a = a + b;
}
constexpr kernel auto& operator+=(numN auto& a, num auto b) {
  return a = a + b;
}
constexpr kernel auto& operator-=(numN auto& a, numN auto const& b) {
  return a = a - b;
}
constexpr kernel auto& operator-=(numN auto& a, num auto b) {
  return a = a - b;
}
constexpr kernel auto& operator*=(numN auto& a, numN auto const& b) {
  return a = a * b;
}
constexpr kernel auto& operator*=(numN auto& a, num auto b) {
  return a = a * b;
}
constexpr kernel auto& operator/=(numN auto& a, numN auto const& b) {
  return a = a / b;
}
constexpr kernel auto& operator/=(numN auto& a, num auto b) {
  return a = a / b;
}
constexpr kernel auto& operator%=(inteN auto& a, inteN auto const& b) {
  return a = a % b;
}
constexpr kernel auto& operator%=(inteN auto& a, integral auto b) {
  return a = a % b;
}

// Vector products and lengths.
constexpr kernel auto dot(numN auto const& a, numN auto const& b) {
  return sum(a * b);
}
constexpr kernel auto cross(num2 auto const& a, num2 auto const& b) {
  return a.x() * b.y() - a.y() * b.x();
}
constexpr kernel auto cross(num3 auto const& a, num3 auto const& b) {
  return vec{a.y() * b.z() - a.z() * b.y(), a.z() * b.x() - a.x() * b.z(),
      a.x() * b.y() - a.y() * b.x()};
}
template <typename T1, typename T2, typename T = common_type<T1, T2>>
constexpr kernel auto angle(numN auto const& a, numN auto const& b) {
  return acos(clamp(dot(normalize(a), normalize(b)), -1, 1));
}

// Orthogonal vectors.
constexpr kernel auto orthogonal(num3 auto const& v) {
  // http://lolengine.net/blog/2013/09/21/picking-orthogonal-vector-combing-coconuts)
  return abs(v.x()) > abs(v.z()) ? vec{-v.y(), v.x(), 0}
                                 : vec{0, -v.z(), v.y()};
}
constexpr kernel auto orthonormalize(numN auto const& a, numN auto const& b) {
  return normalize(a - b * dot(a, b));
}

// Reflected and refracted vector.
constexpr kernel auto reflect(numN auto const& w, numN auto const& n) {
  return -w + 2 * dot(n, w) * n;
}
constexpr kernel auto refract(
    num3 auto const& w, num3 auto const& n, num auto inv_eta) {
  using V     = decltype(w);
  auto cosine = dot(n, w);
  auto k      = 1 + inv_eta * inv_eta * (cosine * cosine - 1);
  if (k < 0) return V{0, 0, 0};
  return -w * inv_eta + (inv_eta * cosine - sqrt(k)) * n;
}

constexpr kernel auto length(numN auto const& a) { return sqrt(dot(a, a)); }
constexpr kernel auto length2(numN auto const& a) { return dot(a, a); }
[[deprecated]] constexpr kernel auto length_squared(numN auto const& a) {
  return dot(a, a);
}
constexpr kernel auto normalize(numN auto const& a) {
  auto l = length(a);
  return (l != 0) ? a / l : a;
}
constexpr kernel auto distance(numN auto const& a, numN auto const& b) {
  return length(a - b);
}
constexpr kernel auto distance2(numN auto const& a, numN auto const& b) {
  return dot(a - b, a - b);
}
[[deprecated]] constexpr kernel auto distance_squared(
    numN auto const& a, numN auto const& b) {
  return dot(a - b, a - b);
}
constexpr kernel auto angle(numN auto const& a, numN auto const& b) {
  return acos(clamp(dot(normalize(a), normalize(b)), -1, 1));
}

constexpr kernel auto slerp(
    num4 auto const& a, num4 auto const& b, num auto u) {
  using T = decltype(u);
  // https://en.wikipedia.org/wiki/Slerp
  auto an = normalize(a), bn = normalize(b);
  auto d = dot(an, bn);
  if (d < 0) {
    bn = -bn;
    d  = -d;
  }
  if (d > (T)0.9995) return normalize(an + u * (bn - an));
  auto th = acos(clamp(d, -1, 1));
  if (th == 0) return an;
  return an * (sin(th * (1 - u)) / sin(th)) + bn * (sin(th * u) / sin(th));
}

// Max element and clamp.
constexpr kernel auto max(numN auto const& a, numN auto const& b) {
  return map(a, b, [](num auto a, num auto b) { return max(a, b); });
}
constexpr kernel auto max(numN auto const& a, num auto b) {
  return map(a, b, [](num auto a, num auto b) { return max(a, b); });
}
constexpr kernel auto min(numN auto const& a, num auto b) {
  return map(a, b, [](num auto a, num auto b) { return min(a, b); });
}
constexpr kernel auto min(numN auto const& a, numN auto const& b) {
  return map(a, b, [](num auto a, num auto b) { return min(a, b); });
}
constexpr kernel auto clamp(
    numN auto const& x, numN auto const& min, numN auto const& max) {
  return map(x, min, max,
      [](num auto a, num auto b, num auto c) { return clamp(a, b, c); });
}
constexpr kernel auto clamp(numN auto const& x, num auto min, num auto max) {
  return map(x, min, max,
      [](num auto a, num auto b, num auto c) { return clamp(a, b, c); });
}
constexpr kernel auto clamp(
    numN auto const& x, numN auto const& min, num auto max) {
  return map(x, min, max,
      [](num auto a, num auto b, num auto c) { return clamp(a, b, c); });
}
constexpr kernel auto clamp(
    numN auto const& x, num auto min, numN auto const& max) {
  return map(x, min, max,
      [](num auto a, num auto b, num auto c) { return clamp(a, b, c); });
}
constexpr kernel auto radians(numN auto const& a) {
  return map(a, [](num auto a) { return radians(a); });
}
constexpr kernel auto degrees(numN auto const& a) {
  return map(a, [](num auto a) { return degrees(a); });
}
constexpr kernel auto lerp(numN auto const& a, numN auto const& b, num auto u) {
  return a * (1 - u) + b * u;
}
constexpr kernel auto lerp(
    numN auto const& a, numN auto const& b, numN auto const& u) {
  return a * (1 - u) + b * u;
}

constexpr kernel auto max(numN auto const& a) {
  return fold(a, [](num auto a, num auto b) { return max(a, b); });
}
constexpr kernel auto min(numN auto const& a) {
  return fold(a, [](num auto a, num auto b) { return min(a, b); });
}
constexpr kernel index_t argmax(numN auto const& a) {
  constexpr auto N = vec_size<decltype(a)>;
  if constexpr (N == 1) {
    return 0;
  } else if constexpr (N == 2) {
    return a[0] >= a[1] ? 0 : 1;
  } else if constexpr (N == 3) {
    return a[0] >= a[1] ? (a[0] >= a[2] ? 0 : 2) : (a[1] >= a[2] ? 1 : 2);
  } else if constexpr (N == 4) {
    if (a[3] >= a[0] && a[3] >= a[1] && a[3] >= a[2])
      return 3;
    else
      return argmax(xyz(a));
  }
}
constexpr kernel index_t argmin(numN auto const& a) {
  constexpr auto N = vec_size<decltype(a)>;
  if constexpr (N == 1) {
    return 0;
  } else if constexpr (N == 2) {
    return a[0] <= a[1] ? 0 : 1;
  } else if constexpr (N == 3) {
    return a[0] <= a[1] ? (a[0] <= a[2] ? 0 : 2) : (a[1] <= a[2] ? 1 : 2);
  } else if constexpr (N == 4) {
    if (a[3] <= a[0] && a[3] <= a[1] && a[3] <= a[2])
      return 3;
    else
      return argmin(xyz(a));
  }
}

constexpr kernel auto sum(numN auto const& a) { return fold_sum(a); }
constexpr kernel auto prod(numN auto const& a) { return fold_prod(a); }
constexpr kernel auto mean(numN auto const& a) { return sum(a) / a.size(); }
constexpr kernel bool all(boolN auto const& a) { return fold_and(a); }
constexpr kernel bool any(boolN auto const& a) { return fold_or(a); }

// Functions applied to vector elements
constexpr kernel auto abs(numN auto const& a) {
  return map(a, [](num auto a) { return abs(a); });
}
constexpr kernel auto sqr(numN auto const& a) {
  return map(a, [](num auto a) { return sqr(a); });
}
constexpr kernel auto sqrt(numN auto const& a) {
  return map(a, [](num auto a) { return sqrt(a); });
}
constexpr kernel auto exp(numN auto const& a) {
  return map(a, [](num auto a) { return exp(a); });
}
constexpr kernel auto log(numN auto const& a) {
  return map(a, [](num auto a) { return log(a); });
}
constexpr kernel auto exp2(numN auto const& a) {
  return map(a, [](num auto a) { return exp2(a); });
}
constexpr kernel auto log2(numN auto const& a) {
  return map(a, [](num auto a) { return log2(a); });
}
constexpr kernel auto pow(numN auto const& a, num auto b) {
  return map(a, b, [](num auto a, num auto b) { return pow(a, b); });
}
constexpr kernel auto pow(numN auto const& a, numN auto const& b) {
  return map(a, b, [](num auto a, num auto b) { return pow(a, b); });
}
constexpr kernel auto round(numN auto const& a) {
  return map(a, [](num auto a) { return round(a); });
}
constexpr kernel auto fmod(numN auto const& a, numN auto const& b) {
  return map(a, b, [](num auto a, num auto b) { return fmod(a, b); });
}
constexpr kernel auto fmod(numN auto const& a, num auto b) {
  return map(a, b, [](num auto a, num auto b) { return fmod(a, b); });
}
constexpr kernel auto mod(numN auto const& a, numN auto const& b) {
  return map(a, b, [](num auto a, num auto b) { return mod(a, b); });
}
constexpr kernel auto mod(numN auto const& a, num auto b) {
  return map(a, b, [](num auto a, num auto b) { return mod(a, b); });
}
constexpr kernel auto bias(numN auto const& a, num auto b) {
  return map(a, b, [](num auto a, num auto b) { return bias(a, b); });
}
constexpr kernel auto gain(numN auto const& a, num auto b) {
  return map(a, b, [](num auto a, num auto b) { return bias(a, b); });
}
constexpr kernel bool isfinite(numN auto const& a) {
  return all(map(a, [](num auto a) { return isfinite(a); }));
}

// Conversion between ranges
constexpr kernel auto unit_to_uv(num auto a) { return (a + 1) / 2; }
constexpr kernel auto unit_to_uv(numN auto const& a) { return (a + 1) / 2; }
constexpr kernel auto uv_to_unit(num auto a) { return a * 2 - 1; }
constexpr kernel auto uv_to_unit(numN auto const& a) { return a * 2 - 1; }

// Conversion between coordinates
constexpr kernel auto cartesian_to_sphericaluv(num3 auto const& w) {
  using T           = vec_etype<decltype(w)>;
  auto [wx, wy, wz] = w;
  auto uv = vec{atan2(wy, wx), acos(clamp(wz, -1, 1))} / vec{2 * (T)pi, (T)pi};
  return mod(uv, 1);
}
constexpr kernel auto cartesiany_to_sphericaluv(num3 auto const& w) {
  using T           = vec_etype<decltype(w)>;
  auto [wx, wy, wz] = w;
  auto uv = vec{atan2(wz, wx), acos(clamp(wy, -1, 1))} / vec{2 * (T)pi, (T)pi};
  return mod(uv, 1);
}
constexpr kernel auto sphericaluv_to_cartesian(num2 auto const& uv) {
  using T           = vec_etype<decltype(uv)>;
  auto [phi, theta] = uv * vec{2 * (T)pi, (T)pi};
  return vec{cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta)};
}
constexpr kernel auto sphericaluv_to_cartesiany(num2 auto const& uv) {
  using T           = vec_etype<decltype(uv)>;
  auto [phi, theta] = uv * vec{2 * (T)pi, (T)pi};
  return vec{cos(phi) * sin(theta), cos(theta), sin(phi) * sin(theta)};
}

// Quaternion operations represented as xi + yj + zk + w
// const auto identity_quat4f = vec4f{0, 0, 0, 1};
constexpr kernel auto quat_mul(num4 auto const a, num auto b) {
  return map(a, b, [](num auto a, num auto b) { return a * b; });
}
constexpr kernel auto quat_mul(num4 auto const a, num4 auto const& b) {
  return vec{a.x() * b.w() + a.w() * b.x() + a.y() * b.w() - a.z() * b.y(),
      a.y() * b.w() + a.w() * b.y() + a.z() * b.x() - a.x() * b.z(),
      a.z() * b.w() + a.w() * b.z() + a.x() * b.y() - a.y() * b.x(),
      a.w() * b.w() - a.x() * b.x() - a.y() * b.y() - a.z() * b.z()};
}
constexpr kernel auto quat_conjugate(num4 auto const& a) {
  return vec{-xyz(a), w(a)};
}
constexpr kernel auto quat_inverse(num4 auto const a) {
  return quat_conjugate(a) / dot(a, a);
}

// Component-wise comparison operations.
constexpr kernel auto component_equal(numN auto const& a, numN auto const& b) {
  return map(a, b, [](num auto a, num auto b) { return a == b; });
}
constexpr kernel auto component_equal(numN auto const& a, num auto b) {
  return map(a, b, [](num auto a, num auto b) { return a == b; });
}
constexpr kernel auto component_not_equal(
    numN auto const& a, numN auto const& b) {
  return map(a, b, [](num auto a, num auto b) { return a != b; });
}
constexpr kernel auto component_not_equal(
    numN auto const& a, const num auto b) {
  return map(a, b, [](num auto a, num auto b) { return a != b; });
}
constexpr kernel auto component_less(numN auto const& a, numN auto const& b) {
  return map(a, b, [](num auto a, num auto b) { return a < b; });
}
constexpr kernel auto component_less(numN auto const& a, const num auto b) {
  return map(a, b, [](num auto a, num auto b) { return a < b; });
}
constexpr kernel auto component_greater(
    numN auto const& a, numN auto const& b) {
  return map(a, b, [](num auto a, num auto b) { return a > b; });
}
constexpr kernel auto component_greater(numN auto const& a, const num auto b) {
  return map(a, b, [](num auto a, num auto b) { return a > b; });
}
constexpr kernel auto component_less_equal(
    numN auto const& a, numN auto const& b) {
  return map(a, b, [](num auto a, num auto b) { return a <= b; });
}
constexpr kernel auto component_less_equal(
    numN auto const& a, const num auto b) {
  return map(a, b, [](num auto a, num auto b) { return a <= b; });
}
constexpr kernel auto component_greater_equal(
    numN auto const& a, numN auto const& b) {
  return map(a, b, [](num auto a, num auto b) { return a >= b; });
}
constexpr kernel auto component_greater_equal(
    numN auto const& a, const num auto b) {
  return map(a, b, [](num auto a, num auto b) { return a >= b; });
}
constexpr kernel auto select(
    boolN auto const& a, numN auto const& b, numN auto const& c) {
  return map(a, b, c, [](bool a, num auto b, num auto c) { return a ? b : c; });
}
constexpr kernel auto select(
    boolN auto const& a, num auto b, numN auto const& c) {
  return map(a, b, c, [](bool a, num auto b, num auto c) { return a ? b : c; });
}
constexpr kernel auto select(
    boolN auto const& a, numN auto const& b, num auto c) {
  return map(a, b, c, [](bool a, num auto b, num auto c) { return a ? b : c; });
}
constexpr kernel auto select(boolN auto const& a, num auto b, num auto c) {
  return map(a, b, c, [](bool a, num auto b, num auto c) { return a ? b : c; });
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// CUSTOMIZATION POINTS FOR STD
// -----------------------------------------------------------------------------
namespace std {

// Structured binding support
template <typename T, size_t N>
struct tuple_size<yocto::vec<T, N>> : std::integral_constant<size_t, N> {};
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
template <typename T, index_t N, index_t M>
struct mat {
  vec<T, N> d[M];

  constexpr kernel mat() {
    for (auto j = 0; j < M; j++) {
      for (auto i = 0; i < N; i++) {
        d[j][i] = i == j ? 1 : 0;
      }
    }
  }
  constexpr kernel mat(T v) {
    for (auto j = 0; j < M; j++) {
      for (auto i = 0; i < N; i++) {
        d[j][i] = i == j ? v : 0;
      }
    }
  }

  // TODO: make generic
  constexpr kernel mat(const vec<T, N>& c1, const vec<T, N>& c2) : d{c1, c2} {}
  constexpr kernel mat(
      const vec<T, N>& c1, const vec<T, N>& c2, const vec<T, N>& c3)
      : d{c1, c2, c3} {}
  constexpr kernel mat(const vec<T, N>& c1, const vec<T, N>& c2,
      const vec<T, N>& c3, const vec<T, N>& c4)
      : d{c1, c2, c3, c4} {}

  constexpr kernel vec<T, N>& operator[](index_t i) { return d[i]; }
  constexpr kernel const vec<T, N>& operator[](index_t i) const { return d[i]; }
};

// Deduction guides
template <numN... Args>
mat(Args...) -> mat<common_type<vec_etype<Args>...>,
    common_size<vec_size<Args>...>, sizeof...(Args)>;

// Matrix type traits
template <typename T, index_t N, index_t M>
struct is_mat_s<mat<T, N, M>> : std::bool_constant<true> {};
template <typename T, index_t N, index_t M>
struct mat_rows_s<mat<T, N, M>> : std::integral_constant<index_t, N> {};
template <typename T, index_t N, index_t M>
struct mat_cols_s<mat<T, N, M>> : std::integral_constant<index_t, M> {};
template <typename T, index_t N, index_t M>
struct mat_etype_s<mat<T, N, M>> {
  using type = T;
};

// Matrix aliases
using mat1x1f = mat<float, 1, 1>;
using mat2x2f = mat<float, 2, 2>;
using mat3x3f = mat<float, 3, 3>;
using mat4x4f = mat<float, 4, 4>;

// Identity matrices constants.
constexpr auto identity2x2f = mat2x2f{{1, 0}, {0, 1}};
constexpr auto identity3x3f = mat3x3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}};
constexpr auto identity4x4f = mat4x4f{
    {1, 0, 0, 0}, {0, 1, 0, 0}, {0, 0, 1, 0}, {0, 0, 0, 1}};

// Data access
template <typename T, index_t N, index_t M>
constexpr kernel T* data(mat<T, N, M>& a) {
  return data(a[0]);
}

// Implementation of operations using map and fold
template <typename T1, index_t N, index_t M, typename Func, index_t... Is,
    typename T = result_t<Func, vec<T1, N>>>
constexpr kernel mat<T, N, M> map(
    index_seq<Is...>, const mat<T1, N, M>& a, Func&& func) {
  return mat{func(a[Is])...};
}
template <typename T1, index_t N, index_t M, typename Func,
    typename T = vec_etype<result_t<Func, vec<T1, N>>>>
constexpr kernel mat<T, N, M> map(const mat<T1, N, M>& a, Func&& func) {
  return map(indices<M>(), a, std::forward<Func>(func));
}
template <typename T1, typename T2, index_t N, index_t M, typename Func,
    index_t... Is,
    typename T = vec_etype<result_t<Func, vec<T1, N>, vec<T2, N>>>>
constexpr kernel mat<T, N, M> map(index_seq<Is...>, const mat<T1, N, M>& a,
    const mat<T2, N, M>& b, Func&& func) {
  return {func(a[Is], b[Is])...};
}
template <typename T1, typename T2, index_t N, index_t M, typename Func,
    typename T = vec_etype<result_t<Func, vec<T2, N>, vec<T2, N>>>>
constexpr kernel mat<T, N, M> map(
    const mat<T1, N, M>& a, const mat<T2, N, M>& b, Func&& func) {
  return map(indices<M>(), a, b, std::forward<Func>(func));
}
template <typename T1, typename T2, index_t N, index_t M, typename Func,
    index_t... Is, typename T = vec_etype<result_t<Func, vec<T1, N>, T2>>>
constexpr kernel mat<T, N, M> map(
    index_seq<Is...>, const mat<T1, N, M>& a, T2 b, Func&& func) {
  return {func(a[Is], b)...};
}
template <typename T1, typename T2, index_t N, index_t M, typename Func,
    typename T = vec_etype<result_t<Func, vec<T2, N>, T2>>>
constexpr kernel mat<T, N, M> map(const mat<T1, N, M>& a, T2 b, Func&& func) {
  return map(indices<M>(), a, b, std::forward<Func>(func));
}
template <typename T1, typename T2, index_t N, index_t M, index_t... Is,
    typename T = common_type<T1, T2>>
constexpr kernel vec<T, N> mul(
    index_seq<Is...>, const mat<T1, N, M>& a, const vec<T2, M>& b) {
  return {((a[Is] * b[Is]) + ...)};
}
template <typename T1, typename T2, index_t N, index_t M,
    typename T = common_type<T1, T2>>
constexpr kernel vec<T, N> mul(const mat<T1, N, M>& a, const vec<T2, M>& b) {
  return mul(indices<M>(), a, b);
}
template <typename T1, typename T2, index_t N, index_t M, index_t... Is,
    typename T = common_type<T1, T2>>
constexpr kernel vec<T, M> mul(
    index_seq<Is...>, const vec<T1, N>& a, const mat<T2, N, M>& b) {
  return {dot(a, b[Is])...};
}
template <typename T1, typename T2, index_t N, index_t M,
    typename T = common_type<T1, T2>>
constexpr kernel vec<T, M> mul(const vec<T1, N>& a, const mat<T2, N, M>& b) {
  return mul(indices<M>(), a, b);
}
template <typename T1, typename T2, index_t N, index_t M, index_t K,
    index_t... Is, typename T = common_type<T1, T2>>
constexpr kernel mat<T, N, M> mul(
    index_seq<Is...>, const mat<T1, N, K>& a, const mat<T2, K, M>& b) {
  return {(a * b[Is])...};
}
template <typename T1, typename T2, index_t N, index_t M, index_t K,
    typename T = common_type<T1, T2>>
constexpr kernel mat<T, N, M> mul(
    const mat<T1, N, K>& a, const mat<T2, K, M>& b) {
  return mul(indices<M>(), a, b);
}

// Matrix comparisons.
constexpr kernel bool operator==(numNxM auto const& a, numNxM auto const& b) {
  return fold_and(
      map(a, b, [](numN auto const& a, numN auto const& b) { return a == b; }));
}
constexpr kernel bool operator!=(numNxM auto const& a, numNxM auto const& b) {
  return fold_or(
      map(a, b, [](numN auto const& a, numN auto const& b) { return a != b; }));
}

// Matrix operations.
constexpr kernel auto operator+(numNxM auto const& a, numNxM auto const& b) {
  return map(
      a, b, [](numN auto const& a, numN auto const& b) { return a + b; });
}
constexpr kernel auto operator*(numNxM auto const& a, num auto b) {
  return map(a, b, [](numN auto const& a, num auto b) { return a * b; });
}
constexpr kernel auto operator*(numNxM auto const& a, numN auto const& b) {
  return mul(a, b);
}
constexpr kernel auto operator*(numN auto const& a, numNxM auto const& b) {
  return mul(a, b);
}
constexpr kernel auto operator*(numNxM auto const& a, numNxM auto const& b) {
  return mul(a, b);
}

// Matrix assignments.
constexpr kernel auto& operator+=(numNxM auto& a, numNxM auto const& b) {
  return a = a + b;
}
constexpr kernel auto& operator*=(numNxM auto& a, numNxM auto const& b) {
  return a = a * b;
}
constexpr kernel auto operator*=(numNxM auto& a, num auto b) {
  return a = a * b;
}

// Matrix diagonals and transposes.
constexpr kernel auto diagonal(numNxM auto const& a) -> numN auto{
  constexpr auto N = mat_rows<decltype(a)>;
  constexpr auto M = mat_cols<decltype(a)>;
  static_assert(N == M && N <= 4);
  if constexpr (N == 1) {
    return vec{a[0][0]};
  } else if constexpr (N == 2) {
    return vec{a[0][0], a[1][1]};
  } else if constexpr (N == 3) {
    return vec{a[0][0], a[1][1], a[2][2]};
  } else if constexpr (N == 4) {
    return vec{a[0][0], a[1][1], a[2][2], a[3][3]};
  }
}
constexpr kernel auto transpose(numNxM auto const& a) -> numNxM auto{
  using T          = mat_etype<decltype(a)>;
  constexpr auto N = mat_rows<decltype(a)>;
  constexpr auto M = mat_cols<decltype(a)>;
  static_assert(N == M);
  if constexpr (N == 1) {
    return mat<T, N, N>{{a[0][0]}};
  } else if constexpr (N == 2) {
    return mat<T, N, N>{{a[0][0], a[1][0]}, {a[0][1], a[1][1]}};
  } else if constexpr (N == 3) {
    return mat<T, N, N>{
        {a[0][0], a[1][0], a[2][0]},
        {a[0][1], a[1][1], a[2][1]},
        {a[0][2], a[1][2], a[2][2]},
    };
  } else if constexpr (N == 4) {
    return mat<T, N, N>{
        {a[0][0], a[1][0], a[2][0], a[3][0]},
        {a[0][1], a[1][1], a[2][1], a[3][1]},
        {a[0][2], a[1][2], a[2][2], a[3][2]},
        {a[0][3], a[1][3], a[2][3], a[3][3]},
    };
  }
}

// Matrix adjoints, determinants and inverses.
constexpr kernel auto determinant(numNxM auto const& a) -> num auto{
  constexpr auto N = mat_rows<decltype(a)>;
  constexpr auto M = mat_cols<decltype(a)>;
  static_assert(N == M);
  if constexpr (N == 1) {
    return a[0];
  } else if constexpr (N == 2) {
    return cross(a[0], a[1]);
  } else if constexpr (N == 3) {
    return dot(a[0], cross(a[1], a[2]));
  } else if constexpr (N == 4) {
    return 0;  // TODO
  }
}
constexpr kernel auto adjoint(numNxM auto const& a) -> numNxM auto{
  using T          = mat_etype<decltype(a)>;
  constexpr auto N = mat_rows<decltype(a)>;
  constexpr auto M = mat_cols<decltype(a)>;
  static_assert(N == M && N < 4);
  if constexpr (N == 1) {
    return mat<T, N, N>{{a[0][0]}};
  } else if constexpr (N == 2) {
    return mat<T, N, N>{{a[1][1], -a[0][1]}, {-a[1][0], a[0][0]}};
  } else if constexpr (N == 3) {
    return transpose(
        mat<T, N, N>{cross(a[1], a[2]), cross(a[2], a[0]), cross(a[0], a[1])});
  }
}
constexpr kernel auto inverse(numNxM auto const& a) -> numNxM auto{
  constexpr auto N = mat_rows<decltype(a)>;
  constexpr auto M = mat_cols<decltype(a)>;
  static_assert(N == M);
  return adjoint(a) * (1 / determinant(a));
}

// Constructs a basis from a direction
constexpr kernel auto basis_fromz(num3 auto const& v) -> numNxM auto{
  using T = vec_etype<decltype(v)>;
  // https://graphics.pixar.com/library/OrthonormalB/paper.pdf
  if constexpr (std::is_same_v<T, float>) {
    auto z    = normalize(v);
    auto sign = copysignf(1.0f, z.z());
    auto a    = -1.0f / (sign + z.z());
    auto b    = z.x() * z.y() * a;
    auto x    = vec<T, 3>{
        1.0f + sign * z.x() * z.x() * a, sign * b, -sign * z.x()};
    auto y = vec<T, 3>{b, sign + z.y() * z.y() * a, -z.y()};
    return mat<T, 3, 3>{x, y, z};
  } else if constexpr (std::is_same_v<T, double>) {
    // TODO: double
    return mat<T, 3, 3>{};
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// RIGID BODY TRANSFORMS/FRAMES
// -----------------------------------------------------------------------------
namespace yocto {

// Rigid frames stored as a column-major affine transform matrix.
template <typename T, index_t N>
struct frame;

// Rigid frames stored as a column-major affine transform matrix.
template <typename T>
struct frame<T, 2> {
  vec<T, 2> d[3];

  constexpr kernel frame() : d{{1, 0}, {0, 1}, {0, 0}} {}
  constexpr kernel frame(
      const vec<T, 2>& x_, const vec<T, 2>& y_, const vec<T, 2>& o_)
      : d{x_, y_, o_} {}
  constexpr kernel frame(const mat<T, 2, 2>& xy_, const vec<T, 2>& o_)
      : d{xy_[0], xy_[1], o_} {}

  explicit constexpr kernel frame(const mat<T, 3, 3>& m)
      : d{xy(m[0]), xy(m[1]), xy(m[2])} {}
  explicit constexpr kernel operator mat<T, 3, 3>() const {
    return {{d[0], 0}, {d[1], 0}, {d[2], 1}};
  }

  constexpr kernel vec<T, 2>& operator[](index_t i) { return d[i]; }
  constexpr kernel const vec<T, 2>& operator[](index_t i) const { return d[i]; }

  constexpr kernel vec<T, 2>& x() { return d[0]; }
  constexpr kernel const vec<T, 2>& x() const { return d[0]; }
  constexpr kernel vec<T, 2>& y() { return d[1]; }
  constexpr kernel const vec<T, 2>& y() const { return d[1]; }
  constexpr kernel vec<T, 2>& o() { return d[2]; }
  constexpr kernel const vec<T, 2>& o() const { return d[2]; }

  constexpr kernel mat<T, 2, 2>& m() { return d[0]; }
  constexpr kernel const mat<T, 2, 2>& m() const { return d[0]; }
  constexpr kernel vec<T, 2>& t() { return d[2]; }
  constexpr kernel const vec<T, 2>& t() const { return d[2]; }
};

// Rigid frames stored as a column-major affine transform matrix.
template <typename T>
struct frame<T, 3> {
  vec<T, 3> d[4];

  constexpr kernel frame() : d{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}} {}
  constexpr kernel frame(const vec<T, 3>& x_, const vec<T, 3>& y_,
      const vec<T, 3>& z_, const vec<T, 3>& o_)
      : d{x_, y_, z_, o_} {}
  constexpr kernel frame(const mat<T, 3, 3>& xyz_, const vec<T, 3>& o_)
      : d{xyz_[0], xyz_[1], xyz_[2], o_} {}

  explicit constexpr kernel frame(const mat<T, 4, 4>& m)
      : d{xyz(m[0]), xyz(m[1]), xyz(m[2]), xyz(m[3])} {}
  explicit constexpr kernel operator mat<T, 4, 4>() const {
    return {{d[0], 0}, {d[1], 0}, {d[2], 0}, {d[3], 1}};
  }

  constexpr kernel vec<T, 3>& operator[](index_t i) { return d[i]; }
  constexpr kernel const vec<T, 3>& operator[](index_t i) const { return d[i]; }

  constexpr kernel vec<T, 3>& x() { return d[0]; }
  constexpr kernel const vec<T, 3>& x() const { return d[0]; }
  constexpr kernel vec<T, 3>& y() { return d[1]; }
  constexpr kernel const vec<T, 3>& y() const { return d[1]; }
  constexpr kernel vec<T, 3>& z() { return d[2]; }
  constexpr kernel const vec<T, 3>& z() const { return d[2]; }
  constexpr kernel vec<T, 3>& o() { return d[3]; }
  constexpr kernel const vec<T, 3>& o() const { return d[3]; }

  constexpr kernel mat<T, 3, 3>& m() { return (mat<T, 3, 3>&)d[0]; }
  constexpr kernel const mat<T, 3, 3>& m() const { return (mat<T, 3, 3>&)d[0]; }
  constexpr kernel vec<T, 3>& t() { return d[3]; }
  constexpr kernel const vec<T, 3>& t() const { return d[3]; }
};

// Deduction guides
template <numN... Args>
frame(Args...) -> frame<common_type<vec_etype<Args>...>, sizeof...(Args) - 1>;
template <numNxM M, numN V>
frame(M, V) -> frame<common_type<mat_etype<M>, vec_etype<V>>, vec_size<V>>;

// Frame type traits
template <typename T, index_t N>
struct is_frame_s<frame<T, N>> : std::bool_constant<true> {};
template <typename T, index_t N>
struct frame_rows_s<frame<T, N>> : std::integral_constant<index_t, N> {};
template <typename T, index_t N>
struct frame_cols_s<frame<T, N>> : std::integral_constant<index_t, N + 1> {};
template <typename T, index_t N>
struct frame_etype_s<frame<T, N>> {
  using type = T;
};

// Frame aliases
using frame2f = frame<float, 2>;
using frame3f = frame<float, 3>;

// Identity frames.
constexpr auto identity2x3f = frame2f{{1, 0}, {0, 1}, {0, 0}};
constexpr auto identity3x4f = frame3f{
    {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 0}};

// Frame properties
constexpr kernel auto rotation(numNxF auto const& a) { return a.m(); }
constexpr kernel auto translation(numNxF auto const& a) { return a.t(); }

// Frame/mat conversion
constexpr kernel auto mat_to_frame(numNxM auto const& m) {
  using T          = mat_etype<decltype(m)>;
  constexpr auto N = mat_rows<decltype(m)>;
  return (frame<T, N - 1>)m;
}
constexpr kernel auto frame_to_mat(numNxF auto const& f) {
  using T          = frame_etype<decltype(f)>;
  constexpr auto N = frame_rows<decltype(f)>;
  return (mat<T, N + 1, N + 1>)f;
}

// Frame comparisons.
template <typename T1, typename T2, index_t N>
constexpr kernel bool operator==(const frame<T1, N>& a, const frame<T2, N>& b) {
  if constexpr (N == 2) {
    return a[0] == b[0] && a[1] == b[1] && a[2] == b[2];
  } else if constexpr (N == 3) {
    return a[0] == b[0] && a[1] == b[1] && a[2] == b[2] && a[3] == b[3];
  }
}
template <typename T1, typename T2, index_t N>
constexpr kernel bool operator!=(const frame<T1, N>& a, const frame<T2, N>& b) {
  return !(a == b);
}

// Frame composition, equivalent to affine matrix product.
constexpr kernel auto operator*(numNxF auto const& a, numNxF auto const& b) {
  return frame{a.m() * b.m(), a.m() * b.t() + a.t()};
}
constexpr kernel auto operator*=(numNxF auto& a, numNxF auto const& b) {
  return a = a * b;
}

// Frame inverse, equivalent to rigid affine inverse.
constexpr kernel auto inverse(numNxF auto const& a, bool non_rigid = false) {
  if (non_rigid) {
    auto minv = inverse(a.m());
    return frame{minv, -(minv * a.t())};
  } else {
    auto minv = transpose(a.m());
    return frame{minv, -(minv * a.t())};
  }
}

// Frame construction from axis.
constexpr kernel auto frame_fromz(num3 auto const& o, num3 auto const& v) {
  using T = vec_etype<common_type<decltype(o), decltype(v)>>;
  // https://graphics.pixar.com/library/OrthonormalB/paper.pdf
  if constexpr (std::is_same_v<T, float>) {
    auto z    = normalize(v);
    auto sign = copysignf((T)1, z.z());
    auto a    = -1 / (sign + z.z());
    auto b    = z.x() * z.y() * a;
    auto x = vec<T, 3>{1 + sign * z.x() * z.x() * a, sign * b, -sign * z.x()};
    auto y = vec<T, 3>{b, sign + z.y() * z.y() * a, -z.y()};
    return frame<T, 3>{x, y, z, o};
  } else if constexpr (std::is_same_v<T, double>) {
    // TODO: double
    return frame<T, 3>{};
  }
}
constexpr kernel auto frame_fromzx(
    num3 auto const& o, num3 auto const& z_, num3 auto const& x_) {
  auto z = normalize(z_);
  auto x = orthonormalize(x_, z);
  auto y = normalize(cross(z, x));
  return frame{x, y, z, o};
}
constexpr kernel auto orthonormalize(num3xF auto const& frame_) {
  auto z = normalize(frame_.z());
  auto x = orthonormalize(frame_.x(), z);
  auto y = normalize(cross(z, x));
  auto o = frame_.o();
  return frame{x, y, z, o};
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// QUATERNIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Quaternions to represent 3D rotations
template <typename T, index_t N>
struct quat;

// Quaternions to represent rotations
template <typename T>
struct quat<T, 4> {
  T d[4];

  constexpr kernel quat() : d{0, 0, 0, 1} {}
  constexpr kernel quat(T x, T y, T z, T w) : d{x, y, z, w} {}

  constexpr kernel T&       operator[](index_t i) { return d[i]; }
  constexpr kernel const T& operator[](index_t i) const { return d[i]; }
};

// Deduction guides
template <num... Args>
quat(Args...) -> quat<common_type<Args...>, sizeof...(Args)>;

// Vector type traits
template <typename T, index_t N>
struct is_quat_s<quat<T, N>> : std::bool_constant<true> {};
template <typename T, index_t N>
struct quat_size_s<quat<T, N>> : std::integral_constant<index_t, N> {};
template <typename T, index_t N>
struct quat_etype_s<quat<T, N>> {
  using type = T;
};

// Quaternion aliases
using quat4f = quat<float, 4>;

[[deprecated]] constexpr auto identity_quat4f = quat4f{0, 0, 0, 1};

template <typename T1, typename T2, typename Func,
    typename T = result_t<Func, T1, T2>>
constexpr kernel quat<T, 4> map(
    const quat<T1, 4>& a, const quat<T2, 4>& b, Func&& func) {
  return {
      func(a[0], b[0]), func(a[1], b[1]), func(a[2], b[2]), func(a[3], b[3])};
}
template <typename T1, typename T2, typename Func,
    typename T = result_t<Func, T1, T2>>
constexpr kernel quat<T, 4> map(const quat<T1, 4>& a, T2 b, Func&& func) {
  return {func(a[0], b), func(a[1], b), func(a[2], b), func(a[3], b)};
}

// Quaternion operations
constexpr kernel auto operator+(numQ auto const& a, numQ auto const& b) {
  return map(a, b, [](num auto a, num auto b) { return a + b; });
}
constexpr kernel auto operator*(numQ auto const& a, num auto b) {
  return map(a, b, [](num auto a, num auto b) { return a * b; });
}
constexpr kernel auto operator/(numQ auto const& a, num auto b) {
  return map(a, b, [](num auto a, num auto b) { return a / b; });
}
constexpr kernel auto operator*(numQ auto const& a, numQ auto const& b) {
  return quat{a.x() * b.w() + a.w() * b.x() + a.y() * b.w() - a.z() * b.y(),
      a.y() * b.w() + a.w() * b.y() + a.z() * b.x() - a.x() * b.z(),
      a.z() * b.w() + a.w() * b.z() + a.x() * b.y() - a.y() * b.x(),
      a.w() * b.w() - a.x() * b.x() - a.y() * b.y() - a.z() * b.z()};
}

// Quaternion operations
constexpr kernel auto dot(numQ auto const& a, numQ auto const& b) {
  return a.x() * b.x() + a.y() * b.y() + a.z() * b.z() + a.w() * b.w();
}
constexpr kernel auto length(numQ auto const& a) { return sqrt(dot(a, a)); }
constexpr kernel auto normalize(numQ auto const& a) {
  auto l = length(a);
  return (l != 0) ? a / l : a;
}
constexpr kernel auto conjugate(numQ auto const& a) {
  return quat{-a.xyz(), a.w()};
}
constexpr kernel auto inverse(numQ auto const& a) {
  return conjugate(a) / dot(a, a);
}
constexpr kernel auto uangle(numQ auto const& a, numQ auto const& b) {
  auto d = dot(a, b);
  return d > 1 ? 0 : acos(d < -1 ? -1 : d);
}
constexpr kernel auto lerp(numQ auto const& a, numQ auto const& b, num auto t) {
  return a * (1 - t) + b * t;
}
constexpr kernel auto nlerp(
    numQ auto const& a, numQ auto const& b, num auto t) {
  return normalize(lerp(a, b, t));
}
constexpr kernel auto slerp(
    numQ auto const& a, numQ auto const& b, num auto t) {
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
constexpr kernel auto transform_point(numNxM auto const& a, numN auto const& b)
  requires(mat_rows<decltype(a)> == vec_size<decltype(b)> + 1)
{
  auto tb = a * vec{b, 1};
  return tb.first() / tb.back();
}
constexpr kernel auto transform_vector(numNxM auto const& a, numN auto const& b)
  requires(mat_rows<decltype(a)> == vec_size<decltype(b)> + 1)
{
  auto tb = a * vec{b, 0};
  return tb.first() / tb.back();
}
constexpr kernel auto transform_direction(
    numNxM auto const& a, numN auto const& b)
  requires(mat_rows<decltype(a)> == vec_size<decltype(b)> + 1)
{
  return normalize(transform_vector(a, b));
}
constexpr kernel auto transform_normal(numNxM auto const& a, numN auto const& b)
  requires(mat_rows<decltype(a)> == vec_size<decltype(b)> + 1)
{
  return normalize(transform_vector(transpose(inverse(a)), b));
}
constexpr kernel auto transform_vector(numNxM auto const& a, numN auto const& b)
  requires(mat_rows<decltype(a)> == vec_size<decltype(b)>)
{
  return a * b;
}
constexpr kernel auto transform_direction(
    numNxM auto const& a, numN auto const& b)
  requires(mat_rows<decltype(a)> == vec_size<decltype(b)>)
{
  return normalize(transform_vector(a, b));
}
constexpr kernel auto transform_normal(numNxM auto const& a, numN auto const& b)
  requires(mat_rows<decltype(a)> == vec_size<decltype(b)>)
{
  return normalize(transform_vector(transpose(inverse(a)), b));
}

// Transforms points, vectors and directions by frames.
constexpr kernel auto transform_point(
    numNxF auto const& a, numN auto const& b) {
  return a.m() * b + a.t();
}
constexpr kernel auto transform_vector(
    numNxF auto const& a, numN auto const& b) {
  return a.m() * b;
}
constexpr kernel auto transform_direction(
    numNxF auto const& a, numN auto const& b) {
  return normalize(transform_vector(a, b));
}
constexpr kernel auto transform_normal(
    numNxF auto const& a, numN auto const& b, bool non_rigid = false) {
  if (non_rigid) {
    return transform_normal(rotation(a), b);
  } else {
    return normalize(transform_vector(a, b));
  }
}

// Transforms points, vectors and directions by frames.
template <typename T, index_t N>
constexpr kernel auto transform_point_inverse(
    numNxF auto const& a, numN auto const& b) {
  return (b - a.t()) * a.m();
}
template <typename T, index_t N>
constexpr kernel auto transform_vector_inverse(
    numNxF auto const& a, numN auto const& b) {
  return b * a.m();
}
template <typename T, index_t N>
constexpr kernel auto transform_direction_inverse(
    numNxF auto const& a, numN auto const& b) {
  return normalize(transform_vector_inverse(a, b));
}

// Translation, scaling and rotations transforms.
constexpr kernel auto translation_frame(num3 auto const& a) {
  return frame{vec{1, 0, 0}, vec{0, 1, 0}, vec{0, 0, 1}, a};
}
constexpr kernel auto scaling_frame(num3 auto const& a) {
  return frame{
      vec{a.x(), 0, 0}, vec{0, a.y(), 0}, vec{0, 0, a.z()}, vec{0, 0, 0}};
}
constexpr kernel auto scaling_frame(num auto a) {
  return scaling_frame(vec{a, a, a});
}
constexpr kernel auto rotation_frame(num3 auto const& axis, num auto angle) {
  auto s = sin(angle), c = cos(angle);
  auto vv = normalize(axis);
  return frame{
      vec{c + (1 - c) * vv.x() * vv.x(), (1 - c) * vv.x() * vv.y() + s * vv.z(),
          (1 - c) * vv.x() * vv.z() - s * vv.y()},
      vec{(1 - c) * vv.x() * vv.y() - s * vv.z(), c + (1 - c) * vv.y() * vv.y(),
          (1 - c) * vv.y() * vv.z() + s * vv.x()},
      vec{(1 - c) * vv.x() * vv.z() + s * vv.y(),
          (1 - c) * vv.y() * vv.z() - s * vv.x(),
          c + (1 - c) * vv.z() * vv.z()},
      vec{0, 0, 0},
  };
}
constexpr kernel auto rotation_frame(num4 auto const& quat) {
  auto v = quat;
  return frame{
      vec{v.w * v.w + v.x() * v.x() - v.y() * v.y() - v.z() * v.z(),
          (v.x() * v.y() + v.z() * v.w) * 2, (v.z() * v.x() - v.y() * v.w) * 2},
      vec{(v.x() * v.y() - v.z() * v.w) * 2,
          v.w * v.w - v.x() * v.x() + v.y() * v.y() - v.z() * v.z(),
          (v.y() * v.z() + v.x() * v.w) * 2},
      vec{(v.z() * v.x() + v.y() * v.w) * 2, (v.y() * v.z() - v.x() * v.w) * 2,
          v.w * v.w - v.x() * v.x() - v.y() * v.y() + v.z() * v.z()},
      vec{0, 0, 0},
  };
}
constexpr kernel auto rotation_frame(numQ auto const& quat) {
  auto v = quat;
  return frame{
      vec{v.w * v.w + v.x() * v.x() - v.y() * v.y() - v.z() * v.z(),
          (v.x() * v.y() + v.z() * v.w) * 2, (v.z() * v.x() - v.y() * v.w) * 2},
      vec{(v.x() * v.y() - v.z() * v.w) * 2,
          v.w * v.w - v.x() * v.x() + v.y() * v.y() - v.z() * v.z(),
          (v.y() * v.z() + v.x() * v.w) * 2},
      vec{(v.z() * v.x() + v.y() * v.w) * 2, (v.y() * v.z() - v.x() * v.w) * 2,
          v.w * v.w - v.x() * v.x() - v.y() * v.y() + v.z() * v.z()},
      vec{0, 0, 0},
  };
}
constexpr kernel auto rotation_frame(numNxM auto const& rot) {
  return frame{rot.x(), rot.y(), rot.z(), {0, 0, 0}};
}

// Lookat frame. Z-axis can be inverted with inv_xz.
constexpr kernel auto lookat_frame(num3 auto const& eye,
    num3 auto const& center, num3 auto const& up, bool inv_xz = false) {
  auto w = normalize(eye - center);
  auto u = normalize(cross(up, w));
  auto v = normalize(cross(w, u));
  if (inv_xz) {
    w = -w;
    u = -u;
  }
  return frame{u, v, w, eye};
}

// OpenGL frustum, ortho and perspecgive matrices.
constexpr kernel auto frustum_mat(
    num auto l, num auto r, num auto b, num auto t, num auto n, num auto f) {
  return mat{
      vec{2 * n / (r - l), 0, 0, 0},
      vec{0, 2 * n / (t - b), 0, 0},
      vec{(r + l) / (r - l), (t + b) / (t - b), -(f + n) / (f - n), -1},
      vec{0, 0, -2 * f * n / (f - n), 0},
  };
}
constexpr kernel auto ortho_mat(
    num auto l, num auto r, num auto b, num auto t, num auto n, num auto f) {
  return mat{
      vec{2 / (r - l), 0, 0, 0},
      vec{0, 2 / (t - b), 0, 0},
      vec{0, 0, -2 / (f - n), 0},
      vec{-(r + l) / (r - l), -(t + b) / (t - b), -(f + n) / (f - n), 1},
  };
}
constexpr kernel auto ortho2d_mat(
    num auto left, num auto right, num auto bottom, num auto top) {
  return ortho_mat(left, right, bottom, top, -1, 1);
}
constexpr kernel auto ortho_mat(
    num auto xmag, num auto ymag, num auto near, num auto far) {
  return mat{
      vec{1 / xmag, 0, 0, 0},
      vec{0, 1 / ymag, 0, 0},
      vec{0, 0, 2 / (near - far), 0},
      vec{0, 0, (far + near) / (near - far), 1},
  };
}
constexpr kernel auto perspective_mat(
    num auto fovy, num auto aspect, num auto near, num auto far) {
  auto tg = tan(fovy / 2);
  return mat{
      vec{1 / (aspect * tg), 0, 0, 0},
      vec{0, 1 / tg, 0, 0},
      vec{0, 0, (far + near) / (near - far), -1},
      vec{0, 0, 2 * far * near / (near - far), 0},
  };
}
constexpr kernel auto perspective_mat(
    num auto fovy, num auto aspect, num auto near) {
  auto tg = tan(fovy / 2);
  return mat{
      vec{1 / (aspect * tg), 0, 0, 0},
      vec{0, 1 / tg, 0, 0},
      vec{0, 0, -1, -1},
      vec{0, 0, 2 * near, 0},
  };
}

// Rotation conversions. Returns axis and angle.
constexpr kernel auto rotation_axisangle(num4 auto const& quat) {
  return pair{normalize(vec{quat.x(), quat.y(), quat.z()}), 2 * acos(quat.w())};
}
constexpr kernel auto rotation_quat(num3 auto const& axis, num auto angle) {
  auto len = length(axis);
  if (len == 0) return vec{0, 0, 0, 1};
  return vec{sin(angle / 2) * axis.x / len, sin(angle / 2) * axis.y / len,
      sin(angle / 2) * axis.z / len, cos(angle / 2)};
}
constexpr kernel auto rotation_quat(num4 auto const& axisangle) {
  return rotation_quat(vec{axisangle.x, axisangle.y, axisangle.z}, axisangle.w);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// USER INTERFACE UTILITIES
// -----------------------------------------------------------------------------
namespace yocto {

// Computes the image uv coordinates corresponding to the view parameters.
// Returns negative coordinates if out of the image.
constexpr kernel auto image_coords(num2 auto const& mouse_pos,
    num2 auto const& center, num auto scale, inte2 auto const& size,
    bool clamped = true) {
  using I = vec_etype<decltype(size)>;
  using T = vec_etype<decltype(center)>;
  auto xy = (mouse_pos - center) / scale;
  auto ij = (vec<I, 2>)round(xy + size / (T)2);
  return clamped ? clamp(ij, 0, size) : ij;
}

// Center image and autofit. Returns center and scale.
constexpr kernel auto camera_imview(num2 auto const& center, num auto scale,
    inte2 auto const& imsize, inte2 auto const& winsize, bool zoom_to_fit) {
  using T = vec_etype<decltype(center)>;
  if (zoom_to_fit) {
    return pair{(vec<T, 2>)winsize / 2, min(winsize / (vec<T, 2>)imsize)};
  } else {
    return pair{select(component_greater_equal(winsize, imsize * scale),
                    (vec<T, 2>)winsize / 2, center),
        scale};
  }
}

// Turntable for UI navigation. Returns from and to.
constexpr kernel auto camera_turntable(num3 auto const& from_,
    num3 auto const& to_, num3 auto const& up, num2 auto const& rotate,
    num auto dolly, num2 auto const& pan) {
  using T = vec_etype<decltype(from_)>;

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
  return pair{from, to};
}

// Turntable for UI navigation. Returns frame and focus.
constexpr kernel auto camera_turntable(num3xF auto const& frame_,
    num auto focus, num2 auto const& rotate, num auto dolly,
    num2 auto const& pan) {
  using T = frame_etype<decltype(frame_)>;

  // copy values
  auto frame = frame_;

  // rotate if necessary
  if (rotate != vec<T, 2>{0, 0}) {
    auto phi   = atan2(frame.z().z(), frame.z().x()) + rotate.x();
    auto theta = acos(frame.z().y()) + rotate.y();
    theta      = clamp(theta, (T)0.001, (T)pi - (T)0.001);
    auto new_z = vec<T, 3>{
        sin(theta) * cos(phi), cos(theta), sin(theta) * sin(phi)};
    auto new_center = frame.o() - frame.z() * focus;
    auto new_o      = new_center + new_z * focus;
    frame           = lookat_frame(new_o, new_center, vec<T, 3>{0, 1, 0});
    focus           = length(new_o - new_center);
  }

  // pan if necessary
  if (dolly != 0) {
    auto c    = frame.o() - frame.z() * focus;
    focus     = max(focus * (1 + dolly), (T)0.001);
    frame.o() = c + frame.z() * focus;
  }

  // pan if necessary
  if (pan != vec<T, 2>{0, 0}) {
    frame.o() += frame.x() * pan.x() + frame.y() * pan.y();
  }

  // done
  return pair{frame, focus};
}

// FPS camera for UI navigation for a frame parametrization. Returns frame.
constexpr kernel auto camera_fpscam(num3xF auto const& frame,
    num3 auto const& transl, num2 auto const& rotate) {
  using T = frame_etype<decltype(frame)>;

  // https://gamedev.stackexchange.com/questions/30644/how-to-keep-my-quaternion-using-fps-camera-from-tilting-and-messing-up
  auto y = vec<T, 3>{0, 1, 0};
  auto z = orthonormalize(frame.z(), y);
  auto x = cross(y, z);

  auto rot = rotation_frame(vec<T, 3>{1, 0, 0}, rotate.y()) *
             yocto::frame<T, 3>{
                 frame.x(), frame.y(), frame.z(), vec<T, 3>{0, 0, 0}} *
             rotation_frame(vec<T, 3>{0, 1, 0}, rotate.x());
  auto pos = frame.o() + transl.x() * x + transl.y() * y + transl.z() * z;

  return yocto::frame{rot.x(), rot.y(), rot.z(), pos};
}

// Computes the image uv coordinates corresponding to the view parameters.
// Returns negative coordinates if out of the image.
constexpr kernel vec2i get_image_coords(num2 auto const& mouse_pos,
    num2 auto const& center, num auto scale, inte2 auto const& size) {
  using T = vec_etype<decltype(center)>;
  auto xy = (mouse_pos - center) / scale;
  return (vec2i)round(xy + size / (T)2);
}

// Center image and autofit.
constexpr kernel void update_imview(inte2 auto& center, num auto& scale,
    inte2 auto const& imsize, inte2 auto const& winsize, bool zoom_to_fit) {
  using T = vec_etype<decltype(center)>;
  if (zoom_to_fit) {
    scale  = min((vec<T, 2>)winsize / imsize);
    center = (vec<T, 2>)winsize / 2;
  } else {
    if (winsize[0] >= imsize[0] * scale) center[0] = (T)winsize[0] / 2;
    if (winsize[1] >= imsize[1] * scale) center[1] = (T)winsize[1] / 2;
  }
}

// Turntable for UI navigation.
constexpr kernel void update_turntable(num3 auto& from, num3 auto& to,
    num3 auto& up, num2 auto const& rotate, num auto dolly,
    num2 auto const& pan) {
  using T = vec_etype<decltype(from)>;

  // rotate if necessary
  if (rotate != vec<T, 2>{0, 0}) {
    auto z     = normalize(to - from);
    auto lz    = length(to - from);
    auto phi   = atan2(z.z(), z.x()) + rotate.x();
    auto theta = acos(z.y()) + rotate.y();
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
    auto t = pan.x() * x + pan.y() * y;
    from += t;
    to += t;
  }
}

// Turntable for UI navigation.
constexpr kernel void update_turntable(num3xF auto& frame, num auto& focus,
    num2 auto const& rotate, num auto dolly, num2 auto const& pan) {
  using T = frame_etype<decltype(frame)>;

  // rotate if necessary
  if (rotate != vec<T, 2>{0, 0}) {
    auto phi   = atan2(frame.z.z(), frame.z.x()) + rotate.x();
    auto theta = acos(frame.z.y()) + rotate.y();
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
    frame.o = c + frame.z() * focus;
  }

  // pan if necessary
  if (pan != vec<T, 2>{0, 0}) {
    frame.o() += frame.x() * pan.x() + frame.y() * pan.y();
  }
}

// FPS camera for UI navigation for a frame parametrization.
constexpr kernel void update_fpscam(
    num3xF auto& frame, num3 auto const& transl, num2 auto const& rotate) {
  // https://gamedev.stackexchange.com/questions/30644/how-to-keep-my-quaternion-using-fps-camera-from-tilting-and-messing-up
  using T = frame_etype<decltype(frame)>;

  auto y = vec<T, 3>{0, 1, 0};
  auto z = orthonormalize(frame.z, y);
  auto x = cross(y, z);

  auto rot = rotation_frame(vec<T, 3>{1, 0, 0}, rotate.y) *
             yocto::frame<T, 3>{frame.x, frame.y, frame.z, vec<T, 3>{0, 0, 0}} *
             rotation_frame(vec<T, 3>{0, 1, 0}, rotate.x());
  auto pos = frame.o + transl.x() * x + transl.y() * y + transl.z() * z;

  frame = {rot, pos};
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

// Python range in 3d. Construct an object that iterates over a 3d integer
// sequence.
template <typename T>
constexpr kernel auto range(vec<T, 3> max);

// Python enumerate
template <typename Sequence, typename T = index_t>
constexpr kernel auto enumerate(const Sequence& sequence, T start = 0);
template <typename Sequence, typename T = index_t>
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

// Python range in 2d. Construct an object that iterates over a 2d integer
// sequence.
template <typename T>
constexpr kernel auto range(vec<T, 3> max) {
  struct range_sentinel {};
  struct range_iterator {
    vec<T, 3>             index, end;
    constexpr kernel void operator++() {
      ++index[0];
      if (index[0] >= end[0]) {
        index[0] = 0;
        index[1]++;
      }
      if (index[1] >= end[1]) {
        index[1] = 0;
        index[2]++;
      }
    }
    constexpr kernel bool operator!=(const range_sentinel&) const {
      return index[2] != end[2];
    }
    constexpr kernel vec<T, 3> operator*() const { return index; }
  };
  struct range_sequence {
    vec<T, 3>                       end_ = {0, 0, 0};
    constexpr kernel range_iterator begin() const { return {{0, 0, 0}, end_}; }
    constexpr kernel range_sentinel end() const { return {}; }
  };
  return range_sequence{max};
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
