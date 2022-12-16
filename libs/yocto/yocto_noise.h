//
// # Yocto/Noise: Noise functions
//
// Yocto/Noise provides a Perlin noise implementation.
// Yocto/Noise is implemented in `yocto_noise.h`.
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

#ifndef _YOCTO_NOISE_H_
#define _YOCTO_NOISE_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_math.h"

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
// HASH FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

inline kernel vec<uint32_t, 3> hash_pcg(const vec<uint32_t, 3>& v_) {
  // clang-format off
  auto v = v_ * 1664525u + 1013904223u;
  v.x += v.y * v.z; v.y += v.z * v.x; v.z += v.x * v.y;
  v = v ^ (v >> 16u);
  v.x += v.y * v.z; v.y += v.z * v.x; v.z += v.x * v.y;
  // clang-format on
  return v;
}
inline kernel vec<uint32_t, 4> hash_pcg(const vec<uint32_t, 4>& v_) {
  // clang-format off
  auto v = v_ * 1664525u + 1013904223u;
  v.x += v.y * v.w; v.y += v.z * v.x; v.z += v.x * v.y; v.w += v.y * v.z;
  v = v ^ (v >> 16u);
  v.x += v.y * v.w; v.y += v.z * v.x; v.z += v.x * v.y; v.w += v.y * v.z;
  // clang-format on
  return v;
}
inline kernel vec<uint32_t, 2> hash_pcg(const vec<uint32_t, 2>& v) {
  return xy(hash_pcg(vec<uint32_t, 3>{v, 0}));
}
template <typename T, size_t N>
inline kernel vec<T, N> hash_pcgr(const vec<uint32_t, N>& v) {
  return hash_pcg(v) * (T)(1.0 / std::numeric_limits<uint32_t>::max());
}
template <typename T, size_t N>
inline kernel vec<T, N> hash_pcgr(const vec<int, N>& v) {
  auto vu = *(vec<uint32_t, N>*)&v;
  return hash_pcg(v) * (T)(1.0 / std::numeric_limits<uint32_t>::max());
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// NOISE FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// the implementation follows ideas from "Building up Perlin Noise"
// http://eastfarthing.com/blog/2015-04-21-noise/

// Interpolants and surflets
template <typename T>
constexpr kernel T _cubic_kernel(T x_) {
  auto x = abs(x_);
  return x >= 1 ? 0 : 1 - (3 - 2 * x) * x * x;
}
template <typename T, size_t N>
constexpr kernel T _cubic_kernel(const vec<T, N>& x) {
  return prod(map(x, [](T x) { return _cubic_kernel(x); }));
}

// Random values
template <typename T, size_t N>
inline kernel vec<T, N> _random_gradient(const vec<int, N>& cell) {
  return normalize(2 * hash_pcgr<T>(cell) - 1);
}
template <typename T, size_t N>
inline kernel vec<T, N> _random_value(const vec<int, N>& cell) {
  return hash_pcgr<T>(cell);
}

// Gradient noise
template <typename T, size_t N>
inline T gradient_noise(const vec<T, N>& p) {
  auto result = (T)0;
  auto cell   = (vec<int, N>)p;
  for (auto offset : range(vec<int, N>(2))) {
    auto index = cell + offset;
    result += _cubic_kernel(p - index) *
              dot(p - index, _random_gradient<T>(index));
  }
  result = (result + 1) / 2;
  return result;
}
// Value noise
template <typename T, size_t N>
inline T value_noise(const vec<T, N>& p) {
  auto result = (T)0;
  auto cell   = (vec<int, N>)p;
  for (auto offset : range(vec<int, N>(2))) {
    auto index = cell + offset;
    result += _cubic_kernel(p - index) * _random_value<T>(index)[0];
  }
  return result;
}

// Fractal noise
template <typename T, size_t N>
inline T fractal_noise(const vec<T, N>& p, int octaves = 6) {
  auto sum = (T)0;
  for (auto octave : range(octaves)) {
    sum += uv_to_unit(gradient_noise(p * pow2(octave))) / pow2(octave);
  }
  return unit_to_uv(sum);
}
// Ridge noise
template <typename T, size_t N>
inline T ridge_noise(const vec<T, N>& p, int octaves = 6) {
  auto sum = (T)0;
  for (auto octave : range(octaves)) {
    sum += (sqr(1 - abs(uv_to_unit(gradient_noise(p * pow2(octave))))) / 2) /
           pow2(octave);
  }
  return sum;
}
// Turbulence noise
template <typename T, size_t N>
inline T turbulence_noise(const vec<T, N>& p, int octaves = 6) {
  auto sum = (T)0;
  for (auto octave : range(octaves)) {
    sum += abs(uv_to_unit(gradient_noise(p * pow2(octave)))) / pow2(octave);
  }
  return sum;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// CUDA SUPPORT
// -----------------------------------------------------------------------------
#ifdef kernel
#undef kernel
#endif

#endif
