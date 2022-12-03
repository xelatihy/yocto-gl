//
// # Yocto/Sampling: Sampling routines
//
// Yocto/Sampling provides many functions to generate points and directions
// useful in path tracing and procedural generation. We also include a random
// number generator suitable for ray tracing.
// This library includes a stand-alone implementation of the PCG32 random number
// generator by M.E. O'Neill.
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
//
// LICENSE OF INCLUDED SOFTWARE for Pcg random number generator
//
// This code also includes a small exerpt from http://www.pcg-random.org/
// licensed as follows
// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)
//

#ifndef _YOCTO_SAMPLING_H_
#define _YOCTO_SAMPLING_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include <algorithm>  // std::upper_bound
#include <array>
#include <utility>
#include <vector>

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// USING DIRECTIVES
// -----------------------------------------------------------------------------
namespace yocto {

// using directives
using std::array;
using std::vector;

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
// RANDOM NUMBER GENERATION
// -----------------------------------------------------------------------------
namespace yocto {

// PCG random numbers from http://www.pcg-random.org/
struct rng_state {
  uint64_t state = 0x853c49e6748fea9bULL;
  uint64_t inc   = 0xda3e39cb94b95bdbULL;

  constexpr kernel rng_state()
      : state{0x853c49e6748fea9bULL}, inc{0xda3e39cb94b95bdbULL} {}
  constexpr kernel rng_state(uint64_t state_, uint64_t inc_)
      : state{state_}, inc{inc_} {}
};

// Next random number, used internally only.
constexpr kernel uint32_t _advance_rng(rng_state& rng) {
  uint64_t oldstate = rng.state;
  rng.state         = oldstate * 6364136223846793005ULL + rng.inc;
  auto xorshifted   = (uint32_t)(((oldstate >> 18u) ^ oldstate) >> 27u);
  auto rot          = (uint32_t)(oldstate >> 59u);
  // return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
  return (xorshifted >> rot) | (xorshifted << ((~rot + 1u) & 31));
}

// Init a random number generator with a state state from the sequence seq.
constexpr kernel rng_state make_rng(uint64_t seed, uint64_t seq = 1) {
  auto rng  = rng_state{};
  rng.state = 0U;
  rng.inc   = (seq << 1u) | 1u;
  _advance_rng(rng);
  rng.state += seed;
  _advance_rng(rng);
  return rng;
}

// Next random numbers: floats in [0,1), ints in [0,n).
constexpr kernel int rand1i(rng_state& rng, int n) {
  return _advance_rng(rng) % n;
}
constexpr kernel float rand1f(rng_state& rng) {
  // union {
  //   uint32_t u;
  //   float    f;
  // } x;
  // x.u = (_advance_rng(rng) >> 9) | 0x3f800000u;
  // return x.f - 1.0f;
  // alternate implementation
  constexpr auto scale = (float)(1.0 / std::numeric_limits<uint32_t>::max());
  return _advance_rng(rng) * scale;
}
constexpr kernel vec2f rand2f(rng_state& rng) {
  // force order of evaluation by using separate assignments.
  auto x = rand1f(rng);
  auto y = rand1f(rng);
  return {x, y};
}
constexpr kernel vec3f rand3f(rng_state& rng) {
  // force order of evaluation by using separate assignments.
  auto x = rand1f(rng);
  auto y = rand1f(rng);
  auto z = rand1f(rng);
  return {x, y, z};
}

// Shuffles a sequence of elements
template <typename T>
constexpr kernel void shuffle(span<T> vals, rng_state& rng) {
  // https://en.wikipedia.org/wiki/Fisher–Yates_shuffle
  for (auto i = (int)vals.size() - 1; i > 0; i--) {
    auto j = rand1i(rng, i + 1);
    std::swap(vals[j], vals[i]);
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// MONTECARLO SAMPLING FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Sample an hemispherical direction with uniform distribution.
template <typename T>
constexpr kernel vec<T, 3> sample_hemisphere(const vec<T, 2>& ruv) {
  auto z   = ruv.y;
  auto r   = sqrt(clamp(1 - z * z, 0.0f, 1.0f));
  auto phi = 2 * pif * ruv.x;
  return {r * cos(phi), r * sin(phi), z};
}
template <typename T>
constexpr kernel float sample_hemisphere_pdf(const vec<T, 3>& direction) {
  return (direction.z <= 0) ? 0 : 1 / (2 * pif);
}

// Sample an hemispherical direction with uniform distribution.
template <typename T>
constexpr kernel vec<T, 3> sample_hemisphere(
    const vec<T, 3>& normal, const vec<T, 2>& ruv) {
  auto z               = ruv.y;
  auto r               = sqrt(clamp(1 - z * z, 0.0f, 1.0f));
  auto phi             = 2 * pif * ruv.x;
  auto local_direction = vec<T, 3>{r * cos(phi), r * sin(phi), z};
  return transform_direction(basis_fromz(normal), local_direction);
}
template <typename T>
constexpr kernel float sample_hemisphere_pdf(
    const vec<T, 3>& normal, const vec<T, 3>& direction) {
  return (dot(normal, direction) <= 0) ? 0 : 1 / (2 * pif);
}

// Sample a spherical direction with uniform distribution.
template <typename T>
constexpr kernel vec<T, 3> sample_sphere(const vec<T, 2>& ruv) {
  auto z   = 2 * ruv.y - 1;
  auto r   = sqrt(clamp(1 - z * z, 0.0f, 1.0f));
  auto phi = 2 * pif * ruv.x;
  return {r * cos(phi), r * sin(phi), z};
}
template <typename T>
constexpr kernel float sample_sphere_pdf(const vec<T, 3>& w) {
  return 1 / (4 * pif);
}

// Sample an hemispherical direction with cosine distribution.
template <typename T>
constexpr kernel vec<T, 3> sample_hemisphere_cos(const vec<T, 2>& ruv) {
  auto z   = sqrt(ruv.y);
  auto r   = sqrt(1 - z * z);
  auto phi = 2 * pif * ruv.x;
  return {r * cos(phi), r * sin(phi), z};
}
template <typename T>
constexpr kernel float sample_hemisphere_cos_pdf(const vec<T, 3>& direction) {
  return (direction.z <= 0) ? 0 : direction.z / pif;
}

// Sample an hemispherical direction with cosine distribution.
template <typename T>
constexpr kernel vec<T, 3> sample_hemisphere_cos(
    const vec<T, 3>& normal, const vec<T, 2>& ruv) {
  auto z               = sqrt(ruv.y);
  auto r               = sqrt(1 - z * z);
  auto phi             = 2 * pif * ruv.x;
  auto local_direction = vec<T, 3>{r * cos(phi), r * sin(phi), z};
  return transform_direction(basis_fromz(normal), local_direction);
}
template <typename T>
constexpr kernel float sample_hemisphere_cos_pdf(
    const vec<T, 3>& normal, const vec<T, 3>& direction) {
  auto cosw = dot(normal, direction);
  return (cosw <= 0) ? 0 : cosw / pif;
}

// Sample an hemispherical direction with cosine power distribution.
template <typename T>
constexpr kernel vec<T, 3> sample_hemisphere_cospower(
    float exponent, const vec<T, 2>& ruv) {
  auto z   = pow(ruv.y, 1 / (exponent + 1));
  auto r   = sqrt(1 - z * z);
  auto phi = 2 * pif * ruv.x;
  return {r * cos(phi), r * sin(phi), z};
}
template <typename T>
constexpr kernel float sample_hemisphere_cospower_pdf(
    float exponent, const vec<T, 3>& direction) {
  return (direction.z <= 0)
             ? 0
             : pow(direction.z, exponent) * (exponent + 1) / (2 * pif);
}

// Sample an hemispherical direction with cosine power distribution.
template <typename T>
constexpr kernel vec<T, 3> sample_hemisphere_cospower(
    float exponent, const vec<T, 3>& normal, const vec<T, 2>& ruv) {
  auto z               = pow(ruv.y, 1 / (exponent + 1));
  auto r               = sqrt(1 - z * z);
  auto phi             = 2 * pif * ruv.x;
  auto local_direction = vec<T, 3>{r * cos(phi), r * sin(phi), z};
  return transform_direction(basis_fromz(normal), local_direction);
}
template <typename T>
constexpr kernel float sample_hemisphere_cospower_pdf(
    float exponent, const vec<T, 3>& normal, const vec<T, 3>& direction) {
  auto cosw = dot(normal, direction);
  return (cosw <= 0) ? 0 : pow(cosw, exponent) * (exponent + 1) / (2 * pif);
}

// Sample a point uniformly on a disk.
template <typename T>
constexpr kernel vec<T, 2> sample_disk(const vec<T, 2>& ruv) {
  auto r   = sqrt(ruv.y);
  auto phi = 2 * pif * ruv.x;
  return {cos(phi) * r, sin(phi) * r};
}
template <typename T>
constexpr kernel float sample_disk_pdf() {
  return 1 / pif;
}

// Sample a point uniformly on a cylinder, without caps.
template <typename T>
constexpr kernel vec<T, 3> sample_cylinder(const vec<T, 2>& ruv) {
  auto phi = 2 * pif * ruv.x;
  return {sin(phi), cos(phi), ruv.y * 2 - 1};
}
template <typename T>
constexpr kernel float sample_cylinder_pdf(const vec<T, 3>& point) {
  return 1 / pif;
}

// Sample a point uniformly on a triangle returning the baricentric coordinates.
template <typename T>
constexpr kernel vec<T, 2> sample_triangle(const vec<T, 2>& ruv) {
  return {1 - sqrt(ruv.x), ruv.y * sqrt(ruv.x)};
}

// Sample a point uniformly on a triangle.
template <typename T>
constexpr kernel vec<T, 3> sample_triangle(const vec<T, 3>& p0,
    const vec<T, 3>& p1, const vec<T, 3>& p2, const vec<T, 2>& ruv) {
  auto uv = sample_triangle(ruv);
  return p0 * (1 - uv.x - uv.y) + p1 * uv.x + p2 * uv.y;
}
// Pdf for uniform triangle sampling, i.e. triangle area.
template <typename T>
constexpr kernel float sample_triangle_pdf(
    const vec<T, 3>& p0, const vec<T, 3>& p1, const vec<T, 3>& p2) {
  return 2 / length(cross(p1 - p0, p2 - p0));
}

// Sample an index with uniform distribution.
template <typename T, typename I>
constexpr kernel I sample_uniform(I size, T r) {
  return clamp((I)(r * size), 0, size - 1);
}
template <typename I, typename T = float>
constexpr kernel T sample_uniform_pdf(I size) {
  return (T)1 / (T)size;
}

// Sample an index with uniform distribution.
template <typename T, typename E>
constexpr kernel E sample_uniform(cspan<E> elements, T r) {
  if (elements.empty()) return E{};
  auto size = elements.size();
  return elements[clamp((size_t)(r * size), (size_t)0, size - 1)];
}
template <typename E, typename T = float>
constexpr kernel T sample_uniform_pdf(cspan<E> elements) {
  if (elements.empty()) return 0;
  return (T)1 / (T)elements.size();
}

// TODO: this should be constexpr, once we understand why
// Sample a discrete distribution represented by its cdf.
template <typename T, typename I = int>
inline kernel I sample_discrete(cspan<T> cdf, T r) {
  r        = clamp(r * cdf.back(), (T)0, cdf.back() - (T)0.00001);
  auto idx = (I)(std::upper_bound(cdf.data(), cdf.data() + cdf.size(), r) -
                 cdf.data());
  return clamp(idx, (I)0, (I)cdf.size() - 1);
}
// Pdf for uniform discrete distribution sampling.
template <typename T, typename I>
inline kernel T sample_discrete_pdf(cspan<T> cdf, I idx) {
  if (idx == 0) return cdf[0];
  return cdf[idx] - cdf[idx - 1];
}

// TODO: eventually remove these
// Sample a discrete distribution represented by its cdf.
template <typename T, typename I = int>
inline kernel int sample_discrete(const vector<T>& cdf, T r) {
  return sample_discrete(cspan<T>{cdf}, r);
}
// Pdf for uniform discrete distribution sampling.
template <typename T, typename I>
inline kernel T sample_discrete_pdf(const vector<T>& cdf, I idx) {
  return sample_discrete_pdf(cspan<T>{cdf}, idx);
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// CUDA SUPPORT
// -----------------------------------------------------------------------------
#ifdef kernel
#undef kernel
#endif

#endif
