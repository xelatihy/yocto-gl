//
// # Yocto/Sampling: Sampling routines
//
// Yocto/Sampling provides many functions to generate points and directions
// useful in path tracing and procedural generation. We also include a random
// number generator suitable for ray tracing.
// This library includes a stand-alone implementaton of the PCG32 random number
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
#ifdef __CUDACC__
#define inline inline __device__ __forceinline__
#endif

// -----------------------------------------------------------------------------
// RANDOM NUMBER GENERATION
// -----------------------------------------------------------------------------
namespace yocto {

// PCG random numbers from http://www.pcg-random.org/
struct rng_state {
  uint64_t state = 0x853c49e6748fea9bULL;
  uint64_t inc   = 0xda3e39cb94b95bdbULL;
};

// Init a random number generator with a state state from the sequence seq.
inline rng_state make_rng(uint64_t seed, uint64_t seq = 1);

// Init a sequence of random number generators.
inline vector<rng_state> make_rngs(
    int num, uint64_t seed1 = 961748941ull, uint64_t seed2 = 1301081ull);

// Next random numbers: floats in [0,1), ints in [0,n).
inline int   rand1i(rng_state& rng, int n);
inline float rand1f(rng_state& rng);
inline vec2f rand2f(rng_state& rng);
inline vec3f rand3f(rng_state& rng);

// Shuffles a sequence of elements
template <typename T>
inline void shuffle(vector<T>& vals, rng_state& rng);

}  // namespace yocto

// -----------------------------------------------------------------------------
// MONETACARLO SAMPLING FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Sample an hemispherical direction with uniform distribution.
inline vec3f sample_hemisphere(vec2f ruv);
inline float sample_hemisphere_pdf(vec3f direction);

// Sample an hemispherical direction with uniform distribution.
inline vec3f sample_hemisphere(vec3f normal, vec2f ruv);
inline float sample_hemisphere_pdf(vec3f normal, vec3f direction);

// Sample a spherical direction with uniform distribution.
inline vec3f sample_sphere(vec2f ruv);
inline float sample_sphere_pdf(vec3f w);

// Sample an hemispherical direction with cosine distribution.
inline vec3f sample_hemisphere_cos(vec2f ruv);
inline float sample_hemisphere_cos_pdf(vec3f direction);

// Sample an hemispherical direction with cosine distribution.
inline vec3f sample_hemisphere_cos(vec3f normal, vec2f ruv);
inline float sample_hemisphere_cos_pdf(vec3f normal, vec3f direction);

// Sample an hemispherical direction with cosine power distribution.
inline vec3f sample_hemisphere_cospower(float exponent, vec2f ruv);
inline float sample_hemisphere_cospower_pdf(float exponent, vec3f direction);

// Sample an hemispherical direction with cosine power distribution.
inline vec3f sample_hemisphere_cospower(
    float exponent, vec3f normal, vec2f ruv);
inline float sample_hemisphere_cospower_pdf(
    float exponent, vec3f normal, vec3f direction);

// Sample a point uniformly on a disk.
inline vec2f sample_disk(vec2f ruv);
inline float sample_disk_pdf(vec2f point);

// Sample a point uniformly on a cylinder, without caps.
inline vec3f sample_cylinder(vec2f ruv);
inline float sample_cylinder_pdf(vec3f point);

// Sample a point uniformly on a triangle returning the baricentric coordinates.
inline vec2f sample_triangle(vec2f ruv);

// Sample a point uniformly on a triangle.
inline vec3f sample_triangle(vec3f p0, vec3f p1, vec3f p2, vec2f ruv);
// Pdf for uniform triangle sampling, i.e. triangle area.
inline float sample_triangle_pdf(vec3f p0, vec3f p1, vec3f p2);

// Sample an index with uniform distribution.
inline int   sample_uniform(int size, float r);
inline float sample_uniform_pdf(int size);

// Sample an index with uniform distribution.
inline float sample_uniform(const vector<float>& elements, float r);
inline float sample_uniform_pdf(const vector<float>& elements);

// Sample a discrete distribution represented by its cdf.
inline int sample_discrete(const vector<float>& cdf, float r);
// Pdf for uniform discrete distribution sampling.
inline float sample_discrete_pdf(const vector<float>& cdf, int idx);

}  // namespace yocto

// -----------------------------------------------------------------------------
// SHAPE SAMPLING
// -----------------------------------------------------------------------------
namespace yocto {

// Pick a point in a point set uniformly.
inline int           sample_points(int npoints, float re);
inline int           sample_points(const vector<float>& cdf, float re);
inline vector<float> sample_points_cdf(int npoints);

// Pick a point on lines uniformly.
inline pair<int, float> sample_lines(
    const vector<float>& cdf, float re, float ru);
inline vector<float> sample_lines_cdf(
    const vector<vec2i>& lines, const vector<vec3f>& positions);

// Pick a point on a triangle mesh uniformly.
inline pair<int, vec2f> sample_triangles(
    const vector<float>& cdf, float re, vec2f ruv);
inline vector<float> sample_triangles_cdf(
    const vector<vec3i>& triangles, const vector<vec3f>& positions);

// Pick a point on a quad mesh uniformly.
inline pair<int, vec2f> sample_quads(
    const vector<float>& cdf, float re, vec2f ruv);
inline vector<float> sample_quads_cdf(
    const vector<vec4i>& quads, const vector<vec3f>& positions);

}  // namespace yocto

// -----------------------------------------------------------------------------
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR RANDOM NUMBER GENERATION
// -----------------------------------------------------------------------------
namespace yocto {

// PCG random numbers from http://www.pcg-random.org/

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

// Init a sequence of random number generators.
inline vector<rng_state> make_rngs(int num, uint64_t seed1, uint64_t seed2) {
  auto rng_ = make_rng(seed2);
  auto rngs = vector<rng_state>(num);
  for (auto& rng : rngs) {
    rng = make_rng(seed1, rand1i(rng_, 1 << 31) / 2 + 1);
  }
  return rngs;
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

// Shuffles a sequence of elements
template <typename T>
inline void shuffle(vector<T>& vals, rng_state& rng) {
  // https://en.wikipedia.org/wiki/Fisher–Yates_shuffle
  for (auto i = (int)vals.size() - 1; i > 0; i--) {
    auto j = rand1i(rng, i + 1);
    std::swap(vals[j], vals[i]);
  }
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF MONTECARLO SAMPLING FUNCTIONS
// -----------------------------------------------------------------------------
namespace yocto {

// Sample an hemispherical direction with uniform distribution.
inline vec3f sample_hemisphere(vec2f ruv) {
  auto z   = ruv.y;
  auto r   = sqrt(clamp(1 - z * z, 0.0f, 1.0f));
  auto phi = 2 * pif * ruv.x;
  return {r * cos(phi), r * sin(phi), z};
}
inline float sample_hemisphere_pdf(vec3f direction) {
  return (direction.z <= 0) ? 0 : 1 / (2 * pif);
}

// Sample an hemispherical direction with uniform distribution.
inline vec3f sample_hemisphere(vec3f normal, vec2f ruv) {
  auto z               = ruv.y;
  auto r               = sqrt(clamp(1 - z * z, 0.0f, 1.0f));
  auto phi             = 2 * pif * ruv.x;
  auto local_direction = vec3f{r * cos(phi), r * sin(phi), z};
  return transform_direction(basis_fromz(normal), local_direction);
}
inline float sample_hemisphere_pdf(vec3f normal, vec3f direction) {
  return (dot(normal, direction) <= 0) ? 0 : 1 / (2 * pif);
}

// Sample a spherical direction with uniform distribution.
inline vec3f sample_sphere(vec2f ruv) {
  auto z   = 2 * ruv.y - 1;
  auto r   = sqrt(clamp(1 - z * z, 0.0f, 1.0f));
  auto phi = 2 * pif * ruv.x;
  return {r * cos(phi), r * sin(phi), z};
}
inline float sample_sphere_pdf(vec3f w) { return 1 / (4 * pif); }

// Sample an hemispherical direction with cosine distribution.
inline vec3f sample_hemisphere_cos(vec2f ruv) {
  auto z   = sqrt(ruv.y);
  auto r   = sqrt(1 - z * z);
  auto phi = 2 * pif * ruv.x;
  return {r * cos(phi), r * sin(phi), z};
}
inline float sample_hemisphere_cos_pdf(vec3f direction) {
  return (direction.z <= 0) ? 0 : direction.z / pif;
}

// Sample an hemispherical direction with cosine distribution.
inline vec3f sample_hemisphere_cos(vec3f normal, vec2f ruv) {
  auto z               = sqrt(ruv.y);
  auto r               = sqrt(1 - z * z);
  auto phi             = 2 * pif * ruv.x;
  auto local_direction = vec3f{r * cos(phi), r * sin(phi), z};
  return transform_direction(basis_fromz(normal), local_direction);
}
inline float sample_hemisphere_cos_pdf(vec3f normal, vec3f direction) {
  auto cosw = dot(normal, direction);
  return (cosw <= 0) ? 0 : cosw / pif;
}

// Sample an hemispherical direction with cosine power distribution.
inline vec3f sample_hemisphere_cospower(float exponent, vec2f ruv) {
  auto z   = pow(ruv.y, 1 / (exponent + 1));
  auto r   = sqrt(1 - z * z);
  auto phi = 2 * pif * ruv.x;
  return {r * cos(phi), r * sin(phi), z};
}
inline float sample_hemisphere_cospower_pdf(float exponent, vec3f direction) {
  return (direction.z <= 0)
             ? 0
             : pow(direction.z, exponent) * (exponent + 1) / (2 * pif);
}

// Sample an hemispherical direction with cosine power distribution.
inline vec3f sample_hemisphere_cospower(
    float exponent, vec3f normal, vec2f ruv) {
  auto z               = pow(ruv.y, 1 / (exponent + 1));
  auto r               = sqrt(1 - z * z);
  auto phi             = 2 * pif * ruv.x;
  auto local_direction = vec3f{r * cos(phi), r * sin(phi), z};
  return transform_direction(basis_fromz(normal), local_direction);
}
inline float sample_hemisphere_cospower_pdf(
    float exponent, vec3f normal, vec3f direction) {
  auto cosw = dot(normal, direction);
  return (cosw <= 0) ? 0 : pow(cosw, exponent) * (exponent + 1) / (2 * pif);
}

// Sample a point uniformly on a disk.
inline vec2f sample_disk(vec2f ruv) {
  auto r   = sqrt(ruv.y);
  auto phi = 2 * pif * ruv.x;
  return {cos(phi) * r, sin(phi) * r};
}
inline float sample_disk_pdf() { return 1 / pif; }

// Sample a point uniformly on a cylinder, without caps.
inline vec3f sample_cylinder(vec2f ruv) {
  auto phi = 2 * pif * ruv.x;
  return {sin(phi), cos(phi), ruv.y * 2 - 1};
}
inline float sample_cylinder_pdf(vec3f point) { return 1 / pif; }

// Sample a point uniformly on a triangle returning the baricentric coordinates.
inline vec2f sample_triangle(vec2f ruv) {
  return {1 - sqrt(ruv.x), ruv.y * sqrt(ruv.x)};
}

// Sample a point uniformly on a triangle.
inline vec3f sample_triangle(vec3f p0, vec3f p1, vec3f p2, vec2f ruv) {
  auto uv = sample_triangle(ruv);
  return p0 * (1 - uv.x - uv.y) + p1 * uv.x + p2 * uv.y;
}
// Pdf for uniform triangle sampling, i.e. triangle area.
inline float sample_triangle_pdf(vec3f p0, vec3f p1, vec3f p2) {
  return 2 / length(cross(p1 - p0, p2 - p0));
}

// Sample an index with uniform distribution.
inline int sample_uniform(int size, float r) {
  return clamp((int)(r * size), 0, size - 1);
}
inline float sample_uniform_pdf(int size) { return (float)1 / (float)size; }

// Sample an index with uniform distribution.
inline float sample_uniform(const vector<float>& elements, float r) {
  if (elements.empty()) return {};
  auto size = (int)elements.size();
  return elements[clamp((int)(r * size), 0, size - 1)];
}
inline float sample_uniform_pdf(const vector<float>& elements) {
  if (elements.empty()) return 0;
  return 1.0f / (int)elements.size();
}

// Sample a discrete distribution represented by its cdf.
inline int sample_discrete(const vector<float>& cdf, float r) {
  r        = clamp(r * cdf.back(), (float)0, cdf.back() - (float)0.00001);
  auto idx = (int)(std::upper_bound(cdf.data(), cdf.data() + cdf.size(), r) -
                   cdf.data());
  return clamp(idx, 0, (int)cdf.size() - 1);
}
// Pdf for uniform discrete distribution sampling.
inline float sample_discrete_pdf(const vector<float>& cdf, int idx) {
  if (idx == 0) return cdf[0];
  return cdf[idx] - cdf[idx - 1];
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHAPE SAMPLING
// -----------------------------------------------------------------------------
namespace yocto {

// Pick a point in a point set uniformly.
inline int sample_points(int npoints, float re) {
  return sample_uniform(npoints, re);
}
inline int sample_points(const vector<float>& cdf, float re) {
  return sample_discrete(cdf, re);
}
inline vector<float> sample_points_cdf(int npoints) {
  auto cdf = vector<float>(npoints);
  for (auto i : range(cdf.size())) cdf[i] = 1 + (i != 0 ? cdf[i - 1] : 0);
  return cdf;
}

// Pick a point on lines uniformly.
inline pair<int, float> sample_lines(
    const vector<float>& cdf, float re, float ru) {
  return {sample_discrete(cdf, re), ru};
}
inline vector<float> sample_lines_cdf(
    const vector<vec2i>& lines, const vector<vec3f>& positions) {
  auto cdf = vector<float>(lines.size());
  for (auto i : range(cdf.size())) {
    auto& l = lines[i];
    auto  w = line_length(positions[l.x], positions[l.y]);
    cdf[i]  = w + (i != 0 ? cdf[i - 1] : 0);
  }
  return cdf;
}

// Pick a point on a triangle mesh uniformly.
inline pair<int, vec2f> sample_triangles(
    const vector<float>& cdf, float re, vec2f ruv) {
  return {sample_discrete(cdf, re), sample_triangle(ruv)};
}
inline vector<float> sample_triangles_cdf(
    const vector<vec3i>& triangles, const vector<vec3f>& positions) {
  auto cdf = vector<float>(triangles.size());
  for (auto i : range(cdf.size())) {
    auto& t = triangles[i];
    auto  w = triangle_area(positions[t.x], positions[t.y], positions[t.z]);
    cdf[i]  = w + (i != 0 ? cdf[i - 1] : 0);
  }
  return cdf;
}

// Pick a point on a quad mesh uniformly.
inline pair<int, vec2f> sample_quads(
    const vector<float>& cdf, float re, vec2f ruv) {
  return {sample_discrete(cdf, re), ruv};
}
inline pair<int, vec2f> sample_quads(
    const vector<vec4i>& quads, const vector<float>& cdf, float re, vec2f ruv) {
  auto element = sample_discrete(cdf, re);
  if (quads[element].z == quads[element].w) {
    return {element, sample_triangle(ruv)};
  } else {
    return {element, ruv};
  }
}
inline vector<float> sample_quads_cdf(
    const vector<vec4i>& quads, const vector<vec3f>& positions) {
  auto cdf = vector<float>(quads.size());
  for (auto i : range(cdf.size())) {
    auto& q = quads[i];
    auto  w = quad_area(
        positions[q.x], positions[q.y], positions[q.z], positions[q.w]);
    cdf[i] = w + (i ? cdf[i - 1] : 0);
  }
  return cdf;
}

}  // namespace yocto

// -----------------------------------------------------------------------------
// CUDA SUPPORT
// -----------------------------------------------------------------------------
#ifdef __CUDACC__
#undef inline
#endif

#endif
