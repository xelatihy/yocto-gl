//
// # Yocto/Noise: Noise functions
//
// Yocto/Noise provides a Perlin noise implementation.
// _This library should to be considered a placeholder since it will grow later
// as a collection of noise functions_ used in procedural modeling.
// For now, the implementation used is the one
// found in the [stb libraries](https://github.com/nothings/stb),
// that are released in the public domain.
// Yocto/Noise is implemented in `yocto_noise.h`.
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

#ifndef _YOCTO_NOISE_H_
#define _YOCTO_NOISE_H_

// -----------------------------------------------------------------------------
// INCLUDES
// -----------------------------------------------------------------------------

#include "yocto_math.h"

// -----------------------------------------------------------------------------
// PERLIN NOISE FUNCTION
// -----------------------------------------------------------------------------
namespace yocto {

// Compute the revised Perlin noise function. Wrap provides a wrapping noise
// but must be power of two (wraps at 256 anyway). For octave based noise,
// good values are obtained with octaves=6 (numerber of noise calls),
// lacunarity=~2.0 (spacing between successive octaves: 2.0 for warpping
// output), gain=0.5 (relative weighting applied to each successive octave),
// offset=1.0 (used to invert the ridges).
inline float perlin_noise(const vec3f& p, const vec3i& wrap = zero3i);
inline float perlin_ridge(const vec3f& p, float lacunarity = 2,
    float gain = 0.5, int octaves = 6, float offset = 1,
    const vec3i& wrap = zero3i);
inline float perlin_fbm(const vec3f& p, float lacunarity = 2, float gain = 0.5,
    int octaves = 6, const vec3i& wrap = zero3i);
inline float perlin_turbulence(const vec3f& p, float lacunarity = 2,
    float gain = 0.5, int octaves = 6, const vec3i& wrap = zero3i);

}  // namespace yocto

// -----------------------------------------------------------------------------
//
//
// IMPLEMENTATION
//
//
// -----------------------------------------------------------------------------

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PERLIN NOISE
// -----------------------------------------------------------------------------
namespace yocto {

// not same permutation table as Perlin's reference to avoid copyright issues;
// Perlin's table can be found at http://mrl.nyu.edu/~perlin/noise/
// @OPTIMIZE: should this be unsigned char instead of int for cache?
inline const unsigned char p[512] = {
    // clang-format off
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
    // clang-format on
};

inline float perlin_noise(float x, float y, int wx, int wy) {
  auto ease   = [](float a) { return ((a * 6 - 15) * a + 10) * a * a * a; };
  auto ifloor = [](float a) -> int {
    int ai = (int)a;
    return (a < ai) ? ai - 1 : ai;
  };
  auto grad = [](int hash, float x, float y) -> float {
    int   h = hash & 7;       // Convert low 3 bits of hash code
    float u = h < 4 ? x : y;  // into 8 simple gradient directions,
    float v = h < 4 ? y : x;  // and compute the dot product with (x,y).
    return ((h & 1) ? -u : u) + ((h & 2) ? -2 * v : 2 * v);
  };
  auto hash = [](int x, int y) -> int { return p[p[x] + y]; };

  uint mx = (wx - 1) & 255, my = (wy - 1) & 255;
  auto ix = ifloor(x), iy = ifloor(y);
  auto fx = x - ix, fy = y - iy;
  auto x0 = ix & mx, y0 = iy & my;
  auto x1 = (ix + 1) & mx, y1 = (iy + 1) & my;

  auto ux = ease(fx), uy = ease(fy);

  auto n00 = grad(hash(x0, y0), fx, fy);
  auto n01 = grad(hash(x0, y1), fx, fy - 1);
  auto n10 = grad(hash(x1, y0), fx - 1, fy);
  auto n11 = grad(hash(x1, y1), fx - 1, fy - 1);

  auto n0 = lerp(n00, n01, uy);
  auto n1 = lerp(n10, n11, uy);

  return lerp(n0, n1, ux);
}

inline float perlin_noise(float x, float y, float z, int wx, int wy, int wz) {
  auto ease   = [](float a) { return ((a * 6 - 15) * a + 10) * a * a * a; };
  auto ifloor = [](float a) -> int {
    int ai = (int)a;
    return (a < ai) ? ai - 1 : ai;
  };
  auto grad = [](int hash, float x, float y, float z) -> float {
    auto h = hash & 15;      // Convert low 4 bits of hash code into 12 simple
    auto u = h < 8 ? x : y;  // gradient directions, and compute dot product.
    auto v = h < 4 ? y
                   : h == 12 || h == 14 ? x : z;  // Fix repeats at h = 12 to 15
    return ((h & 1) ? -u : u) + ((h & 2) ? -v : v);
  };
  auto hash = [](int x, int y, int z) -> int { return p[p[p[x] + y] + z]; };

  uint mx = (wx - 1) & 255, my = (wy - 1) & 255, mz = (wz - 1) & 255;
  auto ix = ifloor(x), iy = ifloor(y), iz = ifloor(z);
  auto fx = x - ix, fy = y - iy, fz = z - iz;
  auto x0 = ix & mx, y0 = iy & my, z0 = iz & mz;
  auto x1 = (ix + 1) & mx, y1 = (iy + 1) & my, z1 = (iz + 1) & mz;

  auto ux = ease(fx), uy = ease(fy), uz = ease(fz);

  auto n000 = grad(hash(x0, y0, z0), fx, fy, fz);
  auto n001 = grad(hash(x0, y0, z1), fx, fy, fz - 1);
  auto n010 = grad(hash(x0, y1, z0), fx, fy - 1, fz);
  auto n011 = grad(hash(x0, y1, z1), fx, fy - 1, fz - 1);
  auto n100 = grad(hash(x1, y0, z0), fx - 1, fy, fz);
  auto n101 = grad(hash(x1, y0, z1), fx - 1, fy, fz - 1);
  auto n110 = grad(hash(x1, y1, z0), fx - 1, fy - 1, fz);
  auto n111 = grad(hash(x1, y1, z1), fx - 1, fy - 1, fz - 1);

  auto n00 = lerp(n000, n001, uz);
  auto n01 = lerp(n010, n011, uz);
  auto n10 = lerp(n100, n101, uz);
  auto n11 = lerp(n110, n111, uz);

  auto n0 = lerp(n00, n01, uy);
  auto n1 = lerp(n10, n11, uy);

  return lerp(n0, n1, ux);
}

// noise
inline float perlin_noise(const vec3f& p, const vec3i& wrap) {
  return perlin_noise(p.x, p.y, p.z, wrap.x, wrap.y, wrap.z);
}

// ridge
inline float perlin_ridge(const vec3f& p, float lacunarity, float gain,
    int octaves, float offset, const vec3i& wrap) {
  auto frequency = 1.0f;
  auto prev      = 1.0f;
  auto amplitude = 0.5f;
  auto sum       = 0.0f;
  for (auto i = 0; i < octaves; i++) {
    auto r = perlin_noise(p * frequency, wrap);
    r      = r < 0 ? -r : r;  // fabs()
    r      = offset - r;
    r      = r * r;
    sum += r * amplitude * prev;
    prev = r;
    frequency *= lacunarity;
    amplitude *= gain;
  }
  return sum;
}

// fmb
inline float perlin_fbm(const vec3f& p, float lacunarity, float gain,
    int octaves, const vec3i& wrap) {
  auto frequency = 1.0f;
  auto amplitude = 1.0f;
  auto sum       = 0.0f;
  for (auto i = 0; i < octaves; i++) {
    sum += perlin_noise(p * frequency, wrap) * amplitude;
    frequency *= lacunarity;
    amplitude *= gain;
  }
  return sum;
}

// turbulence
inline float perlin_turbulence(const vec3f& p, float lacunarity, float gain,
    int octaves, const vec3i& wrap) {
  auto frequency = 1.0f;
  auto amplitude = 1.0f;
  auto sum       = 0.0f;
  for (auto i = 0; i < octaves; i++) {
    auto r = perlin_noise(p * frequency, wrap) * amplitude;
    r      = r < 0 ? -r : r;  // fabs()
    sum += r;
    frequency *= lacunarity;
    amplitude *= gain;
  }
  return sum;
}

}  // namespace yocto

#endif
