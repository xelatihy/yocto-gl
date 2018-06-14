//
// Implementation for Yocto/GL.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2018 Fabio Pellacini
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
// ------------------------------------------------------------------------------
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
// ------------------------------------------------------------------------------

#include "ygl.h"

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PERLIN NOISE
// -----------------------------------------------------------------------------
namespace ygl {

// clang-format off
float stb__perlin_lerp(float a, float b, float t)
{
   return a + (b-a) * t;
}

int stb__perlin_fastfloor(float a)
{
	int ai = (int) a;
	return (a < ai) ? ai-1 : ai;
}

// different grad function from Perlin's, but easy to modify to match reference
float stb__perlin_grad(int hash, float x, float y, float z)
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

float stb_perlin_noise3(float x, float y, float z, int x_wrap, int y_wrap, int z_wrap)
{
    // not same permutation table as Perlin's reference to avoid copyright issues;
    // Perlin's table can be found at http://mrl.nyu.edu/~perlin/noise/
    // @OPTIMIZE: should this be unsigned char instead of int for cache?
    static unsigned char stb__perlin_randtab[512] =
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
   int px = stb__perlin_fastfloor(x);
   int py = stb__perlin_fastfloor(y);
   int pz = stb__perlin_fastfloor(z);
   int x0 = px & x_mask, x1 = (px+1) & x_mask;
   int y0 = py & y_mask, y1 = (py+1) & y_mask;
   int z0 = pz & z_mask, z1 = (pz+1) & z_mask;
   int r0,r1, r00,r01,r10,r11;

   #define stb__perlin_ease(a)   (((a*6-15)*a + 10) * a * a * a)

   x -= px; u = stb__perlin_ease(x);
   y -= py; v = stb__perlin_ease(y);
   z -= pz; w = stb__perlin_ease(z);

   r0 = stb__perlin_randtab[x0];
   r1 = stb__perlin_randtab[x1];

   r00 = stb__perlin_randtab[r0+y0];
   r01 = stb__perlin_randtab[r0+y1];
   r10 = stb__perlin_randtab[r1+y0];
   r11 = stb__perlin_randtab[r1+y1];

   n000 = stb__perlin_grad(stb__perlin_randtab[r00+z0], x  , y  , z   );
   n001 = stb__perlin_grad(stb__perlin_randtab[r00+z1], x  , y  , z-1 );
   n010 = stb__perlin_grad(stb__perlin_randtab[r01+z0], x  , y-1, z   );
   n011 = stb__perlin_grad(stb__perlin_randtab[r01+z1], x  , y-1, z-1 );
   n100 = stb__perlin_grad(stb__perlin_randtab[r10+z0], x-1, y  , z   );
   n101 = stb__perlin_grad(stb__perlin_randtab[r10+z1], x-1, y  , z-1 );
   n110 = stb__perlin_grad(stb__perlin_randtab[r11+z0], x-1, y-1, z   );
   n111 = stb__perlin_grad(stb__perlin_randtab[r11+z1], x-1, y-1, z-1 );

   n00 = stb__perlin_lerp(n000,n001,w);
   n01 = stb__perlin_lerp(n010,n011,w);
   n10 = stb__perlin_lerp(n100,n101,w);
   n11 = stb__perlin_lerp(n110,n111,w);

   n0 = stb__perlin_lerp(n00,n01,v);
   n1 = stb__perlin_lerp(n10,n11,v);

   return stb__perlin_lerp(n0,n1,u);
}

float stb_perlin_ridge_noise3(float x, float y, float z,float lacunarity, float gain, float offset, int octaves,int x_wrap, int y_wrap, int z_wrap)
{
   int i;
   float frequency = 1.0f;
   float prev = 1.0f;
   float amplitude = 0.5f;
   float sum = 0.0f;

   for (i = 0; i < octaves; i++) {
      float r = (float)(stb_perlin_noise3(x*frequency,y*frequency,z*frequency,x_wrap,y_wrap,z_wrap));
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

float stb_perlin_fbm_noise3(float x, float y, float z,float lacunarity, float gain, int octaves,int x_wrap, int y_wrap, int z_wrap)
{
   int i;
   float frequency = 1.0f;
   float amplitude = 1.0f;
   float sum = 0.0f;

   for (i = 0; i < octaves; i++) {
      sum += stb_perlin_noise3(x*frequency,y*frequency,z*frequency,x_wrap,y_wrap,z_wrap)*amplitude;
      frequency *= lacunarity;
      amplitude *= gain;
   }
   return sum;
}

float stb_perlin_turbulence_noise3(float x, float y, float z, float lacunarity, float gain, int octaves,int x_wrap, int y_wrap, int z_wrap)
{
   int i;
   float frequency = 1.0f;
   float amplitude = 1.0f;
   float sum = 0.0f;

   for (i = 0; i < octaves; i++) {
      float r = stb_perlin_noise3(x*frequency,y*frequency,z*frequency,x_wrap,y_wrap,z_wrap)*amplitude;
      r = r<0 ? -r : r; // fabs()
      sum += r;
      frequency *= lacunarity;
      amplitude *= gain;
   }
   return sum;
}
// clang-format on

// adapeted  stb_perlin.h
float perlin_noise(vec3f p, vec3i wrap) {
    return stb_perlin_noise3(p.x, p.y, p.z, wrap.x, wrap.y, wrap.z);
}

// adapeted  stb_perlin.h
float perlin_ridge_noise(vec3f p, float lacunarity, float gain, float offset,
    int octaves, vec3i wrap) {
    return stb_perlin_ridge_noise3(p.x, p.y, p.z, lacunarity, gain, offset,
        octaves, wrap.x, wrap.y, wrap.z);
}

// adapeted  stb_perlin.h
float perlin_fbm_noise(
    vec3f p, float lacunarity, float gain, int octaves, vec3i wrap) {
    return stb_perlin_fbm_noise3(
        p.x, p.y, p.z, lacunarity, gain, octaves, wrap.x, wrap.y, wrap.z);
}

// adapeted  stb_perlin.h
float perlin_turbulence_noise(
    vec3f p, float lacunarity, float gain, int octaves, vec3i wrap) {
    return stb_perlin_turbulence_noise3(
        p.x, p.y, p.z, lacunarity, gain, octaves, wrap.x, wrap.y, wrap.z);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHAPE UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Compute per-vertex tangents for lines.
void compute_tangents(const std::vector<vec2i>& lines,
    const std::vector<vec3f>& pos, std::vector<vec3f>& norm) {
    norm.resize(pos.size());
    for (auto& n : norm) n = zero3f;
    for (auto& l : lines) {
        auto n = line_tangent(pos[l.x], pos[l.y]);
        norm[l.x] += n;
        norm[l.y] += n;
    }
    for (auto& n : norm) n = normalize(n);
}

// Compute per-vertex normals for triangles.
void compute_normals(const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, std::vector<vec3f>& norm) {
    norm.resize(pos.size());
    for (auto& n : norm) n = zero3f;
    for (auto& t : triangles) {
        auto n = triangle_normal(pos[t.x], pos[t.y], pos[t.z]);
        norm[t.x] += n;
        norm[t.y] += n;
        norm[t.z] += n;
    }
    for (auto& n : norm) n = normalize(n);
}

// Compute per-vertex normals for quads.
void compute_normals(const std::vector<vec4i>& quads,
    const std::vector<vec3f>& pos, std::vector<vec3f>& norm) {
    norm.resize(pos.size());
    for (auto& n : norm) n = zero3f;
    for (auto q : quads) {
        auto n = quad_normal(pos[q.x], pos[q.y], pos[q.z], pos[q.w]);
        norm[q.x] += n;
        norm[q.y] += n;
        norm[q.z] += n;
        if (q.z != q.w) norm[q.w] += n;
    }
    for (auto& n : norm) n = normalize(n);
}

// Compute per-vertex tangent frame for triangle meshes.
// Tangent space is defined by a four component vector.
// The first three components are the tangent with respect to the U texcoord.
// The fourth component is the sign of the tangent wrt the V texcoord.
// Tangent frame is useful in normal mapping.
void compute_tangent_space(const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, const std::vector<vec3f>& norm,
    const std::vector<vec2f>& texcoord, std::vector<vec4f>& tangsp,
    bool weighted) {
    auto tangu = std::vector<vec3f>(pos.size(), zero3f);
    auto tangv = std::vector<vec3f>(pos.size(), zero3f);
    for (auto t : triangles) {
        auto tutv = triangle_tangents_fromuv(pos[t.x], pos[t.y], pos[t.z],
            texcoord[t.x], texcoord[t.y], texcoord[t.z]);
        if (!weighted) tutv = {normalize(tutv.first), normalize(tutv.second)};
        for (auto vid : {t.x, t.y, t.z}) tangu[vid] += tutv.first;
        for (auto vid : {t.x, t.y, t.z}) tangv[vid] += tutv.second;
    }
    for (auto& t : tangu) t = normalize(t);
    for (auto& t : tangv) t = normalize(t);
    tangsp.resize(pos.size());
    for (auto& t : tangsp) t = zero4f;
    for (auto i = 0; i < pos.size(); i++) {
        tangu[i] = orthonormalize(tangu[i], norm[i]);
        auto s = (dot(cross(norm[i], tangu[i]), tangv[i]) < 0) ? -1.0f : 1.0f;
        tangsp[i] = {tangu[i].x, tangu[i].y, tangu[i].z, s};
    }
}

// Apply skinning
void compute_skinning(const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec4f>& weights,
    const std::vector<vec4i>& joints, const std::vector<mat4f>& xforms,
    std::vector<vec3f>& skinned_pos, std::vector<vec3f>& skinned_norm) {
    skinned_pos.resize(pos.size());
    skinned_norm.resize(norm.size());
    for (auto i = 0; i < pos.size(); i++) {
        skinned_pos[i] =
            transform_point(xforms[joints[i].x], pos[i]) * weights[i].x +
            transform_point(xforms[joints[i].y], pos[i]) * weights[i].y +
            transform_point(xforms[joints[i].z], pos[i]) * weights[i].z +
            transform_point(xforms[joints[i].w], pos[i]) * weights[i].w;
    }
    for (auto i = 0; i < pos.size(); i++) {
        skinned_norm[i] = normalize(
            transform_direction(xforms[joints[i].x], norm[i]) * weights[i].x +
            transform_direction(xforms[joints[i].y], norm[i]) * weights[i].y +
            transform_direction(xforms[joints[i].z], norm[i]) * weights[i].z +
            transform_direction(xforms[joints[i].w], norm[i]) * weights[i].w);
    }
}

// Apply skinning
void compute_skinning(const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec4f>& weights,
    const std::vector<vec4i>& joints, const std::vector<frame3f>& xforms,
    std::vector<vec3f>& skinned_pos, std::vector<vec3f>& skinned_norm) {
    skinned_pos.resize(pos.size());
    skinned_norm.resize(norm.size());
    for (auto i = 0; i < pos.size(); i++) {
        skinned_pos[i] =
            transform_point(xforms[joints[i].x], pos[i]) * weights[i].x +
            transform_point(xforms[joints[i].y], pos[i]) * weights[i].y +
            transform_point(xforms[joints[i].z], pos[i]) * weights[i].z +
            transform_point(xforms[joints[i].w], pos[i]) * weights[i].w;
    }
    for (auto i = 0; i < pos.size(); i++) {
        skinned_norm[i] = normalize(
            transform_direction(xforms[joints[i].x], norm[i]) * weights[i].x +
            transform_direction(xforms[joints[i].y], norm[i]) * weights[i].y +
            transform_direction(xforms[joints[i].z], norm[i]) * weights[i].z +
            transform_direction(xforms[joints[i].w], norm[i]) * weights[i].w);
    }
}

// Apply skinning as specified in Khronos glTF
void compute_matrix_skinning(const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec4f>& weights,
    const std::vector<vec4i>& joints, const std::vector<mat4f>& xforms,
    std::vector<vec3f>& skinned_pos, std::vector<vec3f>& skinned_norm) {
    skinned_pos.resize(pos.size());
    skinned_norm.resize(norm.size());
    for (auto i = 0; i < pos.size(); i++) {
        auto xform = xforms[joints[i].x] * weights[i].x +
                     xforms[joints[i].y] * weights[i].y +
                     xforms[joints[i].z] * weights[i].z +
                     xforms[joints[i].w] * weights[i].w;
        skinned_pos[i] = transform_point(xform, pos[i]);
        skinned_norm[i] = normalize(transform_direction(xform, norm[i]));
    }
}

// Create an array of edges.
std::vector<vec2i> get_edges(const std::vector<vec3i>& triangles) {
    auto edges = std::vector<vec2i>();
    auto eset = std::unordered_set<vec2i>();
    for (auto t : triangles) {
        for (auto e : {vec2i{t.x, t.y}, vec2i{t.y, t.z}, vec2i{t.z, t.x}}) {
            e = {min(e.x, e.y), max(e.x, e.y)};
            if (!eset.insert(e).second) continue;
            eset.insert({e.y, e.x});
            edges.push_back(e);
        }
    }
    return edges;
}

// Create an array of edges.
std::vector<vec2i> get_edges(const std::vector<vec4i>& quads) {
    auto edges = std::vector<vec2i>();
    auto eset = std::unordered_set<vec2i>();
    for (auto q : quads) {
        for (auto e : {vec2i{q.x, q.y}, vec2i{q.y, q.z}, vec2i{q.z, q.w},
                 vec2i{q.w, q.x}}) {
            if (e.x == e.y) continue;
            e = {min(e.x, e.y), max(e.x, e.y)};
            if (!eset.insert(e).second) continue;
            eset.insert({e.y, e.x});
            edges.push_back(e);
        }
    }
    return edges;
}

// Convert quads to triangles
std::vector<vec3i> convert_quads_to_triangles(const std::vector<vec4i>& quads) {
    auto triangles = std::vector<vec3i>();
    triangles.reserve(quads.size() * 2);
    for (auto q : quads) {
        triangles.push_back({q.x, q.y, q.w});
        if (q.z != q.w) triangles.push_back({q.z, q.w, q.y});
    }
    return triangles;
}

// Convert quads to triangles with a diamond-like topology.
// Quads have to be consecutive one row after another.
std::vector<vec3i> convert_quads_to_triangles(
    const std::vector<vec4i>& quads, int row_length) {
    auto triangles = std::vector<vec3i>();
    triangles.reserve(quads.size() * 2);
    for (auto q : quads) {
        triangles.push_back({q.x, q.y, q.w});
        if (q.z != q.w) triangles.push_back({q.z, q.w, q.y});
    }
    return triangles;
#if 0
        triangles.resize(usteps * vsteps * 2);
        for (auto j = 0; j < vsteps; j++) {
            for (auto i = 0; i < usteps; i++) {
                auto f1 = triangles[(j * usteps + i) * 2 + 0];
                auto f2 = triangles[(j * usteps + i) * 2 + 1];
                if ((i + j) % 2) {
                    f1 = {vid(i, j), vid(i + 1, j), vid(i + 1, j + 1)};
                    f2 = {vid(i + 1, j + 1), vid(i, j + 1), vid(i, j)};
                } else {
                    f1 = {vid(i, j), vid(i + 1, j), vid(i, j + 1)};
                    f2 = {vid(i + 1, j + 1), vid(i, j + 1), vid(i + 1, j)};
                }
            }
        }
#endif
    return triangles;
}

// Convert beziers to lines using 3 lines for each bezier.
std::vector<vec2i> convert_bezier_to_lines(const std::vector<vec4i>& beziers) {
    auto lines = std::vector<vec2i>();
    lines.reserve(beziers.size() * 3);
    for (auto b : beziers) {
        lines.push_back({b.x, b.y});
        lines.push_back({b.y, b.z});
        lines.push_back({b.z, b.w});
    }
    return lines;
}

// Convert face varying data to single primitives. Returns the quads indices
// and filled vectors for pos, norm and texcoord.
void convert_face_varying(std::vector<vec4i>& qquads, std::vector<vec3f>& qpos,
    std::vector<vec3f>& qnorm, std::vector<vec2f>& qtexcoord,
    std::vector<vec4f>& qcolor, const std::vector<vec4i>& quads_pos,
    const std::vector<vec4i>& quads_norm,
    const std::vector<vec4i>& quads_texcoord,
    const std::vector<vec4i>& quads_color, const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec2f>& texcoord,
    const std::vector<vec4f>& color) {
    // make faces unique
    std::unordered_map<vec4i, int> vert_map;
    qquads = std::vector<vec4i>(quads_pos.size());
    for (auto fid = 0; fid < quads_pos.size(); fid++) {
        for (auto c = 0; c < 4; c++) {
            auto v = vec4i{
                (&quads_pos[fid].x)[c],
                (!quads_norm.empty()) ? (&quads_norm[fid].x)[c] : -1,
                (!quads_texcoord.empty()) ? (&quads_texcoord[fid].x)[c] : -1,
                (!quads_color.empty()) ? (&quads_color[fid].x)[c] : -1,
            };
            auto it = vert_map.find(v);
            if (it == vert_map.end()) {
                auto s = (int)vert_map.size();
                vert_map.insert(it, {v, s});
                (&qquads[fid].x)[c] = s;
            } else {
                (&qquads[fid].x)[c] = it->second;
            }
        }
    }

    // fill vert data
    qpos.clear();
    if (!pos.empty()) {
        qpos.resize(vert_map.size());
        for (auto kv : vert_map) { qpos[kv.second] = pos[kv.first.x]; }
    }
    qnorm.clear();
    if (!norm.empty()) {
        qnorm.resize(vert_map.size());
        for (auto kv : vert_map) { qnorm[kv.second] = norm[kv.first.y]; }
    }
    qtexcoord.clear();
    if (!texcoord.empty()) {
        qtexcoord.resize(vert_map.size());
        for (auto kv : vert_map) {
            qtexcoord[kv.second] = texcoord[kv.first.z];
        }
    }
}

// Subdivide lines.
template <typename T>
void subdivide_lines(
    std::vector<vec2i>& lines, std::vector<T>& vert, int level) {
    if (lines.empty() || vert.empty()) return;

    for (auto l = 0; l < level; l++) {
        auto tvert = vert;
        auto tlines = std::vector<vec2i>();
        for (auto l : lines) {
            tlines.push_back({l.x, (int)tvert.size()});
            tlines.push_back({(int)tvert.size(), l.y});
            tvert.push_back((vert[l.x] + vert[l.y]) / 2);
        }

        std::swap(vert, tvert);
        std::swap(lines, tlines);
    }
}

template void subdivide_lines(std::vector<vec2i>&, std::vector<float>&, int);
template void subdivide_lines(std::vector<vec2i>&, std::vector<vec2f>&, int);
template void subdivide_lines(std::vector<vec2i>&, std::vector<vec3f>&, int);
template void subdivide_lines(std::vector<vec2i>&, std::vector<vec4f>&, int);

// Subdivide triangle.
template <typename T>
void subdivide_triangles(
    std::vector<vec3i>& triangles, std::vector<T>& vert, int level) {
    if (triangles.empty() || vert.empty()) return;

    for (auto l = 0; l < level; l++) {
        auto tvert = vert;
        auto ttriangles = std::vector<vec3i>();
        auto emap = std::unordered_map<vec2i, int>();
        for (auto t : triangles) {
            for (auto e : {vec2i{t.x, t.y}, vec2i{t.y, t.z}, vec2i{t.z, t.x}}) {
                if (emap.find(e) != emap.end()) continue;
                emap[{e.x, e.y}] = (int)tvert.size();
                emap[{e.y, e.x}] = (int)tvert.size();
                tvert.push_back((vert[e.x] + vert[e.y]) / 2);
            }
            ttriangles.push_back(
                {t.x, emap.at({t.x, t.y}), emap.at({t.z, t.x})});
            ttriangles.push_back(
                {t.y, emap.at({t.y, t.z}), emap.at({t.x, t.y})});
            ttriangles.push_back(
                {t.z, emap.at({t.z, t.x}), emap.at({t.y, t.z})});
            ttriangles.push_back({emap.at({t.x, t.y}), emap.at({t.y, t.z}),
                emap.at({t.z, t.x})});
        }

        std::swap(vert, tvert);
        std::swap(triangles, ttriangles);
    }
}

template void subdivide_triangles(
    std::vector<vec3i>&, std::vector<float>&, int);
template void subdivide_triangles(
    std::vector<vec3i>&, std::vector<vec2f>&, int);
template void subdivide_triangles(
    std::vector<vec3i>&, std::vector<vec3f>&, int);
template void subdivide_triangles(
    std::vector<vec3i>&, std::vector<vec4f>&, int);

// Subdivide quads.
template <typename T>
void subdivide_quads(
    std::vector<vec4i>& quads, std::vector<T>& vert, int level) {
    if (quads.empty() || vert.empty()) return;

    for (auto l = 0; l < level; l++) {
        auto tvert = vert;
        auto tquads = std::vector<vec4i>();
        auto emap = std::unordered_map<vec2i, int>();
        for (auto q : quads) {
            for (auto e : {vec2i{q.x, q.y}, vec2i{q.y, q.z}, vec2i{q.z, q.w},
                     vec2i{q.w, q.x}}) {
                if (e.x == e.y) continue;
                if (emap.find(e) != emap.end()) continue;
                emap[{e.x, e.y}] = (int)tvert.size();
                emap[{e.y, e.x}] = (int)tvert.size();
                tvert.push_back((vert[e.x] + vert[e.y]) / 2);
            }
            if (q.z != q.w) {
                tquads.push_back({q.x, emap.at({q.x, q.y}), (int)tvert.size(),
                    emap.at({q.w, q.x})});
                tquads.push_back({q.y, emap.at({q.y, q.z}), (int)tvert.size(),
                    emap.at({q.x, q.y})});
                tquads.push_back({q.z, emap.at({q.z, q.w}), (int)tvert.size(),
                    emap.at({q.y, q.z})});
                tquads.push_back({q.w, emap.at({q.w, q.x}), (int)tvert.size(),
                    emap.at({q.z, q.w})});
                tvert.push_back(
                    (vert[q.x] + vert[q.y] + vert[q.z] + vert[q.w]) / 4);
            } else {
                tquads.push_back({q.x, emap.at({q.x, q.y}), (int)tvert.size(),
                    emap.at({q.z, q.x})});
                tquads.push_back({q.y, emap.at({q.y, q.z}), (int)tvert.size(),
                    emap.at({q.x, q.y})});
                tquads.push_back({q.z, emap.at({q.z, q.x}), (int)tvert.size(),
                    emap.at({q.y, q.z})});
                tvert.push_back((vert[q.x] + vert[q.y] + vert[q.y]) / 3);
            }
        }

        std::swap(vert, tvert);
        std::swap(quads, tquads);
    }
}

template void subdivide_quads(std::vector<vec4i>&, std::vector<float>&, int);
template void subdivide_quads(std::vector<vec4i>&, std::vector<vec2f>&, int);
template void subdivide_quads(std::vector<vec4i>&, std::vector<vec3f>&, int);
template void subdivide_quads(std::vector<vec4i>&, std::vector<vec4f>&, int);

// Subdivide beziers.
template <typename T>
void subdivide_beziers(
    std::vector<vec4i>& beziers, std::vector<T>& vert, int level) {
    if (beziers.empty() || vert.empty()) return;

    for (auto l = 0; l < level; l++) {
        auto vmap = std::unordered_map<int, int>();
        auto tvert = std::vector<T>();
        auto tbeziers = std::vector<vec4i>();
        for (auto b : beziers) {
            if (vmap.find(b.x) == vmap.end()) {
                vmap[b.x] = (int)tvert.size();
                tvert.push_back(vert[b.x]);
            }
            if (vmap.find(b.w) == vmap.end()) {
                vmap[b.w] = (int)tvert.size();
                tvert.push_back(vert[b.w]);
            }
            auto bo = (int)tvert.size();
            tbeziers.push_back({vmap.at(b.x), bo + 0, bo + 1, bo + 2});
            tbeziers.push_back({bo + 2, bo + 3, bo + 4, vmap.at(b.w)});
            tvert.push_back(vert[b.x] / 2 + vert[b.y] / 2);
            tvert.push_back(vert[b.x] / 4 + vert[b.y] / 2 + vert[b.z] / 4);
            tvert.push_back(vert[b.x] / 8 + 3 * vert[b.y] / 8 +
                            3 * vert[b.z] / 8 + vert[b.w] / 8);
            tvert.push_back(vert[b.y] / 4 + vert[b.z] / 2 + vert[b.w] / 4);
            tvert.push_back(vert[b.z] / 2 + vert[b.w] / 2);
        }

        std::swap(vert, tvert);
        std::swap(beziers, tbeziers);
    }
}

template void subdivide_beziers(std::vector<vec4i>&, std::vector<float>&, int);
template void subdivide_beziers(std::vector<vec4i>&, std::vector<vec2f>&, int);
template void subdivide_beziers(std::vector<vec4i>&, std::vector<vec3f>&, int);
template void subdivide_beziers(std::vector<vec4i>&, std::vector<vec4f>&, int);

// Subdivide catmullclark.
template <typename T>
void subdivide_catmullclark(std::vector<vec4i>& quads, std::vector<T>& vert,
    int level, bool lock_boundary) {
    if (quads.empty() || vert.empty()) return;

    for (auto l = 0; l < level; l++) {
        // split elements ------------------------------------
        auto tvert = vert;
        auto tquads = std::vector<vec4i>();
        auto emap = std::unordered_map<vec2i, int>();
        for (auto q : quads) {
            for (auto e : {vec2i{q.x, q.y}, vec2i{q.y, q.z}, vec2i{q.z, q.w},
                     vec2i{q.w, q.x}}) {
                if (e.x == e.y) continue;
                if (emap.find(e) != emap.end()) continue;
                emap[{e.x, e.y}] = (int)tvert.size();
                emap[{e.y, e.x}] = (int)tvert.size();
                tvert.push_back((vert[e.x] + vert[e.y]) / 2);
            }
            if (q.z != q.w) {
                tquads.push_back({q.x, emap.at({q.x, q.y}), (int)tvert.size(),
                    emap.at({q.w, q.x})});
                tquads.push_back({q.y, emap.at({q.y, q.z}), (int)tvert.size(),
                    emap.at({q.x, q.y})});
                tquads.push_back({q.z, emap.at({q.z, q.w}), (int)tvert.size(),
                    emap.at({q.y, q.z})});
                tquads.push_back({q.w, emap.at({q.w, q.x}), (int)tvert.size(),
                    emap.at({q.z, q.w})});
                tvert.push_back(
                    (vert[q.x] + vert[q.y] + vert[q.z] + vert[q.w]) / 4);
            } else {
                tquads.push_back({q.x, emap.at({q.x, q.y}), (int)tvert.size(),
                    emap.at({q.z, q.x})});
                tquads.push_back({q.y, emap.at({q.y, q.z}), (int)tvert.size(),
                    emap.at({q.x, q.y})});
                tquads.push_back({q.z, emap.at({q.z, q.x}), (int)tvert.size(),
                    emap.at({q.y, q.z})});
                tvert.push_back((vert[q.x] + vert[q.y] + vert[q.y]) / 3);
            }
        }

        // find boundary -------------------------------------
        auto eset = std::unordered_set<vec2i>();
        for (auto q : quads) {
            for (auto e : {vec2i{q.x, q.y}, vec2i{q.y, q.z}, vec2i{q.z, q.w},
                     vec2i{q.w, q.x}}) {
                if (e.x == e.y) continue;
                eset.insert({e.x, e.y});
            }
        }
        auto boundary = std::vector<vec2i>();
        for (auto e : eset) {
            if (eset.find({e.y, e.x}) != eset.end()) continue;
            boundary.push_back(e);
        }

        // split boundary
        auto tboundary = std::vector<vec2i>();
        for (auto e : boundary) {
            tboundary.push_back({e.x, emap.at(e)});
            tboundary.push_back({emap.at(e), e.y});
        }

        // setup creases -----------------------------------
        auto tcrease_edges = std::vector<vec2i>();
        auto tcrease_verts = std::vector<int>();
        if (lock_boundary) {
            for (auto b : tboundary) {
                tcrease_verts.push_back(b.x);
                tcrease_verts.push_back(b.y);
            }
        } else {
            for (auto b : tboundary) tcrease_edges.push_back(b);
        }

        // define vertex valence ---------------------------
        auto tvert_val = std::vector<int>(tvert.size(), 2);
        for (auto e : tboundary) {
            tvert_val[e.x] = (lock_boundary) ? 0 : 1;
            tvert_val[e.y] = (lock_boundary) ? 0 : 1;
        }

        // averaging pass ----------------------------------
        auto avert = std::vector<T>(tvert.size(), T());
        auto acount = std::vector<int>(tvert.size(), 0);
        for (auto p : tcrease_verts) {
            if (tvert_val[p] != 0) continue;
            avert[p] += tvert[p];
            acount[p] += 1;
        }
        for (auto e : tcrease_edges) {
            auto c = (tvert[e.x] + tvert[e.y]) / 2.0f;
            for (auto vid : {e.x, e.y}) {
                if (tvert_val[vid] != 1) continue;
                avert[vid] += c;
                acount[vid] += 1;
            }
        }
        for (auto q : tquads) {
            auto c = (tvert[q.x] + tvert[q.y] + tvert[q.z] + tvert[q.w]) / 4.0f;
            for (auto vid : {q.x, q.y, q.z, q.w}) {
                if (tvert_val[vid] != 2) continue;
                avert[vid] += c;
                acount[vid] += 1;
            }
        }
        for (auto i = 0; i < tvert.size(); i++) avert[i] /= (float)acount[i];

        // correction pass ----------------------------------
        // p = p + (avg_p - p) * (4/avg_count)
        for (auto i = 0; i < tvert.size(); i++) {
            if (tvert_val[i] != 2) continue;
            avert[i] = tvert[i] + (avert[i] - tvert[i]) * (4.0f / acount[i]);
        }
        tvert = avert;

        std::swap(vert, tvert);
        std::swap(quads, tquads);
    }
}

template void subdivide_catmullclark(
    std::vector<vec4i>&, std::vector<float>&, int, bool);
template void subdivide_catmullclark(
    std::vector<vec4i>&, std::vector<vec2f>&, int, bool);
template void subdivide_catmullclark(
    std::vector<vec4i>&, std::vector<vec3f>&, int, bool);
template void subdivide_catmullclark(
    std::vector<vec4i>&, std::vector<vec4f>&, int, bool);

// Merge lines between shapes.
void merge_lines(std::vector<vec2i>& lines, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    const std::vector<vec2i>& lines1, const std::vector<vec3f>& pos1,
    const std::vector<vec3f>& norm1, const std::vector<vec2f>& texcoord1) {
    lines.reserve(lines.size() + lines1.size());
    auto nverts = (int)pos.size();
    for (auto l : lines1) lines.push_back({l.x + nverts, l.y + nverts});
    for (auto v : pos1) pos.push_back(v);
    for (auto v : norm1) norm.push_back(v);
    for (auto v : texcoord1) texcoord.push_back(v);
}

// Merge triangles between shapes.
void merge_triangles(std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    const std::vector<vec3i>& triangles1, const std::vector<vec3f>& pos1,
    const std::vector<vec3f>& norm1, const std::vector<vec2f>& texcoord1) {
    triangles.reserve(triangles.size() + triangles1.size());
    auto nverts = (int)pos.size();
    for (auto t : triangles1)
        triangles.push_back({t.x + nverts, t.y + nverts, t.z + nverts});
    for (auto v : pos1) pos.push_back(v);
    for (auto v : norm1) norm.push_back(v);
    for (auto v : texcoord1) texcoord.push_back(v);
}

// Merge quads between shapes.
void merge_quads(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    const std::vector<vec4i>& quads1, const std::vector<vec3f>& pos1,
    const std::vector<vec3f>& norm1, const std::vector<vec2f>& texcoord1) {
    quads.reserve(quads.size() + quads1.size());
    auto nverts = (int)pos.size();
    for (auto q : quads1)
        quads.push_back(
            {q.x + nverts, q.y + nverts, q.z + nverts, q.w + nverts});
    for (auto v : pos1) pos.push_back(v);
    for (auto v : norm1) norm.push_back(v);
    for (auto v : texcoord1) texcoord.push_back(v);
}

// Weld vertices within a threshold. For noe the implementation is O(n^2).
std::vector<int> weld_vertices(std::vector<vec3f>& pos, float threshold) {
    auto vid = std::vector<int>(pos.size());
    auto wpos = std::vector<vec3f>();
    for (auto i = 0; i < pos.size(); i++) {
        vid[i] = (int)wpos.size();
        for (auto j = 0; j < wpos.size(); j++) {
            if (length(pos[i] - wpos[j]) < threshold) {
                vid[i] = j;
                break;
            }
        }
        if (vid[i] == (int)wpos.size()) wpos.push_back(pos[i]);
    }
    std::swap(pos, wpos);
    return vid;
}
void weld_triangles(
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos, float threshold) {
    auto vid = weld_vertices(pos, threshold);
    for (auto t : triangles) {
        t.x = vid[t.x];
        t.y = vid[t.y];
        t.z = vid[t.z];
    }
}
void weld_quads(
    std::vector<vec4i>& quads, std::vector<vec3f>& pos, float threshold) {
    auto vid = weld_vertices(pos, threshold);
    for (auto& q : quads) {
        q.x = vid[q.x];
        q.y = vid[q.y];
        q.z = vid[q.z];
        q.w = vid[q.w];
    }
}

// Samples a set of points over a triangle mesh uniformly. The rng function
// takes the point index and returns vec3f numbers uniform directibuted in
// [0,1]^3. unorm and texcoord are optional.
std::tuple<std::vector<vec3f>, std::vector<vec3f>, std::vector<vec2f>>
sample_triangles_points(const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, const std::vector<vec3f>& norm,
    const std::vector<vec2f>& texcoord, int npoints, int seed) {
    auto sampled_pos = std::vector<vec3f>(npoints);
    auto sampled_norm = std::vector<vec3f>(npoints);
    auto sampled_texcoord = std::vector<vec2f>(npoints);
    auto cdf = sample_triangles_cdf(triangles, pos);
    auto rng = make_rng(seed);
    for (auto i = 0; i < npoints; i++) {
        auto ei = 0;
        auto uv = zero2f;
        std::tie(ei, uv) =
            sample_triangles(cdf, rand1f(rng), {rand1f(rng), rand1f(rng)});
        auto t = triangles[ei];
        sampled_pos[i] = interpolate_triangle(pos[t.x], pos[t.y], pos[t.z], uv);
        if (!sampled_norm.empty()) {
            sampled_norm[i] = normalize(
                interpolate_triangle(norm[t.x], norm[t.y], norm[t.z], uv));
        } else {
            sampled_norm[i] = triangle_normal(pos[t.x], pos[t.y], pos[t.z]);
        }
        if (!sampled_texcoord.empty()) {
            sampled_texcoord[i] = interpolate_triangle(
                texcoord[t.x], texcoord[t.y], texcoord[t.z], uv);
        } else {
            sampled_texcoord[i] = zero2f;
        }
    }

    return {sampled_pos, sampled_norm, sampled_texcoord};
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENRTATION OF RAY-PRIMITIVE INTERSECTION FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Intersect a ray with a point (approximate)
bool intersect_point(
    const ray3f& ray, vec3f p, float r, float& dist, vec2f& uv) {
    // find parameter for line-point minimum distance
    auto w = p - ray.o;
    auto t = dot(w, ray.d) / dot(ray.d, ray.d);

    // exit if not within bounds
    if (t < ray.tmin || t > ray.tmax) return false;

    // test for line-point distance vs point radius
    auto rp = ray.o + ray.d * t;
    auto prp = p - rp;
    if (dot(prp, prp) > r * r) return false;

    // intersection occurred: set params and exit
    dist = t;
    uv = {0, 0};

    return true;
}

// Intersect a ray with a line
bool intersect_line(const ray3f& ray, const vec3f& v0, const vec3f& v1,
    float r0, float r1, float& dist, vec2f& uv) {
    // setup intersection params
    auto u = ray.d;
    auto v = v1 - v0;
    auto w = ray.o - v0;

    // compute values to solve a linear system
    auto a = dot(u, u);
    auto b = dot(u, v);
    auto c = dot(v, v);
    auto d = dot(u, w);
    auto e = dot(v, w);
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
    s = clamp(s, (float)0, (float)1);

    // compute segment-segment distance on the closest points
    auto p0 = ray.o + ray.d * t;
    auto p1 = v0 + (v1 - v0) * s;
    auto p01 = p0 - p1;

    // check with the line radius at the same point
    auto d2 = dot(p01, p01);
    auto r = r0 * (1 - s) + r1 * s;
    if (d2 > r * r) return false;

    // intersection occurred: set params and exit
    dist = t;
    uv = {s, sqrt(d2) / r};

    return true;
}

// Intersect a ray with a triangle
bool intersect_triangle(const ray3f& ray, const vec3f& v0, const vec3f& v1,
    const vec3f& v2, float& dist, vec2f& uv) {
    // compute triangle edges
    auto edge1 = v1 - v0;
    auto edge2 = v2 - v0;

    // compute determinant to solve a linear system
    auto pvec = cross(ray.d, edge2);
    auto det = dot(edge1, pvec);

    // check determinant and exit if triangle and ray are parallel
    // (could use EPSILONS if desired)
    if (det == 0) return false;
    auto inv_det = 1.0f / det;

    // compute and check first bricentric coordinated
    auto tvec = ray.o - v0;
    auto u = dot(tvec, pvec) * inv_det;
    if (u < 0 || u > 1) return false;

    // compute and check second bricentric coordinated
    auto qvec = cross(tvec, edge1);
    auto v = dot(ray.d, qvec) * inv_det;
    if (v < 0 || u + v > 1) return false;

    // compute and check ray parameter
    auto t = dot(edge2, qvec) * inv_det;
    if (t < ray.tmin || t > ray.tmax) return false;

    // intersection occurred: set params and exit
    dist = t;
    uv = {u, v};

    return true;
}

// Intersect a ray with a quad.
bool intersect_quad(const ray3f& ray, const vec3f& v0, const vec3f& v1,
    const vec3f& v2, const vec3f& v3, float& dist, vec2f& uv) {
    auto hit = false;
    auto tray = ray;
    if (intersect_triangle(tray, v0, v1, v3, dist, uv)) {
        tray.tmax = dist;
        hit = true;
    }
    if (intersect_triangle(tray, v2, v3, v1, dist, uv)) {
        uv = {1 - uv.x, 1 - uv.y};
        tray.tmax = dist;
        hit = true;
    }
    return hit;
}

// Min/max used in BVH traversal. Copied here since the traversal code
// relies on the specific behaviour wrt NaNs.
static inline const float& _safemin(const float& a, const float& b) {
    return (a < b) ? a : b;
}
// Min/max used in BVH traversal. Copied here since the traversal code
// relies on the specific behaviour wrt NaNs.
static inline const float& _safemax(const float& a, const float& b) {
    return (a > b) ? a : b;
}

// Intersect a ray with a axis-aligned bounding box
bool intersect_bbox(const ray3f& ray, const bbox3f& bbox) {
    // determine intersection ranges
    auto invd = 1.0f / ray.d;
    auto t0 = (bbox.min - ray.o) * invd;
    auto t1 = (bbox.max - ray.o) * invd;
    // flip based on range directions
    if (invd.x < 0.0f) std::swap(t0.x, t1.x);
    if (invd.y < 0.0f) std::swap(t0.y, t1.y);
    if (invd.z < 0.0f) std::swap(t0.z, t1.z);
    auto tmin = _safemax(t0.z, _safemax(t0.y, _safemax(t0.x, ray.tmin)));
    auto tmax = _safemin(t1.z, _safemin(t1.y, _safemin(t1.x, ray.tmax)));
    tmax *= 1.00000024f;  // for double: 1.0000000000000004
    return tmin <= tmax;
}

// Intersect a ray with a axis-aligned bounding box
bool intersect_bbox(
    const ray3f& ray, vec3f ray_dinv, vec3i ray_dsign, const bbox3f& bbox_) {
    auto bbox = &bbox_.min;
    auto txmin = (bbox[ray_dsign.x].x - ray.o.x) * ray_dinv.x;
    auto txmax = (bbox[1 - ray_dsign.x].x - ray.o.x) * ray_dinv.x;
    auto tymin = (bbox[ray_dsign.y].y - ray.o.y) * ray_dinv.y;
    auto tymax = (bbox[1 - ray_dsign.y].y - ray.o.y) * ray_dinv.y;
    auto tzmin = (bbox[ray_dsign.z].z - ray.o.z) * ray_dinv.z;
    auto tzmax = (bbox[1 - ray_dsign.z].z - ray.o.z) * ray_dinv.z;
    auto tmin = _safemax(tzmin, _safemax(tymin, _safemax(txmin, ray.tmin)));
    auto tmax = _safemin(tzmax, _safemin(tymax, _safemin(txmax, ray.tmax)));
    tmax *= 1.00000024f;  // for double: 1.0000000000000004
    return tmin <= tmax;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENRTATION OF POINT-PRIMITIVE DISTANCE FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// TODO: documentation
bool overlap_point(const vec3f& pos, float dist_max, const vec3f& p, float r,
    float& dist, vec2f& uv) {
    auto d2 = dot(pos - p, pos - p);
    if (d2 > (dist_max + r) * (dist_max + r)) return false;
    dist = sqrt(d2);
    uv = {0, 0};
    return true;
}

// TODO: documentation
float closestuv_line(const vec3f& pos, const vec3f& v0, const vec3f& v1) {
    auto ab = v1 - v0;
    auto d = dot(ab, ab);
    // Project c onto ab, computing parameterized position d(t) = a + t*(b –
    // a)
    auto u = dot(pos - v0, ab) / d;
    u = clamp(u, (float)0, (float)1);
    return u;
}

// TODO: documentation
bool overlap_line(const vec3f& pos, float dist_max, const vec3f& v0,
    const vec3f& v1, float r0, float r1, float& dist, vec2f& uv) {
    auto u = closestuv_line(pos, v0, v1);
    // Compute projected position from the clamped t d = a + t * ab;
    auto p = v0 + (v1 - v0) * u;
    auto r = r0 + (r1 - r0) * u;
    auto d2 = dot(pos - p, pos - p);
    // check distance
    if (d2 > (dist_max + r) * (dist_max + r)) return false;
    // done
    dist = sqrt(d2);
    uv = {u, 0};
    return true;
}

// TODO: documentation
// this is a complicated test -> I probably "--"+prefix to use a sequence of
// test (triangle body, and 3 edges)
vec2f closestuv_triangle(
    const vec3f& pos, const vec3f& v0, const vec3f& v1, const vec3f& v2) {
    auto ab = v1 - v0;
    auto ac = v2 - v0;
    auto ap = pos - v0;

    auto d1 = dot(ab, ap);
    auto d2 = dot(ac, ap);

    // corner and edge cases
    if (d1 <= 0 && d2 <= 0) return {0, 0};

    auto bp = pos - v1;
    auto d3 = dot(ab, bp);
    auto d4 = dot(ac, bp);
    if (d3 >= 0 && d4 <= d3) return {1, 0};

    auto vc = d1 * d4 - d3 * d2;
    if ((vc <= 0) && (d1 >= 0) && (d3 <= 0)) return {d1 / (d1 - d3), 0};

    auto cp = pos - v2;
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
    auto u = vb * denom;
    auto v = vc * denom;
    return {u, v};
}

// TODO: documentation
bool overlap_triangle(const vec3f& pos, float dist_max, const vec3f& v0,
    const vec3f& v1, const vec3f& v2, float r0, float r1, float r2, float& dist,
    vec2f& uv) {
    uv = closestuv_triangle(pos, v0, v1, v2);
    auto p = interpolate_triangle(v0, v1, v2, uv);
    auto r = interpolate_triangle(r0, r1, r2, uv);
    auto dd = dot(p - pos, p - pos);
    if (dd > (dist_max + r) * (dist_max + r)) return false;
    dist = sqrt(dd);
    return true;
}

// TODO: documentation
bool overlap_quad(const vec3f& pos, float dist_max, const vec3f& v0,
    const vec3f& v1, const vec3f& v2, const vec3f& v3, float r0, float r1,
    float r2, float r3, float& dist, vec2f& uv) {
    auto hit = false;
    if (overlap_triangle(pos, dist_max, v0, v1, v3, r0, r1, r3, dist, uv)) {
        dist_max = dist;
        hit = true;
    }
    if (overlap_triangle(pos, dist_max, v2, v3, v1, r2, r3, r1, dist, uv)) {
        // dist_max = dist;
        uv = {1 - uv.x, 1 - uv.y};
        hit = true;
    }
    return hit;
}

// TODO: documentation
bool overlap_tetrahedron(const vec3f& pos, const vec3f& v0, const vec3f& v1,
    const vec3f& v2, const vec3f& v3, vec4f& uv) {
    // TODO: fix uv
    auto vol = dot(v3 - v0, cross(v3 - v1, v3 - v0));
    if (vol == 0) return false;
    auto u = dot(v3 - v0, cross(v3 - v1, v3 - v0)) / vol;
    if (u < 0 || u > 1) return false;
    auto v = dot(v3 - v0, cross(v3 - v1, v3 - v0)) / vol;
    if (v < 0 || v > 1 || u + v > 1) return false;
    auto w = dot(v3 - v0, cross(v3 - v1, v3 - v0)) / vol;
    if (w < 0 || w > 1 || u + v + w > 1) return false;
    uv = {u, v, w, 1 - u - v - w};
    return true;
}

// TODO: documentation
bool overlap_tetrahedron(const vec3f& pos, float dist_max, const vec3f& v0,
    const vec3f& v1, const vec3f& v2, const vec3f& v3, float r0, float r1,
    float r2, float r3, float& dist, vec4f& uv) {
    // TODO: FIX UVs
    // check interior
    if (overlap_tetrahedron(pos, v0, v1, v2, v3, uv)) {
        dist = 0;
        return true;
    }

    // check faces
    auto hit = false;
    auto tuv = zero2f;
    if (overlap_triangle(pos, dist_max, v0, v1, v2, r0, r1, r2, dist, tuv)) {
        hit = true;
        dist_max = dist;
    }
    if (overlap_triangle(pos, dist_max, v0, v1, v3, r0, r1, r3, dist, tuv)) {
        hit = true;
        dist_max = dist;
    }
    if (overlap_triangle(pos, dist_max, v0, v2, v3, r0, r2, r3, dist, tuv)) {
        hit = true;
        dist_max = dist;
    }
    if (overlap_triangle(pos, dist_max, v1, v2, v3, r1, r2, r3, dist, tuv)) {
        hit = true;
        // dist_max = dist;
    }

    return hit;
}

// TODO: documentation
bool distance_check_bbox(const vec3f& pos, float dist_max, const bbox3f& bbox) {
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

// TODO: doc
bool overlap_bbox(const bbox3f& bbox1, const bbox3f& bbox2) {
    if (bbox1.max.x < bbox2.min.x || bbox1.min.x > bbox2.max.x) return false;
    if (bbox1.max.y < bbox2.min.y || bbox1.min.y > bbox2.max.y) return false;
    if (bbox1.max.z < bbox2.min.z || bbox1.min.z > bbox2.max.z) return false;
    return true;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR BVH
// -----------------------------------------------------------------------------
namespace ygl {

// BVH primitive with its bbox, its center and the index to the primitive
struct bvh_prim {
    bbox3f bbox = invalid_bbox3f;
    vec3f center = zero3f;
    int primid = 0;
};

// Initializes the BVH node node that contains the primitives sorted_prims
// from start to end, by either splitting it into two other nodes,
// or initializing it as a leaf. When splitting, the heuristic heuristic is
// used and nodes added sequentially in the preallocated nodes array and
// the number of nodes nnodes is updated.
int make_bvh_node(std::vector<bvh_node>& nodes, std::vector<bvh_prim>& prims,
    int start, int end, bvh_node_type type, bool equal_size) {
    // add a new node
    auto nodeid = (int)nodes.size();
    nodes.push_back({});
    auto& node = nodes.back();

    // compute bounds
    node.bbox = invalid_bbox3f;
    for (auto i = start; i < end; i++) node.bbox += prims[i].bbox;

    // split into two children
    if (end - start > bvh_max_prims) {
        // initialize split axis and position
        auto split_axis = 0;
        auto mid = (start + end) / 2;

        // compute primintive bounds and size
        auto cbbox = invalid_bbox3f;
        for (auto i = start; i < end; i++) cbbox += prims[i].center;
        auto csize = cbbox.max - cbbox.min;

        // choose the split axis and position
        if (csize != zero3f) {
            // split along largest
            auto largest_axis = 0;
            if (csize.x >= csize.y && csize.x >= csize.z) largest_axis = 0;
            if (csize.y >= csize.x && csize.y >= csize.z) largest_axis = 1;
            if (csize.z >= csize.x && csize.z >= csize.y) largest_axis = 2;

            // check heuristic
            if (equal_size) {
                // split the space in the middle along the largest axis
                split_axis = largest_axis;
                auto csize = (cbbox.max + cbbox.min) / 2;
                auto middle = (&csize.x)[largest_axis];
                mid = (int)(std::partition(prims.data() + start,
                                prims.data() + end,
                                [split_axis, middle](auto& a) {
                                    return (&a.center.x)[split_axis] < middle;
                                }) -
                            prims.data());
            } else {
                // balanced tree split: find the largest axis of the bounding
                // box and split along this one right in the middle
                split_axis = largest_axis;
                mid = (start + end) / 2;
                std::nth_element(prims.data() + start, prims.data() + mid,
                    prims.data() + end, [split_axis](auto& a, auto& b) {
                        return (&a.center.x)[split_axis] <
                               (&b.center.x)[split_axis];
                    });
            }

            // if we were able to split, just break the primitives in half
            if (mid == start || mid == end) {
                split_axis = 0;
                mid = (start + end) / 2;
            }
        }

        // make an internal node
        node.type = bvh_node_type::internal;
        node.split_axis = split_axis;
        node.count = 2;
        node.prims[0] =
            make_bvh_node(nodes, prims, start, mid, type, equal_size);
        node.prims[1] = make_bvh_node(nodes, prims, mid, end, type, equal_size);
    } else {
        // Make a leaf node
        node.type = type;
        node.count = end - start;
        for (auto i = 0; i < node.count; i++)
            node.prims[i] = prims[start + i].primid;
    }

    // return nodeid
    return nodeid;
}

// Build a BVH from a set of primitives.
void build_bvh(const std::shared_ptr<bvh_tree>& bvh, bool equal_size) {
    // get the number of primitives and the primitive type
    auto prims = std::vector<bvh_prim>();
    auto type = bvh_node_type::internal;
    if (!bvh->points.empty()) {
        for (auto& p : bvh->points) {
            prims.push_back({point_bbox(bvh->pos[p], bvh->radius[p])});
        }
        type = bvh_node_type::point;
    } else if (!bvh->lines.empty()) {
        for (auto& l : bvh->lines) {
            prims.push_back({line_bbox(bvh->pos[l.x], bvh->pos[l.y],
                bvh->radius[l.x], bvh->radius[l.y])});
        }
        type = bvh_node_type::line;
    } else if (!bvh->triangles.empty()) {
        for (auto& t : bvh->triangles) {
            prims.push_back(
                {triangle_bbox(bvh->pos[t.x], bvh->pos[t.y], bvh->pos[t.z])});
        }
        type = bvh_node_type::triangle;
    } else if (!bvh->pos.empty()) {
        for (auto i = 0; i < bvh->pos.size(); i++) {
            prims.push_back({point_bbox(bvh->pos[i], bvh->radius[i])});
        }
        type = bvh_node_type::vertex;
    } else if (!bvh->ist_bvhs.empty()) {
        for (auto i = 0; i < bvh->ist_bvhs.size(); i++) {
            prims.push_back({transform_bbox(
                bvh->ist_frames[i], bvh->ist_bvhs[i]->nodes[0].bbox)});
        }
        type = bvh_node_type::instance;
    }

    // create an array of primitives to sort
    for (auto i = 0; i < prims.size(); i++) {
        prims[i].center = (prims[i].bbox.min + prims[i].bbox.max) / 2;
        prims[i].primid = i;
    }

    // build nodes
    bvh->nodes.clear();
    bvh->nodes.reserve(prims.size() * 2);
    make_bvh_node(bvh->nodes, prims, 0, (int)prims.size(), type, equal_size);
    bvh->nodes.shrink_to_fit();
}

// Recursively recomputes the node bounds for a shape bvh
void refit_bvh(const std::shared_ptr<bvh_tree>& bvh, int nodeid) {
    // refit
    auto& node = bvh->nodes[nodeid];
    node.bbox = invalid_bbox3f;
    switch (node.type) {
        case bvh_node_type::internal: {
            for (auto i = 0; i < 2; i++) {
                refit_bvh(bvh, node.prims[i]);
                node.bbox += bvh->nodes[node.prims[i]].bbox;
            }
        } break;
        case bvh_node_type::triangle: {
            for (auto i = 0; i < node.count; i++) {
                auto& t = bvh->triangles[node.prims[i]];
                node.bbox +=
                    triangle_bbox(bvh->pos[t.x], bvh->pos[t.y], bvh->pos[t.z]);
            }
        } break;
        case bvh_node_type::line: {
            for (auto i = 0; i < node.count; i++) {
                auto& l = bvh->lines[node.prims[i]];
                node.bbox += line_bbox(bvh->pos[l.x], bvh->pos[l.y],
                    bvh->radius[l.x], bvh->radius[l.y]);
            }
        } break;
        case bvh_node_type::point: {
            for (auto i = 0; i < node.count; i++) {
                auto& p = bvh->points[node.prims[i]];
                node.bbox += point_bbox(bvh->pos[p], bvh->radius[p]);
            }
        } break;
        case bvh_node_type::vertex: {
            for (auto i = 0; i < node.count; i++) {
                auto idx = node.prims[i];
                node.bbox += point_bbox(bvh->pos[idx], bvh->radius[idx]);
            }
        } break;
        case bvh_node_type::instance: {
            for (auto i = 0; i < node.count; i++) {
                auto idx = node.prims[i];
                node.bbox += transform_bbox(
                    bvh->ist_frames[idx], bvh->ist_bvhs[idx]->nodes[0].bbox);
            }
        } break;
    }
}

// Recursively recomputes the node bounds for a shape bvh
void refit_bvh(const std::shared_ptr<bvh_tree>& bvh) { refit_bvh(bvh, 0); }

// Intersect ray with a bvh.
bool intersect_bvh(const std::shared_ptr<bvh_tree>& bvh, const ray3f& ray_,
    bool find_any, float& dist, int& iid, int& eid, vec2f& uv) {
    // node stack
    int node_stack[128];
    auto node_cur = 0;
    node_stack[node_cur++] = 0;

    // shared variables
    auto hit = false;

    // copy ray to modify it
    auto ray = ray_;

    // prepare ray for fast queries
    auto ray_dinv = vec3f{1 / ray.d.x, 1 / ray.d.y, 1 / ray.d.z};
    auto ray_dsign = vec3i{(ray_dinv.x < 0) ? 1 : 0, (ray_dinv.y < 0) ? 1 : 0,
        (ray_dinv.z < 0) ? 1 : 0};

    // walking stack
    while (node_cur) {
        // grab node
        auto& node = bvh->nodes[node_stack[--node_cur]];

        // intersect bbox
        if (!intersect_bbox(ray, ray_dinv, ray_dsign, node.bbox)) continue;

        // intersect node, switching based on node type
        // for each type, iterate over the the primitive list
        switch (node.type) {
            case bvh_node_type::internal: {
                // for internal nodes, attempts to proceed along the
                // split axis from smallest to largest nodes
                if ((&ray_dsign.x)[node.split_axis]) {
                    node_stack[node_cur++] = node.prims[0];
                    node_stack[node_cur++] = node.prims[1];
                } else {
                    node_stack[node_cur++] = node.prims[1];
                    node_stack[node_cur++] = node.prims[0];
                }
            } break;
            case bvh_node_type::triangle: {
                for (auto i = 0; i < node.count; i++) {
                    auto& t = bvh->triangles[node.prims[i]];
                    if (intersect_triangle(ray, bvh->pos[t.x], bvh->pos[t.y],
                            bvh->pos[t.z], dist, uv)) {
                        hit = true;
                        ray.tmax = dist;
                        eid = node.prims[i];
                    }
                }
            } break;
            case bvh_node_type::line: {
                for (auto i = 0; i < node.count; i++) {
                    auto& l = bvh->lines[node.prims[i]];
                    if (intersect_line(ray, bvh->pos[l.x], bvh->pos[l.y],
                            bvh->radius[l.x], bvh->radius[l.y], dist, uv)) {
                        hit = true;
                        ray.tmax = dist;
                        eid = node.prims[i];
                    }
                }
            } break;
            case bvh_node_type::point: {
                for (auto i = 0; i < node.count; i++) {
                    auto& p = bvh->points[node.prims[i]];
                    if (intersect_point(
                            ray, bvh->pos[p], bvh->radius[p], dist, uv)) {
                        hit = true;
                        ray.tmax = dist;
                        eid = node.prims[i];
                    }
                }
            } break;
            case bvh_node_type::vertex: {
                for (auto i = 0; i < node.count; i++) {
                    auto idx = node.prims[i];
                    if (intersect_point(
                            ray, bvh->pos[idx], bvh->radius[idx], dist, uv)) {
                        hit = true;
                        ray.tmax = dist;
                        eid = node.prims[i];
                    }
                }
            } break;
            case bvh_node_type::instance: {
                for (auto i = 0; i < node.count; i++) {
                    auto idx = node.prims[i];
                    if (intersect_bvh(bvh->ist_bvhs[idx],
                            transform_ray(bvh->ist_inv_frames[idx], ray),
                            find_any, dist, iid, eid, uv)) {
                        hit = true;
                        ray.tmax = dist;
                        iid = node.prims[i];
                    }
                }
            } break;
        }

        // check for early exit
        if (find_any && hit) return true;
    }

    return hit;
}

// Finds the closest element with a bvh.
bool overlap_bvh(const std::shared_ptr<bvh_tree>& bvh, const vec3f& pos,
    float max_dist, bool find_any, float& dist, int& iid, int& eid, vec2f& uv) {
    // node stack
    int node_stack[64];
    auto node_cur = 0;
    node_stack[node_cur++] = 0;

    // hit
    auto hit = false;

    // walking stack
    while (node_cur) {
        // grab node
        auto node = bvh->nodes[node_stack[--node_cur]];

        // intersect bbox
        if (!distance_check_bbox(pos, max_dist, node.bbox)) continue;

        // intersect node, switching based on node type
        // for each type, iterate over the the primitive list
        switch (node.type) {
            case bvh_node_type::internal: {
                // internal node
                node_stack[node_cur++] = node.prims[0];
                node_stack[node_cur++] = node.prims[1];
            } break;
            case bvh_node_type::triangle: {
                for (auto i = 0; i < node.count; i++) {
                    auto& t = bvh->triangles[node.prims[i]];
                    if (overlap_triangle(pos, max_dist, bvh->pos[t.x],
                            bvh->pos[t.y], bvh->pos[t.z], bvh->radius[t.x],
                            bvh->radius[t.y], bvh->radius[t.z], dist, uv)) {
                        hit = true;
                        max_dist = dist;
                        eid = node.prims[i];
                    }
                }
            } break;
            case bvh_node_type::line: {
                for (auto i = 0; i < node.count; i++) {
                    auto& l = bvh->lines[node.prims[i]];
                    if (overlap_line(pos, max_dist, bvh->pos[l.x],
                            bvh->pos[l.y], bvh->radius[l.x], bvh->radius[l.y],
                            dist, uv)) {
                        hit = true;
                        max_dist = dist;
                        eid = node.prims[i];
                    }
                }
            } break;
            case bvh_node_type::point: {
                for (auto i = 0; i < node.count; i++) {
                    auto& p = bvh->points[node.prims[i]];
                    if (overlap_point(pos, max_dist, bvh->pos[p],
                            bvh->radius[p], dist, uv)) {
                        hit = true;
                        max_dist = dist;
                        eid = node.prims[i];
                    }
                }
            } break;
            case bvh_node_type::vertex: {
                for (auto i = 0; i < node.count; i++) {
                    auto idx = node.prims[i];
                    if (overlap_point(pos, max_dist, bvh->pos[idx],
                            bvh->radius[idx], dist, uv)) {
                        hit = true;
                        max_dist = dist;
                        eid = node.prims[i];
                    }
                }
            } break;
            case bvh_node_type::instance: {
                for (auto i = 0; i < node.count; i++) {
                    auto idx = node.prims[i];
                    if (overlap_bvh(bvh->ist_bvhs[idx],
                            transform_point(bvh->ist_inv_frames[idx], pos),
                            max_dist, find_any, dist, iid, eid, uv)) {
                        hit = true;
                        max_dist = dist;
                        iid = node.prims[i];
                    }
                }
            } break;
        }

        // check for early exit
        if (find_any && hit) return true;
    }

    return hit;
}

#if 0
    // Finds the overlap between BVH leaf nodes.
    template <typename OverlapElem>
    void overlap_bvh_elems(const std::shared_ptr<bvh_tree>& bvh1, const std::shared_ptr<bvh_tree>& bvh2,
                           bool skip_duplicates, bool skip_self, std::vector<vec2i>& overlaps,
                           const OverlapElem& overlap_elems) {
        // node stack
        vec2i node_stack[128];
        auto node_cur = 0;
        node_stack[node_cur++] = {0, 0};

        // walking stack
        while (node_cur) {
            // grab node
            auto node_idx = node_stack[--node_cur];
            const auto node1 = bvh1->nodes[node_idx.x];
            const auto node2 = bvh2->nodes[node_idx.y];

            // intersect bbox
            if (!overlap_bbox(node1.bbox, node2.bbox)) continue;

            // check for leaves
            if (node1.isleaf && node2.isleaf) {
                // collide primitives
                for (auto i1 = node1.start; i1 < node1.start + node1.count; i1++) {
                    for (auto i2 = node2.start; i2 < node2.start + node2.count;
                         i2++) {
                        auto idx1 = bvh1->sorted_prim[i1];
                        auto idx2 = bvh2->sorted_prim[i2];
                        if (skip_duplicates && idx1 > idx2) continue;
                        if (skip_self && idx1 == idx2) continue;
                        if (overlap_elems(idx1, idx2))
                            overlaps.push_back({idx1, idx2});
                    }
                }
            } else {
                // descend
                if (node1.isleaf) {
                    for (auto idx2 = node2.start; idx2 < node2.start + node2.count;
                         idx2++) {
                        node_stack[node_cur++] = {node_idx.x, (int)idx2};
                    }
                } else if (node2.isleaf) {
                    for (auto idx1 = node1.start; idx1 < node1.start + node1.count;
                         idx1++) {
                        node_stack[node_cur++] = {(int)idx1, node_idx.y};
                    }
                } else {
                    for (auto idx2 = node2.start; idx2 < node2.start + node2.count;
                         idx2++) {
                        for (auto idx1 = node1.start;
                             idx1 < node1.start + node1.count; idx1++) {
                            node_stack[node_cur++] = {(int)idx1, (int)idx2};
                        }
                    }
                }
            }
        }
    }
#endif

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHAPE EXAMPLES
// -----------------------------------------------------------------------------
namespace ygl {

// Make a quad.
void make_quad(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec2i& steps,
    const vec2f& size, const vec2f& uvsize) {
    auto vid = [steps](int i, int j) { return j * (steps.x + 1) + i; };

    pos.resize((steps.x + 1) * (steps.y + 1));
    norm.resize((steps.x + 1) * (steps.y + 1));
    texcoord.resize((steps.x + 1) * (steps.y + 1));
    for (auto j = 0; j <= steps.y; j++) {
        for (auto i = 0; i <= steps.x; i++) {
            auto uv = vec2f{i / (float)steps.x, j / (float)steps.y};
            pos[vid(i, j)] = {
                (uv.x - 0.5f) * size.x, (uv.y - 0.5f) * size.y, 0};
            norm[vid(i, j)] = {0, 0, 1};
            texcoord[vid(i, j)] = uv * uvsize;
        }
    }

    quads.resize(steps.x * steps.y);
    for (auto j = 0; j < steps.y; j++) {
        for (auto i = 0; i < steps.x; i++) {
            quads[j * steps.x + i] = {
                vid(i, j), vid(i + 1, j), vid(i + 1, j + 1), vid(i, j + 1)};
        }
    }
}

// Make a stack of quads
void make_quad_stack(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec3i& steps,
    const vec3f& size, const vec2f& uvsize) {
    std::vector<vec4i> qquads;
    std::vector<vec3f> qpos;
    std::vector<vec3f> qnorm;
    std::vector<vec2f> qtexcoord;

    quads.clear();
    pos.clear();
    norm.clear();
    texcoord.clear();

    for (auto i = 0; i <= steps.z; i++) {
        make_quad(qquads, qpos, qnorm, qtexcoord, {steps.x, steps.y},
            {size.x, size.y}, uvsize);
        for (auto& p : qpos) p.z = (-0.5f + (float)i / steps.z) * size.z;
        merge_quads(quads, pos, norm, texcoord, qquads, qpos, qnorm, qtexcoord);
    }
}

// Make a cube.
void make_cube(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec3i& steps,
    const vec3f& size, const vec3f& uvsize) {
    std::vector<vec4i> qquads;
    std::vector<vec3f> qpos;
    std::vector<vec3f> qnorm;
    std::vector<vec2f> qtexcoord;

    quads.clear();
    pos.clear();
    norm.clear();
    texcoord.clear();

    // +z
    make_quad(qquads, qpos, qnorm, qtexcoord, {steps.x, steps.y},
        {size.x, size.y}, {uvsize.x, uvsize.y});
    for (auto i = 0; i < qpos.size(); i++) {
        qpos[i] = {qpos[i].x, qpos[i].y, size.z / 2};
        qnorm[i] = {0, 0, 1};
    }
    merge_quads(quads, pos, norm, texcoord, qquads, qpos, qnorm, qtexcoord);
    // -z
    make_quad(qquads, qpos, qnorm, qtexcoord, {steps.x, steps.y},
        {size.x, size.y}, {uvsize.x, uvsize.y});
    for (auto i = 0; i < qpos.size(); i++) {
        qpos[i] = {-qpos[i].x, qpos[i].y, -size.z / 2};
        qnorm[i] = {0, 0, -1};
    }
    merge_quads(quads, pos, norm, texcoord, qquads, qpos, qnorm, qtexcoord);
    // +x
    make_quad(qquads, qpos, qnorm, qtexcoord, {steps.z, steps.y},
        {size.z, size.y}, {uvsize.z, uvsize.y});
    for (auto i = 0; i < qpos.size(); i++) {
        qpos[i] = {size.x / 2, qpos[i].y, -qpos[i].x};
        qnorm[i] = {1, 0, 0};
    }
    merge_quads(quads, pos, norm, texcoord, qquads, qpos, qnorm, qtexcoord);
    // -x
    make_quad(qquads, qpos, qnorm, qtexcoord, {steps.z, steps.y},
        {size.z, size.y}, {uvsize.z, uvsize.y});
    for (auto i = 0; i < qpos.size(); i++) {
        qpos[i] = {-size.x / 2, qpos[i].y, qpos[i].x};
        qnorm[i] = {-1, 0, 0};
    }
    merge_quads(quads, pos, norm, texcoord, qquads, qpos, qnorm, qtexcoord);
    // +y
    make_quad(qquads, qpos, qnorm, qtexcoord, {steps.x, steps.y},
        {size.x, size.y}, {uvsize.x, uvsize.y});
    for (auto i = 0; i < qpos.size(); i++) {
        qpos[i] = {qpos[i].x, size.y / 2, -qpos[i].y};
        qnorm[i] = {0, 1, 0};
    }
    merge_quads(quads, pos, norm, texcoord, qquads, qpos, qnorm, qtexcoord);
    // +y
    make_quad(qquads, qpos, qnorm, qtexcoord, {steps.x, steps.y},
        {size.x, size.y}, {uvsize.x, uvsize.y});
    for (auto i = 0; i < qpos.size(); i++) {
        qpos[i] = {qpos[i].x, -size.y / 2, qpos[i].y};
        qnorm[i] = {0, -1, 0};
    }
    merge_quads(quads, pos, norm, texcoord, qquads, qpos, qnorm, qtexcoord);
}

// Make a rounded cube.
void make_cube_rounded(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec3i& steps,
    const vec3f& size, const vec3f& uvsize, float radius) {
    make_cube(quads, pos, norm, texcoord, steps, size, uvsize);
    auto c = size / 2 - vec3f{radius, radius, radius};
    for (auto i = 0; i < pos.size(); i++) {
        auto pc = vec3f{fabs(pos[i].x), fabs(pos[i].y), fabs(pos[i].z)};
        auto ps = vec3f{pos[i].x < 0 ? -1.0f : 1.0f,
            pos[i].y < 0 ? -1.0f : 1.0f, pos[i].z < 0 ? -1.0f : 1.0f};
        if (pc.x >= c.x && pc.y >= c.y && pc.z >= c.z) {
            auto pn = normalize(pc - c);
            pos[i] = c + radius * pn;
            norm[i] = pn;
        } else if (pc.x >= c.x && pc.y >= c.y) {
            auto pn = normalize((pc - c) * vec3f{1, 1, 0});
            pos[i] = {c.x + radius * pn.x, c.y + radius * pn.y, pc.z};
            norm[i] = pn;
        } else if (pc.x >= c.x && pc.z >= c.z) {
            auto pn = normalize((pc - c) * vec3f{1, 0, 1});
            pos[i] = {c.x + radius * pn.x, pc.y, c.z + radius * pn.z};
            norm[i] = pn;
        } else if (pc.y >= c.y && pc.z >= c.z) {
            auto pn = normalize((pc - c) * vec3f{0, 1, 1});
            pos[i] = {pc.x, c.y + radius * pn.y, c.z + radius * pn.z};
            norm[i] = pn;
        } else {
            continue;
        }
        pos[i] *= ps;
        norm[i] *= ps;
    }
}

// Make a sphere.
void make_sphere(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec2i& steps,
    float size, const vec2f& uvsize) {
    make_quad(quads, pos, norm, texcoord, steps, {1, 1}, {1, 1});
    for (auto i = 0; i < pos.size(); i++) {
        auto uv = texcoord[i];
        auto a = vec2f{2 * pi * uv.x, pi * (1 - uv.y)};
        auto p = vec3f{cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
        pos[i] = p * (size / 2);
        norm[i] = normalize(p);
        texcoord[i] = uv * uvsize;
    }
}

// Make a spherecube.
void make_sphere_cube(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, int steps,
    float size, float uvsize) {
    make_cube(quads, pos, norm, texcoord, {steps, steps, steps}, {1, 1, 1},
        {uvsize, uvsize, uvsize});
    for (auto i = 0; i < pos.size(); i++) {
        auto p = pos[i];
        pos[i] = normalize(p) * (size / 2);
        norm[i] = normalize(p);
    }
}

// Make a flipped sphere. This is not watertight.
void make_sphere_flipcap(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec2i& steps,
    float size, const vec2f& uvsize, const vec2f& zflip) {
    make_sphere(quads, pos, norm, texcoord, steps, size, uvsize);
    for (auto i = 0; i < pos.size(); i++) {
        if (pos[i].z > zflip.y) {
            pos[i].z = 2 * zflip.y - pos[i].z;
            norm[i].x = -norm[i].x;
            norm[i].y = -norm[i].y;
        } else if (pos[i].z < zflip.x) {
            pos[i].z = 2 * zflip.x - pos[i].z;
            norm[i].x = -norm[i].x;
            norm[i].y = -norm[i].y;
        }
    }
}

// Make a disk.
void make_disk(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec2i& steps,
    float size, const vec2f& uvsize) {
    make_quad(quads, pos, norm, texcoord, steps, {1, 1}, {1, 1});
    for (auto i = 0; i < pos.size(); i++) {
        auto uv = texcoord[i];
        auto phi = 2 * pi * uv.x;
        pos[i] = {cos(phi) * uv.y * size / 2, sin(phi) * uv.y * size / 2, 0};
        norm[i] = {0, 0, 1};
        texcoord[i] = uv * uvsize;
    }
}

// Make a disk from a quad.
void make_disk_quad(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, int steps,
    float size, float uvsize) {
    make_quad(
        quads, pos, norm, texcoord, {steps, steps}, {2, 2}, {uvsize, uvsize});
    for (auto i = 0; i < pos.size(); i++) {
        // Analytical Methods for Squaring the Disc, by C. Fong
        // https://arxiv.org/abs/1509.06344
        auto xy = vec2f{pos[i].x, pos[i].y};
        auto uv = vec2f{
            xy.x * sqrt(1 - xy.y * xy.y / 2), xy.y * sqrt(1 - xy.x * xy.x / 2)};
        pos[i] = {uv.x * size / 2, uv.y * size / 2, 0};
    }
}

// Make a bulged disk from a quad.
void make_disk_bulged(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, int steps,
    float size, float uvsize, float height) {
    make_disk_quad(quads, pos, norm, texcoord, steps, size, uvsize);
    if (height == 0) return;
    auto radius = (size * size / 4 + height * height) / (2 * height);
    auto center = vec3f{0, 0, -radius + height};
    for (auto i = 0; i < pos.size(); i++) {
        auto pn = normalize(pos[i] - center);
        pos[i] = center + pn * radius;
        norm[i] = pn;
    }
}

// Make a cylinder (side-only).
void make_cylinder_side(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec2i& steps,
    const vec2f& size, const vec2f& uvsize) {
    make_quad(quads, pos, norm, texcoord, steps, {1, 1}, {1, 1});
    for (auto i = 0; i < pos.size(); i++) {
        auto uv = texcoord[i];
        auto phi = 2 * pi * uv.x;
        pos[i] = {cos(phi) * size.x / 2, sin(phi) * size.x / 2,
            (uv.y - 0.5f) * size.y};
        norm[i] = {cos(phi), sin(phi), 0};
        texcoord[i] = uv * uvsize;
    }
}

// Make a cylinder.
void make_cylinder(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec3i& steps,
    const vec2f& size, const vec3f& uvsize) {
    std::vector<vec4i> qquads;
    std::vector<vec3f> qpos;
    std::vector<vec3f> qnorm;
    std::vector<vec2f> qtexcoord;

    quads.clear();
    pos.clear();
    norm.clear();
    texcoord.clear();

    // side
    make_cylinder_side(qquads, qpos, qnorm, qtexcoord, {steps.x, steps.y},
        {size.x, size.y}, {uvsize.x, uvsize.y});
    merge_quads(quads, pos, norm, texcoord, qquads, qpos, qnorm, qtexcoord);
    // top
    make_disk(qquads, qpos, qnorm, qtexcoord, {steps.x, steps.z}, size.x,
        {uvsize.x, uvsize.z});
    for (auto i = 0; i < qpos.size(); i++) { qpos[i].z = size.y / 2; }
    merge_quads(quads, pos, norm, texcoord, qquads, qpos, qnorm, qtexcoord);
    // bottom
    make_disk(qquads, qpos, qnorm, qtexcoord, {steps.x, steps.z}, size.x,
        {uvsize.x, uvsize.z});
    for (auto i = 0; i < qpos.size(); i++) {
        qpos[i].z = -size.y / 2;
        qnorm[i] = -qnorm[i];
    }
    for (auto i = 0; i < qquads.size(); i++)
        std::swap(qquads[i].x, qquads[i].z);
    merge_quads(quads, pos, norm, texcoord, qquads, qpos, qnorm, qtexcoord);
}

// Make a rounded cylinder.
void make_cylinder_rounded(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec3i& steps,
    const vec2f& size, const vec3f& uvsize, float radius) {
    auto c = size / 2 - vec2f{radius, radius};
    for (auto i = 0; i < pos.size(); i++) {
        auto phi = atan2(pos[i].y, pos[i].x);
        auto r = length(vec2f{pos[i].x, pos[i].y});
        auto z = pos[i].z;
        auto pc = vec2f{r, fabs(z)};
        auto ps = (z < 0) ? -1.0f : 1.0f;
        if (pc.x >= c.x && pc.y >= c.y) {
            auto pn = normalize(pc - c);
            pos[i] = {cos(phi) * c.x + radius * pn.x,
                sin(phi) * c.x + radius * pn.x, ps * (c.y + radius * pn.y)};
            norm[i] = {cos(phi) * pn.x, sin(phi) * pn.x, ps * pn.y};
        } else {
            continue;
        }
    }
}

// Make a geodesic sphere.
void make_geodesic_sphere(std::vector<vec3i>& triangles,
    std::vector<vec3f>& pos, int tesselation, float size) {
    // https://stackoverflow.com/questions/17705621/algorithm-for-a-geodesic-sphere
    const float X = 0.525731112119133606f;
    const float Z = 0.850650808352039932f;
    pos = std::vector<vec3f>{{-X, 0.0, Z}, {X, 0.0, Z}, {-X, 0.0, -Z},
        {X, 0.0, -Z}, {0.0, Z, X}, {0.0, Z, -X}, {0.0, -Z, X}, {0.0, -Z, -X},
        {Z, X, 0.0}, {-Z, X, 0.0}, {Z, -X, 0.0}, {-Z, -X, 0.0}};
    triangles = std::vector<vec3i>{{0, 1, 4}, {0, 4, 9}, {9, 4, 5}, {4, 8, 5},
        {4, 1, 8}, {8, 1, 10}, {8, 10, 3}, {5, 8, 3}, {5, 3, 2}, {2, 3, 7},
        {7, 3, 10}, {7, 10, 6}, {7, 6, 11}, {11, 6, 0}, {0, 6, 1}, {6, 10, 1},
        {9, 11, 0}, {9, 2, 11}, {9, 5, 2}, {7, 11, 2}};
    subdivide_triangles(triangles, pos, max(0, tesselation - 2));
    for (auto& p : pos) p = normalize(p) * size / 2;
}

// Make a facevarying cube with unique vertices but different texture
// coordinates.
void make_fvcube(std::vector<vec4i>& quads_pos, std::vector<vec3f>& pos,
    std::vector<vec4i>& quads_norm, std::vector<vec3f>& norm,
    std::vector<vec4i>& quads_texcoord, std::vector<vec2f>& texcoord,
    vec3i steps, const vec3f& size, const vec3f& uvsize) {
    make_cube(quads_pos, pos, norm, texcoord, steps, size, uvsize);
    quads_norm = quads_pos;
    quads_texcoord = quads_pos;
    weld_quads(quads_pos, pos,
        min(0.1f * size /
            vec3f{(float)steps.x, (float)steps.y, (float)steps.z}));
}

// Make a suzanne monkey model for testing. Note that some quads are
// degenerate.
void make_suzanne(std::vector<vec4i>& quads, std::vector<vec3f>& pos) {
    static auto suzanne_pos = std::vector<vec3f>{{0.4375, 0.1640625, 0.765625},
        {-0.4375, 0.1640625, 0.765625}, {0.5, 0.09375, 0.6875},
        {-0.5, 0.09375, 0.6875}, {0.546875, 0.0546875, 0.578125},
        {-0.546875, 0.0546875, 0.578125}, {0.3515625, -0.0234375, 0.6171875},
        {-0.3515625, -0.0234375, 0.6171875}, {0.3515625, 0.03125, 0.71875},
        {-0.3515625, 0.03125, 0.71875}, {0.3515625, 0.1328125, 0.78125},
        {-0.3515625, 0.1328125, 0.78125}, {0.2734375, 0.1640625, 0.796875},
        {-0.2734375, 0.1640625, 0.796875}, {0.203125, 0.09375, 0.7421875},
        {-0.203125, 0.09375, 0.7421875}, {0.15625, 0.0546875, 0.6484375},
        {-0.15625, 0.0546875, 0.6484375}, {0.078125, 0.2421875, 0.65625},
        {-0.078125, 0.2421875, 0.65625}, {0.140625, 0.2421875, 0.7421875},
        {-0.140625, 0.2421875, 0.7421875}, {0.2421875, 0.2421875, 0.796875},
        {-0.2421875, 0.2421875, 0.796875}, {0.2734375, 0.328125, 0.796875},
        {-0.2734375, 0.328125, 0.796875}, {0.203125, 0.390625, 0.7421875},
        {-0.203125, 0.390625, 0.7421875}, {0.15625, 0.4375, 0.6484375},
        {-0.15625, 0.4375, 0.6484375}, {0.3515625, 0.515625, 0.6171875},
        {-0.3515625, 0.515625, 0.6171875}, {0.3515625, 0.453125, 0.71875},
        {-0.3515625, 0.453125, 0.71875}, {0.3515625, 0.359375, 0.78125},
        {-0.3515625, 0.359375, 0.78125}, {0.4375, 0.328125, 0.765625},
        {-0.4375, 0.328125, 0.765625}, {0.5, 0.390625, 0.6875},
        {-0.5, 0.390625, 0.6875}, {0.546875, 0.4375, 0.578125},
        {-0.546875, 0.4375, 0.578125}, {0.625, 0.2421875, 0.5625},
        {-0.625, 0.2421875, 0.5625}, {0.5625, 0.2421875, 0.671875},
        {-0.5625, 0.2421875, 0.671875}, {0.46875, 0.2421875, 0.7578125},
        {-0.46875, 0.2421875, 0.7578125}, {0.4765625, 0.2421875, 0.7734375},
        {-0.4765625, 0.2421875, 0.7734375}, {0.4453125, 0.3359375, 0.78125},
        {-0.4453125, 0.3359375, 0.78125}, {0.3515625, 0.375, 0.8046875},
        {-0.3515625, 0.375, 0.8046875}, {0.265625, 0.3359375, 0.8203125},
        {-0.265625, 0.3359375, 0.8203125}, {0.2265625, 0.2421875, 0.8203125},
        {-0.2265625, 0.2421875, 0.8203125}, {0.265625, 0.15625, 0.8203125},
        {-0.265625, 0.15625, 0.8203125}, {0.3515625, 0.2421875, 0.828125},
        {-0.3515625, 0.2421875, 0.828125}, {0.3515625, 0.1171875, 0.8046875},
        {-0.3515625, 0.1171875, 0.8046875}, {0.4453125, 0.15625, 0.78125},
        {-0.4453125, 0.15625, 0.78125}, {0.0, 0.4296875, 0.7421875},
        {0.0, 0.3515625, 0.8203125}, {0.0, -0.6796875, 0.734375},
        {0.0, -0.3203125, 0.78125}, {0.0, -0.1875, 0.796875},
        {0.0, -0.7734375, 0.71875}, {0.0, 0.40625, 0.6015625},
        {0.0, 0.5703125, 0.5703125}, {0.0, 0.8984375, -0.546875},
        {0.0, 0.5625, -0.8515625}, {0.0, 0.0703125, -0.828125},
        {0.0, -0.3828125, -0.3515625}, {0.203125, -0.1875, 0.5625},
        {-0.203125, -0.1875, 0.5625}, {0.3125, -0.4375, 0.5703125},
        {-0.3125, -0.4375, 0.5703125}, {0.3515625, -0.6953125, 0.5703125},
        {-0.3515625, -0.6953125, 0.5703125}, {0.3671875, -0.890625, 0.53125},
        {-0.3671875, -0.890625, 0.53125}, {0.328125, -0.9453125, 0.5234375},
        {-0.328125, -0.9453125, 0.5234375}, {0.1796875, -0.96875, 0.5546875},
        {-0.1796875, -0.96875, 0.5546875}, {0.0, -0.984375, 0.578125},
        {0.4375, -0.140625, 0.53125}, {-0.4375, -0.140625, 0.53125},
        {0.6328125, -0.0390625, 0.5390625}, {-0.6328125, -0.0390625, 0.5390625},
        {0.828125, 0.1484375, 0.4453125}, {-0.828125, 0.1484375, 0.4453125},
        {0.859375, 0.4296875, 0.59375}, {-0.859375, 0.4296875, 0.59375},
        {0.7109375, 0.484375, 0.625}, {-0.7109375, 0.484375, 0.625},
        {0.4921875, 0.6015625, 0.6875}, {-0.4921875, 0.6015625, 0.6875},
        {0.3203125, 0.7578125, 0.734375}, {-0.3203125, 0.7578125, 0.734375},
        {0.15625, 0.71875, 0.7578125}, {-0.15625, 0.71875, 0.7578125},
        {0.0625, 0.4921875, 0.75}, {-0.0625, 0.4921875, 0.75},
        {0.1640625, 0.4140625, 0.7734375}, {-0.1640625, 0.4140625, 0.7734375},
        {0.125, 0.3046875, 0.765625}, {-0.125, 0.3046875, 0.765625},
        {0.203125, 0.09375, 0.7421875}, {-0.203125, 0.09375, 0.7421875},
        {0.375, 0.015625, 0.703125}, {-0.375, 0.015625, 0.703125},
        {0.4921875, 0.0625, 0.671875}, {-0.4921875, 0.0625, 0.671875},
        {0.625, 0.1875, 0.6484375}, {-0.625, 0.1875, 0.6484375},
        {0.640625, 0.296875, 0.6484375}, {-0.640625, 0.296875, 0.6484375},
        {0.6015625, 0.375, 0.6640625}, {-0.6015625, 0.375, 0.6640625},
        {0.4296875, 0.4375, 0.71875}, {-0.4296875, 0.4375, 0.71875},
        {0.25, 0.46875, 0.7578125}, {-0.25, 0.46875, 0.7578125},
        {0.0, -0.765625, 0.734375}, {0.109375, -0.71875, 0.734375},
        {-0.109375, -0.71875, 0.734375}, {0.1171875, -0.8359375, 0.7109375},
        {-0.1171875, -0.8359375, 0.7109375}, {0.0625, -0.8828125, 0.6953125},
        {-0.0625, -0.8828125, 0.6953125}, {0.0, -0.890625, 0.6875},
        {0.0, -0.1953125, 0.75}, {0.0, -0.140625, 0.7421875},
        {0.1015625, -0.1484375, 0.7421875}, {-0.1015625, -0.1484375, 0.7421875},
        {0.125, -0.2265625, 0.75}, {-0.125, -0.2265625, 0.75},
        {0.0859375, -0.2890625, 0.7421875}, {-0.0859375, -0.2890625, 0.7421875},
        {0.3984375, -0.046875, 0.671875}, {-0.3984375, -0.046875, 0.671875},
        {0.6171875, 0.0546875, 0.625}, {-0.6171875, 0.0546875, 0.625},
        {0.7265625, 0.203125, 0.6015625}, {-0.7265625, 0.203125, 0.6015625},
        {0.7421875, 0.375, 0.65625}, {-0.7421875, 0.375, 0.65625},
        {0.6875, 0.4140625, 0.7265625}, {-0.6875, 0.4140625, 0.7265625},
        {0.4375, 0.546875, 0.796875}, {-0.4375, 0.546875, 0.796875},
        {0.3125, 0.640625, 0.8359375}, {-0.3125, 0.640625, 0.8359375},
        {0.203125, 0.6171875, 0.8515625}, {-0.203125, 0.6171875, 0.8515625},
        {0.1015625, 0.4296875, 0.84375}, {-0.1015625, 0.4296875, 0.84375},
        {0.125, -0.1015625, 0.8125}, {-0.125, -0.1015625, 0.8125},
        {0.2109375, -0.4453125, 0.7109375}, {-0.2109375, -0.4453125, 0.7109375},
        {0.25, -0.703125, 0.6875}, {-0.25, -0.703125, 0.6875},
        {0.265625, -0.8203125, 0.6640625}, {-0.265625, -0.8203125, 0.6640625},
        {0.234375, -0.9140625, 0.6328125}, {-0.234375, -0.9140625, 0.6328125},
        {0.1640625, -0.9296875, 0.6328125}, {-0.1640625, -0.9296875, 0.6328125},
        {0.0, -0.9453125, 0.640625}, {0.0, 0.046875, 0.7265625},
        {0.0, 0.2109375, 0.765625}, {0.328125, 0.4765625, 0.7421875},
        {-0.328125, 0.4765625, 0.7421875}, {0.1640625, 0.140625, 0.75},
        {-0.1640625, 0.140625, 0.75}, {0.1328125, 0.2109375, 0.7578125},
        {-0.1328125, 0.2109375, 0.7578125}, {0.1171875, -0.6875, 0.734375},
        {-0.1171875, -0.6875, 0.734375}, {0.078125, -0.4453125, 0.75},
        {-0.078125, -0.4453125, 0.75}, {0.0, -0.4453125, 0.75},
        {0.0, -0.328125, 0.7421875}, {0.09375, -0.2734375, 0.78125},
        {-0.09375, -0.2734375, 0.78125}, {0.1328125, -0.2265625, 0.796875},
        {-0.1328125, -0.2265625, 0.796875}, {0.109375, -0.1328125, 0.78125},
        {-0.109375, -0.1328125, 0.78125}, {0.0390625, -0.125, 0.78125},
        {-0.0390625, -0.125, 0.78125}, {0.0, -0.203125, 0.828125},
        {0.046875, -0.1484375, 0.8125}, {-0.046875, -0.1484375, 0.8125},
        {0.09375, -0.15625, 0.8125}, {-0.09375, -0.15625, 0.8125},
        {0.109375, -0.2265625, 0.828125}, {-0.109375, -0.2265625, 0.828125},
        {0.078125, -0.25, 0.8046875}, {-0.078125, -0.25, 0.8046875},
        {0.0, -0.2890625, 0.8046875}, {0.2578125, -0.3125, 0.5546875},
        {-0.2578125, -0.3125, 0.5546875}, {0.1640625, -0.2421875, 0.7109375},
        {-0.1640625, -0.2421875, 0.7109375}, {0.1796875, -0.3125, 0.7109375},
        {-0.1796875, -0.3125, 0.7109375}, {0.234375, -0.25, 0.5546875},
        {-0.234375, -0.25, 0.5546875}, {0.0, -0.875, 0.6875},
        {0.046875, -0.8671875, 0.6875}, {-0.046875, -0.8671875, 0.6875},
        {0.09375, -0.8203125, 0.7109375}, {-0.09375, -0.8203125, 0.7109375},
        {0.09375, -0.7421875, 0.7265625}, {-0.09375, -0.7421875, 0.7265625},
        {0.0, -0.78125, 0.65625}, {0.09375, -0.75, 0.6640625},
        {-0.09375, -0.75, 0.6640625}, {0.09375, -0.8125, 0.640625},
        {-0.09375, -0.8125, 0.640625}, {0.046875, -0.8515625, 0.6328125},
        {-0.046875, -0.8515625, 0.6328125}, {0.0, -0.859375, 0.6328125},
        {0.171875, 0.21875, 0.78125}, {-0.171875, 0.21875, 0.78125},
        {0.1875, 0.15625, 0.7734375}, {-0.1875, 0.15625, 0.7734375},
        {0.3359375, 0.4296875, 0.7578125}, {-0.3359375, 0.4296875, 0.7578125},
        {0.2734375, 0.421875, 0.7734375}, {-0.2734375, 0.421875, 0.7734375},
        {0.421875, 0.3984375, 0.7734375}, {-0.421875, 0.3984375, 0.7734375},
        {0.5625, 0.3515625, 0.6953125}, {-0.5625, 0.3515625, 0.6953125},
        {0.5859375, 0.2890625, 0.6875}, {-0.5859375, 0.2890625, 0.6875},
        {0.578125, 0.1953125, 0.6796875}, {-0.578125, 0.1953125, 0.6796875},
        {0.4765625, 0.1015625, 0.71875}, {-0.4765625, 0.1015625, 0.71875},
        {0.375, 0.0625, 0.7421875}, {-0.375, 0.0625, 0.7421875},
        {0.2265625, 0.109375, 0.78125}, {-0.2265625, 0.109375, 0.78125},
        {0.1796875, 0.296875, 0.78125}, {-0.1796875, 0.296875, 0.78125},
        {0.2109375, 0.375, 0.78125}, {-0.2109375, 0.375, 0.78125},
        {0.234375, 0.359375, 0.7578125}, {-0.234375, 0.359375, 0.7578125},
        {0.1953125, 0.296875, 0.7578125}, {-0.1953125, 0.296875, 0.7578125},
        {0.2421875, 0.125, 0.7578125}, {-0.2421875, 0.125, 0.7578125},
        {0.375, 0.0859375, 0.7265625}, {-0.375, 0.0859375, 0.7265625},
        {0.4609375, 0.1171875, 0.703125}, {-0.4609375, 0.1171875, 0.703125},
        {0.546875, 0.2109375, 0.671875}, {-0.546875, 0.2109375, 0.671875},
        {0.5546875, 0.28125, 0.671875}, {-0.5546875, 0.28125, 0.671875},
        {0.53125, 0.3359375, 0.6796875}, {-0.53125, 0.3359375, 0.6796875},
        {0.4140625, 0.390625, 0.75}, {-0.4140625, 0.390625, 0.75},
        {0.28125, 0.3984375, 0.765625}, {-0.28125, 0.3984375, 0.765625},
        {0.3359375, 0.40625, 0.75}, {-0.3359375, 0.40625, 0.75},
        {0.203125, 0.171875, 0.75}, {-0.203125, 0.171875, 0.75},
        {0.1953125, 0.2265625, 0.75}, {-0.1953125, 0.2265625, 0.75},
        {0.109375, 0.4609375, 0.609375}, {-0.109375, 0.4609375, 0.609375},
        {0.1953125, 0.6640625, 0.6171875}, {-0.1953125, 0.6640625, 0.6171875},
        {0.3359375, 0.6875, 0.59375}, {-0.3359375, 0.6875, 0.59375},
        {0.484375, 0.5546875, 0.5546875}, {-0.484375, 0.5546875, 0.5546875},
        {0.6796875, 0.453125, 0.4921875}, {-0.6796875, 0.453125, 0.4921875},
        {0.796875, 0.40625, 0.4609375}, {-0.796875, 0.40625, 0.4609375},
        {0.7734375, 0.1640625, 0.375}, {-0.7734375, 0.1640625, 0.375},
        {0.6015625, 0.0, 0.4140625}, {-0.6015625, 0.0, 0.4140625},
        {0.4375, -0.09375, 0.46875}, {-0.4375, -0.09375, 0.46875},
        {0.0, 0.8984375, 0.2890625}, {0.0, 0.984375, -0.078125},
        {0.0, -0.1953125, -0.671875}, {0.0, -0.4609375, 0.1875},
        {0.0, -0.9765625, 0.4609375}, {0.0, -0.8046875, 0.34375},
        {0.0, -0.5703125, 0.3203125}, {0.0, -0.484375, 0.28125},
        {0.8515625, 0.234375, 0.0546875}, {-0.8515625, 0.234375, 0.0546875},
        {0.859375, 0.3203125, -0.046875}, {-0.859375, 0.3203125, -0.046875},
        {0.7734375, 0.265625, -0.4375}, {-0.7734375, 0.265625, -0.4375},
        {0.4609375, 0.4375, -0.703125}, {-0.4609375, 0.4375, -0.703125},
        {0.734375, -0.046875, 0.0703125}, {-0.734375, -0.046875, 0.0703125},
        {0.59375, -0.125, -0.1640625}, {-0.59375, -0.125, -0.1640625},
        {0.640625, -0.0078125, -0.4296875}, {-0.640625, -0.0078125, -0.4296875},
        {0.3359375, 0.0546875, -0.6640625}, {-0.3359375, 0.0546875, -0.6640625},
        {0.234375, -0.3515625, 0.40625}, {-0.234375, -0.3515625, 0.40625},
        {0.1796875, -0.4140625, 0.2578125}, {-0.1796875, -0.4140625, 0.2578125},
        {0.2890625, -0.7109375, 0.3828125}, {-0.2890625, -0.7109375, 0.3828125},
        {0.25, -0.5, 0.390625}, {-0.25, -0.5, 0.390625},
        {0.328125, -0.9140625, 0.3984375}, {-0.328125, -0.9140625, 0.3984375},
        {0.140625, -0.7578125, 0.3671875}, {-0.140625, -0.7578125, 0.3671875},
        {0.125, -0.5390625, 0.359375}, {-0.125, -0.5390625, 0.359375},
        {0.1640625, -0.9453125, 0.4375}, {-0.1640625, -0.9453125, 0.4375},
        {0.21875, -0.28125, 0.4296875}, {-0.21875, -0.28125, 0.4296875},
        {0.2109375, -0.2265625, 0.46875}, {-0.2109375, -0.2265625, 0.46875},
        {0.203125, -0.171875, 0.5}, {-0.203125, -0.171875, 0.5},
        {0.2109375, -0.390625, 0.1640625}, {-0.2109375, -0.390625, 0.1640625},
        {0.296875, -0.3125, -0.265625}, {-0.296875, -0.3125, -0.265625},
        {0.34375, -0.1484375, -0.5390625}, {-0.34375, -0.1484375, -0.5390625},
        {0.453125, 0.8671875, -0.3828125}, {-0.453125, 0.8671875, -0.3828125},
        {0.453125, 0.9296875, -0.0703125}, {-0.453125, 0.9296875, -0.0703125},
        {0.453125, 0.8515625, 0.234375}, {-0.453125, 0.8515625, 0.234375},
        {0.4609375, 0.5234375, 0.4296875}, {-0.4609375, 0.5234375, 0.4296875},
        {0.7265625, 0.40625, 0.3359375}, {-0.7265625, 0.40625, 0.3359375},
        {0.6328125, 0.453125, 0.28125}, {-0.6328125, 0.453125, 0.28125},
        {0.640625, 0.703125, 0.0546875}, {-0.640625, 0.703125, 0.0546875},
        {0.796875, 0.5625, 0.125}, {-0.796875, 0.5625, 0.125},
        {0.796875, 0.6171875, -0.1171875}, {-0.796875, 0.6171875, -0.1171875},
        {0.640625, 0.75, -0.1953125}, {-0.640625, 0.75, -0.1953125},
        {0.640625, 0.6796875, -0.4453125}, {-0.640625, 0.6796875, -0.4453125},
        {0.796875, 0.5390625, -0.359375}, {-0.796875, 0.5390625, -0.359375},
        {0.6171875, 0.328125, -0.5859375}, {-0.6171875, 0.328125, -0.5859375},
        {0.484375, 0.0234375, -0.546875}, {-0.484375, 0.0234375, -0.546875},
        {0.8203125, 0.328125, -0.203125}, {-0.8203125, 0.328125, -0.203125},
        {0.40625, -0.171875, 0.1484375}, {-0.40625, -0.171875, 0.1484375},
        {0.4296875, -0.1953125, -0.2109375},
        {-0.4296875, -0.1953125, -0.2109375}, {0.890625, 0.40625, -0.234375},
        {-0.890625, 0.40625, -0.234375}, {0.7734375, -0.140625, -0.125},
        {-0.7734375, -0.140625, -0.125}, {1.0390625, -0.1015625, -0.328125},
        {-1.0390625, -0.1015625, -0.328125}, {1.28125, 0.0546875, -0.4296875},
        {-1.28125, 0.0546875, -0.4296875}, {1.3515625, 0.3203125, -0.421875},
        {-1.3515625, 0.3203125, -0.421875}, {1.234375, 0.5078125, -0.421875},
        {-1.234375, 0.5078125, -0.421875}, {1.0234375, 0.4765625, -0.3125},
        {-1.0234375, 0.4765625, -0.3125}, {1.015625, 0.4140625, -0.2890625},
        {-1.015625, 0.4140625, -0.2890625}, {1.1875, 0.4375, -0.390625},
        {-1.1875, 0.4375, -0.390625}, {1.265625, 0.2890625, -0.40625},
        {-1.265625, 0.2890625, -0.40625}, {1.2109375, 0.078125, -0.40625},
        {-1.2109375, 0.078125, -0.40625}, {1.03125, -0.0390625, -0.3046875},
        {-1.03125, -0.0390625, -0.3046875}, {0.828125, -0.0703125, -0.1328125},
        {-0.828125, -0.0703125, -0.1328125}, {0.921875, 0.359375, -0.21875},
        {-0.921875, 0.359375, -0.21875}, {0.9453125, 0.3046875, -0.2890625},
        {-0.9453125, 0.3046875, -0.2890625},
        {0.8828125, -0.0234375, -0.2109375},
        {-0.8828125, -0.0234375, -0.2109375}, {1.0390625, 0.0, -0.3671875},
        {-1.0390625, 0.0, -0.3671875}, {1.1875, 0.09375, -0.4453125},
        {-1.1875, 0.09375, -0.4453125}, {1.234375, 0.25, -0.4453125},
        {-1.234375, 0.25, -0.4453125}, {1.171875, 0.359375, -0.4375},
        {-1.171875, 0.359375, -0.4375}, {1.0234375, 0.34375, -0.359375},
        {-1.0234375, 0.34375, -0.359375}, {0.84375, 0.2890625, -0.2109375},
        {-0.84375, 0.2890625, -0.2109375}, {0.8359375, 0.171875, -0.2734375},
        {-0.8359375, 0.171875, -0.2734375}, {0.7578125, 0.09375, -0.2734375},
        {-0.7578125, 0.09375, -0.2734375}, {0.8203125, 0.0859375, -0.2734375},
        {-0.8203125, 0.0859375, -0.2734375}, {0.84375, 0.015625, -0.2734375},
        {-0.84375, 0.015625, -0.2734375}, {0.8125, -0.015625, -0.2734375},
        {-0.8125, -0.015625, -0.2734375}, {0.7265625, 0.0, -0.0703125},
        {-0.7265625, 0.0, -0.0703125}, {0.71875, -0.0234375, -0.171875},
        {-0.71875, -0.0234375, -0.171875}, {0.71875, 0.0390625, -0.1875},
        {-0.71875, 0.0390625, -0.1875}, {0.796875, 0.203125, -0.2109375},
        {-0.796875, 0.203125, -0.2109375}, {0.890625, 0.2421875, -0.265625},
        {-0.890625, 0.2421875, -0.265625}, {0.890625, 0.234375, -0.3203125},
        {-0.890625, 0.234375, -0.3203125}, {0.8125, -0.015625, -0.3203125},
        {-0.8125, -0.015625, -0.3203125}, {0.8515625, 0.015625, -0.3203125},
        {-0.8515625, 0.015625, -0.3203125}, {0.828125, 0.078125, -0.3203125},
        {-0.828125, 0.078125, -0.3203125}, {0.765625, 0.09375, -0.3203125},
        {-0.765625, 0.09375, -0.3203125}, {0.84375, 0.171875, -0.3203125},
        {-0.84375, 0.171875, -0.3203125}, {1.0390625, 0.328125, -0.4140625},
        {-1.0390625, 0.328125, -0.4140625}, {1.1875, 0.34375, -0.484375},
        {-1.1875, 0.34375, -0.484375}, {1.2578125, 0.2421875, -0.4921875},
        {-1.2578125, 0.2421875, -0.4921875}, {1.2109375, 0.0859375, -0.484375},
        {-1.2109375, 0.0859375, -0.484375}, {1.046875, 0.0, -0.421875},
        {-1.046875, 0.0, -0.421875}, {0.8828125, -0.015625, -0.265625},
        {-0.8828125, -0.015625, -0.265625}, {0.953125, 0.2890625, -0.34375},
        {-0.953125, 0.2890625, -0.34375}, {0.890625, 0.109375, -0.328125},
        {-0.890625, 0.109375, -0.328125}, {0.9375, 0.0625, -0.3359375},
        {-0.9375, 0.0625, -0.3359375}, {1.0, 0.125, -0.3671875},
        {-1.0, 0.125, -0.3671875}, {0.9609375, 0.171875, -0.3515625},
        {-0.9609375, 0.171875, -0.3515625}, {1.015625, 0.234375, -0.375},
        {-1.015625, 0.234375, -0.375}, {1.0546875, 0.1875, -0.3828125},
        {-1.0546875, 0.1875, -0.3828125}, {1.109375, 0.2109375, -0.390625},
        {-1.109375, 0.2109375, -0.390625}, {1.0859375, 0.2734375, -0.390625},
        {-1.0859375, 0.2734375, -0.390625}, {1.0234375, 0.4375, -0.484375},
        {-1.0234375, 0.4375, -0.484375}, {1.25, 0.46875, -0.546875},
        {-1.25, 0.46875, -0.546875}, {1.3671875, 0.296875, -0.5},
        {-1.3671875, 0.296875, -0.5}, {1.3125, 0.0546875, -0.53125},
        {-1.3125, 0.0546875, -0.53125}, {1.0390625, -0.0859375, -0.4921875},
        {-1.0390625, -0.0859375, -0.4921875}, {0.7890625, -0.125, -0.328125},
        {-0.7890625, -0.125, -0.328125}, {0.859375, 0.3828125, -0.3828125},
        {-0.859375, 0.3828125, -0.3828125}};
    static auto suzanne_triangles = std::vector<vec3i>{{60, 64, 48},
        {49, 65, 61}, {62, 64, 60}, {61, 65, 63}, {60, 58, 62}, {63, 59, 61},
        {60, 56, 58}, {59, 57, 61}, {60, 54, 56}, {57, 55, 61}, {60, 52, 54},
        {55, 53, 61}, {60, 50, 52}, {53, 51, 61}, {60, 48, 50}, {51, 49, 61},
        {224, 228, 226}, {227, 229, 225}, {72, 283, 73}, {73, 284, 72},
        {341, 347, 383}, {384, 348, 342}, {299, 345, 343}, {344, 346, 300},
        {323, 379, 351}, {352, 380, 324}, {441, 443, 445}, {446, 444, 442},
        {463, 491, 465}, {466, 492, 464}, {495, 497, 499}, {500, 498, 496}};
    static auto suzanne_quads = std::vector<vec4i>{{46, 0, 2, 44},
        {3, 1, 47, 45}, {44, 2, 4, 42}, {5, 3, 45, 43}, {2, 8, 6, 4},
        {7, 9, 3, 5}, {0, 10, 8, 2}, {9, 11, 1, 3}, {10, 12, 14, 8},
        {15, 13, 11, 9}, {8, 14, 16, 6}, {17, 15, 9, 7}, {14, 20, 18, 16},
        {19, 21, 15, 17}, {12, 22, 20, 14}, {21, 23, 13, 15}, {22, 24, 26, 20},
        {27, 25, 23, 21}, {20, 26, 28, 18}, {29, 27, 21, 19}, {26, 32, 30, 28},
        {31, 33, 27, 29}, {24, 34, 32, 26}, {33, 35, 25, 27}, {34, 36, 38, 32},
        {39, 37, 35, 33}, {32, 38, 40, 30}, {41, 39, 33, 31}, {38, 44, 42, 40},
        {43, 45, 39, 41}, {36, 46, 44, 38}, {45, 47, 37, 39}, {46, 36, 50, 48},
        {51, 37, 47, 49}, {36, 34, 52, 50}, {53, 35, 37, 51}, {34, 24, 54, 52},
        {55, 25, 35, 53}, {24, 22, 56, 54}, {57, 23, 25, 55}, {22, 12, 58, 56},
        {59, 13, 23, 57}, {12, 10, 62, 58}, {63, 11, 13, 59}, {10, 0, 64, 62},
        {65, 1, 11, 63}, {0, 46, 48, 64}, {49, 47, 1, 65}, {88, 173, 175, 90},
        {175, 174, 89, 90}, {86, 171, 173, 88}, {174, 172, 87, 89},
        {84, 169, 171, 86}, {172, 170, 85, 87}, {82, 167, 169, 84},
        {170, 168, 83, 85}, {80, 165, 167, 82}, {168, 166, 81, 83},
        {78, 91, 145, 163}, {146, 92, 79, 164}, {91, 93, 147, 145},
        {148, 94, 92, 146}, {93, 95, 149, 147}, {150, 96, 94, 148},
        {95, 97, 151, 149}, {152, 98, 96, 150}, {97, 99, 153, 151},
        {154, 100, 98, 152}, {99, 101, 155, 153}, {156, 102, 100, 154},
        {101, 103, 157, 155}, {158, 104, 102, 156}, {103, 105, 159, 157},
        {160, 106, 104, 158}, {105, 107, 161, 159}, {162, 108, 106, 160},
        {107, 66, 67, 161}, {67, 66, 108, 162}, {109, 127, 159, 161},
        {160, 128, 110, 162}, {127, 178, 157, 159}, {158, 179, 128, 160},
        {125, 155, 157, 178}, {158, 156, 126, 179}, {123, 153, 155, 125},
        {156, 154, 124, 126}, {121, 151, 153, 123}, {154, 152, 122, 124},
        {119, 149, 151, 121}, {152, 150, 120, 122}, {117, 147, 149, 119},
        {150, 148, 118, 120}, {115, 145, 147, 117}, {148, 146, 116, 118},
        {113, 163, 145, 115}, {146, 164, 114, 116}, {113, 180, 176, 163},
        {176, 181, 114, 164}, {109, 161, 67, 111}, {67, 162, 110, 112},
        {111, 67, 177, 182}, {177, 67, 112, 183}, {176, 180, 182, 177},
        {183, 181, 176, 177}, {134, 136, 175, 173}, {175, 136, 135, 174},
        {132, 134, 173, 171}, {174, 135, 133, 172}, {130, 132, 171, 169},
        {172, 133, 131, 170}, {165, 186, 184, 167}, {185, 187, 166, 168},
        {130, 169, 167, 184}, {168, 170, 131, 185}, {143, 189, 188, 186},
        {188, 189, 144, 187}, {184, 186, 188, 68}, {188, 187, 185, 68},
        {129, 130, 184, 68}, {185, 131, 129, 68}, {141, 192, 190, 143},
        {191, 193, 142, 144}, {139, 194, 192, 141}, {193, 195, 140, 142},
        {138, 196, 194, 139}, {195, 197, 138, 140}, {137, 70, 196, 138},
        {197, 70, 137, 138}, {189, 143, 190, 69}, {191, 144, 189, 69},
        {69, 190, 205, 207}, {206, 191, 69, 207}, {70, 198, 199, 196},
        {200, 198, 70, 197}, {196, 199, 201, 194}, {202, 200, 197, 195},
        {194, 201, 203, 192}, {204, 202, 195, 193}, {192, 203, 205, 190},
        {206, 204, 193, 191}, {198, 203, 201, 199}, {202, 204, 198, 200},
        {198, 207, 205, 203}, {206, 207, 198, 204}, {138, 139, 163, 176},
        {164, 140, 138, 176}, {139, 141, 210, 163}, {211, 142, 140, 164},
        {141, 143, 212, 210}, {213, 144, 142, 211}, {143, 186, 165, 212},
        {166, 187, 144, 213}, {80, 208, 212, 165}, {213, 209, 81, 166},
        {208, 214, 210, 212}, {211, 215, 209, 213}, {78, 163, 210, 214},
        {211, 164, 79, 215}, {130, 129, 71, 221}, {71, 129, 131, 222},
        {132, 130, 221, 219}, {222, 131, 133, 220}, {134, 132, 219, 217},
        {220, 133, 135, 218}, {136, 134, 217, 216}, {218, 135, 136, 216},
        {216, 217, 228, 230}, {229, 218, 216, 230}, {217, 219, 226, 228},
        {227, 220, 218, 229}, {219, 221, 224, 226}, {225, 222, 220, 227},
        {221, 71, 223, 224}, {223, 71, 222, 225}, {223, 230, 228, 224},
        {229, 230, 223, 225}, {182, 180, 233, 231}, {234, 181, 183, 232},
        {111, 182, 231, 253}, {232, 183, 112, 254}, {109, 111, 253, 255},
        {254, 112, 110, 256}, {180, 113, 251, 233}, {252, 114, 181, 234},
        {113, 115, 249, 251}, {250, 116, 114, 252}, {115, 117, 247, 249},
        {248, 118, 116, 250}, {117, 119, 245, 247}, {246, 120, 118, 248},
        {119, 121, 243, 245}, {244, 122, 120, 246}, {121, 123, 241, 243},
        {242, 124, 122, 244}, {123, 125, 239, 241}, {240, 126, 124, 242},
        {125, 178, 235, 239}, {236, 179, 126, 240}, {178, 127, 237, 235},
        {238, 128, 179, 236}, {127, 109, 255, 237}, {256, 110, 128, 238},
        {237, 255, 257, 275}, {258, 256, 238, 276}, {235, 237, 275, 277},
        {276, 238, 236, 278}, {239, 235, 277, 273}, {278, 236, 240, 274},
        {241, 239, 273, 271}, {274, 240, 242, 272}, {243, 241, 271, 269},
        {272, 242, 244, 270}, {245, 243, 269, 267}, {270, 244, 246, 268},
        {247, 245, 267, 265}, {268, 246, 248, 266}, {249, 247, 265, 263},
        {266, 248, 250, 264}, {251, 249, 263, 261}, {264, 250, 252, 262},
        {233, 251, 261, 279}, {262, 252, 234, 280}, {255, 253, 259, 257},
        {260, 254, 256, 258}, {253, 231, 281, 259}, {282, 232, 254, 260},
        {231, 233, 279, 281}, {280, 234, 232, 282}, {66, 107, 283, 72},
        {284, 108, 66, 72}, {107, 105, 285, 283}, {286, 106, 108, 284},
        {105, 103, 287, 285}, {288, 104, 106, 286}, {103, 101, 289, 287},
        {290, 102, 104, 288}, {101, 99, 291, 289}, {292, 100, 102, 290},
        {99, 97, 293, 291}, {294, 98, 100, 292}, {97, 95, 295, 293},
        {296, 96, 98, 294}, {95, 93, 297, 295}, {298, 94, 96, 296},
        {93, 91, 299, 297}, {300, 92, 94, 298}, {307, 308, 327, 337},
        {328, 308, 307, 338}, {306, 307, 337, 335}, {338, 307, 306, 336},
        {305, 306, 335, 339}, {336, 306, 305, 340}, {88, 90, 305, 339},
        {305, 90, 89, 340}, {86, 88, 339, 333}, {340, 89, 87, 334},
        {84, 86, 333, 329}, {334, 87, 85, 330}, {82, 84, 329, 331},
        {330, 85, 83, 332}, {329, 335, 337, 331}, {338, 336, 330, 332},
        {329, 333, 339, 335}, {340, 334, 330, 336}, {325, 331, 337, 327},
        {338, 332, 326, 328}, {80, 82, 331, 325}, {332, 83, 81, 326},
        {208, 341, 343, 214}, {344, 342, 209, 215}, {80, 325, 341, 208},
        {342, 326, 81, 209}, {78, 214, 343, 345}, {344, 215, 79, 346},
        {78, 345, 299, 91}, {300, 346, 79, 92}, {76, 323, 351, 303},
        {352, 324, 76, 303}, {303, 351, 349, 77}, {350, 352, 303, 77},
        {77, 349, 347, 304}, {348, 350, 77, 304}, {304, 347, 327, 308},
        {328, 348, 304, 308}, {325, 327, 347, 341}, {348, 328, 326, 342},
        {295, 297, 317, 309}, {318, 298, 296, 310}, {75, 315, 323, 76},
        {324, 316, 75, 76}, {301, 357, 355, 302}, {356, 358, 301, 302},
        {302, 355, 353, 74}, {354, 356, 302, 74}, {74, 353, 315, 75},
        {316, 354, 74, 75}, {291, 293, 361, 363}, {362, 294, 292, 364},
        {363, 361, 367, 365}, {368, 362, 364, 366}, {365, 367, 369, 371},
        {370, 368, 366, 372}, {371, 369, 375, 373}, {376, 370, 372, 374},
        {313, 377, 373, 375}, {374, 378, 314, 376}, {315, 353, 373, 377},
        {374, 354, 316, 378}, {353, 355, 371, 373}, {372, 356, 354, 374},
        {355, 357, 365, 371}, {366, 358, 356, 372}, {357, 359, 363, 365},
        {364, 360, 358, 366}, {289, 291, 363, 359}, {364, 292, 290, 360},
        {73, 359, 357, 301}, {358, 360, 73, 301}, {283, 285, 287, 289},
        {288, 286, 284, 290}, {283, 289, 359, 73}, {360, 290, 284, 73},
        {293, 295, 309, 361}, {310, 296, 294, 362}, {309, 311, 367, 361},
        {368, 312, 310, 362}, {311, 381, 369, 367}, {370, 382, 312, 368},
        {313, 375, 369, 381}, {370, 376, 314, 382}, {347, 349, 385, 383},
        {386, 350, 348, 384}, {317, 383, 385, 319}, {386, 384, 318, 320},
        {297, 299, 383, 317}, {384, 300, 298, 318}, {299, 343, 341, 383},
        {342, 344, 300, 384}, {313, 321, 379, 377}, {380, 322, 314, 378},
        {315, 377, 379, 323}, {380, 378, 316, 324}, {319, 385, 379, 321},
        {380, 386, 320, 322}, {349, 351, 379, 385}, {380, 352, 350, 386},
        {399, 387, 413, 401}, {414, 388, 400, 402}, {399, 401, 403, 397},
        {404, 402, 400, 398}, {397, 403, 405, 395}, {406, 404, 398, 396},
        {395, 405, 407, 393}, {408, 406, 396, 394}, {393, 407, 409, 391},
        {410, 408, 394, 392}, {391, 409, 411, 389}, {412, 410, 392, 390},
        {409, 419, 417, 411}, {418, 420, 410, 412}, {407, 421, 419, 409},
        {420, 422, 408, 410}, {405, 423, 421, 407}, {422, 424, 406, 408},
        {403, 425, 423, 405}, {424, 426, 404, 406}, {401, 427, 425, 403},
        {426, 428, 402, 404}, {401, 413, 415, 427}, {416, 414, 402, 428},
        {317, 319, 443, 441}, {444, 320, 318, 442}, {319, 389, 411, 443},
        {412, 390, 320, 444}, {309, 317, 441, 311}, {442, 318, 310, 312},
        {381, 429, 413, 387}, {414, 430, 382, 388}, {411, 417, 439, 443},
        {440, 418, 412, 444}, {437, 445, 443, 439}, {444, 446, 438, 440},
        {433, 445, 437, 435}, {438, 446, 434, 436}, {431, 447, 445, 433},
        {446, 448, 432, 434}, {429, 447, 431, 449}, {432, 448, 430, 450},
        {413, 429, 449, 415}, {450, 430, 414, 416}, {311, 447, 429, 381},
        {430, 448, 312, 382}, {311, 441, 445, 447}, {446, 442, 312, 448},
        {415, 449, 451, 475}, {452, 450, 416, 476}, {449, 431, 461, 451},
        {462, 432, 450, 452}, {431, 433, 459, 461}, {460, 434, 432, 462},
        {433, 435, 457, 459}, {458, 436, 434, 460}, {435, 437, 455, 457},
        {456, 438, 436, 458}, {437, 439, 453, 455}, {454, 440, 438, 456},
        {439, 417, 473, 453}, {474, 418, 440, 454}, {427, 415, 475, 463},
        {476, 416, 428, 464}, {425, 427, 463, 465}, {464, 428, 426, 466},
        {423, 425, 465, 467}, {466, 426, 424, 468}, {421, 423, 467, 469},
        {468, 424, 422, 470}, {419, 421, 469, 471}, {470, 422, 420, 472},
        {417, 419, 471, 473}, {472, 420, 418, 474}, {457, 455, 479, 477},
        {480, 456, 458, 478}, {477, 479, 481, 483}, {482, 480, 478, 484},
        {483, 481, 487, 485}, {488, 482, 484, 486}, {485, 487, 489, 491},
        {490, 488, 486, 492}, {463, 475, 485, 491}, {486, 476, 464, 492},
        {451, 483, 485, 475}, {486, 484, 452, 476}, {451, 461, 477, 483},
        {478, 462, 452, 484}, {457, 477, 461, 459}, {462, 478, 458, 460},
        {453, 473, 479, 455}, {480, 474, 454, 456}, {471, 481, 479, 473},
        {480, 482, 472, 474}, {469, 487, 481, 471}, {482, 488, 470, 472},
        {467, 489, 487, 469}, {488, 490, 468, 470}, {465, 491, 489, 467},
        {490, 492, 466, 468}, {391, 389, 503, 501}, {504, 390, 392, 502},
        {393, 391, 501, 499}, {502, 392, 394, 500}, {395, 393, 499, 497},
        {500, 394, 396, 498}, {397, 395, 497, 495}, {498, 396, 398, 496},
        {399, 397, 495, 493}, {496, 398, 400, 494}, {387, 399, 493, 505},
        {494, 400, 388, 506}, {493, 501, 503, 505}, {504, 502, 494, 506},
        {493, 495, 499, 501}, {500, 496, 494, 502}, {313, 381, 387, 505},
        {388, 382, 314, 506}, {313, 505, 503, 321}, {504, 506, 314, 322},
        {319, 321, 503, 389}, {504, 322, 320, 390}};

    pos = suzanne_pos;
    quads = suzanne_quads;
    for (auto t : suzanne_triangles) { quads.push_back({t.x, t.y, t.z, t.z}); }
}

// Watertight cube
void make_cube(std::vector<vec4i>& quads, std::vector<vec3f>& pos) {
    static auto cube_pos = std::vector<vec3f>{{-1, -1, -1}, {-1, +1, -1},
        {+1, +1, -1}, {+1, -1, -1}, {-1, -1, +1}, {-1, +1, +1}, {+1, +1, +1},
        {+1, -1, +1}};
    static auto cube_quads = std::vector<vec4i>{{0, 1, 2, 3}, {7, 6, 5, 4},
        {4, 5, 1, 0}, {6, 7, 3, 2}, {2, 1, 5, 6}, {0, 3, 7, 4}};
    static auto cube_quad_uv =
        std::vector<vec2f>{{0, 0}, {1, 0}, {1, 1}, {0, 1}};
    quads = cube_quads;
    pos = cube_pos;
}

// Generate lines set along a quad.
void make_lines(std::vector<vec2i>& lines, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<float>& radius, const vec2i& steps, const vec2f& size,
    const vec2f& uvsize, vec2f line_radius) {
    auto nverts = (steps.x + 1) * steps.y;
    auto nlines = steps.x * steps.y;
    auto vid = [steps](int i, int j) { return j * (steps.x + 1) + i; };
    auto fid = [steps](int i, int j) { return j * steps.x + i; };

    pos.resize(nverts);
    norm.resize(nverts);
    texcoord.resize(nverts);
    radius.resize(nverts);
    if (steps.y > 1) {
        for (auto j = 0; j < steps.y; j++) {
            for (auto i = 0; i <= steps.x; i++) {
                auto uv = vec2f{i / (float)steps.x,
                    j / (float)(steps.y > 1 ? steps.y - 1 : 1)};
                pos[vid(i, j)] = {
                    (uv.x - 0.5f) * size.x, (uv.y - 0.5f) * size.y, 0};
                norm[vid(i, j)] = {1, 0, 0};
                texcoord[vid(i, j)] = uv * uvsize;
            }
        }
    } else {
        for (auto i = 0; i <= steps.x; i++) {
            auto uv = vec2f{i / (float)steps.x, 0};
            pos[vid(i, 0)] = {(uv.x - 0.5f) * size.x, 0, 0};
            norm[vid(i, 0)] = {1, 0, 0};
            texcoord[vid(i, 0)] = uv * uvsize;
        }
    }

    lines.resize(nlines);
    for (int j = 0; j < steps.y; j++) {
        for (int i = 0; i < steps.x; i++) {
            lines[fid(i, j)] = {vid(i, j), vid(i + 1, j)};
        }
    }
}

// Generate a point set with points placed at the origin with texcoords varying
// along u.
void make_points(std::vector<int>& points, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<float>& radius, int num, float uvsize, float point_radius) {
    points.resize(num);
    for (auto i = 0; i < num; i++) points[i] = i;
    pos.assign(num, {0, 0, 0});
    norm.assign(num, {0, 0, 1});
    texcoord.assign(num, {0, 0});
    radius.assign(num, point_radius);
    for (auto i = 0; i < texcoord.size(); i++)
        texcoord[i] = {(float)i / (float)num, 0};
}

// Generate a point set.
void make_random_points(std::vector<int>& points, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<float>& radius, int num, const vec3f& size, float uvsize,
    float point_radius, uint64_t seed) {
    make_points(points, pos, norm, texcoord, radius, num, uvsize, point_radius);
    auto rng = make_rng(seed);
    for (auto i = 0; i < pos.size(); i++) {
        pos[i] = (rand3f(rng) - vec3f{0.5f, 0.5f, 0.5f}) * size;
    }
}

// Make a point.
void make_point(std::vector<int>& points, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<float>& radius, float point_radius) {
    points = {0};
    pos = {{0, 0, 0}};
    norm = {{0, 0, 1}};
    texcoord = {{0, 0}};
    radius = {point_radius};
}

// Make a bezier circle. Returns bezier, pos.
void make_bezier_circle(std::vector<vec4i>& beziers, std::vector<vec3f>& pos) {
    // constant from http://spencermortensen.com/articles/bezier-circle/
    auto c = 0.551915024494f;
    pos = std::vector<vec3f>{{1, 0, 0}, {1, c, 0}, {c, 1, 0}, {0, 1, 0},
        {-c, 1, 0}, {-1, c, 0}, {-1, 0, 0}, {-1, -c, 0}, {-c, -1, 0},
        {0, -1, 0}, {c, -1, 0}, {1, -c, 0}};
    beziers = std::vector<vec4i>{
        {0, 1, 2, 3}, {3, 4, 5, 6}, {6, 7, 8, 9}, {9, 10, 11, 0}};
}

// Make a hair ball around a shape
void make_hair(std::vector<vec2i>& lines, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<float>& radius, const vec2i& steps,
    const std::vector<vec3i>& striangles, const std::vector<vec3f>& spos,
    const std::vector<vec3f>& snorm, const std::vector<vec2f>& stexcoord,
    const vec2f& len, const vec2f& rad, const vec2f& noise, const vec2f& clump,
    const vec2f& rotation, int seed) {
    std::vector<vec3f> bpos;
    std::vector<vec3f> bnorm;
    std::vector<vec2f> btexcoord;
    std::tie(bpos, bnorm, btexcoord) = sample_triangles_points(
        striangles, spos, snorm, stexcoord, steps.y, seed);

    auto rng = make_rng(seed, 3);
    auto blen = std::vector<float>(bpos.size());
    for (auto& l : blen) l = lerp(len.x, len.y, rand1f(rng));

    auto cidx = std::vector<int>();
    if (clump.x > 0) {
        for (auto bidx = 0; bidx < bpos.size(); bidx++) {
            cidx.push_back(0);
            auto cdist = flt_max;
            for (auto c = 0; c < clump.y; c++) {
                auto d = length(bpos[bidx] - bpos[c]);
                if (d < cdist) {
                    cdist = d;
                    cidx.back() = c;
                }
            }
        }
    }

    make_lines(lines, pos, norm, texcoord, radius, steps, {1, 1}, {1, 1});
    for (auto i = 0; i < pos.size(); i++) {
        auto u = texcoord[i].x;
        auto bidx = i / (steps.x + 1);
        pos[i] = bpos[bidx] + bnorm[bidx] * u * blen[bidx];
        norm[i] = bnorm[bidx];
        radius[i] = lerp(rad.x, rad.y, u);
        if (clump.x > 0) {
            pos[i] = pos[i] +
                     (pos[i + (cidx[bidx] - bidx) * (steps.x + 1)] - pos[i]) *
                         u * clump.x;
        }
        if (noise.x > 0) {
            auto nx = perlin_noise(pos[i] * noise.y + vec3f{0, 0, 0}) * noise.x;
            auto ny =
                perlin_noise(pos[i] * noise.y + vec3f{3, 7, 11}) * noise.x;
            auto nz =
                perlin_noise(pos[i] * noise.y + vec3f{13, 17, 19}) * noise.x;
            pos[i] += {nx, ny, nz};
        }
    }

    if (clump.x > 0 || noise.x > 0 || rotation.x > 0)
        compute_tangents(lines, pos, norm);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR COLOR UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Convert between CIE XYZ and xyY
vec3f xyz_to_xyY(const vec3f& xyz) {
    if (xyz == zero3f) return zero3f;
    return {xyz.x / (xyz.x + xyz.y + xyz.z), xyz.y / (xyz.x + xyz.y + xyz.z),
        xyz.y};
}
// Convert between CIE XYZ and xyY
vec3f xyY_to_xyz(const vec3f& xyY) {
    if (xyY.y == 0) return zero3f;
    return {xyY.x * xyY.z / xyY.y, xyY.z, (1 - xyY.x - xyY.y) * xyY.z / xyY.y};
}
// Convert between CIE XYZ and RGB
vec3f xyz_to_rgb(const vec3f& xyz) {
    // from http://www.brucelindbloom.com/index.html?Eqn_RGB_to_XYZ.html
    if (xyz == zero3f) return zero3f;
    return {+3.2404542f * xyz.x - 1.5371385f * xyz.y - 0.4985314f * xyz.z,
        -0.9692660f * xyz.x + 1.8760108f * xyz.y + 0.0415560f * xyz.z,
        +0.0556434f * xyz.x - 0.2040259f * xyz.y + 1.0572252f * xyz.z};
}
// Convert between CIE XYZ and RGB
vec3f rgb_to_xyz(const vec3f& rgb) {
    // from http://www.brucelindbloom.com/index.html?Eqn_RGB_to_XYZ.html
    if (rgb == zero3f) return zero3f;
    return {0.4124564f * rgb.x + 0.3575761f * rgb.y + 0.1804375f * rgb.z,
        0.2126729f * rgb.x + 0.7151522f * rgb.y + 0.0721750f * rgb.z,
        0.0193339f * rgb.x + 0.1191920f * rgb.y + 0.9503041f * rgb.z};
}

// Convert HSV to RGB
vec3f hsv_to_rgb(const vec3f& hsv) {
    // from Imgui.cpp
    auto h = hsv.x, s = hsv.y, v = hsv.z;
    if (hsv.y == 0.0f) return {v, v, v};

    h = fmodf(h, 1.0f) / (60.0f / 360.0f);
    int i = (int)h;
    float f = h - (float)i;
    float p = v * (1.0f - s);
    float q = v * (1.0f - s * f);
    float t = v * (1.0f - s * (1.0f - f));

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
vec3f rgb_to_hsv(vec3f rgb) {
    // from Imgui.cpp
    auto r = rgb.x, g = rgb.y, b = rgb.z;
    float K = 0.f;
    if (g < b) {
        std::swap(g, b);
        K = -1.f;
    }
    if (r < g) {
        std::swap(r, g);
        K = -2.f / 6.f - K;
    }

    float chroma = r - (g < b ? g : b);
    return {
        fabsf(K + (g - b) / (6.f * chroma + 1e-20f)), chroma / (r + 1e-20f), r};
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR IMAGE UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Conversion between linear and gamma-encoded images.
image4f gamma_to_linear(const image4f& srgb, float gamma) {
    auto lin = image4f{srgb.width(), srgb.height()};
    for (auto i = 0; i < srgb.size(); i++) {
        xyz(lin[i]) = gamma_to_linear(xyz(srgb[i]), gamma);
        lin[i].w = srgb[i].w;
    }
    return lin;
}
image4f linear_to_gamma(const image4f& lin, float gamma) {
    auto srgb = image4f{lin.width(), lin.height()};
    for (auto i = 0; i < lin.size(); i++) {
        xyz(srgb[i]) = linear_to_gamma(xyz(lin[i]), gamma);
        srgb[i].w = lin[i].w;
    }
    return srgb;
}

// Conversion from/to floats.
image4f byte_to_float(const image4b& bt) {
    auto fl = image4f{bt.width(), bt.height()};
    for (auto i = 0; i < bt.size(); i++) fl[i] = byte_to_float(bt[i]);
    return fl;
}
image4b float_to_byte(const image4f& fl) {
    auto bt = image4b{fl.width(), fl.height()};
    for (auto i = 0; i < fl.size(); i++) bt[i] = float_to_byte(fl[i]);
    return bt;
}

vec3f filmic_tonemap(const vec3f& hdr) {
    // https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
    auto x = hdr;
    // x *= 0.6; // brings it back to ACES range
    x = (x * (2.51f * x + vec3f{0.03f, 0.03f, 0.03f})) /
        (x * (2.43f * x + vec3f{0.59f, 0.59f, 0.59f}) +
            vec3f{0.14f, 0.14f, 0.14f});
    return x;
}

// Tone mapping HDR to LDR images.
image4f expose_image(const image4f& hdr, float exposure) {
    if (!exposure) return hdr;
    auto ldr = image4f{hdr.width(), hdr.height()};
    auto scale = pow(2.0f, exposure);
    for (auto i = 0; i < hdr.size(); i++) {
        xyz(ldr[i]) = xyz(hdr[i]) * scale;
        ldr[i].w = hdr[i].w;
    }
    return ldr;
}

// Tone mapping HDR to LDR images.
image4f filmic_tonemap_image(const image4f& hdr) {
    auto ldr = image4f{hdr.width(), hdr.height()};
    for (auto i = 0; i < hdr.size(); i++) {
        xyz(ldr[i]) = filmic_tonemap(xyz(hdr[i]));
        ldr[i].w = hdr[i].w;
    }
    return ldr;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR IMAGE EXAMPLES
// -----------------------------------------------------------------------------
namespace ygl {

// Make a grid image
image4f make_grid_image(
    int width, int height, int tiles, const vec4f& c0, const vec4f& c1) {
    auto img = image4f{width, height};
    auto tile = width / tiles;
    for (int j = 0; j < width; j++) {
        for (int i = 0; i < height; i++) {
            auto c = i % tile == 0 || i % tile == tile - 1 || j % tile == 0 ||
                     j % tile == tile - 1;
            img[{i, j}] = (c) ? c0 : c1;
        }
    }
    return img;
}

// Make a checkerboard image
image4f make_checker_image(
    int width, int height, int tiles, const vec4f& c0, const vec4f& c1) {
    auto img = image4f{width, height};
    auto tile = width / tiles;
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto c = (i / tile + j / tile) % 2 == 0;
            img[{i, j}] = (c) ? c0 : c1;
        }
    }
    return img;
}

// Make an image with bumps and dimples.
image4f make_bumpdimple_image(int width, int height, int tiles) {
    auto img = image4f{width, height};
    auto tile = width / tiles;
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto c = (i / tile + j / tile) % 2 == 0;
            auto ii = i % tile - tile / 2, jj = j % tile - tile / 2;
            auto r =
                sqrt(float(ii * ii + jj * jj)) / sqrt(float(tile * tile) / 4);
            auto h = 0.5f;
            if (r < 0.5f) { h += (c) ? (0.5f - r) : -(0.5f - r); }
            img[{i, j}] = {h, h, h, 1};
        }
    }
    return img;
}

// Make a uv colored grid
image4f make_ramp_image(
    int width, int height, const vec4f& c0, const vec4f& c1) {
    auto img = image4f{width, height};
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto u = (float)i / (float)width;
            img[{i, j}] = c0 * (1 - u) + c1 * u;
        }
    }
    return img;
}

// Make a gamma ramp image
image4f make_gammaramp_imagef(int width, int height) {
    auto img = image4f{width, height};
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto u = j / float(height - 1);
            if (i < width / 3) u = pow(u, 2.2f);
            if (i > (width * 2) / 3) u = pow(u, 1 / 2.2f);
            img[{i, j}] = {u, u, u, 1};
        }
    }
    return img;
}

// Make an image color with red/green in the [0,1] range. Helpful to visualize
// uv texture coordinate application.
image4f make_uv_image(int width, int height) {
    auto img = image4f{width, height};
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            img[{i, j}] = {
                i / (float)(width - 1), j / (float)(height - 1), 0, 1};
        }
    }
    return img;
}

// Make a uv colored grid
image4f make_uvgrid_image(int width, int height, int tiles, bool colored) {
    auto img = image4f{width, height};
    auto tile = width / tiles;
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto ii = i / tile, jj = j / tile;
            auto ww = width / tile, hh = height / tile;
            auto ph =
                (((256 / (ww * hh)) * (ii + jj * ww) - 64 + 256) % 256) / 360.f;
            auto pv = 0.5f;
            auto ps = 0.8f;
            if (i % (tile / 2) && j % (tile / 2)) {
                if ((i / tile + j / tile) % 2)
                    pv += 0.05f;
                else
                    pv -= 0.05f;
            } else {
                pv = 0.8f;
                ps = 0.2f;
            }
            auto rgb = (colored) ? hsv_to_rgb({ph, ps, pv}) : vec3f{pv, pv, pv};
            img[{i, height - j - 1}] = {rgb.x, rgb.y, rgb.z, 1};
        }
    }
    return img;
}

// Comvert a bump map to a normal map.
image4f bump_to_normal_map(const image4f& img, float scale) {
    auto norm = image4f{img.width(), img.height()};
    for (int j = 0; j < img.height(); j++) {
        for (int i = 0; i < img.width(); i++) {
            auto i1 = (i + 1) % img.width(), j1 = (j + 1) % img.height();
            auto p00 = img[{i, j}], p10 = img[{i1, j}], p01 = img[{i, j1}];
            auto g00 = (float(p00.x) + float(p00.y) + float(p00.z)) / (3 * 255);
            auto g01 = (float(p01.x) + float(p01.y) + float(p01.z)) / (3 * 255);
            auto g10 = (float(p10.x) + float(p10.y) + float(p10.z)) / (3 * 255);
            auto n = vec3f{scale * (g00 - g10), scale * (g00 - g01), 1.0f};
            n.y = -n.y;  // make green pointing up, even if y axis points down
            n = normalize(n) * 0.5f + vec3f{0.5f, 0.5f, 0.5f};
            norm[{i, j}] = {n.x, n.y, n.z, 1};
        }
    }
    return norm;
}

// Implementation of sunsky modified heavily from pbrt
image4f make_sunsky_image(int width, int height, float thetaSun,
    float turbidity, bool has_sun, vec3f ground_albedo) {
    auto wSun = vec3f{0, cos(thetaSun), sin(thetaSun)};

    // sunSpectralRad =  ComputeAttenuatedSunlight(thetaS, turbidity);
    auto sunAngularRadius = 9.35e-03f / 2;  // Wikipedia
    auto thetaS = thetaSun;

    auto t1 = thetaSun, t2 = thetaSun * thetaSun,
         t3 = thetaSun * thetaSun * thetaSun;
    auto T = turbidity;
    auto T2 = turbidity * turbidity;

    auto zenith_xyY = vec3f{
        (+0.00165f * t3 - 0.00374f * t2 + 0.00208f * t1 + 0) * T2 +
            (-0.02902f * t3 + 0.06377f * t2 - 0.03202f * t1 + 0.00394f) * T +
            (+0.11693f * t3 - 0.21196f * t2 + 0.06052f * t1 + 0.25885f),
        (+0.00275f * t3 - 0.00610f * t2 + 0.00316f * t1 + 0) * T2 +
            (-0.04214f * t3 + 0.08970f * t2 - 0.04153f * t1 + 0.00515f) * T +
            (+0.15346f * t3 - 0.26756f * t2 + 0.06669f * t1 + 0.26688f),
        1000 * (4.0453f * T - 4.9710f) *
                tan((4.0f / 9.0f - T / 120.0f) * (pi - 2 * t1)) -
            .2155f * T + 2.4192f};

    auto perez_A_xyY = vec3f{-0.01925f * T - 0.25922f, -0.01669f * T - 0.26078f,
        +0.17872f * T - 1.46303f};
    auto perez_B_xyY = vec3f{-0.06651f * T + 0.00081f, -0.09495f * T + 0.00921f,
        -0.35540f * T + 0.42749f};
    auto perez_C_xyY = vec3f{-0.00041f * T + 0.21247f, -0.00792f * T + 0.21023f,
        -0.02266f * T + 5.32505f};
    auto perez_D_xyY = vec3f{-0.06409f * T - 0.89887f, -0.04405f * T - 1.65369f,
        +0.12064f * T - 2.57705f};
    auto perez_E_xyY = vec3f{-0.00325f * T + 0.04517f, -0.01092f * T + 0.05291f,
        -0.06696f * T + 0.37027f};

    auto perez_f = [thetaS](float A, float B, float C, float D, float E,
                       float theta, float gamma, float zenith) -> float {
        auto den = ((1 + A * exp(B)) *
                    (1 + C * exp(D * thetaS) + E * cos(thetaS) * cos(thetaS)));
        auto num = ((1 + A * exp(B / cos(theta))) *
                    (1 + C * exp(D * gamma) + E * cos(gamma) * cos(gamma)));
        return zenith * num / den;
    };

    auto sky = [&perez_f, perez_A_xyY, perez_B_xyY, perez_C_xyY, perez_D_xyY,
                   perez_E_xyY, zenith_xyY](auto theta, auto gamma) -> vec3f {
        auto x = perez_f(perez_A_xyY.x, perez_B_xyY.x, perez_C_xyY.x,
            perez_D_xyY.x, perez_E_xyY.x, theta, gamma, zenith_xyY.x);
        auto y = perez_f(perez_A_xyY.y, perez_B_xyY.y, perez_C_xyY.y,
            perez_D_xyY.y, perez_E_xyY.y, theta, gamma, zenith_xyY.y);
        auto Y = perez_f(perez_A_xyY.z, perez_B_xyY.z, perez_C_xyY.z,
            perez_D_xyY.z, perez_E_xyY.z, theta, gamma, zenith_xyY.z);
        return xyz_to_rgb(xyY_to_xyz({x, y, Y})) / 10000.0f;
    };

    // compute sun luminance
    // TODO: how this relates to zenith intensity?
    auto sun_ko = vec3f{0.48f, 0.75f, 0.14f};
    auto sun_kg = vec3f{0.1f, 0.0f, 0.0f};
    auto sun_kwa = vec3f{0.02f, 0.0f, 0.0f};
    auto sun_sol = vec3f{20000.0f, 27000.0f, 30000.0f};
    auto sun_lambda = vec3f{680, 530, 480};
    auto sun_beta = 0.04608365822050f * turbidity - 0.04586025928522f;
    auto sun_m =
        1.0f / (cos(thetaSun) + 0.000940f * pow(1.6386f - thetaSun, -1.253f));

    auto sun_le = zero3f;
    for (auto i = 0; i < 3; i++) {
        auto tauR =
            exp(-sun_m * 0.008735f * pow((&sun_lambda.x)[i] / 1000, -4.08f));
        auto tauA =
            exp(-sun_m * sun_beta * pow((&sun_lambda.x)[i] / 1000, -1.3f));
        auto tauO = exp(-sun_m * (&sun_ko.x)[i] * .35f);
        auto tauG = exp(-1.41f * (&sun_kg.x)[i] * sun_m /
                        pow(1 + 118.93f * (&sun_kg.x)[i] * sun_m, 0.45f));
        auto tauWA =
            exp(-0.2385f * (&sun_kwa.x)[i] * 2.0f * sun_m /
                pow(1 + 20.07f * (&sun_kwa.x)[i] * 2.0f * sun_m, 0.45f));
        (&sun_le.x)[i] = (&sun_sol.x)[i] * tauR * tauA * tauO * tauG * tauWA;
    }

    auto sun = [has_sun, sunAngularRadius, sun_le](auto theta, auto gamma) {
        return (has_sun && gamma < sunAngularRadius) ? sun_le / 10000.0f :
                                                       zero3f;
    };

    auto img = image4f{width, height, {0, 0, 0, 1}};
    for (auto j = 0; j < height / 2; j++) {
        auto theta = pi * ((j + 0.5f) / height);
        theta = clamp(theta, 0.0f, pi / 2 - flt_eps);
        for (int i = 0; i < width; i++) {
            auto phi = 2 * pi * (float(i + 0.5f) / width);
            auto w =
                vec3f{cos(phi) * sin(theta), cos(theta), sin(phi) * sin(theta)};
            auto gamma = acos(clamp(dot(w, wSun), -1.0f, 1.0f));
            auto col = sky(theta, gamma) + sun(theta, gamma);
            img[{i, j}] = {col.x, col.y, col.z, 1};
        }
    }

    if (ground_albedo != zero3f) {
        auto ground = zero3f;
        for (auto j = 0; j < height / 2; j++) {
            auto theta = pi * ((j + 0.5f) / height);
            for (int i = 0; i < width; i++) {
                auto pxl = img[{i, j}];
                auto le = vec3f{pxl.x, pxl.y, pxl.z};
                auto angle = sin(theta) * 4 * pi / (width * height);
                ground += le * (ground_albedo / pi) * cos(theta) * angle;
            }
        }
        for (auto j = height / 2; j < height; j++) {
            for (int i = 0; i < width; i++) {
                img[{i, j}] = {ground.x, ground.y, ground.z, 1};
            }
        }
    }

    return img;
}

// Make an image of multiple lights.
image4f make_lights_image(int width, int height, vec3f le, int nlights,
    float langle, float lwidth, float lheight) {
    auto img = image4f{width, height, {0, 0, 0, 1}};
    for (auto j = 0; j < height / 2; j++) {
        auto theta = pi * ((j + 0.5f) / height);
        theta = clamp(theta, 0.0f, pi / 2 - flt_eps);
        if (fabs(theta - langle) > lheight / 2) continue;
        for (int i = 0; i < width; i++) {
            auto phi = 2 * pi * (float(i + 0.5f) / width);
            auto inlight = false;
            for (auto l = 0; l < nlights; l++) {
                auto lphi = 2 * pi * (l + 0.5f) / nlights;
                inlight = inlight || fabs(phi - lphi) < lwidth / 2;
            }
            img[{i, j}] = {le.x, le.y, le.z, 1};
        }
    }
    return img;
}

// Make a noise image. Wrap works only if both width and height are powers of
// two.
image4f make_noise_image(int width, int height, float scale, bool wrap) {
    auto wrap3i = (wrap) ? vec3i{width, height, 2} : zero3i;
    auto img = image4f{width, height};
    for (auto j = 0; j < height; j++) {
        for (auto i = 0; i < width; i++) {
            auto p = vec3f{i / (float)width, j / (float)height, 0.5f} * scale;
            auto g = perlin_noise(p, wrap3i);
            g = clamp(0.5f + 0.5f * g, 0.0f, 1.0f);
            img[{i, j}] = {g, g, g, 1};
        }
    }
    return img;
}

// Make a noise image. Wrap works only if both width and height are powers of
// two.
image4f make_fbm_image(int width, int height, float scale, float lacunarity,
    float gain, int octaves, bool wrap) {
    auto wrap3i = (wrap) ? vec3i{width, height, 2} : zero3i;
    auto img = image4f{width, height};
    for (auto j = 0; j < height; j++) {
        for (auto i = 0; i < width; i++) {
            auto p = vec3f{i / (float)width, j / (float)height, 0.5f} * scale;
            auto g = perlin_fbm_noise(p, lacunarity, gain, octaves, wrap3i);
            g = clamp(0.5f + 0.5f * g, 0.0f, 1.0f);
            img[{i, j}] = {g, g, g, 1};
        }
    }
    return img;
}

// Make a noise image. Wrap works only if both width and height are powers of
// two.
image4f make_ridge_image(int width, int height, float scale, float lacunarity,
    float gain, float offset, int octaves, bool wrap) {
    auto wrap3i = (wrap) ? vec3i{width, height, 2} : zero3i;
    auto img = image4f{width, height};
    for (auto j = 0; j < height; j++) {
        for (auto i = 0; i < width; i++) {
            auto p = vec3f{i / (float)width, j / (float)height, 0.5f} * scale;
            auto g = perlin_ridge_noise(
                p, lacunarity, gain, offset, octaves, wrap3i);
            g = clamp(g, 0.0f, 1.0f);
            img[{i, j}] = {g, g, g, 1};
        }
    }
    return img;
}

// Make a noise image. Wrap works only if both width and height are powers of
// two.
image4f make_turbulence_image(int width, int height, float scale,
    float lacunarity, float gain, int octaves, bool wrap) {
    auto wrap3i = (wrap) ? vec3i{width, height, 2} : zero3i;
    auto img = image4f{width, height};
    for (auto j = 0; j < height; j++) {
        for (auto i = 0; i < width; i++) {
            auto p = vec3f{i / (float)width, j / (float)height, 0.5f} * scale;
            auto g =
                perlin_turbulence_noise(p, lacunarity, gain, octaves, wrap3i);
            g = clamp(g, 0.0f, 1.0f);
            img[{i, j}] = {g, g, g, 1};
        }
    }
    return img;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SCENE UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Compute shape normals.
void update_normals(const std::shared_ptr<shape>& shp) {
    if (!shp->triangles.empty()) {
        compute_normals(shp->triangles, shp->pos, shp->norm);
    } else if (!shp->lines.empty()) {
        compute_tangents(shp->lines, shp->pos, shp->norm);
    } else {
        shp->norm.assign(shp->pos.size(), {0, 0, 1});
    }
}

// Computes a shape bounding box.
void update_bbox(const std::shared_ptr<shape>& shp) {
    shp->bbox = invalid_bbox3f;
    for (auto p : shp->pos) shp->bbox += p;
}

// Updates the scene and scene's instances bounding boxes
void update_bbox(const std::shared_ptr<scene>& scn, bool do_shapes) {
    if (do_shapes) {
        for (auto shp : scn->shapes) update_bbox(shp);
    }
    scn->bbox = invalid_bbox3f;
    for (auto ist : scn->instances) {
        ist->bbox = transform_bbox(ist->frame, ist->shp->bbox);
        scn->bbox += ist->bbox;
    }
}

// Updates tesselation.
void update_tesselation(
    const std::shared_ptr<subdiv>& sbd, std::shared_ptr<shape> shp) {
    shp->name = sbd->name;
    auto quads_pos = sbd->quads_pos;
    auto quads_texcoord = sbd->quads_texcoord;
    auto quads_color = sbd->quads_color;
    auto pos = sbd->pos;
    auto texcoord = sbd->texcoord;
    auto color = sbd->color;
    subdivide_catmullclark(quads_pos, pos, sbd->level);
    subdivide_catmullclark(quads_texcoord, texcoord, sbd->level, true);
    subdivide_catmullclark(quads_color, color, sbd->level);
    auto norm = std::vector<vec3f>();
    if (sbd->compute_normals) compute_normals(quads_pos, pos, norm);
    auto quads = quads_pos;
    convert_face_varying(quads, shp->pos, shp->norm, shp->texcoord, shp->color,
        quads_pos, quads_pos, quads_texcoord, quads_color, pos, norm, texcoord,
        color);
    shp->triangles = convert_quads_to_triangles(quads);
    update_bbox(shp);
}
void update_tesselation(const std::shared_ptr<scene>& scn) {
    for (auto ist : scn->instances) {
        if (!ist->sbd) continue;
        update_tesselation(ist->sbd, ist->shp);
    }
}

// Update animation transforms
void update_transforms(const std::shared_ptr<animation>& anm, float time,
    const std::string& anim_group) {
    if (anim_group != "" && anim_group != anm->group) return;

    if (!anm->translation.empty()) {
        auto val = vec3f{0, 0, 0};
        switch (anm->type) {
            case animation_type::step:
                val = eval_keyframed_step(anm->times, anm->translation, time);
                break;
            case animation_type::linear:
                val = eval_keyframed_linear(anm->times, anm->translation, time);
                break;
            case animation_type::bezier:
                val = eval_keyframed_bezier(anm->times, anm->translation, time);
                break;
            default: throw std::runtime_error("should not have been here");
        }
        for (auto target : anm->targets) target->translation = val;
    }
    if (!anm->rotation.empty()) {
        auto val = vec4f{0, 0, 0, 1};
        switch (anm->type) {
            case animation_type::step:
                val = eval_keyframed_step(anm->times, anm->rotation, time);
                break;
            case animation_type::linear:
                val = eval_keyframed_linear(anm->times, anm->rotation, time);
                break;
            case animation_type::bezier:
                val = eval_keyframed_bezier(anm->times, anm->rotation, time);
                break;
        }
        for (auto target : anm->targets) target->rotation = val;
    }
    if (!anm->scale.empty()) {
        auto val = vec3f{1, 1, 1};
        switch (anm->type) {
            case animation_type::step:
                val = eval_keyframed_step(anm->times, anm->scale, time);
                break;
            case animation_type::linear:
                val = eval_keyframed_linear(anm->times, anm->scale, time);
                break;
            case animation_type::bezier:
                val = eval_keyframed_bezier(anm->times, anm->scale, time);
                break;
        }
        for (auto target : anm->targets) target->scale = val;
    }
}

// Update node transforms
void update_transforms(const std::shared_ptr<node>& nde,
    const frame3f& parent = identity_frame3f) {
    auto frame = parent * nde->frame * translation_frame(nde->translation) *
                 rotation_frame(nde->rotation) * scaling_frame(nde->scale);
    if (nde->ist) nde->ist->frame = frame;
    if (nde->cam) nde->cam->frame = frame;
    if (nde->env) nde->env->frame = frame;
    for (auto& child : nde->children) update_transforms(child.lock(), frame);
}

// Update node transforms
void update_transforms(const std::shared_ptr<scene>& scn, float time,
    const std::string& anim_group) {
    for (auto& agr : scn->animations) update_transforms(agr, time, anim_group);
    for (auto& nde : scn->nodes) nde->children.clear();
    for (auto& nde : scn->nodes)
        if (nde->parent) nde->parent->children.push_back(nde);
    for (auto& nde : scn->nodes)
        if (!nde->parent) update_transforms(nde);
}

// Compute animation range
vec2f compute_animation_range(
    const std::shared_ptr<scene>& scn, const std::string& anim_group) {
    if (scn->animations.empty()) return zero2f;
    auto range = vec2f{+flt_max, -flt_max};
    for (auto anm : scn->animations) {
        if (anim_group != "" && anm->group != anim_group) continue;
        range.x = min(range.x, anm->times.front());
        range.y = max(range.y, anm->times.back());
    }
    if (range.y < range.x) return zero2f;
    return range;
}

// Update lights.
void update_lights(
    const std::shared_ptr<scene>& scn, bool do_shapes, bool do_environments) {
    if (do_shapes) {
        for (auto shp : scn->shapes) shp->elem_cdf.clear();
    }
    if (do_environments) {
        for (auto env : scn->environments) env->elem_cdf.clear();
    }
    scn->lights.clear();

    for (auto ist : scn->instances) {
        if (!ist->mat || ist->mat->ke == zero3f) continue;
        if (ist->shp->triangles.empty()) continue;
        scn->lights.push_back(ist);
        if (ist->shp->elem_cdf.empty()) update_shape_cdf(ist->shp);
    }

    for (auto env : scn->environments) {
        if (env->ke == zero3f) continue;
        if (env->elem_cdf.empty()) update_environment_cdf(env);
    }
}

// Generate a distribution for sampling a shape uniformly based
// on area/length.
void update_shape_cdf(const std::shared_ptr<shape>& shp) {
    shp->elem_cdf.clear();
    if (!shp->triangles.empty()) {
        shp->elem_cdf = sample_triangles_cdf(shp->triangles, shp->pos);
    } else if (!shp->lines.empty()) {
        shp->elem_cdf = sample_lines_cdf(shp->lines, shp->pos);
    } else if (!shp->pos.empty()) {
        shp->elem_cdf = sample_points_cdf(shp->pos.size());
    } else {
        throw std::runtime_error("empty shape not supported");
    }
}

// Update environment CDF for sampling.
void update_environment_cdf(std::shared_ptr<environment> env) {
    env->elem_cdf.clear();
    auto txt = env->ke_txt;
    if (!txt) return;
    env->elem_cdf.resize(txt->img.width() * txt->img.height());
    if (!txt->img.empty()) {
        for (auto i = 0; i < env->elem_cdf.size(); i++) {
            auto th = (i / txt->img.width() + 0.5f) * pi / txt->img.height();
            env->elem_cdf[i] = max(xyz(txt->img[i])) * sin(th);
            if (i) env->elem_cdf[i] += env->elem_cdf[i - 1];
        }
    } else {
        throw std::runtime_error("empty texture");
    }
}

// Build a shape BVH
void update_bvh(const std::shared_ptr<shape>& shp, bool equalsize) {
    if (!shp->bvh) shp->bvh = std::make_shared<bvh_tree>();
    shp->bvh->pos = shp->pos;
    shp->bvh->radius = shp->radius;
    if (shp->bvh->radius.empty())
        shp->bvh->radius.assign(shp->bvh->radius.size(), 0.001f);
    shp->bvh->points = shp->points;
    shp->bvh->lines = shp->lines;
    shp->bvh->triangles = shp->triangles;
    build_bvh(shp->bvh, equalsize);
}

// Build a scene BVH
void update_bvh(
    const std::shared_ptr<scene>& scn, bool do_shapes, bool equalsize) {
    if (do_shapes) {
        for (auto shp : scn->shapes) update_bvh(shp, equalsize);
    }

    // tree bvh
    if (!scn->bvh) scn->bvh = std::make_shared<bvh_tree>();
    scn->bvh->ist_frames.resize(scn->instances.size());
    scn->bvh->ist_inv_frames.resize(scn->instances.size());
    scn->bvh->ist_bvhs.resize(scn->instances.size());
    for (auto i = 0; i < scn->instances.size(); i++) {
        auto ist = scn->instances[i];
        scn->bvh->ist_frames[i] = ist->frame;
        scn->bvh->ist_inv_frames[i] = inverse(ist->frame, false);
        scn->bvh->ist_bvhs[i] = ist->shp->bvh;
    }
    build_bvh(scn->bvh, equalsize);
}

// Refits a scene BVH
void refit_bvh(const std::shared_ptr<shape>& shp) {
    shp->bvh->pos = shp->pos;
    shp->bvh->radius = shp->radius;
    if (shp->bvh->radius.empty())
        shp->bvh->radius.assign(shp->bvh->radius.size(), 0.001f);
    refit_bvh(shp->bvh);
}

// Refits a scene BVH
void refit_bvh(const std::shared_ptr<scene>& scn, bool do_shapes) {
    if (do_shapes) {
        for (auto shp : scn->shapes) refit_bvh(shp);
    }
    scn->bvh->ist_frames.resize(scn->instances.size());
    scn->bvh->ist_inv_frames.resize(scn->instances.size());
    scn->bvh->ist_bvhs.resize(scn->instances.size());
    for (auto i = 0; i < scn->instances.size(); i++) {
        auto ist = scn->instances[i];
        scn->bvh->ist_frames[i] = ist->frame;
        scn->bvh->ist_inv_frames[i] = inverse(ist->frame);
        scn->bvh->ist_bvhs[i] = ist->shp->bvh;
    }
    refit_bvh(scn->bvh);
}

// Add missing names and resolve duplicated names.
void add_missing_names(const std::shared_ptr<scene>& scn) {
    auto fix_names = [](auto& vals, const std::string& base) {
        auto nmap = std::unordered_map<std::string, int>();
        for (auto val : vals) {
            if (val->name == "") val->name = base;
            if (nmap.find(val->name) == nmap.end()) {
                nmap[val->name] = 0;
            } else {
                nmap[val->name] += 1;
                val->name = val->name + "_" + std::to_string(nmap[val->name]);
            }
        }
    };
    fix_names(scn->cameras, "cam");
    fix_names(scn->shapes, "shp");
    fix_names(scn->textures, "txt");
    fix_names(scn->materials, "mat");
    fix_names(scn->environments, "env");
    fix_names(scn->nodes, "nde");
    fix_names(scn->animations, "anm");
}

// Add missing normals.
void add_missing_normals(const std::shared_ptr<scene>& scn) {
    for (auto shp : scn->shapes) {
        if (shp->norm.empty()) update_normals(shp);
    }
}

// Add missing tangent space if needed.
void add_missing_tangent_space(const std::shared_ptr<scene>& scn) {
    for (auto ist : scn->instances) {
        if (!ist->shp->tangsp.empty() || ist->shp->texcoord.empty()) continue;
        if (!ist->mat || (!ist->mat->norm_txt && !ist->mat->bump_txt)) continue;
        if (!ist->shp->triangles.empty()) {
            if (ist->shp->norm.empty()) update_normals(ist->shp);
            compute_tangent_space(ist->shp->triangles, ist->shp->pos,
                ist->shp->norm, ist->shp->texcoord, ist->shp->tangsp);
        } else {
            throw std::runtime_error("type not supported");
        }
    }
}

// Add missing materials.
void add_missing_materials(const std::shared_ptr<scene>& scn) {
    auto mat = std::shared_ptr<material>();
    for (auto ist : scn->instances) {
        if (ist->mat) continue;
        if (!mat) {
            mat = make_default_material("<default>");
            scn->materials.push_back(mat);
        }
        ist->mat = mat;
    }
}

// Checks for validity of the scene.
std::vector<std::string> validate(
    const std::shared_ptr<scene>& scn, bool skip_textures) {
    auto errs = std::vector<std::string>();
    auto check_names = [&errs](const auto& vals, const std::string& base) {
        auto used = std::unordered_map<std::string, int>();
        for (auto val : vals) used[val->name] += 1;
        for (auto& kv : used) {
            if (kv.first == "")
                errs.push_back("empty " + base + " name");
            else if (kv.second > 1)
                errs.push_back("duplicated " + base + " name " + kv.first);
        }
    };
    auto check_empty_textures =
        [&errs](const std::vector<std::shared_ptr<texture>>& vals) {
            for (auto val : vals) {
                if (val->img.empty())
                    errs.push_back("empty texture " + val->name);
            }
        };

    check_names(scn->cameras, "camera");
    check_names(scn->shapes, "shape");
    check_names(scn->textures, "texture");
    check_names(scn->materials, "material");
    check_names(scn->environments, "environment");
    check_names(scn->nodes, "node");
    check_names(scn->animations, "animation");
    if (!skip_textures) check_empty_textures(scn->textures);

    return errs;
}

// add missing camera
std::shared_ptr<camera> make_bbox_camera(const std::string& name,
    const bbox3f& bbox, float width, float height, float focal) {
    auto bbox_center = (bbox.max + bbox.min) / 2.0f;
    auto bbox_size = bbox.max - bbox.min;
    auto bbox_msize = max(bbox_size.x, max(bbox_size.y, bbox_size.z));
    auto cam = std::make_shared<camera>();
    cam->name = name;
    auto camera_dir = vec3f{1, 0.4f, 1};
    auto from = camera_dir * bbox_msize + bbox_center;
    auto to = bbox_center;
    auto up = vec3f{0, 1, 0};
    cam->frame = lookat_frame(from, to, up);
    cam->ortho = false;
    cam->width = width;
    cam->height = height;
    cam->focal = focal;
    cam->focus = flt_max;
    cam->aperture = 0;
    return cam;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR EVAL AND SAMPLING FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Scene intersection.
scene_intersection intersect_ray(
    const std::shared_ptr<scene>& scn, const ray3f& ray, bool find_any) {
    auto iid = 0;
    auto isec = scene_intersection();
    if (!intersect_bvh(
            scn->bvh, ray, find_any, isec.dist, iid, isec.ei, isec.uv))
        return {};
    isec.ist = scn->instances[iid];
    return isec;
}

// Shape element normal.
vec3f eval_elem_norm(const std::shared_ptr<shape>& shp, int ei) {
    auto norm = zero3f;
    if (!shp->triangles.empty()) {
        auto t = shp->triangles[ei];
        norm = triangle_normal(shp->pos[t.x], shp->pos[t.y], shp->pos[t.z]);
    } else if (!shp->lines.empty()) {
        auto l = shp->lines[ei];
        norm = line_tangent(shp->pos[l.x], shp->pos[l.y]);
    } else {
        norm = {0, 0, 1};
    }
    return norm;
}

// Shape element normal.
vec4f eval_elem_tangsp(const std::shared_ptr<shape>& shp, int ei) {
    auto tangsp = zero4f;
    if (!shp->triangles.empty()) {
        auto t = shp->triangles[ei];
        auto norm =
            triangle_normal(shp->pos[t.x], shp->pos[t.y], shp->pos[t.z]);
        auto txty = std::pair<vec3f, vec3f>();
        if (shp->texcoord.empty()) {
            txty = triangle_tangents_fromuv(shp->pos[t.x], shp->pos[t.y],
                shp->pos[t.z], {0, 0}, {1, 0}, {0, 1});
        } else {
            txty = triangle_tangents_fromuv(shp->pos[t.x], shp->pos[t.y],
                shp->pos[t.z], shp->texcoord[t.x], shp->texcoord[t.y],
                shp->texcoord[t.z]);
        }
        auto tx = txty.first, ty = txty.second;
        tx = orthonormalize(tx, norm);
        auto s = (dot(cross(norm, tx), ty) < 0) ? -1.0f : 1.0f;
        tangsp = {tx.x, tx.y, tx.z, s};
    }
    return tangsp;
}

// Shape value interpolated using barycentric coordinates
template <typename T>
T eval_elem(const std::shared_ptr<shape>& shp, const std::vector<T>& vals,
    int ei, const vec2f& uv) {
    if (vals.empty()) return {};
    if (!shp->triangles.empty()) {
        return interpolate_triangle(vals, shp->triangles[ei], uv);
    } else if (!shp->lines.empty()) {
        return interpolate_line(vals, shp->lines[ei], uv.x);
    } else if (!shp->points.empty()) {
        return vals[shp->points[ei]];
    } else if (!shp->pos.empty()) {
        return vals[ei];
    } else {
        return {};
    }
}

// Shape values interpolated using barycentric coordinates
vec3f eval_pos(const std::shared_ptr<shape>& shp, int ei, const vec2f& uv) {
    return eval_elem(shp, shp->pos, ei, uv);
}
vec3f eval_norm(const std::shared_ptr<shape>& shp, int ei, const vec2f& uv) {
    if (shp->norm.empty()) return eval_elem_norm(shp, ei);
    return normalize(eval_elem(shp, shp->norm, ei, uv));
}
vec2f eval_texcoord(
    const std::shared_ptr<shape>& shp, int ei, const vec2f& uv) {
    if (shp->texcoord.empty()) return uv;
    return eval_elem(shp, shp->texcoord, ei, uv);
}
vec4f eval_color(const std::shared_ptr<shape>& shp, int ei, const vec2f& uv) {
    if (shp->color.empty()) return {1, 1, 1, 1};
    return eval_elem(shp, shp->color, ei, uv);
}
float eval_radius(const std::shared_ptr<shape>& shp, int ei, const vec2f& uv) {
    if (shp->radius.empty()) return 0.001f;
    return eval_elem(shp, shp->radius, ei, uv);
}
vec4f eval_tangsp(const std::shared_ptr<shape>& shp, int ei, const vec2f& uv) {
    if (shp->tangsp.empty()) return eval_elem_tangsp(shp, ei);
    return eval_elem(shp, shp->tangsp, ei, uv);
}
vec3f eval_tangsp(const std::shared_ptr<shape>& shp, int ei, const vec2f& uv,
    bool& left_handed) {
    auto tangsp = (shp->tangsp.empty()) ? eval_elem_tangsp(shp, ei) :
                                          eval_elem(shp, shp->tangsp, ei, uv);
    left_handed = tangsp.w < 0;
    return {tangsp.x, tangsp.y, tangsp.z};
}

// Instance values interpolated using barycentric coordinates.
vec3f eval_pos(const std::shared_ptr<instance>& ist, int ei, const vec2f& uv) {
    return transform_point(ist->frame, eval_pos(ist->shp, ei, uv));
}
vec3f eval_norm(const std::shared_ptr<instance>& ist, int ei, const vec2f& uv) {
    return transform_direction(ist->frame, eval_norm(ist->shp, ei, uv));
}
vec2f eval_texcoord(
    const std::shared_ptr<instance>& ist, int ei, const vec2f& uv) {
    return eval_texcoord(ist->shp, ei, uv);
}
vec4f eval_color(
    const std::shared_ptr<instance>& ist, int ei, const vec2f& uv) {
    return eval_color(ist->shp, ei, uv);
}
float eval_radius(
    const std::shared_ptr<instance>& ist, int ei, const vec2f& uv) {
    return eval_radius(ist->shp, ei, uv);
}
vec3f eval_tangsp(const std::shared_ptr<instance>& ist, int ei, const vec2f& uv,
    bool& left_handed) {
    return transform_direction(
        ist->frame, eval_tangsp(ist->shp, ei, uv, left_handed));
}
// Instance element values.
vec3f eval_elem_norm(const std::shared_ptr<instance>& ist, int ei) {
    return transform_direction(ist->frame, eval_elem_norm(ist->shp, ei));
}
// Shading normals including material perturbations.
vec3f eval_shading_norm(
    const std::shared_ptr<instance>& ist, int ei, const vec2f& uv, vec3f o) {
    if (!ist->shp->triangles.empty()) {
        auto n = eval_norm(ist, ei, uv);
        if (ist->mat && ist->mat->norm_txt) {
            auto texcoord = eval_texcoord(ist, ei, uv);
            auto left_handed = false;
            auto txt = xyz(eval_texture(ist->mat->norm_txt, texcoord));
            txt = txt * 2 - vec3f{1, 1, 1};
            txt.y = -txt.y;  // flip vertical axis to align green with image up
            auto tu = orthonormalize(eval_tangsp(ist, ei, uv, left_handed), n);
            auto tv = normalize(cross(n, tu) * (left_handed ? -1.0f : 1.0f));
            n = normalize(txt.x * tu + txt.y * tv + txt.z * n);
        }
        if (ist->mat && ist->mat->double_sided && dot(n, o) < 0) n = -n;
        return n;
    } else if (!ist->shp->lines.empty()) {
        return orthonormalize(o, eval_norm(ist, ei, uv));
    } else {
        return o;
    }
}

// Environment texture coordinates from the direction.
vec2f eval_texcoord(const std::shared_ptr<environment>& env, vec3f w) {
    auto wl = transform_direction_inverse(env->frame, w);
    auto uv = vec2f{
        atan2(wl.z, wl.x) / (2 * pi), acos(clamp(wl.y, -1.0f, 1.0f)) / pi};
    if (uv.x < 0) uv.x += 1;
    return uv;
}
// Evaluate the environment direction.
vec3f eval_direction(const std::shared_ptr<environment>& env, const vec2f& uv) {
    return transform_direction(
        env->frame, {cos(uv.x * 2 * pi) * sin(uv.y * pi), cos(uv.y * pi),
                        sin(uv.x * 2 * pi) * sin(uv.y * pi)});
}
// Evaluate the environment color.
vec3f eval_environment(const std::shared_ptr<environment>& env, vec3f w) {
    auto ke = env->ke;
    if (env->ke_txt) {
        ke *= xyz(eval_texture(env->ke_txt, eval_texcoord(env, w)));
    }
    return ke;
}

// Evaluate a texture
vec4f eval_texture(const std::shared_ptr<texture>& txt, const vec2f& texcoord) {
    if (!txt || txt->img.empty()) return {1, 1, 1, 1};

    // get image width/height
    auto w = txt->img.width(), h = txt->img.height();

    // get coordinates normalized for tiling
    auto s = 0.0f, t = 0.0f;
    if (txt->clamp) {
        s = clamp(texcoord.x, 0.0f, 1.0f) * w;
        t = clamp(texcoord.y, 0.0f, 1.0f) * h;
    } else {
        s = std::fmod(texcoord.x, 1.0f) * w;
        if (s < 0) s += w;
        t = std::fmod(texcoord.y, 1.0f) * h;
        if (t < 0) t += h;
    }

    // get image coordinates and residuals
    auto i = clamp((int)s, 0, w - 1), j = clamp((int)t, 0, h - 1);
    auto ii = (i + 1) % w, jj = (j + 1) % h;
    auto u = s - i, v = t - j;

    // handle interpolation
    return txt->img[{i, j}] * (1 - u) * (1 - v) +
           txt->img[{i, jj}] * (1 - u) * v + txt->img[{ii, j}] * u * (1 - v) +
           txt->img[{ii, jj}] * u * v;
}

// Set and evaluate camera parameters. Setters take zeros as default values.
float eval_camera_fovy(const std::shared_ptr<camera>& cam) {
    return 2 * std::atan(cam->height / (2 * cam->focal));
}
void set_camera_fovy(
    const std::shared_ptr<camera>& cam, float fovy, float aspect, float width) {
    cam->width = width;
    cam->height = width / aspect;
    cam->focal = cam->height / (2 * std::tan(fovy / 2));
}

// Generates a ray from a camera for image plane coordinate uv and
// the lens coordinates luv.
ray3f eval_camera_ray(
    const std::shared_ptr<camera>& cam, const vec2f& uv, const vec2f& luv) {
    auto dist = cam->focal;
    if (cam->focus < flt_max) {
        dist = cam->focal * cam->focus / (cam->focus - cam->focal);
    }
    auto e = vec3f{luv.x * cam->aperture, luv.y * cam->aperture, 0};
    // auto q = vec3f{cam->width * (uv.x - 0.5f),
    //     cam->height * (uv.y - 0.5f), dist};
    // X flipped for mirror
    auto q =
        vec3f{cam->width * (0.5f - uv.x), cam->height * (uv.y - 0.5f), dist};
    auto ray = make_ray(transform_point(cam->frame, e),
        transform_direction(cam->frame, normalize(e - q)));
    return ray;
}

// Generates a ray from a camera for pixel coordinates `ij`, the
// resolution `res`, the sub-pixel coordinates `puv` and the lens
// coordinates `luv` and the image resolution `res`.
ray3f eval_camera_ray(const std::shared_ptr<camera>& cam, const vec2i& ij,
    const vec2i& imsize, const vec2f& puv, const vec2f& luv) {
    auto uv = vec2f{(ij.x + puv.x) / imsize.x, (ij.y + puv.y) / imsize.y};
    return eval_camera_ray(cam, uv, luv);
}

// Evaluates material parameters.
vec3f eval_emission(
    const std::shared_ptr<instance>& ist, int ei, const vec2f& uv) {
    if (!ist || !ist->mat) return zero3f;
    return ist->mat->ke * xyz(eval_color(ist, ei, uv)) *
           xyz(eval_texture(ist->mat->ke_txt, eval_texcoord(ist, ei, uv)));
}
vec3f eval_diffuse(
    const std::shared_ptr<instance>& ist, int ei, const vec2f& uv) {
    if (!ist || !ist->mat) return zero3f;
    if (!ist->mat->base_metallic) {
        return ist->mat->kd * xyz(eval_color(ist, ei, uv)) *
               xyz(eval_texture(ist->mat->kd_txt, eval_texcoord(ist, ei, uv)));
    } else {
        auto kb =
            ist->mat->kd * xyz(eval_color(ist, ei, uv)) *
            xyz(eval_texture(ist->mat->kd_txt, eval_texcoord(ist, ei, uv)));
        auto km = ist->mat->ks.x *
                  eval_texture(ist->mat->ks_txt, eval_texcoord(ist, ei, uv)).z;
        return kb * (1 - km);
    }
}
vec3f eval_specular(
    const std::shared_ptr<instance>& ist, int ei, const vec2f& uv) {
    if (!ist || !ist->mat) return zero3f;
    if (!ist->mat->base_metallic) {
        return ist->mat->ks * xyz(eval_color(ist, ei, uv)) *
               xyz(eval_texture(ist->mat->ks_txt, eval_texcoord(ist, ei, uv)));
    } else {
        auto kb =
            ist->mat->kd * xyz(eval_color(ist, ei, uv)) *
            xyz(eval_texture(ist->mat->kd_txt, eval_texcoord(ist, ei, uv)));
        auto km = ist->mat->ks.x *
                  eval_texture(ist->mat->ks_txt, eval_texcoord(ist, ei, uv)).z;
        return kb * km + vec3f{0.04f, 0.04f, 0.04f} * (1 - km);
    }
}
float eval_roughness(
    const std::shared_ptr<instance>& ist, int ei, const vec2f& uv) {
    if (!ist || !ist->mat) return 1;
    if (!ist->mat->base_metallic) {
        if (!ist->mat->gltf_textures) {
            auto rs =
                ist->mat->rs *
                eval_texture(ist->mat->rs_txt, eval_texcoord(ist, ei, uv)).x;
            return rs * rs;
        } else {
            auto gs =
                (1 - ist->mat->rs) *
                eval_texture(ist->mat->rs_txt, eval_texcoord(ist, ei, uv)).w;
            auto rs = 1 - gs;
            return rs * rs;
        }
    } else {
        auto rs = ist->mat->rs *
                  eval_texture(ist->mat->rs_txt, eval_texcoord(ist, ei, uv)).y;
        return rs * rs;
    }
}
vec3f eval_transmission(
    const std::shared_ptr<instance>& ist, int ei, const vec2f& uv) {
    if (!ist || !ist->mat) return zero3f;
    return ist->mat->kt * xyz(eval_color(ist, ei, uv)) *
           xyz(eval_texture(ist->mat->kt_txt, eval_texcoord(ist, ei, uv)));
}
float eval_opacity(
    const std::shared_ptr<instance>& ist, int ei, const vec2f& uv) {
    if (!ist || !ist->mat) return 1;
    return ist->mat->op * eval_color(ist->shp, ei, uv).w *
           eval_texture(ist->mat->kd_txt, eval_texcoord(ist, ei, uv)).w *
           eval_texture(ist->mat->op_txt, eval_texcoord(ist, ei, uv)).x;
}

// Evaluates the bsdf at a location.
bsdf eval_bsdf(const std::shared_ptr<instance>& ist, int ei, const vec2f& uv) {
    auto f = bsdf();
    f.kd = eval_diffuse(ist, ei, uv);
    f.ks = eval_specular(ist, ei, uv);
    f.kt = eval_transmission(ist, ei, uv);
    f.rs = eval_roughness(ist, ei, uv);
    f.refract = (ist && ist->mat) ? ist->mat->refract : false;
    if (f.kd != zero3f) {
        f.rs = clamp(f.rs, 0.03f * 0.03f, 1.0f);
    } else if (f.rs <= 0.03f * 0.03f)
        f.rs = 0;
    return f;
}
bool is_delta_bsdf(const bsdf& f) { return f.rs == 0 && f.kd == zero3f; }

// Sample a shape based on a distribution.
std::pair<int, vec2f> sample_shape(
    const std::shared_ptr<shape>& shp, float re, const vec2f& ruv) {
    // TODO: implement sampling without cdf
    if (shp->elem_cdf.empty()) return {};
    if (!shp->triangles.empty()) {
        return sample_triangles(shp->elem_cdf, re, ruv);
    } else if (!shp->lines.empty()) {
        return {sample_lines(shp->elem_cdf, re, ruv.x).first, ruv};
    } else if (!shp->pos.empty()) {
        return {sample_points(shp->elem_cdf, re), ruv};
    } else {
        return {0, zero2f};
    }
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SCENE UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Merge scene into one another
void merge_into(const std::shared_ptr<scene>& merge_into,
    std::shared_ptr<scene> merge_from) {
    auto merge = [](auto& v1, auto& v2) {
        v1.insert(v1.end(), v2.begin(), v2.end());
        v2.clear();
    };
    merge(merge_into->cameras, merge_from->cameras);
    merge(merge_into->textures, merge_from->textures);
    merge(merge_into->materials, merge_from->materials);
    merge(merge_into->shapes, merge_from->shapes);
    merge(merge_into->environments, merge_from->environments);
    merge(merge_into->nodes, merge_from->nodes);
    merge(merge_into->animations, merge_from->animations);
}

void print_stats(const std::shared_ptr<scene>& scn) {
    uint64_t num_cameras = 0;
    uint64_t num_shape_groups = 0;
    uint64_t num_shapes = 0;
    uint64_t num_instances = 0;
    uint64_t num_materials = 0;
    uint64_t num_textures = 0;
    uint64_t num_environments = 0;
    uint64_t num_nodes = 0;
    uint64_t num_animations = 0;

    uint64_t elem_lines = 0;
    uint64_t elem_triangles = 0;
    uint64_t vert_pos = 0;
    uint64_t vert_norm = 0;
    uint64_t vert_texcoord = 0;
    uint64_t vert_color = 0;
    uint64_t vert_radius = 0;
    uint64_t vert_tangsp = 0;

    uint64_t texel_hdr = 0;
    uint64_t texel_ldr = 0;

    uint64_t memory_imgs = 0;
    uint64_t memory_elems = 0;
    uint64_t memory_verts = 0;

    update_bbox(scn, true);
    auto bbox = scn->bbox;

    num_cameras = scn->cameras.size();
    num_shapes = scn->shapes.size();
    num_materials = scn->materials.size();
    num_textures = scn->textures.size();
    num_environments = scn->environments.size();
    num_instances = scn->instances.size();
    num_nodes = scn->nodes.size();
    num_animations = scn->animations.size();

    for (auto shp : scn->shapes) {
        elem_lines += shp->lines.size();
        elem_triangles += shp->triangles.size();
        vert_pos += shp->pos.size();
        vert_norm += shp->norm.size();
        vert_texcoord += shp->texcoord.size();
        vert_color += shp->color.size();
        vert_radius += shp->radius.size();
        vert_tangsp += shp->tangsp.size();
    }
    memory_elems = elem_lines * sizeof(vec2i) + elem_triangles * sizeof(vec3i);
    memory_verts = vert_pos * sizeof(vec3f) + vert_norm * sizeof(vec3f) +
                   vert_texcoord * sizeof(vec3f) + vert_color * sizeof(vec4f) +
                   vert_tangsp * sizeof(vec4f) + vert_radius * sizeof(float);

    for (auto txt : scn->textures) { texel_hdr += txt->img.size(); }
    memory_imgs = texel_hdr * sizeof(vec4f) + texel_ldr * sizeof(vec4b);

    std::cout << "num_cameras: " << num_cameras << "\n";
    std::cout << "num_shape_groups: " << num_shape_groups << "\n";
    std::cout << "num_shapes: " << num_shapes << "\n";
    std::cout << "num_instances: " << num_instances << "\n";
    std::cout << "num_materials: " << num_materials << "\n";
    std::cout << "num_textures: " << num_textures << "\n";
    std::cout << "num_environments: " << num_environments << "\n";
    std::cout << "num_nodes: " << num_nodes << "\n";
    std::cout << "num_animations: " << num_animations << "\n";
    std::cout << "elem_lines: " << elem_lines << "\n";
    std::cout << "elem_triangles: " << elem_triangles << "\n";
    std::cout << "vert_pos: " << vert_pos << "\n";
    std::cout << "vert_norm: " << vert_norm << "\n";
    std::cout << "vert_texcoord: " << vert_texcoord << "\n";
    std::cout << "vert_color: " << vert_color << "\n";
    std::cout << "vert_radius: " << vert_radius << "\n";
    std::cout << "vert_tangsp: " << vert_tangsp << "\n";
    std::cout << "texel_hdr: " << texel_hdr << "\n";
    std::cout << "texel_ldr: " << texel_ldr << "\n";
    std::cout << "memory_imgs: " << memory_imgs << "\n";
    std::cout << "memory_elems: " << memory_elems << "\n";
    std::cout << "memory_verts: " << memory_verts << "\n";
    std::cout << "bbox: " << bbox << "\n";
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PATH TRACING
// -----------------------------------------------------------------------------
namespace ygl {

// Intersect a scene handling opacity.
scene_intersection intersect_ray_cutout(const std::shared_ptr<scene>& scn,
    const ray3f& ray_, rng_state& rng, int nbounces) {
    auto ray = ray_;
    for (auto b = 0; b < nbounces; b++) {
        auto isec = intersect_ray(scn, ray);
        if (!isec.ist) return isec;
        auto op = eval_opacity(isec.ist, isec.ei, isec.uv);
        if (op > 0.999f) return isec;
        if (rand1f(rng) < op) return isec;
        ray = make_ray(eval_pos(isec.ist, isec.ei, isec.uv), ray.d);
    }
    return {};
}

// Check if we are near the mirror direction.
inline bool check_near_mirror(const vec3f& n, const vec3f& o, const vec3f& i) {
    return fabs(dot(i, normalize(n * 2.0f * dot(o, n) - o)) - 1) < 0.001f;
}

// Schlick approximation of the Fresnel term
vec3f fresnel_schlick(const vec3f& ks, const vec3f& h, const vec3f& i) {
    if (ks == zero3f) return zero3f;
    return ks + (vec3f{1, 1, 1} - ks) *
                    pow(clamp(1.0f - fabs(dot(h, i)), 0.0f, 1.0f), 5.0f);
}
vec3f fresnel_schlick(
    const vec3f& ks, const vec3f& h, const vec3f& i, float rs) {
    if (ks == zero3f) return zero3f;
    auto fks = fresnel_schlick(ks, fabs(dot(h, i)));
    return ks + (fks - ks) * (1 - sqrt(clamp(rs, 0.0f, 1.0f)));
}

// Evaluates the GGX distribution and geometric term
float eval_ggx_dist(float rs, const vec3f& n, const vec3f& h) {
    auto di = (dot(n, h) * dot(n, h)) * (rs * rs - 1) + 1;
    return rs * rs / (pi * di * di);
}
float eval_ggx_sm(float rs, const vec3f& n, const vec3f& o, const vec3f& i) {
#if 0
    // evaluate G from Heitz
    auto lambda_o = (-1 + sqrt(1 + alpha2 * (1 - ndo * ndo) / (ndo * ndo))) / 2;
    auto lambda_i = (-1 + sqrt(1 + alpha2 * (1 - ndi * ndi) / (ndi * ndi))) / 2;
    auto g = 1 / (1 + lambda_o + lambda_i);
#else
    auto Go = (2 * fabs(dot(n, o))) /
              (fabs(dot(n, o)) +
                  sqrt(rs * rs + (1 - rs * rs) * dot(n, o) * dot(n, o)));
    auto Gi = (2 * fabs(dot(n, i))) /
              (fabs(dot(n, i)) +
                  sqrt(rs * rs + (1 - rs * rs) * dot(n, i) * dot(n, i)));
    return Go * Gi;
#endif
}

// Evaluates the BRDF scaled by the cosine of the incoming direction.
// - ggx from [Heitz 2014] and [Walter 2007] and [Lagarde 2014]
// "Understanding the Masking-Shadowing Function in Microfacet-Based
// BRDFs" http://jcgt.org/published/0003/02/03/
// - "Microfacet Models for Refraction through Rough Surfaces" EGSR 07
// https://www.cs.cornell.edu/~srm/publications/EGSR07-btdf.pdf
vec3f eval_bsdf(const bsdf& f, const vec3f& n, const vec3f& o, const vec3f& i) {
    if (is_delta_bsdf(f)) return zero3f;
    auto bsdf = zero3f;

    // diffuse
    if (f.kd != zero3f && dot(n, o) * dot(n, i) > 0) {
        auto h = normalize(i + o);
        auto F = fresnel_schlick(f.ks, h, o);
        bsdf += f.kd * (vec3f{1, 1, 1} - F) / pi;
    }

    // specular
    if (f.ks != zero3f && dot(n, o) * dot(n, i) > 0) {
        auto h = normalize(i + o);
        auto F = fresnel_schlick(f.ks, h, o);
        auto D = eval_ggx_dist(f.rs, n, h);
        auto G = eval_ggx_sm(f.rs, n, o, i);
        bsdf += F * D * G / (4 * fabs(dot(n, o)) * fabs(dot(n, i)));
    }

    // transmission (thin sheet)
    if (f.kt != zero3f && dot(n, o) * dot(n, i) < 0) {
        auto ir = (dot(n, o) >= 0) ? reflect(-i, n) : reflect(-i, -n);
        auto h = normalize(ir + o);
        auto F = fresnel_schlick(f.ks, h, o);
        auto D = eval_ggx_dist(f.rs, n, h);
        auto G = eval_ggx_sm(f.rs, n, o, ir);
        bsdf += f.kt * (vec3f{1, 1, 1} - F) * D * G /
                (4 * fabs(dot(n, o)) * fabs(dot(n, ir)));
    }

    return bsdf;
}

// Evaluates the BRDF assuming that it is called only from the directions
// generated by sample_brdf.
vec3f eval_delta_brdf(
    const bsdf& f, const vec3f& n, const vec3f& o, const vec3f& i) {
    if (!is_delta_bsdf(f)) return zero3f;
    auto bsdf = zero3f;

    // specular
    if (f.ks != zero3f && dot(n, o) * dot(n, i) > 0) {
        auto F = fresnel_schlick(f.ks, n, o);
        bsdf += F / fabs(dot(n, i));
    }

    // transmission (thin sheet)
    if (f.kt != zero3f && dot(n, o) * dot(n, i) < 0) {
        auto F = fresnel_schlick(f.ks, n, o);
        bsdf += f.kt * (vec3f{1, 1, 1} - F) / fabs(dot(n, i));
    }

    return bsdf;
}

// Picks a direction based on the BRDF
vec3f sample_brdf(
    const bsdf& f, const vec3f& n, const vec3f& o, float rnl, const vec2f& rn) {
    if (is_delta_bsdf(f)) return zero3f;
    auto F = fresnel_schlick(f.ks, n, o);
    auto prob = vec3f{max(f.kd * (vec3f{1, 1, 1} - F)), max(F),
        max(f.kt * (vec3f{1, 1, 1} - F))};
    if (prob == zero3f) return zero3f;
    prob /= prob.x + prob.y + prob.z;

    // sample according to diffuse
    if (f.kd != zero3f && rnl < prob.x) {
        auto rz = sqrtf(rn.y), rr = sqrtf(1 - rz * rz), rphi = 2 * pi * rn.x;
        auto il = vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
        auto fp = dot(n, o) >= 0 ? make_frame_fromz(zero3f, n) :
                                   make_frame_fromz(zero3f, -n);
        return transform_direction(fp, il);
    }
    // sample according to specular GGX
    else if (f.ks != zero3f && rnl < prob.x + prob.y) {
        auto hl = sample_ggx(f.rs, rn);
        auto fp = dot(n, o) >= 0 ? make_frame_fromz(zero3f, n) :
                                   make_frame_fromz(zero3f, -n);
        auto h = transform_direction(fp, hl);
        return reflect(o, h);
    }
    // transmission hack
    else if (f.kt != zero3f && rnl < prob.x + prob.y + prob.z) {
        auto hl = sample_ggx(f.rs, rn);
        auto fp = dot(n, o) >= 0 ? make_frame_fromz(zero3f, n) :
                                   make_frame_fromz(zero3f, -n);
        auto h = transform_direction(fp, hl);
        auto ir = reflect(o, h);
        return dot(n, o) >= 0 ? reflect(-ir, -n) : reflect(-ir, n);
    } else {
        return zero3f;
    }
}

// Picks a direction based on the BRDF
vec3f sample_delta_brdf(
    const bsdf& f, const vec3f& n, const vec3f& o, float rnl, const vec2f& rn) {
    if (!is_delta_bsdf(f)) return zero3f;
    auto F = fresnel_schlick(f.ks, n, o);
    auto prob = vec3f{0, max(F), max(f.kt * (vec3f{1, 1, 1} - F))};
    if (prob == zero3f) return zero3f;
    prob /= prob.x + prob.y + prob.z;

    // sample according to specular mirror
    if (f.ks != zero3f && rnl < prob.x + prob.y) {
        return reflect(o, dot(n, o) >= 0 ? n : -n);
    }
    // sample according to transmission
    else if (f.kt != zero3f && !f.refract && rnl < prob.x + prob.y + prob.z) {
        return -o;
    }
    // sample according to transmission
    else if (f.kt != zero3f && f.refract && rnl < prob.x + prob.y + prob.z) {
        if (dot(n, o) >= 0) {
            return refract(o, n, 1 / specular_to_eta(f.ks));
        } else {
            return refract(o, -n, specular_to_eta(f.ks));
        }
    }
    // no sampling
    else {
        return zero3f;
    }
}

// Compute the weight for sampling the BRDF
float sample_brdf_pdf(
    const bsdf& f, const vec3f& n, const vec3f& o, const vec3f& i) {
    if (is_delta_bsdf(f)) return 0;
    auto F = fresnel_schlick(f.ks, n, o);
    auto prob = vec3f{max(f.kd * (vec3f{1, 1, 1} - F)), max(F),
        max(f.kt * (vec3f{1, 1, 1} - F))};
    if (prob == zero3f) return 0;
    prob /= prob.x + prob.y + prob.z;

    auto pdf = 0.0f;

    if (f.kd != zero3f && dot(n, o) * dot(n, i) > 0) {
        pdf += prob.x * fabs(dot(n, i)) / pi;
    }
    if (f.ks != zero3f && dot(n, o) * dot(n, i) > 0) {
        auto h = normalize(i + o);
        auto d = sample_ggx_pdf(f.rs, fabs(dot(n, h)));
        pdf += prob.y * d / (4 * fabs(dot(o, h)));
    }
    if (f.kt != zero3f && dot(n, o) * dot(n, i) < 0) {
        auto ir = (dot(n, o) >= 0) ? reflect(-i, n) : reflect(-i, -n);
        auto h = normalize(ir + o);
        auto d = sample_ggx_pdf(f.rs, fabs(dot(n, h)));
        pdf += prob.z * d / (4 * fabs(dot(o, h)));
    }

    return pdf;
}

// Compute the weight for sampling the BRDF
float sample_delta_brdf_pdf(
    const bsdf& f, const vec3f& n, const vec3f& o, const vec3f& i) {
    if (!is_delta_bsdf(f)) return 0;
    auto F = fresnel_schlick(f.ks, n, o);
    auto prob = vec3f{0, max(F), max(f.kt * (vec3f{1, 1, 1} - F))};
    if (prob == zero3f) return 0;
    prob /= prob.x + prob.y + prob.z;

    auto pdf = 0.0f;

    if (f.ks != zero3f && dot(n, o) * dot(n, i) > 0) { return prob.y; }
    if (f.kt != zero3f && dot(n, o) * dot(n, i) < 0) { return prob.z; }

    return pdf;
}

// Sample pdf for an environment.
float sample_environment_pdf(
    const std::shared_ptr<environment> env, const vec3f& i) {
    auto txt = env->ke_txt;
    if (!env->elem_cdf.empty() && txt) {
        auto texcoord = eval_texcoord(env, i);
        auto i = (int)(texcoord.x * txt->img.width());
        auto j = (int)(texcoord.y * txt->img.height());
        auto idx = j * txt->img.width() + i;
        auto prob =
            sample_discrete_pdf(env->elem_cdf, idx) / env->elem_cdf.back();
        auto angle = (2 * pi / txt->img.width()) * (pi / txt->img.height()) *
                     sin(pi * (j + 0.5f) / txt->img.height());
        return prob / angle;
    } else {
        return 1 / (4 * pi);
    }
}

// Picks a point on an environment.
vec3f sample_environment(
    const std::shared_ptr<environment> env, float rel, const vec2f& ruv) {
    auto txt = env->ke_txt;
    if (!env->elem_cdf.empty() && txt) {
        auto idx = sample_discrete(env->elem_cdf, rel);
        auto u = (idx % txt->img.width() + 0.5f) / txt->img.width();
        auto v = (idx / txt->img.width() + 0.5f) / txt->img.height();
        return eval_direction(env, {u, v});
    } else {
        return sample_sphere(ruv);
    }
}

// Picks a point on a light.
vec3f sample_light(const std::shared_ptr<instance> ist, const vec3f& p,
    float rel, const vec2f& ruv) {
    auto sample = sample_shape(ist->shp, rel, ruv);
    return normalize(eval_pos(ist, sample.first, sample.second) - p);
}

// Sample pdf for a light point.
float sample_light_pdf(const std::shared_ptr<instance> ist, const vec3f& p,
    const vec3f& i, const vec3f& lp, const vec3f& ln) {
    if (ist->mat->ke == zero3f) return 0;
    // prob triangle * area triangle = area triangle mesh
    auto area = ist->shp->elem_cdf.back();
    return dot(lp - p, lp - p) / (fabs(dot(ln, i)) * area);
}

// Test occlusion.
vec3f eval_transmission(const std::shared_ptr<scene>& scn, const vec3f& from,
    const vec3f& to, int nbounces) {
    auto weight = vec3f{1, 1, 1};
    auto p = from;
    for (auto bounce = 0; bounce < nbounces; bounce++) {
        auto ray = make_segment(p, to);
        auto isec = intersect_ray(scn, ray);
        if (!isec.ist) break;
        auto f = eval_bsdf(isec.ist, isec.ei, isec.uv);
        auto op = eval_opacity(isec.ist, isec.ei, isec.uv);
        weight *= f.kt + vec3f{1 - op, 1 - op, 1 - op};
        if (weight == zero3f) break;
        p = eval_pos(isec.ist, isec.ei, isec.uv);
    }
    return weight;
}

// Recursive path tracing.
vec3f trace_path(const std::shared_ptr<scene>& scn, const ray3f& ray_,
    rng_state& rng, int nbounces, bool* hit) {
    if (scn->lights.empty() && scn->environments.empty()) return zero3f;

    // initialize
    auto l = zero3f;
    auto weight = vec3f{1, 1, 1};
    auto emission = true;
    auto ray = ray_;

    // trace  path
    for (auto bounce = 0; bounce < nbounces; bounce++) {
        // intersect ray
        auto isec = intersect_ray_cutout(scn, ray, rng, nbounces);
        if (!isec.ist) {
            if (emission) {
                for (auto env : scn->environments)
                    l += weight * eval_environment(env, ray.d);
            }
            break;
        }
        if (hit) *hit = true;

        // point
        auto o = -ray.d;
        auto p = eval_pos(isec.ist, isec.ei, isec.uv);
        auto n = eval_shading_norm(isec.ist, isec.ei, isec.uv, o);
        auto f = eval_bsdf(isec.ist, isec.ei, isec.uv);

        // emission
        if (emission) l += weight * eval_emission(isec.ist, isec.ei, isec.uv);

        // early exit and russian roulette
        if (f.kd + f.ks + f.kt == zero3f || bounce >= nbounces - 1) break;
        if (bounce > 2) {
            auto rrprob = 1.0f - min(max(weight), 0.95f);
            if (rand1f(rng) < rrprob) break;
            weight *= 1 / (1 - rrprob);
        }

        // direct
        if (!is_delta_bsdf(f) &&
            (!scn->lights.empty() || !scn->environments.empty())) {
            auto i = zero3f;
            auto nlights = (int)(scn->lights.size() + scn->environments.size());
            if (rand1f(rng) < 0.5f) {
                auto idx = sample_index(nlights, rand1f(rng));
                if (idx < scn->lights.size()) {
                    auto lgt = scn->lights[idx];
                    i = sample_light(lgt, p, rand1f(rng), rand2f(rng));
                } else {
                    auto env = scn->environments[idx - scn->lights.size()];
                    i = sample_environment(env, rand1f(rng), rand2f(rng));
                }
            } else {
                i = sample_brdf(f, n, o, rand1f(rng), rand2f(rng));
            }
            auto isec =
                intersect_ray_cutout(scn, make_ray(p, i), rng, nbounces);
            auto pdf = 0.5f * sample_brdf_pdf(f, n, o, i);
            auto le = zero3f;
            if (isec.ist) {
                auto lp = eval_pos(isec.ist, isec.ei, isec.uv);
                auto ln = eval_shading_norm(isec.ist, isec.ei, isec.uv, -i);
                pdf += 0.5f * sample_light_pdf(isec.ist, p, i, lp, ln) *
                       sample_index_pdf<float>(nlights);
                le += eval_emission(isec.ist, isec.ei, isec.uv);
            } else {
                for (auto env : scn->environments) {
                    pdf += 0.5f * sample_environment_pdf(env, i) *
                           sample_index_pdf<float>(nlights);
                    le += eval_environment(env, i);
                }
            }
            auto brdfcos = eval_bsdf(f, n, o, i) * fabs(dot(n, i));
            if (pdf != 0) l += weight * le * brdfcos / pdf;
        }

        // continue path
        auto i = zero3f, brdfcos = zero3f;
        auto pdf = 0.0f;
        if (!is_delta_bsdf(f)) {
            i = sample_brdf(f, n, o, rand1f(rng), rand2f(rng));
            brdfcos = eval_bsdf(f, n, o, i) * fabs(dot(n, i));
            pdf = sample_brdf_pdf(f, n, o, i);
        } else {
            i = sample_delta_brdf(f, n, o, rand1f(rng), rand2f(rng));
            brdfcos = eval_delta_brdf(f, n, o, i) * fabs(dot(n, i));
            pdf = sample_delta_brdf_pdf(f, n, o, i);
        }

        // accumulate weight
        if (pdf == 0) break;
        weight *= brdfcos / pdf;
        if (weight == zero3f) break;

        // setup next ray
        ray = make_ray(p, i);
        emission = is_delta_bsdf(f);
    }

    return l;
}

// Recursive path tracing.
vec3f trace_path_naive(const std::shared_ptr<scene>& scn, const ray3f& ray_,
    rng_state& rng, int nbounces, bool* hit) {
    if (scn->lights.empty() && scn->environments.empty()) return zero3f;

    // initialize
    auto l = zero3f;
    auto weight = vec3f{1, 1, 1};
    auto ray = ray_;

    // trace  path
    for (auto bounce = 0; bounce < nbounces; bounce++) {
        // intersect ray
        auto isec = intersect_ray_cutout(scn, ray, rng, nbounces);
        if (!isec.ist) {
            for (auto env : scn->environments)
                l += weight * eval_environment(env, ray.d);
            break;
        }
        if (hit) *hit = true;

        // point
        auto o = -ray.d;
        auto p = eval_pos(isec.ist, isec.ei, isec.uv);
        auto n = eval_shading_norm(isec.ist, isec.ei, isec.uv, o);
        auto f = eval_bsdf(isec.ist, isec.ei, isec.uv);

        // emission
        l += weight * eval_emission(isec.ist, isec.ei, isec.uv);

        // early exit and russian roulette
        if (f.kd + f.ks + f.kt == zero3f || bounce >= nbounces - 1) break;
        if (bounce > 2) {
            auto rrprob = 1.0f - min(max(weight), 0.95f);
            if (rand1f(rng) < rrprob) break;
            weight *= 1 / (1 - rrprob);
        }

        // continue path
        auto i = zero3f, brdfcos = zero3f;
        auto pdf = 0.0f;
        if (!is_delta_bsdf(f)) {
            i = sample_brdf(f, n, o, rand1f(rng), rand2f(rng));
            brdfcos = eval_bsdf(f, n, o, i) * fabs(dot(n, i));
            pdf = sample_brdf_pdf(f, n, o, i);
        } else {
            i = sample_delta_brdf(f, n, o, rand1f(rng), rand2f(rng));
            brdfcos = eval_delta_brdf(f, n, o, i) * fabs(dot(n, i));
            pdf = sample_delta_brdf_pdf(f, n, o, i);
        }

        // accumulate weight
        if (pdf == 0) break;
        weight *= brdfcos / pdf;
        if (weight == zero3f) break;

        // setup next ray
        ray = make_ray(p, i);
    }

    return l;
}

// Recursive path tracing.
vec3f trace_path_nomis(const std::shared_ptr<scene>& scn, const ray3f& ray_,
    rng_state& rng, int nbounces, bool* hit) {
    if (scn->lights.empty() && scn->environments.empty()) return zero3f;

    // initialize
    auto l = zero3f;
    auto weight = vec3f{1, 1, 1};
    auto emission = true;
    auto ray = ray_;

    // trace  path
    for (auto bounce = 0; bounce < nbounces; bounce++) {
        // intersect ray
        auto isec = intersect_ray_cutout(scn, ray, rng, nbounces);
        if (!isec.ist) {
            for (auto env : scn->environments)
                l += weight * eval_environment(env, ray.d);
            break;
        }
        if (hit) *hit = true;

        // point
        auto o = -ray.d;
        auto p = eval_pos(isec.ist, isec.ei, isec.uv);
        auto n = eval_shading_norm(isec.ist, isec.ei, isec.uv, o);
        auto f = eval_bsdf(isec.ist, isec.ei, isec.uv);

        // emission
        if (emission) l += weight * eval_emission(isec.ist, isec.ei, isec.uv);

        // early exit and russian roulette
        if (f.kd + f.ks + f.kt == zero3f || bounce >= nbounces - 1) break;
        if (bounce > 2) {
            auto rrprob = 1.0f - min(max(weight), 0.95f);
            if (rand1f(rng) < rrprob) break;
            weight *= 1 / (1 - rrprob);
        }

        // direct
        if (!is_delta_bsdf(f) && !scn->lights.empty()) {
            auto lgt =
                scn->lights[sample_index(scn->lights.size(), rand1f(rng))];
            auto i = sample_light(lgt, p, rand1f(rng), rand2f(rng));
            auto isec =
                intersect_ray_cutout(scn, make_ray(p, i), rng, nbounces);
            if (isec.ist && isec.ist->mat->ke != zero3f) {
                auto lp = eval_pos(isec.ist, isec.ei, isec.uv);
                auto ln = eval_shading_norm(isec.ist, isec.ei, isec.uv, -i);
                auto pdf = sample_light_pdf(isec.ist, p, i, lp, ln) *
                           sample_index_pdf<float>(scn->lights.size());
                auto le = eval_emission(isec.ist, isec.ei, isec.uv);
                auto brdfcos = eval_bsdf(f, n, o, i) * fabs(dot(n, i));
                if (pdf != 0) l += weight * le * brdfcos / pdf;
            }
        }

        // continue path
        auto i = zero3f, brdfcos = zero3f;
        auto pdf = 0.0f;
        if (!is_delta_bsdf(f)) {
            i = sample_brdf(f, n, o, rand1f(rng), rand2f(rng));
            brdfcos = eval_bsdf(f, n, o, i) * fabs(dot(n, i));
            pdf = sample_brdf_pdf(f, n, o, i);
        } else {
            i = sample_delta_brdf(f, n, o, rand1f(rng), rand2f(rng));
            brdfcos = eval_delta_brdf(f, n, o, i) * fabs(dot(n, i));
            pdf = sample_delta_brdf_pdf(f, n, o, i);
        }

        // accumulate weight
        if (pdf == 0) break;
        weight *= brdfcos / pdf;
        if (weight == zero3f) break;

        // setup next ray
        ray = make_ray(p, i);
        emission = is_delta_bsdf(f);
    }

    return l;
}

// Direct illumination.
vec3f trace_direct(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit) {
    if (scn->lights.empty() && scn->environments.empty()) return zero3f;

    // intersect scene
    auto isec = intersect_ray(scn, ray);
    auto l = zero3f;

    // handle environment
    if (!isec.ist) {
        for (auto env : scn->environments) l += eval_environment(env, ray.d);
        return l;
    }
    if (hit) *hit = true;

    // point
    auto o = -ray.d;
    auto p = eval_pos(isec.ist, isec.ei, isec.uv);
    auto n = eval_shading_norm(isec.ist, isec.ei, isec.uv, o);
    auto f = eval_bsdf(isec.ist, isec.ei, isec.uv);

    // emission
    l += eval_emission(isec.ist, isec.ei, isec.uv);

    // direct lights
    for (auto lgt : scn->lights) {
        auto i = zero3f;
        if (rand1f(rng) < 0.5f) {
            i = sample_light(lgt, p, rand1f(rng), rand2f(rng));
        } else {
            i = sample_brdf(f, n, o, rand1f(rng), rand2f(rng));
        }
        auto isec = intersect_ray(scn, make_ray(p, i));
        if (lgt != isec.ist) continue;
        auto lp = eval_pos(isec.ist, isec.ei, isec.uv);
        auto ln = eval_shading_norm(isec.ist, isec.ei, isec.uv, -i);
        auto pdf = 0.5f * sample_light_pdf(isec.ist, p, i, lp, ln) +
                   0.5f * sample_brdf_pdf(f, n, o, i);
        auto le = eval_emission(isec.ist, isec.ei, isec.uv);
        auto brdfcos = eval_bsdf(f, n, o, i) * fabs(dot(n, i));
        if (pdf != 0) l += le * brdfcos / pdf;
    }

    // direct environments
    for (auto env : scn->environments) {
        auto i = zero3f;
        if (rand1f(rng) < 0.5f) {
            i = sample_environment(env, rand1f(rng), rand2f(rng));
        } else {
            i = sample_brdf(f, n, o, rand1f(rng), rand2f(rng));
        }
        auto isec = intersect_ray(scn, make_ray(p, i));
        if (isec.ist) continue;
        auto pdf = 0.5f * sample_environment_pdf(env, i) +
                   0.5f * sample_brdf_pdf(f, n, o, i);
        auto le = eval_environment(env, i);
        auto brdfcos = eval_bsdf(f, n, o, i) * fabs(dot(n, i));
        if (pdf != 0) l += le * brdfcos / pdf;
    }

    // exit if needed
    if (nbounces <= 0) return l;

    // reflection
    if (f.ks != zero3f && !f.rs) {
        auto i = reflect(o, n);
        l += f.ks * trace_direct(scn, make_ray(p, i), rng, nbounces - 1, hit);
    }

    // refraction
    if (f.kt != zero3f) {
        l += f.kt * trace_direct(scn, make_ray(p, -o), rng, nbounces - 1, hit);
    }

    // opacity
    auto op = eval_opacity(isec.ist, isec.ei, isec.uv);
    if (op != 1) {
        l = op * l + (1 - op) * trace_direct(scn, make_ray(p, -o), rng,
                                    nbounces - 1, hit);
    }

    // done
    return l;
}

// Direct illumination.
vec3f trace_direct_nomis(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit) {
    if (scn->lights.empty() && scn->environments.empty()) return zero3f;

    // intersect scene
    auto isec = intersect_ray(scn, ray);
    auto l = zero3f;

    // handle environment
    if (!isec.ist) {
        for (auto env : scn->environments) l += eval_environment(env, ray.d);
        return l;
    }
    if (hit) *hit = true;

    // point
    auto o = -ray.d;
    auto p = eval_pos(isec.ist, isec.ei, isec.uv);
    auto n = eval_shading_norm(isec.ist, isec.ei, isec.uv, o);
    auto f = eval_bsdf(isec.ist, isec.ei, isec.uv);

    // emission
    l += eval_emission(isec.ist, isec.ei, isec.uv);

    // direct lights
    for (auto lgt : scn->lights) {
        auto i = sample_light(lgt, p, rand1f(rng), rand2f(rng));
        auto isec = intersect_ray(scn, make_ray(p, i));
        if (lgt != isec.ist) continue;
        auto lp = eval_pos(isec.ist, isec.ei, isec.uv);
        auto ln = eval_shading_norm(isec.ist, isec.ei, isec.uv, -i);
        auto pdf = sample_light_pdf(isec.ist, p, i, lp, ln);
        auto le = eval_emission(isec.ist, isec.ei, isec.uv);
        auto brdfcos = eval_bsdf(f, n, o, i) * fabs(dot(n, i));
        if (pdf != 0) l += le * brdfcos / pdf;
    }

    // direct environments
    for (auto env : scn->environments) {
        auto i = sample_environment(env, rand1f(rng), rand2f(rng));
        auto isec = intersect_ray(scn, make_ray(p, i));
        if (isec.ist) continue;
        auto pdf = sample_environment_pdf(env, i);
        auto le = eval_environment(env, i);
        auto brdfcos = eval_bsdf(f, n, o, i) * fabs(dot(n, i));
        if (pdf != 0) l += le * brdfcos / pdf;
    }

    // exit if needed
    if (nbounces <= 0) return l;

    // reflection
    if (f.ks != zero3f && !f.rs) {
        auto i = reflect(o, n);
        l += f.ks * trace_direct(scn, make_ray(p, i), rng, nbounces - 1);
    }

    // opacity
    if (f.kt != zero3f) {
        l += f.kt * trace_direct(scn, make_ray(p, -o), rng, nbounces - 1);
    }

    // opacity
    auto op = eval_opacity(isec.ist, isec.ei, isec.uv);
    if (op != 1) {
        l = op * l +
            (1 - op) * trace_direct(scn, make_ray(p, -o), rng, nbounces - 1);
    }

    // done
    return l;
}

// Environment illumination only with no shadows.
vec3f trace_environment(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit) {
    if (scn->environments.empty()) return zero3f;

    // intersect scene
    auto isec = intersect_ray(scn, ray);
    auto l = zero3f;

    // handle environment
    if (!isec.ist) {
        for (auto env : scn->environments) l += eval_environment(env, ray.d);
        return l;
    }
    if (hit) *hit = true;

    // point
    auto o = -ray.d;
    auto p = eval_pos(isec.ist, isec.ei, isec.uv);
    auto n = eval_shading_norm(isec.ist, isec.ei, isec.uv, o);
    auto f = eval_bsdf(isec.ist, isec.ei, isec.uv);

    // emission
    l += eval_emission(isec.ist, isec.ei, isec.uv);

    // pick indirect direction
    auto i = zero3f, brdfcos = zero3f;
    auto pdf = 0.0f;
    if (!is_delta_bsdf(f)) {
        i = sample_brdf(f, n, o, rand1f(rng), rand2f(rng));
        brdfcos = eval_bsdf(f, n, o, i) * fabs(dot(n, i));
        pdf = sample_brdf_pdf(f, n, o, i);
    } else {
        i = sample_delta_brdf(f, n, o, rand1f(rng), rand2f(rng));
        brdfcos = eval_delta_brdf(f, n, o, i) * fabs(dot(n, i));
        pdf = sample_delta_brdf_pdf(f, n, o, i);
    }

    // accumulate environment illumination
    if (pdf != 0 && brdfcos != zero3f) {
        for (auto env : scn->environments)
            l += brdfcos * eval_environment(env, i) / pdf;
    }

    // exit if needed
    if (nbounces <= 0) return l;

    // opacity
    auto op = eval_opacity(isec.ist, isec.ei, isec.uv);
    if (op != 1) {
        l = op * l + (1 - op) * trace_direct(scn, make_ray(p, -o), rng,
                                    nbounces - 1, hit);
    }

    // done
    return l;
}

// Eyelight for quick previewing.
vec3f trace_eyelight(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit) {
    // intersect scene
    auto isec = intersect_ray(scn, ray);
    auto l = zero3f;

    // handle environment
    if (!isec.ist) {
        for (auto env : scn->environments) l += eval_environment(env, ray.d);
        return l;
    }
    if (hit) *hit = true;

    // point
    auto o = -ray.d;
    auto p = eval_pos(isec.ist, isec.ei, isec.uv);
    auto n = eval_shading_norm(isec.ist, isec.ei, isec.uv, o);
    auto f = eval_bsdf(isec.ist, isec.ei, isec.uv);

    // emission
    l += eval_emission(isec.ist, isec.ei, isec.uv);

    // bsdf*light
    l += eval_bsdf(f, n, o, o) * fabs(dot(n, o)) * pi;

    // opacity
    if (nbounces <= 0) return l;
    if (f.kt != zero3f) {
        l += f.kt * trace_eyelight(scn, make_ray(p, -o), rng, nbounces - 1);
    }
    auto op = eval_opacity(isec.ist, isec.ei, isec.uv);
    if (op != 1) {
        l = op * l +
            (1 - op) * trace_eyelight(scn, make_ray(p, -o), rng, nbounces - 1);
    }

    // done
    return l;
}

// Debug previewing.
vec3f trace_debug_normal(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit) {
    // intersect scene
    auto isec = intersect_ray(scn, ray);
    if (!isec.ist) return zero3f;
    if (hit) *hit = true;

    // point
    auto o = -ray.d;
    auto n = eval_shading_norm(isec.ist, isec.ei, isec.uv, o);

    // shade
    return n * 0.5f + vec3f{0.5f, 0.5f, 0.5f};
}

// Debug frontfacing.
vec3f trace_debug_frontfacing(const std::shared_ptr<scene>& scn,
    const ray3f& ray, rng_state& rng, int nbounces, bool* hit) {
    // intersect scene
    auto isec = intersect_ray(scn, ray);
    if (!isec.ist) return zero3f;
    if (hit) *hit = true;

    // point
    auto o = -ray.d;
    auto n = eval_shading_norm(isec.ist, isec.ei, isec.uv, o);

    // shade
    return dot(n, o) > 0 ? vec3f{0, 1, 0} : vec3f{1, 0, 0};
}

// Debug previewing.
vec3f trace_debug_albedo(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit) {
    // intersect scene
    auto isec = intersect_ray(scn, ray);
    if (!isec.ist) return zero3f;
    if (hit) *hit = true;

    // point
    auto f = eval_bsdf(isec.ist, isec.ei, isec.uv);

    // shade
    return f.kd + f.ks + f.kt;
}

// Debug previewing.
vec3f trace_debug_diffuse(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit) {
    // intersect scene
    auto isec = intersect_ray(scn, ray);
    if (!isec.ist) return zero3f;
    if (hit) *hit = true;

    // point
    auto f = eval_bsdf(isec.ist, isec.ei, isec.uv);

    // shade
    return f.kd;
}

// Debug previewing.
vec3f trace_debug_specular(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit) {
    // intersect scene
    auto isec = intersect_ray(scn, ray);
    if (!isec.ist) return zero3f;
    if (hit) *hit = true;

    // point
    auto f = eval_bsdf(isec.ist, isec.ei, isec.uv);

    // shade
    return f.ks;
}

// Debug previewing.
vec3f trace_debug_roughness(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit) {
    // intersect scene
    auto isec = intersect_ray(scn, ray);
    if (!isec.ist) return zero3f;
    if (hit) *hit = true;

    // point
    auto f = eval_bsdf(isec.ist, isec.ei, isec.uv);

    // shade
    return {f.rs, f.rs, f.rs};
}

// Debug previewing.
vec3f trace_debug_texcoord(const std::shared_ptr<scene>& scn, const ray3f& ray,
    rng_state& rng, int nbounces, bool* hit) {
    // intersect scene
    auto isec = intersect_ray(scn, ray);
    if (!isec.ist) return zero3f;
    if (hit) *hit = true;

    // point
    auto texcoord = eval_texcoord(isec.ist, isec.ei, isec.uv);

    // shade
    return {texcoord.x, texcoord.y, 0};
}

// Trace a single sample
vec4f trace_sample(const std::shared_ptr<scene>& scn,
    const std::shared_ptr<camera>& cam, int i, int j, int width, int height,
    rng_state& rng, trace_func tracer, int nbounces, float pixel_clamp = 100) {
    auto ray =
        eval_camera_ray(cam, {i, j}, {width, height}, rand2f(rng), rand2f(rng));
    auto hit = false;
    auto l = tracer(scn, ray, rng, nbounces, &hit);
    if (!isfinite(l.x) || !isfinite(l.y) || !isfinite(l.z)) {
        std::cout << "NaN detected\n";
        l = zero3f;
    }
    if (max(l) > pixel_clamp) l = l * (pixel_clamp / max(l));
    return {l.x, l.y, l.z, (hit || !scn->environments.empty()) ? 1.0f : 0.0f};
}

// Trace the next nsamples.
void trace_samples(const std::shared_ptr<scene>& scn,
    const std::shared_ptr<camera>& cam, image4f& img, image<rng_state>& rngs,
    int sample, int nsamples, trace_func tracer, int nbounces,
    float pixel_clamp) {
    for (auto j = 0; j < img.height(); j++) {
        for (auto i = 0; i < img.width(); i++) {
            img[{i, j}] *= sample;
            for (auto s = 0; s < nsamples; s++)
                img[{i, j}] += trace_sample(scn, cam, i, j, img.width(),
                    img.height(), rngs[{i, j}], tracer, nbounces, pixel_clamp);
            img[{i, j}] /= sample + nsamples;
        }
    }
}

// Trace the next nsamples.
void trace_samples_mt(const std::shared_ptr<scene>& scn,
    const std::shared_ptr<camera>& cam, image4f& img, image<rng_state>& rngs,
    int sample, int nsamples, trace_func tracer, int nbounces,
    float pixel_clamp) {
    auto nthreads = std::thread::hardware_concurrency();
    auto threads = std::vector<std::thread>();
    for (auto tid = 0; tid < std::thread::hardware_concurrency(); tid++) {
        threads.push_back(std::thread([=, &img, &rngs]() {
            for (auto j = tid; j < img.height(); j += nthreads) {
                for (auto i = 0; i < img.width(); i++) {
                    img[{i, j}] *= sample;
                    for (auto s = 0; s < nsamples; s++)
                        img[{i, j}] += trace_sample(scn, cam, i, j, img.width(),
                            img.height(), rngs[{i, j}], tracer, nbounces,
                            pixel_clamp);
                    img[{i, j}] /= sample + nsamples;
                }
            }
        }));
    }
    for (auto& t : threads) t.join();
    threads.clear();
}

// Starts an anyncrhounous renderer.
void trace_async_start(const std::shared_ptr<scene>& scn,
    const std::shared_ptr<camera>& cam, image4f& img, image<rng_state>& rngs,
    int nsamples, trace_func tracer, int nbounces,
    std::vector<std::thread>& threads, bool& stop_flag, int& cur_sample,
    float pixel_clamp) {
    auto nthreads = std::thread::hardware_concurrency();
    for (auto tid = 0; tid < nthreads; tid++) {
        threads.push_back(std::thread([=, &img, &rngs, &stop_flag,
                                          &cur_sample]() {
            for (auto sample = 0; sample < nsamples; sample++) {
                if (!tid) cur_sample = sample;
                for (auto j = tid; j < img.height(); j += nthreads) {
                    for (auto i = 0; i < img.width(); i++) {
                        if (stop_flag) return;
                        auto l = trace_sample(scn, cam, i, j, img.width(),
                            img.height(), rngs[{i, j}], tracer, nbounces);
                        img[{i, j}] = (img[{i, j}] * sample + l) / (sample + 1);
                    }
                }
                if (!tid) cur_sample = nsamples;
            }
        }));
    }
}

// Stop the asynchronous renderer.
void trace_async_stop(std::vector<std::thread>& threads, bool& stop_flag) {
    stop_flag = true;
    for (auto& t : threads) t.join();
    stop_flag = false;
    threads.clear();
}

}  // namespace ygl
// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PATH TRACING SUPPORT FUNCTION
// -----------------------------------------------------------------------------
namespace ygl {

// Phong exponent to roughness.
float specular_exponent_to_roughness(float n) { return sqrtf(2 / (n + 2)); }

// Specular to fresnel eta.
void specular_fresnel_from_ks(const vec3f& ks, vec3f& es, vec3f& esk) {
    es = {(1 + sqrt(ks.x)) / (1 - sqrt(ks.x)),
        (1 + sqrt(ks.y)) / (1 - sqrt(ks.y)),
        (1 + sqrt(ks.z)) / (1 - sqrt(ks.z))};
    esk = {0, 0, 0};
}

// Specular to  eta.
float specular_to_eta(const vec3f& ks) {
    auto f0 = (ks.x + ks.y + ks.z) / 3;
    return (1 + sqrt(f0)) / (1 - sqrt(f0));
}

// Compute the fresnel term for dielectrics. Implementation from
// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
vec3f fresnel_dielectric(float cosw, const vec3f& eta_) {
    auto eta = eta_;
    if (cosw < 0) {
        eta = 1.0f / eta;
        cosw = -cosw;
    }

    auto sin2 = 1 - cosw * cosw;
    auto eta2 = eta * eta;

    auto cos2t = vec3f{1, 1, 1} - sin2 / eta2;
    if (cos2t.x < 0 || cos2t.y < 0 || cos2t.z < 0)
        return vec3f{1, 1, 1};  // tir

    auto t0 = vec3f{sqrt(cos2t.x), sqrt(cos2t.y), sqrt(cos2t.z)};
    auto t1 = eta * t0;
    auto t2 = eta * cosw;

    auto rs = (vec3f{cosw, cosw, cosw} - t1) / (vec3f{cosw, cosw, cosw} + t1);
    auto rp = (t0 - t2) / (t0 + t2);

    return (rs * rs + rp * rp) / 2.0f;
}

// Compute the fresnel term for metals. Implementation from
// https://seblagarde.wordpress.com/2013/04/29/memo-on-fresnel-equations/
vec3f fresnel_metal(float cosw, const vec3f& eta, const vec3f& etak) {
    if (etak == zero3f) return fresnel_dielectric(cosw, eta);

    cosw = clamp(cosw, (float)-1, (float)1);
    auto cos2 = cosw * cosw;
    auto sin2 = clamp(1 - cos2, (float)0, (float)1);
    auto eta2 = eta * eta;
    auto etak2 = etak * etak;

    auto t0 = eta2 - etak2 - vec3f{sin2, sin2, sin2};
    auto a2plusb2_2 = t0 * t0 + 4.0f * eta2 * etak2;
    auto a2plusb2 =
        vec3f{sqrt(a2plusb2_2.x), sqrt(a2plusb2_2.y), sqrt(a2plusb2_2.z)};
    auto t1 = a2plusb2 + vec3f{cos2, cos2, cos2};
    auto a_2 = (a2plusb2 + t0) / 2.0f;
    auto a = vec3f{sqrt(a_2.x), sqrt(a_2.y), sqrt(a_2.z)};
    auto t2 = 2.0f * a * cosw;
    auto rs = (t1 - t2) / (t1 + t2);

    auto t3 = vec3f{cos2, cos2, cos2} * a2plusb2 +
              vec3f{sin2, sin2, sin2} * vec3f{sin2, sin2, sin2};
    auto t4 = t2 * sin2;
    auto rp = rs * (t3 - t4) / (t3 + t4);

    return (rp + rs) / 2.0f;
}

// Schlick approximation of the Fresnel term
vec3f fresnel_schlick(const vec3f& ks, float cosw) {
    if (ks == zero3f) return zero3f;
    return ks +
           (vec3f{1, 1, 1} - ks) * pow(clamp(1.0f - cosw, 0.0f, 1.0f), 5.0f);
}
vec3f fresnel_schlick(const vec3f& ks, float cosw, float rs) {
    if (ks == zero3f) return zero3f;
    auto fks = fresnel_schlick(ks, cosw);
    return ks + (fks - ks) * (1 - sqrt(clamp(rs, 0.0f, 1.0f)));
}

// Evaluates the GGX pdf
float sample_ggx_pdf(float rs, float ndh) {
    auto alpha2 = rs * rs;
    auto di = (ndh * ndh) * (alpha2 - 1) + 1;
    auto d = alpha2 / (pi * di * di);
    return d * ndh;
}

// Sample the GGX distribution
vec3f sample_ggx(float rs, const vec2f& rn) {
    auto tan2 = rs * rs * rn.y / (1 - rn.y);
    auto rz = sqrt(1 / (tan2 + 1));
    auto rr = sqrt(1 - rz * rz);
    auto rphi = 2 * pi * rn.x;
    // set to wh
    auto wh_local = vec3f{rr * cos(rphi), rr * sin(rphi), rz};
    return wh_local;
}

}  // namespace ygl
