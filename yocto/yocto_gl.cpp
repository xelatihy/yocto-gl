//
// Implementation for Yocto/GL. See yocto_gl.h for documentation.
//

//
// LICENSE:
//
// Copyright (c) 2016 -- 2017 Fabio Pellacini
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
//
// LICENSE OF INCLUDED CODE FOR PERLIN NOISE (stb_perlin.h)
//
// Copyright (c) 2017 Sean Barrett
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions: The above copyright
// notice and this permission notice shall be included in all copies or
// substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS",
// WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
// TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
// FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
// TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR
// THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//
//
// LICENSE OF INCLUDED CODE FOR BASE64 (base64.h, base64.cpp)
//
// Copyright (C) 2004-2008 René Nyffenegger
//
// This source code is provided 'as-is', without any express or implied
// warranty. In no event will the author be held liable for any damages
// arising from the use of this software.
//
// Permission is granted to anyone to use this software for any purpose,
// including commercial applications, and to alter it and redistribute it
// freely, subject to the following restrictions:
//
// 1. The origin of this source code must not be misrepresented; you must not
// claim that you wrote the original source code. If you use this source code
// in a product, an acknowledgment in the product documentation would be
// appreciated but is not required.
//
// 2. Altered source versions must be plainly marked as such, and must not be
// misrepresented as being the original source code.
//
// 3. This notice may not be removed or altered from any source distribution.
//
// René Nyffenegger rene.nyffenegger@adp-gmbh.ch
//
//

#include "yocto_gl.h"

#include <fstream>

#if YGL_IMAGEIO
#include "ext/stb_image.h"
#include "ext/stb_image_resize.h"
#include "ext/stb_image_write.h"
#include "ext/tinyexr.h"
#endif

#if YGL_GLTF
#include "ext/json.hpp"
#endif

#if YGL_OPENGL
#ifdef __APPLE__
#include <OpenGL/gl3.h>
#define GLFW_INCLUDE_GLCOREARB
#else
#include <GL/glew.h>
#endif
#include <GLFW/glfw3.h>
#include "ext/imgui/imgui.h"
#include "ext/imgui/imgui_impl_glfw_gl3.h"
unsigned int imgui_extrafont_compressed_size();
const unsigned int* imgui_extrafont_compressed_data();
#endif

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR UI UTILITIES
// -----------------------------------------------------------------------------
namespace ygl {

// Turntable for UI navigation.
void camera_turntable(vec3f& from, vec3f& to, vec3f& up, const vec2f& rotate,
    float dolly, const vec2f& pan) {
    // rotate if necessary
    if (rotate.x || rotate.y) {
        auto z = normalize(to - from);
        auto lz = length(to - from);
        auto phi = atan2(z.z, z.x) + rotate.x;
        auto theta = acos(z.y) + rotate.y;
        theta = clamp(theta, 0.001f, pi - 0.001f);
        auto nz = vec3f{sin(theta) * cos(phi) * lz, cos(theta) * lz,
            sin(theta) * sin(phi) * lz};
        from = to - nz;
    }

    // dolly if necessary
    if (dolly) {
        auto z = normalize(to - from);
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
void camera_turntable(frame3f& frame, float& focus, const vec2f& rotate,
    float dolly, const vec2f& pan) {
    // rotate if necessary
    if (rotate != zero2f) {
        auto phi = atan2(frame.z.z, frame.z.x) + rotate.x;
        auto theta = acos(frame.z.y) + rotate.y;
        theta = clamp(theta, 0.001f, pi - 0.001f);
        auto new_z =
            vec3f{sin(theta) * cos(phi), cos(theta), sin(theta) * sin(phi)};
        auto new_center = frame.o - frame.z * focus;
        auto new_o = new_center + new_z * focus;
        frame = lookat_frame(new_o, new_center, {0, 1, 0});
        focus = length(new_o - new_center);
    }

    // pan if necessary
    if (dolly) {
        auto c = frame.o - frame.z * focus;
        focus = max(focus * (1 + dolly), 0.001f);
        frame.o = c + frame.z * focus;
    }

    // pan if necessary
    if (pan.x || pan.y) { frame.o += frame.x * pan.x + frame.y * pan.y; }
}

// FPS camera for UI navigation for a frame parametrization.
void camera_fps(frame3f& frame, const vec3f& transl, const vec2f& rotate) {
    // https://gamedev.stackexchange.com/questions/30644/how-to-keep-my-quaternion-using-fps-camera-from-tilting-and-messing-up
    auto y = vec3f{0, 1, 0};
    auto z = orthonormalize(frame.z, y);
    auto x = cross(y, z);

    auto rot = rotation_frame(vec3f{1, 0, 0}, rotate.y) *
               frame3f{frame.x, frame.y, frame.z, zero3f} *
               rotation_frame(vec3f{0, 1, 0}, rotate.x);
    auto pos = frame.o + transl.x * x + transl.y * y + transl.z * z;

    frame = {rot.x, rot.y, rot.z, pos};
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PERLIN NOISE
// -----------------------------------------------------------------------------
namespace ygl {

// clang-format off
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

static float stb__perlin_lerp(float a, float b, float t)
{
   return a + (b-a) * t;
}

static int stb__perlin_fastfloor(float a)
{
	int ai = (int) a;
	return (a < ai) ? ai-1 : ai;
}

// different grad function from Perlin's, but easy to modify to match reference
static float stb__perlin_grad(int hash, float x, float y, float z)
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
float perlin_noise(const vec3f& p, const vec3i& wrap) {
    return stb_perlin_noise3(p.x, p.y, p.z, wrap.x, wrap.y, wrap.z);
}

// adapeted  stb_perlin.h
float perlin_ridge_noise(const vec3f& p, float lacunarity, float gain,
    float offset, int octaves, const vec3i& wrap) {
    return stb_perlin_ridge_noise3(p.x, p.y, p.z, lacunarity, gain, offset,
        octaves, wrap.x, wrap.y, wrap.z);
}

// adapeted  stb_perlin.h
float perlin_fbm_noise(const vec3f& p, float lacunarity, float gain,
    int octaves, const vec3i& wrap) {
    return stb_perlin_fbm_noise3(
        p.x, p.y, p.z, lacunarity, gain, octaves, wrap.x, wrap.y, wrap.z);
}

// adapeted  stb_perlin.h
float perlin_turbulence_noise(const vec3f& p, float lacunarity, float gain,
    int octaves, const vec3i& wrap) {
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
    const std::vector<vec3f>& pos, std::vector<vec3f>& norm, bool weighted) {
    norm.resize(pos.size());
    for (auto& n : norm) n = zero3f;
    for (auto& l : lines) {
        auto n = pos[l.y] - pos[l.x];
        if (!weighted) n = normalize(n);
        norm[l.x] += n;
        norm[l.y] += n;
    }
    for (auto& n : norm) n = normalize(n);
}

// Compute per-vertex normals for triangles.
void compute_normals(const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, std::vector<vec3f>& norm, bool weighted) {
    norm.resize(pos.size());
    for (auto& n : norm) n = zero3f;
    for (auto& t : triangles) {
        auto n = cross(pos[t.y] - pos[t.x], pos[t.z] - pos[t.x]);
        if (!weighted) n = normalize(n);
        for (auto vid : {t.x, t.y, t.z}) norm[vid] += n;
    }
    for (auto& n : norm) n = normalize(n);
}

// Compute per-vertex normals for quads.
void compute_normals(const std::vector<vec4i>& quads,
    const std::vector<vec3f>& pos, std::vector<vec3f>& norm, bool weighted) {
    norm.resize(pos.size());
    for (auto& n : norm) n = zero3f;
    for (auto& q : quads) {
        auto n = cross(pos[q.y] - pos[q.x], pos[q.w] - pos[q.x]) +
                 cross(pos[q.w] - pos[q.z], pos[q.x] - pos[q.z]);
        if (!weighted) n = normalize(n);
        for (auto vid : {q.x, q.y, q.z, q.w}) norm[vid] += n;
    }
    for (auto& n : norm) n = normalize(n);
}

// Compute per-vertex tangent frame for triangle meshes.
// Tangent space is defined by a four component vector.
// The first three components are the tangent with respect to the U texcoord.
// The fourth component is the sign of the tangent wrt the V texcoord.
// Tangent frame is useful in normal mapping.
void compute_tangent_frames(const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, const std::vector<vec3f>& norm,
    const std::vector<vec2f>& texcoord, std::vector<vec4f>& tangsp,
    bool weighted) {
    auto tangu = std::vector<vec3f>(pos.size(), zero3f);
    auto tangv = std::vector<vec3f>(pos.size(), zero3f);
    for (auto& t : triangles) {
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
std::vector<vec2i> get_edges(const std::vector<vec2i>& lines,
    const std::vector<vec3i>& triangles, const std::vector<vec4i>& quads) {
    auto edges = std::vector<vec2i>();
    auto eset = std::unordered_set<vec2i>();
    for (auto e : lines) {
        e = {min(e.x, e.y), max(e.x, e.y)};
        if (!eset.insert(e).second) continue;
        eset.insert({e.y, e.x});
        edges.push_back(e);
    }
    for (auto& t : triangles) {
        for (auto e : {vec2i{t.x, t.y}, vec2i{t.y, t.z}, vec2i{t.z, t.x}}) {
            e = {min(e.x, e.y), max(e.x, e.y)};
            if (!eset.insert(e).second) continue;
            eset.insert({e.y, e.x});
            edges.push_back(e);
        }
    }
    for (auto& q : quads) {
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
    for (auto& q : quads) {
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
    for (auto& q : quads) {
        triangles.push_back({q.x, q.y, q.w});
        if (q.z != q.w) triangles.push_back({q.z, q.w, q.y});
    }
    return triangles;
#if 0
        triangles.resize(usteps * vsteps * 2);
        for (auto j = 0; j < vsteps; j++) {
            for (auto i = 0; i < usteps; i++) {
                auto& f1 = triangles[(j * usteps + i) * 2 + 0];
                auto& f2 = triangles[(j * usteps + i) * 2 + 1];
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
    for (auto& b : beziers) {
        lines.push_back({b.x, b.y});
        lines.push_back({b.y, b.z});
        lines.push_back({b.z, b.w});
    }
    return lines;
}

// Convert face varying data to single primitives. Returns the quads indices
// and filled vectors for pos, norm and texcoord.
std::tuple<std::vector<vec4i>, std::vector<vec3f>, std::vector<vec3f>,
    std::vector<vec2f>>
convert_face_varying(const std::vector<vec4i>& quads_pos,
    const std::vector<vec4i>& quads_norm,
    const std::vector<vec4i>& quads_texcoord, const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec2f>& texcoord) {
    // make faces unique
    std::unordered_map<vec3i, int> vert_map;
    auto quads = std::vector<vec4i>(quads_pos.size());
    for (auto fid = 0; fid < quads_pos.size(); fid++) {
        for (auto c = 0; c < 4; c++) {
            auto v = vec3i{
                (&quads_pos[fid].x)[c],
                (!quads_norm.empty()) ? (&quads_norm[fid].x)[c] : -1,
                (!quads_texcoord.empty()) ? (&quads_texcoord[fid].x)[c] : -1,
            };
            if (vert_map.find(v) == vert_map.end()) {
                auto s = (int)vert_map.size();
                vert_map[v] = s;
            }
            (&quads[fid].x)[c] = vert_map.at(v);
        }
    }

    // fill vert data
    auto qpos = std::vector<vec3f>();
    if (!pos.empty()) {
        qpos.resize(vert_map.size());
        for (auto& kv : vert_map) { qpos[kv.second] = pos[kv.first.x]; }
    }
    auto qnorm = std::vector<vec3f>();
    if (!norm.empty()) {
        qnorm.resize(vert_map.size());
        for (auto& kv : vert_map) { qnorm[kv.second] = norm[kv.first.y]; }
    }
    auto qtexcoord = std::vector<vec2f>();
    if (!texcoord.empty()) {
        qtexcoord.resize(vert_map.size());
        for (auto& kv : vert_map) {
            qtexcoord[kv.second] = texcoord[kv.first.z];
        }
    }

    // done
    return {quads, qpos, qnorm, qtexcoord};
}

// Subdivide lines.
template <typename T>
void subdivide_lines(
    std::vector<vec2i>& lines, std::vector<T>& vert, bool update_lines) {
    if (lines.empty() || vert.empty()) return;

    auto tvert = vert;
    auto tlines = std::vector<vec2i>();
    for (auto& l : lines) {
        tlines.push_back({l.x, (int)tvert.size()});
        tlines.push_back({(int)tvert.size(), l.y});
        tvert.push_back((vert[l.x] + vert[l.y]) / 2);
    }

    std::swap(vert, tvert);
    if (update_lines) std::swap(lines, tlines);
}

// Subdivide triangle.
template <typename T>
void subdivide_triangles(std::vector<vec3i>& triangles, std::vector<T>& vert,
    bool update_triangles) {
    if (triangles.empty() || vert.empty()) return;

    auto tvert = vert;
    auto ttriangles = std::vector<vec3i>();
    auto emap = std::unordered_map<vec2i, int>();
    for (auto& t : triangles) {
        for (auto e : {vec2i{t.x, t.y}, vec2i{t.y, t.z}, vec2i{t.z, t.x}}) {
            if (emap.find(e) != emap.end()) continue;
            emap[{e.x, e.y}] = (int)tvert.size();
            emap[{e.y, e.x}] = (int)tvert.size();
            tvert.push_back((vert[e.x] + vert[e.y]) / 2);
        }
        ttriangles.push_back({t.x, emap.at({t.x, t.y}), emap.at({t.z, t.x})});
        ttriangles.push_back({t.y, emap.at({t.y, t.z}), emap.at({t.x, t.y})});
        ttriangles.push_back({t.z, emap.at({t.z, t.x}), emap.at({t.y, t.z})});
        ttriangles.push_back(
            {emap.at({t.x, t.y}), emap.at({t.y, t.z}), emap.at({t.z, t.x})});
    }

    std::swap(vert, tvert);
    if (update_triangles) std::swap(triangles, ttriangles);
}

// Subdivide quads.
template <typename T>
void subdivide_quads(
    std::vector<vec4i>& quads, std::vector<T>& vert, bool update_quads) {
    if (quads.empty() || vert.empty()) return;

    auto tvert = vert;
    auto tquads = std::vector<vec4i>();
    auto emap = std::unordered_map<vec2i, int>();
    for (auto& q : quads) {
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
    if (update_quads) std::swap(quads, tquads);
}

// Subdivide beziers.
template <typename T>
void subdivide_beziers(
    std::vector<vec4i>& beziers, std::vector<T>& vert, bool update_beziers) {
    if (beziers.empty() || vert.empty()) return;

    auto vmap = std::unordered_map<int, int>();
    auto tvert = std::vector<T>();
    auto tbeziers = std::vector<vec4i>();
    for (auto& b : beziers) {
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
        tvert.push_back(vert[b.x] / 8 + 3 * vert[b.y] / 8 + 3 * vert[b.z] / 8 +
                        vert[b.w] / 8);
        tvert.push_back(vert[b.y] / 4 + vert[b.z] / 2 + vert[b.w] / 4);
        tvert.push_back(vert[b.z] / 2 + vert[b.w] / 2);
    }

    std::swap(vert, tvert);
    if (update_beziers) std::swap(beziers, tbeziers);
}

// Subdivide catmullclark.
template <typename T>
void subdivide_catmullclark(
    std::vector<vec4i>& quads, std::vector<T>& vert, bool update_quads) {
    if (quads.empty() || vert.empty()) return;

    auto tvert = vert;
    auto tquads = std::vector<vec4i>();
    auto emap = std::unordered_map<vec2i, int>();
    for (auto& q : quads) {
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

    auto tboundary = std::vector<vec2i>();
    for (auto e_kv : emap) {
        auto e = e_kv.first;
        auto v = e_kv.second;
        if (emap.find({e.y, e.x}) != emap.end()) continue;
        tboundary.push_back({e.x, v});
        tboundary.push_back({v, e.y});
    }

    // setup
    auto tcrease_edges = tboundary;
    auto tcrease_verts = std::vector<int>();

    // define vertex valence ---------------------------
    auto tvert_val = std::vector<int>(tvert.size(), 2);
    for (auto e : tcrease_edges) {
        tvert_val[e.x] = 1;
        tvert_val[e.y] = 1;
    }
    for (auto vid : tcrease_verts) tvert_val[vid] = 0;

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
    for (auto& q : tquads) {
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
    if (update_quads) std::swap(quads, tquads);
}

// Subdivide lines by splitting each line in half.
void subdivide_lines(std::vector<vec2i>& lines, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<vec4f>& color, std::vector<float>& radius) {
    if (lines.empty()) return;
    subdivide_lines(lines, norm, false);
    for (auto& n : norm) n = normalize(n);
    subdivide_lines(lines, texcoord, false);
    subdivide_lines(lines, color, false);
    subdivide_lines(lines, radius, false);
    subdivide_lines(lines, pos);
}

// Subdivide triangles.
void subdivide_triangles(std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<vec4f>& color, std::vector<float>& radius) {
    if (triangles.empty()) return;
    subdivide_triangles(triangles, norm, false);
    for (auto& n : norm) n = normalize(n);
    subdivide_triangles(triangles, texcoord, false);
    subdivide_triangles(triangles, color, false);
    subdivide_triangles(triangles, radius, false);
    subdivide_triangles(triangles, pos);
}

// Subdivide quads.
void subdivide_quads(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<vec4f>& color, std::vector<float>& radius) {
    if (quads.empty()) return;
    subdivide_quads(quads, norm, false);
    for (auto& n : norm) n = normalize(n);
    subdivide_quads(quads, texcoord, false);
    subdivide_quads(quads, color, false);
    subdivide_quads(quads, radius, false);
    subdivide_quads(quads, pos);
}

// Subdivide beziers.
void subdivide_beziers(std::vector<vec4i>& beziers, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<vec4f>& color, std::vector<float>& radius) {
    if (beziers.empty()) return;
    subdivide_beziers(beziers, norm, false);
    for (auto& n : norm) n = normalize(n);
    subdivide_beziers(beziers, texcoord, false);
    subdivide_beziers(beziers, color, false);
    subdivide_beziers(beziers, radius, false);
    subdivide_beziers(beziers, pos);
}

// Subdivide quads.
void subdivide_catmullclark(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<vec4f>& color, std::vector<float>& radius) {
    if (quads.empty()) return;
    subdivide_catmullclark(quads, norm, false);
    for (auto& n : norm) n = normalize(n);
    subdivide_catmullclark(quads, texcoord, false);
    subdivide_catmullclark(quads, color, false);
    subdivide_catmullclark(quads, radius, false);
    subdivide_catmullclark(quads, pos);
}

// Merge lines between shapes.
void merge_lines(std::vector<vec2i>& lines, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    const std::vector<vec2i>& lines1, const std::vector<vec3f>& pos1,
    const std::vector<vec3f>& norm1, const std::vector<vec2f>& texcoord1) {
    lines.reserve(lines.size() + lines1.size());
    auto nverts = (int)pos.size();
    for (auto& l : lines1) lines.push_back({l.x + nverts, l.y + nverts});
    for (auto& v : pos1) pos.push_back(v);
    for (auto& v : norm1) norm.push_back(v);
    ;
    for (auto& v : texcoord1) texcoord.push_back(v);
}

// Merge triangles between shapes.
void merge_triangles(std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    const std::vector<vec3i>& triangles1, const std::vector<vec3f>& pos1,
    const std::vector<vec3f>& norm1, const std::vector<vec2f>& texcoord1) {
    triangles.reserve(triangles.size() + triangles1.size());
    auto nverts = (int)pos.size();
    for (auto& t : triangles1)
        triangles.push_back({t.x + nverts, t.y + nverts, t.z + nverts});
    for (auto& v : pos1) pos.push_back(v);
    for (auto& v : norm1) norm.push_back(v);
    for (auto& v : texcoord1) texcoord.push_back(v);
}

// Merge quads between shapes.
void merge_quads(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    const std::vector<vec4i>& quads1, const std::vector<vec3f>& pos1,
    const std::vector<vec3f>& norm1, const std::vector<vec2f>& texcoord1) {
    quads.reserve(quads.size() + quads1.size());
    auto nverts = (int)pos.size();
    for (auto& q : quads1)
        quads.push_back(
            {q.x + nverts, q.y + nverts, q.z + nverts, q.w + nverts});
    for (auto& v : pos1) pos.push_back(v);
    for (auto& v : norm1) norm.push_back(v);
    for (auto& v : texcoord1) texcoord.push_back(v);
}

// Samples a set of points over a triangle mesh uniformly. The rng function
// takes the point index and returns vec3f numbers uniform directibuted in
// [0,1]^3. unorm and texcoord are optional.
std::tuple<std::vector<vec3f>, std::vector<vec3f>, std::vector<vec2f>>
sample_triangles_points(const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, const std::vector<vec3f>& norm,
    const std::vector<vec2f>& texcoord, int npoints, uint64_t seed) {
    auto sampled_pos = std::vector<vec3f>(npoints);
    auto sampled_norm = std::vector<vec3f>(npoints);
    auto sampled_texcoord = std::vector<vec2f>(npoints);
    auto cdf = sample_triangles_cdf(triangles, pos);
    auto rng = make_rng(seed);
    for (auto i = 0; i < npoints; i++) {
        auto eid = 0;
        auto euv = zero2f;
        std::tie(eid, euv) = sample_triangles(
            cdf, next_rand1f(rng), {next_rand1f(rng), next_rand1f(rng)});
        auto t = triangles[eid];
        sampled_pos[i] =
            interpolate_triangle(pos[t.x], pos[t.y], pos[t.z], euv);
        if (!sampled_norm.empty()) {
            sampled_norm[i] = normalize(
                interpolate_triangle(norm[t.x], norm[t.y], norm[t.z], euv));
        } else {
            sampled_norm[i] = triangle_normal(pos[t.x], pos[t.y], pos[t.z]);
        }
        if (!sampled_texcoord.empty()) {
            sampled_texcoord[i] = interpolate_triangle(
                texcoord[t.x], texcoord[t.y], texcoord[t.z], euv);
        } else {
            sampled_texcoord[i] = zero2f;
        }
    }

    return {sampled_pos, sampled_norm, sampled_texcoord};
}

// Make a quad.
void make_quad(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec2i& steps,
    const vec2f& size, const vec2f& uvsize) {
    auto nverts = (steps.x + 1) * (steps.y + 1);
    auto nfaces = steps.x * steps.y;
    auto vid = [steps](int i, int j) { return j * (steps.x + 1) + i; };
    auto fid = [steps](int i, int j) { return j * steps.x + i; };

    pos.resize(nverts);
    norm.resize(nverts);
    texcoord.resize(nverts);
    for (auto j = 0; j <= steps.y; j++) {
        for (auto i = 0; i <= steps.x; i++) {
            auto uv = vec2f{i / (float)steps.x, j / (float)steps.y};
            pos[vid(i, j)] = {
                (uv.x - 0.5f) * size.x, (uv.y - 0.5f) * size.y, 0};
            norm[vid(i, j)] = {0, 0, 1};
            texcoord[vid(i, j)] = uv * uvsize;
        }
    }

    quads.resize(nfaces);
    for (auto j = 0; j < steps.y; j++) {
        for (auto i = 0; i < steps.x; i++) {
            quads[fid(i, j)] = {
                vid(i, j), vid(i + 1, j), vid(i + 1, j + 1), vid(i, j + 1)};
        }
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
        auto uv = cartesian_to_elliptical(vec2f{pos[i].x, pos[i].y}) * size / 2;
        pos[i] = {uv.x, uv.y, 0};
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
    make_cylinder(quads, pos, norm, texcoord, steps, size, uvsize);
    auto c = size / 2 - vec2f{radius, radius};
    for (auto i = 0; i < pos.size(); i++) {
        auto pcyl = cartesian_to_cylindrical(pos[i]);
        auto pc = vec2f{pcyl.y, fabs(pcyl.z)};
        auto pp = pcyl.x;
        auto ps = (pcyl.z < 0) ? -1.0f : 1.0f;
        if (pc.x >= c.x && pc.y >= c.y) {
            auto pn = normalize(pc - c);
            pos[i] = cylindrical_to_cartesian(
                vec3f{pp, c.x + radius * pn.x, ps * (c.y + radius * pn.y)});
            norm[i] = cylindrical_to_cartesian(vec3f{pp, pn.x, ps * pn.y});
        } else {
            continue;
        }
    }
}

// Make a geodesic sphere.
void make_geodesic_sphere(
    std::vector<vec3i>& triangles, std::vector<vec3f>& pos, int tesselation) {
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
    for (auto l = 0; l < max(0, tesselation - 2); l++)
        subdivide_triangles(triangles, pos);
    for (auto& p : pos) p = normalize(p);
}

// Make a cube with unique vertices. This is watertight but has no
// texture coordinates or normals.
void make_cube(
    std::vector<vec4i>& quads, std::vector<vec3f>& pos, int tesselation) {
    pos = std::vector<vec3f>{{-1, -1, -1}, {-1, +1, -1}, {+1, +1, -1},
        {+1, -1, -1}, {-1, -1, +1}, {-1, +1, +1}, {+1, +1, +1}, {+1, -1, +1}};
    quads = std::vector<vec4i>{{0, 1, 2, 3}, {7, 6, 5, 4}, {4, 5, 1, 0},
        {6, 7, 3, 2}, {2, 1, 5, 6}, {0, 3, 7, 4}};

    for (auto l = 0; l < tesselation; l++) subdivide_quads(quads, pos);
}

// Make a facevarying cube with unique vertices but different texture
// coordinates.
void make_fvcube(std::vector<vec4i>& quads_pos, std::vector<vec3f>& pos,
    std::vector<vec4i>& quads_norm, std::vector<vec3f>& norm,
    std::vector<vec4i>& quads_texcoord, std::vector<vec2f>& texcoord,
    int tesselation) {
    pos = std::vector<vec3f>{{-1, -1, -1}, {-1, +1, -1}, {+1, +1, -1},
        {+1, -1, -1}, {-1, -1, +1}, {-1, +1, +1}, {+1, +1, +1}, {+1, -1, +1}};
    quads_pos = std::vector<vec4i>{{0, 1, 2, 3}, {7, 6, 5, 4}, {4, 5, 1, 0},
        {6, 7, 3, 2}, {2, 1, 5, 6}, {0, 3, 7, 4}};
    norm = std::vector<vec3f>{{0, 0, -1}, {0, 0, -1}, {0, 0, -1}, {0, 0, -1},
        {0, 0, +1}, {0, 0, +1}, {0, 0, +1}, {0, 0, +1}, {-1, 0, 0}, {-1, 0, 0},
        {-1, 0, 0}, {-1, 0, 0}, {+1, 0, 0}, {+1, 0, 0}, {+1, 0, 0}, {+1, 0, 0},
        {0, +1, 0}, {0, +1, 0}, {0, +1, 0}, {0, +1, 0}, {0, -1, 0}, {0, -1, 0},
        {0, -1, 0}, {0, -1, 0}};
    quads_norm = std::vector<vec4i>{{0, 1, 2, 3}, {4, 5, 6, 7}, {8, 9, 10, 11},
        {12, 13, 14, 15}, {16, 17, 18, 19}, {20, 21, 22, 23}};
    texcoord = std::vector<vec2f>{{0, 0}, {1, 0}, {1, 1}, {0, 1}, {0, 0},
        {1, 0}, {1, 1}, {0, 1}, {0, 0}, {1, 0}, {1, 1}, {0, 1}, {0, 0}, {1, 0},
        {1, 1}, {0, 1}, {0, 0}, {1, 0}, {1, 1}, {0, 1}, {0, 0}, {1, 0}, {1, 1},
        {0, 1}};
    quads_texcoord = std::vector<vec4i>{{0, 1, 2, 3}, {4, 5, 6, 7},
        {8, 9, 10, 11}, {12, 13, 14, 15}, {16, 17, 18, 19}, {20, 21, 22, 23}};

    for (auto l = 0; l < tesselation; l++) {
        subdivide_quads(quads_pos, pos);
        subdivide_quads(quads_norm, norm);
        subdivide_quads(quads_texcoord, texcoord);
    }
}

// Make a facevarying sphere with unique vertices but different texture
// coordinates.
void make_fvsphere(std::vector<vec4i>& quads_pos, std::vector<vec3f>& pos,
    std::vector<vec4i>& quads_norm, std::vector<vec3f>& norm,
    std::vector<vec4i>& quads_texcoord, std::vector<vec2f>& texcoord,
    int tesselation) {
    log_error("fix implementation");
#if 0
    make_quads(quads_pos, pos, pow2(tesselation + 2), pow2(tesselation + 1),
        [](auto uv) {
            auto a = vec2f{2 * pi * uv.x, pi * (1 - uv.y)};
            return vec3f{cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
        },
        true, false, true, true);
    make_quads(quads_norm, norm, pow2(tesselation + 2), pow2(tesselation + 1),
        [](auto uv) {
            auto a = vec2f{2 * pi * uv.x, pi * (1 - uv.y)};
            return vec3f{cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
        },
        true, false, true, true);
    make_quads(quads_texcoord, texcoord, pow2(tesselation + 2),
        pow2(tesselation + 1), [](auto uv) { return uv; });
#endif
}

// Make a suzanne monkey model for testing. Note that some quads are
// degenerate.
void make_suzanne(
    std::vector<vec4i>& quads, std::vector<vec3f>& pos, int tesselation) {
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
    for (auto& t : suzanne_triangles) { quads.push_back({t.x, t.y, t.z, t.z}); }

    for (auto l = 0; l < tesselation; l++) subdivide_quads(quads, pos);
}

// Generate lines set along a quad.
void make_lines(std::vector<vec2i>& lines, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<float>& radius, const vec2i& steps, const vec2f& size,
    const vec2f& uvsize, const vec2f& line_radius) {
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
        pos[i] = (next_rand3f(rng) - vec3f{0.5f, 0.5f, 0.5f}) * size;
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
    const std::vector<vec3i>& striangles, const std::vector<vec4i>& squads,
    const std::vector<vec3f>& spos, const std::vector<vec3f>& snorm,
    const std::vector<vec2f>& stexcoord, const make_hair_params& params) {
    std::vector<vec3f> bpos;
    std::vector<vec3f> bnorm;
    std::vector<vec2f> btexcoord;
    auto all_triangles = convert_quads_to_triangles(squads);
    all_triangles.insert(
        all_triangles.end(), striangles.begin(), striangles.end());
    std::tie(bpos, bnorm, btexcoord) = sample_triangles_points(
        all_triangles, spos, snorm, stexcoord, steps.y, params.seed);

    auto rng = make_rng(params.seed, 3);
    auto blen = std::vector<float>(bpos.size());
    for (auto& l : blen)
        l = lerp(params.length.x, params.length.y, next_rand1f(rng));

    auto cidx = std::vector<int>();
    if (params.clump.x > 0) {
        for (auto bidx = 0; bidx < bpos.size(); bidx++) {
            cidx.push_back(0);
            auto cdist = flt_max;
            for (auto c = 0; c < params.clump.y; c++) {
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
        radius[i] = lerp(params.radius.x, params.radius.y, u);
        if (params.clump.x > 0) {
            pos[i] = pos[i] +
                     (pos[i + (cidx[bidx] - bidx) * (steps.x + 1)] - pos[i]) *
                         u * params.clump.x;
        }
        if (params.noise.x > 0) {
            auto nx = perlin_noise(pos[i] * params.noise.y + vec3f{0, 0, 0}) *
                      params.noise.x;
            auto ny = perlin_noise(pos[i] * params.noise.y + vec3f{3, 7, 11}) *
                      params.noise.x;
            auto nz =
                perlin_noise(pos[i] * params.noise.y + vec3f{13, 17, 19}) *
                params.noise.x;
            pos[i] += {nx, ny, nz};
        }
    }

    if (params.clump.x > 0 || params.noise.x > 0 || params.rotation.x > 0)
        compute_tangents(lines, pos, norm);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR IMAGEIO
// -----------------------------------------------------------------------------
namespace ygl {

// Pfm load
float* load_pfm(const char* filename, int* w, int* h, int* nc, int req) {
    auto split = [](const std::string& str) {
        auto ret = std::vector<std::string>();
        if (str.empty()) return ret;
        auto lpos = (size_t)0;
        while (lpos != str.npos) {
            auto pos = str.find_first_of(" \t\n\r", lpos);
            if (pos != str.npos) {
                if (pos > lpos) ret.push_back(str.substr(lpos, pos - lpos));
                lpos = pos + 1;
            } else {
                if (lpos < str.size()) ret.push_back(str.substr(lpos));
                lpos = pos;
            }
        }
        return ret;
    };

    auto f = fopen(filename, "rb");
    if (!f) return nullptr;

    // buffer
    char buf[256];
    auto toks = std::vector<std::string>();

    // read magic
    if (!fgets(buf, 256, f)) return nullptr;
    toks = split(buf);
    if (toks[0] == "Pf")
        *nc = 1;
    else if (toks[0] == "PF")
        *nc = 3;
    else
        return nullptr;

    // read w, h
    if (!fgets(buf, 256, f)) return nullptr;
    toks = split(buf);
    *w = atoi(toks[0].c_str());
    *h = atoi(toks[1].c_str());

    // read scale
    if (!fgets(buf, 256, f)) return nullptr;
    toks = split(buf);
    auto s = atof(toks[0].c_str());

    // read the data (flip y)
    auto npixels = (*w) * (*h);
    auto nvalues = (*w) * (*h) * (*nc);
    auto nrow = (*w) * (*nc);
    auto pixels = std::unique_ptr<float[]>(new float[nvalues]);
    for (auto j = *h - 1; j >= 0; j--) {
        if (fread(pixels.get() + j * nrow, sizeof(float), nrow, f) != nrow)
            return nullptr;
    }

    // done reading
    fclose(f);

    // endian conversion
    if (s > 0) {
        for (auto i = 0; i < nvalues; ++i) {
            auto dta = (uint8_t*)(pixels.get() + i);
            std::swap(dta[0], dta[3]);
            std::swap(dta[1], dta[2]);
        }
    }

    // scale
    auto scl = (s > 0) ? s : -s;
    if (scl != 1) {
        for (auto i = 0; i < nvalues; i++) pixels[i] *= scl;
    }

    // proper number of channels
    if (!req || *nc == req) return pixels.release();

    // pack into channels
    if (req < 0 || req > 4) return nullptr;
    auto cpixels = new float[req * npixels];
    for (auto i = 0; i < npixels; i++) {
        auto vp = pixels.get() + i * (*nc);
        auto cp = cpixels + i * req;
        if (*nc == 1) {
            switch (req) {
                case 1: cp[0] = vp[0]; break;
                case 2:
                    cp[0] = vp[0];
                    cp[1] = vp[0];
                    break;
                case 3:
                    cp[0] = vp[0];
                    cp[1] = vp[0];
                    cp[2] = vp[0];
                    break;
                case 4:
                    cp[0] = vp[0];
                    cp[1] = vp[0];
                    cp[2] = vp[0];
                    cp[3] = 1;
                    break;
            }
        } else {
            switch (req) {
                case 1: cp[0] = vp[0]; break;
                case 2:
                    cp[0] = vp[0];
                    cp[1] = vp[1];
                    break;
                case 3:
                    cp[0] = vp[0];
                    cp[1] = vp[1];
                    cp[2] = vp[2];
                    break;
                case 4:
                    cp[0] = vp[0];
                    cp[1] = vp[1];
                    cp[2] = vp[2];
                    cp[3] = 1;
                    break;
            }
        }
    }
    return cpixels;
}

// save pfm
bool save_pfm(const char* filename, int w, int h, int nc, const float* pixels) {
    auto f = fopen(filename, "wb");
    if (!f) return false;

    fprintf(f, "%s\n", (nc == 1) ? "Pf" : "PF");
    fprintf(f, "%d %d\n", w, h);
    fprintf(f, "-1\n");
    if (nc == 1 || nc == 3) {
        fwrite(pixels, sizeof(float), w * h * nc, f);
    } else {
        for (auto i = 0; i < w * h; i++) {
            auto vz = 0.0f;
            auto v = pixels + i * nc;
            fwrite(v + 0, sizeof(float), 1, f);
            fwrite(v + 1, sizeof(float), 1, f);
            if (nc == 2)
                fwrite(&vz, sizeof(float), 1, f);
            else
                fwrite(v + 2, sizeof(float), 1, f);
        }
    }

    fclose(f);

    return true;
}

// check hdr extensions
bool is_hdr_filename(const std::string& filename) {
    auto ext = path_extension(filename);
    return ext == ".hdr" || ext == ".exr" || ext == ".pfm";
}

// Loads an ldr image.
image4b load_image4b(const std::string& filename) {
    auto w = 0, h = 0, c = 0;
    auto pixels =
        std::unique_ptr<byte>(stbi_load(filename.c_str(), &w, &h, &c, 4));
    if (!pixels) return {};
    return make_image4b(w, h, (vec4b*)pixels.get());
}

// Loads an hdr image.
image4f load_image4f(const std::string& filename) {
    auto ext = path_extension(filename);
    auto w = 0, h = 0, c = 0;
    auto pixels = std::unique_ptr<float>(nullptr);
    if (ext == ".exr") {
        auto pixels_ = (float*)nullptr;
        if (LoadEXR(&pixels_, &w, &h, filename.c_str(), nullptr) < 0) return {};
        pixels = std::unique_ptr<float>(pixels_);
        c = 4;
    } else if (ext == ".pfm") {
        pixels =
            std::unique_ptr<float>(load_pfm(filename.c_str(), &w, &h, &c, 4));
    } else {
        pixels =
            std::unique_ptr<float>(stbi_loadf(filename.c_str(), &w, &h, &c, 4));
    }
    if (!pixels) return {};
    return make_image4f(w, h, (vec4f*)pixels.get());
}

// Saves an ldr image.
bool save_image4b(const std::string& filename, const image4b& img) {
    if (path_extension(filename) == ".png") {
        return stbi_write_png(filename.c_str(), img.width, img.height, 4,
            (byte*)img.pixels.data(), img.width * 4);
    } else if (path_extension(filename) == ".jpg") {
        return stbi_write_jpg(filename.c_str(), img.width, img.height, 4,
            (byte*)img.pixels.data(), 75);
    } else {
        return false;
    }
}

// Saves an hdr image.
bool save_image4f(const std::string& filename, const image4f& img) {
    if (path_extension(filename) == ".hdr") {
        return stbi_write_hdr(filename.c_str(), img.width, img.height, 4,
            (float*)img.pixels.data());
    } else if (path_extension(filename) == ".pfm") {
        return save_pfm(filename.c_str(), img.width, img.height, 4,
            (float*)img.pixels.data());
    } else if (path_extension(filename) == ".exr") {
        return !SaveEXR((float*)img.pixels.data(), img.width, img.height, 4,
            filename.c_str());
    } else {
        return false;
    }
}

// Loads an image
std::vector<float> load_imagef(
    const std::string& filename, int& width, int& height, int& ncomp) {
    auto ext = path_extension(filename);
    auto pixels = (float*)nullptr;
    if (ext == ".exr") {
        if (LoadEXR(&pixels, &width, &height, filename.c_str(), nullptr) < 0)
            return {};
        ncomp = 4;
    } else if (ext == ".pfm") {
        pixels = load_pfm(filename.c_str(), &width, &height, &ncomp, 0);
    } else {
        pixels = stbi_loadf(filename.c_str(), &width, &height, &ncomp, 0);
    }
    if (!pixels) return {};
    auto ret = std::vector<float>(pixels, pixels + width * height * ncomp);
    free(pixels);
    return ret;
}

// Loads an image
std::vector<byte> load_imageb(
    const std::string& filename, int& width, int& height, int& ncomp) {
    auto pixels = stbi_load(filename.c_str(), &width, &height, &ncomp, 0);
    if (!pixels) return {};
    auto ret = std::vector<byte>(pixels, pixels + width * height * ncomp);
    free(pixels);
    return ret;
}

// Loads an image from memory.
std::vector<float> load_imagef_from_memory(const std::string& filename,
    const byte* data, int length, int& width, int& height, int& ncomp) {
    auto pixels =
        stbi_loadf_from_memory(data, length, &width, &height, &ncomp, 0);
    if (!pixels) return {};
    auto ret = std::vector<float>(pixels, pixels + width * height * ncomp);
    free(pixels);
    return ret;
}

// Loads an image from memory.
std::vector<byte> load_imageb_from_memory(const std::string& filename,
    const byte* data, int length, int& width, int& height, int& ncomp) {
    auto pixels =
        stbi_load_from_memory(data, length, &width, &height, &ncomp, 0);
    if (!pixels) return {};
    auto ret = std::vector<byte>(pixels, pixels + width * height * ncomp);
    free(pixels);
    return ret;
}

// Saves an image
bool save_imagef(const std::string& filename, int width, int height, int ncomp,
    const float* hdr) {
    if (path_extension(filename) == ".hdr") {
        return stbi_write_hdr(filename.c_str(), width, height, ncomp, hdr);
    } else if (path_extension(filename) == ".pfm") {
        return save_pfm(filename.c_str(), width, height, ncomp, hdr);
    } else {
        return false;
    }
}

// Saves an image
bool save_imageb(const std::string& filename, int width, int height, int ncomp,
    const byte* ldr) {
    if (path_extension(filename) == ".png") {
        return stbi_write_png(
            filename.c_str(), width, height, ncomp, ldr, width * ncomp);
    } else if (path_extension(filename) == ".jpg") {
        return stbi_write_jpg(filename.c_str(), width, height, ncomp, ldr, 75);
    } else {
        return false;
    }
}

// Save an HDR or LDR image with tonemapping based on filename
bool save_image(const std::string& filename, const image4f& hdr,
    tonemap_type tonemapper, float exposure) {
    if (is_hdr_filename(filename)) {
        return save_image4f(filename, hdr);
    } else {
        auto ldr = tonemap_image(hdr, tonemapper, exposure);
        return save_image4b(filename, ldr);
    }
}

// Resize image.
void resize_image(const image4f& img, image4f& res_img, resize_filter filter,
    resize_edge edge, bool premultiplied_alpha) {
    static const auto filter_map = std::map<resize_filter, stbir_filter>{
        {resize_filter::def, STBIR_FILTER_DEFAULT},
        {resize_filter::box, STBIR_FILTER_BOX},
        {resize_filter::triangle, STBIR_FILTER_TRIANGLE},
        {resize_filter::cubic_spline, STBIR_FILTER_CUBICBSPLINE},
        {resize_filter::catmull_rom, STBIR_FILTER_CATMULLROM},
        {resize_filter::mitchell, STBIR_FILTER_MITCHELL}};

    static const auto edge_map =
        std::map<resize_edge, stbir_edge>{{resize_edge::def, STBIR_EDGE_CLAMP},
            {resize_edge::clamp, STBIR_EDGE_CLAMP},
            {resize_edge::reflect, STBIR_EDGE_REFLECT},
            {resize_edge::wrap, STBIR_EDGE_WRAP},
            {resize_edge::zero, STBIR_EDGE_ZERO}};

    stbir_resize_float_generic((float*)img.pixels.data(), img.width, img.height,
        sizeof(vec4f) * img.width, (float*)res_img.pixels.data(), res_img.width,
        res_img.height, sizeof(vec4f) * res_img.width, 4, 3,
        (premultiplied_alpha) ? STBIR_FLAG_ALPHA_PREMULTIPLIED : 0,
        edge_map.at(edge), filter_map.at(filter), STBIR_COLORSPACE_LINEAR,
        nullptr);
}

// Resize image.
void resize_image(const image4b& img, image4b& res_img, resize_filter filter,
    resize_edge edge, bool premultiplied_alpha) {
    static const auto filter_map = std::map<resize_filter, stbir_filter>{
        {resize_filter::def, STBIR_FILTER_DEFAULT},
        {resize_filter::box, STBIR_FILTER_BOX},
        {resize_filter::triangle, STBIR_FILTER_TRIANGLE},
        {resize_filter::cubic_spline, STBIR_FILTER_CUBICBSPLINE},
        {resize_filter::catmull_rom, STBIR_FILTER_CATMULLROM},
        {resize_filter::mitchell, STBIR_FILTER_MITCHELL}};

    static const auto edge_map =
        std::map<resize_edge, stbir_edge>{{resize_edge::def, STBIR_EDGE_CLAMP},
            {resize_edge::clamp, STBIR_EDGE_CLAMP},
            {resize_edge::reflect, STBIR_EDGE_REFLECT},
            {resize_edge::wrap, STBIR_EDGE_WRAP},
            {resize_edge::zero, STBIR_EDGE_ZERO}};

    stbir_resize_uint8_generic((unsigned char*)img.pixels.data(), img.width,
        img.height, sizeof(vec4b) * img.width,
        (unsigned char*)res_img.pixels.data(), res_img.width, res_img.height,
        sizeof(vec4b) * res_img.width, 4, 3,
        (premultiplied_alpha) ? STBIR_FLAG_ALPHA_PREMULTIPLIED : 0,
        edge_map.at(edge), filter_map.at(filter), STBIR_COLORSPACE_LINEAR,
        nullptr);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF IMAGE OPERATIONS
// -----------------------------------------------------------------------------
namespace ygl {

inline vec3f tonemap_gamma(const vec3f& x) {
    return {pow(x.x, 1 / 2.2f), pow(x.y, 1 / 2.2f), pow(x.z, 1 / 2.2f)};
}

inline float tonemap_srgb(float x) {
    if (x <= 0.0031308f) return 12.92f * x;
    return 1.055f * pow(x, 1 / 2.4f) - 0.055f;
}

inline vec3f tonemap_srgb(const vec3f& x) {
    return {tonemap_srgb(x.x), tonemap_srgb(x.y), tonemap_srgb(x.z)};
}

inline vec3f tonemap_filmic1(const vec3f& hdr) {
    // http://filmicworlds.com/blog/filmic-tonemapping-operators/
    auto x = vec3f{max(0.0f, hdr.x - 0.004f), max(0.0f, hdr.y - 0.004f),
        max(0.0f, hdr.z - 0.004f)};
    return (x * (6.2f * x + vec3f{0.5f, 0.5f, 0.5f})) /
           (x * (6.2f * x + vec3f{1.7f, 1.7f, 1.7f}) +
               vec3f{0.06f, 0.06f, 0.06f});
}

inline vec3f tonemap_filmic2(const vec3f& hdr) {
    // https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
    auto x = hdr;
    // x *= 0.6; // brings it back to ACES range
    x = (x * (2.51f * x + vec3f{0.03f, 0.03f, 0.03f})) /
        (x * (2.43f * x + vec3f{0.59f, 0.59f, 0.59f}) +
            vec3f{0.14f, 0.14f, 0.14f});
    return tonemap_gamma(x);
}

inline vec3f tonemap_filmic3(const vec3f& hdr) {
    // https://github.com/TheRealMJP/BakingLab/blob/master/BakingLab/ACES.hlsl

    // sRGB => XYZ => D65_2_D60 => AP1 => RRT_SAT
    static const mat3f ACESInputMat = transpose(mat3f{
        vec3f{0.59719, 0.35458, 0.04823}, vec3f{0.07600, 0.90834, 0.01566},
        vec3f{0.02840, 0.13383, 0.83777}});

    // ODT_SAT => XYZ => D60_2_D65 => sRGB
    static const mat3f ACESOutputMat = transpose(mat3f{
        vec3f{1.60475, -0.53108, -0.07367}, vec3f{-0.10208, 1.10813, -0.00605},
        vec3f{-0.00327, -0.07276, 1.07602}});

    auto x = hdr;
    x = 2 * x;  // matches standard range
    x = ACESInputMat * x;
    // Apply RRT and ODT
    vec3f a = x * (x + vec3f{0.0245786f, 0.0245786f, 0.0245786f}) -
              vec3f{0.000090537f, 0.000090537f, 0.000090537f};
    vec3f b = x * (0.983729f * x + vec3f{0.4329510f, 0.4329510f, 0.4329510f}) +
              vec3f{0.238081f, 0.238081f, 0.238081f};
    x = a / b;
    x = ACESOutputMat * x;
    return tonemap_gamma(x);
}

// Tone mapping HDR to LDR images.
image4b tonemap_image(
    const image4f& hdr, tonemap_type tonemapper, float exposure) {
    auto ldr = make_image4b(hdr.width, hdr.height);
    auto scale = pow(2.0f, exposure);
    for (auto j = 0; j < hdr.height; j++) {
        for (auto i = 0; i < hdr.width; i++) {
            auto h4 = hdr.at(i, j);
            auto h = vec3f{h4.x, h4.y, h4.z} * scale;
            auto a = h4.w;
            switch (tonemapper) {
                case tonemap_type::linear: break;
                case tonemap_type::gamma: h = tonemap_gamma(h); break;
                case tonemap_type::srgb: h = tonemap_srgb(h); break;
                case tonemap_type::filmic1: h = tonemap_filmic1(h); break;
                case tonemap_type::filmic2: h = tonemap_filmic2(h); break;
                case tonemap_type::filmic3: h = tonemap_filmic3(h); break;
            }
            ldr.at(i, j) = float_to_byte({h.x, h.y, h.z, a});
        }
    }
    return ldr;
}

// Image over operator
void image_over(
    vec4f* img, int width, int height, int nlayers, vec4f** layers) {
    for (auto i = 0; i < width * height; i++) {
        img[i] = {0, 0, 0, 0};
        auto weight = 1.0f;
        for (auto l = 0; l < nlayers; l++) {
            img[i].x += layers[l][i].x * layers[l][i].w * weight;
            img[i].y += layers[l][i].y * layers[l][i].w * weight;
            img[i].z += layers[l][i].z * layers[l][i].w * weight;
            img[i].w += layers[l][i].w * weight;
            weight *= (1 - layers[l][i].w);
        }
        if (img[i].w) {
            img[i].x /= img[i].w;
            img[i].y /= img[i].w;
            img[i].z /= img[i].w;
        }
    }
}

// Image over operator
void image_over(
    vec4b* img, int width, int height, int nlayers, vec4b** layers) {
    for (auto i = 0; i < width * height; i++) {
        auto comp = zero4f;
        auto weight = 1.0f;
        for (auto l = 0; l < nlayers && weight > 0; l++) {
            auto w = layers[l][i].w / 255.0f;
            comp.x += layers[l][i].x / 255.0f * w * weight;
            comp.y += layers[l][i].y / 255.0f * w * weight;
            comp.z += layers[l][i].z / 255.0f * w * weight;
            comp.w += w * weight;
            weight *= (1 - w);
        }
        if (comp.w) {
            img[i] = float_to_byte(
                {comp.x / comp.w, comp.y / comp.w, comp.z / comp.w, comp.w});
        } else {
            img[i] = {0, 0, 0, 0};
        }
    }
}

// Convert HSV to RGB
// Implementatkion from
// http://stackoverflow.com/questions/3018313/algorithm-to-convert-rgb-to-hsv-and-hsv-to-rgb-in-range-0-255-for-both
vec4b hsv_to_rgb(const vec4b& hsv) {
    vec4b rgb = {0, 0, 0, hsv.w};
    byte region, remainder, p, q, t;

    byte h = hsv.x, s = hsv.y, v = hsv.z;

    if (s == 0) {
        rgb.x = v;
        rgb.y = v;
        rgb.z = v;
        return rgb;
    }

    region = h / 43;
    remainder = (h - (region * 43)) * 6;

    p = (v * (255 - s)) >> 8;
    q = (v * (255 - ((s * remainder) >> 8))) >> 8;
    t = (v * (255 - ((s * (255 - remainder)) >> 8))) >> 8;

    switch (region) {
        case 0:
            rgb.x = v;
            rgb.y = t;
            rgb.z = p;
            break;
        case 1:
            rgb.x = q;
            rgb.y = v;
            rgb.z = p;
            break;
        case 2:
            rgb.x = p;
            rgb.y = v;
            rgb.z = t;
            break;
        case 3:
            rgb.x = p;
            rgb.y = q;
            rgb.z = v;
            break;
        case 4:
            rgb.x = t;
            rgb.y = p;
            rgb.z = v;
            break;
        default:
            rgb.x = v;
            rgb.y = p;
            rgb.z = q;
            break;
    }

    return rgb;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR IMAGE EXAMPLES
// -----------------------------------------------------------------------------
namespace ygl {

// Make a grid image
image4b make_grid_image(
    int width, int height, int tile, const vec4b& c0, const vec4b& c1) {
    auto img = make_image4b(width, height);
    for (int j = 0; j < width; j++) {
        for (int i = 0; i < height; i++) {
            auto c = i % tile == 0 || i % tile == tile - 1 || j % tile == 0 ||
                     j % tile == tile - 1;
            img.at(i, j) = (c) ? c0 : c1;
        }
    }
    return img;
}

// Make a checkerboard image
image4b make_checker_image(
    int width, int height, int tile, const vec4b& c0, const vec4b& c1) {
    auto img = make_image4b(width, height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto c = (i / tile + j / tile) % 2 == 0;
            img.at(i, j) = (c) ? c0 : c1;
        }
    }
    return img;
}

// Make an image with bumps and dimples.
image4b make_bumpdimple_image(int width, int height, int tile) {
    auto img = make_image4b(width, height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto c = (i / tile + j / tile) % 2 == 0;
            auto ii = i % tile - tile / 2, jj = j % tile - tile / 2;
            auto r =
                sqrt(float(ii * ii + jj * jj)) / sqrt(float(tile * tile) / 4);
            auto h = 0.5f;
            if (r < 0.5f) { h += (c) ? (0.5f - r) : -(0.5f - r); }
            img.at(i, j) = float_to_byte({h, h, h, 1});
        }
    }
    return img;
}

// Make a uv colored grid
image4b make_ramp_image(
    int width, int height, const vec4b& c0, const vec4b& c1, bool srgb) {
    auto img = make_image4b(width, height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto u = (float)i / (float)width;
            if (srgb) {
                img.at(i, j) = linear_to_srgb(
                    srgb_to_linear(c0) * (1 - u) + srgb_to_linear(c1) * u);
            } else {
                img.at(i, j) = float_to_byte(
                    byte_to_float(c0) * (1 - u) + byte_to_float(c1) * u);
            }
        }
    }
    return img;
}

// Make a gamma ramp image
image4b make_gammaramp_image(int width, int height) {
    auto img = make_image4b(width, height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto u = j / float(height - 1);
            if (i < width / 3) u = pow(u, 2.2f);
            if (i > (width * 2) / 3) u = pow(u, 1 / 2.2f);
            auto c = (unsigned char)(u * 255);
            img.at(i, j) = {c, c, c, 255};
        }
    }
    return img;
}

// Make a gamma ramp image
image4f make_gammaramp_imagef(int width, int height) {
    auto img = make_image4f(width, height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto u = j / float(height - 1);
            if (i < width / 3) u = pow(u, 2.2f);
            if (i > (width * 2) / 3) u = pow(u, 1 / 2.2f);
            img.at(i, j) = {u, u, u, 1};
        }
    }
    return img;
}

// Make an image color with red/green in the [0,1] range. Helpful to visualize
// uv texture coordinate application.
image4b make_uv_image(int width, int height) {
    auto img = make_image4b(width, height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            img.at(i, j) = float_to_byte(
                {i / (float)(width - 1), j / (float)(height - 1), 0, 255});
        }
    }
    return img;
}

// Make a uv colored grid
image4b make_uvgrid_image(int width, int height, int tile, bool colored) {
    auto img = make_image4b(width, height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            byte ph = 32 * (i / (height / 8));
            byte pv = 128;
            byte ps = 64 + 16 * (7 - j / (height / 8));
            if (i % (tile / 2) && j % (tile / 2)) {
                if ((i / tile + j / tile) % 2)
                    pv += 16;
                else
                    pv -= 16;
            } else {
                pv = 196;
                ps = 32;
            }
            img.at(i, j) = (colored) ? hsv_to_rgb({ph, ps, pv, 255}) :
                                       vec4b{pv, pv, pv, 255};
        }
    }
    return img;
}

// Make a uv recusive colored grid
image4b make_recuvgrid_image(int width, int height, int tile, bool colored) {
    auto img = make_image4b(width, height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            byte ph = 32 * (i / (height / 8));
            byte pv = 128;
            byte ps = 64 + 16 * (7 - j / (height / 8));
            if (i % (tile / 2) && j % (tile / 2)) {
                if ((i / tile + j / tile) % 2)
                    pv += 16;
                else
                    pv -= 16;
                if ((i / (tile / 4) + j / (tile / 4)) % 2)
                    pv += 4;
                else
                    pv -= 4;
                if ((i / (tile / 8) + j / (tile / 8)) % 2)
                    pv += 1;
                else
                    pv -= 1;
            } else {
                pv = 196;
                ps = 32;
            }
            img.at(i, j) = (colored) ? hsv_to_rgb({ph, ps, pv, 255}) :
                                       vec4b{pv, pv, pv, 255};
        }
    }
    return img;
}

// Comvert a bump map to a normal map.
image4b bump_to_normal_map(const image4b& img, float scale) {
    auto norm = make_image4b(img.width, img.height);
    for (int j = 0; j < img.height; j++) {
        for (int i = 0; i < img.width; i++) {
            auto i1 = (i + 1) % img.width, j1 = (j + 1) % img.height;
            auto p00 = img.at(i, j), p10 = img.at(i1, j), p01 = img.at(i, j1);
            auto g00 = (float(p00.x) + float(p00.y) + float(p00.z)) / (3 * 255);
            auto g01 = (float(p01.x) + float(p01.y) + float(p01.z)) / (3 * 255);
            auto g10 = (float(p10.x) + float(p10.y) + float(p10.z)) / (3 * 255);
            auto n = vec3f{scale * (g00 - g10), scale * (g00 - g01), 1.0f};
            n = normalize(n) * 0.5f + vec3f{0.5f, 0.5f, 0.5f};
            auto c =
                vec4b{byte(n.x * 255), byte(n.y * 255), byte(n.z * 255), 255};
            norm.at(i, j) = c;
        }
    }
    return norm;
}

// Implementation of sunsky modified heavily from pbrt
image4f make_sunsky_image(
    int res, float thetaSun, float turbidity, bool has_sun, bool has_ground) {
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

    auto img = make_image4f(2 * res, res);
    for (auto j = 0; j < img.height; j++) {
        if (!has_ground && j >= img.height / 2) continue;
        auto theta = pi * ((j + 0.5f) / img.height);
        theta = clamp(theta, 0.0f, pi / 2 - flt_eps);
        for (int i = 0; i < img.width; i++) {
            auto phi = 2 * pi * (float(i + 0.5f) / img.width);
            auto w =
                vec3f{cos(phi) * sin(theta), cos(theta), sin(phi) * sin(theta)};
            auto gamma = acos(clamp(dot(w, wSun), -1.0f, 1.0f));
            auto col = sky(theta, gamma) + sun(theta, gamma);
            img.at(i, j) = {col.x, col.y, col.z, 1};
        }
    }

    return img;
}

// Make a noise image. Wrap works only if both resx and resy are powers of two.
image4b make_noise_image(int resx, int resy, float scale, bool wrap) {
    auto wrap3i = (wrap) ? vec3i{resx, resy, 2} : zero3i;
    auto img = make_image4b(resx, resy);
    for (auto j = 0; j < resy; j++) {
        for (auto i = 0; i < resx; i++) {
            auto p = vec3f{i / (float)resx, j / (float)resy, 0.5f} * scale;
            auto g = perlin_noise(p, wrap3i);
            g = clamp(0.5f + 0.5f * g, 0.0f, 1.0f);
            img.at(i, j) = float_to_byte({g, g, g, 1});
        }
    }
    return img;
}

// Make a noise image. Wrap works only if both resx and resy are powers of two.
image4b make_fbm_image(int resx, int resy, float scale, float lacunarity,
    float gain, int octaves, bool wrap) {
    auto wrap3i = (wrap) ? vec3i{resx, resy, 2} : zero3i;
    auto img = make_image4b(resx, resy);
    for (auto j = 0; j < resy; j++) {
        for (auto i = 0; i < resx; i++) {
            auto p = vec3f{i / (float)resx, j / (float)resy, 0.5f} * scale;
            auto g = perlin_fbm_noise(p, lacunarity, gain, octaves, wrap3i);
            g = clamp(0.5f + 0.5f * g, 0.0f, 1.0f);
            img.at(i, j) = float_to_byte({g, g, g, 1});
        }
    }
    return img;
}

// Make a noise image. Wrap works only if both resx and resy are powers of two.
image4b make_ridge_image(int resx, int resy, float scale, float lacunarity,
    float gain, float offset, int octaves, bool wrap) {
    auto wrap3i = (wrap) ? vec3i{resx, resy, 2} : zero3i;
    auto img = make_image4b(resx, resy);
    for (auto j = 0; j < resy; j++) {
        for (auto i = 0; i < resx; i++) {
            auto p = vec3f{i / (float)resx, j / (float)resy, 0.5f} * scale;
            auto g = perlin_ridge_noise(
                p, lacunarity, gain, offset, octaves, wrap3i);
            g = clamp(g, 0.0f, 1.0f);
            img.at(i, j) = float_to_byte({g, g, g, 1});
        }
    }
    return img;
}

// Make a noise image. Wrap works only if both resx and resy are powers of two.
image4b make_turbulence_image(int resx, int resy, float scale, float lacunarity,
    float gain, int octaves, bool wrap) {
    auto wrap3i = (wrap) ? vec3i{resx, resy, 2} : zero3i;
    auto img = make_image4b(resx, resy);
    for (auto j = 0; j < resy; j++) {
        for (auto i = 0; i < resx; i++) {
            auto p = vec3f{i / (float)resx, j / (float)resy, 0.5f} * scale;
            auto g =
                perlin_turbulence_noise(p, lacunarity, gain, octaves, wrap3i);
            g = clamp(g, 0.0f, 1.0f);
            img.at(i, j) = float_to_byte({g, g, g, 1});
        }
    }
    return img;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENRTATION OF RAY-PRIMITIVE INTERSECTION FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Intersect a ray with a point (approximate)
bool intersect_point(const ray3f& ray, const vec3f& p, float r, float& ray_t) {
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
    ray_t = t;

    return true;
}

// Intersect a ray with a line
bool intersect_line(const ray3f& ray, const vec3f& v0, const vec3f& v1,
    float r0, float r1, float& ray_t, vec2f& euv) {
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
    ray_t = t;
    euv = {s, sqrt(d2) / r};

    return true;
}

// Intersect a ray with a triangle
bool intersect_triangle(const ray3f& ray, const vec3f& v0, const vec3f& v1,
    const vec3f& v2, float& ray_t, vec2f& euv) {
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
    ray_t = t;
    euv = {u, v};

    return true;
}

// Intersect a ray with a quad.
bool intersect_quad(const ray3f& ray, const vec3f& v0, const vec3f& v1,
    const vec3f& v2, const vec3f& v3, float& ray_t, vec2f& euv) {
    auto hit = false;
    auto tray = ray;
    if (intersect_triangle(tray, v0, v1, v3, ray_t, euv)) {
        tray.tmax = ray_t;
        hit = true;
    }
    if (intersect_triangle(tray, v2, v3, v1, ray_t, euv)) {
        euv = {1 - euv.x, 1 - euv.y};
        tray.tmax = ray_t;
        hit = true;
    }
    return hit;
}

// Intersect a ray with a tetrahedron.
bool intersect_tetrahedron(const ray3f& ray_, const vec3f& v0, const vec3f& v1,
    const vec3f& v2, const vec3f& v3, float& ray_t, vec4f& euv) {
    // check intersction for each face
    auto hit = false;
    auto ray = ray_;
    auto tuv = zero2f;
    // TODO: fix uvs
    if (intersect_triangle(ray, v0, v1, v2, ray_t, tuv)) {
        hit = true;
        ray.tmax = ray_t;
    }
    if (intersect_triangle(ray, v0, v1, v3, ray_t, tuv)) {
        hit = true;
        ray.tmax = ray_t;
    }
    if (intersect_triangle(ray, v0, v2, v3, ray_t, tuv)) {
        hit = true;
        ray.tmax = ray_t;
    }
    if (intersect_triangle(ray, v1, v2, v3, ray_t, tuv)) {
        hit = true;
        ray.tmax = ray_t;
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
bool intersect_check_bbox(const ray3f& ray, const bbox3f& bbox) {
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
bool intersect_check_bbox(const ray3f& ray, const vec3f& ray_dinv,
    const vec3i& ray_dsign, const bbox3f& bbox_) {
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
bool overlap_point(
    const vec3f& pos, float dist_max, const vec3f& p, float r, float& dist) {
    auto d2 = dot(pos - p, pos - p);
    if (d2 > (dist_max + r) * (dist_max + r)) return false;
    dist = sqrt(d2);
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
    const vec3f& v1, float r0, float r1, float& dist, vec2f& euv) {
    auto u = closestuv_line(pos, v0, v1);
    // Compute projected position from the clamped t d = a + t * ab;
    auto p = v0 + (v1 - v0) * u;
    auto r = r0 + (r1 - r0) * u;
    auto d2 = dot(pos - p, pos - p);
    // check distance
    if (d2 > (dist_max + r) * (dist_max + r)) return false;
    // done
    dist = sqrt(d2);
    euv = {u, 0};
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
    vec2f& euv) {
    auto uv = closestuv_triangle(pos, v0, v1, v2);
    auto p = interpolate_triangle(v0, v1, v2, uv);
    auto r = interpolate_triangle(r0, r1, r2, uv);
    auto dd = dot(p - pos, p - pos);
    if (dd > (dist_max + r) * (dist_max + r)) return false;
    dist = sqrt(dd);
    euv = uv;
    return true;
}

// TODO: documentation
bool overlap_quad(const vec3f& pos, float dist_max, const vec3f& v0,
    const vec3f& v1, const vec3f& v2, const vec3f& v3, float r0, float r1,
    float r2, float r3, float& dist, vec2f& euv) {
    auto hit = false;
    if (overlap_triangle(pos, dist_max, v0, v1, v3, r0, r1, r3, dist, euv)) {
        dist_max = dist;
        hit = true;
    }
    if (overlap_triangle(pos, dist_max, v2, v3, v1, r2, r3, r1, dist, euv)) {
        // dist_max = dist;
        euv = {1 - euv.x, 1 - euv.y};
        hit = true;
    }
    return hit;
}

// TODO: documentation
bool overlap_tetrahedron(const vec3f& pos, const vec3f& v0, const vec3f& v1,
    const vec3f& v2, const vec3f& v3, vec4f& euv) {
    // TODO: fix uv
    auto vol = dot(v3 - v0, cross(v3 - v1, v3 - v0));
    if (vol == 0) return false;
    auto u = dot(v3 - v0, cross(v3 - v1, v3 - v0)) / vol;
    if (u < 0 || u > 1) return false;
    auto v = dot(v3 - v0, cross(v3 - v1, v3 - v0)) / vol;
    if (v < 0 || v > 1 || u + v > 1) return false;
    auto w = dot(v3 - v0, cross(v3 - v1, v3 - v0)) / vol;
    if (w < 0 || w > 1 || u + v + w > 1) return false;
    euv = {u, v, w, 1 - u - v - w};
    return true;
}

// TODO: documentation
bool overlap_tetrahedron(const vec3f& pos, float dist_max, const vec3f& v0,
    const vec3f& v1, const vec3f& v2, const vec3f& v3, float r0, float r1,
    float r2, float r3, float& dist, vec4f& euv) {
    // TODO: FIX UVs
    // check interior
    if (overlap_tetrahedron(pos, v0, v1, v2, v3, euv)) {
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

// Build a BVH node list and sorted primitive array
std::vector<bvh_node> make_bvh_nodes(
    const std::vector<bbox3f>& bboxes, bvh_node_type type, bool equal_size) {
    // create an array of primitives to sort
    auto prims = std::vector<bvh_prim>(bboxes.size());
    for (auto i = 0; i < bboxes.size(); i++) {
        prims[i].bbox = bboxes[i];
        prims[i].center = (bboxes[i].min + bboxes[i].max) / 2;
        prims[i].primid = i;
    }

    // allocate nodes (over-allocate now then shrink)
    auto nodes = std::vector<bvh_node>();
    nodes.reserve(prims.size() * 2);

    // start recursive splitting
    make_bvh_node(nodes, prims, 0, (int)prims.size(), type, equal_size);

    // shrink back
    nodes.shrink_to_fit();

    // done
    return nodes;
}

// Build a BVH from the data already set
void make_bvh_nodes(bvh_tree* bvh, bool equal_size) {
    // get the number of primitives and the primitive type
    auto bboxes = std::vector<bbox3f>();
    if (!bvh->points.empty()) {
        for (auto& p : bvh->points) {
            bboxes.push_back(point_bbox(bvh->pos[p], bvh->radius[p]));
        }
        bvh->type = bvh_node_type::point;
    } else if (!bvh->lines.empty()) {
        for (auto& l : bvh->lines) {
            bboxes.push_back(line_bbox(bvh->pos[l.x], bvh->pos[l.y],
                bvh->radius[l.x], bvh->radius[l.y]));
        }
        bvh->type = bvh_node_type::line;
    } else if (!bvh->triangles.empty()) {
        for (auto& t : bvh->triangles) {
            bboxes.push_back(
                triangle_bbox(bvh->pos[t.x], bvh->pos[t.y], bvh->pos[t.z]));
        }
        bvh->type = bvh_node_type::triangle;
    } else if (!bvh->quads.empty()) {
        for (auto& q : bvh->quads) {
            bboxes.push_back(quad_bbox(
                bvh->pos[q.x], bvh->pos[q.y], bvh->pos[q.z], bvh->pos[q.w]));
        }
        bvh->type = bvh_node_type::quad;
    } else if (!bvh->pos.empty()) {
        for (auto i = 0; i < bvh->pos.size(); i++) {
            bboxes.push_back(point_bbox(bvh->pos[i], bvh->radius[i]));
        }
        bvh->type = bvh_node_type::vertex;
    } else if (!bvh->instances.empty()) {
        for (auto& ist : bvh->instances) {
            auto sbvh = bvh->shape_bvhs[ist.shape_id];
            bboxes.push_back(transform_bbox(ist.frame, sbvh->nodes[0].bbox));
        }
        bvh->type = bvh_node_type::instance;
    }

    // make node bvh
    bvh->nodes = make_bvh_nodes(bboxes, bvh->type, equal_size);
}

// Build a BVH from a set of primitives.
bvh_tree* make_bvh(const std::vector<int>& points,
    const std::vector<vec2i>& lines, const std::vector<vec3i>& triangles,
    const std::vector<vec4i>& quads, const std::vector<vec3f>& pos,
    const std::vector<float>& radius, float def_radius, bool equal_size) {
    // allocate the bvh
    auto bvh = new bvh_tree();

    // set values
    bvh->points = points;
    bvh->lines = lines;
    bvh->triangles = triangles;
    bvh->quads = quads;
    bvh->pos = pos;
    bvh->radius =
        (radius.empty()) ? std::vector<float>(pos.size(), def_radius) : radius;

    // make bvh nodes
    make_bvh_nodes(bvh, equal_size);

    // done
    return bvh;
}

// Build a BVH from a set of shape instances.
bvh_tree* make_bvh(const std::vector<bvh_instance>& instances,
    const std::vector<bvh_tree*>& shape_bvhs, bool equal_size,
    bool own_shape_bvhs) {
    // allocate the bvh
    auto bvh = new bvh_tree();

    // set values
    bvh->instances = instances;
    bvh->shape_bvhs = shape_bvhs;
    bvh->own_shape_bvhs = own_shape_bvhs;

    // make bvh nodes
    make_bvh_nodes(bvh, equal_size);

    // done
    return bvh;
}

// Recursively recomputes the node bounds for a shape bvh
void refit_bvh(bvh_tree* bvh, int nodeid) {
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
        case bvh_node_type::point: {
            for (auto i = 0; i < node.count; i++) {
                auto& p = bvh->points[node.prims[i]];
                node.bbox += point_bbox(bvh->pos[p], bvh->radius[p]);
            }
        } break;
        case bvh_node_type::line: {
            for (auto i = 0; i < node.count; i++) {
                auto& l = bvh->lines[node.prims[i]];
                node.bbox += line_bbox(bvh->pos[l.x], bvh->pos[l.y],
                    bvh->radius[l.x], bvh->radius[l.y]);
            }
        } break;
        case bvh_node_type::triangle: {
            for (auto i = 0; i < node.count; i++) {
                auto& t = bvh->triangles[node.prims[i]];
                node.bbox +=
                    triangle_bbox(bvh->pos[t.x], bvh->pos[t.y], bvh->pos[t.z]);
            }
        } break;
        case bvh_node_type::quad: {
            for (auto i = 0; i < node.count; i++) {
                auto& q = bvh->quads[node.prims[i]];
                node.bbox += quad_bbox(
                    bvh->pos[q.x], bvh->pos[q.y], bvh->pos[q.z], bvh->pos[q.w]);
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
                auto& ist = bvh->instances[node.prims[node.prims[i]]];
                auto sbvh = bvh->shape_bvhs[ist.shape_id];
                node.bbox += transform_bbox(ist.frame, sbvh->nodes[0].bbox);
            }
        } break;
    }
}

// Recursively recomputes the node bounds for a shape bvh
void refit_bvh(bvh_tree* bvh, const std::vector<vec3f>& pos,
    const std::vector<float>& radius, float def_radius) {
    bvh->pos = pos;
    bvh->radius =
        (radius.empty()) ? std::vector<float>(pos.size(), def_radius) : radius;
    refit_bvh(bvh, 0);
}

// Recursively recomputes the node bounds for a scene bvh
void refit_bvh(bvh_tree* bvh, const std::vector<frame3f>& frames,
    const std::vector<frame3f>& frames_inv) {
    for (auto i = 0; i < frames.size(); i++) {
        bvh->instances[i].frame = frames[i];
        bvh->instances[i].frame_inv = frames_inv[i];
    }
    refit_bvh(bvh, 0);
}

// Intersect ray with a bvh.
bool intersect_bvh(const bvh_tree* bvh, const ray3f& ray_, bool find_any,
    float& ray_t, int& iid, int& eid, vec2f& euv) {
    // node stack
    int node_stack[128];
    auto node_cur = 0;
    node_stack[node_cur++] = 0;

    // shared variables
    auto hit = false;

    // copy ray to modify it
    auto ray = ray_;

    // prepare ray for fast queries
    auto ray_dinv = vec3f{1, 1, 1} / ray.d;
    auto ray_dsign = vec3i{(ray_dinv.x < 0) ? 1 : 0, (ray_dinv.y < 0) ? 1 : 0,
        (ray_dinv.z < 0) ? 1 : 0};
    auto ray_reverse = std::array<bool, 4>{
        {(bool)ray_dsign.x, (bool)ray_dsign.y, (bool)ray_dsign.z, false}};

    // walking stack
    while (node_cur) {
        // grab node
        auto& node = bvh->nodes[node_stack[--node_cur]];

        // intersect bbox
        if (!intersect_check_bbox(ray, ray_dinv, ray_dsign, node.bbox))
            continue;

        // intersect node, switching based on node type
        // for each type, iterate over the the primitive list
        switch (node.type) {
            case bvh_node_type::internal: {
                // for internal nodes, attempts to proceed along the
                // split axis from smallest to largest nodes
                if (ray_reverse[node.split_axis]) {
                    node_stack[node_cur++] = node.prims[0];
                    node_stack[node_cur++] = node.prims[1];
                } else {
                    node_stack[node_cur++] = node.prims[1];
                    node_stack[node_cur++] = node.prims[0];
                }
            } break;
            case bvh_node_type::point: {
                for (auto i = 0; i < node.count; i++) {
                    auto& p = bvh->points[node.prims[i]];
                    if (intersect_point(
                            ray, bvh->pos[p], bvh->radius[p], ray_t)) {
                        hit = true;
                        ray.tmax = ray_t;
                        eid = node.prims[i];
                        euv = {1, 0};
                    }
                }
            } break;
            case bvh_node_type::line: {
                for (auto i = 0; i < node.count; i++) {
                    auto& l = bvh->lines[node.prims[i]];
                    if (intersect_line(ray, bvh->pos[l.x], bvh->pos[l.y],
                            bvh->radius[l.x], bvh->radius[l.y], ray_t, euv)) {
                        hit = true;
                        ray.tmax = ray_t;
                        eid = node.prims[i];
                    }
                }
            } break;
            case bvh_node_type::triangle: {
                for (auto i = 0; i < node.count; i++) {
                    auto& t = bvh->triangles[node.prims[i]];
                    if (intersect_triangle(ray, bvh->pos[t.x], bvh->pos[t.y],
                            bvh->pos[t.z], ray_t, euv)) {
                        hit = true;
                        ray.tmax = ray_t;
                        eid = node.prims[i];
                    }
                }
            } break;
            case bvh_node_type::quad: {
                for (auto i = 0; i < node.count; i++) {
                    auto& q = bvh->quads[node.prims[i]];
                    if (intersect_quad(ray, bvh->pos[q.x], bvh->pos[q.y],
                            bvh->pos[q.z], bvh->pos[q.w], ray_t, euv)) {
                        hit = true;
                        ray.tmax = ray_t;
                        eid = node.prims[i];
                    }
                }
            } break;
            case bvh_node_type::vertex: {
                for (auto i = 0; i < node.count; i++) {
                    auto idx = node.prims[i];
                    if (intersect_point(
                            ray, bvh->pos[idx], bvh->radius[idx], ray_t)) {
                        hit = true;
                        ray.tmax = ray_t;
                        eid = node.prims[i];
                        euv = {1, 0};
                    }
                }
            } break;
            case bvh_node_type::instance: {
                for (auto i = 0; i < node.count; i++) {
                    auto& ist = bvh->instances[node.prims[i]];
                    auto sbvh = bvh->shape_bvhs[ist.shape_id];
                    if (intersect_bvh(sbvh, transform_ray(ist.frame_inv, ray),
                            find_any, ray_t, iid, eid, euv)) {
                        hit = true;
                        ray.tmax = ray_t;
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
bool overlap_bvh(const bvh_tree* bvh, const vec3f& pos, float max_dist,
    bool find_any, float& dist, int& iid, int& eid, vec2f& euv) {
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
            case bvh_node_type::point: {
                for (auto i = 0; i < node.count; i++) {
                    auto& p = bvh->points[node.prims[i]];
                    if (overlap_point(
                            pos, max_dist, bvh->pos[p], bvh->radius[p], dist)) {
                        hit = true;
                        max_dist = dist;
                        eid = node.prims[i];
                        euv = {1, 0};
                    }
                }
            } break;
            case bvh_node_type::line: {
                for (auto i = 0; i < node.count; i++) {
                    auto& l = bvh->lines[node.prims[i]];
                    if (overlap_line(pos, max_dist, bvh->pos[l.x],
                            bvh->pos[l.y], bvh->radius[l.x], bvh->radius[l.y],
                            dist, euv)) {
                        hit = true;
                        max_dist = dist;
                        eid = node.prims[i];
                    }
                }
            } break;
            case bvh_node_type::triangle: {
                for (auto i = 0; i < node.count; i++) {
                    auto& t = bvh->triangles[node.prims[i]];
                    if (overlap_triangle(pos, max_dist, bvh->pos[t.x],
                            bvh->pos[t.y], bvh->pos[t.z], bvh->radius[t.x],
                            bvh->radius[t.y], bvh->radius[t.z], dist, euv)) {
                        hit = true;
                        max_dist = dist;
                        eid = node.prims[i];
                    }
                }
            } break;
            case bvh_node_type::quad: {
                for (auto i = 0; i < node.count; i++) {
                    auto& q = bvh->quads[node.prims[i]];
                    if (overlap_quad(pos, max_dist, bvh->pos[q.x],
                            bvh->pos[q.y], bvh->pos[q.z], bvh->pos[q.w],
                            bvh->radius[q.x], bvh->radius[q.y],
                            bvh->radius[q.z], bvh->radius[q.w], dist, euv)) {
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
                            bvh->radius[idx], dist)) {
                        hit = true;
                        max_dist = dist;
                        eid = node.prims[i];
                        euv = {1, 0};
                    }
                }
            } break;
            case bvh_node_type::instance: {
                for (auto i = 0; i < node.count; i++) {
                    auto& ist = bvh->instances[node.prims[i]];
                    auto sbvh = bvh->shape_bvhs[ist.shape_id];
                    if (overlap_bvh(sbvh, transform_point(ist.frame_inv, pos),
                            max_dist, find_any, dist, iid, eid, euv)) {
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

// Intersect ray with a bvh (convenience wrapper).
intersection_point intersect_bvh(
    const bvh_tree* bvh, const ray3f& ray, bool find_any) {
    auto isec = intersection_point();
    if (!intersect_bvh(
            bvh, ray, find_any, isec.dist, isec.iid, isec.eid, isec.euv))
        return {};
    return isec;
}

// Finds the closest element with a bvh (convenience wrapper).
intersection_point overlap_bvh(
    const bvh_tree* bvh, const vec3f& pos, float max_dist, bool find_any) {
    auto isec = intersection_point();
    if (!overlap_bvh(bvh, pos, max_dist, find_any, isec.dist, isec.iid,
            isec.eid, isec.euv))
        return {};
    return isec;
}

#if 0
    // Finds the overlap between BVH leaf nodes.
    template <typename OverlapElem>
    void overlap_bvh_elems(const bvh_tree* bvh1, const bvh_tree* bvh2,
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
                        assert(node_cur < 128);
                    }
                } else if (node2.isleaf) {
                    for (auto idx1 = node1.start; idx1 < node1.start + node1.count;
                         idx1++) {
                        node_stack[node_cur++] = {(int)idx1, node_idx.y};
                        assert(node_cur < 128);
                    }
                } else {
                    for (auto idx2 = node2.start; idx2 < node2.start + node2.count;
                         idx2++) {
                        for (auto idx1 = node1.start;
                             idx1 < node1.start + node1.count; idx1++) {
                            node_stack[node_cur++] = {(int)idx1, (int)idx2};
                            assert(node_cur < 128);
                        }
                    }
                }
            }
        }
    }
#endif

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR SIMPLE SCENE
// -----------------------------------------------------------------------------
namespace ygl {

// Synchronizes shape element type.
shape_elem_type get_shape_type(const shape* shp) {
    if (!shp->triangles.empty()) {
        return shape_elem_type::triangles;
    } else if (!shp->lines.empty()) {
        return shape_elem_type::lines;
    } else if (!shp->points.empty()) {
        return shape_elem_type::points;
    } else if (!shp->quads.empty()) {
        return shape_elem_type::quads;
    } else if (!shp->beziers.empty()) {
        return shape_elem_type::beziers;
    } else if (!shp->quads_pos.empty()) {
        return shape_elem_type::facevarying;
    } else if (!shp->pos.empty()) {
        return shape_elem_type::vertices;
    } else {
        return shape_elem_type::none;
    }
}

// Shape element normal.
vec3f eval_elem_norm(const shape* shp, int eid) {
    auto norm = zero3f;
    switch (get_shape_type(shp)) {
        case shape_elem_type::none: {
            throw std::runtime_error("type not supported");
        } break;
        case shape_elem_type::triangles: {
            auto t = shp->triangles[eid];
            norm = triangle_normal(shp->pos[t.x], shp->pos[t.y], shp->pos[t.z]);
        } break;
        case shape_elem_type::lines: {
            auto l = shp->lines[eid];
            norm = line_tangent(shp->pos[l.x], shp->pos[l.y]);
        } break;
        case shape_elem_type::points: {
            norm = {0, 0, 1};
        } break;
        case shape_elem_type::quads: {
            auto q = shp->quads[eid];
            norm = quad_normal(
                shp->pos[q.x], shp->pos[q.y], shp->pos[q.z], shp->pos[q.w]);
        } break;
        case shape_elem_type::beziers: {
            auto l = shp->beziers[eid];
            norm = line_tangent(shp->pos[l.x], shp->pos[l.w]);
        } break;
        case shape_elem_type::vertices: {
            norm = {0, 0, 1};
        } break;
        case shape_elem_type::facevarying: {
            auto q = (shp->quads_norm.empty()) ? shp->quads_pos[eid] :
                                                 shp->quads_norm[eid];
            norm = quad_normal(
                shp->pos[q.x], shp->pos[q.y], shp->pos[q.z], shp->pos[q.w]);
        } break;
    }
    return norm;
}

// Shape value interpolated using barycentric coordinates
template <typename T>
T eval_elem(const shape* shp, const std::vector<T>& vals, int eid,
    const vec2f& euv, const T& def) {
    if (vals.empty()) return def;
    switch (get_shape_type(shp)) {
        case shape_elem_type::none: {
            throw std::runtime_error("type not supported");
        } break;
        case shape_elem_type::triangles: {
            return interpolate_triangle(vals, shp->triangles[eid], euv);
        } break;
        case shape_elem_type::lines: {
            return interpolate_line(vals, shp->lines[eid], euv.x);
        } break;
        case shape_elem_type::points: {
            return interpolate_point(vals, shp->points[eid]);
        } break;
        case shape_elem_type::quads: {
            return interpolate_quad(vals, shp->quads[eid], euv);
        } break;
        case shape_elem_type::beziers: {
            return interpolate_bezier(vals, shp->beziers[eid], euv.x);
        } break;
        case shape_elem_type::vertices: {
            return vals[eid];
        } break;
        case shape_elem_type::facevarying: {
            throw std::runtime_error("should not have gotten here");
            return def;
        } break;
    }
}

// Shape position interpolated using barycentric coordinates
vec3f eval_pos(const shape* shp, int eid, const vec2f& euv) {
    return eval_elem(shp, shp->pos, eid, euv, {0, 0, 0});
}
// Shape normal interpolated using barycentric coordinates
vec3f eval_norm(const shape* shp, int eid, const vec2f& euv) {
    if (shp->norm.empty()) return eval_elem_norm(shp, eid);
    return normalize(eval_elem(shp, shp->norm, eid, euv, {0, 0, 1}));
}
// Shape texcoord interpolated using barycentric coordinates
vec2f eval_texcoord(const shape* shp, int eid, const vec2f& euv) {
    return eval_elem(shp, shp->texcoord, eid, euv, {0, 0});
}
// Shape color interpolated using barycentric coordinates
vec4f eval_color(const shape* shp, int eid, const vec2f& euv) {
    return eval_elem(shp, shp->color, eid, euv, {1, 1, 1, 1});
}
// Shape radius interpolated using barycentric coordinates
float eval_radius(const shape* shp, int eid, const vec2f& euv) {
    return eval_elem(shp, shp->radius, eid, euv, 0.0f);
}
// Shape tangent space interpolated using barycentric coordinates
vec4f eval_tangsp(const shape* shp, int eid, const vec2f& euv) {
    return eval_elem(shp, shp->tangsp, eid, euv, {0, 0, 0, 1});
}

// Environment position interpolated using uv parametrization.
vec3f eval_pos(const environment* env, const vec2f& uv, bool transformed) {
    auto pos = sphericaly_to_cartesian(
        vec3f{uv.x * pi * 2, uv.y * pi, environment_distance});
    if (transformed) pos = transform_point(env->frame, pos);
    return pos;
}
// Environment normal interpolated using uv parametrization.
vec3f eval_norm(const environment* env, const vec2f& uv, bool transformed) {
    auto norm = normalize(
        -sphericaly_to_cartesian(vec3f{uv.x * pi * 2, uv.y * pi, 1.0f}));
    if (transformed) norm = transform_direction(env->frame, norm);
    return norm;
}
// Environment texture coordinates from uv parametrization.
vec2f eval_texcoord(const environment* env, const vec2f& uv) { return uv; }
// Evaluate uv parameters for environment.
vec2f eval_uv(const environment* env, const vec3f& w, bool transformed) {
    auto wl = w;
    if (transformed) wl = transform_direction_inverse(env->frame, wl);
    auto sh = cartesian_to_sphericaly(wl);
    return {sh.x / (2 * pi), sh.y / pi};
}

// Evaluate a texture
vec4f eval_texture(const texture_info& info, const vec2f& texcoord, bool srgb,
    const vec4f& def) {
    auto txt = info.txt;
    if (!txt || (txt->hdr.pixels.empty() && txt->ldr.pixels.empty()))
        return def;

    auto lookup = [&def, &txt, &srgb](int i, int j) {
        if (!txt->ldr.pixels.empty())
            return (srgb) ? srgb_to_linear(txt->ldr.at(i, j)) :
                            byte_to_float(txt->ldr.at(i, j));
        else if (!txt->hdr.pixels.empty())
            return txt->hdr.at(i, j);
        else
            return def;
    };

    // get image width/height
    auto w = (!txt->ldr.pixels.empty()) ? txt->ldr.width : txt->hdr.width,
         h = (!txt->ldr.pixels.empty()) ? txt->ldr.height : txt->hdr.height;

    // get coordinates normalized for tiling
    auto s = 0.0f, t = 0.0f;
    if (!info.wrap_s) {
        s = clamp(texcoord.x, 0.0f, 1.0f) * w;
    } else {
        s = std::fmod(texcoord.x, 1.0f) * w;
        if (s < 0) s += w;
    }
    if (!info.wrap_t) {
        t = clamp(texcoord.y, 0.0f, 1.0f) * h;
    } else {
        t = std::fmod(texcoord.y, 1.0f) * h;
        if (t < 0) t += h;
    }

    // get image coordinates and residuals
    auto i = clamp((int)s, 0, w - 1), j = clamp((int)t, 0, h - 1);
    auto ii = (i + 1) % w, jj = (j + 1) % h;
    auto u = s - i, v = t - j;

    // nearest lookup
    if (!info.linear) return lookup(i, j);

    // handle interpolation
    return lookup(i, j) * (1 - u) * (1 - v) + lookup(i, jj) * (1 - u) * v +
           lookup(ii, j) * u * (1 - v) + lookup(ii, jj) * u * v;
}

// Generates a ray from a camera for image plane coordinate uv and the
// lens coordinates luv.
ray3f eval_camera_ray(const camera* cam, const vec2f& uv, const vec2f& luv) {
    auto h = 2 * tan(cam->yfov / 2);
    auto w = h * cam->aspect;
    auto o = vec3f{luv.x * cam->aperture, luv.y * cam->aperture, 0};
    auto q = vec3f{w * cam->focus * (uv.x - 0.5f),
        h * cam->focus * (uv.y - 0.5f), -cam->focus};
    return make_ray(transform_point(cam->frame, o),
        transform_direction(cam->frame, normalize(q - o)));
}

// Generates a ray from a camera for pixel coordinates `ij`, the resolution
// `res`, the sub-pixel coordinates `puv` and the lens coordinates `luv` and the
// image resolution `res`.
ray3f eval_camera_ray(const camera* cam, const vec2i& ij, int res,
    const vec2f& puv, const vec2f& luv) {
    auto uv =
        vec2f{(ij.x + puv.x) / (cam->aspect * res), 1 - (ij.y - puv.y) / res};
    return eval_camera_ray(cam, uv, luv);
}

// Generate a distribution for sampling a shape uniformly based on area/length.
std::vector<float> sample_shape_cdf(const shape* shp) {
    switch (get_shape_type(shp)) {
        case shape_elem_type::none:
            throw std::runtime_error("type not supported");
        case shape_elem_type::triangles:
            return sample_triangles_cdf(shp->triangles, shp->pos);
        case shape_elem_type::lines:
            return sample_lines_cdf(shp->lines, shp->pos);
        case shape_elem_type::points:
            return sample_points_cdf(shp->points.size());
        case shape_elem_type::quads:
            return sample_quads_cdf(shp->quads, shp->pos);
        case shape_elem_type::beziers:
            throw std::runtime_error("type not supported");
        case shape_elem_type::vertices:
            return sample_points_cdf(shp->pos.size());
        case shape_elem_type::facevarying:
            return sample_quads_cdf(shp->quads_pos, shp->pos);
    }
}

// Sample a shape based on a distribution.
std::pair<int, vec2f> sample_shape(const shape* shp,
    const std::vector<float>& cdf, float re, const vec2f& ruv) {
    switch (get_shape_type(shp)) {
        case shape_elem_type::none:
            throw std::runtime_error("type not supported");
        case shape_elem_type::triangles: return sample_triangles(cdf, re, ruv);
        case shape_elem_type::lines:
            return {sample_lines(cdf, re, ruv.x).first, ruv};
        case shape_elem_type::points: return {sample_points(cdf, re), ruv};
        case shape_elem_type::quads: return sample_quads(cdf, re, ruv);
        case shape_elem_type::beziers:
            throw std::runtime_error("type not supported");
        case shape_elem_type::vertices: return {sample_points(cdf, re), ruv};
        case shape_elem_type::facevarying: return sample_quads(cdf, re, ruv);
    }
}

// Subdivides shape elements.
void subdivide_shape_once(shape* shp, bool subdiv) {
    switch (get_shape_type(shp)) {
        case shape_elem_type::none: break;
        case shape_elem_type::points: break;
        case shape_elem_type::vertices: break;
        case shape_elem_type::triangles: {
            subdivide_triangles(shp->triangles, shp->pos, shp->norm,
                shp->texcoord, shp->color, shp->radius);
        } break;
        case shape_elem_type::lines: {
            subdivide_lines(shp->lines, shp->pos, shp->norm, shp->texcoord,
                shp->color, shp->radius);
        } break;
        case shape_elem_type::quads: {
            if (!subdiv) {
                subdivide_quads(shp->quads, shp->pos, shp->norm, shp->texcoord,
                    shp->color, shp->radius);
            } else {
                subdivide_catmullclark(shp->quads, shp->pos, shp->norm,
                    shp->texcoord, shp->color, shp->radius);
            }
        } break;
        case shape_elem_type::beziers: {
            subdivide_beziers(shp->beziers, shp->pos, shp->norm, shp->texcoord,
                shp->color, shp->radius);
        } break;
        case shape_elem_type::facevarying: {
            if (!subdiv) {
                subdivide_quads(shp->quads_pos, shp->pos);
                subdivide_quads(shp->quads_norm, shp->norm);
                subdivide_quads(shp->quads_texcoord, shp->texcoord);
            } else {
                subdivide_catmullclark(shp->quads_pos, shp->pos);
                subdivide_catmullclark(shp->quads_norm, shp->norm);
                subdivide_catmullclark(shp->quads_texcoord, shp->texcoord);
            }
        } break;
    }
}

// Compute shape normals. Supports only non-facevarying shapes.
void compute_normals(shape* shp) {
    switch (get_shape_type(shp)) {
        case shape_elem_type::none: break;
        case shape_elem_type::points: {
            shp->norm.assign(shp->pos.size(), {0, 0, 1});
        } break;
        case shape_elem_type::vertices: {
            shp->norm.assign(shp->pos.size(), {0, 0, 1});
        } break;
        case shape_elem_type::triangles: {
            compute_normals(shp->triangles, shp->pos, shp->norm);
        } break;
        case shape_elem_type::lines: {
            compute_tangents(shp->lines, shp->pos, shp->norm);
        } break;
        case shape_elem_type::quads: {
            compute_normals(shp->quads, shp->pos, shp->norm);
        } break;
        case shape_elem_type::beziers: {
            throw std::runtime_error("type not supported");
        } break;
        case shape_elem_type::facevarying: {
            throw std::runtime_error("type not supported");
        } break;
    }
}

// Tesselate a shape into basic primitives
void tesselate_shape(shape* shp, bool subdivide,
    bool facevarying_to_sharedvertex, bool quads_to_triangles,
    bool bezier_to_lines) {
    if (subdivide && shp->subdivision) {
        for (auto l = 0; l < shp->subdivision; l++) {
            subdivide_shape_once(shp, shp->catmullclark);
        }
    }
    auto type = get_shape_type(shp);
    if (facevarying_to_sharedvertex && type == shape_elem_type::facevarying) {
        std::tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
            convert_face_varying(shp->quads_pos, shp->quads_norm,
                shp->quads_texcoord, shp->pos, shp->norm, shp->texcoord);
        shp->quads_pos = {};
        shp->quads_norm = {};
        shp->quads_texcoord = {};
        type = get_shape_type(shp);
    }
    if (quads_to_triangles && type == shape_elem_type::quads) {
        shp->triangles = convert_quads_to_triangles(shp->quads);
        shp->quads = {};
        type = get_shape_type(shp);
    }
    if (bezier_to_lines && type == shape_elem_type::beziers) {
        shp->lines = convert_bezier_to_lines(shp->beziers);
        shp->beziers = {};
        type = get_shape_type(shp);
    }
}

// Tesselate scene shapes and update pointers
void tesselate_shapes(scene* scn, bool subdivide,
    bool facevarying_to_sharedvertex, bool quads_to_triangles,
    bool bezier_to_lines) {
    for (auto shp : scn->shapes) {
        tesselate_shape(shp, subdivide, facevarying_to_sharedvertex,
            quads_to_triangles, bezier_to_lines);
    }
}

// Update animation transforms
void update_transforms(
    const animation* anm, float time, const std::string& anim_group) {
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
void update_transforms(node* nde, const frame3f& parent = identity_frame3f) {
    auto frame = parent * nde->frame * translation_frame(nde->translation) *
                 rotation_frame(nde->rotation) * scaling_frame(nde->scale);
    if (nde->ist) nde->ist->frame = frame;
    if (nde->cam) nde->cam->frame = frame;
    if (nde->env) nde->env->frame = frame;
    for (auto child : nde->children_) update_transforms(child, frame);
}

// Update node transforms
void update_transforms(scene* scn, float time, const std::string& anim_group) {
    for (auto agr : scn->animations) update_transforms(agr, time, anim_group);
    for (auto nde : scn->nodes) nde->children_.clear();
    for (auto nde : scn->nodes)
        if (nde->parent) nde->parent->children_.push_back(nde);
    for (auto nde : scn->nodes)
        if (!nde->parent) update_transforms(nde);
}

// Compute animation range
vec2f compute_animation_range(const scene* scn, const std::string& anim_group) {
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

// Add missing names and resolve duplicated names.
void add_names(scene* scn) {
    auto fix_names = [](auto& vals, const std::string& base) {
        auto nmap = std::map<std::string, int>();
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
void add_normals(scene* scn) {
    for (auto shp : scn->shapes) {
        if (!shp->norm.empty()) continue;
        compute_normals(shp);
    }
}

// Add missing tangent space if needed.
void add_tangent_space(scene* scn) {
    for (auto ist : scn->instances) {
        if (!ist->shp->tangsp.empty() || ist->shp->texcoord.empty()) continue;
        if (!ist->mat || (!ist->mat->norm_txt.txt && !ist->mat->bump_txt.txt))
            continue;
        auto type = get_shape_type(ist->shp);
        if (type == shape_elem_type::triangles) {
            compute_tangent_frames(ist->shp->triangles, ist->shp->pos,
                ist->shp->norm, ist->shp->texcoord, ist->shp->tangsp);
        } else if (type == shape_elem_type::quads) {
            auto triangles = convert_quads_to_triangles(ist->shp->quads);
            compute_tangent_frames(triangles, ist->shp->pos, ist->shp->norm,
                ist->shp->texcoord, ist->shp->tangsp);
        } else {
            throw std::runtime_error("type not supported");
        }
    }
}

// Checks for validity of the scene.
std::vector<std::string> validate(
    const scene* scn, bool skip_missing, bool log_as_warning) {
    auto errs = std::vector<std::string>();
    auto check_names = [&errs](const auto& vals, const std::string& base) {
        auto used = std::map<std::string, int>();
        for (auto val : vals) used[val->name] += 1;
        for (auto& kv : used) {
            if (kv.first == "")
                errs.push_back("empty " + base + " name");
            else if (kv.second > 1)
                errs.push_back("duplicated " + base + " name " + kv.first);
        }
    };
    auto check_empty_textures = [&errs](const std::vector<texture*>& vals) {
        for (auto val : vals) {
            if (val->ldr.pixels.empty() && val->hdr.pixels.empty())
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
    if (!skip_missing) check_empty_textures(scn->textures);

    if (log_as_warning) {
        for (auto& err : errs) log_warning(err);
    }

    return errs;
}

// Make a view camera either copying a given one or building a
// default one.
camera* make_view_camera(const scene* scn, int camera_id) {
    if (scn->cameras.empty()) {
        auto bbox = compute_bbox(scn);
        auto bbox_center = (bbox.max + bbox.min) / 2.0f;
        auto bbox_size = bbox.max - bbox.min;
        auto bbox_msize = max(bbox_size.x, max(bbox_size.y, bbox_size.z));
        // set up camera
        auto cam = new camera();
        cam->name = "<view>";
        auto camera_dir = vec3f{1, 0.4f, 1};
        auto from = camera_dir * bbox_msize + bbox_center;
        auto to = bbox_center;
        auto up = vec3f{0, 1, 0};
        cam->frame = lookat_frame(from, to, up);
        cam->ortho = false;
        cam->aspect = 16.0f / 9.0f;
        cam->yfov = 2 * atanf(0.5f);
        cam->aperture = 0;
        cam->focus = length(to - from);
        return cam;
    } else {
        camera_id = clamp(camera_id, 0, (int)scn->cameras.size());
        auto cam = new camera(*scn->cameras[camera_id]);
        cam->name = "<view>";
        return cam;
    }
}

// Merge scene into one another
void merge_into(scene* merge_into, scene* merge_from) {
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

// Computes a shape bounding box (quick computation that ignores
// radius)
bbox3f compute_bbox(const shape* shp) {
    auto bbox = invalid_bbox3f;
    for (auto p : shp->pos) bbox += p;
    return bbox;
}

// Updates the scene and scene's instances bounding boxes
bbox3f compute_bbox(const scene* scn, bool skip_emitting) {
    auto shape_bboxes = std::unordered_map<shape*, bbox3f>();
    for (auto shp : scn->shapes) { shape_bboxes[shp] += compute_bbox(shp); }
    auto bbox = invalid_bbox3f;
    for (auto ist : scn->instances) {
        if (!ist->shp) continue;
        if (skip_emitting && ist->mat && ist->mat->ke != zero3f) continue;
        bbox += transform_bbox(ist->frame, shape_bboxes.at(ist->shp));
    }
    return bbox;
}

// Build a shape BVH
bvh_tree* make_bvh(const shape* shp, float def_radius, bool equalsize) {
    return make_bvh(shp->points, shp->lines, shp->triangles, shp->quads,
        shp->pos, shp->radius, def_radius, equalsize);
}

// Build a scene BVH
bvh_tree* make_bvh(const scene* scn, float def_radius, bool equalsize) {
    // do shapes
    auto shape_bvhs = std::vector<bvh_tree*>();
    auto smap = std::unordered_map<shape*, int>();
    for (auto shp : scn->shapes) {
        shape_bvhs.push_back(make_bvh(shp, def_radius, equalsize));
        smap[shp] = (int)shape_bvhs.size() - 1;
    }

    // tree bvh
    auto bists = std::vector<bvh_instance>();
    for (auto ist : scn->instances) {
        auto bist = bvh_instance();
        bist.frame = ist->frame;
        bist.frame_inv = inverse(ist->frame);
        bist.shape_id = smap.at(ist->shp);
        bists.push_back(bist);
    }
    return make_bvh(bists, shape_bvhs, equalsize, true);
}

// Refits a scene BVH
void refit_bvh(bvh_tree* bvh, const shape* shp, float def_radius) {
    refit_bvh(bvh, shp->pos, shp->radius, def_radius);
}

// Refits a scene BVH
void refit_bvh(
    bvh_tree* bvh, const scene* scn, bool do_shapes, float def_radius) {
    if (do_shapes) {
        auto sid = 0;
        for (auto shp : scn->shapes) {
            refit_bvh(get_shape_bvhs(bvh).at(sid++), shp->pos, shp->radius,
                def_radius);
        }
    }
    auto ist_frames = std::vector<frame3f>();
    auto ist_frames_inv = std::vector<frame3f>();
    for (auto ist : scn->instances) {
        ist_frames.push_back(ist->frame);
        ist_frames_inv.push_back(inverse(ist->frame));
    }
    refit_bvh(bvh, ist_frames, ist_frames_inv);
}

void print_stats(const scene* scn) {
    uint64_t num_cameras = 0;
    uint64_t num_shape_groups = 0;
    uint64_t num_shapes = 0;
    uint64_t num_instances = 0;
    uint64_t num_materials = 0;
    uint64_t num_textures = 0;
    uint64_t num_environments = 0;
    uint64_t num_nodes = 0;
    uint64_t num_animations = 0;

    uint64_t elem_points = 0;
    uint64_t elem_lines = 0;
    uint64_t elem_triangles = 0;
    uint64_t elem_quads = 0;
    uint64_t vert_pos = 0;
    uint64_t vert_norm = 0;
    uint64_t vert_texcoord = 0;
    uint64_t vert_color = 0;
    uint64_t vert_radius = 0;
    uint64_t vert_tangsp = 0;

    uint64_t texel_ldrs = 0;
    uint64_t texel_hdrs = 0;

    uint64_t memory_ldrs = 0;
    uint64_t memory_hdrs = 0;
    uint64_t memory_elems = 0;
    uint64_t memory_verts = 0;

    bbox3f bbox_scn = invalid_bbox3f;
    bbox3f bbox_nolights = invalid_bbox3f;

    num_cameras = scn->cameras.size();
    num_shapes = scn->shapes.size();
    num_materials = scn->materials.size();
    num_textures = scn->textures.size();
    num_environments = scn->environments.size();
    num_instances = scn->instances.size();
    num_nodes = scn->nodes.size();
    num_animations = scn->animations.size();

    for (auto shp : scn->shapes) {
        elem_points += shp->points.size();
        elem_lines += shp->lines.size();
        elem_triangles += shp->triangles.size();
        elem_quads += shp->quads.size();
        vert_pos += shp->pos.size();
        vert_norm += shp->norm.size();
        vert_texcoord += shp->texcoord.size();
        vert_color += shp->color.size();
        vert_radius += shp->radius.size();
        vert_tangsp += shp->tangsp.size();
    }
    memory_elems = elem_points * sizeof(int) + elem_lines * sizeof(vec2i) +
                   elem_triangles * sizeof(vec3i) + elem_quads * sizeof(vec4i);
    memory_verts = vert_pos * sizeof(vec3f) + vert_norm * sizeof(vec3f) +
                   vert_texcoord * sizeof(vec3f) + vert_color * sizeof(vec4f) +
                   vert_tangsp * sizeof(vec4f) + vert_radius * sizeof(float);

    for (auto txt : scn->textures) {
        texel_ldrs = txt->ldr.width * txt->ldr.height;
        texel_hdrs = txt->hdr.width * txt->hdr.height;
    }
    memory_ldrs = texel_ldrs * sizeof(vec4b);
    memory_hdrs = texel_hdrs * sizeof(vec4f);

    bbox_scn = compute_bbox(scn);
    bbox_nolights = compute_bbox(scn, true);

    println("num_cameras: {}", num_cameras);
    println("num_shape_groups: {}", num_shape_groups);
    println("num_shapes: {}", num_shapes);
    println("num_instances: {}", num_instances);
    println("num_materials: {}", num_materials);
    println("num_textures: {}", num_textures);
    println("num_environments: {}", num_environments);
    println("num_nodes: {}", num_nodes);
    println("num_animations: {}", num_animations);
    println("elem_points: {}", elem_points);
    println("elem_lines: {}", elem_lines);
    println("elem_triangles: {}", elem_triangles);
    println("elem_quads: {}", elem_quads);
    println("vert_pos: {}", vert_pos);
    println("vert_norm: {}", vert_norm);
    println("vert_texcoord: {}", vert_texcoord);
    println("vert_color: {}", vert_color);
    println("vert_radius: {}", vert_radius);
    println("vert_tangsp: {}", vert_tangsp);
    println("texel_ldrs: {}", texel_ldrs);
    println("texel_hdrs: {}", texel_hdrs);
    println("memory_ldrs: {}", memory_ldrs);
    println("memory_hdrs: {}", memory_hdrs);
    println("memory_elems: {}", memory_elems);
    println("memory_verts: {}", memory_verts);
    println("bbox_scn: {} {}", bbox_scn.min, bbox_scn.max);
    println("bbox_nolights: {} {}", bbox_nolights.min, bbox_nolights.max);
    println("bbox_min   : {}", bbox_scn.min);
    println("bbox_max   : {}", bbox_scn.max);
    println("bbox_size  : {}", bbox_scn.max - bbox_scn.min);
    println("bbox_center: {}", (bbox_scn.max + bbox_scn.min) / 2);
}

// Flattens an scene
scene* obj_to_scene(const obj_scene* obj, const load_options& opts) {
    // clear scene
    auto scn = new scene();

    struct obj_vertex_hash {
        std::hash<int> Th;
        size_t operator()(const obj_vertex& vv) const {
            auto v = (const int*)&vv;
            size_t h = 0;
            for (auto i = 0; i < sizeof(obj_vertex) / sizeof(int); i++) {
                // embads hash_combine below
                h ^= (Th(v[i]) + 0x9e3779b9 + (h << 6) + (h >> 2));
            }
            return h;
        }
    };

    // convert textures
    auto tmap = std::unordered_map<std::string, texture*>{{"", nullptr}};
    for (auto otxt : obj->textures) {
        auto txt = new texture();
        txt->name = otxt->path;
        txt->path = otxt->path;
        if (!otxt->datab.empty()) {
            txt->ldr = make_image4b(otxt->width, otxt->height, otxt->ncomp,
                otxt->datab.data(), vec4b{0, 0, 0, 255});
        } else if (!otxt->dataf.empty()) {
            txt->hdr = make_image4f(otxt->width, otxt->height, otxt->ncomp,
                otxt->dataf.data(), vec4f{0, 0, 0, 1});
        }
        scn->textures.push_back(txt);
        tmap[txt->path] = txt;
    }

    auto make_texture_info = [&tmap](const obj_texture_info& oinfo) {
        auto info = texture_info();
        if (oinfo.path == "") return info;
        info.txt = tmap.at(oinfo.path);
        info.wrap_s = !oinfo.clamp;
        info.wrap_t = !oinfo.clamp;
        info.scale = oinfo.scale;
        return info;
    };

    // convert materials and build textures
    auto mmap = std::unordered_map<std::string, material*>{{"", nullptr}};
    for (auto omat : obj->materials) {
        auto mat = new material();
        mat->name = omat->name;
        mat->type = material_type::specular_roughness;
        mat->ke = omat->ke;
        mat->kd = omat->kd;
        mat->ks = omat->ks;
        mat->kr = omat->kr;
        mat->kt = omat->kt;
        if (omat->ns >= 1e6f)
            mat->rs = 0;
        else if (omat->ns < 1)
            mat->rs = 1;
        else
            mat->rs = pow(2 / (omat->ns + 2), 1 / 4.0f);
        mat->op = omat->op;
        mat->ke_txt = make_texture_info(omat->ke_txt);
        mat->kd_txt = make_texture_info(omat->kd_txt);
        mat->ks_txt = make_texture_info(omat->ks_txt);
        mat->kr_txt = make_texture_info(omat->kr_txt);
        mat->kt_txt = make_texture_info(omat->kt_txt);
        mat->rs_txt = make_texture_info(omat->ns_txt);
        mat->op_txt = make_texture_info(omat->op_txt);
        mat->norm_txt = make_texture_info(omat->norm_txt);
        mat->bump_txt = make_texture_info(omat->bump_txt);
        mat->disp_txt = make_texture_info(omat->disp_txt);
        scn->materials.push_back(mat);
        mmap[mat->name] = mat;
    }

    // convert meshes
    auto omap = std::unordered_map<std::string,
        std::vector<std::pair<shape*, material*>>>{{"", {}}};
    for (auto omsh : obj->objects) {
        if (omsh->verts.empty()) continue;
        if (omsh->elems.empty()) continue;
        for (auto gid = 0; gid < omsh->groups.size(); gid++) {
            auto ogrp = omsh->groups.at(gid);
            auto shp = new shape();
            shp->name = omsh->name + ((gid) ? std::to_string(gid) : ""s);

            // check to see if this shuold be face-varying or flat
            // quads
            auto as_facevarying = false, as_quads = false;
            if (opts.obj_preserve_quads || opts.obj_preserve_facevarying) {
                auto face_max = 0;
                for (auto& elem : omsh->elems) {
                    if (elem.type != obj_element_type::face) {
                        face_max = 0;
                        break;
                    } else {
                        face_max = max(face_max, (int)elem.size);
                    }
                }
                as_quads = opts.obj_preserve_quads && face_max > 3;
                as_facevarying = opts.obj_preserve_facevarying && face_max > 2;
                // in case of facevarying, check if there is really
                // need for it
                if (as_facevarying) {
                    auto need_facevarying = false;
                    for (auto& elem : omsh->elems) {
                        for (auto i = elem.start; i < elem.start + elem.size;
                             i++) {
                            auto& v = omsh->verts[i];
                            if ((v.norm >= 0 && v.pos != v.norm) ||
                                (v.texcoord >= 0 && v.pos != v.texcoord) ||
                                (v.norm >= 0 && v.texcoord >= 0 &&
                                    v.norm != v.texcoord))
                                need_facevarying = true;
                            if (v.color >= 0 || v.radius >= 0)
                                as_facevarying = false;
                        }
                        if (!as_facevarying) break;
                    }
                    as_facevarying = need_facevarying;
                }
            }

            if (!as_facevarying) {
                // insert all vertices
                std::unordered_map<obj_vertex, int, obj_vertex_hash> vert_map;
                std::vector<int> vert_ids;
                for (auto& vert : omsh->verts) {
                    if (vert_map.find(vert) == vert_map.end()) {
                        auto s = (int)vert_map.size();
                        vert_map[vert] = s;
                    }
                    vert_ids.push_back(vert_map.at(vert));
                }

                // convert elements
                for (auto& elem : omsh->elems) {
                    if (elem.groupid != gid) continue;
                    switch (elem.type) {
                        case obj_element_type::point: {
                            for (auto i = elem.start;
                                 i < elem.start + elem.size; i++) {
                                shp->points.push_back(vert_ids[i]);
                            }
                        } break;
                        case obj_element_type::line: {
                            for (auto i = elem.start;
                                 i < elem.start + elem.size - 1; i++) {
                                shp->lines.push_back(
                                    {vert_ids[i], vert_ids[i + 1]});
                            }
                        } break;
                        case obj_element_type::face: {
                            if (as_quads && elem.size == 4) {
                                shp->quads.push_back({vert_ids[elem.start + 0],
                                    vert_ids[elem.start + 1],
                                    vert_ids[elem.start + 2],
                                    vert_ids[elem.start + 3]});
                            } else if (as_quads && elem.size != 4) {
                                for (auto i = elem.start + 2;
                                     i < elem.start + elem.size; i++) {
                                    shp->quads.push_back(
                                        {vert_ids[elem.start], vert_ids[i - 1],
                                            vert_ids[i], vert_ids[i]});
                                }
                            } else {
                                for (auto i = elem.start + 2;
                                     i < elem.start + elem.size; i++) {
                                    shp->triangles.push_back(
                                        {vert_ids[elem.start], vert_ids[i - 1],
                                            vert_ids[i]});
                                }
                            }
                        } break;
                        case obj_element_type::bezier: {
                            if ((elem.size - 1) % 3)
                                throw std::runtime_error("bad obj bezier");
                            for (auto i = elem.start + 1;
                                 i < elem.start + elem.size; i += 3) {
                                shp->beziers.push_back(
                                    {vert_ids[i - 1], vert_ids[i],
                                        vert_ids[i + 1], vert_ids[i + 2]});
                            }
                        } break;
                        default: { assert(false); }
                    }
                }

                // copy vertex data
                auto v = omsh->verts[0];
                if (v.pos >= 0) shp->pos.resize(vert_map.size());
                if (v.texcoord >= 0) shp->texcoord.resize(vert_map.size());
                if (v.norm >= 0) shp->norm.resize(vert_map.size());
                if (v.color >= 0) shp->color.resize(vert_map.size());
                if (v.radius >= 0) shp->radius.resize(vert_map.size());
                for (auto& kv : vert_map) {
                    auto idx = kv.second;
                    auto vert = kv.first;
                    if (v.pos >= 0 && vert.pos >= 0)
                        shp->pos[idx] = obj->pos[vert.pos];
                    if (v.texcoord >= 0 && vert.texcoord >= 0)
                        shp->texcoord[idx] = obj->texcoord[vert.texcoord];
                    if (v.norm >= 0 && vert.norm >= 0)
                        shp->norm[idx] = obj->norm[vert.norm];
                    if (v.color >= 0 && vert.color >= 0)
                        shp->color[idx] = obj->color[vert.color];
                    if (v.radius >= 0 && vert.radius >= 0)
                        shp->radius[idx] = obj->radius[vert.radius];
                }
            } else {
                // insert all vertices
                std::unordered_map<int, int> pos_map, norm_map, texcoord_map;
                std::vector<int> pos_ids, norm_ids, texcoord_ids;
                for (auto& vert : omsh->verts) {
                    if (vert.pos >= 0) {
                        if (pos_map.find(vert.pos) == pos_map.end()) {
                            auto s = (int)pos_map.size();
                            pos_map[vert.pos] = s;
                        }
                        pos_ids.push_back(pos_map.at(vert.pos));
                    } else {
                        if (!pos_ids.empty())
                            throw std::runtime_error("malformed obj");
                    }
                    if (vert.norm >= 0) {
                        if (norm_map.find(vert.norm) == norm_map.end()) {
                            auto s = (int)norm_map.size();
                            norm_map[vert.norm] = s;
                        }
                        norm_ids.push_back(norm_map.at(vert.norm));
                    } else {
                        if (!norm_ids.empty())
                            throw std::runtime_error("malformed obj");
                    }
                    if (vert.texcoord >= 0) {
                        if (texcoord_map.find(vert.texcoord) ==
                            texcoord_map.end()) {
                            auto s = (int)texcoord_map.size();
                            texcoord_map[vert.texcoord] = s;
                        }
                        texcoord_ids.push_back(texcoord_map.at(vert.texcoord));
                    } else {
                        if (!texcoord_ids.empty())
                            throw std::runtime_error("malformed obj");
                    }
                }

                // convert elements
                for (auto elem : omsh->elems) {
                    if (elem.groupid != gid) continue;
                    if (elem.size == 4) {
                        if (!pos_ids.empty()) {
                            shp->quads_pos.push_back({pos_ids[elem.start + 0],
                                pos_ids[elem.start + 1],
                                pos_ids[elem.start + 2],
                                pos_ids[elem.start + 3]});
                        }
                        if (!texcoord_ids.empty()) {
                            shp->quads_texcoord.push_back(
                                {texcoord_ids[elem.start + 0],
                                    texcoord_ids[elem.start + 1],
                                    texcoord_ids[elem.start + 2],
                                    texcoord_ids[elem.start + 3]});
                        }
                        if (!norm_ids.empty()) {
                            shp->quads_norm.push_back({norm_ids[elem.start + 0],
                                norm_ids[elem.start + 1],
                                norm_ids[elem.start + 2],
                                norm_ids[elem.start + 3]});
                        }
                    } else {
                        if (!pos_ids.empty()) {
                            for (auto i = elem.start + 2;
                                 i < elem.start + elem.size; i++) {
                                shp->quads_pos.push_back({pos_ids[elem.start],
                                    pos_ids[i - 1], pos_ids[i], pos_ids[i]});
                            }
                        }
                        if (!texcoord_ids.empty()) {
                            for (auto i = elem.start + 2;
                                 i < elem.start + elem.size; i++) {
                                shp->quads_texcoord.push_back(
                                    {texcoord_ids[elem.start],
                                        texcoord_ids[i - 1], texcoord_ids[i],
                                        texcoord_ids[i]});
                            }
                        }
                        if (!norm_ids.empty()) {
                            for (auto i = elem.start + 2;
                                 i < elem.start + elem.size; i++) {
                                shp->quads_norm.push_back({norm_ids[elem.start],
                                    norm_ids[i - 1], norm_ids[i], norm_ids[i]});
                            }
                        }
                    }
                }

                // copy vertex data
                shp->pos.resize(pos_map.size());
                shp->texcoord.resize(texcoord_map.size());
                shp->norm.resize(norm_map.size());
                for (auto& kv : pos_map) {
                    shp->pos[kv.second] = obj->pos[kv.first];
                }
                for (auto& kv : texcoord_map) {
                    shp->texcoord[kv.second] = obj->texcoord[kv.first];
                }
                for (auto& kv : norm_map) {
                    shp->norm[kv.second] = obj->norm[kv.first];
                }
            }
            scn->shapes.push_back(shp);
            omap[omsh->name].push_back({shp, mmap[ogrp.matname]});
        }
    }

    // convert cameras
    auto cmap = std::unordered_map<std::string, camera*>{{"", nullptr}};
    for (auto ocam : obj->cameras) {
        auto cam = new camera();
        cam->name = ocam->name;
        cam->ortho = ocam->ortho;
        cam->yfov = ocam->yfov;
        cam->aspect = ocam->aspect;
        cam->aperture = ocam->aperture;
        cam->focus = ocam->focus;
        cam->frame = ocam->frame;
        scn->cameras.push_back(cam);
        cmap[cam->name] = cam;
    }

    // convert envs
    std::unordered_set<material*> env_mat;
    auto emap = std::unordered_map<std::string, environment*>{{"", nullptr}};
    for (auto oenv : obj->environments) {
        auto env = new environment();
        env->name = oenv->name;
        env->ke = oenv->ke;
        env->ke_txt = make_texture_info(oenv->ke_txt);
        env->frame = oenv->frame;
        scn->environments.push_back(env);
        emap[env->name] = env;
    }

    // convert nodes
    if (!obj->nodes.empty()) {
        for (auto onde : obj->nodes) {
            auto nde = new node();
            nde->name = onde->name;
            nde->cam = cmap.at(onde->camname);
            nde->env = emap.at(onde->envname);
            nde->translation = onde->translation;
            nde->rotation = onde->rotation;
            nde->scale = onde->scale;
            nde->frame = onde->frame;
            scn->nodes.push_back(nde);
        }

        // set up parent pointers
        for (auto nid = 0; nid < obj->nodes.size(); nid++) {
            auto onde = obj->nodes[nid];
            if (onde->parent.empty()) continue;
            auto nde = scn->nodes[nid];
            for (auto parent : scn->nodes) {
                if (parent->name == onde->parent) {
                    nde->parent = parent;
                    break;
                }
            }
        }

        // set up instances
        for (auto nid = 0; nid < obj->nodes.size(); nid++) {
            auto onde = obj->nodes[nid];
            if (onde->objname.empty()) continue;
            auto nde = scn->nodes[nid];
            auto& shps = omap.at(onde->objname);
            if (shps.empty()) continue;
            if (shps.size() == 1) {
                nde->ist = new instance();
                nde->ist->name = nde->name;
                nde->ist->shp = shps[0].first;
                nde->ist->mat = shps[0].second;
                scn->instances.push_back(nde->ist);
            } else {
                for (auto shp : shps) {
                    auto child = new node();
                    child->name = nde->name + "_" + shp.first->name;
                    child->parent = nde;
                    child->ist = new instance();
                    child->ist->name = child->name;
                    child->ist->shp = shp.first;
                    child->ist->mat = shp.second;
                    scn->instances.push_back(child->ist);
                }
            }
        }
    } else {
        for (auto& shps : omap) {
            for (auto shp : shps.second) {
                auto ist = new instance();
                ist->name = shp.first->name;
                ist->shp = shp.first;
                ist->mat = shp.second;
                scn->instances.push_back(ist);
            }
        }
    }

    // update transforms
    update_transforms(scn);

    // done
    return scn;
}

// Load an obj scene
scene* load_obj_scene(const std::string& filename, const load_options& opts) {
    auto oscn = load_obj(filename, opts.obj_split_shapes, opts.load_textures,
        opts.skip_missing, opts.obj_flip_texcoord, opts.obj_flip_tr);
    auto scn = obj_to_scene(oscn, opts);
    delete oscn;
    return scn;
}

// Save an scene
obj_scene* scene_to_obj(const scene* scn, const save_options& opts) {
    auto obj = new obj_scene();

    auto make_texture_info = [](const texture_info& info, bool bump = false) {
        auto oinfo = obj_texture_info();
        if (!info.txt) return oinfo;
        oinfo.path = info.txt->path;
        oinfo.clamp = !info.wrap_s && !info.wrap_t;
        if (bump) oinfo.scale = info.scale;
        return oinfo;
    };

    // convert textures
    for (auto txt : scn->textures) {
        auto otxt = new obj_texture();
        otxt->path = txt->path;
        if (!txt->hdr.pixels.empty()) {
            otxt->width = txt->hdr.width;
            otxt->height = txt->hdr.height;
            otxt->ncomp = 4;
            otxt->dataf.assign((float*)txt->hdr.pixels.data(),
                (float*)txt->hdr.pixels.data() +
                    txt->hdr.width * txt->hdr.height * 4);
        }
        if (!txt->ldr.pixels.empty()) {
            otxt->width = txt->ldr.width;
            otxt->height = txt->ldr.height;
            otxt->ncomp = 4;
            otxt->datab.assign((uint8_t*)txt->ldr.pixels.data(),
                (uint8_t*)txt->ldr.pixels.data() +
                    txt->ldr.width * txt->ldr.height * 4);
        }
        obj->textures.push_back(otxt);
    }

    // convert materials
    for (auto mat : scn->materials) {
        auto omat = new obj_material();
        omat->name = mat->name;
        omat->ke = {mat->ke.x, mat->ke.y, mat->ke.z};
        omat->ke_txt = make_texture_info(mat->ke_txt);
        switch (mat->type) {
            case material_type::specular_roughness: {
                omat->kd = {mat->kd.x, mat->kd.y, mat->kd.z};
                omat->ks = {mat->ks.x, mat->ks.y, mat->ks.z};
                omat->kr = {mat->kr.x, mat->kr.y, mat->kr.z};
                omat->kt = {mat->kt.x, mat->kt.y, mat->kt.z};
                if (mat->rs <= 0)
                    omat->ns = 1e6f;
                else if (mat->rs >= 1)
                    omat->ns = 0;
                else
                    omat->ns = 2 / pow(mat->rs, 4.0f) - 2;
                omat->op = mat->op;
                omat->kd_txt = make_texture_info(mat->kd_txt);
                omat->ks_txt = make_texture_info(mat->ks_txt);
                omat->kr_txt = make_texture_info(mat->kr_txt);
                omat->kt_txt = make_texture_info(mat->kt_txt);
                omat->op_txt = make_texture_info(mat->op_txt);
            } break;
            case material_type::metallic_roughness: {
                if (mat->rs >= 1 && mat->ks.x == 0) {
                    omat->kd = mat->kd;
                    omat->ks = {0, 0, 0};
                    omat->ns = 1;
                } else {
                    auto kd = mat->kd * (1 - 0.04f) * (1 - mat->ks.x);
                    auto ks = mat->kd * mat->ks.x +
                              vec3f{0.04f, 0.04f, 0.04f} * (1 - mat->ks.x);
                    omat->kd = {kd.x, kd.y, kd.z};
                    omat->ks = {ks.x, ks.y, ks.z};
                    if (mat->rs <= 0)
                        omat->ns = 1e6f;
                    else if (mat->rs >= 1)
                        omat->ns = 0;
                    else
                        omat->ns = 2 / pow(mat->rs, 4.0f) - 2;
                }
                omat->op = mat->op;
                if (mat->ks.x < 0.5f) {
                    omat->kd_txt = make_texture_info(mat->kd_txt);
                } else {
                    omat->ks_txt = make_texture_info(mat->ks_txt);
                }
            } break;
            case material_type::specular_glossiness: {
                omat->kd = {mat->kd.x, mat->kd.y, mat->kd.z};
                omat->ks = {mat->ks.x, mat->ks.y, mat->ks.z};
                if (mat->rs <= 0)
                    omat->ns = 1e6f;
                else if (mat->rs >= 1)
                    omat->ns = 0;
                else
                    omat->ns = 2 / pow(mat->rs, 4.0f) - 2;
                omat->op = mat->op;
                omat->kd_txt = make_texture_info(mat->kd_txt);
                omat->ks_txt = make_texture_info(mat->ks_txt);
            } break;
        }
        omat->bump_txt = make_texture_info(mat->bump_txt, true);
        omat->disp_txt = make_texture_info(mat->disp_txt, true);
        omat->norm_txt = make_texture_info(mat->norm_txt, true);
        if (mat->op < 1 || mat->kt != zero3f) {
            omat->illum = 4;
        } else {
            omat->illum = 2;
        }
        obj->materials.push_back(omat);
    }

    // add elem
    auto add_elem = [](auto& shp, auto& oobj, auto etype, auto esize,
                        auto eid) {
        auto elem = obj_element();
        elem.start = (uint32_t)oobj->verts.size();
        elem.type = etype;
        elem.size = (uint16_t)esize;
        elem.groupid = 0;
        oobj->elems.push_back(elem);
    };
    // add vertex
    auto add_vert = [](auto& shp, auto& oobj, auto& offset, int vid) {
        auto vert = obj_vertex{-1, -1, -1, -1, -1};
        if (!shp->pos.empty()) vert.pos = offset.pos + vid;
        if (!shp->texcoord.empty()) vert.texcoord = offset.texcoord + vid;
        if (!shp->norm.empty()) vert.norm = offset.norm + vid;
        if (!shp->color.empty()) vert.color = offset.color + vid;
        if (!shp->radius.empty()) vert.radius = offset.radius + vid;
        oobj->verts.push_back(vert);
    };

    // flatten instances if necessary
    auto flatten_instances = !opts.obj_save_instances && scn->nodes.empty() &&
                             !scn->instances.empty();
    auto shapes = std::vector<std::pair<shape*, frame3f>>();
    if (flatten_instances) {
        for (auto ist : scn->instances)
            shapes.push_back({ist->shp, ist->frame});
    } else {
        for (auto shp : scn->shapes) shapes.push_back({shp, identity_frame3f});
    }
    auto shape_mats = std::map<shape*, material*>();
    for (auto shp : scn->shapes) shape_mats[shp] = nullptr;
    for (auto ist : scn->instances) {
        if (shape_mats.at(ist->shp) && shape_mats.at(ist->shp) != ist->mat)
            log_error("shapes can only have one material associated");
        else
            shape_mats[ist->shp] = ist->mat;
    }

    // convert shapes
    for (auto& shp_frame : shapes) {
        auto shp = shp_frame.first;
        auto frame = shp_frame.second;
        auto mat = shape_mats.at(shp);
        auto oobj = new obj_object();
        oobj->name = shp->name;
        if (shp->subdivision)
            oobj->props["subdivision"].push_back(
                std::to_string(shp->subdivision));
        if (shp->catmullclark) oobj->props["catmullclark"].push_back("1");
        auto offset = obj_vertex{(int)obj->pos.size(),
            (int)obj->texcoord.size(), (int)obj->norm.size(),
            (int)obj->color.size(), (int)obj->radius.size()};
        if (frame != identity_frame3f) {
            for (auto& v : shp->pos)
                obj->pos.push_back(transform_point(frame, v));
            for (auto& v : shp->norm)
                obj->norm.push_back(transform_direction(frame, v));
        } else {
            for (auto& v : shp->pos) obj->pos.push_back(v);
            for (auto& v : shp->norm) obj->norm.push_back(v);
        }
        for (auto& v : shp->texcoord) obj->texcoord.push_back(v);
        for (auto& v : shp->color) obj->color.push_back({v.x, v.y, v.z, v.w});
        for (auto& v : shp->radius) obj->radius.push_back(v);
        oobj->groups.push_back(
            {"", (mat) ? mat->name : ""s, shp->norm.empty()});
        for (auto eid = 0; eid < shp->points.size(); eid++) {
            auto vid = shp->points[eid];
            add_elem(shp, oobj, obj_element_type::point, 1, eid);
            add_vert(shp, oobj, offset, vid);
        }
        for (auto eid = 0; eid < shp->lines.size(); eid++) {
            auto l = shp->lines[eid];
            add_elem(shp, oobj, obj_element_type::line, 2, eid);
            for (auto vid : {l.x, l.y}) add_vert(shp, oobj, offset, vid);
        }
        for (auto eid = 0; eid < shp->triangles.size(); eid++) {
            auto t = shp->triangles[eid];
            add_elem(shp, oobj, obj_element_type::face, 3, eid);
            for (auto vid : {t.x, t.y, t.z}) add_vert(shp, oobj, offset, vid);
        }
        for (auto eid = 0; eid < shp->quads.size(); eid++) {
            auto q = shp->quads[eid];
            add_elem(shp, oobj, obj_element_type::face,
                (uint16_t)((q.z == q.w) ? 3 : 4), eid);
            if (oobj->elems.back().size == 3) {
                for (auto vid : {q.x, q.y, q.z})
                    add_vert(shp, oobj, offset, vid);
            } else {
                for (auto vid : {q.x, q.y, q.z, q.w})
                    add_vert(shp, oobj, offset, vid);
            }
        }
        for (auto eid = 0; eid < shp->beziers.size(); eid++) {
            auto b = shp->beziers[eid];
            add_elem(shp, oobj, obj_element_type::bezier, 4, eid);
            for (auto vid : {b.x, b.y, b.z, b.w})
                add_vert(shp, oobj, offset, vid);
        }
        for (auto eid = 0; eid < shp->quads_pos.size(); eid++) {
            add_elem(shp, oobj, obj_element_type::face, 4, eid);
            auto last_vid = -1;
            for (auto i = 0; i < 4; i++) {
                if (last_vid == (&shp->quads_pos[eid].x)[i]) continue;
                auto vert = obj_vertex{-1, -1, -1, -1, -1};
                if (!shp->pos.empty() && !shp->quads_pos.empty())
                    vert.pos = offset.pos + (&shp->quads_pos[eid].x)[i];
                if (!shp->texcoord.empty() && !shp->quads_texcoord.empty())
                    vert.texcoord =
                        offset.texcoord + (&shp->quads_texcoord[eid].x)[i];
                if (!shp->norm.empty() && !shp->quads_norm.empty())
                    vert.norm = offset.norm + (&shp->quads_norm[eid].x)[i];
                oobj->verts.push_back(vert);
                last_vid = (&shp->quads_pos[eid].x)[i];
            }
        }
        obj->objects.push_back(oobj);
    }

    // convert cameras
    for (auto cam : scn->cameras) {
        auto ocam = new obj_camera();
        ocam->name = cam->name;
        ocam->ortho = cam->ortho;
        ocam->yfov = cam->yfov;
        ocam->aspect = cam->aspect;
        ocam->focus = cam->focus;
        ocam->aperture = cam->aperture;
        ocam->frame = cam->frame;
        obj->cameras.push_back(ocam);
    }

    // convert envs
    for (auto env : scn->environments) {
        auto oenv = new obj_environment();
        oenv->name = env->name;
        oenv->ke = env->ke;
        oenv->ke_txt = make_texture_info(env->ke_txt);
        oenv->frame = env->frame;
        obj->environments.push_back(oenv);
    }

    // convert hierarchy
    if (!scn->nodes.empty()) {
        for (auto nde : scn->nodes) {
            auto onde = new obj_node();
            onde->name = nde->name;
            if (nde->cam) onde->camname = nde->cam->name;
            if (nde->ist) onde->objname = nde->ist->shp->name;
            if (nde->env) onde->envname = nde->env->name;
            onde->frame = nde->frame;
            onde->translation = nde->translation;
            onde->rotation = nde->rotation;
            onde->scale = nde->scale;
            obj->nodes.push_back(onde);
        }

        // parent
        for (auto idx = 0; idx < scn->nodes.size(); idx++) {
            auto nde = scn->nodes.at(idx);
            if (!nde->parent) continue;
            auto onde = obj->nodes.at(idx);
            onde->parent = nde->parent->name;
        }
    } else if (!flatten_instances) {
        for (auto ist : scn->instances) {
            auto onde = new obj_node();
            onde->name = ist->name;
            onde->objname = ist->shp->name;
            onde->frame = ist->frame;
            obj->nodes.push_back(onde);
        }
    }

    return obj;
}

// Save an obj scene
void save_obj_scene(
    const std::string& filename, const scene* scn, const save_options& opts) {
    auto oscn = scene_to_obj(scn, opts);
    save_obj(filename, oscn, opts.save_textures, opts.skip_missing,
        opts.obj_flip_texcoord, opts.obj_flip_tr);
    delete oscn;
}

#if YGL_GLTF

// Flattens a gltf file into a flattened asset.
scene* gltf_to_scene(const glTF* gltf, const load_options& opts) {
    auto startswith = [](const std::string& str, const std::string& substr) {
        if (str.length() < substr.length()) return false;
        for (auto i = 0; i < substr.length(); i++)
            if (str[i] != substr[i]) return false;
        return true;
    };

    // clear asset
    auto scn = new scene();

    // convert images
    for (auto gtxt : gltf->images) {
        auto txt = new texture();
        txt->name = gtxt->name;
        txt->path = (startswith(gtxt->uri, "data:")) ? std::string("inlines") :
                                                       gtxt->uri;
        if (!gtxt->data.datab.empty()) {
            txt->ldr = make_image4b(gtxt->data.width, gtxt->data.height,
                gtxt->data.ncomp, gtxt->data.datab.data(), vec4b{0, 0, 0, 255});
        } else if (!gtxt->data.dataf.empty()) {
            txt->hdr = make_image4f(gtxt->data.width, gtxt->data.height,
                gtxt->data.ncomp, gtxt->data.dataf.data(), vec4f{0, 0, 0, 1});
        }
        scn->textures.push_back(txt);
    }

    // add a texture
    auto make_texture_info = [gltf, scn](const glTFTextureInfo* ginfo,
                                 bool normal = false, bool occlusion = false) {
        auto info = texture_info();
        if (!ginfo) return info;
        auto gtxt = gltf->get(ginfo->index);
        if (!gtxt || !gtxt->source) return info;
        info.txt = scn->textures.at((int)gtxt->source);
        if (!info.txt) return info;
        auto gsmp = gltf->get(gtxt->sampler);
        if (gsmp) {
            info.linear = gsmp->magFilter != glTFSamplerMagFilter::Nearest;
            info.mipmap = gsmp->minFilter != glTFSamplerMinFilter::Linear &&
                          gsmp->minFilter != glTFSamplerMinFilter::Nearest;
            info.wrap_s = gsmp->wrapS != glTFSamplerWrapS::ClampToEdge;
            info.wrap_t = gsmp->wrapT != glTFSamplerWrapT::ClampToEdge;
        }
        if (normal) {
            auto ninfo = (glTFMaterialNormalTextureInfo*)ginfo;
            info.scale = ninfo->scale;
        }
        if (occlusion) {
            auto ninfo = (glTFMaterialOcclusionTextureInfo*)ginfo;
            info.scale = ninfo->strength;
        }
        return info;
    };

    // convert materials
    for (auto gmat : gltf->materials) {
        auto mat = new material();
        mat->name = gmat->name;
        mat->ke = gmat->emissiveFactor;
        mat->ke_txt = make_texture_info(gmat->emissiveTexture);
        if (gmat->pbrMetallicRoughness) {
            mat->type = material_type::metallic_roughness;
            auto gmr = gmat->pbrMetallicRoughness;
            mat->kd = {gmr->baseColorFactor.x, gmr->baseColorFactor.y,
                gmr->baseColorFactor.z};
            mat->op = gmr->baseColorFactor.w;
            mat->ks = {
                gmr->metallicFactor, gmr->metallicFactor, gmr->metallicFactor};
            mat->rs = gmr->roughnessFactor;
            mat->kd_txt = make_texture_info(gmr->baseColorTexture);
            mat->ks_txt = make_texture_info(gmr->metallicRoughnessTexture);
        }
        if (gmat->pbrSpecularGlossiness) {
            mat->type = material_type::specular_glossiness;
            auto gsg = gmat->pbrSpecularGlossiness;
            mat->kd = {gsg->diffuseFactor.x, gsg->diffuseFactor.y,
                gsg->diffuseFactor.z};
            mat->op = gsg->diffuseFactor.w;
            mat->ks = gsg->specularFactor;
            mat->rs = gsg->glossinessFactor;
            mat->kd_txt = make_texture_info(gsg->diffuseTexture);
            mat->ks_txt = make_texture_info(gsg->specularGlossinessTexture);
        }
        mat->norm_txt = make_texture_info(gmat->normalTexture, true, false);
        mat->occ_txt = make_texture_info(gmat->occlusionTexture, false, true);
        mat->double_sided = gmat->doubleSided;
        scn->materials.push_back(mat);
    }

    // convert meshes
    auto meshes = std::vector<std::vector<std::pair<shape*, material*>>>();
    for (auto gmesh : gltf->meshes) {
        meshes.push_back({});
        auto sid = 0;
        for (auto gprim : gmesh->primitives) {
            auto shp = new shape();
            shp->name = gmesh->name + ((sid) ? std::to_string(sid) : ""s);
            sid++;
            // vertex data
            for (auto gattr : gprim->attributes) {
                auto semantic = gattr.first;
                auto vals = accessor_view(gltf, gltf->get(gattr.second));
                if (semantic == "POSITION") {
                    shp->pos.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        shp->pos.push_back(vals.getv3f(i));
                } else if (semantic == "NORMAL") {
                    shp->norm.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        shp->norm.push_back(vals.getv3f(i));
                } else if (semantic == "TEXCOORD" || semantic == "TEXCOORD_0") {
                    shp->texcoord.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        shp->texcoord.push_back(vals.getv2f(i));
                } else if (semantic == "TEXCOORD_1") {
                    shp->texcoord1.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        shp->texcoord1.push_back(vals.getv2f(i));
                } else if (semantic == "COLOR" || semantic == "COLOR_0") {
                    shp->color.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        shp->color.push_back(vals.getv4f(i, {0, 0, 0, 1}));
                } else if (semantic == "TANGENT") {
                    shp->tangsp.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        shp->tangsp.push_back(vals.getv4f(i));
                } else if (semantic == "RADIUS") {
                    shp->radius.reserve(vals.size());
                    for (auto i = 0; i < vals.size(); i++)
                        shp->radius.push_back(vals.get(i, 0));
                } else {
                    // ignore
                }
            }
            // indices
            if (!gprim->indices) {
                switch (gprim->mode) {
                    case glTFMeshPrimitiveMode::Triangles: {
                        shp->triangles.reserve(shp->pos.size() / 3);
                        for (auto i = 0; i < shp->pos.size() / 3; i++) {
                            shp->triangles.push_back(
                                {i * 3 + 0, i * 3 + 1, i * 3 + 2});
                        }
                    } break;
                    case glTFMeshPrimitiveMode::TriangleFan: {
                        shp->triangles.reserve(shp->pos.size() - 2);
                        for (auto i = 2; i < shp->pos.size(); i++) {
                            shp->triangles.push_back({0, i - 1, i});
                        }
                    } break;
                    case glTFMeshPrimitiveMode::TriangleStrip: {
                        shp->triangles.reserve(shp->pos.size() - 2);
                        for (auto i = 2; i < shp->pos.size(); i++) {
                            shp->triangles.push_back({i - 2, i - 1, i});
                        }
                    } break;
                    case glTFMeshPrimitiveMode::Lines: {
                        shp->lines.reserve(shp->pos.size() / 2);
                        for (auto i = 0; i < shp->pos.size() / 2; i++) {
                            shp->lines.push_back({i * 2 + 0, i * 2 + 1});
                        }
                    } break;
                    case glTFMeshPrimitiveMode::LineLoop: {
                        shp->lines.reserve(shp->pos.size());
                        for (auto i = 1; i < shp->pos.size(); i++) {
                            shp->lines.push_back({i - 1, i});
                        }
                        shp->lines.back() = {(int)shp->pos.size() - 1, 0};
                    } break;
                    case glTFMeshPrimitiveMode::LineStrip: {
                        shp->lines.reserve(shp->pos.size() - 1);
                        for (auto i = 1; i < shp->pos.size(); i++) {
                            shp->lines.push_back({i - 1, i});
                        }
                    } break;
                    case glTFMeshPrimitiveMode::NotSet:
                    case glTFMeshPrimitiveMode::Points: {
                        shp->points.reserve(shp->pos.size());
                        for (auto i = 0; i < shp->pos.size(); i++) {
                            shp->points.push_back(i);
                        }
                    } break;
                }
            } else {
                auto indices = accessor_view(gltf, gltf->get(gprim->indices));
                switch (gprim->mode) {
                    case glTFMeshPrimitiveMode::Triangles: {
                        shp->triangles.reserve(indices.size());
                        for (auto i = 0; i < indices.size() / 3; i++) {
                            shp->triangles.push_back({indices.geti(i * 3 + 0),
                                indices.geti(i * 3 + 1),
                                indices.geti(i * 3 + 2)});
                        }
                    } break;
                    case glTFMeshPrimitiveMode::TriangleFan: {
                        shp->triangles.reserve(indices.size() - 2);
                        for (auto i = 2; i < indices.size(); i++) {
                            shp->triangles.push_back({indices.geti(0),
                                indices.geti(i - 1), indices.geti(i)});
                        }
                    } break;
                    case glTFMeshPrimitiveMode::TriangleStrip: {
                        shp->triangles.reserve(indices.size() - 2);
                        for (auto i = 2; i < indices.size(); i++) {
                            shp->triangles.push_back({indices.geti(i - 2),
                                indices.geti(i - 1), indices.geti(i)});
                        }
                    } break;
                    case glTFMeshPrimitiveMode::Lines: {
                        shp->lines.reserve(indices.size() / 2);
                        for (auto i = 0; i < indices.size() / 2; i++) {
                            shp->lines.push_back({indices.geti(i * 2 + 0),
                                indices.geti(i * 2 + 1)});
                        }
                    } break;
                    case glTFMeshPrimitiveMode::LineLoop: {
                        shp->lines.reserve(indices.size());
                        for (auto i = 1; i < indices.size(); i++) {
                            shp->lines.push_back(
                                {indices.geti(i - 1), indices.geti(i)});
                        }
                        shp->lines.back() = {
                            indices.geti(indices.size() - 1), indices.geti(0)};
                    } break;
                    case glTFMeshPrimitiveMode::LineStrip: {
                        shp->lines.reserve(indices.size() - 1);
                        for (auto i = 1; i < indices.size(); i++) {
                            shp->lines.push_back(
                                {indices.geti(i - 1), indices.geti(i)});
                        }
                    } break;
                    case glTFMeshPrimitiveMode::NotSet:
                    case glTFMeshPrimitiveMode::Points: {
                        shp->points.reserve(indices.size());
                        for (auto i = 0; i < indices.size(); i++) {
                            shp->points.push_back(indices.geti(i));
                        }
                    } break;
                }
            }
            meshes.back().push_back(
                {shp, scn->materials[(int)gprim->material]});
            scn->shapes.push_back(shp);
        }
    }

    // convert cameras
    for (auto gcam : gltf->cameras) {
        auto cam = new camera();
        cam->name = gcam->name;
        cam->ortho = gcam->type == glTFCameraType::Orthographic;
        if (cam->ortho) {
            auto ortho = gcam->orthographic;
            cam->yfov = ortho->ymag;
            cam->aspect = ortho->xmag / ortho->ymag;
            cam->near = ortho->znear;
            cam->far = ortho->zfar;
        } else {
            auto persp = gcam->perspective;
            cam->yfov = persp->yfov;
            cam->aspect = persp->aspectRatio;
            if (!cam->aspect) cam->aspect = 16.0f / 9.0f;
            cam->near = persp->znear;
            cam->far = persp->zfar;
        }
        scn->cameras.push_back(cam);
    }

    // convert nodes
    for (auto gnde : gltf->nodes) {
        auto nde = new node();
        nde->name = gnde->name;
        if (gnde->camera) nde->cam = scn->cameras[(int)gnde->camera];
        nde->translation = gnde->translation;
        nde->rotation = gnde->rotation;
        nde->scale = gnde->scale;
        nde->frame = mat_to_frame(gnde->matrix);
        scn->nodes.push_back(nde);
    }

    // set up parent pointers
    for (auto nid = 0; nid < gltf->nodes.size(); nid++) {
        auto gnde = gltf->nodes[nid];
        auto nde = scn->nodes[nid];
        for (auto cid : gnde->children) scn->nodes[(int)cid]->parent = nde;
    }

    // set up instances
    for (auto nid = 0; nid < gltf->nodes.size(); nid++) {
        auto gnde = gltf->nodes[nid];
        if (!gnde->mesh) continue;
        auto nde = scn->nodes[nid];
        auto& shps = meshes.at((int)gnde->mesh);
        if (shps.empty()) continue;
        if (shps.size() == 1) {
            nde->ist = new instance();
            nde->ist->name = nde->name;
            nde->ist->shp = shps[0].first;
            nde->ist->mat = shps[0].second;
            scn->instances.push_back(nde->ist);
        } else {
            for (auto shp : shps) {
                auto child = new node();
                child->name = nde->name + "_" + shp.first->name;
                child->parent = nde;
                child->ist = new instance();
                child->ist->name = child->name;
                child->ist->shp = shp.first;
                child->ist->mat = shp.second;
                scn->instances.push_back(child->ist);
            }
        }
    }

    // keyframe type conversion
    static auto keyframe_types =
        std::unordered_map<glTFAnimationSamplerInterpolation, animation_type>{
            {glTFAnimationSamplerInterpolation::NotSet, animation_type::linear},
            {glTFAnimationSamplerInterpolation::Linear, animation_type::linear},
            {glTFAnimationSamplerInterpolation::Step, animation_type::step},
            {glTFAnimationSamplerInterpolation::CubicSpline,
                animation_type::bezier},
        };

    // convert animations
    for (auto ganm : gltf->animations) {
        auto aid = 0;
        auto sampler_map = std::unordered_map<vec2i, int>();
        for (auto gchannel : ganm->channels) {
            if (sampler_map.find({(int)gchannel->sampler,
                    (int)gchannel->target->path}) == sampler_map.end()) {
                auto gsampler = ganm->get(gchannel->sampler);
                auto anm = new animation();
                anm->name = ((ganm->name != "") ? ganm->name : "anim"s) +
                            std::to_string(aid++);
                anm->group = ganm->name;
                auto input_view =
                    accessor_view(gltf, gltf->get(gsampler->input));
                anm->times.resize(input_view.size());
                for (auto i = 0; i < input_view.size(); i++)
                    anm->times[i] = input_view.get(i);
                anm->type = keyframe_types.at(gsampler->interpolation);
                auto output_view =
                    accessor_view(gltf, gltf->get(gsampler->output));
                switch (gchannel->target->path) {
                    case glTFAnimationChannelTargetPath::Translation: {
                        anm->translation.reserve(output_view.size());
                        for (auto i = 0; i < output_view.size(); i++)
                            anm->translation.push_back(output_view.getv3f(i));
                    } break;
                    case glTFAnimationChannelTargetPath::Rotation: {
                        anm->rotation.reserve(output_view.size());
                        for (auto i = 0; i < output_view.size(); i++)
                            anm->rotation.push_back(output_view.getv4f(i));
                    } break;
                    case glTFAnimationChannelTargetPath::Scale: {
                        anm->scale.reserve(output_view.size());
                        for (auto i = 0; i < output_view.size(); i++)
                            anm->scale.push_back(output_view.getv3f(i));
                    } break;
                    case glTFAnimationChannelTargetPath::Weights: {
                        // get a node that it refers to
                        auto ncomp = 0;
                        auto gnode = gltf->get(gchannel->target->node);
                        auto gmesh = gltf->get(gnode->mesh);
                        if (gmesh) {
                            for (auto gshp : gmesh->primitives) {
                                ncomp = max((int)gshp->targets.size(), ncomp);
                            }
                        }
                        if (ncomp) {
                            auto values = std::vector<float>();
                            values.reserve(output_view.size());
                            for (auto i = 0; i < output_view.size(); i++)
                                values.push_back(output_view.get(i));
                            anm->weights.resize(values.size() / ncomp);
                            for (auto i = 0; i < anm->weights.size(); i++) {
                                anm->weights[i].resize(ncomp);
                                for (auto j = 0; j < ncomp; j++)
                                    anm->weights[i][j] = values[i * ncomp + j];
                            }
                        }
                    } break;
                    default: {
                        throw std::runtime_error("should not have gotten here");
                    }
                }
                sampler_map[{(int)gchannel->sampler,
                    (int)gchannel->target->path}] = (int)scn->animations.size();
                scn->animations.push_back(anm);
            }
            scn->animations[sampler_map.at({(int)gchannel->sampler,
                                (int)gchannel->target->path})]
                ->targets.push_back(scn->nodes[(int)gchannel->target->node]);
        }
    }

    // compute transforms
    update_transforms(scn, 0);

    return scn;
}

// Load an gltf scene
scene* load_gltf_scene(const std::string& filename, const load_options& opts) {
    auto gscn =
        load_gltf(filename, true, opts.load_textures, opts.skip_missing);
    auto scn = gltf_to_scene(gscn, opts);
    delete gscn;
    if (!scn) {
        throw std::runtime_error("could not convert gltf scene");
        return nullptr;
    }
    return scn;
}

// Unflattnes gltf
glTF* scene_to_gltf(
    const scene* scn, const std::string& buffer_uri, bool separate_buffers) {
    auto gltf = new glTF();

    // add asset info
    gltf->asset = new glTFAsset();
    gltf->asset->generator = "Yocto/gltf";
    gltf->asset->version = "2.0";

    // convert cameras
    for (auto cam : scn->cameras) {
        auto gcam = new glTFCamera();
        gcam->name = cam->name;
        gcam->type = (cam->ortho) ? glTFCameraType::Orthographic :
                                    glTFCameraType::Perspective;
        if (cam->ortho) {
            auto ortho = new glTFCameraOrthographic();
            ortho->ymag = cam->yfov;
            ortho->xmag = cam->aspect * cam->yfov;
            ortho->znear = cam->near;
            ortho->znear = cam->far;
            gcam->orthographic = ortho;
        } else {
            auto persp = new glTFCameraPerspective();
            persp->yfov = cam->yfov;
            persp->aspectRatio = cam->aspect;
            persp->znear = cam->near;
            persp->zfar = cam->far;
            gcam->perspective = persp;
        }
        gltf->cameras.push_back(gcam);
    }

    // convert images
    for (auto txt : scn->textures) {
        auto gimg = new glTFImage();
        gimg->uri = txt->path;
        if (!txt->hdr.pixels.empty()) {
            gimg->data.width = txt->hdr.width;
            gimg->data.height = txt->hdr.height;
            gimg->data.ncomp = 4;
            gimg->data.dataf.assign((float*)txt->hdr.pixels.data(),
                (float*)txt->hdr.pixels.data() +
                    txt->hdr.width * txt->hdr.height * 4);
        }
        if (!txt->ldr.pixels.empty()) {
            gimg->data.width = txt->ldr.width;
            gimg->data.height = txt->ldr.height;
            gimg->data.ncomp = 4;
            gimg->data.datab.assign((uint8_t*)txt->ldr.pixels.data(),
                (uint8_t*)txt->ldr.pixels.data() +
                    txt->ldr.width * txt->ldr.height * 4);
        }
        gltf->images.push_back(gimg);
    }

    // index of an object
    auto index = [](const auto& vec, auto& val) -> int {
        auto pos = find(vec.begin(), vec.end(), val);
        if (pos == vec.end()) return -1;
        return (int)(pos - vec.begin());
    };

    // add a texture and sampler
    auto add_texture_info = [&gltf, &index, scn](const texture_info& info,
                                bool norm = false, bool occ = false) {
        if (!info.txt) return (glTFTextureInfo*)nullptr;
        auto gtxt = new glTFTexture();
        gtxt->name = info.txt->name;
        gtxt->source = glTFid<glTFImage>(index(scn->textures, info.txt));

        // check if it is default
        auto is_default =
            info.wrap_s && info.wrap_t && info.linear && info.mipmap;
        if (!is_default) {
            auto gsmp = new glTFSampler();
            gsmp->wrapS = (info.wrap_s) ? glTFSamplerWrapS::Repeat :
                                          glTFSamplerWrapS::ClampToEdge;
            gsmp->wrapT = (info.wrap_t) ? glTFSamplerWrapT::Repeat :
                                          glTFSamplerWrapT::ClampToEdge;
            gsmp->minFilter = (info.mipmap) ?
                                  glTFSamplerMinFilter::LinearMipmapLinear :
                                  glTFSamplerMinFilter::Nearest;
            gsmp->magFilter = (info.linear) ? glTFSamplerMagFilter::Linear :
                                              glTFSamplerMagFilter::Nearest;
            gtxt->sampler = glTFid<glTFSampler>((int)gltf->samplers.size());
            gltf->samplers.push_back(gsmp);
        }
        gltf->textures.push_back(gtxt);
        if (norm) {
            auto ginfo = new glTFMaterialNormalTextureInfo();
            ginfo->index = glTFid<glTFTexture>{(int)gltf->textures.size() - 1};
            ginfo->scale = info.scale;
            return (glTFTextureInfo*)ginfo;
        } else if (occ) {
            auto ginfo = new glTFMaterialOcclusionTextureInfo();
            ginfo->index = glTFid<glTFTexture>{(int)gltf->textures.size() - 1};
            ginfo->strength = info.scale;
            return (glTFTextureInfo*)ginfo;
        } else {
            auto ginfo = new glTFTextureInfo();
            ginfo->index = glTFid<glTFTexture>{(int)gltf->textures.size() - 1};
            return ginfo;
        }
    };

    // convert materials
    for (auto mat : scn->materials) {
        auto gmat = new glTFMaterial();
        gmat->name = mat->name;
        gmat->emissiveFactor = mat->ke;
        gmat->emissiveTexture = add_texture_info(mat->ke_txt);
        switch (mat->type) {
            case material_type::specular_roughness: {
                gmat->pbrSpecularGlossiness =
                    new glTFMaterialPbrSpecularGlossiness();
                auto gsg = gmat->pbrSpecularGlossiness;
                gsg->diffuseFactor = {mat->kd.x, mat->kd.y, mat->kd.z, mat->op};
                gsg->specularFactor = mat->ks;
                gsg->glossinessFactor = 1 - mat->rs;
                gsg->diffuseTexture = add_texture_info(mat->kd_txt);
                gsg->specularGlossinessTexture = add_texture_info(mat->ks_txt);
            } break;
            case material_type::metallic_roughness: {
                gmat->pbrMetallicRoughness =
                    new glTFMaterialPbrMetallicRoughness();
                auto gmr = gmat->pbrMetallicRoughness;
                gmr->baseColorFactor = {
                    mat->kd.x, mat->kd.y, mat->kd.z, mat->op};
                gmr->metallicFactor = mat->ks.x;
                gmr->roughnessFactor = mat->rs;
                gmr->baseColorTexture = add_texture_info(mat->kd_txt);
                gmr->metallicRoughnessTexture = add_texture_info(mat->ks_txt);
            } break;
            case material_type::specular_glossiness: {
                gmat->pbrSpecularGlossiness =
                    new glTFMaterialPbrSpecularGlossiness();
                auto gsg = gmat->pbrSpecularGlossiness;
                gsg->diffuseFactor = {mat->kd.x, mat->kd.y, mat->kd.z, mat->op};
                gsg->specularFactor = mat->ks;
                gsg->glossinessFactor = mat->rs;
                gsg->diffuseTexture = add_texture_info(mat->kd_txt);
                gsg->specularGlossinessTexture = add_texture_info(mat->ks_txt);
            } break;
        }
        gmat->normalTexture = (glTFMaterialNormalTextureInfo*)(add_texture_info(
            mat->norm_txt, true, false));
        gmat->occlusionTexture =
            (glTFMaterialOcclusionTextureInfo*)(add_texture_info(
                mat->occ_txt, false, true));
        gmat->doubleSided = mat->double_sided;
        gltf->materials.push_back(gmat);
    }

    // add buffer
    auto add_buffer = [&gltf](const std::string& buffer_uri) {
        auto gbuffer = new glTFBuffer();
        gltf->buffers.push_back(gbuffer);
        gbuffer->uri = buffer_uri;
        return gbuffer;
    };

    // init buffers
    auto gbuffer_global = add_buffer(buffer_uri);

    // add an optional buffer
    auto add_opt_buffer = [&gbuffer_global, buffer_uri, &add_buffer,
                              separate_buffers](const std::string& uri) {
        if (separate_buffers && uri != "") {
            return add_buffer(uri);
        } else {
            if (!gbuffer_global) gbuffer_global = add_buffer(buffer_uri);
            return gbuffer_global;
        }
    };

    // attribute handling
    auto add_accessor = [&gltf, &index](glTFBuffer* gbuffer,
                            const std::string& name, glTFAccessorType type,
                            glTFAccessorComponentType ctype, int count,
                            int csize, const void* data, bool save_min_max) {
        gltf->bufferViews.push_back(new glTFBufferView());
        auto bufferView = gltf->bufferViews.back();
        bufferView->buffer = glTFid<glTFBuffer>(index(gltf->buffers, gbuffer));
        bufferView->byteOffset = (int)gbuffer->data.size();
        bufferView->byteStride = 0;
        bufferView->byteLength = count * csize;
        gbuffer->data.resize(gbuffer->data.size() + bufferView->byteLength);
        gbuffer->byteLength += bufferView->byteLength;
        auto ptr = gbuffer->data.data() + gbuffer->data.size() -
                   bufferView->byteLength;
        bufferView->target = glTFBufferViewTarget::ArrayBuffer;
        memcpy(ptr, data, bufferView->byteLength);
        gltf->accessors.push_back(new glTFAccessor());
        auto accessor = gltf->accessors.back();
        accessor->bufferView =
            glTFid<glTFBufferView>((int)gltf->bufferViews.size() - 1);
        accessor->byteOffset = 0;
        accessor->componentType = ctype;
        accessor->count = count;
        accessor->type = type;
        if (save_min_max && count &&
            ctype == glTFAccessorComponentType::Float) {
            float dmin[4] = {flt_max, flt_max, flt_max, flt_max};
            float dmax[4] = {flt_min, flt_min, flt_min, flt_min};
            auto d = (float*)data;
            auto nc = 0;
            switch (type) {
                case glTFAccessorType::Scalar: nc = 1; break;
                case glTFAccessorType::Vec2: nc = 2; break;
                case glTFAccessorType::Vec3: nc = 3; break;
                case glTFAccessorType::Vec4: nc = 4; break;
                default: break;
            }
            for (auto i = 0; i < count; i++) {
                for (auto c = 0; c < nc; c++) {
                    dmin[c] = min(dmin[c], d[i * nc + c]);
                    dmax[c] = max(dmax[c], d[i * nc + c]);
                }
            }
            for (auto c = 0; c < nc; c++) {
                accessor->min.push_back(dmin[c]);
                accessor->max.push_back(dmax[c]);
            }
        }
        return glTFid<glTFAccessor>((int)gltf->accessors.size() - 1);
    };

    // instances
    auto shape_mats = std::map<shape*, material*>();
    for (auto shp : scn->shapes) shape_mats[shp] = nullptr;
    for (auto ist : scn->instances) {
        if (shape_mats.at(ist->shp) && shape_mats.at(ist->shp) != ist->mat)
            log_error("shapes can only have one material associated");
        else
            shape_mats[ist->shp] = ist->mat;
    }

    // convert meshes
    for (auto shp : scn->shapes) {
        auto gmesh = new glTFMesh();
        gmesh->name = shp->name;
        auto gbuffer = add_opt_buffer(shp->path);
        auto gprim = new glTFMeshPrimitive();
        gprim->material =
            glTFid<glTFMaterial>(index(scn->materials, shape_mats.at(shp)));
        if (!shp->pos.empty())
            gprim->attributes["POSITION"] =
                add_accessor(gbuffer, shp->name + "_pos",
                    glTFAccessorType::Vec3, glTFAccessorComponentType::Float,
                    (int)shp->pos.size(), sizeof(vec3f), shp->pos.data(), true);
        if (!shp->norm.empty())
            gprim->attributes["NORMAL"] = add_accessor(gbuffer,
                shp->name + "_norm", glTFAccessorType::Vec3,
                glTFAccessorComponentType::Float, (int)shp->norm.size(),
                sizeof(vec3f), shp->norm.data(), false);
        if (!shp->texcoord.empty())
            gprim->attributes["TEXCOORD_0"] = add_accessor(gbuffer,
                shp->name + "_texcoord", glTFAccessorType::Vec2,
                glTFAccessorComponentType::Float, (int)shp->texcoord.size(),
                sizeof(vec2f), shp->texcoord.data(), false);
        if (!shp->texcoord1.empty())
            gprim->attributes["TEXCOORD_1"] = add_accessor(gbuffer,
                shp->name + "_texcoord1", glTFAccessorType::Vec2,
                glTFAccessorComponentType::Float, (int)shp->texcoord1.size(),
                sizeof(vec2f), shp->texcoord1.data(), false);
        if (!shp->color.empty())
            gprim->attributes["COLOR_0"] = add_accessor(gbuffer,
                shp->name + "_color", glTFAccessorType::Vec4,
                glTFAccessorComponentType::Float, (int)shp->color.size(),
                sizeof(vec4f), shp->color.data(), false);
        if (!shp->radius.empty())
            gprim->attributes["RADIUS"] = add_accessor(gbuffer,
                shp->name + "_radius", glTFAccessorType::Scalar,
                glTFAccessorComponentType::Float, (int)shp->radius.size(),
                sizeof(float), shp->radius.data(), false);
        // auto elem_as_uint = shp->pos.size() >
        // numeric_limits<unsigned short>::max();
        if (!shp->points.empty()) {
            gprim->indices = add_accessor(gbuffer, shp->name + "_points",
                glTFAccessorType::Scalar,
                glTFAccessorComponentType::UnsignedInt, (int)shp->points.size(),
                sizeof(int), (int*)shp->points.data(), false);
            gprim->mode = glTFMeshPrimitiveMode::Points;
        } else if (!shp->lines.empty()) {
            gprim->indices = add_accessor(gbuffer, shp->name + "_lines",
                glTFAccessorType::Scalar,
                glTFAccessorComponentType::UnsignedInt,
                (int)shp->lines.size() * 2, sizeof(int),
                (int*)shp->lines.data(), false);
            gprim->mode = glTFMeshPrimitiveMode::Lines;
        } else if (!shp->triangles.empty()) {
            gprim->indices = add_accessor(gbuffer, shp->name + "_triangles",
                glTFAccessorType::Scalar,
                glTFAccessorComponentType::UnsignedInt,
                (int)shp->triangles.size() * 3, sizeof(int),
                (int*)shp->triangles.data(), false);
            gprim->mode = glTFMeshPrimitiveMode::Triangles;
        } else if (!shp->quads.empty()) {
            auto triangles = convert_quads_to_triangles(shp->quads);
            gprim->indices = add_accessor(gbuffer, shp->name + "_quads",
                glTFAccessorType::Scalar,
                glTFAccessorComponentType::UnsignedInt,
                (int)triangles.size() * 3, sizeof(int), (int*)triangles.data(),
                false);
            gprim->mode = glTFMeshPrimitiveMode::Triangles;
        } else if (!shp->quads_pos.empty()) {
            throw std::runtime_error("face varying not supported in glTF");
        } else {
            throw std::runtime_error("empty mesh");
        }
        gmesh->primitives.push_back(gprim);
        gltf->meshes.push_back(gmesh);
    }

    // hierarchy
    if (scn->nodes.empty()) {
        // shapes
        for (auto ist : scn->instances) {
            auto gnode = new glTFNode();
            gnode->name = ist->name;
            gnode->mesh = glTFid<glTFMesh>(index(scn->shapes, ist->shp));
            gnode->matrix = frame_to_mat(ist->frame);
            gltf->nodes.push_back(gnode);
        }

        // cameras
        for (auto cam : scn->cameras) {
            auto gnode = new glTFNode();
            gnode->name = cam->name;
            gnode->camera = glTFid<glTFCamera>(index(scn->cameras, cam));
            gnode->matrix = frame_to_mat(cam->frame);
            gltf->nodes.push_back(gnode);
        }

        // scenes
        if (!gltf->nodes.empty()) {
            auto gscene = new glTFScene();
            gscene->name = "scene";
            for (auto i = 0; i < gltf->nodes.size(); i++) {
                gscene->nodes.push_back(glTFid<glTFNode>(i));
            }
            gltf->scenes.push_back(gscene);
            gltf->scene = glTFid<glTFScene>(0);
        }
    } else {
        for (auto nde : scn->nodes) {
            auto gnode = new glTFNode();
            gnode->name = nde->name;
            if (nde->cam) {
                gnode->camera =
                    glTFid<glTFCamera>(index(scn->cameras, nde->cam));
            }
            if (nde->ist) {
                gnode->mesh =
                    glTFid<glTFMesh>(index(scn->shapes, nde->ist->shp));
            }
            gnode->matrix = frame_to_mat(nde->frame);
            gnode->translation = nde->translation;
            gnode->rotation = nde->rotation;
            gnode->scale = nde->scale;
            gltf->nodes.push_back(gnode);
        }

        // children
        for (auto idx = 0; idx < scn->nodes.size(); idx++) {
            auto nde = scn->nodes.at(idx);
            if (!nde->parent) continue;
            auto gnde = gltf->nodes.at(index(scn->nodes, nde->parent));
            gnde->children.push_back(glTFid<glTFNode>(idx));
        }

        // root nodes
        auto is_root = std::vector<bool>(gltf->nodes.size(), true);
        for (auto idx = 0; idx < gltf->nodes.size(); idx++) {
            auto gnde = gltf->nodes.at(idx);
            for (auto idx1 = 0; idx1 < gnde->children.size(); idx1++) {
                is_root[(int)gnde->children.at(idx1)] = false;
            }
        }

        // scene with root nodes
        auto gscene = new glTFScene();
        gscene->name = "scene";
        for (auto idx = 0; idx < gltf->nodes.size(); idx++) {
            if (is_root[idx]) gscene->nodes.push_back(glTFid<glTFNode>(idx));
        }
        gltf->scenes.push_back(gscene);
        gltf->scene = glTFid<glTFScene>(0);
    }

    // interpolation map
    static const auto interpolation_map =
        std::map<animation_type, glTFAnimationSamplerInterpolation>{
            {animation_type::step, glTFAnimationSamplerInterpolation::Step},
            {animation_type::linear, glTFAnimationSamplerInterpolation::Linear},
            {animation_type::bezier,
                glTFAnimationSamplerInterpolation::CubicSpline},
        };

    // gruop animations
    struct anim_group {
        std::string path;
        std::string name;
        std::vector<animation*> animations;
    };

    std::map<std::string, anim_group> anim_groups;
    for (auto anm : scn->animations) {
        auto agr = &anim_groups[anm->group];
        if (agr->path == "") agr->path = anm->path;
        agr->animations.push_back(anm);
    }

    // animation
    for (auto& agr_kv : anim_groups) {
        auto agr = &agr_kv.second;
        auto ganm = new glTFAnimation();
        ganm->name = agr->name;
        auto gbuffer = add_opt_buffer(agr->path);
        auto count = 0;
        for (auto anm : scn->animations) {
            auto aid = ganm->name + "_" + std::to_string(count++);
            auto gsmp = new glTFAnimationSampler();
            gsmp->input =
                add_accessor(gbuffer, aid + "_time", glTFAccessorType::Scalar,
                    glTFAccessorComponentType::Float, (int)anm->times.size(),
                    sizeof(float), anm->times.data(), false);
            auto path = glTFAnimationChannelTargetPath::NotSet;
            if (!anm->translation.empty()) {
                gsmp->output = add_accessor(gbuffer, aid + "_translation",
                    glTFAccessorType::Vec3, glTFAccessorComponentType::Float,
                    (int)anm->translation.size(), sizeof(vec3f),
                    anm->translation.data(), false);
                path = glTFAnimationChannelTargetPath::Translation;
            } else if (!anm->rotation.empty()) {
                gsmp->output = add_accessor(gbuffer, aid + "_rotation",
                    glTFAccessorType::Vec4, glTFAccessorComponentType::Float,
                    (int)anm->rotation.size(), sizeof(vec4f),
                    anm->rotation.data(), false);
                path = glTFAnimationChannelTargetPath::Rotation;
            } else if (!anm->scale.empty()) {
                gsmp->output = add_accessor(gbuffer, aid + "_scale",
                    glTFAccessorType::Vec3, glTFAccessorComponentType::Float,
                    (int)anm->scale.size(), sizeof(vec3f), anm->scale.data(),
                    false);
                path = glTFAnimationChannelTargetPath::Scale;
            } else if (!anm->weights.empty()) {
                auto values = std::vector<float>();
                values.reserve(anm->weights.size() * anm->weights[0].size());
                for (auto i = 0; i < anm->weights.size(); i++) {
                    values.insert(values.end(), anm->weights[i].begin(),
                        anm->weights[i].end());
                }
                gsmp->output = add_accessor(gbuffer, aid + "_weights",
                    glTFAccessorType::Scalar, glTFAccessorComponentType::Float,
                    (int)values.size(), sizeof(float), values.data(), false);
                path = glTFAnimationChannelTargetPath::Weights;
            } else {
                throw std::runtime_error("should not have gotten here");
            }
            gsmp->interpolation = interpolation_map.at(anm->type);
            for (auto target : anm->targets) {
                auto gchan = new glTFAnimationChannel();
                gchan->sampler =
                    glTFid<glTFAnimationSampler>{(int)ganm->samplers.size()};
                gchan->target = new glTFAnimationChannelTarget();
                gchan->target->node =
                    glTFid<glTFNode>{index(scn->nodes, target)};
                gchan->target->path = path;
                ganm->channels.push_back(gchan);
            }
            ganm->samplers.push_back(gsmp);
        }

        gltf->animations.push_back(ganm);
    }

    // done
    return gltf;
}

// Save a gltf scene
void save_gltf_scene(
    const std::string& filename, const scene* scn, const save_options& opts) {
    auto buffer_uri = path_basename(filename) + ".bin";
    auto gscn = scene_to_gltf(scn, buffer_uri, opts.gltf_separate_buffers);
    save_gltf(filename, gscn, true, opts.save_textures);
}

#endif

// Load a scene
scene* load_scene(const std::string& filename, const load_options& opts) {
    auto ext = path_extension(filename);
    if (ext == ".obj" || ext == ".OBJ") return load_obj_scene(filename, opts);
#if YGL_GLTF
    if (ext == ".gltf" || ext == ".GLTF")
        return load_gltf_scene(filename, opts);
#endif
#if YGL_SVG
    if (ext == ".svg" || ext == ".SVG") return load_svg_scene(filename, opts);
#endif
    throw std::runtime_error("unsupported extension " + ext);
    return nullptr;
}

// Save a scene
void save_scene(
    const std::string& filename, const scene* scn, const save_options& opts) {
    auto ext = path_extension(filename);
    if (ext == ".obj" || ext == ".OBJ")
        return save_obj_scene(filename, scn, opts);
#if YGL_GLTF
    if (ext == ".gltf" || ext == ".GLTF")
        return save_gltf_scene(filename, scn, opts);
#endif
    throw std::runtime_error("unsupported extension " + ext);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR EXAMPLE SCENES
// -----------------------------------------------------------------------------
namespace ygl {

camera* make_camera(const std::string& name, const vec3f& from, const vec3f& to,
    float yfov, float aspect) {
    auto cam = new camera();
    cam->name = name;
    cam->frame = lookat_frame(from, to, vec3f{0, 1, 0});
    cam->yfov = yfov;
    cam->aspect = aspect;
    cam->near = 0.01f;
    cam->far = 10000;
    cam->aperture = 0;
    cam->focus = length(from - to);
    return cam;
};

material* make_material(
    const std::string& name, const vec3f& kd, const vec3f& ks, float rs) {
    auto mat = new material();
    mat->name = name;
    mat->kd = kd;
    mat->ks = ks;
    mat->rs = rs;
    return mat;
}

texture* make_texture(const std::string& name, const std::string& path,
    const image4b& ldr, const image4f& hdr) {
    auto txt = new texture();
    txt->name = name;
    txt->path = path;
    txt->ldr = ldr;
    txt->hdr = hdr;
    return txt;
}

camera* make_proc_camera(const std::string& name, const std::string& type) {
    if (type == "") {
        return make_camera(name, {0, 0, 1}, {0, 0, 0}, 75 * pi / 180, 1);
    } else if (type == "cam1") {
        return make_camera(name, {0, 4, 10}, {0, 1, 0}, 15 * pi / 180, 1);
    } else if (type == "cam1") {
        return make_camera(
            name, {0, 4, 10}, {0, 1, 0}, 15 * pi / 180, 16.0f / 9.0f);
    } else if (type == "cam1") {
        return make_camera(
            name, {0, 6, 24}, {0, 1, 0}, 7.5f * pi / 180, 2.35f / 1.0f);
    } else {
        throw std::runtime_error("unknown camera type");
    }
}

// Makes/updates a test texture
texture* make_proc_texture(const std::string& name, const std::string& type,
    int res, float scale, float sky_sunangle, float bump_scale) {
    auto txt = new texture();

    txt->name = name;
    txt->path = "";

    if (type == "") {
    } else if (type == "grid") {
        auto ts = (int)(res / scale + 0.5f);
        txt->ldr = make_grid_image(res, res, ts);
    } else if (type == "checker") {
        auto ts = (int)(res / scale + 0.5f);
        txt->ldr = make_checker_image(res, res, ts);
    } else if (type == "colored") {
        auto ts = (int)(res / scale + 0.5f);
        txt->ldr = make_uvgrid_image(res, res, ts);
    } else if (type == "rcolored") {
        auto ts = (int)(res / scale + 0.5f);
        txt->ldr = make_recuvgrid_image(res, res, ts);
    } else if (type == "bump") {
        auto ts = (int)(res / (2 * scale) + 0.5f);
        txt->ldr = make_bumpdimple_image(res, res, ts);
    } else if (type == "uv") {
        txt->ldr = make_uv_image(res, res);
    } else if (type == "gamma") {
        txt->ldr = make_gammaramp_image(res, res);
    } else if (type == "noise") {
        txt->ldr = make_noise_image(res, res, scale);
    } else if (type == "ridge") {
        txt->ldr = make_ridge_image(res, res, scale);
    } else if (type == "fbm") {
        txt->ldr = make_fbm_image(res, res, scale);
    } else if (type == "turbulence") {
        txt->ldr = make_turbulence_image(res, res, scale);
    } else if (type == "grid_norm") {
        auto ts = (int)(res / scale + 0.5f);
        txt->ldr = make_grid_image(res, res, ts);
        txt->ldr = bump_to_normal_map(txt->ldr, bump_scale);
    } else if (type == "bump_norm") {
        auto ts = (int)(res / (2 * scale) + 0.5f);
        txt->ldr = make_bumpdimple_image(res, res, ts);
        txt->ldr = bump_to_normal_map(txt->ldr, bump_scale);
    } else if (type == "gammaf") {
        txt->hdr = make_gammaramp_imagef(res, res);
    } else if (type == "sky") {
        txt->hdr = make_sunsky_image(res, sky_sunangle);
    } else {
        throw std::runtime_error("unknown texture type " + type);
    }

    if (!txt->ldr.pixels.empty()) txt->path = name + ".png";
    if (!txt->hdr.pixels.empty()) txt->path = name + ".hdr";

    return txt;
}

std::map<std::string, texture*>& texture_presets() {
    static auto presets = std::map<std::string, texture*>();
    if (!presets.empty()) return presets;

    auto res = 512;

    presets["grid"] = make_proc_texture("grid", "grid", res);
    presets["checker"] = make_proc_texture("checker", "checker", res);
    presets["colored"] = make_proc_texture("colored", "colored", res);
    presets["rcolored"] = make_proc_texture("rcolored", "rcolored", res);
    presets["bump"] = make_proc_texture("bump", "bump", res);
    presets["uv"] = make_proc_texture("uv", "uv", res);
    presets["gamma"] = make_proc_texture("gamma", "gamma", res);
    presets["grid_norm"] = make_proc_texture("grid_norm", "grid_norm");
    presets["bump_norm"] = make_proc_texture("bump_norm", "bump_norm");
    presets["noise"] = make_proc_texture("noise", "noise", res);
    presets["ridge"] = make_proc_texture("ridge", "ridge", res);
    presets["fbm"] = make_proc_texture("fbm", "fbm", res);
    presets["turbulence"] = make_proc_texture("turbulence", "turbulence", res);

    presets["gammaf"] = make_proc_texture("gammaf", "gammaf", res);
    presets["sky1"] = make_proc_texture("sky1", "sky", res, 0, pi / 4);
    presets["sky2"] = make_proc_texture("sky2", "sky", res, 0, pi / 2);

    return presets;
}

// Makes/updates a test material
material* make_proc_material(const std::string& name, const std::string& type_,
    const vec3f& col_, float rs, const std::string& txt_type_,
    const std::string& norm_type_) {
    auto erase = [](std::string& str, const std::string& substr) {
        auto pos = str.find(substr);
        if (pos == str.npos) return false;
        str = str.erase(pos, substr.size());
        return true;
    };

    auto mat = new material();

    mat->name = name;

    auto col = col_;
    auto type = type_, txt_type = txt_type_, norm_type = norm_type_;
    for (auto txt_name : {"grid_norm", "bump_norm"})
        if (erase(type, "_"s + txt_name)) norm_type = txt_name;
    for (auto txt_name : {"grid", "checker", "colored", "rcolored", "bump",
             "uv", "gamma", "noise", "ridge", "fbm", "turbulence"})
        if (erase(type, "_"s + txt_name)) txt_type = txt_name;

    if (erase(type, "_sharp")) rs = 0.05f;
    if (erase(type, "_rough")) rs = 0.25f;
    if (erase(type, "_mirror")) rs = 0;

    if (erase(type, "_zero")) col = {0, 0, 0};
    if (erase(type, "_one")) col = {1, 1, 1};
    if (erase(type, "_black")) col = {0.01f, 0.01f, 0.01f};
    if (erase(type, "_gray")) col = {0.2f, 0.2f, 0.2f};
    if (erase(type, "_lgray")) col = {0.5f, 0.5f, 0.5f};
    if (erase(type, "_red")) col = {0.5f, 0.2f, 0.2f};
    if (erase(type, "_green")) col = {0.2f, 0.5f, 0.2f};
    if (erase(type, "_blue")) col = {0.2f, 0.2f, 0.5f};
    if (erase(type, "_white")) col = {0.9f, 0.9f, 0.9f};
    if (erase(type, "_gold")) col = {0.66f, 0.45f, 0.34f};
    if (erase(type, "_silver")) col = {0.7f, 0.7f, 0.7f};

    auto txt = (texture*)nullptr, norm = (texture*)nullptr;
    if (txt_type != "") txt = make_proc_texture(txt_type, txt_type, 512);
    if (norm_type != "") norm = make_proc_texture(norm_type, norm_type, 512);

    auto def_rs = rs == 1;

    if (type == "") {
    } else if (type == "emission") {
        mat->ke = col;
        mat->ke_txt.txt = txt;
    } else if (type == "matte") {
        mat->kd = col;
        mat->kd_txt.txt = txt;
        mat->rs = 1;
    } else if (type == "plastic") {
        mat->kd = col;
        mat->ks = {0.04f, 0.04f, 0.04f};
        mat->rs = (def_rs) ? 0.1f : rs;
        mat->kd_txt.txt = txt;
    } else if (type == "metal") {
        mat->ks = col;
        mat->rs = (def_rs) ? 0.1f : rs;
        mat->ks_txt.txt = txt;
    } else if (type == "carpaint") {
        mat->ks = col;
        mat->rs = (def_rs) ? 0.1f : rs;
        mat->ks_txt.txt = txt;
        mat->kr = {0.04f, 0.04f, 0.04f};
    } else if (type == "glass") {
        mat->ks = {0.04f, 0.04f, 0.04f};
        mat->kt = col;
        mat->rs = (def_rs) ? 0.1f : rs;
        mat->kt_txt.txt = txt;
    } else if (type == "transparent") {
        mat->kd = col;
        mat->op = (def_rs) ? 0.5f : rs;
        mat->kd_txt.txt = txt;
    } else {
        throw std::runtime_error("unknown material type " + type);
    }

    mat->norm_txt.txt = norm;

    return mat;
}

// Makes/updates a test shape
shape* make_proc_shape(const std::string& name, const std::string& type_,
    const vec3i& tesselation_, const vec3f& size_, const vec3f& uvsize_,
    float rounded, float radius, const make_hair_params& hair_params_) {
    auto shp = new shape();

    auto tesselation = tesselation_;
    auto def_tesselation = (tesselation == zero3i);

    auto type = type_;
    auto hair_params = hair_params_;
    if (type == "hairball_noise") {
        hair_params.noise = {0.5f, 8};
        type = "hairball";
    }
    if (type == "hairball_clump") {
        hair_params.clump = {0.5f, 128};
        type = "hairball";
    }

    auto size = size_;
    auto def_size = (size == zero3f);
    if (def_size) size = {2, 2, 2};

    auto uvsize = uvsize_;
    auto def_uvsize = (uvsize == zero3f);
    if (def_size) uvsize = {1, 1, 1};

    shp->name = name;
    if (type == "") {
    } else if (type == "floor") {
        if (def_size) size = {40, 40, 40};
        if (def_uvsize) uvsize = {20, 20, 20};
        make_quad(shp->quads, shp->pos, shp->norm, shp->texcoord,
            {1 << tesselation.x, 1 << tesselation.y}, {size.x, size.y},
            {uvsize.x, uvsize.y});
        for (auto& p : shp->pos) p = {-p.x, p.z, p.y};
        for (auto& n : shp->norm) n = {n.x, n.z, n.y};
    } else if (type == "quad") {
        make_quad(shp->quads, shp->pos, shp->norm, shp->texcoord,
            {1 << tesselation.x, 1 << tesselation.y}, {size.x, size.y},
            {uvsize.x, uvsize.y});
    } else if (type == "cube") {
        make_cube(shp->quads, shp->pos, shp->norm, shp->texcoord,
            {1 << tesselation.x, 1 << tesselation.y, 1 << tesselation.z}, size,
            uvsize);
    } else if (type == "cube_rounded") {
        if (def_tesselation) tesselation = {5, 5, 5};
        make_cube_rounded(shp->quads, shp->pos, shp->norm, shp->texcoord,
            {1 << tesselation.x, 1 << tesselation.y, 1 << tesselation.z}, size,
            uvsize, min(size) * (1 - rounded) / 2);
    } else if (type == "sphere") {
        if (def_tesselation) tesselation = {6, 5, 5};
        make_sphere(shp->quads, shp->pos, shp->norm, shp->texcoord,
            {1 << tesselation.x, 1 << tesselation.y}, size.x,
            {uvsize.x, uvsize.y});
    } else if (type == "sphere_flipcap") {
        if (def_tesselation) tesselation = {6, 5, 5};
        make_sphere_flipcap(shp->quads, shp->pos, shp->norm, shp->texcoord,
            {1 << tesselation.x, 1 << tesselation.y}, size.x,
            {uvsize.x, uvsize.y}, {-rounded, rounded});
    } else if (type == "sphere_cube") {
        if (def_tesselation) tesselation = {5, 5, 5};
        make_sphere_cube(shp->quads, shp->pos, shp->norm, shp->texcoord,
            1 << tesselation.x, size.x, uvsize.x);
    } else if (type == "disk") {
        if (def_tesselation) tesselation = {5, 5, 5};
        make_disk(shp->quads, shp->pos, shp->norm, shp->texcoord,
            {1 << tesselation.x, 1 << tesselation.y}, size.x,
            {uvsize.x, uvsize.y});
    } else if (type == "disk_quad") {
        if (def_tesselation) tesselation = {5, 5, 5};
        make_disk_quad(shp->quads, shp->pos, shp->norm, shp->texcoord,
            1 << tesselation.x, size.x, uvsize.x);
    } else if (type == "disk_bulged") {
        if (def_tesselation) tesselation = {5, 5, 5};
        make_disk_bulged(shp->quads, shp->pos, shp->norm, shp->texcoord,
            1 << tesselation.x, size.x, uvsize.x, rounded);
    } else if (type == "cylinder") {
        if (def_tesselation) tesselation = {5, 5, 5};
        make_cylinder(shp->quads, shp->pos, shp->norm, shp->texcoord,
            {1 << tesselation.x, 1 << tesselation.y, 1 << tesselation.z},
            {size.x, size.y}, uvsize);
    } else if (type == "cylinder_rounded") {
        if (def_tesselation) tesselation = {5, 5, 5};
        make_cylinder_rounded(shp->quads, shp->pos, shp->norm, shp->texcoord,
            {1 << tesselation.x, 1 << tesselation.y, 1 << tesselation.z},
            {size.x, size.y}, uvsize, rounded);
    } else if (type == "cylindery") {
        if (def_tesselation) tesselation = {5, 5, 5};
        make_cylinder(shp->quads, shp->pos, shp->norm, shp->texcoord,
            {1 << tesselation.x, 1 << tesselation.y, 1 << tesselation.z},
            {size.x, size.y}, uvsize);
        for (auto& p : shp->pos) std::swap(p.y, p.z);
        for (auto& n : shp->norm) std::swap(n.y, n.z);
    } else if (type == "cylindery_rounded") {
        if (def_tesselation) tesselation = {5, 5, 5};
        make_cylinder_rounded(shp->quads, shp->pos, shp->norm, shp->texcoord,
            {1 << tesselation.x, 1 << tesselation.y, 1 << tesselation.z},
            {size.x, size.y}, uvsize, rounded);
        for (auto& p : shp->pos) std::swap(p.y, p.z);
        for (auto& n : shp->norm) std::swap(n.y, n.z);
    } else if (type == "geodesic_sphere") {
        if (def_tesselation) tesselation = {4, 0, 0};
        make_geodesic_sphere(shp->triangles, shp->pos, tesselation.x);
        shp->norm = shp->pos;
        shp->texcoord.assign(shp->pos.size(), {0, 0});
    } else if (type == "suzanne") {
        make_suzanne(shp->quads, shp->pos, tesselation.x);
    } else if (type == "cube_subdiv") {
        make_cube(shp->quads, shp->pos, 0);
        for (auto i = 0; i < tesselation.x; i++) {
            subdivide_catmullclark(shp->quads, shp->pos);
        }
        if (tesselation.x) compute_normals(shp->quads, shp->pos, shp->norm);
    } else if (type == "suzanne_subdiv") {
        make_suzanne(shp->quads, shp->pos, 0);
        for (auto i = 0; i < tesselation.x; i++) {
            subdivide_catmullclark(shp->quads, shp->pos);
        }
        if (tesselation.x) compute_normals(shp->quads, shp->pos, shp->norm);
    } else if (type == "fvcube_subdiv") {
        make_fvcube(shp->quads_pos, shp->pos, shp->quads_norm, shp->norm,
            shp->quads_texcoord, shp->texcoord, 0);
        for (auto i = 0; i < tesselation.x; i++) {
            subdivide_catmullclark(shp->quads_pos, shp->pos);
            subdivide_catmullclark(shp->quads_norm, shp->pos);
            subdivide_catmullclark(shp->quads_texcoord, shp->pos);
        }
    } else if (type == "matball") {
        if (def_tesselation) tesselation = {5, 4, 0};
        make_sphere_flipcap(shp->quads, shp->pos, shp->norm, shp->texcoord,
            {1 << tesselation.x, 1 << tesselation.y}, size.x,
            {uvsize.x, uvsize.y}, {-rounded, rounded});
    } else if (type == "point") {
        make_point(
            shp->points, shp->pos, shp->norm, shp->texcoord, shp->radius);
    } else if (type == "pointscube") {
        if (def_tesselation) tesselation = {16, 4, 0};
        make_random_points(shp->points, shp->pos, shp->norm, shp->texcoord,
            shp->radius, 1 << tesselation.x, size, uvsize.x, radius);
    } else if (type == "hairball") {
        if (def_tesselation) tesselation = {16, 4, 0};
        auto shp1 = new shape();
        shp1->name = "interior";
        make_sphere_cube(shp1->quads, shp1->pos, shp1->norm, shp1->texcoord,
            1 << 5, size.x * 0.8f, 1);
        shp1->radius.assign(shp1->pos.size(), 0);
        make_hair(shp->lines, shp->pos, shp->norm, shp->texcoord, shp->radius,
            {1 << tesselation.x, 1 << tesselation.y}, {}, shp1->quads,
            shp1->pos, shp1->norm, shp1->texcoord, hair_params);
        delete shp1;
    } else if (type == "beziercircle") {
        make_bezier_circle(shp->beziers, shp->pos);
        shp->subdivision = 2;
    } else {
        throw std::runtime_error("unknown shape type " + type);
    }

    return shp;
}

// makes the cornell box scene
// http://graphics.cs.williams.edu/data
// http://www.graphics.cornell.edu/online/box/data.html
scene* make_cornell_box_scene() {
    auto add_camera = [](scene* scn, std::string name, vec3f from, vec3f to,
                          float yfov, float aperture,
                          float aspect = 16.0f / 9.0f) {
        auto cam = new camera();
        cam->name = name;
        cam->frame = lookat_frame(from, to, {0, 1, 0});
        cam->aperture = aperture;
        cam->focus = length(from - to);
        cam->yfov = yfov * pi / 180;
        cam->aspect = aspect;
        scn->cameras.push_back(cam);
    };

    auto add_quad = [](scene* scn, std::string name, material* mat, vec3f pos,
                        vec3f rot = {0, 0, 0}, vec2f size = {2, 2}) {
        auto shp = new shape();
        shp->name = name;
        ygl::make_quad(shp->quads, shp->pos, shp->norm, shp->texcoord, {1, 1},
            size, {1, 1});
        scn->shapes.push_back(shp);
        auto ist = new instance();
        ist->name = name;
        ist->shp = shp;
        ist->mat = mat;
        ist->frame = translation_frame(pos);
        if (rot != zero3f) {
            ist->frame = ist->frame *
                         rotation_frame(vec3f{0, 0, 1}, rot.z * pi / 180) *
                         rotation_frame(vec3f{0, 1, 0}, rot.y * pi / 180) *
                         rotation_frame(vec3f{1, 0, 0}, rot.x * pi / 180);
        }
        scn->instances.push_back(ist);
    };

    auto add_box = [](scene* scn, std::string name, material* mat, vec3f pos,
                       vec3f rot = {0, 0, 0}, vec3f size = {2, 2, 2}) {
        auto shp = new shape();
        shp->name = name;
        shp->name = name;
        make_cube(shp->quads, shp->pos, shp->norm, shp->texcoord, {1, 1, 1},
            size, {1, 1, 1});
        scn->shapes.push_back(shp);
        auto ist = new instance();
        ist->name = name;
        ist->shp = shp;
        ist->mat = mat;
        ist->frame = translation_frame(pos);
        if (rot != zero3f) {
            ist->frame = ist->frame *
                         rotation_frame(vec3f{0, 0, 1}, rot.z * pi / 180) *
                         rotation_frame(vec3f{0, 1, 0}, rot.y * pi / 180) *
                         rotation_frame(vec3f{1, 0, 0}, rot.x * pi / 180);
        }
        scn->instances.push_back(ist);
    };

    auto add_material = [](scene* scn, std::string name, vec3f kd,
                            vec3f ke = {0, 0, 0}) {
        auto mat = new material();
        mat->type = material_type::specular_roughness;
        mat->name = name;
        mat->ke = ke;
        mat->kd = kd;
        mat->ks = zero3f;
        mat->rs = 1;
        scn->materials.push_back(mat);
        return mat;
    };

    auto scn = new scene();
    add_camera(scn, "cb_cam", {0, 1, 5.15f}, {0, 1, 0}, 27, 0, 1);
    add_material(scn, "cb_white", {0.725f, 0.71f, 0.68f});
    add_material(scn, "cb_red", {0.63f, 0.065f, 0.05f});
    add_material(scn, "cb_green", {0.14f, 0.45f, 0.091f});
    add_material(scn, "cb_light", zero3f, {17, 12, 4});
    add_quad(scn, "cb_floor", scn->materials[0], {0, 0, 0}, {-90, 0, 0});
    add_quad(scn, "cb_ceiling", scn->materials[0], {0, 2, 0}, {90, 0, 0});
    add_quad(scn, "cb_back", scn->materials[0], {0, 1, -1}, {0, 0, 0});
    add_quad(scn, "cb_left", scn->materials[2], {+1, 1, 0}, {0, -90, 0});
    add_quad(scn, "cb_right", scn->materials[1], {-1, 1, 0}, {0, 90, 0});
    add_box(scn, "cb_tallbox", scn->materials[0], {-0.33f, 0.6f, -0.29f},
        {0, 15, 0}, {0.6f, 1.2f, 0.6f});
    add_box(scn, "cb_shortbox", scn->materials[0], {0.33f, 0.3f, 0.33f},
        {0, -15, 0}, {0.6f, 0.6f, 0.6f});
    add_quad(scn, "cb_light", scn->materials[3], {0, 1.999f, 0}, {90, 0, 0},
        {0.5f, 0.5f});
    return scn;
}

instance* make_instance(
    const std::string& name, shape* shp, material* mat, const frame3f& frame) {
    auto ist = new instance();
    ist->name = name;
    ist->shp = shp;
    ist->mat = mat;
    ist->frame = frame;
    return ist;
}

instance* make_proc_instance(const std::string& name, const std::string& stype,
    const std::string& mtype, const frame3f& frame) {
    return make_instance(name, make_proc_shape(name, stype),
        make_proc_material(name, mtype), frame);
}

node* make_node(const std::string& name, camera* cam, instance* ist,
    environment* env, const frame3f& frame) {
    auto nde = new node();
    nde->name = name;
    nde->cam = cam;
    nde->ist = ist;
    nde->env = env;
    nde->frame = frame;
    return nde;
}

environment* make_environment(const std::string& name, const vec3f& ke,
    texture* ke_txt, const frame3f& frame) {
    auto env = new environment();
    env->name = name;
    env->ke = ke;
    env->ke_txt.txt = ke_txt;
    env->frame = frame;
    return env;
}

// make simple scene
scene* make_simple_scene(const std::vector<std::string>& shapes,
    const std::vector<std::string>& mats, const std::string& lights, bool nodes,
    const std::vector<std::string>& animations, const std::string& floor_mat) {
    auto frames =
        std::vector<frame3f>{{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {-2.50f, 1, 0}},
            {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 1, 0}},
            {{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {+2.50f, 1, 0}}};
    auto scn = new scene();
    scn->cameras.push_back(make_camera(
        "cam", {0, 6, 24}, {0, 1, 0}, 7.5f * pi / 180, 2.35f / 1.0f));
    if (floor_mat != "") {
        scn->instances.push_back(
            make_proc_instance("floor", "floor", floor_mat));
    }
    for (auto i = 0; i < shapes.size(); i++) {
        scn->instances.push_back(make_proc_instance(
            "obj" + std::to_string(i + 1), shapes[i], mats[i], frames[i]));
    }

    if (lights == "pointlights") {
        auto pos = std::vector<vec3f>{{-2, 10, 8}, {+2, 10, 8}};
        for (auto i = 0; i < 2; i++) {
            scn->instances.push_back(
                make_proc_instance("light" + std::to_string(i + 1), "point",
                    "emission", translation_frame(pos[i])));
            scn->instances.back()->mat->ke *= 120;
        }
    }
    if (lights == "arealights") {
        auto pos = std::vector<vec3f>{{0, 16, 0}, {0, 16, 16}};
        for (auto i = 0; i < 2; i++) {
            scn->instances.push_back(make_proc_instance(
                "light" + std::to_string(i + 1), "quad", "emission",
                lookat_frame(pos[i], {0, 1, 0}, {0, 0, 1}, true)));
            scn->instances.back()->mat->ke *= 8;
            for (auto& p : scn->instances.back()->shp->pos) p *= 8;
        }
    }
    if (lights == "arealights1") {
        auto pos = std::vector<vec3f>{{-4, 5, 8}, {+4, 5, 8}};
        for (auto i = 0; i < 2; i++) {
            scn->instances.push_back(make_proc_instance(
                "light" + std::to_string(i + 1), "quad", "emission",
                lookat_frame(pos[i], {0, 1, 0}, {0, 0, 1}, true)));
            scn->instances.back()->mat->ke *= 20;
            for (auto& p : scn->instances.back()->shp->pos) p *= 4;
        }
    }
    if (lights == "envlights") {
        scn->textures.push_back(make_proc_texture("sky", "sky", 512));
        scn->environments.push_back(
            make_environment("sky", {1, 1, 1}, scn->textures.back()));
    }
    for (auto ist : scn->instances) {
        scn->shapes.push_back(ist->shp);
        scn->materials.push_back(ist->mat);
    }
    if (nodes || !animations.empty()) {
        for (auto cam : scn->cameras) {
            scn->nodes.push_back(
                make_node(cam->name, cam, nullptr, nullptr, cam->frame));
        }
        for (auto ist : scn->instances) {
            scn->nodes.push_back(
                make_node(ist->name, nullptr, ist, nullptr, ist->frame));
        }
        for (auto env : scn->environments) {
            scn->nodes.push_back(
                make_node(env->name, nullptr, nullptr, env, env->frame));
        }
    }
    /*
    if (!animations.empty()) {
        for (auto i = 0; i < shapes.size(); i++) {
            auto name = "obj" + std::to_string(i + 1);
            auto nde = (node*)nullptr;
            for (auto n : scn->nodes)
                if (n->name == name) nde = n;
            nde->translation = nde->frame.o;
            nde->frame = identity_frame3f;
            scn->animations.push_back(new animation(*animations[i]));
            scn->animations.back()->name = name;
            scn->animations.back()->targets = {nde};
        }
    }
    */

    auto textures = std::map<std::string, texture*>();
    auto fix_texture = [&textures, scn](texture_info& tinfo) {
        if (!tinfo.txt) return;
        if (textures.find(tinfo.txt->name) != textures.end()) {
            delete tinfo.txt;
            tinfo.txt = textures.at(tinfo.txt->name);
        } else {
            scn->textures.push_back(tinfo.txt);
            textures[tinfo.txt->name] = tinfo.txt;
        }
    };
    for (auto mat : scn->materials) {
        fix_texture(mat->ke_txt);
        fix_texture(mat->kd_txt);
        fix_texture(mat->ks_txt);
        fix_texture(mat->norm_txt);
    }
    for (auto env : scn->environments) { fix_texture(env->ke_txt); }
    return scn;
}

// instances shared functions
scene* make_random_scene(const vec2i& num, const bbox3f& bbox, uint64_t seed) {
    auto rscale = 0.9f * 0.25f *
                  min((bbox.max.x - bbox.min.x) / num.x,
                      (bbox.max.x - bbox.min.x) / num.y);

    auto scn = new scene();
    scn->cameras.push_back(make_camera(
        "cam", {0, 6, 24}, {0, 1, 0}, 7.5f * pi / 180, 2.35f / 1.0f));
    {
        auto name = "floor"s;
        scn->shapes.push_back(make_proc_shape(name, "floor"));
        scn->materials.push_back(make_material(name, {0.2f, 0.2f, 0.2f}));
        scn->instances.push_back(make_instance(
            name, scn->shapes.back(), scn->materials.back(), identity_frame3f));
    }

    auto cols =
        std::vector<vec3f>{{0.5, 0.2, 0.2}, {0.2, 0.5, 0.2}, {0.2, 0.2, 0.5}};
    auto shps =
        std::vector<std::string>{"sphere", "sphere_flipcap", "sphere_flipcap"};
    auto shapes = std::vector<shape*>();
    auto materials = std::vector<material*>();
    auto sid = 0;
    for (auto shp : shps) {
        for (auto col : cols) {
            auto name = "shp" + std::to_string(sid++);
            scn->materials.push_back(make_material(name, col));
            materials.push_back(scn->materials.back());
            scn->shapes.push_back(make_proc_shape(
                name, shp, {32, 32, 32}, {2 * rscale, 2 * rscale, 2 * rscale}));
            shapes.push_back(scn->shapes.back());
        }
    }

    auto rng = make_rng(seed, 7);
    auto iid = 0;
    for (auto j = 0; j < num.y; j++) {
        for (auto i = 0; i < num.x; i++) {
            auto rpos = next_rand2f(rng);
            auto pos = vec3f{
                bbox.min.x + (bbox.max.x - bbox.min.x) *
                                 (i + 0.45f + 0.1f * rpos.x) / num.x,
                rscale,
                bbox.min.y + (bbox.max.y - bbox.min.y) *
                                 (j + 0.45f + 0.1f * rpos.y) / num.y,
            };
            auto idx = next_rand1i(rng, (int)shapes.size());
            scn->instances.push_back(
                make_instance("ist" + std::to_string(iid++), shapes[idx],
                    materials[idx], translation_frame(pos)));
        }
    }

    {
        auto pos = std::vector<vec3f>{{-2, 10, 8}, {+2, 10, 8}};
        for (auto i = 0; i < 2; i++) {
            auto name = "light" + std::to_string(i + 1);
            scn->shapes.push_back(make_proc_shape(name, "point"));
            scn->materials.push_back(
                make_proc_material(name, "emission", {80, 80, 80}));
            scn->instances.push_back(make_instance(name, scn->shapes.back(),
                scn->materials.back(), translation_frame(pos[i])));
        }
    }

    return scn;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PATH TRACING SUPPORT FUNCTION
// -----------------------------------------------------------------------------
namespace ygl {

// Phong exponent to roughness. Public API, see above.
float specular_exponent_to_roughness(float n) { return sqrtf(2 / (n + 2)); }

// Specular to fresnel eta. Public API, see above.
void specular_fresnel_from_ks(const vec3f& ks, vec3f& es, vec3f& esk) {
    es = {(1 + sqrt(ks.x)) / (1 - sqrt(ks.x)),
        (1 + sqrt(ks.y)) / (1 - sqrt(ks.y)),
        (1 + sqrt(ks.z)) / (1 - sqrt(ks.z))};
    esk = {0, 0, 0};
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

// Schlick approximation of Fresnel term
vec3f fresnel_schlick(const vec3f& ks, float cosw) {
    if (ks == zero3f) return zero3f;
    return ks +
           (vec3f{1, 1, 1} - ks) * pow(clamp(1.0f - cosw, 0.0f, 1.0f), 5.0f);
}

// Schlick approximation of Fresnel term weighted by roughness.
// This is a hack, but works better than not doing it.
vec3f fresnel_schlick(const vec3f& ks, float cosw, float rs) {
    if (ks == zero3f) return zero3f;
    auto fks = fresnel_schlick(ks, cosw);
    return ks + (fks - ks) * (1 - sqrt(clamp(rs, 0.0f, 1.0f)));
}

// Evaluates the GGX distribution and geometric term
float eval_ggx(float rs, float ndh, float ndi, float ndo) {
    // evaluate D
    auto alpha2 = rs * rs;
    auto di = (ndh * ndh) * (alpha2 - 1) + 1;
    auto d = alpha2 / (pi * di * di);
#if 0
    // evaluate G from Heitz
    auto lambda_o = (-1 + sqrt(1 + alpha2 * (1 - ndo * ndo) / (ndo * ndo))) / 2;
    auto lambda_i = (-1 + sqrt(1 + alpha2 * (1 - ndi * ndi) / (ndi * ndi))) / 2;
    auto g = 1 / (1 + lambda_o + lambda_i);
#else
    // evaluate G from Smith
    auto go = (2 * ndo) / (ndo + sqrt(alpha2 + (1 - alpha2) * ndo * ndo));
    auto gi = (2 * ndi) / (ndi + sqrt(alpha2 + (1 - alpha2) * ndi * ndi));
    auto g = go * gi;
#endif
    return d * g;
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

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR PATH TRACING
// -----------------------------------------------------------------------------
namespace ygl {

// Generates a 1-dimensional sample.
float sample_next1f(trace_pixel& pxl) {
    return clamp(next_rand1f(pxl.rng), 0.0f, 1 - flt_eps);
}

// Generates a 2-dimensional sample.
vec2f sample_next2f(trace_pixel& pxl) {
    return {next_rand1f(pxl.rng), next_rand1f(pxl.rng)};
}

// Type of point.
enum struct trace_point_type {
    none = 0,         // unitialized point
    point = 1,        // point
    curve = 2,        // curve
    surface = 3,      // surface
    environment = 4,  // environment
};

// Surface point with geometry and material data. Supports point on
// envmap too. This is the key data manipulated in the path tracer.
struct trace_point {
    const instance* ist = nullptr;                   // shape instance
    const environment* env = nullptr;                // environment
    trace_point_type type = trace_point_type::none;  // type
    vec3f pos = zero3f;                              // pos
    vec3f norm = {0, 0, 1};                          // norm
    vec2f texcoord = zero2f;                         // texcoord
    vec3f ke = zero3f;                               // emission
    vec3f kd = {0, 0, 0};                            // diffuse
    vec3f ks = {0, 0, 0};                            // specular
    float rs = 0;                                    // specular roughness
    vec3f kt = {0, 0, 0};                            // thin glass transmission
    float op = 1;                                    // opacity
    bool double_sided = false;                       // double sided
};

// Create a point for an environment map. Resolves material with textures.
trace_point eval_point(const environment* env, const vec2f& uv) {
    auto pt = trace_point();
    pt.env = env;
    pt.type = trace_point_type::environment;

    pt.pos = eval_pos(env, uv);
    pt.norm = eval_norm(env, uv);
    pt.texcoord = eval_texcoord(env, uv);

    pt.pos = transform_point(env->frame, pt.pos);
    pt.norm = transform_direction(env->frame, pt.norm);

    pt.ke = env->ke;
    if (env->ke_txt.txt) {
        auto txt = eval_texture(env->ke_txt, pt.texcoord);
        pt.ke *= {txt.x, txt.y, txt.z};
    }
    return pt;
}

// Create a point for a shape. Resolves geometry and material with textures.
trace_point eval_point(
    const instance* ist, int eid, const vec2f& euv, bool force_double_sided) {
    // default material
    static auto def_material = (material*)nullptr;
    if (!def_material) {
        def_material = new material();
        def_material->kd = {0.2f, 0.2f, 0.2f};
        def_material->rs = 1;
    }

    // point
    auto pt = trace_point();

    // shape
    pt.ist = ist;
    switch (get_shape_type(ist->shp)) {
        case shape_elem_type::none: {
            return pt;
        } break;
        case shape_elem_type::points:
        case shape_elem_type::vertices: {
            pt.type = trace_point_type::point;
        } break;
        case shape_elem_type::lines:
        case shape_elem_type::beziers: {
            pt.type = trace_point_type::curve;
        } break;
        case shape_elem_type::triangles:
        case shape_elem_type::quads:
        case shape_elem_type::facevarying: {
            pt.type = trace_point_type::surface;
        } break;
    }

    // shape values
    pt.pos = eval_pos(ist->shp, eid, euv);
    pt.norm = eval_norm(ist->shp, eid, euv);
    pt.texcoord = eval_texcoord(ist->shp, eid, euv);

    // shortcuts
    auto mat = (ist->mat) ? ist->mat : def_material;

    // handle normal map
    if (mat->norm_txt.txt) {
        auto tangsp = eval_tangsp(ist->shp, eid, euv);
        auto txt = eval_texture(mat->norm_txt, pt.texcoord, false) * 2.0f -
                   vec4f{1, 1, 1, 1};
        auto ntxt = normalize(vec3f{txt.x, -txt.y, txt.z});
        auto frame = make_frame_fromzx(
            {0, 0, 0}, pt.norm, {tangsp.x, tangsp.y, tangsp.z});
        frame.y *= tangsp.w;
        pt.norm = transform_direction(frame, ntxt);
    }

    // move to world coordinates
    pt.pos = transform_point(ist->frame, pt.pos);
    pt.norm = transform_direction(ist->frame, pt.norm);

    // double-sided
    pt.double_sided = mat->double_sided || force_double_sided;

    // initialized material values
    auto kx = vec3f{1, 1, 1};
    auto op = 1.0f;
    if (!ist->shp->color.empty()) {
        auto col = eval_color(ist->shp, eid, euv);
        kx *= {col.x, col.y, col.z};
        op *= col.w;
    }

    // handle occlusion
    if (mat->occ_txt.txt) {
        auto txt = eval_texture(mat->occ_txt, pt.texcoord);
        kx *= {txt.x, txt.y, txt.z};
    }

    // handle opacity
    pt.op = mat->op * op;
    if (mat->op_txt.txt) {
        auto txt = eval_texture(mat->op_txt, pt.texcoord);
        pt.op *= (txt.x + txt.y + txt.z) / 3;
    }

    // sample emission
    pt.ke = mat->ke * kx;
    if (mat->ke_txt.txt) {
        auto txt = eval_texture(mat->ke_txt, pt.texcoord);
        pt.ke *= {txt.x, txt.y, txt.z};
    }

    // sample reflectance
    switch (mat->type) {
        case material_type::specular_roughness: {
            pt.kd = mat->kd * kx;
            if (mat->kd_txt.txt) {
                auto txt = eval_texture(mat->kd_txt, pt.texcoord);
                pt.kd *= {txt.x, txt.y, txt.z};
                pt.op *= txt.w;
            }
            pt.ks = mat->ks * kx;
            pt.rs = mat->rs;
            if (mat->ks_txt.txt) {
                auto txt = eval_texture(mat->ks_txt, pt.texcoord);
                pt.ks *= {txt.x, txt.y, txt.z};
            }
            pt.kt = mat->kt * kx;
            if (mat->kt_txt.txt) {
                auto txt = eval_texture(mat->kt_txt, pt.texcoord);
                pt.kt *= {txt.x, txt.y, txt.z};
            }
        } break;
        case material_type::metallic_roughness: {
            auto kb = mat->kd * kx;
            if (mat->kd_txt.txt) {
                auto txt = eval_texture(mat->kd_txt, pt.texcoord);
                kb *= {txt.x, txt.y, txt.z};
                pt.op *= txt.w;
            }
            auto km = mat->ks.x;
            pt.rs = mat->rs;
            if (mat->ks_txt.txt) {
                auto txt = eval_texture(mat->ks_txt, pt.texcoord);
                km *= txt.y;
                pt.rs *= txt.z;
            }
            pt.kd = kb * (1 - km);
            pt.ks = kb * km + vec3f{0.04f, 0.04f, 0.04f} * (1 - km);
        } break;
        case material_type::specular_glossiness: {
            pt.kd = mat->kd * kx;
            if (mat->kd_txt.txt) {
                auto txt = eval_texture(mat->kd_txt, pt.texcoord);
                pt.kd *= {txt.x, txt.y, txt.z};
                pt.op *= txt.w;
            }
            pt.ks = mat->ks * kx;
            pt.rs = mat->rs;
            if (mat->ks_txt.txt) {
                auto txt = eval_texture(mat->ks_txt, pt.texcoord);
                pt.ks *= {txt.x, txt.y, txt.z};
                pt.rs *= txt.w;
            }
            pt.rs = 1 - pt.rs;  // glossiness -> roughnes
            pt.kt = mat->kt * kx;
            if (mat->kt_txt.txt) {
                auto txt = eval_texture(mat->kt_txt, pt.texcoord);
                pt.kt *= {txt.x, txt.y, txt.z};
            }
        } break;
    }

    // set up final values
    if (pt.ks != zero3f && pt.rs < 0.9999f) {
        pt.rs = pt.rs * pt.rs;
        if (pt.rs < 0.03f * 0.03f)
            pt.rs = 0;
        else
            pt.rs = clamp(pt.rs, 0.03f * 0.03f, 1.0f);
    } else {
        pt.ks = zero3f;
        pt.rs = 1;
    }
    if (pt.kt != zero3f) pt.double_sided = true;
    if (pt.op > 0.999f) pt.op = 1;

    // done
    return pt;
}

// Check if it has reflectance.
bool has_reflectance(const trace_point& pt) {
    return pt.kd + pt.ks + pt.kt != zero3f;
}

// Evaluates emission.
vec3f eval_emission(const trace_point& pt, const vec3f& wo) {
    if (pt.type != trace_point_type::surface || pt.double_sided ||
        dot(pt.norm, wo) >= 0)
        return pt.ke;
    return zero3f;
}

// Check if we are near the mirror direction.
inline bool check_near_mirror(
    const vec3f& wn, const vec3f& wo, const vec3f& wi) {
    return fabs(dot(wi, normalize(wn * 2.0f * dot(wo, wn) - wo)) - 1) < 0.001f;
}

// Evaluates the BRDF scaled by the cosine of the incoming direction.
// - ggx from [Heitz 2014] and [Walter 2007] and [Lagarde 2014]
// "Understanding the Masking-Shadowing Function in Microfacet-Based
// BRDFs" http://jcgt.org/published/0003/02/03/
// - "Microfacet Models for Refraction through Rough Surfaces" EGSR 07
// https://www.cs.cornell.edu/~srm/publications/EGSR07-btdf.pdf
vec3f eval_surface_brdfcos(
    const trace_point& pt, const vec3f& wo, const vec3f& wi) {
    auto brdfcos = zero3f;
    auto wn = pt.norm;
    if (pt.double_sided && dot(wo, wn) < 0) wn = -wn;

    auto ndo = dot(wn, wo), ndi = dot(wn, wi);

    if (pt.kd != zero3f && ndi > 0 && ndo > 0) { brdfcos += pt.kd * ndi / pi; }

    if (pt.ks != zero3f && pt.rs && ndi > 0 && ndo > 0) {
        auto wh = normalize(wo + wi);
        auto ndh = clamp(dot(wh, wn), -1.0f, 1.0f);
        auto dg = eval_ggx(pt.rs, ndh, ndi, ndo);
        auto odh = clamp(dot(wo, wh), 0.0f, 1.0f);
        auto ks = fresnel_schlick(pt.ks, odh, pt.rs);
        brdfcos += ks * ndi * dg / (4 * ndi * ndo);
    }

    if (pt.kt != zero3f && pt.rs && ndo > 0 && ndi < 0) {
        auto wir = wi - 2 * dot(wi, wn) * wn;
        auto wh = normalize(wo + wir);
        auto ndh = clamp(dot(wh, wn), -1.0f, 1.0f);
        auto dg = eval_ggx(pt.rs, ndh, -ndi, ndo);
        auto odh = clamp(dot(wo, wh), 0.0f, 1.0f);
        auto kt = pt.kt * (vec3f{1, 1, 1} - fresnel_schlick(pt.ks, odh, pt.rs));
        brdfcos += kt * ndi * dg / (4 * ndi * ndo);
    }

    assert(isfinite(brdfcos.x) && isfinite(brdfcos.y) && isfinite(brdfcos.z));
    return brdfcos;
}

// Evaluates the BRDF scaled by the cosine of the incoming direction.
vec3f eval_surface_delta_brdfcos(
    const trace_point& pt, const vec3f& wo, const vec3f& wi) {
    auto brdfcos = zero3f;
    auto wn = pt.norm;
    if (pt.double_sided && dot(wo, wn) < 0) wn = -wn;

    auto ndo = dot(wn, wo), ndi = dot(wn, wi);

    if (pt.ks != zero3f && !pt.rs && ndo > 0 && check_near_mirror(wn, wo, wi)) {
        auto ks = fresnel_schlick(pt.ks, ndo);
        brdfcos += ks;
    }

    if (pt.kt != zero3f && !pt.rs && wo == -wi) {
        auto kt = pt.kt * (vec3f{1, 1, 1} - fresnel_schlick(pt.ks, ndo));
        brdfcos += kt;
    }

    assert(isfinite(brdfcos.x) && isfinite(brdfcos.y) && isfinite(brdfcos.z));
    return brdfcos;
}

// Evaluates the BRDF scaled by the cosine of the incoming direction.
// - uses Kajiya-Kay for hair
vec3f eval_curve_brdfcos(
    const trace_point& pt, const vec3f& wo, const vec3f& wi) {
    auto brdfcos = zero3f;
    auto wh = normalize(wo + wi);
    auto wn = pt.norm;

    auto ndo = dot(wn, wo), ndi = dot(wn, wi),
         ndh = clamp(dot(wn, wh), 0.0f, 1.0f);
    auto so = sqrt(clamp(1 - ndo * ndo, 0.0001f, 1.0f)),
         si = sqrt(clamp(1 - ndi * ndi, 0.0001f, 1.0f)),
         sh = sqrt(clamp(1 - ndh * ndh, 0.0001f, 1.0f));

    if (pt.kd != zero3f) { brdfcos += pt.kd * si / pi; }

    if (pt.ks != zero3f && pt.rs) {
        auto ns = 2 / (pt.rs * pt.rs) - 2;
        auto d = (ns + 2) * pow(sh, ns) / (2 + pi);
        brdfcos += pt.ks * si * d / (4.0f * si * so);
    }

    assert(isfinite(brdfcos.x) && isfinite(brdfcos.y) && isfinite(brdfcos.z));
    return brdfcos;
}

// Evaluates the BRDF scaled by the cosine of the incoming direction.
// - uses a hack for points
vec3f eval_point_brdfcos(
    const trace_point& pt, const vec3f& wo, const vec3f& wi) {
    auto brdfcos = zero3f;

    auto ido = dot(wo, wi);
    brdfcos += pt.kd * (2 * ido + 1) / (2 * pi);

    assert(isfinite(brdfcos.x) && isfinite(brdfcos.y) && isfinite(brdfcos.z));
    return brdfcos;
}

// Evaluates the BRDF scaled by the cosine of the incoming direction.
vec3f eval_brdfcos(const trace_point& pt, const vec3f& wo, const vec3f& wi,
    bool delta = false) {
    if (!has_reflectance(pt)) return zero3f;
    switch (pt.type) {
        case trace_point_type::none: return zero3f;
        case trace_point_type::surface:
            return (!delta) ? eval_surface_brdfcos(pt, wo, wi) :
                              eval_surface_delta_brdfcos(pt, wo, wi);
        case trace_point_type::curve: return eval_curve_brdfcos(pt, wo, wi);
        case trace_point_type::point: return eval_point_brdfcos(pt, wo, wi);
        case trace_point_type::environment: return zero3f;
    }
}

// Compute the weight for sampling the BRDF
float weight_surface_brdf(
    const trace_point& pt, const vec3f& wo, const vec3f& wi) {
    auto prob_kd = max(pt.kd), prob_ks = max(pt.ks), prob_kt = max(pt.kt);
    auto prob_sum = prob_kd + prob_ks + prob_kt;
    if (prob_sum == 0) return 0;
    prob_kd /= prob_sum;
    prob_ks /= prob_sum;
    prob_kt /= prob_sum;

    // normal
    auto wn = pt.norm;
    if (pt.double_sided && dot(wo, wn) < 0) wn = -wn;

    auto ndo = dot(wn, wo), ndi = dot(wn, wi);

    auto pdf = 0.0f;

    if (prob_kd && ndo > 0 && ndi > 0) { pdf += prob_kd * ndi / pi; }

    if (prob_ks && pt.rs && ndo > 0 && ndi > 0) {
        auto wh = normalize(wi + wo);
        auto ndh = dot(wn, wh);
        auto d = sample_ggx_pdf(pt.rs, ndh);
        auto hdo = dot(wo, wh);
        pdf += prob_ks * d / (4 * hdo);
    }

    if (prob_kt && pt.rs && ndo > 0 && ndi < 0) {
        auto wir = wi - 2 * dot(wi, wn) * wn;
        auto wh = normalize(wo + wir);
        auto ndh = dot(wn, wh);
        auto d = sample_ggx_pdf(pt.rs, ndh);
        auto hdo = dot(wo, wh);
        pdf += prob_kt * d / (4 * hdo);
    }

    assert(isfinite(pdf));
    if (!pdf) return 0;
    return 1 / pdf;
}

// Compute the weight for sampling the BRDF
float weight_surface_delta_brdf(
    const trace_point& pt, const vec3f& wo, const vec3f& wi) {
    if (pt.rs) return 0;
    auto prob_ks = max(pt.ks), prob_kt = max(pt.kt);
    auto prob_sum = prob_ks + prob_kt;
    if (prob_sum == 0) return 0;
    prob_ks /= prob_sum;
    prob_kt /= prob_sum;

    // normal
    auto wn = pt.norm;
    if (pt.double_sided && dot(wo, wn) < 0) wn = -wn;

    auto ndo = dot(wn, wo), ndi = dot(wn, wi);

    auto pdf = 0.0f;

    if (prob_ks && !pt.rs && ndo > 0 && check_near_mirror(wn, wo, wi)) {
        pdf += prob_ks;
    }

    if (prob_kt && !pt.rs && wo == -wi) { pdf += prob_kt; }

    assert(isfinite(pdf));
    if (!pdf) return 0;
    return 1 / pdf;
}

// Compute the weight for sampling the BRDF
float weight_curve_brdf(
    const trace_point& pt, const vec3f& wo, const vec3f& wi) {
    auto pdf = 1 / (4 * pi);
    return 1 / pdf;
}

// Compute the weight for sampling the BRDF
float weight_point_brdf(
    const trace_point& pt, const vec3f& wo, const vec3f& wi) {
    auto pdf = 1 / (4 * pi);
    return 1 / pdf;
}

// Compute the weight for sampling the BRDF
float weight_brdf(const trace_point& pt, const vec3f& wo, const vec3f& wi,
    bool delta = false) {
    if (!has_reflectance(pt)) return 0;
    switch (pt.type) {
        case trace_point_type::none: return 0;
        case trace_point_type::surface:
            return (!delta) ? weight_surface_brdf(pt, wo, wi) :
                              weight_surface_delta_brdf(pt, wo, wi);
        case trace_point_type::curve: return weight_curve_brdf(pt, wo, wi);
        case trace_point_type::point: return weight_point_brdf(pt, wo, wi);
        case trace_point_type::environment: return 0;
    }
}

// Picks a direction based on the BRDF
vec3f sample_surface_brdf(
    const trace_point& pt, const vec3f& wo, float rnl, const vec2f& rn) {
    auto prob_kd = max(pt.kd), prob_ks = (pt.rs) ? max(pt.ks) : 0,
         prob_kt = (pt.rs) ? max(pt.kt) : 0;
    auto prob_sum = prob_kd + prob_ks + prob_kt;
    if (prob_sum == 0) return zero3f;
    prob_kd /= prob_sum;
    prob_ks /= prob_sum;
    prob_kt /= prob_sum;

    auto wn = pt.norm;
    if (pt.double_sided && dot(wo, wn) < 0) wn = -wn;

    auto ndo = dot(wn, wo);
    if (ndo <= 0) return zero3f;

    // sample according to diffuse
    if (rnl < prob_kd) {
        auto fp = make_frame_fromz(pt.pos, wn);
        auto rz = sqrtf(rn.y), rr = sqrtf(1 - rz * rz), rphi = 2 * pi * rn.x;
        auto wi_local = vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
        return transform_direction(fp, wi_local);
    }
    // sample according to specular GGX
    else if (rnl < prob_kd + prob_ks && pt.rs) {
        auto fp = make_frame_fromz(pt.pos, wn);
        auto wh_local = sample_ggx(pt.rs, rn);
        auto wh = transform_direction(fp, wh_local);
        return normalize(wh * 2.0f * dot(wo, wh) - wo);
    }
    // transmission hack
    else if (rnl < prob_kd + prob_ks + prob_kt && pt.rs) {
        auto fp = make_frame_fromz(pt.pos, wn);
        auto wh_local = sample_ggx(pt.rs, rn);
        auto wh = transform_direction(fp, wh_local);
        auto wi = normalize(wh * 2.0f * dot(wo, wh) - wo);
        return normalize(wi - 2 * dot(wi, wn) * wn);
    } else {
        return zero3f;
    }
}

// Picks a direction based on the BRDF
vec3f sample_surface_delta_brdf(
    const trace_point& pt, const vec3f& wo, float rnl, const vec2f& rn) {
    if (pt.rs) return zero3f;
    auto prob_ks = max(pt.ks), prob_kt = max(pt.kt);
    auto prob_sum = prob_ks + prob_kt;
    if (prob_sum == 0) return zero3f;
    prob_ks /= prob_sum;
    prob_kt /= prob_sum;

    auto wn = pt.norm;
    if (pt.double_sided && dot(wo, wn) < 0) wn = -wn;

    auto ndo = dot(wn, wo);
    if (ndo <= 0) return zero3f;

    // sample according to specular mirror
    if (rnl < prob_ks && !pt.rs) {
        return normalize(wn * 2.0f * dot(wo, wn) - wo);
    }
    // transmission hack
    else if (rnl < prob_ks + prob_kt && !pt.rs) {
        return -wo;
    } else {
        return zero3f;
    }
}

// Picks a direction based on the BRDF
vec3f sample_curve_brdf(
    const trace_point& pt, const vec3f& wo, float rnl, const vec2f& rn) {
    auto wn = pt.norm;
    auto fp = make_frame_fromz(pt.pos, wn);
    auto rz = 2 * rn.y - 1, rr = sqrtf(1 - rz * rz), rphi = 2 * pi * rn.x;
    auto wi_local = vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
    return transform_direction(fp, wi_local);
}

// Picks a direction based on the BRDF
vec3f sample_point_brdf(
    const trace_point& pt, const vec3f& wo, float rnl, const vec2f& rn) {
    auto wn = pt.norm;
    auto fp = make_frame_fromz(pt.pos, wn);
    auto rz = 2 * rn.y - 1, rr = sqrtf(1 - rz * rz), rphi = 2 * pi * rn.x;
    auto wi_local = vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
    return transform_direction(fp, wi_local);
}

// Picks a direction based on the BRDF
vec3f sample_brdf(const trace_point& pt, const vec3f& wo, float rnl,
    const vec2f& rn, bool delta = false) {
    if (!has_reflectance(pt)) return zero3f;
    switch (pt.type) {
        case trace_point_type::none: return zero3f;
        case trace_point_type::surface:
            return (!delta) ? sample_surface_brdf(pt, wo, rnl, rn) :
                              sample_surface_delta_brdf(pt, wo, rnl, rn);
        case trace_point_type::curve: return sample_curve_brdf(pt, wo, rnl, rn);
        case trace_point_type::point: return sample_point_brdf(pt, wo, rnl, rn);
        case trace_point_type::environment: return zero3f;
    }
}

float sample_delta_prob(const trace_point& pt, const vec3f& wo) {
    if (pt.type != trace_point_type::surface) return 0;
    if (pt.rs) return 0;
    auto dw = max(pt.ks) + max(pt.kt);
    return dw / (dw + max(pt.kd));
}

// Sample weight for a light point.
float weight_light(
    const trace_lights& lights, const trace_point& lpt, const trace_point& pt) {
    switch (lpt.type) {
        case trace_point_type::none: {
            throw std::runtime_error("should not have gotten here");
        } break;
        case trace_point_type::point: {
            auto area = lights.shape_distribs.at(lpt.ist->shp).back();
            auto dist = length(lpt.pos - pt.pos);
            return area / (dist * dist);
        } break;
        case trace_point_type::curve: {
            throw std::runtime_error("not implemented yet");
        } break;
        case trace_point_type::surface: {
            auto area = lights.shape_distribs.at(lpt.ist->shp).back();
            auto dist = length(lpt.pos - pt.pos);
            return area * fabs(dot(lpt.norm, normalize(lpt.pos - pt.pos))) /
                   (dist * dist);
        } break;
        case trace_point_type::environment: {
            return 4 * pi;
        } break;
    }
    return 0;
}

// Sample weight for a light point.
float weight_lights(
    const trace_lights& lights, const trace_point& lpt, const trace_point& pt) {
    return lights.lights.size() * weight_light(lights, lpt, pt);
}

// Picks a point on a light.
trace_point sample_light(const trace_lights& lights, const trace_light& lgt,
    const trace_point& pt, float rel, const vec2f& ruv,
    const trace_params& params) {
    if (lgt.ist) {
        auto& dst = lights.shape_distribs.at(lgt.ist->shp);
        auto sample = sample_shape(lgt.ist->shp, dst, rel, ruv);
        return eval_point(
            lgt.ist, sample.first, sample.second, params.double_sided);
    } else if (lgt.env) {
        // BUG: this is not uniform sampling
        return eval_point(lgt.env, ruv);
    } else {
        throw std::runtime_error("should not have gotten here");
    }
    return {};
}

// Picks a point on a light.
trace_point sample_lights(const trace_lights& lights, const trace_point& pt,
    float rnl, float rne, const vec2f& ruv, const trace_params& params) {
    auto lidx = sample_discrete(lights.light_distrib, rnl);
    auto& lgt = lights.lights.at(lidx);
    return sample_light(lights, lgt, pt, rne, ruv, params);
}

// Intersects a ray with the scn and return the point (or env point).
trace_point intersect_scene(const scene* scn, const bvh_tree* bvh,
    const ray3f& ray, const trace_params& params) {
    auto iid = 0, eid = 0;
    auto euv = zero2f;
    auto ray_t = 0.0f;
    if (intersect_bvh(bvh, ray, false, ray_t, iid, eid, euv)) {
        return eval_point(scn->instances[iid], eid, euv, params.double_sided);
    } else if (!scn->environments.empty()) {
        return eval_point(
            scn->environments[0], eval_uv(scn->environments[0],
                                      transform_direction_inverse(
                                          scn->environments[0]->frame, ray.d)));
    } else {
        return {};
    }
}

// Test occlusion.
vec3f eval_transmission(const scene* scn, const bvh_tree* bvh,
    const trace_point& pt, const trace_point& lpt, const trace_params& params) {
    if (params.notransmission) {
        auto ray = make_segment(pt.pos, lpt.pos);
        return (intersect_bvh(bvh, ray, true)) ? zero3f : vec3f{1, 1, 1};
    } else {
        auto cpt = pt;
        auto weight = vec3f{1, 1, 1};
        for (auto bounce = 0; bounce < params.max_depth; bounce++) {
            auto ray = make_segment(cpt.pos, lpt.pos);
            cpt = intersect_scene(scn, bvh, ray, params);
            if (!cpt.ist) break;
            weight *= cpt.kt + vec3f{1 - cpt.op, 1 - cpt.op, 1 - cpt.op};
            if (weight == zero3f) break;
        }
        return weight;
    }
}

// Mis weight.
float weight_mis(float w0, float w1) {
    if (w0 == 0) return w1;
    if (w1 == 0) return w0;
    return 1 / (1 / w0 + 1 / w1);
}

// Recursive path tracing.
vec3f trace_path(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt_, const vec3f& wo_,
    rng_state& rng, const trace_params& params) {
    if (lights.lights.empty()) return zero3f;
    auto pt = pt_;
    auto wo = wo_;

    // initialize
    auto l = eval_emission(pt, wo);
    auto weight = vec3f{1, 1, 1};

    // trace  path
    for (auto bounce = 0; bounce < params.max_depth; bounce++) {
        // opacity
        if (pt.op != 1) {
            if (next_rand1f(rng) < 1 - pt.op) {
                pt = intersect_scene(scn, bvh, make_ray(pt.pos, -wo), params);
                l += weight * eval_emission(pt, wo);
                continue;
            }
        }

        // early exit
        if (!has_reflectance(pt)) break;

        // direct – light
        auto& lgt = lights.lights[sample_discrete(
            lights.light_distrib, next_rand1f(rng))];
        auto lpt = sample_light(
            lights, lgt, pt, next_rand1f(rng), next_rand2f(rng), params);
        auto lw = weight_light(lights, lpt, pt) * lights.light_distrib.back();
        auto lwi = normalize(lpt.pos - pt.pos);
        auto lke = eval_emission(lpt, -lwi);
        auto lbc = eval_brdfcos(pt, wo, lwi);
        auto lld = lke * lbc;
        if (lld != zero3f) {
            l += weight * lld * eval_transmission(scn, bvh, pt, lpt, params) *
                 weight_mis(lw, weight_brdf(pt, wo, lwi));
        }

        // direct – brdf
        auto bwi = sample_brdf(pt, wo, next_rand1f(rng), next_rand2f(rng));
        auto bpt = intersect_scene(scn, bvh, make_ray(pt.pos, bwi), params);
        auto bw = weight_brdf(pt, wo, bwi);
        auto bke = eval_emission(bpt, -bwi);
        auto bbc = eval_brdfcos(pt, wo, bwi);
        auto bld = bke * bbc;
        if (bld != zero3f) {
            l += weight * bld *
                 weight_mis(bw, weight_light(lights, bpt, pt) *
                                    lights.light_distrib.back());
        }

        // skip recursion if path ends
        if (bounce == params.max_depth - 1) break;

        // deltas
        auto delta_prob = sample_delta_prob(pt, wo);
        auto delta = (delta_prob && next_rand1f(rng) < delta_prob);
        weight *= (delta) ? 1 / delta_prob : 1 / (1 - delta_prob);
        if (delta) {
            bwi =
                sample_brdf(pt, wo, next_rand1f(rng), next_rand2f(rng), delta);
            bpt = intersect_scene(scn, bvh, make_ray(pt.pos, bwi), params);
            l += weight * eval_emission(bpt, -bwi) *
                 eval_brdfcos(pt, wo, bwi, delta) *
                 weight_brdf(pt, wo, bwi, delta);
        }

        // continue path
        weight *=
            eval_brdfcos(pt, wo, bwi, delta) * weight_brdf(pt, wo, bwi, delta);
        if (weight == zero3f) break;

        // roussian roulette
        if (bounce > 2) {
            auto rrprob = 1.0f - min(max(weight), 0.95f);
            if (next_rand1f(rng) < rrprob) break;
            weight *= 1 / (1 - rrprob);
        }

        // continue path
        pt = bpt;
        wo = -bwi;
    }

    return l;
}

// Recursive path tracing.
vec3f trace_path_nomis(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt_, const vec3f& wo_,
    rng_state& rng, const trace_params& params) {
    if (lights.lights.empty()) return zero3f;
    auto pt = pt_;
    auto wo = wo_;

    // initialize
    auto l = eval_emission(pt, wo);
    auto weight = vec3f{1, 1, 1};

    // trace  path
    for (auto bounce = 0; bounce < params.max_depth; bounce++) {
        // opacity
        if (pt.op != 1) {
            if (next_rand1f(rng) < 1 - pt.op) {
                pt = intersect_scene(scn, bvh, make_ray(pt.pos, -wo), params);
                l += weight * eval_emission(pt, wo);
                continue;
            }
        }

        // early exit
        if (!has_reflectance(pt)) break;

        // direct
        if (!lights.lights.empty()) {
            auto& lgt = lights.lights[sample_discrete(
                lights.light_distrib, next_rand1f(rng))];
            auto lpt = sample_light(
                lights, lgt, pt, next_rand1f(rng), next_rand2f(rng), params);
            auto lwi = normalize(lpt.pos - pt.pos);
            auto ld = eval_emission(lpt, -lwi) * eval_brdfcos(pt, wo, lwi) *
                      weight_light(lights, lpt, pt) *
                      lights.light_distrib.back();
            if (ld != zero3f) {
                l += weight * ld * eval_transmission(scn, bvh, pt, lpt, params);
            }
        }

        // skip recursion if path ends
        if (bounce == params.max_depth - 1) break;

        // choose delta
        auto delta_prob = sample_delta_prob(pt, wo);
        auto delta = (delta_prob && next_rand1f(rng) < delta_prob);
        weight *= (delta) ? 1 / delta_prob : 1 / (1 - delta_prob);

        // continue path
        auto bwi =
            sample_brdf(pt, wo, next_rand1f(rng), next_rand2f(rng), delta);
        weight *=
            eval_brdfcos(pt, wo, bwi, delta) * weight_brdf(pt, wo, bwi, delta);
        if (weight == zero3f) break;

        // roussian roulette
        if (bounce > 2) {
            auto rrprob = 1.0f - min(max(weight), 0.95f);
            if (next_rand1f(rng) < rrprob) break;
            weight *= 1 / (1 - rrprob);
        }

        pt = intersect_scene(scn, bvh, make_ray(pt.pos, bwi), params);
        wo = -bwi;
        if (delta) l += weight * eval_emission(pt, wo);
    }

    return l;
}

// Direct illumination.
vec3f trace_direct(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt, const vec3f& wo,
    int bounce, rng_state& rng, const trace_params& params) {
    // emission
    auto l = eval_emission(pt, wo);

    // ambient
    l += params.ambient * pt.kd;

    // direct
    for (auto& lgt : lights.lights) {
        auto lpt = sample_light(
            lights, lgt, pt, next_rand1f(rng), next_rand2f(rng), params);
        auto lwi = normalize(lpt.pos - pt.pos);
        auto ld = eval_emission(lpt, -lwi) * eval_brdfcos(pt, wo, lwi) *
                  weight_light(lights, lpt, pt);
        if (ld == zero3f) continue;
        l += ld * eval_transmission(scn, bvh, pt, lpt, params);
    }

    // exit if needed
    if (bounce >= params.max_depth) return l;

    // reflection
    if (pt.ks != zero3f && !pt.rs) {
        auto wi = reflect(wo, pt.norm);
        auto rpt = intersect_scene(scn, bvh, make_ray(pt.pos, wi), params);
        l += pt.ks *
             trace_direct(scn, bvh, lights, rpt, -wi, bounce + 1, rng, params);
    }

    // opacity
    if (pt.kt != zero3f) {
        auto opt = intersect_scene(scn, bvh, make_ray(pt.pos, -wo), params);
        l += pt.kt *
             trace_direct(scn, bvh, lights, opt, wo, bounce + 1, rng, params);
    }

    // opacity
    if (pt.op != 1) {
        auto opt = intersect_scene(scn, bvh, make_ray(pt.pos, -wo), params);
        l = pt.op * l + (1 - pt.op) * trace_direct(scn, bvh, lights, opt, wo,
                                          bounce + 1, rng, params);
    }

    // done
    return l;
}

// Direct illumination.
vec3f trace_direct(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt, const vec3f& wo,
    rng_state& rng, const trace_params& params) {
    return trace_direct(scn, bvh, lights, pt, wo, 0, rng, params);
}

// Eyelight for quick previewing.
vec3f trace_eyelight(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt, const vec3f& wo,
    int bounce, rng_state& rng, const trace_params& params) {
    // emission
    auto l = eval_emission(pt, wo);

    // brdf*light
    l += eval_brdfcos(pt, wo, wo) * pi;

    // opacity
    if (bounce >= params.max_depth) return l;
    if (pt.kt != zero3f) {
        auto opt = intersect_scene(scn, bvh, make_ray(pt.pos, -wo), params);
        l += pt.kt *
             trace_eyelight(scn, bvh, lights, opt, wo, bounce + 1, rng, params);
    }
    if (pt.op != 1) {
        auto opt = intersect_scene(scn, bvh, make_ray(pt.pos, -wo), params);
        l = pt.op * l + (1 - pt.op) * trace_eyelight(scn, bvh, lights, opt, wo,
                                          bounce + 1, rng, params);
    }

    // done
    return l;
}

// Eyelight for quick previewing.
vec3f trace_eyelight(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt, const vec3f& wo,
    rng_state& rng, const trace_params& params) {
    return trace_eyelight(scn, bvh, lights, pt, wo, 0, rng, params);
}

// Debug previewing.
vec3f trace_debug_normal(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt, const vec3f& wo,
    rng_state& rng, const trace_params& params) {
    auto wn = pt.norm;
    if (pt.double_sided && dot(wn, wo) < 0) wn = -wn;
    return wn * 0.5f + vec3f{0.5f, 0.5f, 0.5f};
}

// Debug frontfacing.
vec3f trace_debug_frontfacing(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt, const vec3f& wo,
    rng_state& rng, const trace_params& params) {
    auto wn = pt.norm;
    if (pt.double_sided && dot(wn, wo) < 0) wn = -wn;
    return dot(wn, wo) > 0 ? vec3f{0, 1, 0} : vec3f{1, 0, 0};
}

// Debug previewing.
vec3f trace_debug_albedo(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt, const vec3f& wo,
    rng_state& rng, const trace_params& params) {
    return pt.kd + pt.ks + pt.kt;
}

// Debug previewing.
vec3f trace_debug_texcoord(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt, const vec3f& wo,
    rng_state& rng, const trace_params& params) {
    return {pt.texcoord.x, pt.texcoord.y, 0};
}

// Trace function
vec3f trace_func(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt, const vec3f& wo,
    rng_state& rng, const trace_params& params) {
    switch (params.tracer) {
        case trace_type::eyelight:
            return trace_eyelight(scn, bvh, lights, pt, wo, rng, params);
        case trace_type::direct:
            return trace_direct(scn, bvh, lights, pt, wo, rng, params);
        case trace_type::pathtrace:
            return trace_path(scn, bvh, lights, pt, wo, rng, params);
        case trace_type::pathtrace_nomis:
            return trace_path_nomis(scn, bvh, lights, pt, wo, rng, params);
        case trace_type::debug_albedo:
            return trace_debug_albedo(scn, bvh, lights, pt, wo, rng, params);
        case trace_type::debug_normal:
            return trace_debug_normal(scn, bvh, lights, pt, wo, rng, params);
        case trace_type::debug_frontfacing:
            return trace_debug_frontfacing(
                scn, bvh, lights, pt, wo, rng, params);
        case trace_type::debug_texcoord:
            return trace_debug_texcoord(scn, bvh, lights, pt, wo, rng, params);
    }
}

// Trace a single sample
void trace_sample(const scene* scn, const camera* cam, const bvh_tree* bvh,
    const trace_lights& lights, trace_pixel& pxl, const trace_params& params) {
    auto crn = next_rand2f(pxl.rng);
    auto lrn = next_rand2f(pxl.rng);
    auto uv = vec2f{(pxl.i + crn.x) / (cam->aspect * params.resolution),
        1 - (pxl.j + crn.y) / params.resolution};
    auto ray = eval_camera_ray(cam, uv, lrn);
    auto pt = intersect_scene(scn, bvh, ray, params);
    if (!pt.ist && params.envmap_invisible) return;
    auto l = trace_func(scn, bvh, lights, pt, -ray.d, pxl.rng, params);
    if (!isfinite(l.x) || !isfinite(l.y) || !isfinite(l.z)) {
        log_error("NaN detected");
        return;
    }
    if (max(l) > params.pixel_clamp) l = l * (params.pixel_clamp / max(l));
    pxl.acc += {l.x, l.y, l.z, 1};
    pxl.sample += 1;
}

// Trace the next nsamples.
void trace_samples(const scene* scn, const camera* cam, const bvh_tree* bvh,
    const trace_lights& lights, image4f& img, std::vector<trace_pixel>& pixels,
    int nsamples, const trace_params& params) {
    if (params.parallel) {
        auto nthreads = std::thread::hardware_concurrency();
        auto threads = std::vector<std::thread>();
        for (auto tid = 0; tid < std::thread::hardware_concurrency(); tid++) {
            threads.push_back(std::thread([=, &img, &pixels, &params]() {
                for (auto pid = tid; pid < pixels.size(); pid += nthreads) {
                    auto& pxl = pixels.at(pid);
                    for (auto s = 0; s < nsamples; s++)
                        trace_sample(scn, cam, bvh, lights, pxl, params);
                    img.at(pxl.i, pxl.j) = pxl.acc / pxl.sample;
                }
            }));
        }
        for (auto& t : threads) t.join();
        threads.clear();
    } else {
        for (auto& pxl : pixels) {
            for (auto s = 0; s < params.nsamples; s++)
                trace_sample(scn, cam, bvh, lights, pxl, params);
            img.at(pxl.i, pxl.j) = pxl.acc / pxl.sample;
        }
    }
}

// Starts an anyncrhounous renderer.
void trace_async_start(const scene* scn, const camera* cam, const bvh_tree* bvh,
    const trace_lights& lights, image4f& img, std::vector<trace_pixel>& pixels,
    std::vector<std::thread>& threads, bool& stop_flag,
    const trace_params& params, const std::function<void(int, int)>& callback) {
    pixels = make_trace_pixels(img, params);

    // render preview
    if (params.preview_resolution) {
        auto pparams = params;
        pparams.resolution = params.preview_resolution;
        pparams.nsamples = 1;
        auto pimg =
            make_image4f((int)std::round(cam->aspect * pparams.resolution),
                pparams.resolution);
        auto ppixels = make_trace_pixels(pimg, pparams);
        trace_samples(scn, cam, bvh, lights, pimg, ppixels, 1, pparams);
        resize_image(pimg, img, ygl::resize_filter::box);
    } else {
        for (auto& p : img.pixels) p = zero4f;
    }
    if (callback) callback(0, 0);

    // start rendering
    auto nthreads = std::thread::hardware_concurrency();
    for (auto tid = 0; tid < std::thread::hardware_concurrency(); tid++) {
        threads.push_back(std::thread([=, &img, &pixels, &stop_flag]() {
            for (auto s = 0; s < params.nsamples; s++) {
                for (auto pid = tid; pid < pixels.size(); pid += nthreads) {
                    if (stop_flag) return;
                    auto& pxl = pixels.at(pid);
                    trace_sample(scn, cam, bvh, lights, pxl, params);
                    img.at(pxl.i, pxl.j) = pxl.acc / pxl.sample;
                    if (!pxl.i && callback) callback(s, pxl.j);
                }
            }
            if (!tid && callback) callback(params.nsamples, 0);
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

// Initialize trace lights
trace_lights make_trace_lights(const scene* scn) {
    auto lights = trace_lights();

    for (auto ist : scn->instances) {
        if (!ist->mat || ist->mat->ke == zero3f) continue;
        auto lgt = trace_light();
        lgt.ist = ist;
        lights.lights.push_back(lgt);
        if (lights.shape_distribs.find(ist->shp) == lights.shape_distribs.end())
            lights.shape_distribs[ist->shp] = sample_shape_cdf(ist->shp);
    }

    for (auto env : scn->environments) {
        if (env->ke == zero3f) continue;
        auto lgt = trace_light();
        lgt.env = env;
        lights.lights.push_back(lgt);
    }

    lights.light_distrib = std::vector<float>(lights.lights.size());
    for (auto i = 0; i < lights.lights.size(); i++)
        lights.light_distrib[i] = i + 1;

    return lights;
}

// Initialize a rendering state
std::vector<trace_pixel> make_trace_pixels(
    const image4f& img, const trace_params& params) {
    auto pixels = std::vector<trace_pixel>(img.width * img.height);
    for (auto j = 0; j < img.height; j++) {
        for (auto i = 0; i < img.width; i++) {
            auto pxl = trace_pixel();
            pxl.i = i;
            pxl.j = j;
            pxl.rng = make_rng(params.seed, (j * img.width + i) * 2 + 1);
            pixels[j * img.width + i] = pxl;
        }
    }
    return pixels;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR WAVEFRONT OBJ
// -----------------------------------------------------------------------------
namespace ygl {

// skip whitespace
inline void obj_skipws(char*& s) {
    while (*s == ' ') s++;
}

// skip a string if matched
inline bool obj_streq(const char* s, const char* str) {
    while (*s == *str && *s && *str) {
        s++;
        str++;
    }
    return *s == *str;
}

#if YGL_FASTOBJ

// parse base value
inline void obj_parse_base(char*& s, int& val) {
    val = 0;
    auto sn = (*s == '-') ? -1 : 1;
    if (*s == '-' || *s == '+') s++;
    while (*s >= '0' && *s <= '9') val = val * 10 + (*s++ - '0');
    val *= sn;
}

// parse base value
inline void obj_parse_base(char*& s, float& val) {
    //    auto ss = s; auto sss = ss;
    auto mantissa = 0, fractional = 0, fractional_length = 0, exponent = 0;
    auto sn = (*s == '-') ? -1 : 1;
    if (*s == '-' || *s == '+') s++;
    while (*s >= '0' && *s <= '9') mantissa = mantissa * 10 + (*s++ - '0');
    if (*s == '.') {
        s++;
        while (*s >= '0' && *s <= '9') {
            fractional = fractional * 10 + (*s++ - '0');
            fractional_length++;
        }
    }
    mantissa *= sn;
    fractional *= sn;
    if (*s == 'e' || *s == 'E') {
        s++;
        auto en = (*s == '-') ? -1 : 1;
        if (*s == '-' || *s == '+') s++;
        while (*s >= '0' && *s <= '9') exponent = exponent * 10 + (*s++ - '0');
        exponent *= en;
    }
    auto dval = (double)mantissa;
    if (fractional)
        dval += fractional * std::pow(10.0, -(double)fractional_length);
    if (exponent) dval *= std::pow(10.0, (double)exponent);
    val = (float)dval;
#if 0
    auto cval = val;
    sscanf(ss, "%f", &cval);
    if(abs(val - cval) > 0.01f) {
        printf("- %g %g %s\n", val, cval, sss);
    }
    auto len = 0;
    sscanf(s, "%f%n", &val, &len);
    s += len;
#endif
}

// parse base value
inline void obj_parse_base(char*& s, char* val) {
    while (*s && *s != ' ') *val++ = *s++;
    *val = 0;
}

// parse base value
inline void obj_parse_base(char*& s, std::string& val) {
    char buf[4096];
    obj_parse_base(s, buf);
    val = buf;
}

#else

// parse base value
inline void obj_parse_base(char*& s, int& val) {
    auto len = 0;
    sscanf(s, "%d%n", &val, &len);
    s += len;
}

// parse base value
inline void obj_parse_base(char*& s, float& val) {
    auto len = 0;
    sscanf(s, "%f%n", &val, &len);
    s += len;
}

// parse base value
inline void obj_parse_base(char*& s, std::string& val) {
    char buf[4096];
    auto len = 0;
    sscanf(s, "%s%n", buf, &len);
    if (len) {
        s += len;
        val = buf;
    } else {
        val = "";
    }
}

// parse base value
inline void obj_parse_base(char*& s, char* val) {
    auto len = 0;
    sscanf(s, "%s%n", val, &len);
    s += len;
}

#endif

// parse value
inline void obj_parse(char*& s, int& val) {
    obj_skipws(s);
    obj_parse_base(s, val);
}
inline void obj_parse(char*& s, float& val) {
    obj_skipws(s);
    obj_parse_base(s, val);
}
inline void obj_parse(char*& s, bool& val) {
    auto i = 0;
    obj_parse(s, i);
    val = i;
}
inline void obj_parse(char*& s, std::string& val) {
    obj_skipws(s);
    obj_parse_base(s, val);
}
inline void obj_parse(char*& s, char* val) {
    obj_skipws(s);
    obj_parse_base(s, val);
}
inline void obj_parse(char*& s, vec2f& val) {
    for (auto i = 0; i < 2; i++) obj_parse(s, (&val.x)[i]);
}
inline void obj_parse(char*& s, vec3f& val) {
    for (auto i = 0; i < 3; i++) obj_parse(s, (&val.x)[i]);
}
inline void obj_parse(char*& s, vec4f& val) {
    for (auto i = 0; i < 4; i++) obj_parse(s, (&val.x)[i]);
}
inline void obj_parse(char*& s, frame3f& val) {
    for (auto i = 0; i < 12; i++) obj_parse(s, (&val.x.x)[i]);
}
inline void obj_parse(char*& s, obj_vertex& val, const obj_vertex& vert_size) {
    char buf[1024];
    obj_skipws(s);
    obj_parse_base(s, buf);
    val = obj_vertex{-1, -1, -1, -1, -1};
    auto v = &val.pos;
    auto vs = &vert_size.pos;
    auto i = 0;
    auto sb = buf;
    while (i < 5 && *sb) {
        obj_parse_base(sb, v[i]);
        v[i] = (v[i] < 0) ? vs[i] + v[i] : v[i] - 1;
        if (*sb != '/') break;
        while (*sb == '/') {
            sb++;
            i++;
        }
    }
}

// clear the whitespace
inline void obj_convertws(char* s) {
    while (*s) {
        if (*s == '\t' || *s == '\r' || *s == '\n') *s = ' ';
        s++;
    }
}

// Parse texture options and name
inline void obj_parse(char*& s, obj_texture_info& info) {
    // initialize
    info = obj_texture_info();

    // get tokens
    auto tokens = std::vector<std::string>();
    obj_skipws(s);
    while (*s) {
        tokens.push_back("");
        obj_parse(s, tokens.back());
        obj_skipws(s);
    }

    // exit if no tokens
    if (tokens.empty()) return;

    // texture name
    info.path = tokens.back();
    for (auto& c : info.path)
        if (c == '\\') c = '/';

    // texture options
    auto last = std::string();
    for (auto& tok : tokens) {
        if (tok == tokens.back()) break;
        if (tok[0] == '-') {
            last = tok;
            info.props[last] = {};
        } else {
            info.props[last].push_back(tok);
        }
    }

    // clamp
    if (info.props.find("-clamp") != info.props.end() &&
        info.props.at("-clamp").size() > 0) {
        auto& clamp_vec = info.props.at("-clamp");
        auto clamp_str = (clamp_vec.empty()) ? "" : clamp_vec.front();
        info.clamp = clamp_str == "on" || clamp_str == "1";
        info.props.erase("-clamp");
    }

    if (info.props.find("-bm") != info.props.end() &&
        info.props.at("-bm").size() > 0) {
        auto& bm_vec = info.props.at("-bm");
        auto bm_str = (bm_vec.empty()) ? "" : bm_vec.front();
        info.scale = std::atof(bm_str.c_str());
        info.props.erase("-bm");
    }
}

// Load MTL
std::vector<obj_material*> load_mtl(const std::string& filename, bool flip_tr,
    std::vector<std::string>& textures) {
    // clear materials
    auto materials = std::vector<obj_material*>();

    // open file
    auto fs = fopen(filename.c_str(), "rt");
    if (!fs) throw std::runtime_error("cannot open filename " + filename);

    // add a material preemptively to avoid crashes
    materials.push_back(new obj_material());
    auto mat = materials.back();

    // read the file line by line
    char line[4096];
    char cmd[1024];
    auto linenum = 0;
    while (fgets(line, sizeof(line), fs)) {
        // prepare to parse
        linenum += 1;
        auto ss = line;
        obj_convertws(ss);
        obj_skipws(ss);

        // skip empty and comments
        if (!ss[0] || ss[0] == '#') continue;

        // get command
        obj_parse(ss, cmd);

        // possible token values
        if (obj_streq(cmd, "newmtl")) {
            materials.push_back(new obj_material());
            mat = materials.back();
            obj_parse(ss, mat->name);
        } else if (obj_streq(cmd, "illum")) {
            obj_parse(ss, mat->illum);
        } else if (obj_streq(cmd, "Ke")) {
            obj_parse(ss, mat->ke);
        } else if (obj_streq(cmd, "Ka")) {
            obj_parse(ss, mat->ka);
        } else if (obj_streq(cmd, "Kd")) {
            obj_parse(ss, mat->kd);
        } else if (obj_streq(cmd, "Ks")) {
            obj_parse(ss, mat->ks);
        } else if (obj_streq(cmd, "Kr")) {
            obj_parse(ss, mat->kr);
        } else if (obj_streq(cmd, "Kt")) {
            obj_parse(ss, mat->kt);
        } else if (obj_streq(cmd, "Tf")) {
            auto nchan = 0;
            obj_skipws(ss);
            while (*ss && nchan < 3) {
                obj_parse(ss, (&mat->kt.x)[nchan++]);
                obj_skipws(ss);
            }
            if (nchan < 3) mat->kt = {mat->kt.x, mat->kt.x, mat->kt.x};
            if (flip_tr)
                materials.back()->kt = vec3f{1, 1, 1} - materials.back()->kt;
        } else if (obj_streq(cmd, "Tr")) {
            auto nchan = 0;
            auto tr = zero3f;
            obj_skipws(ss);
            while (*ss && nchan < 3) {
                obj_parse(ss, (&tr.x)[nchan++]);
                obj_skipws(ss);
            }
            if (nchan < 3) tr = {tr.x, tr.x, tr.x};
            materials.back()->op = (tr.x + tr.y + tr.z) / 3;
            if (flip_tr) materials.back()->op = 1 - materials.back()->op;
        } else if (obj_streq(cmd, "Ns")) {
            obj_parse(ss, mat->ns);
        } else if (obj_streq(cmd, "d")) {
            obj_parse(ss, mat->op);
        } else if (obj_streq(cmd, "Ni")) {
            obj_parse(ss, mat->ior);
        } else if (obj_streq(cmd, "map_Ke")) {
            obj_parse(ss, mat->ke_txt);
        } else if (obj_streq(cmd, "map_Ka")) {
            obj_parse(ss, mat->ka_txt);
        } else if (obj_streq(cmd, "map_Kd")) {
            obj_parse(ss, mat->kd_txt);
        } else if (obj_streq(cmd, "map_Ks")) {
            obj_parse(ss, mat->ks_txt);
        } else if (obj_streq(cmd, "map_Kr")) {
            obj_parse(ss, mat->kr_txt);
        } else if (obj_streq(cmd, "map_Tr")) {
            obj_parse(ss, mat->kt_txt);
        } else if (obj_streq(cmd, "map_Ns")) {
            obj_parse(ss, mat->ns_txt);
        } else if (obj_streq(cmd, "map_d")) {
            obj_parse(ss, mat->op_txt);
        } else if (obj_streq(cmd, "map_Ni")) {
            obj_parse(ss, mat->ior_txt);
        } else if (obj_streq(cmd, "map_bump") || obj_streq(cmd, "bump")) {
            obj_parse(ss, mat->bump_txt);
        } else if (obj_streq(cmd, "map_disp") || obj_streq(cmd, "disp")) {
            obj_parse(ss, mat->disp_txt);
        } else if (obj_streq(cmd, "map_norm") || obj_streq(cmd, "norm")) {
            obj_parse(ss, mat->norm_txt);
        } else {
            // copy into strings
            obj_skipws(ss);
            while (*ss) {
                mat->props[cmd].push_back("");
                obj_parse(ss, mat->props[cmd].back());
                obj_skipws(ss);
            }
        }
    }

    // remove first fake material
    materials.erase(materials.begin());

    // create texture array
    textures = {};
    auto texture_set = std::unordered_set<std::string>();
    auto add_texture = [&texture_set, &textures](const obj_texture_info& info) {
        if (info.path == "") return;
        if (texture_set.find(info.path) != texture_set.end()) return;
        texture_set.insert(info.path);
        textures.push_back(info.path);
    };
    for (auto mat : materials) {
        add_texture(mat->ke_txt);
        add_texture(mat->ka_txt);
        add_texture(mat->kd_txt);
        add_texture(mat->ks_txt);
        add_texture(mat->kr_txt);
        add_texture(mat->kt_txt);
        add_texture(mat->ns_txt);
        add_texture(mat->op_txt);
        add_texture(mat->ior_txt);
        add_texture(mat->bump_txt);
        add_texture(mat->bump_txt);
        add_texture(mat->disp_txt);
        add_texture(mat->norm_txt);
    }

    // clone
    fclose(fs);

    // done
    return materials;
}

// Loads textures for an scene.
void load_textures(
    const obj_scene* obj, const std::string& dirname, bool skip_missing) {
    for (auto txt : obj->textures) {
        auto filename = dirname + txt->path;
        for (auto& c : filename)
            if (c == '\\') c = '/';
#if YGL_IMAGEIO
        if (is_hdr_filename(filename)) {
            txt->dataf =
                load_imagef(filename, txt->width, txt->height, txt->ncomp);
        } else {
            txt->datab =
                load_imageb(filename, txt->width, txt->height, txt->ncomp);
        }
#endif
        if (txt->datab.empty() && txt->dataf.empty()) {
            if (skip_missing) continue;
            throw std::runtime_error("cannot laod image " + filename);
        }
    }
}

// Loads an OBJ
obj_scene* load_obj(const std::string& filename, bool split_shapes,
    bool load_txt, bool skip_missing, bool flip_texcoord, bool flip_tr) {
    // clear obj
    auto obj = new obj_scene();

    // splitting policy
    auto split_material = split_shapes;
    auto split_group = split_shapes;
    auto split_smoothing = split_shapes;

    // open file
    auto fs = fopen(filename.c_str(), "rt");
    if (!fs) throw std::runtime_error("cannot open filename " + filename);

    // current parsing values
    auto mtllibs = std::vector<std::string>();
    auto oobj = (obj_object*)nullptr;
    auto matname = ""s;
    auto oname = ""s;
    auto faceted = false;
    auto elems = std::vector<obj_vertex>();

    // initializing obj
    obj->objects.push_back(new obj_object());
    oobj = obj->objects.back();
    oobj->groups.push_back({"", "", false});

    // keep track of array lengths
    auto vert_size = obj_vertex{0, 0, 0, 0, 0};

    // elem type map
    static auto elem_type_map =
        std::unordered_map<std::string, obj_element_type>{
            {"f", obj_element_type::face}, {"l", obj_element_type::line},
            {"p", obj_element_type::point}, {"b", obj_element_type::bezier}};

    // read the file line by line
    char line[4096];
    char cmd[1024];
    auto linenum = 0;
    while (fgets(line, sizeof(line), fs)) {
        // prepare to parse
        linenum += 1;
        auto ss = line;
        obj_convertws(ss);
        obj_skipws(ss);

        // skip empty and comments
        if (!ss[0] || ss[0] == '#') continue;

        // get command
        obj_parse(ss, cmd);

        // possible token values
        if (obj_streq(cmd, "v")) {
            vert_size.pos += 1;
            obj->pos.push_back(zero3f);
            obj_parse(ss, obj->pos.back());
        } else if (obj_streq(cmd, "vn")) {
            vert_size.norm += 1;
            obj->norm.push_back(zero3f);
            obj_parse(ss, obj->norm.back());
        } else if (obj_streq(cmd, "vt")) {
            vert_size.texcoord += 1;
            obj->texcoord.push_back(zero2f);
            obj_parse(ss, obj->texcoord.back());
            if (flip_texcoord)
                obj->texcoord.back().y = 1 - obj->texcoord.back().y;
        } else if (obj_streq(cmd, "vc")) {
            vert_size.color += 1;
            obj->color.push_back(vec4f{0, 0, 0, 1});
            obj_parse(ss, obj->color.back());
        } else if (obj_streq(cmd, "vr")) {
            vert_size.radius += 1;
            obj->radius.push_back(0);
            obj_parse(ss, obj->radius.back());
        } else if (obj_streq(cmd, "f") || obj_streq(cmd, "l") ||
                   obj_streq(cmd, "p") || obj_streq(cmd, "b")) {
            auto elem = obj_element();
            elem.type = elem_type_map.at(cmd);
            elem.start = (uint32_t)oobj->verts.size();
            elem.size = 0;
            elem.groupid = (int)oobj->groups.size() - 1;
            oobj->elems.push_back(elem);
            obj_skipws(ss);
            while (*ss) {
                auto vert = obj_vertex();
                obj_parse(ss, vert, vert_size);
                obj_skipws(ss);
                oobj->verts.push_back(vert);
                oobj->elems.back().size += 1;
            }
        } else if (obj_streq(cmd, "o")) {
            auto name = ""s;
            obj_parse(ss, name);
            obj->objects.push_back(new obj_object());
            oobj = obj->objects.back();
            oobj->name = name;
            oobj->groups.push_back({"", "", false});
            matname = "";
            oname = name;
        } else if (obj_streq(cmd, "usemtl")) {
            obj_parse(ss, matname);
            if (split_material) {
                obj->objects.push_back(new obj_object());
                oobj = obj->objects.back();
                oobj->name = oname;
                oobj->groups.push_back({"", matname, faceted});
            } else {
                oobj->groups.push_back({oobj->groups.back().name, matname,
                    oobj->groups.back().faceted});
            }
        } else if (obj_streq(cmd, "g")) {
            auto name = ""s;
            obj_parse(ss, name);
            if (split_group) {
                obj->objects.push_back(new obj_object());
                oobj = obj->objects.back();
                oobj->name = oname + name;
                oobj->groups.push_back({name, matname, faceted});
            } else {
                oobj->groups.push_back({name, oobj->groups.back().matname,
                    oobj->groups.back().faceted});
            }
        } else if (obj_streq(cmd, "s")) {
            auto name = std::string();
            obj_parse(ss, name);
            faceted = (name != "on");
            if (split_smoothing) {
                obj->objects.push_back(new obj_object());
                oobj = obj->objects.back();
                oobj->name = oname;
                oobj->groups.push_back({"", matname, faceted});
            } else {
                oobj->groups.push_back({oobj->groups.back().name,
                    oobj->groups.back().matname, faceted});
            }
        } else if (obj_streq(cmd, "op")) {
            auto name = std::string();
            obj_parse(ss, name);
            obj_skipws(ss);
            while (*ss) {
                auto tok = std::string();
                obj_parse(ss, tok);
                obj_skipws(ss);
                oobj->props[name].push_back(tok);
            }
        } else if (obj_streq(cmd, "mtllib")) {
            mtllibs.push_back("");
            obj_parse(ss, mtllibs.back());
        } else if (obj_streq(cmd, "c")) {
            auto cam = new obj_camera();
            obj_parse(ss, cam->name);
            obj_parse(ss, cam->ortho);
            obj_parse(ss, cam->yfov);
            obj_parse(ss, cam->aspect);
            obj_parse(ss, cam->aperture);
            obj_parse(ss, cam->focus);
            obj_parse(ss, cam->frame);
            obj->cameras.push_back(cam);
        } else if (obj_streq(cmd, "e")) {
            auto env = new obj_environment();
            obj_parse(ss, env->name);
            obj_parse(ss, env->ke);
            obj_parse(ss, env->ke_txt.path);
            if (env->ke_txt.path == "\"\"") env->ke_txt.path = "";
            obj_parse(ss, env->frame);
            obj->environments.push_back(env);
        } else if (obj_streq(cmd, "n")) {
            auto nde = new obj_node();
            obj_parse(ss, nde->name);
            obj_parse(ss, nde->parent);
            obj_parse(ss, nde->camname);
            obj_parse(ss, nde->objname);
            obj_parse(ss, nde->envname);
            obj_parse(ss, nde->frame);
            obj_parse(ss, nde->translation);
            obj_parse(ss, nde->rotation);
            obj_parse(ss, nde->scale);
            if (nde->parent == "\"\"") nde->parent = "";
            if (nde->camname == "\"\"") nde->camname = "";
            if (nde->objname == "\"\"") nde->objname = "";
            if (nde->envname == "\"\"") nde->envname = "";
            obj->nodes.push_back(nde);
        } else {
            // unused
        }
    }

    // cleanup empty
    for (auto idx = 0; idx < obj->objects.size(); idx++) {
        if (!obj->objects[idx]->elems.empty()) continue;
        if (!obj->objects[idx]->verts.empty()) continue;
        delete obj->objects[idx];
        obj->objects.erase(obj->objects.begin() + idx);
        idx--;
    }

    // cleanup unused
    for (auto oobj : obj->objects) {
        if (oobj->groups.empty()) continue;
        auto used = std::vector<bool>(oobj->groups.size());
        for (auto& elem : oobj->elems) used[elem.groupid] = true;
        auto emap = std::vector<uint16_t>(oobj->groups.size());
        auto valid = 0;
        for (auto idx = 0; idx < oobj->groups.size(); idx++) {
            emap[idx] = valid;
            if (used[idx]) valid++;
        }
        for (auto& elem : oobj->elems) elem.groupid = emap[elem.groupid];
        for (auto idx = 0; idx < used.size(); idx++) {
            if (used[idx]) continue;
            used.erase(used.begin() + idx);
            oobj->groups.erase(oobj->groups.begin() + idx);
            idx--;
        }
        if (oobj->groups.size() == 1 && oobj->groups[0].name == "" &&
            oobj->groups[0].matname == "" && oobj->groups[0].faceted == false)
            oobj->groups.clear();
    }

    auto end = std::remove_if(obj->objects.begin(), obj->objects.end(),
        [](const obj_object* x) { return !x; });
    obj->objects.erase(end, obj->objects.end());

    // parse materials
    auto mtllibs_set =
        std::unordered_set<std::string>(mtllibs.begin(), mtllibs.end());
    mtllibs = std::vector<std::string>{mtllibs_set.begin(), mtllibs_set.end()};
    auto dirname = path_dirname(filename);
    std::unordered_set<std::string> texture_set;
    for (auto mtllib : mtllibs) {
        auto mtlname = dirname + mtllib;
        std::vector<std::string> textures;
        auto materials = load_mtl(mtlname, flip_tr, textures);
        obj->materials.insert(
            obj->materials.end(), materials.begin(), materials.end());
        for (auto& txt : textures) {
            if (texture_set.find(txt) != texture_set.end()) continue;
            obj->textures.push_back(new obj_texture());
            obj->textures.back()->path = txt;
            texture_set.insert(txt);
        }
    }
    for (auto oenv : obj->environments) {
        auto txt = oenv->ke_txt.path;
        if (txt == "") continue;
        if (texture_set.find(txt) != texture_set.end()) continue;
        obj->textures.push_back(new obj_texture());
        obj->textures.back()->path = txt;
        texture_set.insert(txt);
    }

    // load textures
    if (load_txt) load_textures(obj, dirname, skip_missing);

    // close file
    fclose(fs);

    // done
    return obj;
}

// Dumps a value
inline void obj_dump(char*& s, char* val) {
    while (*val) *s++ = *val++;
}
inline void obj_dump(char*& s, const char* val) {
    while (*val) *s++ = *val++;
}
inline void obj_dump(char*& s, const std::string& val) {
    auto val_ = val.c_str();
    while (*val_) *s++ = *val_++;
}
inline void obj_dump(char*& s, int val) { s += sprintf(s, "%d", val); }
inline void obj_dump(char*& s, float val) { s += sprintf(s, "%g", val); }
inline void obj_dump(char*& s, const float* val, int num) {
    for (auto i = 0; i < num; i++) {
        if (i) *s++ = ' ';
        obj_dump(s, val[i]);
    }
}
inline void obj_dump(char*& s, const vec2f& val) { obj_dump(s, &val.x, 2); }
inline void obj_dump(char*& s, const vec3f& val) { obj_dump(s, &val.x, 3); }
inline void obj_dump(char*& s, const vec4f& val) { obj_dump(s, &val.x, 4); }
inline void obj_dump(char*& s, const frame3f& val) {
    obj_dump(s, &val.x.x, 12);
}
inline void obj_dump(char*& s, const obj_vertex& val) {
    auto vert_ptr = &val.pos;
    auto nto_write = 0;
    for (auto i = 0; i < 5; i++) {
        if (vert_ptr[i] >= 0) nto_write = i + 1;
    }
    for (auto i = 0; i < nto_write; i++) {
        if (i) *s++ = '/';
        if (vert_ptr[i] >= 0) s += sprintf(s, "%d", vert_ptr[i] + 1);
    }
}
inline void obj_dump(char*& s, const std::vector<std::string>& vals) {
    for (auto i = 0; i < vals.size(); i++) {
        if (i) *s++ = ' ';
        obj_dump(s, vals[i]);
    }
}
inline void obj_dump(char*& s, const obj_texture_info& v) {
    for (auto&& kv : v.props) {
        obj_dump(s, kv.first);
        *s++ = ' ';
        for (auto&& vv : kv.second) {
            obj_dump(s, vv);
            *s++ = ' ';
        }
    }
    if (v.clamp) obj_dump(s, "-clamp on ");
    obj_dump(s, v.path);
}

// Dumps a line
template <typename T>
inline void obj_dump_line(FILE* fs, const char* lbl, const T& val) {
    char buf[4096];
    buf[0] = 0;
    auto s = buf;
    obj_dump(s, lbl);
    *s++ = ' ';
    obj_dump(s, val);
    *s++ = '\n';
    *s = 0;
    fputs(buf, fs);
}

// Dumps a line
inline void obj_dump_line(
    FILE* fs, const char* lbl, const obj_vertex* vals, int count) {
    char buf[4096];
    buf[0] = 0;
    auto s = buf;
    obj_dump(s, lbl);
    for (auto i = 0; i < count; i++) {
        *s++ = ' ';
        obj_dump(s, vals[i]);
    }
    *s++ = '\n';
    *s = 0;
    fputs(buf, fs);
}

// Dumps a line
template <typename T>
inline void obj_dump_sp(FILE* fs, const T& val) {
    char buf[4096];
    buf[0] = 0;
    auto s = buf;
    obj_dump(s, val);
    *s++ = ' ';
    *s = 0;
    fputs(buf, fs);
}

// Save an MTL file
void save_mtl(const std::string& filename,
    const std::vector<obj_material*>& materials, bool flip_tr) {
    // open file
    auto fs = fopen(filename.c_str(), "wt");
    if (!fs) throw std::runtime_error("cannot open filename " + filename);

    // for each material, dump all the values
    for (auto mat : materials) {
        obj_dump_line(fs, "newmtl", mat->name);
        obj_dump_line(fs, "  illum", mat->illum);
        if (mat->ke != zero3f) obj_dump_line(fs, "  Ke", mat->ke);
        if (mat->ka != zero3f) obj_dump_line(fs, "  Ka", mat->ka);
        if (mat->kd != zero3f) obj_dump_line(fs, "  Kd", mat->kd);
        if (mat->ks != zero3f) obj_dump_line(fs, "  Ks", mat->ks);
        if (mat->kr != zero3f) obj_dump_line(fs, "  Kr", mat->kr);
        if (mat->kt != zero3f) obj_dump_line(fs, "  Kt", mat->kt);
        // if (mat->kt != zero3f) obj_dump_line(fs, "  Tf", mat->kt);
        if (mat->ns != 0.0f)
            obj_dump_line(fs, "  Ns", (int)clamp(mat->ns, 0.0f, 1000000000.0f));
        if (mat->op != 1.0f) obj_dump_line(fs, "  d", mat->op);
        if (mat->ior != 1.0f) obj_dump_line(fs, "  Ni", mat->ior);
        if (mat->ke_txt.path != "") obj_dump_line(fs, "  map_Ke", mat->ke_txt);
        if (mat->ka_txt.path != "") obj_dump_line(fs, "  map_Ka", mat->ka_txt);
        if (mat->kd_txt.path != "") obj_dump_line(fs, "  map_Kd", mat->kd_txt);
        if (mat->ks_txt.path != "") obj_dump_line(fs, "  map_Ks", mat->ks_txt);
        if (mat->kr_txt.path != "") obj_dump_line(fs, "  map_Kr", mat->kr_txt);
        if (mat->kt_txt.path != "") obj_dump_line(fs, "  map_Kt", mat->kt_txt);
        if (mat->ns_txt.path != "") obj_dump_line(fs, "  map_Ns", mat->ns_txt);
        if (mat->op_txt.path != "") obj_dump_line(fs, "  map_d ", mat->op_txt);
        if (mat->ior_txt.path != "")
            obj_dump_line(fs, "  map_Ni", mat->ior_txt);
        if (mat->bump_txt.path != "")
            obj_dump_line(fs, "  map_bump", mat->bump_txt);
        if (mat->disp_txt.path != "")
            obj_dump_line(fs, "  map_disp", mat->disp_txt);
        if (mat->norm_txt.path != "")
            obj_dump_line(fs, "  map_norm", mat->norm_txt);
        for (auto&& kv : mat->props) {
            obj_dump_line(fs, ("  " + kv.first + ' ').c_str(), kv.second);
        }
        fputs("\n", fs);
    }

    fclose(fs);
}

// Loads textures for an scene.
void save_textures(
    const obj_scene* obj, const std::string& dirname, bool skip_missing) {
    for (auto txt : obj->textures) {
        if (txt->datab.empty() && txt->dataf.empty()) continue;
        auto filename = dirname + txt->path;
        for (auto& c : filename)
            if (c == '\\') c = '/';
        auto ok = false;
#if YGL_IMAGEIO
        if (!txt->datab.empty()) {
            ok = save_imageb(filename, txt->width, txt->height, txt->ncomp,
                txt->datab.data());
        }
        if (!txt->dataf.empty()) {
            ok = save_imagef(filename, txt->width, txt->height, txt->ncomp,
                txt->dataf.data());
        }
#endif
        if (!ok) {
            if (skip_missing) continue;
            throw std::runtime_error("cannot save image " + filename);
        }
    }
}

// Save an OBJ
void save_obj(const std::string& filename, const obj_scene* obj, bool save_txt,
    bool skip_missing, bool flip_texcoord, bool flip_tr) {
    // open file
    auto fs = fopen(filename.c_str(), "wt");
    if (!fs) throw std::runtime_error("cannot open filename " + filename);

    // linkup to mtl
    auto dirname = path_dirname(filename);
    auto basename = filename.substr(dirname.length());
    basename = basename.substr(0, basename.length() - 4);
    if (!obj->materials.empty()) {
        obj_dump_line(fs, "mtllib", basename + ".mtl");
    }

    // save cameras
    for (auto cam : obj->cameras) {
        obj_dump_sp(fs, "c");
        obj_dump_sp(fs, cam->name);
        obj_dump_sp(fs, cam->ortho);
        obj_dump_sp(fs, cam->yfov);
        obj_dump_sp(fs, cam->aspect);
        obj_dump_sp(fs, cam->aperture);
        obj_dump_sp(fs, cam->focus);
        obj_dump_sp(fs, cam->frame);
        obj_dump_sp(fs, "\n");
    }

    // save envs
    for (auto env : obj->environments) {
        obj_dump_sp(fs, "e");
        obj_dump_sp(fs, env->name);
        obj_dump_sp(fs, env->ke);
        obj_dump_sp(fs, (env->ke_txt.path != "") ? env->ke_txt.path : "\"\""s);
        obj_dump_sp(fs, env->frame);
        obj_dump_sp(fs, "\n");
    }

    // save nodes
    for (auto nde : obj->nodes) {
        obj_dump_sp(fs, "n");
        obj_dump_sp(fs, nde->name);
        obj_dump_sp(fs, (nde->parent.empty()) ? "\"\""s : nde->parent);
        obj_dump_sp(fs, (nde->camname.empty()) ? "\"\""s : nde->camname);
        obj_dump_sp(fs, (nde->objname.empty()) ? "\"\""s : nde->objname);
        obj_dump_sp(fs, (nde->envname.empty()) ? "\"\""s : nde->envname);
        obj_dump_sp(fs, nde->frame);
        obj_dump_sp(fs, nde->translation);
        obj_dump_sp(fs, nde->rotation);
        obj_dump_sp(fs, nde->scale);
        obj_dump_sp(fs, "\n");
    }

    // save all vertex data
    for (auto& v : obj->pos) obj_dump_line(fs, "v", v);
    if (flip_texcoord) {
        for (auto& v : obj->texcoord)
            obj_dump_line(fs, "vt", vec2f{v.x, 1 - v.y});
    } else {
        for (auto& v : obj->texcoord) obj_dump_line(fs, "vt", v);
    }
    for (auto& v : obj->norm) obj_dump_line(fs, "vn", v);
    for (auto& v : obj->color) obj_dump_line(fs, "vc", v);
    for (auto& v : obj->radius) obj_dump_line(fs, "vr", v);

    // save element data
    static auto elem_labels = std::unordered_map<obj_element_type, std::string>{
        {obj_element_type::point, "p"}, {obj_element_type::line, "l"},
        {obj_element_type::face, "f"}, {obj_element_type::bezier, "b"}};
    for (auto oobj : obj->objects) {
        obj_dump_line(fs, "o", oobj->name);
        for (auto& kv : oobj->props) {
            auto nv = kv.second;
            nv.insert(nv.begin(), kv.first);
            obj_dump_line(fs, "op", nv);
        }
        if (!oobj->groups.empty()) {
            auto& ogrp = oobj->groups[0];
            if (ogrp.name != "") obj_dump_line(fs, "g", ogrp.name);
            if (ogrp.matname != "") obj_dump_line(fs, "usemtl", ogrp.matname);
            if (ogrp.faceted) obj_dump_line(fs, "s", "off");
        }
        auto last_groupid = 0;
        for (auto& elem : oobj->elems) {
            if (last_groupid != elem.groupid) {
                auto& lgrp = oobj->groups[last_groupid];
                auto& ogrp = oobj->groups[elem.groupid];
                if (ogrp.name != lgrp.name) obj_dump_line(fs, "g", ogrp.name);
                if (ogrp.matname != "" && ogrp.matname != lgrp.matname)
                    obj_dump_line(fs, "usemtl", ogrp.matname);
                if (ogrp.faceted != lgrp.faceted)
                    obj_dump_line(fs, "s", ogrp.faceted ? "off" : "on");
                last_groupid = elem.groupid;
            }
            obj_dump_line(fs, elem_labels.at(elem.type).c_str(),
                oobj->verts.data() + elem.start, elem.size);
        }
    }

    fclose(fs);

    // save materials
    if (!obj->materials.empty())
        save_mtl(dirname + basename + ".mtl", obj->materials, flip_tr);

    // save textures
    if (save_txt) save_textures(obj, dirname, skip_missing);
}

}  // namespace ygl

#if YGL_GLTF

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR KHRONOS GLTF
// -----------------------------------------------------------------------------
namespace ygl {

// Json alias
using json = nlohmann::json;

// Parse int function.
void serialize(int& val, json& js, bool reading) {
    if (reading) {
        if (!js.is_number_integer())
            throw std::runtime_error("integer expected");
        val = js;
    } else {
        js = val;
    }
}

// Parse float function.
void serialize(float& val, json& js, bool reading) {
    if (reading) {
        if (!js.is_number()) throw std::runtime_error("number expected");
        val = js;
    } else {
        js = val;
    }
}

// Parse bool function.
void serialize(bool& val, json& js, bool reading) {
    if (reading) {
        if (!js.is_boolean()) throw std::runtime_error("bool expected");
        val = js;
    } else {
        js = val;
    }
}

// Parse string function.
void serialize(std::string& val, json& js, bool reading) {
    if (reading) {
        if (!js.is_string()) throw std::runtime_error("string expected");
        val = js;
    } else {
        js = val;
    }
}

// Parse json function.
void serialize(json& val, json& js, bool reading) {
    if (reading) {
        val = js;
    } else {
        js = val;
    }
}

// Parse support function.
template <typename T>
void serialize(T*& val, json& js, bool reading) {
    if (reading) {
        if (js.is_null()) {
            if (val) delete val;
            val = nullptr;
            return;
        }
        if (!js.is_object()) throw std::runtime_error("object expected");
        if (!val) val = new T();
        serialize(*val, js, reading);
    } else {
        if (!val) {
            js = nullptr;
            return;
        }
        if (!js.is_object()) js = json::object();
        serialize(*val, js, reading);
    }
}

// Parse support function.
template <typename T>
void serialize(std::vector<T>& vals, json& js, bool reading) {
    if (reading) {
        if (!js.is_array()) throw std::runtime_error("array expected");
        vals.resize(js.size());
        for (auto i = 0; i < js.size(); i++) {
            // this is contrived to support for std::vector<bool>
            auto v = T();
            serialize(v, js[i], reading);
            vals[i] = v;
        }
    } else {
        js = json::array();
        for (auto i = 0; i < vals.size(); i++)
            serialize(vals[i], js[i], reading);
    }
}

// Parse support function.
template <typename T, size_t N>
void serialize(std::array<T, N>& vals, json& js, bool reading) {
    if (reading) {
        if (!js.is_array()) throw std::runtime_error("array expected");
        if (N != js.size()) throw std::runtime_error("wrong array size");
        for (auto i = 0; i < N; i++) serialize(vals[i], js.at(i), reading);
    } else {
        js = json::array();
        for (auto i = 0; i < N; i++) serialize(vals[i], js[i], reading);
    }
}

// Parse support function.
template <typename T>
void serialize(std::map<std::string, T>& vals, json& js, bool reading) {
    if (reading) {
        if (!js.is_object()) throw std::runtime_error("object expected");
        for (auto kv = js.begin(); kv != js.end(); ++kv) {
            serialize(vals[kv.key()], kv.value(), reading);
        }
    } else {
        js = json::object();
        for (auto& kv : vals) serialize(kv.second, js[kv.first], reading);
    }
}

// Parse support function.
template <typename T, typename T1>
void serialize(T& val, json& js, bool reading,
    const std::vector<std::pair<T1, T>>& table) {
    if (reading) {
        auto v = T1();
        serialize(v, js, reading);
        auto found = false;
        for (auto& kv : table) {
            if (kv.first == v) {
                val = kv.second;
                found = true;
                break;
            }
        }
        if (!found) throw std::runtime_error("bad enum value");
    } else {
        auto found = false;
        auto v = T1();
        for (auto& kv : table) {
            if (kv.second == val) {
                v = kv.first;
                found = true;
                break;
            }
        }
        if (!found) throw std::runtime_error("invalid value");
        serialize(v, js, reading);
    }
}

// Parse support function.
template <typename T, typename T1>
void serialize(T& val, json& js, bool reading, const std::map<T, T1>& table) {
    if (reading) {
        auto v = T1();
        serialize(v, js, reading);
        auto found = false;
        for (auto& kv : table) {
            if (kv.second == v) {
                val = kv.first;
                found = true;
                break;
            }
        }
        if (!found) throw std::runtime_error("bad enum value");
    } else {
        auto found = false;
        auto v = T1();
        for (auto& kv : table) {
            if (kv.first == val) {
                v = kv.second;
                found = true;
                break;
            }
        }
        if (!found) throw std::runtime_error("invalid value");
        serialize(v, js, reading);
    }
}

// Parse support function.
void serialize(vec2f& vals, json& js, bool reading) {
    serialize((std::array<float, 2>&)vals, js, reading);
}
void serialize(vec3f& vals, json& js, bool reading) {
    serialize((std::array<float, 3>&)vals, js, reading);
}
void serialize(vec4f& vals, json& js, bool reading) {
    serialize((std::array<float, 4>&)vals, js, reading);
}
void serialize(vec2i& vals, json& js, bool reading) {
    serialize((std::array<int, 2>&)vals, js, reading);
}
void serialize(vec3i& vals, json& js, bool reading) {
    serialize((std::array<int, 3>&)vals, js, reading);
}
void serialize(vec4i& vals, json& js, bool reading) {
    serialize((std::array<int, 4>&)vals, js, reading);
}
void serialize(mat3f& vals, json& js, bool reading) {
    serialize((std::array<float, 9>&)vals, js, reading);
}
void serialize(mat4f& vals, json& js, bool reading) {
    serialize((std::array<float, 16>&)vals, js, reading);
}
void serialize(frame3f& vals, json& js, bool reading) {
    serialize((std::array<float, 12>&)vals, js, reading);
}

// Parse support function.
void serialize_obj(json& js, bool reading) {
    if (reading) {
        if (!js.is_object()) throw std::runtime_error("object expected");
    } else {
        if (!js.is_object()) js = json::object();
    }
}

// Parse support function.
template <typename T>
void serialize_attr(T& val, json& js, const char* name, bool reading,
    bool required = true, const T& def = {}) {
    if (reading) {
        if (required) {
            if (!js.count(name)) throw std::runtime_error("missing value");
            serialize(val, js.at(name), reading);
        } else {
            if (!js.count(name))
                val = def;
            else
                serialize(val, js.at(name), reading);
        }
    } else {
        if (required || val != def) serialize(val, js[name], reading);
    }
}

// Dump support function.
template <typename T>
void serialize_attr(std::vector<T>& val, json& js, const char* name,
    bool reading, bool required = true, const std::vector<T>& def = {}) {
    if (reading) {
        if (required) {
            if (!js.count(name)) throw std::runtime_error("missing value");
            serialize(val, js.at(name), reading);
        } else {
            if (!js.count(name))
                val = def;
            else
                serialize(val, js.at(name), reading);
        }
    } else {
        if (required || !val.empty()) serialize(val, js[name], reading);
    }
}

// Dump support function.
template <typename T>
void serialize_attr(std::map<std::string, T>& val, json& js, const char* name,
    bool reading, bool required = true,
    const std::map<std::string, T>& def = {}) {
    if (reading) {
        if (required) {
            if (!js.count(name)) throw std::runtime_error("missing value");
            serialize(val, js.at(name), reading);
        } else {
            if (!js.count(name))
                val = def;
            else
                serialize(val, js.at(name), reading);
        }
    } else {
        if (required || !val.empty()) serialize(val, js[name], reading);
    }
}

// #codegen begin gltf-func

// Check for default value
template <typename T>
bool operator==(const glTFid<T>& a, const glTFid<T>& b) {
    return (int)a == (int)b;
}

// Check for default value
template <typename T>
bool operator!=(const glTFid<T>& a, const glTFid<T>& b) {
    return (int)a != (int)b;
}

// Parse id function.
template <typename T>
void serialize(glTFid<T>& val, json& js, bool reading) {
    if (reading) {
        if (!js.is_number_integer()) throw std::runtime_error("int expected");
        val = glTFid<T>((int)js);
    } else {
        js = (int)val;
    }
}

// Parses a glTFProperty object
void serialize(glTFProperty& val, json& js, bool reading) {
    if (reading) {
        if (!js.is_object()) throw std::runtime_error("object expected");
#if YGL_GLTFJSON
        if (js.count("extensions"))
            serialize(val.extensions, js.at("extensions"), reading);
        if (js.count("extras")) serialize(val.extras, js.at("extras"), reading);
#endif
    } else {
        if (!js.is_object()) js = json::object();
#if YGL_GLTFJSON
        if (!val.extensions.empty())
            serialize(val.extensions, js["extensions"], reading);
        if (!val.extras.is_null()) serialize(val.extras, js["extras"], reading);
#endif
    }
}

// Parses a glTFChildOfRootProperty object
void serialize(glTFChildOfRootProperty& val, json& js, bool reading) {
    static auto def = glTFChildOfRootProperty();
    serialize_obj(js, reading);
    serialize((glTFProperty&)val, js, reading);
    serialize_attr(val.name, js, "name", reading, false, def.name);
}
// Parse a glTFAccessorSparseIndicesComponentType enum
void serialize(
    glTFAccessorSparseIndicesComponentType& val, json& js, bool reading) {
    static std::vector<std::pair<int, glTFAccessorSparseIndicesComponentType>>
        table = {
            {5121, glTFAccessorSparseIndicesComponentType::UnsignedByte},
            {5123, glTFAccessorSparseIndicesComponentType::UnsignedShort},
            {5125, glTFAccessorSparseIndicesComponentType::UnsignedInt},
        };
    serialize(val, js, reading, table);
}

// Parses a glTFAccessorSparseIndices object
void serialize(glTFAccessorSparseIndices& val, json& js, bool reading) {
    static auto def = glTFAccessorSparseIndices();
    serialize_obj(js, reading);
    serialize((glTFProperty&)val, js, reading);
    serialize_attr(
        val.bufferView, js, "bufferView", reading, true, def.bufferView);
    serialize_attr(
        val.byteOffset, js, "byteOffset", reading, false, def.byteOffset);
    serialize_attr(val.componentType, js, "componentType", reading, true,
        def.componentType);
}

// Parses a glTFAccessorSparseValues object
void serialize(glTFAccessorSparseValues& val, json& js, bool reading) {
    static auto def = glTFAccessorSparseValues();
    serialize_obj(js, reading);
    serialize((glTFProperty&)val, js, reading);
    serialize_attr(
        val.bufferView, js, "bufferView", reading, true, def.bufferView);
    serialize_attr(
        val.byteOffset, js, "byteOffset", reading, false, def.byteOffset);
}

// Parses a glTFAccessorSparse object
void serialize(glTFAccessorSparse& val, json& js, bool reading) {
    static auto def = glTFAccessorSparse();
    serialize_obj(js, reading);
    serialize((glTFProperty&)val, js, reading);
    serialize_attr(val.count, js, "count", reading, true, def.count);
    serialize_attr(val.indices, js, "indices", reading, true, def.indices);
    serialize_attr(val.values, js, "values", reading, true, def.values);
}
// Parse a glTFAccessorComponentType enum
void serialize(glTFAccessorComponentType& val, json& js, bool reading) {
    static std::vector<std::pair<int, glTFAccessorComponentType>> table = {
        {5120, glTFAccessorComponentType::Byte},
        {5121, glTFAccessorComponentType::UnsignedByte},
        {5122, glTFAccessorComponentType::Short},
        {5123, glTFAccessorComponentType::UnsignedShort},
        {5125, glTFAccessorComponentType::UnsignedInt},
        {5126, glTFAccessorComponentType::Float},
    };
    serialize(val, js, reading, table);
}

// Parse a glTFAccessorType enum
void serialize(glTFAccessorType& val, json& js, bool reading) {
    static std::vector<std::pair<std::string, glTFAccessorType>> table = {
        {"SCALAR", glTFAccessorType::Scalar},
        {"VEC2", glTFAccessorType::Vec2},
        {"VEC3", glTFAccessorType::Vec3},
        {"VEC4", glTFAccessorType::Vec4},
        {"MAT2", glTFAccessorType::Mat2},
        {"MAT3", glTFAccessorType::Mat3},
        {"MAT4", glTFAccessorType::Mat4},
    };
    serialize(val, js, reading, table);
}

// Parses a glTFAccessor object
void serialize(glTFAccessor& val, json& js, bool reading) {
    static auto def = glTFAccessor();
    serialize_obj(js, reading);
    serialize((glTFChildOfRootProperty&)val, js, reading);
    serialize_attr(
        val.bufferView, js, "bufferView", reading, false, def.bufferView);
    serialize_attr(
        val.byteOffset, js, "byteOffset", reading, false, def.byteOffset);
    serialize_attr(val.componentType, js, "componentType", reading, true,
        def.componentType);
    serialize_attr(
        val.normalized, js, "normalized", reading, false, def.normalized);
    serialize_attr(val.count, js, "count", reading, true, def.count);
    serialize_attr(val.type, js, "type", reading, true, def.type);
    serialize_attr(val.max, js, "max", reading, false, def.max);
    serialize_attr(val.min, js, "min", reading, false, def.min);
    serialize_attr(val.sparse, js, "sparse", reading, false, def.sparse);
}
// Parse a glTFAnimationChannelTargetPath enum
void serialize(glTFAnimationChannelTargetPath& val, json& js, bool reading) {
    static std::vector<std::pair<std::string, glTFAnimationChannelTargetPath>>
        table = {
            {"translation", glTFAnimationChannelTargetPath::Translation},
            {"rotation", glTFAnimationChannelTargetPath::Rotation},
            {"scale", glTFAnimationChannelTargetPath::Scale},
            {"weights", glTFAnimationChannelTargetPath::Weights},
        };
    serialize(val, js, reading, table);
}

// Parses a glTFAnimationChannelTarget object
void serialize(glTFAnimationChannelTarget& val, json& js, bool reading) {
    static auto def = glTFAnimationChannelTarget();
    serialize_obj(js, reading);
    serialize((glTFProperty&)val, js, reading);
    serialize_attr(val.node, js, "node", reading, true, def.node);
    serialize_attr(val.path, js, "path", reading, true, def.path);
}

// Parses a glTFAnimationChannel object
void serialize(glTFAnimationChannel& val, json& js, bool reading) {
    static auto def = glTFAnimationChannel();
    serialize_obj(js, reading);
    serialize((glTFProperty&)val, js, reading);
    serialize_attr(val.sampler, js, "sampler", reading, true, def.sampler);
    serialize_attr(val.target, js, "target", reading, true, def.target);
}
// Parse a glTFAnimationSamplerInterpolation enum
void serialize(glTFAnimationSamplerInterpolation& val, json& js, bool reading) {
    static std::vector<
        std::pair<std::string, glTFAnimationSamplerInterpolation>>
        table = {
            {"LINEAR", glTFAnimationSamplerInterpolation::Linear},
            {"STEP", glTFAnimationSamplerInterpolation::Step},
            {"CUBICSPLINE", glTFAnimationSamplerInterpolation::CubicSpline},
        };
    serialize(val, js, reading, table);
}

// Parses a glTFAnimationSampler object
void serialize(glTFAnimationSampler& val, json& js, bool reading) {
    static auto def = glTFAnimationSampler();
    serialize_obj(js, reading);
    serialize((glTFProperty&)val, js, reading);
    serialize_attr(val.input, js, "input", reading, true, def.input);
    serialize_attr(val.interpolation, js, "interpolation", reading, false,
        def.interpolation);
    serialize_attr(val.output, js, "output", reading, true, def.output);
}

// Parses a glTFAnimation object
void serialize(glTFAnimation& val, json& js, bool reading) {
    static auto def = glTFAnimation();
    serialize_obj(js, reading);
    serialize((glTFChildOfRootProperty&)val, js, reading);
    serialize_attr(val.channels, js, "channels", reading, true, def.channels);
    serialize_attr(val.samplers, js, "samplers", reading, true, def.samplers);
}

// Parses a glTFAsset object
void serialize(glTFAsset& val, json& js, bool reading) {
    static auto def = glTFAsset();
    serialize_obj(js, reading);
    serialize((glTFProperty&)val, js, reading);
    serialize_attr(
        val.copyright, js, "copyright", reading, false, def.copyright);
    serialize_attr(
        val.generator, js, "generator", reading, false, def.generator);
    serialize_attr(val.version, js, "version", reading, true, def.version);
    serialize_attr(
        val.minVersion, js, "minVersion", reading, false, def.minVersion);
}

// Parses a glTFBuffer object
void serialize(glTFBuffer& val, json& js, bool reading) {
    static auto def = glTFBuffer();
    serialize_obj(js, reading);
    serialize((glTFChildOfRootProperty&)val, js, reading);
    serialize_attr(val.uri, js, "uri", reading, false, def.uri);
    serialize_attr(
        val.byteLength, js, "byteLength", reading, true, def.byteLength);
}
// Parse a glTFBufferViewTarget enum
void serialize(glTFBufferViewTarget& val, json& js, bool reading) {
    static std::vector<std::pair<int, glTFBufferViewTarget>> table = {
        {34962, glTFBufferViewTarget::ArrayBuffer},
        {34963, glTFBufferViewTarget::ElementArrayBuffer},
    };
    serialize(val, js, reading, table);
}

// Parses a glTFBufferView object
void serialize(glTFBufferView& val, json& js, bool reading) {
    static auto def = glTFBufferView();
    serialize_obj(js, reading);
    serialize((glTFChildOfRootProperty&)val, js, reading);
    serialize_attr(val.buffer, js, "buffer", reading, true, def.buffer);
    serialize_attr(
        val.byteOffset, js, "byteOffset", reading, false, def.byteOffset);
    serialize_attr(
        val.byteLength, js, "byteLength", reading, true, def.byteLength);
    serialize_attr(
        val.byteStride, js, "byteStride", reading, false, def.byteStride);
    serialize_attr(val.target, js, "target", reading, false, def.target);
}

// Parses a glTFCameraOrthographic object
void serialize(glTFCameraOrthographic& val, json& js, bool reading) {
    static auto def = glTFCameraOrthographic();
    serialize_obj(js, reading);
    serialize((glTFProperty&)val, js, reading);
    serialize_attr(val.xmag, js, "xmag", reading, true, def.xmag);
    serialize_attr(val.ymag, js, "ymag", reading, true, def.ymag);
    serialize_attr(val.zfar, js, "zfar", reading, true, def.zfar);
    serialize_attr(val.znear, js, "znear", reading, true, def.znear);
}

// Parses a glTFCameraPerspective object
void serialize(glTFCameraPerspective& val, json& js, bool reading) {
    static auto def = glTFCameraPerspective();
    serialize_obj(js, reading);
    serialize((glTFProperty&)val, js, reading);
    serialize_attr(
        val.aspectRatio, js, "aspectRatio", reading, false, def.aspectRatio);
    serialize_attr(val.yfov, js, "yfov", reading, true, def.yfov);
    serialize_attr(val.zfar, js, "zfar", reading, false, def.zfar);
    serialize_attr(val.znear, js, "znear", reading, true, def.znear);
}
// Parse a glTFCameraType enum
void serialize(glTFCameraType& val, json& js, bool reading) {
    static std::vector<std::pair<std::string, glTFCameraType>> table = {
        {"perspective", glTFCameraType::Perspective},
        {"orthographic", glTFCameraType::Orthographic},
    };
    serialize(val, js, reading, table);
}

// Parses a glTFCamera object
void serialize(glTFCamera& val, json& js, bool reading) {
    static auto def = glTFCamera();
    serialize_obj(js, reading);
    serialize((glTFChildOfRootProperty&)val, js, reading);
    serialize_attr(
        val.orthographic, js, "orthographic", reading, false, def.orthographic);
    serialize_attr(
        val.perspective, js, "perspective", reading, false, def.perspective);
    serialize_attr(val.type, js, "type", reading, true, def.type);
}
// Parse a glTFImageMimeType enum
void serialize(glTFImageMimeType& val, json& js, bool reading) {
    static std::vector<std::pair<std::string, glTFImageMimeType>> table = {
        {"image/jpeg", glTFImageMimeType::ImageJpeg},
        {"image/png", glTFImageMimeType::ImagePng},
    };
    serialize(val, js, reading, table);
}

// Parses a glTFImage object
void serialize(glTFImage& val, json& js, bool reading) {
    static auto def = glTFImage();
    serialize_obj(js, reading);
    serialize((glTFChildOfRootProperty&)val, js, reading);
    serialize_attr(val.uri, js, "uri", reading, false, def.uri);
    serialize_attr(val.mimeType, js, "mimeType", reading, false, def.mimeType);
    serialize_attr(
        val.bufferView, js, "bufferView", reading, false, def.bufferView);
}

// Parses a glTFTextureInfo object
void serialize(glTFTextureInfo& val, json& js, bool reading) {
    static auto def = glTFTextureInfo();
    serialize_obj(js, reading);
    serialize((glTFProperty&)val, js, reading);
    serialize_attr(val.index, js, "index", reading, true, def.index);
    serialize_attr(val.texCoord, js, "texCoord", reading, false, def.texCoord);
}

// Parses a glTFTexture object
void serialize(glTFTexture& val, json& js, bool reading) {
    static auto def = glTFTexture();
    serialize_obj(js, reading);
    serialize((glTFChildOfRootProperty&)val, js, reading);
    serialize_attr(val.sampler, js, "sampler", reading, false, def.sampler);
    serialize_attr(val.source, js, "source", reading, false, def.source);
}

// Parses a glTFMaterialNormalTextureInfo object
void serialize(glTFMaterialNormalTextureInfo& val, json& js, bool reading) {
    static auto def = glTFMaterialNormalTextureInfo();
    serialize_obj(js, reading);
    serialize((glTFTextureInfo&)val, js, reading);
    serialize_attr(val.scale, js, "scale", reading, false, def.scale);
}

// Parses a glTFMaterialOcclusionTextureInfo object
void serialize(glTFMaterialOcclusionTextureInfo& val, json& js, bool reading) {
    static auto def = glTFMaterialOcclusionTextureInfo();
    serialize_obj(js, reading);
    serialize((glTFTextureInfo&)val, js, reading);
    serialize_attr(val.strength, js, "strength", reading, false, def.strength);
}

// Parses a glTFMaterialPbrMetallicRoughness object
void serialize(glTFMaterialPbrMetallicRoughness& val, json& js, bool reading) {
    static auto def = glTFMaterialPbrMetallicRoughness();
    serialize_obj(js, reading);
    serialize((glTFProperty&)val, js, reading);
    serialize_attr(val.baseColorFactor, js, "baseColorFactor", reading, false,
        def.baseColorFactor);
    serialize_attr(val.baseColorTexture, js, "baseColorTexture", reading, false,
        def.baseColorTexture);
    serialize_attr(val.metallicFactor, js, "metallicFactor", reading, false,
        def.metallicFactor);
    serialize_attr(val.roughnessFactor, js, "roughnessFactor", reading, false,
        def.roughnessFactor);
    serialize_attr(val.metallicRoughnessTexture, js, "metallicRoughnessTexture",
        reading, false, def.metallicRoughnessTexture);
}

// Parses a glTFMaterialPbrSpecularGlossiness object
void serialize(glTFMaterialPbrSpecularGlossiness& val, json& js, bool reading) {
    static auto def = glTFMaterialPbrSpecularGlossiness();
    serialize_obj(js, reading);
    serialize((glTFProperty&)val, js, reading);
    serialize_attr(val.diffuseFactor, js, "diffuseFactor", reading, false,
        def.diffuseFactor);
    serialize_attr(val.diffuseTexture, js, "diffuseTexture", reading, false,
        def.diffuseTexture);
    serialize_attr(val.specularFactor, js, "specularFactor", reading, false,
        def.specularFactor);
    serialize_attr(val.glossinessFactor, js, "glossinessFactor", reading, false,
        def.glossinessFactor);
    serialize_attr(val.specularGlossinessTexture, js,
        "specularGlossinessTexture", reading, false,
        def.specularGlossinessTexture);
}
// Parse a glTFMaterialAlphaMode enum
void serialize(glTFMaterialAlphaMode& val, json& js, bool reading) {
    static std::vector<std::pair<std::string, glTFMaterialAlphaMode>> table = {
        {"OPAQUE", glTFMaterialAlphaMode::Opaque},
        {"MASK", glTFMaterialAlphaMode::Mask},
        {"BLEND", glTFMaterialAlphaMode::Blend},
    };
    serialize(val, js, reading, table);
}

// Parses a glTFMaterial object
void serialize(glTFMaterial& val, json& js, bool reading) {
    static auto def = glTFMaterial();
    serialize_obj(js, reading);
    serialize((glTFChildOfRootProperty&)val, js, reading);
    serialize_attr(val.pbrMetallicRoughness, js, "pbrMetallicRoughness",
        reading, false, def.pbrMetallicRoughness);
    serialize_attr(val.normalTexture, js, "normalTexture", reading, false,
        def.normalTexture);
    serialize_attr(val.occlusionTexture, js, "occlusionTexture", reading, false,
        def.occlusionTexture);
    serialize_attr(val.emissiveTexture, js, "emissiveTexture", reading, false,
        def.emissiveTexture);
    serialize_attr(val.emissiveFactor, js, "emissiveFactor", reading, false,
        def.emissiveFactor);
    serialize_attr(
        val.alphaMode, js, "alphaMode", reading, false, def.alphaMode);
    serialize_attr(
        val.alphaCutoff, js, "alphaCutoff", reading, false, def.alphaCutoff);
    serialize_attr(
        val.doubleSided, js, "doubleSided", reading, false, def.doubleSided);
    if (reading) {
        if (js.count("extensions")) {
            auto& js_ext = js["extensions"];
            serialize_attr(val.pbrSpecularGlossiness, js_ext,
                "KHR_materials_pbrSpecularGlossiness", reading, false,
                def.pbrSpecularGlossiness);
        }
    } else {
        if (val.pbrSpecularGlossiness != nullptr) {
            auto& js_ext = js["extensions"];
            serialize_attr(val.pbrSpecularGlossiness, js_ext,
                "KHR_materials_pbrSpecularGlossiness", reading, false,
                def.pbrSpecularGlossiness);
        }
    }
}
// Parse a glTFMeshPrimitiveMode enum
void serialize(glTFMeshPrimitiveMode& val, json& js, bool reading) {
    static std::vector<std::pair<int, glTFMeshPrimitiveMode>> table = {
        {0, glTFMeshPrimitiveMode::Points},
        {1, glTFMeshPrimitiveMode::Lines},
        {2, glTFMeshPrimitiveMode::LineLoop},
        {3, glTFMeshPrimitiveMode::LineStrip},
        {4, glTFMeshPrimitiveMode::Triangles},
        {5, glTFMeshPrimitiveMode::TriangleStrip},
        {6, glTFMeshPrimitiveMode::TriangleFan},
    };
    serialize(val, js, reading, table);
}

// Parses a glTFMeshPrimitive object
void serialize(glTFMeshPrimitive& val, json& js, bool reading) {
    static auto def = glTFMeshPrimitive();
    serialize_obj(js, reading);
    serialize((glTFProperty&)val, js, reading);
    serialize_attr(
        val.attributes, js, "attributes", reading, true, def.attributes);
    serialize_attr(val.indices, js, "indices", reading, false, def.indices);
    serialize_attr(val.material, js, "material", reading, false, def.material);
    serialize_attr(val.mode, js, "mode", reading, false, def.mode);
    serialize_attr(val.targets, js, "targets", reading, false, def.targets);
}

// Parses a glTFMesh object
void serialize(glTFMesh& val, json& js, bool reading) {
    static auto def = glTFMesh();
    serialize_obj(js, reading);
    serialize((glTFChildOfRootProperty&)val, js, reading);
    serialize_attr(
        val.primitives, js, "primitives", reading, true, def.primitives);
    serialize_attr(val.weights, js, "weights", reading, false, def.weights);
}

// Parses a glTFNode object
void serialize(glTFNode& val, json& js, bool reading) {
    static auto def = glTFNode();
    serialize_obj(js, reading);
    serialize((glTFChildOfRootProperty&)val, js, reading);
    serialize_attr(val.camera, js, "camera", reading, false, def.camera);
    serialize_attr(val.children, js, "children", reading, false, def.children);
    serialize_attr(val.skin, js, "skin", reading, false, def.skin);
    serialize_attr(val.matrix, js, "matrix", reading, false, def.matrix);
    serialize_attr(val.mesh, js, "mesh", reading, false, def.mesh);
    serialize_attr(val.rotation, js, "rotation", reading, false, def.rotation);
    serialize_attr(val.scale, js, "scale", reading, false, def.scale);
    serialize_attr(
        val.translation, js, "translation", reading, false, def.translation);
    serialize_attr(val.weights, js, "weights", reading, false, def.weights);
}
// Parse a glTFSamplerMagFilter enum
void serialize(glTFSamplerMagFilter& val, json& js, bool reading) {
    static std::vector<std::pair<int, glTFSamplerMagFilter>> table = {
        {9728, glTFSamplerMagFilter::Nearest},
        {9729, glTFSamplerMagFilter::Linear},
    };
    serialize(val, js, reading, table);
}

// Parse a glTFSamplerMinFilter enum
void serialize(glTFSamplerMinFilter& val, json& js, bool reading) {
    static std::vector<std::pair<int, glTFSamplerMinFilter>> table = {
        {9728, glTFSamplerMinFilter::Nearest},
        {9729, glTFSamplerMinFilter::Linear},
        {9984, glTFSamplerMinFilter::NearestMipmapNearest},
        {9985, glTFSamplerMinFilter::LinearMipmapNearest},
        {9986, glTFSamplerMinFilter::NearestMipmapLinear},
        {9987, glTFSamplerMinFilter::LinearMipmapLinear},
    };
    serialize(val, js, reading, table);
}

// Parse a glTFSamplerWrapS enum
void serialize(glTFSamplerWrapS& val, json& js, bool reading) {
    static std::vector<std::pair<int, glTFSamplerWrapS>> table = {
        {33071, glTFSamplerWrapS::ClampToEdge},
        {33648, glTFSamplerWrapS::MirroredRepeat},
        {10497, glTFSamplerWrapS::Repeat},
    };
    serialize(val, js, reading, table);
}

// Parse a glTFSamplerWrapT enum
void serialize(glTFSamplerWrapT& val, json& js, bool reading) {
    static std::vector<std::pair<int, glTFSamplerWrapT>> table = {
        {33071, glTFSamplerWrapT::ClampToEdge},
        {33648, glTFSamplerWrapT::MirroredRepeat},
        {10497, glTFSamplerWrapT::Repeat},
    };
    serialize(val, js, reading, table);
}

// Parses a glTFSampler object
void serialize(glTFSampler& val, json& js, bool reading) {
    static auto def = glTFSampler();
    serialize_obj(js, reading);
    serialize((glTFChildOfRootProperty&)val, js, reading);
    serialize_attr(
        val.magFilter, js, "magFilter", reading, false, def.magFilter);
    serialize_attr(
        val.minFilter, js, "minFilter", reading, false, def.minFilter);
    serialize_attr(val.wrapS, js, "wrapS", reading, false, def.wrapS);
    serialize_attr(val.wrapT, js, "wrapT", reading, false, def.wrapT);
}

// Parses a glTFScene object
void serialize(glTFScene& val, json& js, bool reading) {
    static auto def = glTFScene();
    serialize_obj(js, reading);
    serialize((glTFChildOfRootProperty&)val, js, reading);
    serialize_attr(val.nodes, js, "nodes", reading, false, def.nodes);
}

// Parses a glTFSkin object
void serialize(glTFSkin& val, json& js, bool reading) {
    static auto def = glTFSkin();
    serialize_obj(js, reading);
    serialize((glTFChildOfRootProperty&)val, js, reading);
    serialize_attr(val.inverseBindMatrices, js, "inverseBindMatrices", reading,
        false, def.inverseBindMatrices);
    serialize_attr(val.skeleton, js, "skeleton", reading, false, def.skeleton);
    serialize_attr(val.joints, js, "joints", reading, true, def.joints);
}

// Parses a glTF object
void serialize(glTF& val, json& js, bool reading) {
    static auto def = glTF();
    serialize_obj(js, reading);
    serialize((glTFProperty&)val, js, reading);
    serialize_attr(val.extensionsUsed, js, "extensionsUsed", reading, false,
        def.extensionsUsed);
    serialize_attr(val.extensionsRequired, js, "extensionsRequired", reading,
        false, def.extensionsRequired);
    serialize_attr(
        val.accessors, js, "accessors", reading, false, def.accessors);
    serialize_attr(
        val.animations, js, "animations", reading, false, def.animations);
    serialize_attr(val.asset, js, "asset", reading, true, def.asset);
    serialize_attr(val.buffers, js, "buffers", reading, false, def.buffers);
    serialize_attr(
        val.bufferViews, js, "bufferViews", reading, false, def.bufferViews);
    serialize_attr(val.cameras, js, "cameras", reading, false, def.cameras);
    serialize_attr(val.images, js, "images", reading, false, def.images);
    serialize_attr(
        val.materials, js, "materials", reading, false, def.materials);
    serialize_attr(val.meshes, js, "meshes", reading, false, def.meshes);
    serialize_attr(val.nodes, js, "nodes", reading, false, def.nodes);
    serialize_attr(val.samplers, js, "samplers", reading, false, def.samplers);
    serialize_attr(val.scene, js, "scene", reading, false, def.scene);
    serialize_attr(val.scenes, js, "scenes", reading, false, def.scenes);
    serialize_attr(val.skins, js, "skins", reading, false, def.skins);
    serialize_attr(val.textures, js, "textures", reading, false, def.textures);
}

// #codegen end gltf-func

// Encode in base64
std::string base64_encode(
    unsigned char const* bytes_to_encode, unsigned int in_len) {
    static const std::string base64_chars =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz"
        "0123456789+/";

    std::string ret;
    int i = 0;
    int j = 0;
    unsigned char char_array_3[3];
    unsigned char char_array_4[4];

    while (in_len--) {
        char_array_3[i++] = *(bytes_to_encode++);
        if (i == 3) {
            char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
            char_array_4[1] = ((char_array_3[0] & 0x03) << 4) +
                              ((char_array_3[1] & 0xf0) >> 4);
            char_array_4[2] = ((char_array_3[1] & 0x0f) << 2) +
                              ((char_array_3[2] & 0xc0) >> 6);
            char_array_4[3] = char_array_3[2] & 0x3f;

            for (i = 0; (i < 4); i++) ret += base64_chars[char_array_4[i]];
            i = 0;
        }
    }

    if (i) {
        for (j = i; j < 3; j++) char_array_3[j] = '\0';

        char_array_4[0] = (char_array_3[0] & 0xfc) >> 2;
        char_array_4[1] =
            ((char_array_3[0] & 0x03) << 4) + ((char_array_3[1] & 0xf0) >> 4);
        char_array_4[2] =
            ((char_array_3[1] & 0x0f) << 2) + ((char_array_3[2] & 0xc0) >> 6);
        char_array_4[3] = char_array_3[2] & 0x3f;

        for (j = 0; (j < i + 1); j++) ret += base64_chars[char_array_4[j]];

        while ((i++ < 3)) ret += '=';
    }

    return ret;
}

// Decode from base64
std::string base64_decode(std::string const& encoded_string) {
    static const std::string base64_chars =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz"
        "0123456789+/";

    auto is_base64 = [](unsigned char c) -> bool {
        return (isalnum(c) || (c == '+') || (c == '/'));
    };

    int in_len = (int)encoded_string.size();
    int i = 0;
    int j = 0;
    int in_ = 0;
    unsigned char char_array_4[4], char_array_3[3];
    std::string ret;

    while (in_len-- && (encoded_string[in_] != '=') &&
           is_base64(encoded_string[in_])) {
        char_array_4[i++] = encoded_string[in_];
        in_++;
        if (i == 4) {
            for (i = 0; i < 4; i++)
                char_array_4[i] = base64_chars.find(char_array_4[i]);

            char_array_3[0] =
                (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
            char_array_3[1] = ((char_array_4[1] & 0xf) << 4) +
                              ((char_array_4[2] & 0x3c) >> 2);
            char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];

            for (i = 0; (i < 3); i++) ret += char_array_3[i];
            i = 0;
        }
    }

    if (i) {
        for (j = i; j < 4; j++) char_array_4[j] = 0;

        for (j = 0; j < 4; j++)
            char_array_4[j] = base64_chars.find(char_array_4[j]);

        char_array_3[0] =
            (char_array_4[0] << 2) + ((char_array_4[1] & 0x30) >> 4);
        char_array_3[1] =
            ((char_array_4[1] & 0xf) << 4) + ((char_array_4[2] & 0x3c) >> 2);
        char_array_3[2] = ((char_array_4[2] & 0x3) << 6) + char_array_4[3];

        for (j = 0; (j < i - 1); j++) ret += char_array_3[j];
    }

    return ret;
}

// Load buffer data.
void load_buffers(
    const glTF* gltf, const std::string& dirname, bool skip_missing) {
    auto fix_path = [](const std::string& path_) {
        auto path = path_;
        for (auto& c : path)
            if (c == '\\') c = '/';
        return path;
    };
    auto startswith = [](const std::string& str, const std::string& substr) {
        if (str.length() < substr.length()) return false;
        for (auto i = 0; i < substr.length(); i++)
            if (str[i] != substr[i]) return false;
        return true;
    };
    auto load_binary = [](const std::string& filename) {
        // https://stackoverflow.com/questions/174531/easiest-way-to-get-files-contents-in-c
        auto f = fopen(filename.c_str(), "rb");
        if (!f) throw std::runtime_error("cannot read file " + filename);
        fseek(f, 0, SEEK_END);
        auto len = ftell(f);
        fseek(f, 0, SEEK_SET);
        auto buf = std::vector<unsigned char>(len);
        if (fread(buf.data(), 1, len, f) != len)
            throw std::runtime_error("cannot read file " + filename);
        fclose(f);
        return buf;
    };

    for (auto buffer : gltf->buffers) {
        if (buffer->uri == "") continue;
        try {
            if (startswith(buffer->uri, "data:")) {
                // assume it is base64 and find ','
                auto pos = buffer->uri.find(',');
                if (pos == buffer->uri.npos) {
                    if (skip_missing) continue;
                    throw std::runtime_error("could not decode base64 data");
                }
                // decode
                auto data = base64_decode(buffer->uri.substr(pos + 1));
                buffer->data =
                    std::vector<unsigned char>((unsigned char*)data.c_str(),
                        (unsigned char*)data.c_str() + data.length());
            } else {
                buffer->data = load_binary(fix_path(dirname + buffer->uri));
                if (buffer->data.empty()) {
                    if (skip_missing) continue;
                    throw std::runtime_error("could not load binary file " +
                                             fix_path(dirname + buffer->uri));
                }
            }
            if (buffer->byteLength != buffer->data.size()) {
                throw std::runtime_error("mismatched buffer size");
            }
        } catch (const std::exception&) {
            if (skip_missing) continue;
            throw;
        }
    }
}

// Loads images.
void load_images(
    const glTF* gltf, const std::string& dirname, bool skip_missing) {
    auto fix_path = [](const std::string& path_) {
        auto path = path_;
        for (auto& c : path)
            if (c == '\\') c = '/';
        return path;
    };
    auto startswith = [](const std::string& str, const std::string& substr) {
        if (str.length() < substr.length()) return false;
        for (auto i = 0; i < substr.length(); i++)
            if (str[i] != substr[i]) return false;
        return true;
    };

    for (auto image : gltf->images) {
        image->data = image_data();
        auto filename = std::string();
#if YGL_IMAGEIO
        if (image->bufferView || startswith(image->uri, "data:")) {
            auto buffer = std::string();
            auto data = (unsigned char*)nullptr;
            auto data_size = 0;
            if (image->bufferView) {
                auto view = gltf->get(image->bufferView);
                auto buffer = gltf->get(view->buffer);
                if (!view || !buffer || view->byteStride) {
                    if (skip_missing) continue;
                    throw std::runtime_error("invalid image buffer view");
                }
                if (image->mimeType == glTFImageMimeType::ImagePng)
                    filename = "internal_data.png";
                else if (image->mimeType == glTFImageMimeType::ImageJpeg)
                    filename = "internal_data.jpg";
                else {
                    if (skip_missing) continue;
                    throw std::runtime_error("unsupported image format");
                }
                data = buffer->data.data() + view->byteOffset;
                data_size = view->byteLength;
            } else {
                // assume it is base64 and find ','
                auto pos = image->uri.find(',');
                if (pos == image->uri.npos) {
                    if (skip_missing) continue;
                    throw std::runtime_error("could not decode base64 data");
                }
                auto header = image->uri.substr(0, pos);
                for (auto format : {"png", "jpg", "jpeg", "tga", "ppm", "hdr"})
                    if (header.find(format) != header.npos)
                        filename = std::string("fake.") + format;
                if (is_hdr_filename(filename)) {
                    if (skip_missing) continue;
                    throw std::runtime_error(
                        "unsupported embedded image format " +
                        header.substr(0, pos));
                }
                // decode
                buffer = base64_decode(image->uri.substr(pos + 1));
                data_size = (int)buffer.size();
                data = (unsigned char*)buffer.data();
            }
            if (is_hdr_filename(filename)) {
                image->data.dataf = load_imagef_from_memory(filename, data,
                    data_size, image->data.width, image->data.height,
                    image->data.ncomp);
            } else {
                image->data.datab = load_imageb_from_memory(filename, data,
                    data_size, image->data.width, image->data.height,
                    image->data.ncomp);
            }
        } else {
            filename = fix_path(dirname + image->uri);
            if (is_hdr_filename(filename)) {
                image->data.dataf = load_imagef(filename, image->data.width,
                    image->data.height, image->data.ncomp);
            } else {
                image->data.datab = load_imageb(filename, image->data.width,
                    image->data.height, image->data.ncomp);
            }
        }
#endif
        if (image->data.dataf.empty() && image->data.datab.empty()) {
            if (skip_missing) continue;
            throw std::runtime_error("cannot load image " + filename);
        }
    }
}

// Loads a gltf.
glTF* load_gltf(const std::string& filename, bool load_bin, bool load_image,
    bool skip_missing) {
    // clear data
    auto gltf = new glTF();

    // load json
    std::ifstream stream(filename.c_str());
    if (!stream) throw std::runtime_error("could not load json " + filename);
    auto js = json();
    try {
        stream >> js;
    } catch (const std::exception& e) {
        throw std::runtime_error(
            std::string("could not load json with error ") + e.what());
    }

    // parse json
    try {
        serialize(gltf, js, true);
    } catch (const std::exception& e) {
        throw std::runtime_error("error parsing gltf " + std::string(e.what()));
    }

    // load external resources
    auto dirname = path_dirname(filename);
    if (load_bin) load_buffers(gltf, dirname, skip_missing);
    if (load_image) load_images(gltf, dirname, skip_missing);

    // done
    return gltf;
}

// Save buffer data.
void save_buffers(
    const glTF* gltf, const std::string& dirname, bool skip_missing) {
    auto save_binary = [](const std::string& filename,
                           const std::vector<unsigned char>& data) {
        auto f = fopen(filename.c_str(), "wb");
        if (!f) throw std::runtime_error("cannot write file " + filename);
        auto num = fwrite(data.data(), 1, data.size(), f);
        if (num != data.size())
            throw std::runtime_error("cannot write file " + filename);
        fclose(f);
    };

    for (auto buffer : gltf->buffers) {
        auto startswith = [](const std::string& str,
                              const std::string& substr) {
            if (str.length() < substr.length()) return false;
            for (auto i = 0; i < substr.length(); i++)
                if (str[i] != substr[i]) return false;
            return true;
        };

        try {
            if (startswith(buffer->uri, "data:")) {
                throw std::runtime_error(
                    "saving of embedded data not supported");
            }
            save_binary(dirname + buffer->uri, buffer->data);
        } catch (const std::exception&) {
            if (skip_missing) continue;
            throw;
        }
    }
}

// Save images.
void save_images(
    const glTF* gltf, const std::string& dirname, bool skip_missing) {
    auto startswith = [](const std::string& str, const std::string& substr) {
        if (str.length() < substr.length()) return false;
        for (auto i = 0; i < substr.length(); i++)
            if (str[i] != substr[i]) return false;
        return true;
    };

    for (auto image : gltf->images) {
        try {
            if (startswith(image->uri, "data:")) {
                throw std::runtime_error(
                    "saving of embedded data not supported");
            }
            auto filename = dirname + image->uri;
            auto ok = false;
#if YGL_IMAGEIO
            if (!image->data.datab.empty()) {
                ok =
                    save_imageb(filename, image->data.width, image->data.height,
                        image->data.ncomp, image->data.datab.data());
            }
            if (!image->data.dataf.empty()) {
                ok =
                    save_imagef(filename, image->data.width, image->data.height,
                        image->data.ncomp, image->data.dataf.data());
            }
#endif
            if (!ok) {
                throw std::runtime_error("cannot save image " + filename);
            }
        } catch (const std::exception&) {
            if (skip_missing) continue;
            throw;
        }
    }
}

// Saves a gltf.
void save_gltf(const std::string& filename, const glTF* gltf, bool save_bin,
    bool save_image) {
    auto save_text = [](const std::string& filename, const std::string& str) {
        auto f = fopen(filename.c_str(), "wb");
        if (!f) throw std::runtime_error("cannot write file " + filename);
        auto num = fwrite(str.c_str(), 1, str.size(), f);
        if (num != str.size())
            throw std::runtime_error("cannot write file " + filename);
        fclose(f);
    };

    // dumps json
    auto js = json();
    serialize((glTF*&)gltf, js, false);

    // save json
    save_text(filename, js.dump(2));

    // save external resources
    auto dirname = path_dirname(filename);
    if (save_bin) save_buffers(gltf, dirname, false);
    if (save_image) save_images(gltf, dirname, false);
}

// reading shortcut
template <typename T>
void gltf_fread(FILE* f, T* v, int count) {
    if (fread(v, sizeof(T), count, f) != count)
        throw std::runtime_error("could not read binary file");
}

// writing shortcut
template <typename T>
void gltf_fwrite(FILE* f, const T* v, int count) {
    if (fwrite(v, sizeof(T), count, f) != count)
        std::runtime_error("could not write binary file");
}

// Loads a binary gltf.
glTF* load_binary_gltf(const std::string& filename, bool load_bin,
    bool load_image, bool skip_missing) {
    // clear data
    auto gltf = new glTF();

    // opens binary file
    auto f = fopen(filename.c_str(), "rb");
    if (!f) throw std::runtime_error("could not load binary file " + filename);

    // read magic
    uint32_t magic;
    gltf_fread(f, &magic, 1);
    if (magic != 0x46546C67) throw std::runtime_error("corrupted glb format");

    // read version
    uint32_t version;
    gltf_fread(f, &version, 1);
    if (version != 1 && version != 2)
        throw std::runtime_error("unsupported glb version");

    // read length
    uint32_t length;
    gltf_fread(f, &length, 1);

    // data
    auto json_bytes = std::vector<char>();
    auto buffer_bytes = std::vector<unsigned char>();
    uint32_t buffer_length = 0;

    if (version == 1) {
        // read content length and format
        uint32_t json_length, json_format;
        gltf_fread(f, &json_length, 1);
        gltf_fread(f, &json_format, 1);

        // read json bytes
        json_bytes.resize(json_length);
        gltf_fread(f, json_bytes.data(), json_length);

        // read buffer bytes
        if (load_bin) {
            buffer_bytes.resize(length - json_length - 20);
            gltf_fread(f, buffer_bytes.data(), (int)buffer_bytes.size());
            buffer_length = (int)buffer_bytes.size();
        }
    }

    if (version == 2) {
        // read content length and format
        uint32_t json_length, json_format;
        gltf_fread(f, &json_length, 1);
        gltf_fread(f, &json_format, 1);
        if (json_format != 0x4E4F534A) {
            throw std::runtime_error("corrupt binary format");
            return nullptr;
        }

        // read json bytes
        json_bytes.resize(json_length);
        gltf_fread(f, json_bytes.data(), (int)json_bytes.size());

        // read content length and format
        uint32_t buffer_format;
        gltf_fread(f, &buffer_length, 1);
        gltf_fread(f, &buffer_format, 1);
        if (buffer_format != 0x004E4942)
            throw std::runtime_error("corrupt binary format");

        // read buffer bytes
        if (load_bin) {
            buffer_bytes.resize(buffer_length);
            gltf_fread(f, buffer_bytes.data(), (int)buffer_bytes.size());
        }
    }

    // load json
    auto js = json();
    try {
        json_bytes.push_back(0);
        js = json::parse(json_bytes.data());
    } catch (const std::exception& e) {
        throw std::runtime_error(
            std::string("could not load json with error ") + e.what());
    }

    // parse json
    try {
        serialize(gltf, js, true);
    } catch (const std::exception& e) {
        throw std::runtime_error(
            "cannot parse gltf json " + std::string(e.what()));
        return nullptr;
    }

    // fix internal buffer
    auto buffer = gltf->buffers.at(0);
    buffer->byteLength = buffer_length;
    if (version == 2) buffer->uri = "";
    if (load_bin) { buffer->data = buffer_bytes; }

    // load external resources
    auto dirname = path_dirname(filename);
    if (load_bin) load_buffers(gltf, dirname, skip_missing);
    if (load_image) load_images(gltf, dirname, skip_missing);

    // close
    fclose(f);

    // done
    return gltf;
}

// Saves a binary gltf.
void save_binary_gltf(const std::string& filename, const glTF* gltf,
    bool save_bin, bool save_image) {
    // opens binary file
    auto f = fopen(filename.c_str(), "wb");
    if (!f) throw std::runtime_error("could not write binary file");

    // dumps json
    auto js = json();
    serialize((glTF*&)gltf, js, false);

    // fix string
    auto js_str = js.dump(2);
    while (js_str.length() % 4) js_str += " ";
    uint32_t json_length = (uint32_t)js_str.size();

    // internal buffer
    auto buffer = gltf->buffers.at(0);
    uint32_t buffer_length = buffer->byteLength;
    if (buffer_length % 4) buffer_length += 4 - buffer_length % 4;

    // write header
    uint32_t magic = 0x46546C67;
    gltf_fwrite(f, &magic, 1);
    uint32_t version = 2;
    gltf_fwrite(f, &version, 1);
    uint32_t length = 12 + 8 + json_length + 8 + buffer_length;
    gltf_fwrite(f, &length, 1);

    // write json
    uint32_t json_type = 0x4E4F534A;
    gltf_fwrite(f, &json_length, 1);
    gltf_fwrite(f, &json_type, 1);
    gltf_fwrite(f, js_str.data(), (int)json_length);

    if (save_bin) {
        uint32_t buffer_type = 0x004E4942;
        gltf_fwrite(f, &buffer_length, 1);
        gltf_fwrite(f, &buffer_type, 1);
        gltf_fwrite(f, buffer->data.data(), (int)buffer->data.size());
        char pad = 0;
        for (auto i = 0; i < buffer_length - buffer->data.size(); i++)
            gltf_fwrite(f, &pad, 1);
    }

    // close
    fclose(f);

    // save external resources
    auto dirname = path_dirname(filename);
    if (save_bin) save_buffers(gltf, dirname, false);
    if (save_image) save_images(gltf, dirname, false);
}

accessor_view::accessor_view(const glTF* gltf, const glTFAccessor* accessor) {
    _size = accessor->count;
    _ncomp = _num_components(accessor->type);
    _ctype = accessor->componentType;
    _normalize = accessor->normalized;
    auto buffer_view = gltf->get(accessor->bufferView);
    _stride = (buffer_view->byteStride) ? buffer_view->byteStride :
                                          (_ctype_size(_ctype) * _ncomp);
    auto buffer = gltf->get(buffer_view->buffer);
    _data =
        buffer->data.data() + accessor->byteOffset + buffer_view->byteOffset;
    auto remaining_buffer_bytes =
        buffer->data.size() - (_data - buffer->data.data());
    auto view_bytes = _size * _stride;
    _valid = remaining_buffer_bytes >= view_bytes;
    if (!_valid) throw std::runtime_error("corrupted glTF accessor view");
}

float accessor_view::get(int idx, int c) const {
    auto i = min(max(c, 0), ncomp() - 1);
    auto valb = _data + _stride * idx + i * _ctype_size(_ctype);
    // use double for integer conversion to attempt to maintain
    // precision
    if (!_normalize) {
        switch (_ctype) {
            case glTFAccessorComponentType::Float:
                return (float)(*(float*)valb);
            case glTFAccessorComponentType::Byte: return (float)(*(char*)valb);
            case glTFAccessorComponentType::UnsignedByte:
                return (float)(*(unsigned char*)valb);
            case glTFAccessorComponentType::Short:
                return (float)(*(short*)valb);
            case glTFAccessorComponentType::UnsignedShort:
                return (float)(*(unsigned short*)valb);
            case glTFAccessorComponentType::UnsignedInt:
                return (float)(*(unsigned int*)valb);
            case glTFAccessorComponentType::NotSet:
                throw std::runtime_error("bad enum value");
                break;
        }

    } else {
        switch (_ctype) {
            case glTFAccessorComponentType::Float:
                return (float)(*(float*)valb);
            case glTFAccessorComponentType::Byte:
                return (float)max((float)(c / 127.0), -1.0f);
            case glTFAccessorComponentType::UnsignedByte:
                return (float)(c / 255.0);
            case glTFAccessorComponentType::Short:
                return (float)(max((float)(c / 32767.0), -1.0f));
            case glTFAccessorComponentType::UnsignedShort:
                return (float)(c / 65535.0);
            case glTFAccessorComponentType::UnsignedInt:
                return (float)(max((float)(c / 2147483647.0), -1.0f));
            case glTFAccessorComponentType::NotSet:
                throw std::runtime_error("bad enum value");
                break;
        }
    }
    return 0;
}

int accessor_view::geti(int idx, int c) const {
    auto i = min(max(c, 0), ncomp() - 1);
    auto valb = _data + _stride * idx + i * _ctype_size(_ctype);
    // use double for integer conversion to attempt to maintain
    // precision
    switch (_ctype) {
        case glTFAccessorComponentType::Float: return (int)(*(float*)valb);
        case glTFAccessorComponentType::Byte: return (int)(*(char*)valb);
        case glTFAccessorComponentType::UnsignedByte:
            return (int)(*(unsigned char*)valb);
        case glTFAccessorComponentType::Short: return (int)(*(short*)valb);
        case glTFAccessorComponentType::UnsignedShort:
            return (int)(*(unsigned short*)valb);
        case glTFAccessorComponentType::UnsignedInt:
            return (int)(*(unsigned int*)valb);
        case glTFAccessorComponentType::NotSet:
            throw std::runtime_error("bad enum value");
            break;
    }
    return 0;
}

int accessor_view::_num_components(glTFAccessorType type) {
    switch (type) {
        case glTFAccessorType::Scalar: return 1;
        case glTFAccessorType::Vec2: return 2;
        case glTFAccessorType::Vec3: return 3;
        case glTFAccessorType::Vec4: return 4;
        case glTFAccessorType::Mat2: return 4;
        case glTFAccessorType::Mat3: return 9;
        case glTFAccessorType::Mat4: return 16;
        default: assert(false); return 0;
    }
}

int accessor_view::_ctype_size(glTFAccessorComponentType componentType) {
    switch (componentType) {
        case glTFAccessorComponentType::Byte: return 1;
        case glTFAccessorComponentType::UnsignedByte: return 1;
        case glTFAccessorComponentType::Short: return 2;
        case glTFAccessorComponentType::UnsignedShort: return 2;
        case glTFAccessorComponentType::UnsignedInt: return 4;
        case glTFAccessorComponentType::Float: return 4;
        default: assert(false); return 0;
    }
}

}  // namespace ygl

#endif

#if YGL_OPENGL

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR OPENGL
// -----------------------------------------------------------------------------
namespace ygl {
// Checks for GL error and then prints
bool check_glerror(bool print) {
    auto ok = glGetError();
    if (ok == GL_NO_ERROR) return true;
    if (!print) return false;
    switch (ok) {
        case GL_NO_ERROR: printf("GL_NO_ERROR\n"); break;
        case GL_INVALID_ENUM: printf("GL_INVALID_ENUM\n"); break;
        case GL_INVALID_VALUE: printf("GL_INVALID_VALUE\n"); break;
        case GL_INVALID_OPERATION: printf("GL_INVALID_OPERATION\n"); break;
        case GL_INVALID_FRAMEBUFFER_OPERATION:
            printf("GL_INVALID_FRAMEBUFFER_OPERATION\n");
            break;
        case GL_OUT_OF_MEMORY: printf("GL_OUT_OF_MEMORY\n"); break;
        default: printf("<UNKNOWN GL ERROR>\n"); break;
    }
    return false;
}

// Clear window
void clear_glbuffers(const vec4f& background) {
    assert(check_glerror());
    glClearColor(background.x, background.y, background.z, background.w);
    glClearDepth(1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    assert(check_glerror());
}

// Enable/disable depth test
void enable_gldepth_test(bool enabled) {
    assert(check_glerror());
    if (enabled) {
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LEQUAL);
    } else {
        glDisable(GL_DEPTH_TEST);
    }
    assert(check_glerror());
}

// Enable/disable culling
void enable_glculling(bool enabled, bool front, bool back) {
    assert(check_glerror());
    if (enabled && (front || back)) {
        glEnable(GL_CULL_FACE);
        if (front && back)
            glCullFace(GL_FRONT_AND_BACK);
        else if (front)
            glCullFace(GL_FRONT);
        else if (back)
            glCullFace(GL_BACK);
    } else {
        glDisable(GL_CULL_FACE);
        glCullFace(GL_BACK);
    }
    assert(check_glerror());
}

// Enable/disable wireframe
void enable_glwireframe(bool enabled) {
    assert(check_glerror());
    if (enabled)
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    else
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    assert(check_glerror());
}

// Enable/disable blending
void enable_glblending(bool enabled) {
    assert(check_glerror());
    if (enabled) {
        glEnable(GL_BLEND);
    } else {
        glDisable(GL_BLEND);
    }
    assert(check_glerror());
}

// Set blending to over operator
void set_glblend_over() {
    assert(check_glerror());
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    assert(check_glerror());
}

// Line width
void gl_line_width(float w) {
    assert(check_glerror());
    glLineWidth(min(max(w, 0.0f), 1.0f));
    assert(check_glerror());
}

// Set viewport
void set_glviewport(const vec4i& v) {
    assert(check_glerror());
    glViewport(v.x, v.y, v.z, v.w);
    assert(check_glerror());
}

// Set viewport
void set_glviewport(const vec2i& v) {
    assert(check_glerror());
    glViewport(0, 0, v.x, v.y);
    assert(check_glerror());
}

// This is a public API. See above for documentation.
void read_glimagef(float* pixels, int w, int h, int nc) {
    assert(check_glerror());
    int formats[4] = {GL_RED, GL_RG, GL_RGB, GL_RGBA};
    glReadPixels(0, 0, w, h, formats[nc - 1], GL_FLOAT, pixels);
    assert(check_glerror());
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF OPENGL TEXTURE FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Implementation of update_texture.
void update_gltexture(gltexture& txt, int w, int h, int nc, const void* pixels,
    bool floats, bool linear, bool mipmap, bool as_float, bool as_srgb) {
    auto refresh = !txt.tid || txt.width != w || txt.height != h ||
                   txt.ncomp != nc || txt.as_float != as_float ||
                   txt.as_srgb != as_srgb || txt.mipmap != mipmap ||
                   txt.linear != linear;
    txt.width = w;
    txt.height = h;
    txt.ncomp = nc;
    txt.as_float = as_float;
    txt.as_srgb = as_srgb;
    txt.mipmap = mipmap;
    txt.linear = linear;
    assert(!as_srgb || !as_float);
    assert(check_glerror());
    if (w * h) {
        int formats_ub[4] = {GL_RED, GL_RG, GL_RGB, GL_RGBA};
        int formats_sub[4] = {GL_RED, GL_RG, GL_SRGB, GL_SRGB_ALPHA};
        int formats_f[4] = {GL_R32F, GL_RG32F, GL_RGB32F, GL_RGBA32F};
        int* formats =
            (as_float) ? formats_f : ((as_srgb) ? formats_sub : formats_ub);
        assert(check_glerror());
        if (!txt.tid) glGenTextures(1, &txt.tid);
        glBindTexture(GL_TEXTURE_2D, txt.tid);
        if (refresh) {
            glTexImage2D(GL_TEXTURE_2D, 0, formats[nc - 1], w, h, 0,
                formats_ub[nc - 1], (floats) ? GL_FLOAT : GL_UNSIGNED_BYTE,
                pixels);
        } else {
            glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, w, h, formats_ub[nc - 1],
                (floats) ? GL_FLOAT : GL_UNSIGNED_BYTE, pixels);
        }
        if (mipmap) glGenerateMipmap(GL_TEXTURE_2D);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
            (linear) ? GL_LINEAR : GL_NEAREST);
        if (mipmap) {
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
                (linear) ? GL_LINEAR_MIPMAP_LINEAR : GL_NEAREST_MIPMAP_NEAREST);
        } else {
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
                (linear) ? GL_LINEAR : GL_NEAREST);
        }
        glBindTexture(GL_TEXTURE_2D, 0);
    } else {
        if (txt.tid) {
            glBindTexture(GL_TEXTURE_2D, txt.tid);
            glDeleteTextures(1, &txt.tid);
            txt.tid = 0;
            glBindTexture(GL_TEXTURE_2D, 0);
        }
    }
    assert(check_glerror());
}

// Binds a texture to a texture unit
void bind_gltexture(const gltexture& txt, uint unit) {
    glActiveTexture(GL_TEXTURE0 + unit);
    glBindTexture(GL_TEXTURE_2D, txt.tid);
}

// Unbinds
void unbind_gltexture(const gltexture& txt, uint unit) {
    glActiveTexture(GL_TEXTURE0 + unit);
    glBindTexture(GL_TEXTURE_2D, 0);
}

// Destroys the texture tid.
void clear_gltexture(gltexture& txt) {
    assert(check_glerror());
    glDeleteTextures(1, &txt.tid);
    txt.tid = 0;
    assert(check_glerror());
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF OPENGL VERTEX/ELEMENTS BUFFERS
// -----------------------------------------------------------------------------
namespace ygl {

// Updates the buffer with new data.
void update_glbuffer(glvertex_buffer& buf, bool elems, int num, int ncomp,
    const void* values, bool dynamic) {
    auto resize =
        !buf.bid || num * ncomp != buf.num * buf.ncomp || elems != buf.elems;
    buf.num = num;
    buf.ncomp = ncomp;
    buf.elems = elems;
    auto target = (elems) ? GL_ELEMENT_ARRAY_BUFFER : GL_ARRAY_BUFFER;
    auto esize = (elems) ? sizeof(int) : sizeof(float);
    auto dmode = (dynamic) ? GL_DYNAMIC_DRAW : GL_STATIC_DRAW;
    assert(check_glerror());
    if (num) {
        if (!buf.bid) glGenBuffers(1, &buf.bid);
        glBindBuffer(target, buf.bid);
        if (resize) {
            glBufferData(target, buf.num * buf.ncomp * esize, values, dmode);
        } else {
            glBufferSubData(target, 0, buf.num * buf.ncomp * esize, values);
        }
        glBindBuffer(target, 0);
    } else {
        if (buf.bid) {
            glBindBuffer(target, buf.bid);
            glDeleteBuffers(1, &buf.bid);
            buf.bid = 0;
            glBindBuffer(target, 0);
        }
    }
    assert(check_glerror());
}

// Bind the buffer at a particular attribute location
void bind_glbuffer(const glvertex_buffer& buf, uint vattr) {
    if (buf.elems) throw std::runtime_error("Bad OpenGL buffer.");
    glEnableVertexAttribArray(vattr);
    glBindBuffer(GL_ARRAY_BUFFER, buf.bid);
    glVertexAttribPointer(vattr, buf.ncomp, GL_FLOAT, false, 0, 0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

// Unbind the buffer
void unbind_glbuffer(const glvertex_buffer& buf, uint vattr) {
    if (buf.elems) throw std::runtime_error("Bad OpenGL buffer.");
    glDisableVertexAttribArray(vattr);
}

// Draws elements.
void draw_glelems(const glvertex_buffer& buf) {
    if (!buf.elems) throw std::runtime_error("Bad OpenGL buffer.");
    if (!buf.bid) return;
    assert(check_glerror());
    int mode = 0;
    switch (buf.ncomp) {
        case 1: mode = GL_POINTS; break;
        case 2: mode = GL_LINES; break;
        case 3: mode = GL_TRIANGLES; break;
        case 4: mode = GL_QUADS; break;
        default: assert(false);
    };
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf.bid);
    glDrawElements(mode, buf.ncomp * buf.num, GL_UNSIGNED_INT, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    assert(check_glerror());
}

// Unbind the buffer
void unbind_glbuffer(uint vattr) { glDisableVertexAttribArray(vattr); }

// Destroys the buffer
void clear_glbuffer(glvertex_buffer& buf) {
    assert(check_glerror());
    glDeleteBuffers(1, &buf.bid);
    buf.bid = 0;
    assert(check_glerror());
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF OPENGL PROGRAM FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Creates and OpenGL program from vertex and fragment code. Returns the
// program id. Optionally return vertex and fragment shader ids. A VAO
// is created.
glprogram make_glprogram(
    const std::string& vertex, const std::string& fragment) {
    auto prog = glprogram();

    assert(check_glerror());
    glGenVertexArrays(1, &prog.vao);
    glBindVertexArray(prog.vao);
    assert(check_glerror());

    int errflags;
    char errbuf[10000];

    // create vertex
    prog.vid = glCreateShader(GL_VERTEX_SHADER);
    const char* vertex_str = vertex.c_str();
    glShaderSource(prog.vid, 1, &vertex_str, NULL);
    glCompileShader(prog.vid);
    glGetShaderiv(prog.vid, GL_COMPILE_STATUS, &errflags);
    if (!errflags) {
        glGetShaderInfoLog(prog.vid, 10000, 0, errbuf);
        throw std::runtime_error(
            std::string("shader not compiled\n\n") + errbuf);
    }
    assert(glGetError() == GL_NO_ERROR);

    // create fragment
    prog.fid = glCreateShader(GL_FRAGMENT_SHADER);
    const char* fragment_str = fragment.c_str();
    glShaderSource(prog.fid, 1, &fragment_str, NULL);
    glCompileShader(prog.fid);
    glGetShaderiv(prog.fid, GL_COMPILE_STATUS, &errflags);
    if (!errflags) {
        glGetShaderInfoLog(prog.fid, 10000, 0, errbuf);
        throw std::runtime_error(
            std::string("shader not compiled\n\n") + errbuf);
    }
    assert(glGetError() == GL_NO_ERROR);

    // create program
    prog.pid = glCreateProgram();
    glAttachShader(prog.pid, prog.vid);
    glAttachShader(prog.pid, prog.fid);
    glLinkProgram(prog.pid);
    glValidateProgram(prog.pid);
    glGetProgramiv(prog.pid, GL_LINK_STATUS, &errflags);
    if (!errflags) {
        glGetProgramInfoLog(prog.pid, 10000, 0, errbuf);
        throw std::runtime_error(
            std::string("program not linked\n\n") + errbuf);
    }
    glGetProgramiv(prog.pid, GL_VALIDATE_STATUS, &errflags);
    if (!errflags) {
        glGetProgramInfoLog(prog.pid, 10000, 0, errbuf);
        throw std::runtime_error(
            std::string("program not linked\n\n") + errbuf);
    }
    assert(check_glerror());

    glBindVertexArray(0);
    assert(check_glerror());

    return prog;
}

// Destroys the program pid and optionally the sahders vid and fid.
void clear_program(glprogram& prog) {
    assert(check_glerror());
    glDetachShader(prog.pid, prog.vid);
    glDeleteShader(prog.vid);
    prog.vid = 0;
    glDetachShader(prog.pid, prog.fid);
    glDeleteShader(prog.fid);
    prog.fid = 0;
    glDeleteProgram(prog.pid);
    prog.pid = 0;
    glDeleteVertexArrays(1, &prog.vao);
    prog.vao = 0;
    assert(check_glerror());
}

// Get uniform location (simple GL wrapper that avoids GL includes)
int get_gluniform_location(const glprogram& prog, const std::string& name) {
    assert(check_glerror());
    return glGetUniformLocation(prog.pid, name.c_str());
    assert(check_glerror());
}

// Get uniform location (simple GL wrapper that avoids GL includes)
int get_glattrib_location(const glprogram& prog, const std::string& name) {
    assert(check_glerror());
    return glGetAttribLocation(prog.pid, name.c_str());
    assert(check_glerror());
}

// Set uniform values.
void set_gluniform(const glprogram& prog, int pos, bool val) {
    if (pos < 0) throw std::runtime_error("bad OpenGL id");
    glUniform1i(pos, (int)val);
}
void set_gluniform(const glprogram& prog, int pos, int val) {
    if (pos < 0) throw std::runtime_error("bad OpenGL id");
    glUniform1i(pos, val);
}
void set_gluniform(const glprogram& prog, int pos, float val) {
    if (pos < 0) throw std::runtime_error("bad OpenGL id");
    glUniform1f(pos, val);
}
void set_gluniform(const glprogram& prog, int pos, const vec2f& val) {
    if (pos < 0) throw std::runtime_error("bad OpenGL id");
    glUniform2f(pos, val.x, val.y);
}
void set_gluniform(const glprogram& prog, int pos, const vec3f& val) {
    if (pos < 0) throw std::runtime_error("bad OpenGL id");
    glUniform3f(pos, val.x, val.y, val.z);
}
void set_gluniform(const glprogram& prog, int pos, const vec4f& val) {
    if (pos < 0) throw std::runtime_error("bad OpenGL id");
    glUniform4f(pos, val.x, val.y, val.z, val.w);
}
void set_gluniform(const glprogram& prog, int pos, const mat4f& val) {
    if (pos < 0) throw std::runtime_error("bad OpenGL id");
    glUniformMatrix4fv(pos, 1, false, &val.x.x);
}
void set_gluniform(const glprogram& prog, int pos, const frame3f& val) {
    if (pos < 0) throw std::runtime_error("bad OpenGL id");
    glUniformMatrix4x3fv(pos, 1, false, &val.x.x);
}

// Set uniform texture id tid and unit tunit for program pid and
// variable var.
void set_gluniform_texture(
    const glprogram& prog, int pos, const gltexture_info& tinfo, uint tunit) {
    static const auto wrap_mode_map =
        std::map<gltexture_wrap, uint>{{gltexture_wrap::repeat, GL_REPEAT},
            {gltexture_wrap::clamp, GL_CLAMP_TO_EDGE},
            {gltexture_wrap::mirror, GL_MIRRORED_REPEAT}};
    static const auto filter_mode_map = std::map<gltexture_filter, uint>{
        {gltexture_filter::nearest, GL_NEAREST},
        {gltexture_filter::linear, GL_LINEAR},
        {gltexture_filter::nearest_mipmap_nearest, GL_NEAREST_MIPMAP_NEAREST},
        {gltexture_filter::linear_mipmap_nearest, GL_LINEAR_MIPMAP_NEAREST},
        {gltexture_filter::nearest_mipmap_linear, GL_NEAREST_MIPMAP_LINEAR},
        {gltexture_filter::linear_mipmap_linear, GL_LINEAR_MIPMAP_LINEAR}};

    assert(check_glerror());
    if (pos < 0) throw std::runtime_error("bad OpenGL id");
    if (is_gltexture_valid(tinfo.txt)) {
        bind_gltexture(tinfo.txt, tunit);
        if (tinfo.wrap_s != gltexture_wrap::not_set)
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,
                wrap_mode_map.at(tinfo.wrap_s));
        if (tinfo.wrap_t != gltexture_wrap::not_set)
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T,
                wrap_mode_map.at(tinfo.wrap_t));
        if (tinfo.filter_min != gltexture_filter::not_set)
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
                filter_mode_map.at(tinfo.filter_min));
        if (tinfo.filter_mag != gltexture_filter::not_set)
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
                filter_mode_map.at(tinfo.filter_mag));
        glUniform1i(pos, tunit);
    } else {
        unbind_gltexture(tinfo.txt, tunit);
        glUniform1i(pos, tunit);
    }
    assert(check_glerror());
}

// Set uniform texture id tid and unit tunit for program pid and
// variable var.
void set_gluniform_texture(
    const glprogram& prog, int pos, const gltexture& txt, uint tunit) {
    assert(check_glerror());
    if (pos < 0) throw std::runtime_error("bad OpenGL id");
    if (is_gltexture_valid(txt)) {
        bind_gltexture(txt, tunit);
        glUniform1i(pos, tunit);
    } else {
        unbind_gltexture(txt, tunit);
        glUniform1i(pos, tunit);
    }
    assert(check_glerror());
}

// Set uniform texture id tid and unit tunit for program pid and
// variable var.
void set_gluniform_texture(const glprogram& prog, const std::string& name,
    const gltexture& txt, uint tunit) {
    set_gluniform_texture(prog, get_gluniform_location(prog, name), txt, tunit);
}

// Sets a vartex attribute for program pid and variable var.
void set_glattribute(
    const glprogram& prog, int pos, const glvertex_buffer& buf, float def) {
    assert(check_glerror());
    if (pos < 0) throw std::runtime_error("bad OpenGL id");
    if (is_glbuffer_valid(buf)) {
        glEnableVertexAttribArray(pos);
        bind_glbuffer(buf, pos);
    } else {
        glDisableVertexAttribArray(pos);
        unbind_glbuffer(buf, pos);
        glVertexAttrib1f(pos, def);
    }
    assert(check_glerror());
}
// Sets a vartex attribute for program pid and variable var.
void set_glattribute(const glprogram& prog, int pos, const glvertex_buffer& buf,
    const vec2f& def) {
    assert(check_glerror());
    if (pos < 0) throw std::runtime_error("bad OpenGL id");
    if (is_glbuffer_valid(buf)) {
        glEnableVertexAttribArray(pos);
        bind_glbuffer(buf, pos);
    } else {
        glDisableVertexAttribArray(pos);
        unbind_glbuffer(buf, pos);
        glVertexAttrib2f(pos, def.x, def.y);
    }
    assert(check_glerror());
}
// Sets a vartex attribute for program pid and variable var.
void set_glattribute(const glprogram& prog, int pos, const glvertex_buffer& buf,
    const vec3f& def) {
    assert(check_glerror());
    if (pos < 0) throw std::runtime_error("bad OpenGL id");
    if (is_glbuffer_valid(buf)) {
        glEnableVertexAttribArray(pos);
        bind_glbuffer(buf, pos);
    } else {
        glDisableVertexAttribArray(pos);
        unbind_glbuffer(buf, pos);
        glVertexAttrib3f(pos, def.x, def.y, def.z);
    }
    assert(check_glerror());
}
// Sets a vartex attribute for program pid and variable var.
void set_glattribute(const glprogram& prog, int pos, const glvertex_buffer& buf,
    const vec4f& def) {
    assert(check_glerror());
    if (pos < 0) throw std::runtime_error("bad OpenGL id");
    if (is_glbuffer_valid(buf)) {
        glEnableVertexAttribArray(pos);
        bind_glbuffer(buf, pos);
    } else {
        glDisableVertexAttribArray(pos);
        unbind_glbuffer(buf, pos);
        glVertexAttrib4f(pos, def.x, def.y, def.z, def.w);
    }
    assert(check_glerror());
}

// Binds a program
void bind_glprogram(const glprogram& prog) {
    assert(check_glerror());
    if (!prog.pid) return;
    glBindVertexArray(prog.vao);
    glUseProgram(prog.pid);
    assert(check_glerror());
}

// Unbind a program
void unbind_glprogram(const glprogram& prog) {
    assert(check_glerror());
    glUseProgram(0);
    glBindVertexArray(0);
    assert(check_glerror());
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SCENE SHADER FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Init shading
void update_gltexture(const texture* txt, gltexture& gtxt) {
    if (!txt) {
        clear_gltexture(gtxt);
    } else {
        if (!txt->hdr.pixels.empty()) {
            update_gltexture(gtxt, txt->hdr, true, true, true);
        } else if (!txt->ldr.pixels.empty()) {
            update_gltexture(gtxt, txt->ldr, true, true, true);
        } else
            assert(false);
    }
}

template <typename T>
std::vector<std::vector<T>> _split_elems(
    int ngroups, const std::vector<T>& elems, const std::vector<int>& ids) {
    if (ids.empty()) return {elems};
    auto splits = std::vector<std::vector<T>>(ngroups);
    if (elems.empty()) return splits;
    for (auto i = 0; i < elems.size(); i++) {
        splits[ids[i]].push_back(elems[i]);
    }
    return splits;
}

// Update shading
void update_glshape(const shape* shp, glshape& gshp) {
    auto update_vert_buffer = [](auto& buf, const auto& vert) {
        if (vert.empty()) {
            clear_glbuffer(buf);
        } else {
            update_glbuffer(buf, false, vert);
        }
    };
    auto update_elem_buffer = [](auto& buf, const auto& elem) {
        if (elem.empty()) {
            clear_glbuffer(buf);
        } else {
            update_glbuffer(buf, true, elem);
        }
    };
    if (!shp) {
        clear_glshape(gshp);
    } else {
        if (!shp->quads_pos.empty()) {
            auto pos = std::vector<vec3f>();
            auto norm = std::vector<vec3f>();
            auto texcoord = std::vector<vec2f>();
            auto quads = std::vector<vec4i>();
            std::tie(quads, pos, norm, texcoord) =
                convert_face_varying(shp->quads_pos, shp->quads_norm,
                    shp->quads_texcoord, shp->pos, shp->norm, shp->texcoord);
            update_vert_buffer(gshp.pos, pos);
            update_vert_buffer(gshp.norm, norm);
            update_vert_buffer(gshp.texcoord, texcoord);
            update_vert_buffer(gshp.color, std::vector<vec4f>{});
            update_vert_buffer(gshp.tangsp, std::vector<vec4f>{});
            update_elem_buffer(gshp.quads, convert_quads_to_triangles(quads));
            update_elem_buffer(gshp.edges, get_edges({}, {}, quads));
        } else {
            update_vert_buffer(gshp.pos, shp->pos);
            update_vert_buffer(gshp.norm, shp->norm);
            update_vert_buffer(gshp.texcoord, shp->texcoord);
            update_vert_buffer(gshp.color, shp->color);
            update_vert_buffer(gshp.tangsp, shp->tangsp);
            update_elem_buffer(gshp.points, shp->points);
            update_elem_buffer(gshp.lines, shp->lines);
            update_elem_buffer(gshp.triangles, shp->triangles);
            update_elem_buffer(
                gshp.quads, convert_quads_to_triangles(shp->quads));
            update_elem_buffer(
                gshp.beziers, convert_bezier_to_lines(shp->beziers));
            update_elem_buffer(
                gshp.edges, get_edges({}, shp->triangles, shp->quads));
        }
    }
}

// clear glshape
void clear_glshape(glshape& gshp) {
    clear_glbuffer(gshp.pos);
    clear_glbuffer(gshp.norm);
    clear_glbuffer(gshp.texcoord);
    clear_glbuffer(gshp.color);
    clear_glbuffer(gshp.tangsp);
    clear_glbuffer(gshp.points);
    clear_glbuffer(gshp.lines);
    clear_glbuffer(gshp.triangles);
    clear_glbuffer(gshp.quads);
    clear_glbuffer(gshp.beziers);
    clear_glbuffer(gshp.edges);
}

// Add gl lights
void add_gllights(gllights& lights, const instance* ist) {
    if (!ist->mat || ist->mat->ke == zero3f) return;
    if (lights.pos.size() >= 16) return;
    auto shp = ist->shp;
    if (!shp->points.empty()) {
        for (auto p : shp->points) {
            if (lights.pos.size() >= 16) break;
            lights.pos.push_back(transform_point(ist->frame, shp->pos[p]));
            lights.ke.push_back(ist->mat->ke);
            lights.type.push_back(gllight_type::point);
        }
    } else {
        auto bbox = compute_bbox(shp);
        auto pos = (bbox.max + bbox.min) / 2;
        auto area = 0.0f;
        for (auto l : shp->lines)
            area += line_length(shp->pos[l.x], shp->pos[l.y]);
        for (auto t : shp->triangles)
            area += triangle_area(shp->pos[t.x], shp->pos[t.y], shp->pos[t.z]);
        for (auto t : shp->quads)
            area += quad_area(
                shp->pos[t.x], shp->pos[t.y], shp->pos[t.z], shp->pos[t.w]);
        auto ke = ist->mat->ke * area;
        if (lights.pos.size() < 16) {
            lights.pos.push_back(transform_point(ist->frame, pos));
            lights.ke.push_back(ke);
            lights.type.push_back(gllight_type::point);
        }
    }
}

// Initialize gl lights
gllights make_gllights(const scene* scn) {
    auto lights = gllights();
    for (auto ist : scn->instances) add_gllights(lights, ist);
    return lights;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR OPRNGL IMAGE SHADER FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Initialize the program. Call with true only after the GL is initialized.
glimage_program make_glimage_program() {
    std::string _header =
        R"(
        #version 330

        float pi = 3.14159265;

        uniform vec2 offset;
        uniform vec2 win_size;
        uniform float zoom;

        uniform sampler2D img;

        )";

    std::string _vert =
        R"(
        layout(location = 0) in vec2 vert_texcoord;

        out vec2 texcoord;

        void main() {
            vec2 size = textureSize(img, 0).xy;
            texcoord = vert_texcoord.xy;
            vec2 pos = offset + size * vert_texcoord.xy * zoom;
            vec2 upos = 2 * pos / win_size - vec2(1,1);
            upos.y = - upos.y;
            gl_Position = vec4(upos.x, upos.y, 0, 1);
        }

        )";

    std::string _frag_tonemap =
        R"(
        struct Tonemap {
            float exposure;
            int type;
        };
        uniform Tonemap tonemap;

        vec3 eval_gamma(vec3 x) {
            return pow(x, vec3(1/2.2));
        }

        vec3 eval_srgb(vec3 x) {
            float r = (x.x < 0.0031308) ? 12.92 * x.x : (1+0.055) * pow(x.x, 1/2.4) - 0.055;
            float g = (x.y < 0.0031308) ? 12.92 * x.y : (1+0.055) * pow(x.y, 1/2.4) - 0.055;
            float b = (x.z < 0.0031308) ? 12.92 * x.z : (1+0.055) * pow(x.z, 1/2.4) - 0.055;
            return vec3(r,g,b);
        }

        vec3 eval_filmic1(vec3 x) {
            // http://filmicworlds.com/blog/filmic-tonemapping-operators/
            x = max(vec3(0),x-0.004);
            return (x*(6.2*x+.5))/(x*(6.2*x+1.7)+0.06);
        }

        vec3 eval_filmic2(vec3 x) {
            // https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
            // x *= 0.6; // brings it back to ACES range
            return pow(clamp((x*(2.51f*x+0.03f))/(x*(2.43f*x+0.59f)+0.14f),0,1),vec3(1/2.2));
        }

        vec3 eval_filmic3(vec3 x) {
            // https://github.com/TheRealMJP/BakingLab/blob/master/BakingLab/ACES.hlsl

            // sRGB => XYZ => D65_2_D60 => AP1 => RRT_SAT
            mat3 ACESInputMat = transpose(mat3(
                vec3(0.59719, 0.35458, 0.04823),
                vec3(0.07600, 0.90834, 0.01566),
                vec3(0.02840, 0.13383, 0.83777)));

            // ODT_SAT => XYZ => D60_2_D65 => sRGB
            mat3 ACESOutputMat = transpose(mat3(
                vec3( 1.60475, -0.53108, -0.07367),
                vec3(-0.10208,  1.10813, -0.00605),
                vec3(-0.00327, -0.07276,  1.07602)));

            x = 2 * x; // matches standard range
            x = ACESInputMat * x;
            // Apply RRT and ODT
            vec3 a = x * (x + 0.0245786f) - 0.000090537f;
            vec3 b = x * (0.983729f * x + 0.4329510f) + 0.238081f;
            x = a / b;
            x = ACESOutputMat * x;
            x = pow(clamp(x,0,1),vec3(1/2.2));
            return x;
        }

        vec3 eval_tonemap(vec3 c) {
            // final color correction
            c = c*pow(2,tonemap.exposure);
            switch(tonemap.type) {
                case 0: break;
                case 1: c = eval_gamma(c); break;
                case 2: c = eval_srgb(c); break;
                case 3: c = eval_filmic1(c); break;
                case 4: c = eval_filmic2(c); break;
                case 5: c = eval_filmic3(c); break;
                default: c = eval_gamma(c); break;
            }
            return c;
        }

        )";

    std::string _frag_main =
        R"(
        in vec2 texcoord;
        out vec4 color;

        void main() {
            vec4 c = texture(img,texcoord);
            c.xyz = eval_tonemap(c.xyz);
            color = c;
        }
        )";

    auto prog = glimage_program();
    prog.prog =
        make_glprogram(_header + _vert, _header + _frag_tonemap + _frag_main);

    update_glbuffer(
        prog.vbo, false, std::vector<vec2f>{{0, 0}, {0, 1}, {1, 1}, {1, 0}});
    update_glbuffer(prog.ebo, true, std::vector<vec3i>{{0, 1, 2}, {0, 2, 3}});
    return prog;
}

// Draws the stdimage program.
void draw_glimage(const glimage_program& prog, const gltexture& txt,
    const vec2i& win_size, const vec2f& offset, float zoom,
    tonemap_type tonemap, float exposure) {
    assert(is_gltexture_valid(txt));

    bind_glprogram(prog.prog);

    enable_glblending(true);
    set_glblend_over();

    bind_gltexture(txt, 0);
    set_gluniform(prog.prog, "zoom", zoom);
    set_gluniform(
        prog.prog, "win_size", vec2f{(float)win_size.x, (float)win_size.y});
    set_gluniform(prog.prog, "offset", offset);
    set_gluniform(prog.prog, "tonemap.exposure", exposure);
    set_gluniform(prog.prog, "tonemap.type", (int)tonemap);
    set_gluniform_texture(prog.prog, "img", txt, 0);

    set_glattribute(prog.prog, "vert_texcoord", prog.vbo, vec2f{0, 0});
    draw_glelems(prog.ebo);

    unbind_glprogram(prog.prog);

    enable_glblending(false);

    assert(check_glerror());
}

// Computes the image uv coordinates corresponding to the view parameters.
vec2i get_draw_image_coords(
    const vec2f& mouse_pos, const glimage_params& params) {
    auto xy = (mouse_pos - params.offset) / params.zoom;
    return {(int)round(xy.x), (int)round(xy.y)};
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR OPRNGL STANDARD SURFACE SHADER FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Initialize a standard shader. Call with true only after the gl has
// been initialized
glsurface_program make_glsurface_program() {
#ifndef _WIN32
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverlength-strings"
#endif
    std::string _vert =
        R"(
        #version 330

        layout(location = 0) in vec3 vert_pos;            // vertex position (in mesh coordinate frame)
        layout(location = 1) in vec3 vert_norm;           // vertex normal (in mesh coordinate frame)
        layout(location = 2) in vec2 vert_texcoord;       // vertex texcoords
        layout(location = 3) in vec4 vert_color;          // vertex color
        layout(location = 4) in vec4 vert_tangsp;         // vertex tangent space

        uniform mat4 shape_xform;           // shape transform
        uniform float shape_normal_offset;           // shape normal offset

        uniform mat4 cam_xform;          // camera xform
        uniform mat4 cam_xform_inv;      // inverse of the camera frame (as a matrix)
        uniform mat4 cam_proj;           // camera projection

        out vec3 pos;                   // [to fragment shader] vertex position (in world coordinate)
        out vec3 norm;                  // [to fragment shader] vertex normal (in world coordinate)
        out vec2 texcoord;              // [to fragment shader] vertex texture coordinates
        out vec4 color;                 // [to fragment shader] vertex color
        out vec4 tangsp;                // [to fragment shader] vertex tangent space

        // main function
        void main() {
            // copy values
            pos = vert_pos;
            norm = vert_norm;
            tangsp = vert_tangsp;

            // normal offset
            if(shape_normal_offset != 0) {
                pos += shape_normal_offset * norm;
            }

            // world projection
            pos = (shape_xform * vec4(pos,1)).xyz;
            norm = (shape_xform * vec4(norm,0)).xyz;
            tangsp.xyz = (shape_xform * vec4(tangsp.xyz,0)).xyz;

            // copy other vertex properties
            texcoord = vert_texcoord;
            color = vert_color;

            // clip
            gl_Position = cam_proj * cam_xform_inv * vec4(pos,1);
        }
        )";

    std::string _frag_header =
        R"(
        #version 330

        float pi = 3.14159265;

        )";

    std::string _frag_tonemap =
        R"(
        uniform float tm_exposure; 
        uniform int tm_type;

        vec3 eval_gamma(vec3 x) {
            return pow(x, vec3(1/2.2));
        }

        vec3 eval_srgb(vec3 x) {
            float r = (x.x < 0.0031308) ? 12.92 * x.x : (1+0.055) * pow(x.x, 1/2.4) - 0.055;
            float g = (x.y < 0.0031308) ? 12.92 * x.y : (1+0.055) * pow(x.y, 1/2.4) - 0.055;
            float b = (x.z < 0.0031308) ? 12.92 * x.z : (1+0.055) * pow(x.z, 1/2.4) - 0.055;
            return vec3(r,g,b);
        }

        vec3 eval_filmic1(vec3 x) {
            // http://filmicworlds.com/blog/filmic-tonemapping-operators/
            x = max(vec3(0),x-0.004);
            return (x*(6.2*x+.5))/(x*(6.2*x+1.7)+0.06);
        }

        vec3 eval_filmic2(vec3 x) {
            // https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
            // x *= 0.6; // brings it back to ACES range
            return pow(clamp((x*(2.51f*x+0.03f))/(x*(2.43f*x+0.59f)+0.14f),0,1),vec3(1/2.2));
        }

        vec3 eval_filmic3(vec3 x) {
            // https://github.com/TheRealMJP/BakingLab/blob/master/BakingLab/ACES.hlsl

            // sRGB => XYZ => D65_2_D60 => AP1 => RRT_SAT
            mat3 ACESInputMat = transpose(mat3(
                vec3(0.59719, 0.35458, 0.04823),
                vec3(0.07600, 0.90834, 0.01566),
                vec3(0.02840, 0.13383, 0.83777)));

            // ODT_SAT => XYZ => D60_2_D65 => sRGB
            mat3 ACESOutputMat = transpose(mat3(
                vec3( 1.60475, -0.53108, -0.07367),
                vec3(-0.10208,  1.10813, -0.00605),
                vec3(-0.00327, -0.07276,  1.07602)));

            x = 2 * x; // matches standard range
            x = ACESInputMat * x;
            // Apply RRT and ODT
            vec3 a = x * (x + 0.0245786f) - 0.000090537f;
            vec3 b = x * (0.983729f * x + 0.4329510f) + 0.238081f;
            x = a / b;
            x = ACESOutputMat * x;
            x = pow(clamp(x,0,1),vec3(1/2.2));
            return x;
        }

        vec3 eval_tonemap(vec3 c) {
            // final color correction
            c = c*pow(2,tm_exposure);
            switch(tm_type) {
                case 0: break;
                case 1: c = eval_gamma(c); break;
                case 2: c = eval_srgb(c); break;
                case 3: c = eval_filmic1(c); break;
                case 4: c = eval_filmic2(c); break;
                case 5: c = eval_filmic3(c); break;
                default: c = eval_gamma(c); break;
            }
            return c;
        }

        )";

    std::string _frag_lighting =
        R"(
        uniform bool eyelight;        // eyelight shading
        uniform vec3 lamb;             // ambient light
        uniform int lnum;              // number of lights
        uniform int ltype[16];         // light type (0 -> point, 1 -> directional)
        uniform vec3 lpos[16];         // light positions
        uniform vec3 lke[16];          // light intensities

        void eval_light(int lid, vec3 pos, out vec3 cl, out vec3 wi) {
            cl = vec3(0,0,0);
            wi = vec3(0,0,0);
            if(ltype[lid] == 0) {
                // compute point light color at pos
                cl = lke[lid] / pow(length(lpos[lid]-pos),2);
                // compute light direction at pos
                wi = normalize(lpos[lid]-pos);
            }
            else if(ltype[lid] == 1) {
                // compute light color
                cl = lke[lid];
                // compute light direction
                wi = normalize(lpos[lid]);
            }
        }

        )";

    std::string _frag_brdf =
        R"(
        vec3 brdfcos(int etype, vec3 ke, vec3 kd, vec3 ks, float rs, float op,
            vec3 n, vec3 wi, vec3 wo) {
            if(etype == 0) return vec3(0);
            vec3 wh = normalize(wi+wo);
            float ns = 2/(rs*rs)-2;
            float ndi = dot(wi,n), ndo = dot(wo,n), ndh = dot(wh,n);
            if(etype == 1) {
                return ((1+dot(wo,wi))/2) * kd/pi;
            } else if(etype == 2) {
                float si = sqrt(1-ndi*ndi);
                float so = sqrt(1-ndo*ndo);
                float sh = sqrt(1-ndh*ndh);
                if(si <= 0) return vec3(0);
                vec3 diff = si * kd / pi;
                if(sh<=0) return diff;
                float d = ((2+ns)/(2*pi)) * pow(si,ns);
                vec3 spec = si * ks * d / (4*si*so);
                return diff+spec;
            } else if(etype == 3 || etype == 4) {
                if(ndi<=0 || ndo <=0) return vec3(0);
                vec3 diff = ndi * kd / pi;
                if(ndh<=0) return diff;
                if(etype == 4) {
                    float d = ((2+ns)/(2*pi)) * pow(ndh,ns);
                    vec3 spec = ndi * ks * d / (4*ndi*ndo);
                    return diff+spec;
                } else {
                    float cos2 = ndh * ndh;
                    float tan2 = (1 - cos2) / cos2;
                    float alpha2 = rs * rs;
                    float d = alpha2 / (pi * cos2 * cos2 * (alpha2 + tan2) * (alpha2 + tan2));
                    float lambda_o = (-1 + sqrt(1 + (1 - ndo * ndo) / (ndo * ndo))) / 2;
                    float lambda_i = (-1 + sqrt(1 + (1 - ndi * ndi) / (ndi * ndi))) / 2;
                    float g = 1 / (1 + lambda_o + lambda_i);
                    vec3 spec = ndi * ks * d * g / (4*ndi*ndo);
                    return diff+spec;
                }
            }
        }

        )";

    std::string _frag_material =
        R"(
        uniform int elem_type;
        uniform bool elem_faceted;
        uniform vec4 highlight;   // highlighted color

        uniform int mat_type;          // material type
        uniform vec3 mat_ke;           // material ke
        uniform vec3 mat_kd;           // material kd
        uniform vec3 mat_ks;           // material ks
        uniform float mat_rs;          // material rs
        uniform float mat_op;          // material op

        uniform bool mat_txt_ke_on;    // material ke texture on
        uniform sampler2D mat_txt_ke;  // material ke texture
        uniform bool mat_txt_kd_on;    // material kd texture on
        uniform sampler2D mat_txt_kd;  // material kd texture
        uniform bool mat_txt_ks_on;    // material ks texture on
        uniform sampler2D mat_txt_ks;  // material ks texture
        uniform bool mat_txt_rs_on;    // material rs texture on
        uniform sampler2D mat_txt_rs;  // material rs texture

        uniform bool mat_txt_norm_on;    // material norm texture on
        uniform sampler2D mat_txt_norm;  // material norm texture
        uniform float mat_txt_norm_scale;  // material norm scale

        uniform bool mat_txt_occ_on;    // material occ texture on
        uniform sampler2D mat_txt_occ;  // material occ texture
        uniform float mat_txt_occ_scale;  // material occ scale

        uniform bool mat_use_phong;       // material use phong
        uniform bool mat_double_sided;    // material double sided
        uniform bool mat_alpha_cutout;    // material alpha cutout

        bool eval_material(vec2 texcoord, vec4 color, out vec3 ke, 
                           out vec3 kd, out vec3 ks, out float rs, out float op, out bool cutout) {
            if(mat_type == 0) {
                ke = mat_ke;
                kd = vec3(0,0,0);
                ks = vec3(0,0,0);
                op = 1;
                return false;
            }

            ke = color.xyz * mat_ke;
            kd = color.xyz * mat_kd;
            ks = color.xyz * mat_ks;
            rs = mat_rs;
            op = color.w * mat_op;

            vec3 ke_txt = (mat_txt_ke_on) ? texture(mat_txt_ke,texcoord).xyz : vec3(1);
            vec4 kd_txt = (mat_txt_kd_on) ? texture(mat_txt_kd,texcoord) : vec4(1);
            vec4 ks_txt = (mat_txt_ks_on) ? texture(mat_txt_ks,texcoord) : vec4(1);
            float rs_txt = (mat_txt_rs_on) ? texture(mat_txt_rs,texcoord).x : 1;

            // scale common values
            ke *= ke_txt;

            // get material color from textures and adjust values
            if(mat_type == 0) {
            } else if(mat_type == 1) {
                kd *= kd_txt.xyz;
                ks *= ks_txt.xyz;
                rs *= rs_txt;
                rs = rs*rs;
            } else if(mat_type == 2) {
                vec3 kb = kd * kd_txt.xyz;
                float km = ks.x * ks_txt.z;
                kd = kb * (1 - km);
                ks = kb * km + vec3(0.04) * (1 - km);
                rs *= ks_txt.y;
                rs = rs*rs;
                op *= kd_txt.w;
            } else if(mat_type == 3) {
                kd *= kd_txt.xyz;
                ks *= ks_txt.xyz;
                rs *= ks_txt.w;
                rs = (1 - rs) * (1 - rs);
                op *= kd_txt.w;
            }

            cutout = mat_alpha_cutout && op == 0;
            return true;
        }

        vec3 apply_normal_map(vec2 texcoord, vec3 norm, vec4 tangsp) {
            if(!mat_txt_norm_on) return norm;
            vec3 tangu = normalize(tangsp.xyz);
            vec3 tangv = normalize(cross(tangu, norm));
            if(tangsp.w < 0) tangv = -tangv;
            vec3 txt = 2 * pow(texture(mat_txt_norm,texcoord).xyz, vec3(1/2.2)) - 1;
            return normalize( tangu * txt.x + tangv * txt.y + norm * txt.z );
        }

        )";

    std::string _frag_main =
        R"(
        in vec3 pos;                   // [from vertex shader] position in world space
        in vec3 norm;                  // [from vertex shader] normal in world space (need normalization)
        in vec2 texcoord;              // [from vertex shader] texcoord
        in vec4 color;                 // [from vertex shader] color
        in vec4 tangsp;                // [from vertex shader] tangent space

        uniform mat4 cam_xform;          // camera xform
        uniform mat4 cam_xform_inv;      // inverse of the camera frame (as a matrix)
        uniform mat4 cam_proj;           // camera projection

        out vec4 frag_color;        // eyelight shading

        vec3 triangle_normal(vec3 pos) {
            vec3 fdx = dFdx(pos); 
            vec3 fdy = dFdy(pos); 
            return normalize(cross(fdx, fdy));
        }

        // main
        void main() {
            // view vector
            vec3 wo = normalize( (cam_xform*vec4(0,0,0,1)).xyz - pos );

            // prepare normals
            vec3 n;
            if(elem_faceted) {
                n = triangle_normal(pos);
            } else {
                n = normalize(norm);
            }

            // apply normal map
            n = apply_normal_map(texcoord, n, tangsp);

            // use faceforward to ensure the normals points toward us
            if(mat_double_sided) n = faceforward(n,-wo,n);

            // get material color from textures
            vec3 brdf_ke, brdf_kd, brdf_ks; float brdf_rs, brdf_op; bool brdf_cutout;
            bool has_brdf = eval_material(texcoord, color, brdf_ke, brdf_kd, brdf_ks, brdf_rs, brdf_op, brdf_cutout);

            // exit if needed
            if(brdf_cutout) discard;

            // check const color
            if(elem_type == 0) {
                frag_color = vec4(brdf_ke,brdf_op);
                return;
            }

            // emission
            vec3 c = brdf_ke;

            // check early exit
            if(brdf_kd != vec3(0,0,0) || brdf_ks != vec3(0,0,0)) {
                // eyelight shading
                if(eyelight) {
                    vec3 wi = wo;
                    c += pi * brdfcos((has_brdf) ? elem_type : 0, brdf_ke, brdf_kd, brdf_ks, brdf_rs, brdf_op, n,wi,wo);
                } else {
                    // accumulate ambient
                    c += lamb * brdf_kd;
                    // foreach light
                    for(int lid = 0; lid < lnum; lid ++) {
                        vec3 cl = vec3(0,0,0); vec3 wi = vec3(0,0,0);
                        eval_light(lid, pos, cl, wi);
                        c += cl * brdfcos((has_brdf) ? elem_type : 0, brdf_ke, brdf_kd, brdf_ks, brdf_rs, brdf_op, n,wi,wo);
                    }
                }
            }

            // final color correction
            c = eval_tonemap(c);

            // highlighting
            if(highlight.w > 0) {
                if(mod(int(gl_FragCoord.x)/4 + int(gl_FragCoord.y)/4, 2)  == 0)
                    c = highlight.xyz * highlight.w + c * (1-highlight.w);
            }

            // output final color by setting gl_FragColor
            frag_color = vec4(c,brdf_op);
        }
        )";
#ifndef _WIN32
#pragma GCC diagnostic pop
#endif

    assert(check_glerror());
    auto prog = glsurface_program();
    prog.prog =
        make_glprogram(_vert, _frag_header + _frag_tonemap + _frag_lighting +
                                  _frag_brdf + _frag_material + _frag_main);
    assert(check_glerror());
    prog.eyelight_id = get_gluniform_location(prog.prog, "eyelight");
    prog.exposure_id = get_gluniform_location(prog.prog, "tm_exposure");
    prog.tonemap_id = get_gluniform_location(prog.prog, "tm_type");
    prog.cam_xform_id = get_gluniform_location(prog.prog, "cam_xform");
    prog.cam_xform_inv_id = get_gluniform_location(prog.prog, "cam_xform_inv");
    prog.cam_proj_id = get_gluniform_location(prog.prog, "cam_proj");
    prog.lamb_id = get_gluniform_location(prog.prog, "lamb");
    prog.lnum_id = get_gluniform_location(prog.prog, "lnum");
    for (auto i = 0; i < 16; i++) {
        auto is = std::to_string(i);
        prog.lpos_id[i] = get_gluniform_location(prog.prog, "lpos[" + is + "]");
        prog.lke_id[i] = get_gluniform_location(prog.prog, "lke[" + is + "]");
        prog.ltype_id[i] =
            get_gluniform_location(prog.prog, "ltype[" + is + "]");
    }
    prog.shp_xform_id = get_gluniform_location(prog.prog, "shape_xform");
    prog.shp_normal_offset_id =
        get_gluniform_location(prog.prog, "shape_normal_offset");
    prog.highlight_id = get_gluniform_location(prog.prog, "highlight");
    prog.mtype_id = get_gluniform_location(prog.prog, "mat_type");
    prog.ke_id = get_gluniform_location(prog.prog, "mat_ke");
    prog.kd_id = get_gluniform_location(prog.prog, "mat_kd");
    prog.ks_id = get_gluniform_location(prog.prog, "mat_ks");
    prog.rs_id = get_gluniform_location(prog.prog, "mat_rs");
    prog.op_id = get_gluniform_location(prog.prog, "mat_op");
    prog.ke_txt_id = get_gluniform_location(prog.prog, "mat_txt_ke");
    prog.ke_txt_on_id = get_gluniform_location(prog.prog, "mat_txt_ke_on");
    prog.kd_txt_id = get_gluniform_location(prog.prog, "mat_txt_kd");
    prog.kd_txt_on_id = get_gluniform_location(prog.prog, "mat_txt_kd_on");
    prog.ks_txt_id = get_gluniform_location(prog.prog, "mat_txt_ks");
    prog.ks_txt_on_id = get_gluniform_location(prog.prog, "mat_txt_ks_on");
    prog.rs_txt_id = get_gluniform_location(prog.prog, "mat_txt_rs");
    prog.rs_txt_on_id = get_gluniform_location(prog.prog, "mat_txt_rs_on");
    prog.norm_txt_id = get_gluniform_location(prog.prog, "mat_txt_norm");
    prog.norm_txt_on_id = get_gluniform_location(prog.prog, "mat_txt_norm_on");
    prog.occ_txt_id = get_gluniform_location(prog.prog, "mat_txt_occ");
    prog.occ_txt_on_id = get_gluniform_location(prog.prog, "mat_txt_occ_on");
    prog.norm_scale_id = get_gluniform_location(prog.prog, "mat_norm_scale");
    prog.occ_scale_id = get_gluniform_location(prog.prog, "mat_occ_scale");
    prog.use_phong_id = get_gluniform_location(prog.prog, "mat_use_phong");
    prog.double_sided_id =
        get_gluniform_location(prog.prog, "mat_double_sided");
    prog.alpha_cutout_id =
        get_gluniform_location(prog.prog, "mat_alpha_cutout");
    prog.etype_id = get_gluniform_location(prog.prog, "elem_type");
    prog.efaceted_id = get_gluniform_location(prog.prog, "elem_faceted");
    assert(check_glerror());
    prog.pos_id = get_glattrib_location(prog.prog, "vert_pos");
    prog.norm_id = get_glattrib_location(prog.prog, "vert_norm");
    prog.texcoord_id = get_glattrib_location(prog.prog, "vert_texcoord");
    prog.color_id = get_glattrib_location(prog.prog, "vert_color");
    prog.tangsp_id = get_glattrib_location(prog.prog, "vert_tangsp");
    assert(check_glerror());
    return prog;
}

// Starts a frame by setting exposure/gamma values, camera transforms
// and projection. Sets also whether to use full shading or a quick
// eyelight preview.
void begin_glsurface_frame(const glsurface_program& prog, bool shade_eyelight,
    tonemap_type tonemap, float exposure, const mat4f& camera_xform,
    const mat4f& camera_xform_inv, const mat4f& camera_proj) {
    assert(check_glerror());
    bind_glprogram(prog.prog);
    set_gluniform(prog.prog, prog.eyelight_id, shade_eyelight);
    set_gluniform(prog.prog, prog.exposure_id, exposure);
    set_gluniform(prog.prog, prog.tonemap_id, (int)tonemap);
    set_gluniform(prog.prog, prog.cam_xform_id, camera_xform);
    set_gluniform(prog.prog, prog.cam_xform_inv_id, camera_xform_inv);
    set_gluniform(prog.prog, prog.cam_proj_id, camera_proj);
    assert(check_glerror());
    set_glsurface_elems(prog, glelem_type::triangle, false);
}

// Ends a frame.
void end_glsurface_frame(const glsurface_program& prog) {
    assert(check_glerror());
    unbind_glprogram(prog.prog);
    //    glBindVertexArray(0);
    //    glUseProgram(0);
    assert(check_glerror());
}

// Set num lights with position pos, color ke, type ltype. Also set the
// ambient illumination amb.
void set_glsurface_lights(
    const glsurface_program& prog, const vec3f& amb, const gllights& lights) {
    assert(check_glerror());
    set_gluniform(prog.prog, prog.lamb_id, amb);
    set_gluniform(prog.prog, prog.lnum_id, (int)lights.pos.size());
    for (auto i = 0; i < lights.pos.size(); i++) {
        set_gluniform(prog.prog, prog.lpos_id[i], lights.pos[i]);
        set_gluniform(prog.prog, prog.lke_id[i], lights.ke[i]);
        set_gluniform(prog.prog, prog.ltype_id[i], (int)lights.type[i]);
    }
    assert(check_glerror());
}

// Begins drawing a shape with transform xform.
void begin_glsurface_shape(
    const glsurface_program& prog, const mat4f& xform, float normal_offset) {
    assert(check_glerror());
    set_gluniform(prog.prog, prog.shp_xform_id, xform);
    set_gluniform(prog.prog, prog.shp_normal_offset_id, normal_offset);
    assert(check_glerror());
}

// End shade drawing.
void end_glsurface_shape(const glsurface_program& prog) {
    assert(check_glerror());
    for (int i = 0; i < 16; i++) unbind_glbuffer(i);
    assert(check_glerror());
}

// Sets normal offset.
void set_glsurface_normaloffset(
    const glsurface_program& prog, float normal_offset) {
    assert(check_glerror());
    set_gluniform(prog.prog, prog.shp_normal_offset_id, normal_offset);
    assert(check_glerror());
}

// Set the object as highlighted.
void set_glsurface_highlight(
    const glsurface_program& prog, const vec4f& highlight) {
    assert(check_glerror());
    set_gluniform(prog.prog, prog.highlight_id, highlight);
    assert(check_glerror());
}

// Set material values with emission ke, diffuse kd, specular ks and
// specular roughness rs, opacity op. Indicates textures ids with the
// correspoinding XXX_txt variables. Sets also normal and occlusion
// maps. Works for points/lines/triangles (diffuse for points,
// Kajiya-Kay for lines, GGX/Phong for triangles).
// Material type matches the scene material type.
void set_glsurface_material(const glsurface_program& prog, material_type type,
    const vec3f& ke, const vec3f& kd, const vec3f& ks, float rs, float op,
    const gltexture_info& ke_txt, const gltexture_info& kd_txt,
    const gltexture_info& ks_txt, const gltexture_info& rs_txt,
    const gltexture_info& norm_txt, const gltexture_info& occ_txt,
    bool use_phong, bool double_sided, bool alpha_cutout) {
    static auto mtypes = std::unordered_map<material_type, int>{
        {material_type::specular_roughness, 1},
        {material_type::metallic_roughness, 2},
        {material_type::specular_glossiness, 3}};

    assert(check_glerror());
    set_gluniform(prog.prog, prog.mtype_id, mtypes.at(type));
    set_gluniform(prog.prog, prog.ke_id, ke);
    set_gluniform(prog.prog, prog.kd_id, kd);
    set_gluniform(prog.prog, prog.ks_id, ks);
    set_gluniform(prog.prog, prog.rs_id, rs);
    set_gluniform(prog.prog, prog.op_id, op);
    set_gluniform_texture(
        prog.prog, prog.ke_txt_id, prog.ke_txt_on_id, ke_txt, 0);
    set_gluniform_texture(
        prog.prog, prog.kd_txt_id, prog.kd_txt_on_id, kd_txt, 1);
    set_gluniform_texture(
        prog.prog, prog.ks_txt_id, prog.ks_txt_on_id, ks_txt, 2);
    set_gluniform_texture(
        prog.prog, prog.rs_txt_id, prog.rs_txt_on_id, rs_txt, 3);
    set_gluniform_texture(
        prog.prog, prog.norm_txt_id, prog.norm_txt_on_id, norm_txt, 4);
    //    set_gluniform(prog.prog, prog.norm_scale_id, norm_txt.scale);
    //    set_gluniform(prog.prog, prog.occ_scale_id, occ_txt.scale);
    //    set_gluniform(prog.prog, prog.use_phong_id, use_phong);
    set_gluniform(prog.prog, prog.double_sided_id, double_sided);
    set_gluniform(prog.prog, prog.alpha_cutout_id, alpha_cutout);
    assert(check_glerror());
}

// Set constant material values with emission ke.
void set_glsurface_constmaterial(
    glsurface_program& prog, const vec3f& ke, float op) {
    assert(check_glerror());
    set_gluniform(prog.prog, prog.mtype_id, 0);
    set_gluniform(prog.prog, prog.ke_id, ke);
    set_gluniform(prog.prog, prog.op_id, op);
    assert(check_glerror());
}

// Set vertex data with buffers for position pos, normals norm, texture
// coordinates texcoord, per-vertex color color and tangent space
// tangsp.
void set_glsurface_vert(const glsurface_program& prog,
    const glvertex_buffer& pos, const glvertex_buffer& norm,
    const glvertex_buffer& texcoord, const glvertex_buffer& color,
    const glvertex_buffer& tangsp) {
    assert(check_glerror());
    set_glattribute(prog.prog, prog.pos_id, pos, zero3f);
    set_glattribute(prog.prog, prog.norm_id, norm, zero3f);
    set_glattribute(prog.prog, prog.texcoord_id, texcoord, zero2f);
    set_glattribute(prog.prog, prog.color_id, color, vec4f{1, 1, 1, 1});
    set_glattribute(prog.prog, prog.tangsp_id, tangsp, zero4f);
    assert(check_glerror());
}

// Set element properties.
void set_glsurface_elems(
    const glsurface_program& prog, glelem_type etype, bool faceted) {
    assert(check_glerror());
    set_gluniform(prog.prog, prog.etype_id, (int)etype);
    set_gluniform(prog.prog, prog.efaceted_id, (int)faceted);
    assert(check_glerror());
}

// Draw a shape
void draw_stdsurface_shape(const shape* shp, const material* mat,
    const mat4f& xform, bool highlighted, glsurface_program& prog,
    std::unordered_map<shape*, glshape>& shapes,
    std::unordered_map<texture*, gltexture>& textures,
    const glsurface_params& params) {
    static auto default_material = material();
    default_material.kd = {0.2f, 0.2f, 0.2f};

    begin_glsurface_shape(prog, xform);

    auto txt = [&textures](const texture_info& info) -> gltexture_info {
        auto ginfo = gltexture_info();
        if (!info.txt) return ginfo;
        ginfo.txt = textures.at(info.txt);
        return ginfo;
    };

    auto& gshp = shapes.at((shape*)shp);
    if (!mat) mat = &default_material;
    auto faceted = shp->norm.empty();

    set_glsurface_material(prog, mat->type, mat->ke, mat->kd, mat->ks, mat->rs,
        mat->op, txt(mat->ke_txt), txt(mat->kd_txt), txt(mat->ks_txt),
        txt(mat->rs_txt), txt(mat->norm_txt), txt(mat->occ_txt), false,
        mat->double_sided || params.double_sided, params.cutout);

    set_glsurface_vert(
        prog, gshp.pos, gshp.norm, gshp.texcoord, gshp.color, gshp.tangsp);

    if (!is_glbuffer_empty(gshp.points)) {
        set_glsurface_elems(prog, glelem_type::point, faceted);
        draw_glelems(gshp.points);
    }
    if (!is_glbuffer_empty(gshp.lines)) {
        set_glsurface_elems(prog, glelem_type::line, faceted);
        draw_glelems(gshp.lines);
    }
    if (!is_glbuffer_empty(gshp.triangles)) {
        set_glsurface_elems(prog, glelem_type::triangle, faceted);
        draw_glelems(gshp.triangles);
    }
    if (!is_glbuffer_empty(gshp.quads)) {
        set_glsurface_elems(prog, glelem_type::triangle, faceted);
        draw_glelems(gshp.quads);
    }
    if (!is_glbuffer_empty(gshp.beziers)) {
        set_glsurface_elems(prog, glelem_type::line, faceted);
        draw_glelems(gshp.beziers);
    }

    if ((params.edges && !params.wireframe) || highlighted) {
        enable_glculling(false);
        set_glsurface_constmaterial(prog,
            (highlighted) ? params.highlight_color : params.edge_color,
            (highlighted) ? 1 : mat->op);
        set_glsurface_normaloffset(prog, params.edge_offset);
        if (is_glbuffer_empty(gshp.edges)) draw_glelems(gshp.edges);
        enable_glculling(params.cull_backface);
    }

    end_glsurface_shape(prog);
}

// Display a scene
void draw_glsurface_scene(const scene* scn, const camera* cam,
    glsurface_program& prog, std::unordered_map<shape*, glshape>& shapes,
    std::unordered_map<texture*, gltexture>& textures, const gllights& lights,
    const vec2i& viewport_size, const void* highlighted,
    const glsurface_params& params) {
    // begin frame
    enable_gldepth_test(true);
    enable_glculling(params.cull_backface && !params.double_sided);
    enable_glwireframe(params.wireframe);
    set_glviewport(viewport_size);

    auto camera_xform = frame_to_mat(cam->frame);
    auto camera_view = frame_to_mat(inverse(cam->frame));
    auto camera_proj = perspective_mat(cam->yfov,
        (float)viewport_size.x / (float)viewport_size.y, cam->near, cam->far);

    begin_glsurface_frame(prog, params.eyelight, params.tonemapper,
        params.exposure, camera_xform, camera_view, camera_proj);

    if (!params.eyelight) {
        set_glsurface_lights(prog, params.ambient, lights);
    }

    for (auto ist : scn->instances) {
        draw_stdsurface_shape(ist->shp, ist->mat, frame_to_mat(ist->frame),
            ist == highlighted || ist->shp == highlighted ||
                ist->mat == highlighted,
            prog, shapes, textures, params);
    }

    end_glsurface_frame(prog);
    enable_glwireframe(false);
}

glwindow::~glwindow() {
    if (widget_enabled) {
        ImGui_ImplGlfwGL3_Shutdown();
        widget_enabled = false;
    }
    if (gwin) {
        glfwDestroyWindow(gwin);
        glfwTerminate();
        gwin = nullptr;
    }
}

// Support
void _glfw_error_cb(int error, const char* description) {
    log_error("GLFW error: {}\n", description);
}

// Support
void _glfw_text_cb(GLFWwindow* gwin, unsigned key) {
    auto win = (glwindow*)glfwGetWindowUserPointer(gwin);
    if (win->widget_enabled) { ImGui_ImplGlfw_CharCallback(win->gwin, key); }
    if (win->text_cb) win->text_cb(key);
}

// Support
void _glfw_key_cb(
    GLFWwindow* gwin, int key, int scancode, int action, int mods) {
    auto win = (glwindow*)glfwGetWindowUserPointer(gwin);
    if (win->widget_enabled) {
        ImGui_ImplGlfw_KeyCallback(win->gwin, key, scancode, action, mods);
    }
}

// Support
void _glfw_mouse_cb(GLFWwindow* gwin, int button, int action, int mods) {
    auto win = (glwindow*)glfwGetWindowUserPointer(gwin);
    if (win->widget_enabled) {
        ImGui_ImplGlfw_MouseButtonCallback(win->gwin, button, action, mods);
    }
    if (win->mouse_cb) win->mouse_cb(button, action == GLFW_PRESS, mods);
}

// Support
void _glfw_scroll_cb(GLFWwindow* gwin, double xoffset, double yoffset) {
    auto win = (glwindow*)glfwGetWindowUserPointer(gwin);
    if (win->widget_enabled) {
        ImGui_ImplGlfw_ScrollCallback(win->gwin, xoffset, yoffset);
    }
}

// Support
void _glfw_refresh_cb(GLFWwindow* gwin) {
    auto win = (glwindow*)glfwGetWindowUserPointer(gwin);
    if (win->refresh_cb) win->refresh_cb();
}

// Initialize glwindow
glwindow* make_glwindow(
    int width, int height, const std::string& title, bool opengl4) {
    auto win = new glwindow();

    // glwindow
    glfwSetErrorCallback(_glfw_error_cb);
    if (!glfwInit()) throw std::runtime_error("cannot open glwindow");

    // profile creation
    if (opengl4) {
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    } else {
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    }
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#if __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    win->gwin = glfwCreateWindow(width, height, title.c_str(), 0, 0);
    glfwMakeContextCurrent(win->gwin);
    glfwSetWindowUserPointer(win->gwin, win);
    glfwSwapInterval(1);  // Enable vsync

    glfwSetCharCallback(win->gwin, _glfw_text_cb);
    glfwSetKeyCallback(win->gwin, _glfw_key_cb);
    glfwSetMouseButtonCallback(win->gwin, _glfw_mouse_cb);
    glfwSetScrollCallback(win->gwin, _glfw_scroll_cb);

    glfwSetWindowRefreshCallback(win->gwin, _glfw_refresh_cb);

// init gl extensions
#ifndef __APPLE__
    if (glewInit() != GLEW_OK) return nullptr;
#endif
    return win;
}

// Set glwindow callbacks
void set_glwindow_callbacks(glwindow* win, text_glcallback text_cb,
    mouse_glcallback mouse_cb, refresh_glcallback refresh_cb) {
    win->text_cb = text_cb;
    win->mouse_cb = mouse_cb;
    win->refresh_cb = refresh_cb;
    if (win->text_cb) glfwSetCharCallback(win->gwin, _glfw_text_cb);
}

// Set glwindow title
void set_glwindow_title(glwindow* win, const std::string& title) {
    glfwSetWindowTitle(win->gwin, title.c_str());
}

// Wait events
void wait_glwindow_events(glwindow* win) { glfwWaitEvents(); }

// Poll events
void poll_glwindow_events(glwindow* win) { glfwPollEvents(); }

#ifdef __APPLE__

// Wait events
void wait_glwindow_events_timeout(glwindow* win, double timeout_sec) {
    glfwWaitEventsTimeout(timeout_sec);
}

// Wait events
void post_glwindow_event(glwindow* win) { glfwPostEmptyEvent(); }

#endif

// Swap buffers
void swap_glwindow_buffers(glwindow* win) { glfwSwapBuffers(win->gwin); }

// Should close
bool should_glwindow_close(glwindow* win) {
    return glfwWindowShouldClose(win->gwin);
}

// Mouse button
int get_glmouse_button(glwindow* win) {
    auto mouse1 =
        glfwGetMouseButton(win->gwin, GLFW_MOUSE_BUTTON_1) == GLFW_PRESS;
    auto mouse2 =
        glfwGetMouseButton(win->gwin, GLFW_MOUSE_BUTTON_2) == GLFW_PRESS;
    auto mouse3 =
        glfwGetMouseButton(win->gwin, GLFW_MOUSE_BUTTON_3) == GLFW_PRESS;
    if (mouse1) return 1;
    if (mouse2) return 2;
    if (mouse3) return 3;
#if 0
        if (action == GLFW_RELEASE) {
            vparams.mouse_button = 0;
        } else if (button == GLFW_MOUSE_BUTTON_1 && !mods) {
            vparams.mouse_button = 1;
        } else if (button == GLFW_MOUSE_BUTTON_1 && (mods & GLFW_MOD_CONTROL)) {
            vparams.mouse_button = 2;
        } else if (button == GLFW_MOUSE_BUTTON_1 && (mods & GLFW_MOD_SHIFT)) {
            vparams.mouse_button = 3;
        } else if (button == GLFW_MOUSE_BUTTON_2) {
            vparams.mouse_button = 2;
        } else {
            vparams.mouse_button = 0;
        }
#endif
    return 0;
}

// Mouse position
vec2i get_glmouse_pos(glwindow* win) {
    double x, y;
    glfwGetCursorPos(win->gwin, &x, &y);
    return {(int)x, (int)y};
}

// Mouse position
vec2f get_glmouse_posf(glwindow* win) {
    double x, y;
    glfwGetCursorPos(win->gwin, &x, &y);
    return {(float)x, (float)y};
}

// Window size
vec2i get_glwindow_size(glwindow* win) {
    auto ret = vec2i{0, 0};
    glfwGetWindowSize(win->gwin, &ret.x, &ret.y);
    return ret;
}

// Check if a key is pressed (not all keys are supported)
bool get_glkey(glwindow* win, int key) {
    key = std::toupper(key);
    return glfwGetKey(win->gwin, key) == GLFW_PRESS;
}

// Check if the alt key is down
bool get_glalt_key(glwindow* win) {
    return glfwGetKey(win->gwin, GLFW_KEY_LEFT_ALT) == GLFW_PRESS ||
           glfwGetKey(win->gwin, GLFW_KEY_RIGHT_ALT) == GLFW_PRESS;
}

// Check if the alt key is down
bool get_glctrl_key(glwindow* win) {
    return glfwGetKey(win->gwin, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS ||
           glfwGetKey(win->gwin, GLFW_KEY_RIGHT_CONTROL) == GLFW_PRESS;
}

// Check if the alt key is down
bool get_glshift_key(glwindow* win) {
    return glfwGetKey(win->gwin, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
           glfwGetKey(win->gwin, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS;
}

// Framebuffer size
vec2i get_glframebuffer_size(glwindow* win) {
    auto ret = vec2i{0, 0};
    glfwGetFramebufferSize(win->gwin, &ret.x, &ret.y);
    return ret;
}

// Read pixels
image4b take_glscreenshot4b(glwindow* win, bool flipy, bool back) {
    auto wh = get_glframebuffer_size(win);
    auto img = make_image4b(wh.x, wh.y);
    glReadBuffer((back) ? GL_BACK : GL_FRONT);
    glReadPixels(0, 0, img.width, img.height, GL_RGBA, GL_UNSIGNED_BYTE,
        img.pixels.data());
    if (flipy) {
        for (int j = 0; j < img.height / 2; j++) {
            for (auto i = 0; i < img.width; i++) {
                std::swap(img.at(i, j), img.at(i, img.height - 1 - j));
            }
        }
    }
    return img;
}

// Handles camera navigation
bool handle_glcamera_navigation(
    glwindow* win, camera* cam, bool navigation_fps) {
    static auto mouse_last = zero2f;
    auto mouse_pos = get_glmouse_posf(win);
    auto mouse_button = get_glmouse_button(win);
    auto alt_down = get_glalt_key(win);

    // updated
    auto updated = false;

    // handle mouse and keyboard for navigation
    if (mouse_button && alt_down && !get_imgui_active(win)) {
        if (mouse_button == 1 && get_glshift_key(win)) mouse_button = 3;
        if (navigation_fps) {
            auto dolly = 0.0f;
            auto pan = zero2f;
            auto rotate = zero2f;
            switch (mouse_button) {
                case 1: rotate = (mouse_pos - mouse_last) / 100.0f; break;
                case 2: dolly = (mouse_pos.x - mouse_last.x) / 100.0f; break;
                case 3: pan = (mouse_pos - mouse_last) / 100.0f; break;
                default: break;
            }
            camera_fps(cam->frame, {0, 0, 0}, rotate);
            updated = true;
        } else {
            auto dolly = 0.0f;
            auto pan = zero2f;
            auto rotate = zero2f;
            switch (mouse_button) {
                case 1: rotate = (mouse_pos - mouse_last) / 100.0f; break;
                case 2: dolly = (mouse_pos.x - mouse_last.x) / 100.0f; break;
                case 3: pan = (mouse_pos - mouse_last) / 100.0f; break;
                default: break;
            }
            camera_turntable(cam->frame, cam->focus, rotate, dolly, pan);
            updated = true;
        }
    }

    // handle keytboard for navigation
    if (!get_imgui_active(win) && navigation_fps) {
        auto transl = zero3f;
        if (get_glkey(win, 'a')) transl.x -= 1;
        if (get_glkey(win, 'd')) transl.x += 1;
        if (get_glkey(win, 's')) transl.z += 1;
        if (get_glkey(win, 'w')) transl.z -= 1;
        if (get_glkey(win, 'e')) transl.y += 1;
        if (get_glkey(win, 'q')) transl.y -= 1;
        if (transl != zero3f) {
            camera_fps(cam->frame, transl, {0, 0});
            updated = true;
        }
    }

    // record mouse position
    mouse_last = mouse_pos;

    // done
    return updated;
}

// Handle scene selection
bool handle_glscene_selection(glwindow* win, const scene* scn,
    const camera* cam, const bvh_tree* bvh, int res,
    const glimage_params& params, scene_selection& sel) {
    auto mouse_pos = get_glmouse_posf(win);
    auto mouse_button = get_glmouse_button(win);

    if (!(mouse_button == 1 && !get_imgui_active(win))) return false;
    auto ij = get_draw_image_coords(mouse_pos, params);
    if (ij.x < 0 || ij.x >= (int)round(res * cam->aspect) || ij.y < 0 ||
        ij.y >= res)
        return false;
    auto ray = eval_camera_ray(cam, ij, res, {0.5f, 0.5f}, zero2f);
    auto isec = intersect_bvh(bvh, ray, false);
    if (!isec) return false;
    if (scn->nodes.empty()) {
        sel = scn->shapes[isec.iid];
    } else {
        sel = scn->nodes[isec.iid];
    }
    return true;
}

// Initialize widgets
void init_imgui(glwindow* win, bool light_style, bool alt_font) {
    ImGui::CreateContext();
    ImGui_ImplGlfwGL3_Init(win->gwin, false);
    ImGuiIO& io = ImGui::GetIO();
    // io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;  // Enable Keyboard
    // Controls io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;   // Enable
    // Gamepad Controls
    ImGui::GetStyle().WindowRounding = 0;
    io.IniFilename = nullptr;
    auto size = get_glwindow_size(win);
    ImGui::SetNextWindowPos({(float)size.x - 320, 0});
    ImGui::SetNextWindowSize({(float)320, (float)size.y});
    if (light_style) ImGui::StyleColorsLight();
    if (alt_font) {
        io.Fonts->AddFontFromMemoryCompressedTTF(
            imgui_extrafont_compressed_data(),
            imgui_extrafont_compressed_size(), 16);
    } else {
        io.Fonts->AddFontDefault();
    }
    win->widget_enabled = true;
}

// Begin draw widget
bool begin_imgui_frame(glwindow* win, const std::string& title) {
    static bool first_time = true;
    ImGui_ImplGlfwGL3_NewFrame();
    if (first_time) {
        auto size = get_glwindow_size(win);
        ImGui::SetNextWindowPos({(float)size.x - 320, 0});
        ImGui::SetNextWindowSize({(float)320, (float)size.y});
        ImGui::SetNextWindowCollapsed(true);
        first_time = false;
    }
    ImGui::Begin(title.c_str(), nullptr);
    // ImGui::ShowTestWindow();
    // ImGui::ShowStyleEditor();
    return true;
}

// End draw widget
void end_imgui_frame(glwindow* win) {
    ImGui::End();
    ImGui::Render();
    ImGui_ImplGlfwGL3_RenderDrawData(ImGui::GetDrawData());
}

// Whether widget are active
bool get_imgui_active(glwindow* win) {
    if (!win->widget_enabled) return false;
    auto io = &ImGui::GetIO();
    return io->WantTextInput || io->WantCaptureMouse || io->WantCaptureKeyboard;
}

// Horizontal separator
void draw_imgui_separator(glwindow* win) { ImGui::Separator(); }

// Indent widget
void begin_imgui_indent(glwindow* win) { ImGui::Indent(); }

// Indent widget
void end_imgui_indent(glwindow* win) { ImGui::Unindent(); }

// Continue line with next widget
void continue_imgui_line(glwindow* win) { ImGui::SameLine(); }

// Label widget
void draw_imgui_label(
    glwindow* win, const std::string& lbl, const std::string& msg) {
    ImGui::LabelText(lbl.c_str(), "%s", msg.c_str());
}

// Value widget
bool draw_imgui_text(glwindow* win, const std::string& lbl, std::string& str) {
    char buf[4096];
    if (str.length() >= 4096) throw std::runtime_error("bad memory");
    memcpy(buf, str.c_str(), str.size());
    buf[str.size()] = 0;
    auto ret = ImGui::InputText(lbl.c_str(), buf, 4096);
    str = buf;
    return ret;
}

// Value widget
bool draw_imgui_multiline_text(
    glwindow* win, const std::string& lbl, std::string& str) {
    char sbuf[8192];
    std::vector<char> dbuf;
    char* buf = nullptr;
    int buf_size = 0;
    if (str.size() > sizeof(sbuf) / 2) {
        dbuf.resize(str.size() * 2);
        buf = dbuf.data();
        buf_size = dbuf.size();
    } else {
        buf = sbuf;
        buf_size = sizeof(sbuf);
    }
    memcpy(buf, str.c_str(), str.size());
    buf[str.size()] = 0;
    auto ret = ImGui::InputTextMultiline(lbl.c_str(), buf, buf_size);
    str = buf;
    return ret;
}

static float draw_drag_scale = 1 / 100.0f;

// Drag widget.
bool draw_imgui_dragbox(
    glwindow* win, const std::string& lbl, int& val, int min, int max) {
    auto speed = (min == max) ? 1 : draw_drag_scale * (max - min);
    return ImGui::DragInt(lbl.c_str(), &val, speed, min, max);
}
// Drag widget.
bool draw_imgui_dragbox(
    glwindow* win, const std::string& lbl, vec2i& val, int min, int max) {
    auto speed = (min == max) ? 1 : draw_drag_scale * (max - min);
    return ImGui::DragInt2(lbl.c_str(), &val.x, speed, min, max);
}
// Drag widget.
bool draw_imgui_dragbox(
    glwindow* win, const std::string& lbl, vec3i& val, int min, int max) {
    auto speed = (min == max) ? 1 : draw_drag_scale * (max - min);
    return ImGui::DragInt3(lbl.c_str(), &val.x, speed, min, max);
}
// Drag widget.
bool draw_imgui_dragbox(
    glwindow* win, const std::string& lbl, vec4i& val, int min, int max) {
    auto speed = (min == max) ? 1 : draw_drag_scale * (max - min);
    return ImGui::DragInt4(lbl.c_str(), &val.x, speed, min, max);
}
// Drag widget.
bool draw_imgui_dragbox(
    glwindow* win, const std::string& lbl, float& val, float min, float max) {
    auto speed = (min == max) ? 1 : draw_drag_scale * (max - min);
    return ImGui::DragFloat(lbl.c_str(), &val, speed, min, max);
}
// Drag widget.
bool draw_imgui_dragbox(
    glwindow* win, const std::string& lbl, vec2f& val, float min, float max) {
    auto speed = (min == max) ? 1 : draw_drag_scale * (max - min);
    return ImGui::DragFloat2(lbl.c_str(), &val.x, speed, min, max);
}
// Drag widget.
bool draw_imgui_dragbox(
    glwindow* win, const std::string& lbl, vec3f& val, float min, float max) {
    auto speed = (min == max) ? 1 : draw_drag_scale * (max - min);
    return ImGui::DragFloat3(lbl.c_str(), &val.x, speed, min, max);
}
// Drag widget.
bool draw_imgui_dragbox(
    glwindow* win, const std::string& lbl, vec4f& val, float min, float max) {
    auto speed = (min == max) ? 1 : draw_drag_scale * (max - min);
    return ImGui::DragFloat4(lbl.c_str(), &val.x, speed, min, max);
}
// Drag widget.
bool draw_imgui_dragbox(
    glwindow* win, const std::string& lbl, mat4f& val, float min, float max) {
    auto modx = draw_imgui_dragbox(win, lbl + ".x", val.x, min, max);
    auto mody = draw_imgui_dragbox(win, lbl + ".y", val.y, min, max);
    auto modz = draw_imgui_dragbox(win, lbl + ".z", val.z, min, max);
    auto modw = draw_imgui_dragbox(win, lbl + ".w", val.w, min, max);
    return modx || mody || modz || modw;
}
// Drag widget.
bool draw_imgui_dragbox(
    glwindow* win, const std::string& lbl, frame3f& val, float min, float max) {
    auto modx = draw_imgui_dragbox(win, lbl + ".x", val.x, -1, 1);
    auto mody = draw_imgui_dragbox(win, lbl + ".y", val.y, -1, 1);
    auto modz = draw_imgui_dragbox(win, lbl + ".z", val.z, -1, 1);
    auto modo = draw_imgui_dragbox(win, lbl + ".o", val.o, min, max);
    return modx || mody || modz || modo;
}

// Color widget
bool draw_imgui_colorbox(glwindow* win, const std::string& lbl, vec3f& val) {
    auto mod = ImGui::ColorEdit3(
        lbl.c_str(), (float*)&val.x, ImGuiColorEditFlags_Float);
    // fix for bug in ImGui
    if (mod) {
        if (val.x < 0.0001f) val.x = 0;
        if (val.y < 0.0001f) val.y = 0;
        if (val.z < 0.0001f) val.z = 0;
    }
    return mod;
}
// Color widget
bool draw_imgui_colorbox(glwindow* win, const std::string& lbl, vec4f& val) {
    auto mod = ImGui::ColorEdit4(
        lbl.c_str(), (float*)&val.x, ImGuiColorEditFlags_Float);
    // fix for bug in ImGui
    if (mod) {
        if (val.x < 0.0001f) val.x = 0;
        if (val.y < 0.0001f) val.y = 0;
        if (val.z < 0.0001f) val.z = 0;
        if (val.w < 0.0001f) val.w = 0;
    }
    return mod;
}
// Color widget
bool draw_imgui_colorbox(glwindow* win, const std::string& lbl, vec4b& val) {
    auto valf = ImGui::ColorConvertU32ToFloat4(*(uint32_t*)&val);
    if (ImGui::ColorEdit4(lbl.c_str(), &valf.x)) {
        auto valb = ImGui::ColorConvertFloat4ToU32(valf);
        *(uint32_t*)&val = valb;
        return true;
    }
    return false;
}
bool draw_hdr_color_widget(
    glwindow* win, const std::string& lbl, vec3f& val, float max) {
    auto vall = ygl::max(val);
    auto valc = val;
    if (vall > 1) {
        valc /= vall;
    } else {
        vall = 1;
    }
    auto mod1 = draw_imgui_dragbox(win, lbl + "(m)", vall, 0, max);
    auto mod2 = draw_imgui_colorbox(win, lbl, valc);
    if (mod1 || mod2) {
        val = valc * vall;
        return true;
    } else {
        return false;
    }
}

// Bool widget
bool draw_imgui_checkbox(glwindow* win, const std::string& lbl, bool& val) {
    return ImGui::Checkbox(lbl.c_str(), &val);
}

// Combo widget
bool begin_imgui_combobox(
    glwindow* win, const std::string& lbl, const std::string& label) {
    return ImGui::BeginCombo(lbl.c_str(), label.c_str());
}

// Combo widget
bool draw_imgui_item(
    glwindow* win, const std::string& label, int idx, bool selected) {
    ImGui::PushID((void*)(intptr_t)idx);
    auto clicked = ImGui::Selectable(label.c_str(), selected);
    if (selected) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
    return clicked;
}

// Combo widget
void end_imgui_combobox(glwindow* win) { ImGui::EndCombo(); }

// Button widget
bool draw_imgui_button(glwindow* win, const std::string& lbl) {
    return ImGui::Button(lbl.c_str());
}

// Collapsible header
bool draw_imgui_header(glwindow* win, const std::string& lbl) {
    return ImGui::CollapsingHeader(lbl.c_str());
}

// Start tree node
bool begin_imgui_tree(glwindow* win, const std::string& lbl) {
    return ImGui::TreeNode(lbl.c_str());
}

// Collapsible header
void end_imgui_tree(glwindow* win) { ImGui::TreePop(); }

// Start selectable tree node
bool begin_imgui_tree(
    glwindow* win, const std::string& lbl, void*& selection, void* content) {
    ImGuiTreeNodeFlags node_flags =
        ImGuiTreeNodeFlags_OpenOnArrow | ImGuiTreeNodeFlags_OpenOnDoubleClick;
    if (selection == content) node_flags |= ImGuiTreeNodeFlags_Selected;
    auto open = ImGui::TreeNodeEx(content, node_flags, "%s", lbl.c_str());
    if (ImGui::IsItemClicked()) selection = content;
    return open;
}

// End selectable tree node
void end_imgui_tree(glwindow* win, void* content) { ImGui::TreePop(); }

// Selectable tree leaf node
void draw_imgui_tree_leaf(
    glwindow* win, const std::string& lbl, void*& selection, void* content) {
    ImGuiTreeNodeFlags node_flags =
        ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen;
    if (selection == content) node_flags |= ImGuiTreeNodeFlags_Selected;
    ImGui::TreeNodeEx(content, node_flags, "%s", lbl.c_str());
    if (ImGui::IsItemClicked()) selection = content;
}

// Selectable tree leaf node
void draw_imgui_tree_leaf(glwindow* win, const std::string& lbl,
    void*& selection, void* content, const vec4f& col) {
    ImGui::PushStyleColor(ImGuiCol_Text, {col.x, col.y, col.z, col.w});
    draw_imgui_tree_leaf(win, lbl, selection, content);
    ImGui::PopStyleColor();
}

// Image widget
void draw_imgui_imagebox(
    glwindow* win, int tid, const vec2i& size, const vec2i& imsize) {
    auto w = ImGui::GetContentRegionAvailWidth();
    auto s = vec2f{(float)size.x, (float)size.y};
    auto a = (float)imsize.x / (float)imsize.y;
    if (!s.x && !s.y) {
        s.x = w;
        s.y = w / a;
    } else if (s.x && !s.y) {
        s.y = s.x / a;
    } else if (!s.x && s.y) {
        s.x = s.y * a;
    } else {
        auto as = s.x / s.y;
        if (as / a > 1) {
            s.x = s.y * a;
        } else {
            s.y = s.x / a;
        }
    }
    if (s.x > w) {
        s.x = w;
        s.y = w / a;
    }
    ImGui::Image((void*)(size_t)tid, {s.x, s.y});
}

// Image widget
void draw_imgui_imagebox(
    glwindow* win, const gltexture& txt, const vec2i& size) {
    draw_imgui_imagebox(
        win, get_gltexture_id(txt), size, {txt.width, txt.height});
}

// Scroll region
void begin_imgui_scrollarea(
    glwindow* win, const std::string& lbl, int height, bool border) {
    ImGui::BeginChild(lbl.c_str(), ImVec2(0, height), border);
}
// Scroll region
void end_imgui_scrollarea(glwindow* win) { ImGui::EndChild(); }
// Scroll region
void move_imgui_scrollarea(glwindow* win) { ImGui::SetScrollHere(); }

// Group ids
void push_imgui_groupid(glwindow* win, int gid) { ImGui::PushID(gid); }
// Group ids
void push_imgui_groupid(glwindow* win, void* gid) { ImGui::PushID(gid); }
// Group ids
void push_imgui_groupid(glwindow* win, const void* gid) { ImGui::PushID(gid); }
// Group ids
void push_imgui_groupid(glwindow* win, const char* gid) { ImGui::PushID(gid); }
// Group ids
void pop_imgui_groupid(glwindow* win) { ImGui::PopID(); }

// Widget style
void push_imgui_style(glwindow* win, const vec4f& col) {
    ImGui::PushStyleColor(ImGuiCol_Text, {col.x, col.y, col.z, col.w});
}

// Widget style
void pop_imgui_style(glwindow* win) { ImGui::PopStyleColor(); }

// Image inspection widgets.
void draw_imgui_image_inspector(glwindow* win, const std::string& lbl,
    const image4f& hdr, const image4b& ldr, const vec2f& mouse_pos,
    const glimage_params& params) {
    auto ij = get_draw_image_coords(mouse_pos, params);
    auto v4f = zero4f;
    auto v4b = zero4b;
    if (!hdr.pixels.empty() && contains(hdr, ij.x, ij.y)) {
        v4f = hdr.at(ij.x, ij.y);
        v4b = linear_to_srgb(hdr.at(ij.x, ij.y));
    } else if (!ldr.pixels.empty() && contains(ldr, ij.x, ij.y)) {
        v4f = srgb_to_linear(ldr.at(ij.x, ij.y));
        v4b = ldr.at(ij.x, ij.y);
    }
    char buf[1024];
    sprintf(buf, "%5d %5d", ij.x, ij.y);
    draw_imgui_label(win, lbl + "mouse pos", buf);
    sprintf(buf, "%4.4g %4.4g %4.4g %4.4g", v4f.x, v4f.y, v4f.z, v4f.w);
    draw_imgui_label(win, lbl + "hdr val", buf);
    sprintf(
        buf, "%3d %3d %3d %3d", (int)v4b.x, (int)v4b.y, (int)v4b.z, (int)v4b.w);
    draw_imgui_label(win, lbl + "ldr val", buf);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR SIMPLE SCENE UI
// -----------------------------------------------------------------------------
namespace ygl {

// Draws widgets for params.
bool draw_imgui_stdimage_inspector(
    glwindow* win, const std::string& lbl, glimage_params& params) {
    auto edited = 0;
    edited += draw_imgui_dragbox(win, "offset", params.offset, -4096, 4096);
    edited += draw_imgui_dragbox(win, "zoom", params.zoom, 0.01, 10);
    edited += draw_imgui_colorbox(win, "background", params.background);
    edited += draw_imgui_combobox(
        win, "tonemap", params.tonemapper, tonemap_type_names());
    edited += draw_imgui_dragbox(win, "exposure", params.exposure, -5, 5);
    return edited;
}
bool draw_imgui_stdsurface_inspector(
    glwindow* win, const std::string& lbl, glsurface_params& params) {
    auto edited = 0;
    edited +=
        draw_imgui_dragbox(win, "resolution", params.resolution, 256, 4096);
    edited += draw_imgui_checkbox(win, "wireframe", params.wireframe);
    edited += draw_imgui_checkbox(win, "edges", params.edges);
    edited +=
        draw_imgui_dragbox(win, "edge_offset", params.edge_offset, 0, 0.1);
    edited += draw_imgui_checkbox(win, "cutout", params.cutout);
    edited += draw_imgui_checkbox(win, "eyelight", params.eyelight);
    edited += draw_imgui_combobox(
        win, "tonemap", params.tonemapper, tonemap_type_names());
    edited += draw_imgui_dragbox(win, "exposure", params.exposure);
    edited += draw_imgui_colorbox(win, "background", params.background);
    edited += draw_imgui_colorbox(win, "ambient", params.ambient);
    edited +=
        draw_imgui_colorbox(win, "highlight_color", params.highlight_color);
    edited += draw_imgui_colorbox(win, "edge_color", params.edge_color);
    edited += draw_imgui_checkbox(win, "double_sided", params.double_sided);
    edited += draw_imgui_checkbox(win, "cull_backface", params.cull_backface);
    return edited;
}
bool draw_imgui_trace_inspector(
    glwindow* win, const std::string& lbl, trace_params& params) {
    auto edited = 0;
    edited +=
        draw_imgui_dragbox(win, "resolution", params.resolution, 256, 4096);
    edited += draw_imgui_dragbox(win, "nsamples", params.nsamples, 16, 4096);
    edited +=
        draw_imgui_combobox(win, "tracer", params.tracer, trace_type_names());
    edited += draw_imgui_checkbox(win, "notransmission", params.notransmission);
    edited += draw_imgui_checkbox(win, "double_sided", params.double_sided);
    edited += draw_imgui_colorbox(win, "ambient", params.ambient);
    edited +=
        draw_imgui_checkbox(win, "envmap_invisible", params.envmap_invisible);
    edited += draw_imgui_dragbox(win, "min_depth", params.min_depth, 1, 10);
    edited += draw_imgui_dragbox(win, "max_depth", params.max_depth, 1, 10);
    edited +=
        draw_imgui_dragbox(win, "pixel_clamp", params.pixel_clamp, 10, 1000);
    edited +=
        draw_imgui_dragbox(win, "ray_eps", params.ray_eps, 0.0001f, 0.001f);
    edited += draw_imgui_checkbox(win, "parallel", params.parallel);
    edited += draw_imgui_dragbox(win, "seed", (int&)params.seed, 0, 1000);
    edited +=
        draw_imgui_dragbox(win, "preview", params.preview_resolution, 64, 1080);
    return edited;
}

// Implementation of camera selection
bool draw_imgui_camera_inspector(
    glwindow* win, const std::string& lbl, camera* cam) {
    if (!cam) return false;
    auto edited = 0;
    auto from = cam->frame.o;
    auto to = cam->frame.o - cam->frame.z * cam->focus;
    edited += draw_imgui_text(win, lbl + " name", cam->name);
    edited += draw_imgui_dragbox(win, lbl + " from", from, -10, 10);
    edited += draw_imgui_dragbox(win, lbl + " to", to, -10, 10);
    edited += draw_imgui_dragbox(win, lbl + " yfov", cam->yfov, 0.1, 10);
    edited += draw_imgui_dragbox(win, lbl + " aspect", cam->aspect, 1, 3);
    edited += draw_imgui_dragbox(win, lbl + " aperture", cam->aperture, 0, 1);
    edited += draw_imgui_dragbox(win, lbl + " focus", cam->focus, 0.1, 100);
    edited += draw_imgui_dragbox(win, lbl + " near", cam->near, 0.01f, 1);
    edited += draw_imgui_dragbox(win, lbl + " far", cam->far, 1, 1000);
    if (edited) cam->frame = lookat_frame(from, to, {0, 1, 0});
    return edited;
}

static const std::unordered_map<std::string, vec4f>
    draw_visitor_highlight_colors = {{"red", {1, 0.5f, 0.5f, 1}},
        {"green", {0.5f, 1, 0.5f, 1}}, {"blue", {0.5f, 0.5f, 1, 1}}};

vec4f get_highlight_color(
    const std::unordered_map<std::string, std::string>& highlights,
    const std::string& name) {
    if (!highlights.empty() && highlights.find(name) != highlights.end()) {
        return draw_visitor_highlight_colors.at(highlights.at(name));
    }
    return zero4f;
}

template <typename T>
void draw_scene_tree_widgets_rec(glwindow* win, const std::string& lbl_, T* val,
    scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {}

template <typename T>
void draw_imgui_scene_tree(glwindow* win, const std::string& lbl_, T* val,
    scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    if (!val) return;
    auto lbl = val->name;
    if (!lbl_.empty()) lbl = lbl_ + ": " + val->name;
    auto selection = sel.get_raw();
    auto color = get_highlight_color(highlights, val->name);
    if (color != zero4f) push_imgui_style(win, color);
    auto open = begin_imgui_tree(win, lbl, selection, val);
    if (color != zero4f) pop_imgui_style(win);
    if (selection == val) sel = val;
    if (open) {
        draw_scene_tree_widgets_rec(win, lbl_, val, sel, highlights);
        end_imgui_tree(win);
    }
}

template <>
void draw_scene_tree_widgets_rec<instance>(glwindow* win,
    const std::string& lbl_, instance* val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    draw_imgui_scene_tree(win, "shp", val->shp, sel, highlights);
    draw_imgui_scene_tree(win, "mat", val->mat, sel, highlights);
}

template <>
void draw_scene_tree_widgets_rec<material>(glwindow* win,
    const std::string& lbl_, material* val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    draw_imgui_scene_tree(win, "ke", val->ke_txt.txt, sel, highlights);
    draw_imgui_scene_tree(win, "kd", val->kd_txt.txt, sel, highlights);
    draw_imgui_scene_tree(win, "ks", val->ks_txt.txt, sel, highlights);
    draw_imgui_scene_tree(win, "kr", val->kr_txt.txt, sel, highlights);
    draw_imgui_scene_tree(win, "rs", val->rs_txt.txt, sel, highlights);
    draw_imgui_scene_tree(win, "op", val->op_txt.txt, sel, highlights);
    draw_imgui_scene_tree(win, "bump", val->bump_txt.txt, sel, highlights);
    draw_imgui_scene_tree(win, "disp", val->disp_txt.txt, sel, highlights);
    draw_imgui_scene_tree(win, "norm", val->norm_txt.txt, sel, highlights);
    draw_imgui_scene_tree(win, "occ", val->occ_txt.txt, sel, highlights);
}
template <>
void draw_scene_tree_widgets_rec<environment>(glwindow* win,
    const std::string& lbl_, environment* val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    draw_imgui_scene_tree(win, "ke", val->ke_txt.txt, sel, highlights);
}
template <>
void draw_scene_tree_widgets_rec<node>(glwindow* win, const std::string& lbl_,
    node* val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    draw_imgui_scene_tree(win, "ist", val->ist, sel, highlights);
    draw_imgui_scene_tree(win, "cam", val->cam, sel, highlights);
    draw_imgui_scene_tree(win, "env", val->env, sel, highlights);
    draw_imgui_scene_tree(win, "par", val->parent, sel, highlights);
    auto cid = 0;
    for (auto ch : val->children_) {
        draw_imgui_scene_tree(
            win, "ch" + std::to_string(cid++), ch, sel, highlights);
    }
}
template <>
void draw_scene_tree_widgets_rec<animation>(glwindow* win,
    const std::string& lbl_, animation* val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    auto tid = 0;
    for (auto tg : val->targets) {
        draw_imgui_scene_tree(
            win, "tg" + std::to_string(tid++), tg, sel, highlights);
    }
}

void draw_imgui_scene_tree(glwindow* win, const std::string& lbl,
    texture_info& val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    draw_imgui_scene_tree(win, lbl, val.txt, sel, highlights);
}

void draw_imgui_scene_tree(glwindow* win, scene* scn, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    if (!scn->cameras.empty() && begin_imgui_tree(win, "cameras")) {
        for (auto v : scn->cameras)
            draw_imgui_scene_tree(win, "", v, sel, highlights);
        end_imgui_tree(win);
    }
    if (!scn->shapes.empty() && begin_imgui_tree(win, "shapes")) {
        for (auto v : scn->shapes)
            draw_imgui_scene_tree(win, "", v, sel, highlights);
        end_imgui_tree(win);
    }
    if (!scn->instances.empty() && begin_imgui_tree(win, "instances")) {
        for (auto v : scn->instances)
            draw_imgui_scene_tree(win, "", v, sel, highlights);
        end_imgui_tree(win);
    }
    if (!scn->materials.empty() && begin_imgui_tree(win, "materials")) {
        for (auto v : scn->materials)
            draw_imgui_scene_tree(win, "", v, sel, highlights);
        end_imgui_tree(win);
    }
    if (!scn->textures.empty() && begin_imgui_tree(win, "textures")) {
        for (auto v : scn->textures)
            draw_imgui_scene_tree(win, "", v, sel, highlights);
        end_imgui_tree(win);
    }
    if (!scn->environments.empty() && begin_imgui_tree(win, "environments")) {
        for (auto v : scn->environments)
            draw_imgui_scene_tree(win, "", v, sel, highlights);
        end_imgui_tree(win);
    }
    if (!scn->nodes.empty() && begin_imgui_tree(win, "nodes")) {
        for (auto v : scn->nodes)
            draw_imgui_scene_tree(win, "", v, sel, highlights);
        end_imgui_tree(win);
    }
    if (!scn->animations.empty() && begin_imgui_tree(win, "animations")) {
        for (auto v : scn->animations)
            draw_imgui_scene_tree(win, "", v, sel, highlights);
        end_imgui_tree(win);
    }
}

template <typename T>
void draw_imgui_label(glwindow* win, const std::string& lbl,
    const std::vector<T>& val, bool skip_if_empty = true) {
    if (skip_if_empty && val.empty()) return;
    draw_imgui_label(win, lbl, std::to_string(val.size()));
}
void draw_imgui_label(glwindow* win, const std::string& lbl, const image4b& val,
    bool skip_if_empty = true) {
    if (skip_if_empty && val.pixels.empty()) return;
    draw_imgui_label(win, lbl,
        std::to_string(val.width) + " x " + std::to_string(val.height));
}
void draw_imgui_label(glwindow* win, const std::string& lbl, const image4f& val,
    bool skip_if_empty = true) {
    if (skip_if_empty && val.pixels.empty()) return;
    draw_imgui_label(win, lbl,
        std::to_string(val.width) + " x " + std::to_string(val.height));
}

/// Visit struct elements.
bool draw_imgui_scene_inspector(glwindow* win, camera* val, scene* scn) {
    auto edited = 0;
    edited += draw_imgui_text(win, "name", val->name);
    edited += draw_imgui_dragbox(win, "frame", val->frame);
    edited += draw_imgui_checkbox(win, "ortho", val->ortho);
    edited += draw_imgui_dragbox(win, "yfov", val->yfov, 0.1f, 10);
    edited += draw_imgui_dragbox(win, "aspect", val->aspect, 1, 3);
    edited += draw_imgui_dragbox(win, "focus", val->focus, 0.01f, 1000);
    edited += draw_imgui_dragbox(win, "aperture", val->aperture, 0, 5);
    edited += draw_imgui_dragbox(win, "near", val->near, 0.01f, 10);
    edited += draw_imgui_dragbox(win, "far", val->far, 10, 10000);
    return edited;
}

/// Visit struct elements.
bool draw_imgui_scene_inspector(glwindow* win, texture* val, scene* scn) {
    auto edited = 0;
    edited += draw_imgui_text(win, "name", val->name);
    edited += draw_imgui_text(win, "path", val->path);
    draw_imgui_label(win, "ldr", val->ldr);
    draw_imgui_label(win, "hdr", val->hdr);
    return edited;
}

bool draw_imgui_scene_inspector(
    glwindow* win, const std::string& lbl, texture_info& val, scene* scn) {
    auto edited = 0;
    edited += draw_imgui_combobox(win, lbl + " txt", val.txt, scn->textures);
    edited += draw_imgui_checkbox(win, lbl + " wrap_s", val.wrap_s);
    edited += draw_imgui_checkbox(win, lbl + " wrap_t", val.wrap_t);
    edited += draw_imgui_checkbox(win, lbl + " linear", val.linear);
    edited += draw_imgui_checkbox(win, lbl + " mipmap", val.mipmap);
    edited += draw_imgui_dragbox(win, lbl + " scale", val.scale);
    return edited;
}

bool draw_imgui_scene_inspector(glwindow* win, material* val, scene* scn) {
    auto edited = 0;
    edited += draw_imgui_text(win, "name", val->name);
    edited += draw_imgui_checkbox(win, "double_sided", val->double_sided);
    edited +=
        draw_imgui_combobox(win, "type", val->type, material_type_names());
    edited += draw_hdr_color_widget(win, "ke", val->ke);
    edited += draw_imgui_colorbox(win, "kd", val->kd);
    edited += draw_imgui_colorbox(win, "ks", val->ks);
    edited += draw_imgui_colorbox(win, "kr", val->kr);
    edited += draw_imgui_colorbox(win, "kt", val->kt);
    edited += draw_imgui_dragbox(win, "rs", val->rs);
    edited += draw_imgui_dragbox(win, "op", val->op);
    edited += draw_imgui_scene_inspector(win, "ke", val->ke_txt, scn);
    edited += draw_imgui_scene_inspector(win, "kd", val->kd_txt, scn);
    edited += draw_imgui_scene_inspector(win, "ks", val->ks_txt, scn);
    edited += draw_imgui_scene_inspector(win, "kr", val->kr_txt, scn);
    edited += draw_imgui_scene_inspector(win, "kt", val->kt_txt, scn);
    edited += draw_imgui_scene_inspector(win, "rs", val->rs_txt, scn);
    edited += draw_imgui_scene_inspector(win, "op", val->op_txt, scn);
    edited += draw_imgui_scene_inspector(win, "bump", val->bump_txt, scn);
    edited += draw_imgui_scene_inspector(win, "disp", val->disp_txt, scn);
    edited += draw_imgui_scene_inspector(win, "norm", val->norm_txt, scn);
    edited += draw_imgui_scene_inspector(win, "occ", val->occ_txt, scn);
    return edited;
}

bool draw_imgui_scene_inspector(glwindow* win, shape* val, scene* scn) {
    auto edited = 0;
    edited += draw_imgui_text(win, "name", val->name);
    edited += draw_imgui_text(win, "path", val->path);
    draw_imgui_label(win, "points", val->points);
    draw_imgui_label(win, "lines", val->lines);
    draw_imgui_label(win, "triangles", val->triangles);
    draw_imgui_label(win, "quads", val->quads);
    draw_imgui_label(win, "quads_pos", val->quads_pos);
    draw_imgui_label(win, "quads_norm", val->quads_norm);
    draw_imgui_label(win, "quads_texcoord", val->quads_texcoord);
    draw_imgui_label(win, "beziers", val->beziers);
    draw_imgui_label(win, "pos", val->pos);
    draw_imgui_label(win, "norm", val->norm);
    draw_imgui_label(win, "texcoord", val->texcoord);
    draw_imgui_label(win, "texcoord1", val->texcoord1);
    draw_imgui_label(win, "color", val->color);
    draw_imgui_label(win, "radius", val->radius);
    draw_imgui_label(win, "tangsp", val->tangsp);
    draw_imgui_label(win, "subdivision", val->subdivision);
    draw_imgui_label(win, "catmullclark", val->catmullclark);
    return edited;
}

bool draw_imgui_scene_inspector(glwindow* win, instance* val, scene* scn) {
    auto edited = 0;
    edited += draw_imgui_text(win, "name", val->name);
    edited += draw_imgui_dragbox(win, "frame", val->frame);
    edited += draw_imgui_combobox(win, "shp", val->shp, scn->shapes);
    edited += draw_imgui_combobox(win, "mat", val->mat, scn->materials);
    return edited;
}

bool draw_imgui_scene_inspector(glwindow* win, environment* val, scene* scn) {
    auto edited = 0;
    edited += draw_imgui_text(win, "name", val->name);
    edited += draw_imgui_dragbox(win, "frame", val->frame);
    edited += draw_hdr_color_widget(win, "ke", val->ke);
    edited += draw_imgui_scene_inspector(win, "ke", val->ke_txt, scn);
    return edited;
}

bool draw_imgui_scene_inspector(glwindow* win, node* val, scene* scn) {
    auto edited = 0;
    edited += draw_imgui_text(win, "name", val->name);
    edited += draw_imgui_combobox(win, "parent", val->parent, scn->nodes);
    edited += draw_imgui_dragbox(win, "frame", val->frame);
    edited += draw_imgui_dragbox(win, "translation", val->translation);
    edited += draw_imgui_dragbox(win, "rotation", val->rotation, -1, 1);
    edited += draw_imgui_dragbox(win, "scale", val->scale, 0, 10);
    edited += draw_imgui_combobox(win, "cam", val->cam, scn->cameras);
    edited += draw_imgui_combobox(win, "ist", val->ist, scn->instances);
    edited += draw_imgui_combobox(win, "env", val->env, scn->environments);
    return edited;
}

bool draw_imgui_scene_inspector(glwindow* win, animation* val, scene* scn) {
    auto edited = 0;
    edited += draw_imgui_text(win, "name", val->name);
    edited += draw_imgui_text(win, "path", val->path);
    edited += draw_imgui_text(win, "group", val->group);
    edited +=
        draw_imgui_combobox(win, "type", val->type, animation_type_names());
    draw_imgui_label(win, "times", val->times);
    draw_imgui_label(win, "translation", val->translation);
    draw_imgui_label(win, "rotation", val->rotation);
    draw_imgui_label(win, "scale", val->scale);
    draw_imgui_label(win, "weights", val->weights);
    draw_imgui_label(win, "targets", val->targets);
    return edited;
}

bool draw_imgui_scene_tree(glwindow* win, const std::string& lbl, scene* scn,
    scene_selection& sel, std::vector<ygl::scene_selection>& update_list,
    const std::unordered_map<std::string, std::string>& inspector_highlights) {
    if (!scn) return false;
    push_imgui_groupid(win, scn);
    // begin_imgui_scrollarea(win, "scene #$%^!@", 240, false);
    draw_imgui_scene_tree(win, scn, sel, inspector_highlights);
    // end_imgui_scrollarea(win);

    auto update_len = update_list.size();
#if 0
    if (test_scn) {
        draw_add_elem_widgets(
            win, scn, "cam", scn->cameras, test_scn->cameras, sel, update_list);
        draw_add_elem_widgets(win, scn, "txt", scn->textures,
            test_scn->textures, sel, update_list);
        draw_add_elem_widgets(win, scn, "mat", scn->materials,
            test_scn->materials, sel, update_list);
        draw_add_elem_widgets(
            win, scn, "shp", scn->shapes, test_scn->shapes, sel, update_list);
        draw_add_elem_widgets(win, scn, "ist", scn->instances,
            test_scn->instances, sel, update_list);
        draw_add_elem_widgets(
            win, scn, "nde", scn->nodes, test_scn->nodes, sel, update_list);
        draw_add_elem_widgets(win, scn, "env", scn->environments,
            test_scn->environments, sel, update_list);
        draw_add_elem_widgets(win, scn, "anim", scn->animations,
            test_scn->animations, sel, update_list);
    }
#endif

    pop_imgui_groupid(win);
    return update_list.size() != update_len;
}

bool draw_imgui_scene_inspector(glwindow* win, const std::string& lbl,
    scene* scn, scene_selection& sel,
    std::vector<ygl::scene_selection>& update_list,
    const std::unordered_map<std::string, std::string>& inspector_highlights) {
    if (!scn || !sel.get_raw()) return false;
    push_imgui_groupid(win, sel.get_raw());

    auto update_len = update_list.size();

    auto edited = false;
    if (sel.is<camera>())
        edited = draw_imgui_scene_inspector(win, sel.get<camera>(), scn);
    if (sel.is<shape>())
        edited = draw_imgui_scene_inspector(win, sel.get<shape>(), scn);
    if (sel.is<texture>())
        edited = draw_imgui_scene_inspector(win, sel.get<texture>(), scn);
    if (sel.is<material>())
        edited = draw_imgui_scene_inspector(win, sel.get<material>(), scn);
    if (sel.is<environment>())
        edited = draw_imgui_scene_inspector(win, sel.get<environment>(), scn);
    if (sel.is<instance>())
        edited = draw_imgui_scene_inspector(win, sel.get<instance>(), scn);
    if (sel.is<node>())
        edited = draw_imgui_scene_inspector(win, sel.get<node>(), scn);
    if (sel.is<animation>())
        edited = draw_imgui_scene_inspector(win, sel.get<animation>(), scn);
    if (edited) update_list.push_back(sel);

    pop_imgui_groupid(win);
    return update_list.size() != update_len;
}

}  // namespace ygl

#endif
