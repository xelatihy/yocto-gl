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

#if YGL_IMAGEIO
#include "ext/stb_image.h"
#include "ext/stb_image_resize.h"
#include "ext/stb_image_write.h"
#include "ext/tinyexr.h"
#endif

#if YGL_GLTF
#include "ext/json.hpp"
#endif

#if YGL_SVG
#include "ext/nanosvg.h"
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
        theta = clamp(theta, 0.001f, pif - 0.001f);
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
        theta = clamp(theta, 0.001f, pif - 0.001f);
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

    frame = {frame_rot(rotation_frame(vec3f{1, 0, 0}, rotate.y)) *
                 frame_rot(frame) *
                 frame_rot(rotation_frame(vec3f{0, 1, 0}, rotate.x)),
        frame.o + transl.x * x + transl.y * y + transl.z * z};
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
        for (auto vid : l) norm[vid] += n;
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
        for (auto vid : t) norm[vid] += n;
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
        for (auto vid : q) norm[vid] += n;
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
        for (auto vid : t) tangu[vid] += tutv.first;
        for (auto vid : t) tangv[vid] += tutv.second;
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
                quads_pos[fid][c],
                (!quads_norm.empty()) ? quads_norm[fid][c] : -1,
                (!quads_texcoord.empty()) ? quads_texcoord[fid][c] : -1,
            };
            if (vert_map.find(v) == vert_map.end()) {
                auto s = (int)vert_map.size();
                vert_map[v] = s;
            }
            quads[fid][c] = vert_map.at(v);
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

// Generate a rectangular grid of usteps x vsteps uv values for parametric
// surface generation.
void make_quads_uv(std::vector<vec4i>& quads, std::vector<vec2f>& uv,
    int usteps, int vsteps, bool uwrap, bool vwrap, bool vpole0, bool vpole1) {
    auto uvert = (uwrap) ? usteps : usteps + 1;
    auto vvert = (vwrap) ? vsteps : vsteps + 1;
    auto vid = [=](int i, int j) {
        if (uwrap) i = i % usteps;
        if (vwrap) j = j % vsteps;
        return j * uvert + i;
    };

    uv = std::vector<vec2f>(uvert * vvert);
    for (auto j = 0; j < vvert; j++) {
        for (auto i = 0; i < uvert; i++) {
            uv[vid(i, j)] = {i / (float)usteps, j / (float)vsteps};
        }
    }

    quads = std::vector<vec4i>(usteps * vsteps);
    for (auto j = 0; j < vsteps; j++) {
        for (auto i = 0; i < usteps; i++) {
            quads[j * usteps + i] = {
                vid(i, j), vid(i + 1, j), vid(i + 1, j + 1), vid(i, j + 1)};
        }
    }

    if (vpole0) {
        if (vwrap) throw std::runtime_error("cannot have a pole with wrapping");
        uv = std::vector<vec2f>(uv.begin() + uvert, uv.end());
        uv.insert(uv.begin(), {0, 0});
        for (auto& q : quads) {
            for (auto& vid : q) { vid = (vid < usteps) ? 0 : vid - uvert + 1; }
            if (q.x == 0 && q.y == 0) q = {q.z, q.w, q.x, q.y};
        }
    }

    if (vpole1) {
        if (vwrap) throw std::runtime_error("cannot have a pole with wrapping");
        auto pid = (int)uv.size() - uvert;
        uv = std::vector<vec2f>(uv.begin(), uv.end() - uvert);
        uv.insert(uv.end(), {0, 1});
        for (auto& q : quads) {
            for (auto& vid : q) { vid = (vid < pid) ? vid : pid; }
        }
    }
}

// Generate parametric num lines of usteps segments.
void make_lines_uv(
    std::vector<vec2i>& lines, std::vector<vec2f>& uv, int num, int usteps) {
    auto vid = [usteps](int i, int j) { return j * (usteps + 1) + i; };
    uv = std::vector<vec2f>((usteps + 1) * num);
    for (auto j = 0; j < num; j++) {
        for (auto i = 0; i <= usteps; i++) {
            uv[vid(i, j)] = {i / (float)usteps, j / (float)num};
        }
    }

    lines = std::vector<vec2i>(usteps * num);
    for (int j = 0; j < num; j++) {
        for (int i = 0; i < usteps; i++) {
            lines[j * usteps + i] = {vid(i, j), vid(i + 1, j)};
        }
    }
}

// Generate a parametric point set. Mostly here for completeness.
void make_points_uv(std::vector<int>& points, std::vector<vec2f>& uv, int num) {
    uv = std::vector<vec2f>(num);
    for (auto i = 0; i < num; i++) { uv[i] = {i / (float)num, 0}; }

    points = std::vector<int>(num);
    for (auto i = 0; i < num; i++) points[i] = i;
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
            if (contains(emap, e)) continue;
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
            if (contains(emap, e)) continue;
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
        if (!contains(vmap, b.x)) {
            vmap[b.x] = (int)tvert.size();
            tvert.push_back(vert[b.x]);
        }
        if (!contains(vmap, b.w)) {
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
            if (contains(emap, e)) continue;
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
        if (contains(emap, {e.y, e.x})) continue;
        tboundary.push_back({e.x, v});
        tboundary.push_back({v, e.y});
    }

    // setup
    auto tcrease_edges = tboundary;
    auto tcrease_verts = std::vector<int>();

    // define vertex valence ---------------------------
    auto tvert_val = std::vector<int>(tvert.size(), 2);
    for (auto e : tcrease_edges)
        for (auto vid : e) tvert_val[vid] = 1;
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
        for (auto vid : e) {
            if (tvert_val[vid] != 1) continue;
            avert[vid] += c;
            acount[vid] += 1;
        }
    }
    for (auto& q : tquads) {
        auto c = (tvert[q.x] + tvert[q.y] + tvert[q.z] + tvert[q.w]) / 4.0f;
        for (auto vid : q) {
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
inline void merge_lines(
    std::vector<vec2i>& lines, const std::vector<vec2i>& lines1) {
    auto nverts = 0;
    for (auto& l : lines) nverts = max(nverts, max_element(l));
    for (auto& l : lines1) lines.push_back({l.x + nverts, l.y + nverts});
}
// Merge triangles between shapes. The elements are merged by increasing the
// array size of the second array by the number of vertices of the first.
// Vertex data can then be concatenated successfully.
void merge_triangles(
    std::vector<vec3i>& triangles, const std::vector<vec3i>& triangles1) {
    auto nverts = 0;
    for (auto& t : triangles) nverts = max(nverts, max_element(t));
    for (auto& t : triangles1)
        triangles.push_back({t.x + nverts, t.y + nverts, t.z + nverts});
}
// Merge quads between shapes. The elements are merged by increasing the
// array size of the second array by the number of vertices of the first.
// Vertex data can then be concatenated successfully.
void merge_quads(std::vector<vec4i>& quads, const std::vector<vec4i>& quads1) {
    auto nverts = 0;
    for (auto& q : quads) nverts = max(nverts, max_element(q));
    for (auto& q : quads1)
        quads.push_back(
            {q.x + nverts, q.y + nverts, q.z + nverts, q.w + nverts});
}

// Duplicate vertex data for each line index, giving a faceted look.
template <typename T>
void facet_lines(
    std::vector<vec2i>& lines, std::vector<T>& vert, bool update_lines) {
    if (vert.empty()) return;
    auto tvert = vert;
    tvert.resize(lines.size() * 2);
    for (auto i = 0; i < lines.size(); i++) {
        tvert[2 * i + 0] = vert[lines[i].x];
        tvert[2 * i + 1] = vert[lines[i].y];
    }
    std::swap(vert, tvert);

    if (update_lines) {
        auto tlines = std::vector<vec2i>(lines.size());
        for (auto i = 0; i < lines.size(); i++)
            tlines[i] = {i * 2 + 0, i * 2 + 1};
        std::swap(lines, tlines);
    }
}

// Duplicate vertex data for each line index, giving a faceted look.
template <typename T>
void facet_triangles(std::vector<vec3i>& triangles, std::vector<T>& vert,
    bool update_triangles) {
    if (vert.empty()) return;
    auto tvert = vert;
    tvert.resize(triangles.size() * 3);
    for (auto i = 0; i < triangles.size(); i++) {
        tvert[3 * i + 0] = vert[triangles[i].x];
        tvert[3 * i + 1] = vert[triangles[i].y];
        tvert[3 * i + 2] = vert[triangles[i].z];
    }
    std::swap(vert, tvert);

    if (update_triangles) {
        auto ttriangles = std::vector<vec3i>(triangles.size());
        for (auto i = 0; i < triangles.size(); i++)
            ttriangles[i] = {i * 3 + 0, i * 3 + 1, i * 3 + 2};
        std::swap(triangles, ttriangles);
    }
}
// Duplicate vertex data for each quad index, giving a faceted look.
template <typename T>
void facet_quads(
    std::vector<vec4i>& quads, std::vector<T>& vert, bool update_quads) {
    if (vert.empty()) return;
    auto tvert = vert;
    tvert.resize(quads.size() * 4);
    for (auto i = 0; i < quads.size(); i++) {
        tvert[4 * i + 0] = vert[quads[i].x];
        tvert[4 * i + 1] = vert[quads[i].y];
        tvert[4 * i + 2] = vert[quads[i].z];
        tvert[4 * i + 3] = vert[quads[i].z];
    }
    std::swap(vert, tvert);

    if (update_quads) {
        auto tquads = std::vector<vec4i>(quads.size());
        for (auto i = 0; i < quads.size(); i++)
            tquads[i] = {i * 4 + 0, i * 4 + 1, i * 4 + 2, i * 4 + 3};
        std::swap(quads, tquads);
    }
}

// Unshare line data.
void facet_lines(std::vector<vec2i>& lines, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<vec4f>& color, std::vector<float>& radius) {
    facet_lines(lines, norm, false);
    facet_lines(lines, texcoord, false);
    facet_lines(lines, color, false);
    facet_lines(lines, radius, false);
    facet_lines(lines, pos);
}

// Unshare triangle data.
void facet_triangles(std::vector<vec3i>& triangles, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<vec4f>& color, std::vector<float>& radius) {
    facet_triangles(triangles, norm, false);
    facet_triangles(triangles, texcoord, false);
    facet_triangles(triangles, color, false);
    facet_triangles(triangles, radius, false);
    facet_triangles(triangles, pos);
}

// Unshare quad data.
void facet_quads(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    std::vector<vec4f>& color, std::vector<float>& radius) {
    facet_quads(quads, norm, false);
    facet_quads(quads, texcoord, false);
    facet_quads(quads, color, false);
    facet_quads(quads, radius, false);
    facet_quads(quads, pos);
}

// Unshare vertices for faceting
template <typename T>
std::vector<T> facet_vert(
    const std::vector<T>& vert, const std::vector<int>& vmap) {
    if (vert.empty()) return vert;
    auto tvert = std::vector<T>(vmap.size());
    for (auto vid = 0; vid < vmap.size(); vid++) tvert[vid] = vert[vmap[vid]];
    return tvert;
}

// instantations
template std::vector<float> facet_vert<float>(
    const std::vector<float>& vert, const std::vector<int>& vmap);
template std::vector<vec2f> facet_vert<vec2f>(
    const std::vector<vec2f>& vert, const std::vector<int>& vmap);
template std::vector<vec3f> facet_vert<vec3f>(
    const std::vector<vec3f>& vert, const std::vector<int>& vmap);
template std::vector<vec4f> facet_vert<vec4f>(
    const std::vector<vec4f>& vert, const std::vector<int>& vmap);

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHAPE SAMPLING
// -----------------------------------------------------------------------------
namespace ygl {

// Compute a distribution for sampling points uniformly
std::vector<float> sample_points_cdf(int npoints) {
    auto cdf = std::vector<float>(npoints);
    for (auto i = 0; i < npoints; i++) cdf[i] = i + 1;
    return cdf;
}

// Compute a distribution for sampling lines uniformly
std::vector<float> sample_lines_cdf(
    const std::vector<vec2i>& lines, const std::vector<vec3f>& pos) {
    auto cdf = std::vector<float>(lines.size());
    for (auto i = 0; i < lines.size(); i++)
        cdf[i] = length(pos[lines[i].x] - pos[lines[i].y]);
    for (auto i = 1; i < lines.size(); i++) cdf[i] += cdf[i - 1];
    return cdf;
}

// Compute a distribution for sampling triangle meshes uniformly
std::vector<float> sample_triangles_cdf(
    const std::vector<vec3i>& triangles, const std::vector<vec3f>& pos) {
    auto cdf = std::vector<float>(triangles.size());
    for (auto i = 0; i < triangles.size(); i++)
        cdf[i] = triangle_area(
            pos[triangles[i].x], pos[triangles[i].y], pos[triangles[i].z]);
    for (auto i = 1; i < triangles.size(); i++) cdf[i] += cdf[i - 1];
    return cdf;
}

// Compute a distribution for sampling quad meshes uniformly
std::vector<float> sample_quads_cdf(
    const std::vector<vec4i>& quads, const std::vector<vec3f>& pos) {
    auto cdf = std::vector<float>(quads.size());
    for (auto i = 0; i < quads.size(); i++)
        cdf[i] = quad_area(
            pos[quads[i].x], pos[quads[i].y], pos[quads[i].z], pos[quads[i].w]);
    for (auto i = 1; i < quads.size(); i++) cdf[i] += cdf[i - 1];
    return cdf;
}

// Samples a set of points over a triangle mesh uniformly. The rng function
// takes the point index and returns vec3f numbers uniform directibuted in
// [0,1]^3. unorm and texcoord are optional.
std::tuple<std::vector<vec3f>, std::vector<vec3f>, std::vector<vec2f>>
sample_triangles_points(const std::vector<vec3i>& triangles,
    const std::vector<vec3f>& pos, const std::vector<vec3f>& norm,
    const std::vector<vec2f>& texcoord, int npoints, uint64_t seed) {
    auto sampled_pos = std::vector<vec3f>(npoints);
    auto sampled_norm = std::vector<vec3f>(norm.empty() ? 0 : npoints);
    auto sampled_texcoord = std::vector<vec2f>(texcoord.empty() ? 0 : npoints);
    auto cdf = sample_triangles_cdf(triangles, pos);
    auto rng = init_rng(seed);
    for (auto i = 0; i < npoints; i++) {
        auto eid = 0;
        auto euv = zero2f;
        std::tie(eid, euv) = sample_triangles(
            cdf, next_rand1f(rng), {next_rand1f(rng), next_rand1f(rng)});
        auto t = triangles[eid];
        sampled_pos[i] =
            interpolate_triangle(pos[t.x], pos[t.y], pos[t.z], euv);
        if (!sampled_norm.empty())
            sampled_norm[i] = normalize(
                interpolate_triangle(norm[t.x], norm[t.y], norm[t.z], euv));
        if (!sampled_texcoord.empty())
            sampled_texcoord[i] = interpolate_triangle(
                texcoord[t.x], texcoord[t.y], texcoord[t.z], euv);
    }

    return {sampled_pos, sampled_norm, sampled_texcoord};
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR IMAGEIO
// -----------------------------------------------------------------------------
namespace ygl {

// check hdr extensions
bool is_hdr_filename(const std::string& filename) {
    auto ext = path_extension(filename);
    return ext == ".hdr" || ext == ".exr";
}

// Loads an ldr image.
image4b load_image4b(const std::string& filename) {
    auto w = 0, h = 0, c = 0;
    auto pixels =
        std::unique_ptr<byte>(stbi_load(filename.c_str(), &w, &h, &c, 4));
    if (!pixels) return {};
    return make_image(w, h, (vec4b*)pixels.get());
}

// Loads an hdr image.
image4f load_image4f(const std::string& filename) {
    auto ext = path_extension(filename);
    auto w = 0, h = 0, c = 0;
    auto pixels = std::unique_ptr<float>(nullptr);
    if (ext == ".exr") {
        auto pixels_ = (float*)nullptr;
        if (!LoadEXR(&pixels_, &w, &h, filename.c_str(), nullptr))
            pixels = std::unique_ptr<float>(pixels_);
    } else {
        pixels =
            std::unique_ptr<float>(stbi_loadf(filename.c_str(), &w, &h, &c, 4));
    }
    if (!pixels) return {};
    return make_image(w, h, (vec4f*)pixels.get());
}

// Saves an ldr image.
bool save_image4b(const std::string& filename, const image4b& img) {
    if (path_extension(filename) == ".png") {
        return stbi_write_png(filename.c_str(), img.width(), img.height(), 4,
            (byte*)data(img), img.width() * 4);
    } else if (path_extension(filename) == ".jpg") {
        return stbi_write_jpg(filename.c_str(), img.width(), img.height(), 4,
            (byte*)data(img), 75);
    } else {
        return false;
    }
}

// Saves an hdr image.
bool save_image4f(const std::string& filename, const image4f& img) {
    if (path_extension(filename) == ".hdr") {
        return stbi_write_hdr(
            filename.c_str(), img.width(), img.height(), 4, (float*)data(img));
    } else if (path_extension(filename) == ".exr") {
        return !SaveEXR(
            (float*)data(img), img.width(), img.height(), 4, filename.c_str());
    } else {
        return false;
    }
}

// Loads an image
std::vector<float> load_imagef(
    const std::string& filename, int& width, int& height, int& ncomp) {
    auto pixels = stbi_loadf(filename.c_str(), &width, &height, &ncomp, 0);
    if (!pixels) return {};
    auto ret = std::vector<float>(pixels, pixels + width * height * ncomp);
    free(pixels);
    return ret;
}

// Loads an image
std::vector<byte> load_image(
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
std::vector<byte> load_image_from_memory(const std::string& filename,
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
    } else {
        return false;
    }
}

// Saves an image
bool save_image(const std::string& filename, int width, int height, int ncomp,
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
    const tonemap_params& params) {
    if (is_hdr_filename(filename)) {
        return save_image4f(filename, hdr);
    } else {
        auto ldr = tonemap_image(hdr, params);
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

    stbir_resize_float_generic((float*)img.pixels.data(), img.width(),
        img.height(), sizeof(vec4f) * img.width(),
        (float*)res_img.pixels.data(), res_img.width(), res_img.height(),
        sizeof(vec4f) * res_img.width(), 4, 3,
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

    stbir_resize_uint8_generic((unsigned char*)img.pixels.data(), img.width(),
        img.height(), sizeof(vec4b) * img.width(),
        (unsigned char*)res_img.pixels.data(), res_img.width(),
        res_img.height(), sizeof(vec4b) * res_img.width(), 4, 3,
        (premultiplied_alpha) ? STBIR_FLAG_ALPHA_PREMULTIPLIED : 0,
        edge_map.at(edge), filter_map.at(filter), STBIR_COLORSPACE_LINEAR,
        nullptr);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF IMAGE OPERATIONS
// -----------------------------------------------------------------------------
namespace ygl {

#if 1
// Tone map with a fitted filmic curve.
//
// Implementation from
// https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
inline float tonemap_filmic(float hdr) {
    // rescale
    auto x = hdr * 2.05f;
    // fitted values
    float a = 2.51f, b = 0.03f, c = 2.43f, d = 0.59f, e = 0.14f;
    auto y = ((x * (a * x + b)) / (x * (c * x + d) + e));
    return pow(clamp(y, 0.0f, 1.0f), 1 / 2.2f);
}
#else
inline float tonemap_filmic(float x) {
    auto y =
        (x * (x * (x * (x * 2708.7142 + 6801.1525) + 1079.5474) + 1.1614649) -
            0.00004139375) /
        (x * (x * (x * (x * 983.38937 + 4132.0662) + 2881.6522) + 128.35911) +
            1.0);
    return (float)std::max(y, 0.0);
}
#endif

// Tone mapping HDR to LDR images.
image4b tonemap_image(const image4f& hdr, const tonemap_params& params) {
    auto ldr = image4b(hdr.width(), hdr.height());
    auto scale = pow(2.0f, params.exposure);
    auto inv_gamma = 1 / params.gamma;
    for (auto j = 0; j < hdr.height(); j++) {
        for (auto i = 0; i < hdr.width(); i++) {
            auto h = hdr.at(i, j) * vec4f{scale, scale, scale, 1};
            if (params.filmic) {
                h = {tonemap_filmic(h.x), tonemap_filmic(h.y),
                    tonemap_filmic(h.z), h.w};
            } else {
                h = {pow(h.x, inv_gamma), pow(h.y, inv_gamma),
                    pow(h.z, inv_gamma), h.w};
            }
            ldr.at(i, j) = float_to_byte(h);
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
            auto w = byte_to_float(layers[l][i].w);
            comp.x += byte_to_float(layers[l][i].x) * w * weight;
            comp.y += byte_to_float(layers[l][i].y) * w * weight;
            comp.z += byte_to_float(layers[l][i].z) * w * weight;
            comp.w += w * weight;
            weight *= (1 - w);
        }
        if (comp.w) {
            img[i].x = float_to_byte(comp.x / comp.w);
            img[i].y = float_to_byte(comp.y / comp.w);
            img[i].z = float_to_byte(comp.z / comp.w);
            img[i].w = float_to_byte(comp.w);
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
    image4b pixels(width, height);
    for (int j = 0; j < width; j++) {
        for (int i = 0; i < height; i++) {
            auto c = i % tile == 0 || i % tile == tile - 1 || j % tile == 0 ||
                     j % tile == tile - 1;
            pixels.at(i, j) = (c) ? c0 : c1;
        }
    }
    return pixels;
}

// Make a checkerboard image
image4b make_checker_image(
    int width, int height, int tile, const vec4b& c0, const vec4b& c1) {
    image4b pixels(width, height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto c = (i / tile + j / tile) % 2 == 0;
            pixels.at(i, j) = (c) ? c0 : c1;
        }
    }
    return pixels;
}

// Make an image with bumps and dimples.
image4b make_bumpdimple_image(int width, int height, int tile) {
    image4b pixels(width, height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto c = (i / tile + j / tile) % 2 == 0;
            auto ii = i % tile - tile / 2, jj = j % tile - tile / 2;
            auto r =
                sqrt(float(ii * ii + jj * jj)) / sqrt(float(tile * tile) / 4);
            auto h = 0.5f;
            if (r < 0.5f) { h += (c) ? (0.5f - r) : -(0.5f - r); }
            auto g = float_to_byte(h);
            pixels.at(i, j) = vec4b{g, g, g, 255};
        }
    }
    return pixels;
}

// Make a uv colored grid
image4b make_ramp_image(
    int width, int height, const vec4b& c0, const vec4b& c1, bool srgb) {
    image4b pixels(width, height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto u = (float)i / (float)width;
            if (srgb) {
                pixels.at(i, j) = linear_to_srgb(
                    srgb_to_linear(c0) * (1 - u) + srgb_to_linear(c1) * u);
            } else {
                pixels.at(i, j) = float_to_byte(
                    byte_to_float(c0) * (1 - u) + byte_to_float(c1) * u);
            }
        }
    }
    return pixels;
}

// Make a gamma ramp image
image4b make_gammaramp_image(int width, int height) {
    image4b pixels(width, height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto u = j / float(height - 1);
            if (i < width / 3) u = pow(u, 2.2f);
            if (i > (width * 2) / 3) u = pow(u, 1 / 2.2f);
            auto c = (unsigned char)(u * 255);
            pixels.at(i, j) = {c, c, c, 255};
        }
    }
    return pixels;
}

// Make a gamma ramp image
image4f make_gammaramp_imagef(int width, int height) {
    image4f pixels(width, height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto u = j / float(height - 1);
            if (i < width / 3) u = pow(u, 2.2f);
            if (i > (width * 2) / 3) u = pow(u, 1 / 2.2f);
            pixels.at(i, j) = {u, u, u, 1};
        }
    }
    return pixels;
}

// Make an image color with red/green in the [0,1] range. Helpful to visualize
// uv texture coordinate application.
image4b make_uv_image(int width, int height) {
    image4b pixels(width, height);
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            auto r = float_to_byte(i / (float)(width - 1));
            auto g = float_to_byte(j / (float)(height - 1));
            pixels.at(i, j) = vec4b{r, g, 0, 255};
        }
    }
    return pixels;
}

// Make a uv colored grid
image4b make_uvgrid_image(int width, int height, int tile, bool colored) {
    image4b pixels(width, height);
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
            pixels.at(i, j) = (colored) ? hsv_to_rgb({ph, ps, pv, 255}) :
                                          vec4b{pv, pv, pv, 255};
        }
    }
    return pixels;
}

// Make a uv recusive colored grid
image4b make_recuvgrid_image(int width, int height, int tile, bool colored) {
    image4b pixels(width, height);
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
            pixels.at(i, j) = (colored) ? hsv_to_rgb({ph, ps, pv, 255}) :
                                          vec4b{pv, pv, pv, 255};
        }
    }
    return pixels;
}

// Comvert a bump map to a normal map.
image4b bump_to_normal_map(const image4b& img, float scale) {
    image4b norm(img.width(), img.height());
    for (int j = 0; j < img.height(); j++) {
        for (int i = 0; i < img.width(); i++) {
            auto i1 = (i + 1) % img.width(), j1 = (j + 1) % img.height();
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
                tan((4.0f / 9.0f - T / 120.0f) * (pif - 2 * t1)) -
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
        auto tauR = exp(-sun_m * 0.008735f * pow(sun_lambda[i] / 1000, -4.08f));
        auto tauA = exp(-sun_m * sun_beta * pow(sun_lambda[i] / 1000, -1.3f));
        auto tauO = exp(-sun_m * sun_ko[i] * .35f);
        auto tauG = exp(-1.41f * sun_kg[i] * sun_m /
                        pow(1 + 118.93f * sun_kg[i] * sun_m, 0.45f));
        auto tauWA = exp(-0.2385f * sun_kwa[i] * 2.0f * sun_m /
                         pow(1 + 20.07f * sun_kwa[i] * 2.0f * sun_m, 0.45f));
        sun_le[i] = sun_sol[i] * tauR * tauA * tauO * tauG * tauWA;
    }

    auto sun = [has_sun, sunAngularRadius, sun_le](auto theta, auto gamma) {
        return (has_sun && gamma < sunAngularRadius) ? sun_le / 10000.0f :
                                                       zero3f;
    };

    auto img = image4f(2 * res, res);
    for (auto j = 0; j < img.height(); j++) {
        if (!has_ground && j >= img.height() / 2) continue;
        auto theta = pif * ((j + 0.5f) / img.height());
        theta = clamp(theta, 0.0f, pif / 2 - flt_eps);
        for (int i = 0; i < img.width(); i++) {
            auto phi = 2 * pif * (float(i + 0.5f) / img.width());
            auto w =
                vec3f(cos(phi) * sin(theta), cos(theta), sin(phi) * sin(theta));
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
    auto img = image4b(resx, resy);
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
    auto img = image4b(resx, resy);
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
    auto img = image4b(resx, resy);
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
    auto img = image4b(resx, resy);
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

// Intersect a ray with a axis-aligned bounding box
bool intersect_check_bbox(const ray3f& ray, const bbox3f& bbox) {
    // set up convenient pointers for looping over axes
    auto tmin = ray.tmin, tmax = ray.tmax;

    // for each axis, clip intersection against the bounding planes
    for (int i = 0; i < 3; i++) {
        // determine intersection ranges
        auto invd = 1.0f / ray.d[i];
        auto t0 = (bbox.min[i] - ray.o[i]) * invd;
        auto t1 = (bbox.max[i] - ray.o[i]) * invd;
        // flip based on range directions
        if (invd < 0.0f) {
            float a = t0;
            t0 = t1;
            t1 = a;
        }
        // clip intersection
        tmin = t0 > tmin ? t0 : tmin;
        tmax = t1 < tmax ? t1 : tmax;
        // if intersection is empty, exit
        if (tmin > tmax) return false;
    }

    // passed all planes, then intersection occurred
    return true;
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
    auto p = lerp(v0, v1, u);
    auto r = lerp(r0, r1, u);
    auto d2 = dot(pos - p, pos - p);
    // check distance
    if (d2 > (dist_max + r) * (dist_max + r)) return false;
    // done
    dist = sqrt(d2);
    euv = {u, 0};
    return true;
}

// TODO: documentation
// this is a complicated test -> I probably prefer to use a sequence of test
// (triangle body, and 3 edges)
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
    for (int i = 0; i < 3; i++) {
        auto v = pos[i];
        if (v < bbox.min[i]) dd += (bbox.min[i] - v) * (bbox.min[i] - v);
        if (v > bbox.max[i]) dd += (v - bbox.max[i]) * (v - bbox.max[i]);
    }

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

// cleanup
bvh_tree::~bvh_tree() {
    if (!own_shape_bvhs) return;
    for (auto bvh : shape_bvhs) delete bvh;
}

// number of primitives to avoid splitting on
const int bvh_minprims = 4;

// Initializes the BVH node node that contains the primitives sorted_prims
// from start to end, by either splitting it into two other nodes,
// or initializing it as a leaf. When splitting, the heuristic heuristic is
// used and nodes added sequentially in the preallocated nodes array and
// the number of nodes nnodes is updated.
void make_bvh_node(std::vector<bvh_node>& nodes, int nodeid,
    std::vector<int>& sorted_prims, int start, int end,
    const std::vector<bbox3f>& bboxes, bvh_node_type type, bool equal_size) {
    // compute node bounds
    auto& node = nodes.at(nodeid);
    node.bbox = invalid_bbox3f;
    for (auto i = start; i < end; i++) node.bbox += bboxes[sorted_prims[i]];

    // initialize as a leaf
    node.type = type;
    node.start = start;
    node.count = end - start;

    // try to split into two children
    if (end - start > bvh_minprims) {
        // choose the split axis and position
        // init to default values
        auto axis = 0;
        auto mid = (start + end) / 2;

        // compute primintive bounds and size
        auto centroid_bbox = invalid_bbox3f;
        for (auto i = start; i < end; i++)
            centroid_bbox += bbox_center(bboxes[sorted_prims[i]]);
        auto centroid_size = bbox_diagonal(centroid_bbox);

        // check if it is not possible to split
        if (centroid_size != zero3f) {
            // split along largest
            auto largest_axis = max_element(centroid_size);

            // check heuristic
            if (equal_size) {
                // split the space in the middle along the largest axis
                axis = largest_axis;
                auto middle = bbox_center(centroid_bbox)[largest_axis];
                mid =
                    (int)(std::partition(sorted_prims.data() + start,
                              sorted_prims.data() + end,
                              [axis, middle, &bboxes](auto& a) {
                                  return bbox_center(bboxes[a])[axis] < middle;
                              }) -
                          sorted_prims.data());
            } else {
                // balanced tree split: find the largest axis of the bounding
                // box and split along this one right in the middle
                axis = largest_axis;
                mid = (start + end) / 2;
                std::nth_element(sorted_prims.data() + start,
                    sorted_prims.data() + mid, sorted_prims.data() + end,
                    [axis, &bboxes](auto& a, auto& b) {
                        return bbox_center(bboxes[a])[axis] <
                               bbox_center(bboxes[b])[axis];
                    });
            }

            // check correctness
            assert(axis >= 0 && mid > 0);
            assert(mid > start && mid < end);

            // makes an internal node
            node.type = bvh_node_type::internal;
            // perform the splits by preallocating the child nodes and recurring
            node.axis = axis;
            node.start = (int)nodes.size();
            node.count = 2;
            nodes.emplace_back();
            nodes.emplace_back();
            // build child nodes
            make_bvh_node(nodes, node.start, sorted_prims, start, mid, bboxes,
                type, equal_size);
            make_bvh_node(nodes, node.start + 1, sorted_prims, mid, end, bboxes,
                type, equal_size);
        }
    }
}

// Build a BVH node list and sorted primitive array
std::tuple<std::vector<bvh_node>, std::vector<int>> make_bvh_nodes(
    const std::vector<bbox3f>& bboxes, bvh_node_type type, bool equal_size) {
    // create an array of primitives to sort
    auto sorted_prim = std::vector<int>(bboxes.size());
    for (auto i = 0; i < bboxes.size(); i++) sorted_prim[i] = i;

    // allocate nodes (over-allocate now then shrink)
    auto nodes = std::vector<bvh_node>();
    nodes.reserve(sorted_prim.size() * 2);

    // start recursive splitting
    nodes.emplace_back();
    make_bvh_node(nodes, 0, sorted_prim, 0, (int)sorted_prim.size(), bboxes,
        type, equal_size);

    // shrink back
    nodes.shrink_to_fit();

    // done
    return {nodes, sorted_prim};
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
            bboxes.push_back(transform_bbox(ist.frame, ist.bvh->nodes[0].bbox));
        }
        bvh->type = bvh_node_type::instance;
    }

    // make node bvh
    std::tie(bvh->nodes, bvh->sorted_prim) =
        make_bvh_nodes(bboxes, bvh->type, equal_size);

    // sort primitives
    auto sort_prims = [bvh](auto& prims) {
        if (prims.empty()) return;
        auto sprims = prims;
        for (auto i = 0; i < bvh->sorted_prim.size(); i++) {
            prims[i] = sprims[bvh->sorted_prim[i]];
        }
    };
    sort_prims(bvh->points);
    sort_prims(bvh->lines);
    sort_prims(bvh->triangles);
    sort_prims(bvh->quads);
    sort_prims(bvh->instances);
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
    const std::vector<bvh_tree*>& shape_bvhs, bool own_shape_bvhs,
    bool equal_size) {
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
            for (auto i = node.start; i < node.start + node.count; i++) {
                refit_bvh(bvh, i);
                node.bbox += bvh->nodes[i].bbox;
            }
        } break;
        case bvh_node_type::point: {
            for (auto i = node.start; i < node.start + node.count; i++) {
                auto& p = bvh->points[i];
                node.bbox += point_bbox(bvh->pos[p], bvh->radius[p]);
            }
        } break;
        case bvh_node_type::line: {
            for (auto i = node.start; i < node.start + node.count; i++) {
                auto& l = bvh->lines[i];
                node.bbox += line_bbox(bvh->pos[l.x], bvh->pos[l.y],
                    bvh->radius[l.x], bvh->radius[l.y]);
            }
        } break;
        case bvh_node_type::triangle: {
            for (auto i = node.start; i < node.start + node.count; i++) {
                auto& t = bvh->triangles[i];
                node.bbox +=
                    triangle_bbox(bvh->pos[t.x], bvh->pos[t.y], bvh->pos[t.z]);
            }
        } break;
        case bvh_node_type::quad: {
            for (auto i = node.start; i < node.start + node.count; i++) {
                auto& q = bvh->quads[i];
                node.bbox += quad_bbox(
                    bvh->pos[q.x], bvh->pos[q.y], bvh->pos[q.z], bvh->pos[q.w]);
            }
        } break;
        case bvh_node_type::vertex: {
            for (auto i = node.start; i < node.start + node.count; i++) {
                auto idx = bvh->sorted_prim[i];
                node.bbox += point_bbox(bvh->pos[idx], bvh->radius[idx]);
            }
        } break;
        case bvh_node_type::instance: {
            for (auto i = node.start; i < node.start + node.count; i++) {
                auto& ist = bvh->instances[i];
                node.bbox += transform_bbox(ist.frame, ist.bvh->nodes[0].bbox);
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
        bvh->instances[i].frame = frames[bvh->sorted_prim[i]];
        bvh->instances[i].frame_inv = frames_inv[bvh->sorted_prim[i]];
    }
    refit_bvh(bvh, 0);
}

// Intersect ray with a bvh.
bool intersect_bvh(const bvh_tree* bvh, const ray3f& ray_, bool find_any,
    float& ray_t, int& iid, int& sid, int& eid, vec2f& euv) {
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
                if (ray_reverse[node.axis]) {
                    node_stack[node_cur++] = node.start;
                    node_stack[node_cur++] = node.start + 1;
                } else {
                    node_stack[node_cur++] = node.start + 1;
                    node_stack[node_cur++] = node.start;
                }
            } break;
            case bvh_node_type::point: {
                for (auto i = node.start; i < node.start + node.count; i++) {
                    auto& p = bvh->points[i];
                    if (intersect_point(
                            ray, bvh->pos[p], bvh->radius[p], ray_t)) {
                        hit = true;
                        ray.tmax = ray_t;
                        eid = bvh->sorted_prim[i];
                        euv = {1, 0};
                    }
                }
            } break;
            case bvh_node_type::line: {
                for (auto i = node.start; i < node.start + node.count; i++) {
                    auto& l = bvh->lines[i];
                    if (intersect_line(ray, bvh->pos[l.x], bvh->pos[l.y],
                            bvh->radius[l.x], bvh->radius[l.y], ray_t, euv)) {
                        hit = true;
                        ray.tmax = ray_t;
                        eid = bvh->sorted_prim[i];
                    }
                }
            } break;
            case bvh_node_type::triangle: {
                for (auto i = node.start; i < node.start + node.count; i++) {
                    auto& t = bvh->triangles[i];
                    if (intersect_triangle(ray, bvh->pos[t.x], bvh->pos[t.y],
                            bvh->pos[t.z], ray_t, euv)) {
                        hit = true;
                        ray.tmax = ray_t;
                        eid = bvh->sorted_prim[i];
                    }
                }
            } break;
            case bvh_node_type::quad: {
                for (auto i = node.start; i < node.start + node.count; i++) {
                    auto& q = bvh->quads[i];
                    if (intersect_quad(ray, bvh->pos[q.x], bvh->pos[q.y],
                            bvh->pos[q.z], bvh->pos[q.w], ray_t, euv)) {
                        hit = true;
                        ray.tmax = ray_t;
                        eid = bvh->sorted_prim[i];
                    }
                }
            } break;
            case bvh_node_type::vertex: {
                for (auto i = node.start; i < node.start + node.count; i++) {
                    auto idx = bvh->sorted_prim[i];
                    if (intersect_point(
                            ray, bvh->pos[idx], bvh->radius[idx], ray_t)) {
                        hit = true;
                        ray.tmax = ray_t;
                        eid = idx;
                        euv = {1, 0};
                    }
                }
            } break;
            case bvh_node_type::instance: {
                for (auto i = node.start; i < node.start + node.count; i++) {
                    auto& ist = bvh->instances[i];
                    if (intersect_bvh(ist.bvh,
                            transform_ray(ist.frame_inv, ray), find_any, ray_t,
                            iid, sid, eid, euv)) {
                        hit = true;
                        ray.tmax = ray_t;
                        iid = ist.iid;
                        sid = ist.sid;
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
    bool find_any, float& dist, int& iid, int& sid, int& eid, vec2f& euv) {
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
                node_stack[node_cur++] = node.start;
                node_stack[node_cur++] = node.start + 1;
            } break;
            case bvh_node_type::point: {
                for (auto i = node.start; i < node.start + node.count; i++) {
                    auto& p = bvh->points[i];
                    if (overlap_point(
                            pos, max_dist, bvh->pos[p], bvh->radius[p], dist)) {
                        hit = true;
                        max_dist = dist;
                        eid = bvh->sorted_prim[i];
                        euv = {1, 0};
                    }
                }
            } break;
            case bvh_node_type::line: {
                for (auto i = node.start; i < node.start + node.count; i++) {
                    auto& l = bvh->lines[i];
                    if (overlap_line(pos, max_dist, bvh->pos[l.x],
                            bvh->pos[l.y], bvh->radius[l.x], bvh->radius[l.y],
                            dist, euv)) {
                        hit = true;
                        max_dist = dist;
                        eid = bvh->sorted_prim[i];
                    }
                }
            } break;
            case bvh_node_type::triangle: {
                for (auto i = node.start; i < node.start + node.count; i++) {
                    auto& t = bvh->triangles[i];
                    if (overlap_triangle(pos, max_dist, bvh->pos[t.x],
                            bvh->pos[t.y], bvh->pos[t.z], bvh->radius[t.x],
                            bvh->radius[t.y], bvh->radius[t.z], dist, euv)) {
                        hit = true;
                        max_dist = dist;
                        eid = bvh->sorted_prim[i];
                    }
                }
            } break;
            case bvh_node_type::quad: {
                for (auto i = node.start; i < node.start + node.count; i++) {
                    auto& q = bvh->quads[i];
                    if (overlap_quad(pos, max_dist, bvh->pos[q.x],
                            bvh->pos[q.y], bvh->pos[q.z], bvh->pos[q.w],
                            bvh->radius[q.x], bvh->radius[q.y],
                            bvh->radius[q.z], bvh->radius[q.w], dist, euv)) {
                        hit = true;
                        max_dist = dist;
                        eid = bvh->sorted_prim[i];
                    }
                }
            } break;
            case bvh_node_type::vertex: {
                for (auto i = node.start; i < node.start + node.count; i++) {
                    auto idx = bvh->sorted_prim[i];
                    if (overlap_point(pos, max_dist, bvh->pos[idx],
                            bvh->radius[idx], dist)) {
                        hit = true;
                        max_dist = dist;
                        eid = idx;
                        euv = {1, 0};
                    }
                }
            } break;
            case bvh_node_type::instance: {
                for (auto i = node.start; i < node.start + node.count; i++) {
                    auto& ist = bvh->instances[i];
                    if (overlap_bvh(ist.bvh,
                            transform_point(ist.frame_inv, pos), max_dist,
                            find_any, dist, iid, sid, eid, euv)) {
                        hit = true;
                        max_dist = dist;
                        iid = ist.iid;
                        sid = ist.sid;
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
    if (!intersect_bvh(bvh, ray, find_any, isec.dist, isec.iid, isec.sid,
            isec.eid, isec.euv))
        return {};
    return isec;
}

// Finds the closest element with a bvh (convenience wrapper).
intersection_point overlap_bvh(
    const bvh_tree* bvh, const vec3f& pos, float max_dist, bool find_any) {
    auto isec = intersection_point();
    if (!overlap_bvh(bvh, pos, max_dist, find_any, isec.dist, isec.iid,
            isec.sid, isec.eid, isec.euv))
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

// cleanup
material::~material() {
    if (ke_txt_info) delete ke_txt_info;
    if (kd_txt_info) delete kd_txt_info;
    if (ks_txt_info) delete ks_txt_info;
    if (kr_txt_info) delete kr_txt_info;
    if (kt_txt_info) delete kt_txt_info;
    if (rs_txt_info) delete rs_txt_info;
    if (bump_txt_info) delete bump_txt_info;
    if (disp_txt_info) delete disp_txt_info;
    if (norm_txt_info) delete norm_txt_info;
    if (occ_txt_info) delete occ_txt_info;
}

// cleanup
environment::environment() {
    if (ke_txt_info) delete ke_txt_info;
}

// cleanup
shape_group::~shape_group() {
    for (auto v : shapes) delete v;
}

// cleanup
animation_group::~animation_group() {
    for (auto v : animations) delete v;
}

// cleanup
scene::~scene() {
    for (auto v : shapes) delete v;
    for (auto v : instances) delete v;
    for (auto v : materials) delete v;
    for (auto v : textures) delete v;
    for (auto v : cameras) delete v;
    for (auto v : environments) delete v;
    for (auto v : nodes) delete v;
    for (auto v : animations) delete v;
}

// Shape value interpolated using barycentric coordinates
template <typename T>
T eval_elem(const shape* shp, const std::vector<T>& vals, int eid,
    const vec2f& euv, const T& def) {
    if (vals.empty()) return def;
    if (!shp->triangles.empty()) {
        return interpolate_triangle(vals, shp->triangles[eid], euv);
    } else if (!shp->lines.empty()) {
        return interpolate_line(vals, shp->lines[eid], euv.x);
    } else if (!shp->points.empty()) {
        return interpolate_point(vals, shp->points[eid]);
    } else if (!shp->quads.empty()) {
        return interpolate_quad(vals, shp->quads[eid], euv);
    } else {
        return vals[eid];  // points
    }
}

// Shape position interpolated using barycentric coordinates
vec3f eval_pos(const shape* shp, int eid, const vec2f& euv) {
    return eval_elem(shp, shp->pos, eid, euv, {0, 0, 0});
}
// Shape normal interpolated using barycentric coordinates
vec3f eval_norm(const shape* shp, int eid, const vec2f& euv) {
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

// Instance position interpolated using barycentric coordinates
vec3f eval_pos(const instance* ist, int sid, int eid, const vec2f& euv) {
    auto shp = ist->shp->shapes.at(sid);
    return transform_point(
        ist->frame, eval_elem(shp, shp->pos, eid, euv, {0, 0, 0}));
}
// Instance normal interpolated using barycentric coordinates
vec3f eval_norm(const instance* ist, int sid, int eid, const vec2f& euv) {
    auto shp = ist->shp->shapes.at(sid);
    return transform_direction(
        ist->frame, normalize(eval_elem(shp, shp->norm, eid, euv, {0, 0, 1})));
}

// Evaluate a texture
vec4f eval_texture(const texture* txt, const texture_info& info,
    const vec2f& texcoord, bool srgb, const vec4f& def) {
    if (!txt) return def;
    assert(!txt->hdr.empty() || !txt->ldr.empty());

    auto lookup = [&def, &txt, &srgb](int i, int j) {
        if (!txt->ldr.empty())
            return (srgb) ? srgb_to_linear(txt->ldr.at(i, j)) :
                            byte_to_float(txt->ldr.at(i, j));
        else if (!txt->hdr.empty())
            return txt->hdr.at(i, j);
        else
            return def;
    };

    // get image width/height
    auto w = (!txt->ldr.empty()) ? txt->ldr.width() : txt->hdr.width(),
         h = (!txt->ldr.empty()) ? txt->ldr.height() : txt->hdr.height();

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
    return ray3f(transform_point(cam->frame, o),
        transform_direction(cam->frame, normalize(q - o)));
}

// Generates a ray from a camera for pixel coordinates `ij`, the resolution
// `res`, the sub-pixel coordinates `puv` and the lens coordinates `luv` and the
// image resolution `res`.
ray3f eval_camera_ray(const camera* cam, const vec2i& ij, int res,
    const vec2f& puv, const vec2f& luv) {
    auto uv =
        vec2f{(ij.x + puv.x) / (cam->aspect * res), (1 - ij.y - puv.y) / res};
    return eval_camera_ray(cam, uv, luv);
}

// Subdivides shape elements.
void subdivide_shape_once(shape* shp, bool subdiv) {
    if (!shp->lines.empty()) {
        subdivide_lines(shp->lines, shp->pos, shp->norm, shp->texcoord,
            shp->color, shp->radius);
    } else if (!shp->triangles.empty()) {
        subdivide_triangles(shp->triangles, shp->pos, shp->norm, shp->texcoord,
            shp->color, shp->radius);
    } else if (!shp->quads.empty() && !subdiv) {
        subdivide_quads(shp->quads, shp->pos, shp->norm, shp->texcoord,
            shp->color, shp->radius);
    } else if (!shp->quads.empty() && subdiv) {
        subdivide_catmullclark(shp->quads, shp->pos, shp->norm, shp->texcoord,
            shp->color, shp->radius);
    } else if (!shp->quads_pos.empty() && !subdiv) {
        subdivide_quads(shp->quads_pos, shp->pos);
        subdivide_quads(shp->quads_norm, shp->norm);
        subdivide_quads(shp->quads_texcoord, shp->texcoord);
    } else if (!shp->quads_pos.empty() && subdiv) {
        subdivide_catmullclark(shp->quads_pos, shp->pos);
        subdivide_catmullclark(shp->quads_norm, shp->norm);
        subdivide_catmullclark(shp->quads_texcoord, shp->texcoord);
    } else if (!shp->beziers.empty()) {
        subdivide_beziers(shp->beziers, shp->pos, shp->norm, shp->texcoord,
            shp->color, shp->radius);
    }
}

// Compute shape normals. Supports only non-facevarying shapes.
void compute_normals(shape* shp) {
    if (!shp->points.empty()) {
        shp->norm.assign(shp->pos.size(), {0, 0, 1});
    } else if (!shp->lines.empty()) {
        compute_tangents(shp->lines, shp->pos, shp->norm);
    } else if (!shp->triangles.empty()) {
        compute_normals(shp->triangles, shp->pos, shp->norm);
    } else if (!shp->quads.empty()) {
        compute_normals(shp->quads, shp->pos, shp->norm);
    }
}

// Facet a shape. Supports only non-facevarying shapes
void facet_shape(shape* shp, bool recompute_normals) {
    if (!shp->lines.empty()) {
        facet_lines(shp->lines, shp->pos, shp->norm, shp->texcoord, shp->color,
            shp->radius);
        if (recompute_normals) compute_normals(shp);
    } else if (!shp->triangles.empty()) {
        facet_triangles(shp->triangles, shp->pos, shp->norm, shp->texcoord,
            shp->color, shp->radius);
        if (recompute_normals) compute_normals(shp);
    } else if (!shp->quads.empty()) {
        std::vector<int> verts;
        facet_quads(shp->quads, shp->pos, shp->norm, shp->texcoord, shp->color,
            shp->radius);
        if (recompute_normals) compute_normals(shp);
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
    if (facevarying_to_sharedvertex && !shp->quads_pos.empty()) {
        std::tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
            convert_face_varying(shp->quads_pos, shp->quads_norm,
                shp->quads_texcoord, shp->pos, shp->norm, shp->texcoord);
        shp->quads_pos = {};
        shp->quads_norm = {};
        shp->quads_texcoord = {};
    }
    if (quads_to_triangles && !shp->quads.empty()) {
        shp->triangles = convert_quads_to_triangles(shp->quads);
        shp->quads = {};
    }
    if (bezier_to_lines && !shp->beziers.empty()) {
        shp->lines = convert_bezier_to_lines(shp->beziers);
        shp->beziers = {};
    }
}

// Tesselate scene shapes and update pointers
void tesselate_shapes(scene* scn, bool subdivide,
    bool facevarying_to_sharedvertex, bool quads_to_triangles,
    bool bezier_to_lines) {
    for (auto sgr : scn->shapes) {
        for (auto shp : sgr->shapes) {
            tesselate_shape(shp, subdivide, facevarying_to_sharedvertex,
                quads_to_triangles, bezier_to_lines);
        }
    }
}

// Update animation transforms
void update_transforms(animation_group* agr, float time) {
    auto interpolate = [](keyframe_type type, const std::vector<float>& times,
                           const auto& vals, float time) {
        switch (type) {
            case keyframe_type::step:
                return eval_keyframed_step(times, vals, time);
            case keyframe_type::linear:
                return eval_keyframed_linear(times, vals, time);
            case keyframe_type::catmull_rom: return vals.at(0);
            case keyframe_type::bezier:
                return eval_keyframed_bezier(times, vals, time);
            default: throw std::runtime_error("should not have been here");
        }
        return vals.at(0);
    };

    for (auto anm : agr->animations) {
        if (!anm->translation.empty()) {
            auto val =
                interpolate(anm->type, anm->times, anm->translation, time);
            for (auto target : agr->targets)
                if (target.first == anm) target.second->translation = val;
        }
        if (!anm->rotation.empty()) {
            auto val = interpolate(anm->type, anm->times, anm->rotation, time);
            for (auto target : agr->targets)
                if (target.first == anm) target.second->rotation = val;
        }
        if (!anm->scaling.empty()) {
            auto val = interpolate(anm->type, anm->times, anm->scaling, time);
            for (auto target : agr->targets)
                if (target.first == anm) target.second->scaling = val;
        }
    }
}

// Update node transforms
void update_transforms(node* nde, const frame3f& parent = identity_frame3f) {
    auto frame = parent * nde->frame * translation_frame(nde->translation) *
                 rotation_frame(nde->rotation) * scaling_frame(nde->scaling);
    if (nde->ist) nde->ist->frame = frame;
    if (nde->cam) nde->cam->frame = frame;
    if (nde->env) nde->env->frame = frame;
    for (auto child : nde->children_) update_transforms(child, frame);
}

// Update node transforms
void update_transforms(scene* scn, float time) {
    for (auto agr : scn->animations) update_transforms(agr, time);
    for (auto nde : scn->nodes) nde->children_.clear();
    for (auto nde : scn->nodes)
        if (nde->parent) nde->parent->children_.push_back(nde);
    for (auto nde : scn->nodes)
        if (!nde->parent) update_transforms(nde);
}

// Compute animation range
vec2f compute_animation_range(const scene* scn) {
    if (scn->animations.empty()) return zero2f;
    auto range = vec2f{+flt_max, -flt_max};
    for (auto agr : scn->animations) {
        for (auto anm : agr->animations) {
            range.x = min(range.x, anm->times.front());
            range.y = max(range.y, anm->times.back());
        }
    }
    return range;
}

// Add missing values and elements
void add_elements(scene* scn, const add_elements_options& opts) {
    if (opts.smooth_normals) {
        for (auto sgr : scn->shapes) {
            for (auto shp : sgr->shapes) {
                if (!shp->norm.empty()) continue;
                shp->norm.resize(shp->pos.size(), {0, 0, 1});
                if (!shp->lines.empty() || !shp->triangles.empty() ||
                    !shp->quads.empty()) {
                    compute_normals(shp);
                }
                if (!shp->quads_pos.empty()) {
                    if (!shp->quads_norm.empty())
                        throw std::runtime_error("bad normals");
                    shp->quads_norm = shp->quads_pos;
                    compute_normals(shp->quads_pos, shp->pos, shp->norm);
                }
            }
        }
    }

    if (opts.tangent_space) {
        for (auto sgr : scn->shapes) {
            for (auto shp : sgr->shapes) {
                if (!shp->tangsp.empty()) continue;
                if (shp->texcoord.empty()) continue;
                if (!shp->mat) continue;
                if (!(shp->mat->norm_txt || shp->mat->bump_txt)) continue;
                if (!shp->triangles.empty()) {
                    compute_tangent_frames(shp->triangles, shp->pos, shp->norm,
                        shp->texcoord, shp->tangsp);
                } else if (!shp->quads.empty()) {
                    auto triangles = convert_quads_to_triangles(shp->quads);
                    compute_tangent_frames(triangles, shp->pos, shp->norm,
                        shp->texcoord, shp->tangsp);
                }
            }
        }
    }

    if (opts.texture_data) {
        for (auto txt : scn->textures) {
            if (txt->hdr.empty() && txt->ldr.empty()) {
                printf("unable to load texture %s\n", txt->path.c_str());
                txt->ldr = image4b(1, 1, {255, 255, 255, 255});
            }
        }
    }

    if (opts.shape_instances) {
        if (!scn->instances.empty()) return;
        for (auto shp : scn->shapes) {
            auto ist = new instance();
            ist->name = shp->name;
            ist->shp = shp;
            scn->instances.push_back(ist);
        }
    }

    if (opts.default_names || opts.default_paths) {
        auto cid = 0;
        for (auto cam : scn->cameras) {
            if (cam->name.empty())
                cam->name = "unnamed_camera_" + std::to_string(cid);
            cid++;
        }

        auto tid = 0;
        for (auto txt : scn->textures) {
            if (txt->name.empty())
                txt->name = "unnamed_texture_" + std::to_string(tid);
            tid++;
        }

        auto mid = 0;
        for (auto mat : scn->materials) {
            if (mat->name.empty())
                mat->name = "unnamed_material_" + std::to_string(mid);
            mid++;
        }

        auto sid = 0;
        for (auto shp : scn->shapes) {
            if (shp->name.empty())
                shp->name = "unnamed_shape_" + std::to_string(sid);
            sid++;
        }

        auto iid = 0;
        for (auto ist : scn->instances) {
            if (ist->name.empty())
                ist->name = "unnamed_instance_" + std::to_string(iid);
            iid++;
        }

        auto eid = 0;
        for (auto env : scn->environments) {
            if (env->name.empty())
                env->name = "unnamed_environment_" + std::to_string(eid);
            eid++;
        }
    }

    if (opts.default_paths) {
        for (auto txt : scn->textures) {
            if (txt->path != "") continue;
            txt->path = txt->name + ".png";
        }
        for (auto shp : scn->shapes) {
            if (shp->path != "") continue;
            shp->path = shp->name + ".bin";
        }
        for (auto anm : scn->animations) {
            if (anm->path != "") continue;
            anm->path = anm->name + ".bin";
        }
    }
}

// Make a view camera either copying a given one or building a
// default one.
camera* make_view_camera(const scene* scn, int camera_id) {
    if (scn->cameras.empty()) {
        auto bbox = compute_bounds(scn);
        auto bbox_center = (bbox.max + bbox.min) / 2.0f;
        auto bbox_size = bbox.max - bbox.min;
        auto bbox_msize = max(bbox_size[0], max(bbox_size[1], bbox_size[2]));
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
    merge_into->cameras.insert(merge_from->cameras.begin(),
        merge_from->cameras.end(), merge_into->cameras.end());
    merge_from->cameras.clear();
    merge_into->textures.insert(merge_from->textures.begin(),
        merge_from->textures.end(), merge_into->textures.end());
    merge_from->textures.clear();
    merge_into->materials.insert(merge_from->materials.begin(),
        merge_from->materials.end(), merge_into->materials.end());
    merge_from->materials.clear();
    merge_into->shapes.insert(merge_from->shapes.begin(),
        merge_from->shapes.end(), merge_into->shapes.end());
    merge_from->shapes.clear();
    merge_into->instances.insert(merge_from->instances.begin(),
        merge_from->instances.end(), merge_into->instances.end());
    merge_from->instances.clear();
    merge_into->environments.insert(merge_from->environments.begin(),
        merge_from->environments.end(), merge_into->environments.end());
    merge_from->environments.clear();
}

// Computes a shape bounding box (quick computation that ignores
// radius)
bbox3f compute_bounds(const shape* shp) {
    auto bbox = invalid_bbox3f;
    for (auto p : shp->pos) bbox += vec3f(p);
    return bbox;
}

// Updates the scene and scene's instances bounding boxes
bbox3f compute_bounds(const scene* scn) {
    auto shape_bboxes = std::unordered_map<shape_group*, bbox3f>();
    for (auto sgr : scn->shapes) {
        shape_bboxes[sgr] = invalid_bbox3f;
        for (auto shp : sgr->shapes) shape_bboxes[sgr] += compute_bounds(shp);
    }
    auto bbox = invalid_bbox3f;
    if (!scn->instances.empty()) {
        for (auto ist : scn->instances)
            bbox += transform_bbox(ist->frame, shape_bboxes.at(ist->shp));

    } else {
        for (auto shp : scn->shapes) bbox += shape_bboxes.at(shp);
    }
    return bbox;
}

// Flatten scene instances into separate meshes.
void flatten_instances(scene* scn) {
    if (scn->instances.empty()) return;
    auto shapes = scn->shapes;
    scn->shapes.clear();
    auto instances = scn->instances;
    scn->instances.clear();
    for (auto ist : instances) {
        if (!ist->shp) continue;
        auto nsgr = new shape_group();
        nsgr->name = ist->shp->name;
        nsgr->path = "";
        for (auto shp : ist->shp->shapes) {
            auto nshp = new shape(*shp);
            for (auto& p : nshp->pos) p = transform_point(ist->frame, p);
            for (auto& n : nshp->norm) n = transform_direction(ist->frame, n);
            nsgr->shapes.push_back(nshp);
        }
        scn->shapes.push_back(nsgr);
    }
    for (auto e : shapes) delete e;
    for (auto e : instances) delete e;
    for (auto e : scn->nodes) delete e;
    scn->nodes.clear();
    for (auto e : scn->animations) delete e;
    scn->animations.clear();
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
    auto smap = std::unordered_map<shape*, bvh_tree*>();
    for (auto sgr : scn->shapes) {
        for (auto shp : sgr->shapes) {
            shape_bvhs.push_back(make_bvh(shp, def_radius, equalsize));
            smap[shp] = shape_bvhs.back();
        }
    }

    // tree bvh
    auto bists = std::vector<bvh_instance>();
    for (auto iid = 0; iid < scn->instances.size(); iid++) {
        auto ist = scn->instances[iid];
        for (auto sid = 0; sid < ist->shp->shapes.size(); sid++) {
            auto bist = bvh_instance();
            bist.frame = ist->frame;
            bist.frame_inv = inverse(ist->frame);
            bist.iid = iid;
            bist.sid = sid;
            bist.bvh = smap.at(ist->shp->shapes.at(sid));
            bists.push_back(bist);
        }
    }
    return make_bvh(bists, shape_bvhs, true, equalsize);
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
        for (auto sgr : scn->shapes) {
            for (auto shp : sgr->shapes) {
                refit_bvh(get_shape_bvhs(bvh).at(sid), shp->pos, shp->radius,
                    def_radius);
                sid++;
            }
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

// Print scene info (call update bounds bes before)
void print_info(const scene* scn) {
    auto nverts = 0, nnorms = 0, ntexcoords = 0, npoints = 0, nlines = 0,
         ntriangles = 0, nquads = 0;
    for (auto sgr : scn->shapes) {
        for (auto shp : sgr->shapes) {
            nverts += shp->pos.size();
            nnorms += shp->norm.size();
            ntexcoords += shp->texcoord.size();
            npoints += shp->points.size();
            nlines += shp->lines.size();
            ntriangles += shp->triangles.size();
            nquads += shp->quads.size();
        }
    }

    auto bbox = compute_bounds(scn);
    auto bboxc = vec3f{(bbox.max[0] + bbox.min[0]) / 2,
        (bbox.max[1] + bbox.min[1]) / 2, (bbox.max[2] + bbox.min[2]) / 2};
    auto bboxs = vec3f{bbox.max[0] - bbox.min[0], bbox.max[1] - bbox.min[1],
        bbox.max[2] - bbox.min[2]};

    printf("number of cameras:      %d\n", (int)scn->cameras.size());
    printf("number of shapes:       %d\n", (int)scn->shapes.size());
    printf("number of instances:    %d\n", (int)scn->instances.size());
    printf("number of materials:    %d\n", (int)scn->materials.size());
    printf("number of textures:     %d\n", (int)scn->textures.size());
    printf("number of environments: %d\n", (int)scn->environments.size());
    printf("number of vertices:     %d\n", nverts);
    printf("number of normals:      %d\n", nnorms);
    printf("number of texcoords:    %d\n", ntexcoords);
    printf("number of points:       %d\n", npoints);
    printf("number of lines:        %d\n", nlines);
    printf("number of triangles:    %d\n", ntriangles);
    printf("number of quads:        %d\n", nquads);
    printf("\n");
    printf("bbox min:    %g %g %g\n", bbox.min[0], bbox.min[1], bbox.min[2]);
    printf("bbox max:    %g %g %g\n", bbox.max[0], bbox.max[1], bbox.max[2]);
    printf("bbox center: %g %g %g\n", bboxc[0], bboxc[1], bboxc[2]);
    printf("bbox size:   %g %g %g\n", bboxs[0], bboxs[1], bboxs[2]);
    printf("\n");
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
            txt->ldr = make_image(otxt->width, otxt->height, otxt->ncomp,
                otxt->datab.data(), vec4b{0, 0, 0, 255});
        } else if (!otxt->dataf.empty()) {
            txt->hdr = make_image(otxt->width, otxt->height, otxt->ncomp,
                otxt->dataf.data(), vec4f{0, 0, 0, 1});
        }
        scn->textures.push_back(txt);
        tmap[txt->path] = txt;
    }

    auto make_texture_info = [&tmap](const obj_texture_info& oinfo) {
        if (oinfo.path == "") return (texture_info*)nullptr;
        auto info = new texture_info();
        info->wrap_s = !oinfo.clamp;
        info->wrap_t = !oinfo.clamp;
        info->scale = oinfo.scale;
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
        mat->rs = (omat->ns >= 1e6f) ? 0 : pow(2 / (omat->ns + 2), 1 / 4.0f);
        mat->op = omat->op;
        mat->ke_txt = tmap.at(omat->ke_txt.path);
        mat->kd_txt = tmap.at(omat->kd_txt.path);
        mat->ks_txt = tmap.at(omat->ks_txt.path);
        mat->kr_txt = tmap.at(omat->kr_txt.path);
        mat->kt_txt = tmap.at(omat->kt_txt.path);
        mat->rs_txt = tmap.at(omat->ns_txt.path);
        mat->norm_txt = tmap.at(omat->norm_txt.path);
        mat->bump_txt = tmap.at(omat->bump_txt.path);
        mat->disp_txt = tmap.at(omat->disp_txt.path);
        mat->ke_txt_info = make_texture_info(omat->ke_txt);
        mat->kd_txt_info = make_texture_info(omat->kd_txt);
        mat->ks_txt_info = make_texture_info(omat->ks_txt);
        mat->kr_txt_info = make_texture_info(omat->kr_txt);
        mat->kt_txt_info = make_texture_info(omat->kt_txt);
        mat->rs_txt_info = make_texture_info(omat->ns_txt);
        mat->norm_txt_info = make_texture_info(omat->norm_txt);
        mat->bump_txt_info = make_texture_info(omat->bump_txt);
        mat->disp_txt_info = make_texture_info(omat->disp_txt);
        switch (omat->illum) {
            case 0:  // Color on and Ambient off
            case 1:  // Color on and Ambient on
            case 2:  // Highlight on
            case 3:  // Reflection on and Ray trace on
                mat->op = 1;
                mat->kt = {0, 0, 0};
                break;
            case 4:  // Transparency: Glass on
                // Reflection: Ray trace on
                break;
            case 5:  // Reflection: Fresnel on and Ray trace on
                mat->op = 1;
                mat->kt = {0, 0, 0};
                break;
            case 6:  // Transparency: Refraction on
                     // Reflection: Fresnel off and Ray trace on
            case 7:  // Transparency: Refraction on
                // Reflection: Fresnel on and Ray trace on
                break;
            case 8:  // Reflection on and Ray trace off
                mat->op = 1;
                mat->kt = {0, 0, 0};
                break;
            case 9:  // Transparency: Glass on
                // Reflection: Ray trace off
                break;
        }
        scn->materials.push_back(mat);
        mmap[mat->name] = mat;
    }

    // convert meshes
    auto omap = std::unordered_map<std::string, shape_group*>{{"", nullptr}};
    for (auto omsh : obj->objects) {
        auto sgr = new shape_group();
        sgr->name = omsh->name;
        for (auto oshp : omsh->groups) {
            if (oshp->verts.empty()) continue;
            if (oshp->elems.empty()) continue;

            auto shp = new shape();
            shp->name = oshp->groupname;
            shp->mat = mmap[oshp->matname];
            if (oshp->props.find("subdivision") != oshp->props.end()) {
                shp->subdivision =
                    atoi(oshp->props.at("subdivision").at(0).c_str());
            }
            if (oshp->props.find("catmullclark") != oshp->props.end()) {
                shp->catmullclark =
                    (bool)atoi(oshp->props.at("catmullclark").at(1).c_str());
            }

            // check to see if this shuold be face-varying or flat
            // quads
            auto as_facevarying = false, as_quads = false;
            if (opts.preserve_quads || opts.preserve_facevarying) {
                auto face_max = 0;
                for (auto& elem : oshp->elems) {
                    if (elem.type != obj_element_type::face) {
                        face_max = 0;
                        break;
                    } else {
                        face_max = max(face_max, (int)elem.size);
                    }
                }
                as_quads = opts.preserve_quads && face_max > 3;
                as_facevarying = opts.preserve_facevarying && face_max > 2;
                // in case of facevarying, check if there is really
                // need for it
                if (as_facevarying) {
                    auto need_facevarying = false;
                    for (auto& elem : oshp->elems) {
                        for (auto i = elem.start; i < elem.start + elem.size;
                             i++) {
                            auto& v = oshp->verts[i];
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
                for (auto& vert : oshp->verts) {
                    if (vert_map.find(vert) == vert_map.end()) {
                        auto s = (int)vert_map.size();
                        vert_map[vert] = s;
                    }
                    vert_ids.push_back(vert_map.at(vert));
                }

                // convert elements
                for (auto& elem : oshp->elems) {
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
                auto v = oshp->verts[0];
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

                // fix smoothing
                if (!oshp->smoothing && opts.obj_facet_non_smooth) {
                    auto faceted = new shape();
                    faceted->name = shp->name;
                    auto pidx = std::vector<int>();
                    for (auto point : shp->points) {
                        faceted->points.push_back((int)pidx.size());
                        pidx.push_back(point);
                    }
                    for (auto line : shp->lines) {
                        faceted->lines.push_back(
                            {(int)pidx.size() + 0, (int)pidx.size() + 1});
                        pidx.push_back(line.x);
                        pidx.push_back(line.y);
                    }
                    for (auto triangle : shp->triangles) {
                        faceted->triangles.push_back({(int)pidx.size() + 0,
                            (int)pidx.size() + 1, (int)pidx.size() + 2});
                        pidx.push_back(triangle.x);
                        pidx.push_back(triangle.y);
                        pidx.push_back(triangle.z);
                    }
                    for (auto idx : pidx) {
                        if (!shp->pos.empty())
                            faceted->pos.push_back(shp->pos[idx]);
                        if (!shp->norm.empty())
                            faceted->norm.push_back(shp->norm[idx]);
                        if (!shp->texcoord.empty())
                            faceted->texcoord.push_back(shp->texcoord[idx]);
                        if (!shp->color.empty())
                            faceted->color.push_back(shp->color[idx]);
                        if (!shp->radius.empty())
                            faceted->radius.push_back(shp->radius[idx]);
                    }
                    delete shp;
                    shp = faceted;
                }
            } else {
                // insert all vertices
                std::unordered_map<int, int> pos_map, norm_map, texcoord_map;
                std::vector<int> pos_ids, norm_ids, texcoord_ids;
                for (auto& vert : oshp->verts) {
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
                for (auto elem : oshp->elems) {
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

                // fix smoothing
                if (!oshp->smoothing && opts.obj_facet_non_smooth) {}
            }
            sgr->shapes.push_back(shp);
        }
        scn->shapes.push_back(sgr);
        omap[omsh->name] = sgr;
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
        for (auto mat : scn->materials) {
            if (mat->name == oenv->matname) {
                env->ke = mat->ke;
                env->ke_txt = mat->ke_txt;
                env_mat.insert(mat);
            }
        }
        env->frame = oenv->frame;
        scn->environments.push_back(env);
        emap[env->name] = env;
    }

    // remove env materials
    for (auto sgr : scn->shapes)
        for (auto shp : sgr->shapes) env_mat.erase(shp->mat);
    for (auto mat : env_mat) {
        auto end =
            std::remove(scn->materials.begin(), scn->materials.end(), mat);
        scn->materials.erase(end, scn->materials.end());
        delete mat;
    }

    if (!obj->nodes.empty()) {
        // convert nodes
        for (auto onde : obj->nodes) {
            auto nde = new node();
            nde->name = onde->name;
            nde->cam = cmap.at(onde->camname);
            if (!onde->objname.empty()) {
                auto ist = new instance();
                ist->name = onde->name;
                ist->shp = omap.at(onde->objname);
                nde->ist = ist;
                scn->instances.push_back(ist);
            }
            nde->env = emap.at(onde->envname);
            nde->translation = onde->translation;
            nde->rotation = onde->rotation;
            nde->scaling = onde->scaling;
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
    }

    // update transforms
    update_transforms(scn);

    // remove hierarchy if necessary
    if (!opts.preserve_hierarchy) {
        for (auto nde : scn->nodes) delete nde;
        scn->nodes.clear();
        for (auto anm : scn->animations) delete anm;
        scn->animations.clear();
    }

    // done
    return scn;
}

// Load an obj scene
scene* load_obj_scene(const std::string& filename, const load_options& opts) {
    auto oscn =
        std::unique_ptr<obj_scene>(load_obj(filename, opts.load_textures,
            opts.skip_missing, opts.obj_flip_texcoord, opts.obj_flip_tr));
    auto scn = std::unique_ptr<scene>(obj_to_scene(oscn.get(), opts));
    return scn.release();
}

// Save an scene
obj_scene* scene_to_obj(const scene* scn) {
    auto obj = new obj_scene();

    auto make_texture_info = [](const texture* txt, const texture_info* info,
                                 bool bump = false) {
        auto oinfo = obj_texture_info();
        if (!txt) return oinfo;
        oinfo.path = txt->path;
        if (!info) return oinfo;
        oinfo.clamp = !info->wrap_s && !info->wrap_t;
        if (bump) oinfo.scale = info->scale;
        return oinfo;
    };

    // convert textures
    for (auto txt : scn->textures) {
        auto otxt = new obj_texture();
        otxt->path = txt->path;
        if (!txt->hdr.empty()) {
            otxt->width = txt->hdr.width();
            otxt->height = txt->hdr.height();
            otxt->ncomp = 4;
            otxt->dataf.assign((float*)data(txt->hdr),
                (float*)data(txt->hdr) +
                    txt->hdr.width() * txt->hdr.height() * 4);
        }
        if (!txt->ldr.empty()) {
            otxt->width = txt->ldr.width();
            otxt->height = txt->ldr.height();
            otxt->ncomp = 4;
            otxt->datab.assign((uint8_t*)data(txt->ldr),
                (uint8_t*)data(txt->ldr) +
                    txt->ldr.width() * txt->ldr.height() * 4);
        }
        obj->textures.push_back(otxt);
    }

    // convert materials
    for (auto mat : scn->materials) {
        auto omat = new obj_material();
        omat->name = mat->name;
        omat->ke = {mat->ke.x, mat->ke.y, mat->ke.z};
        omat->ke_txt = make_texture_info(mat->ke_txt, mat->ke_txt_info);
        switch (mat->type) {
            case material_type::specular_roughness: {
                omat->kd = {mat->kd.x, mat->kd.y, mat->kd.z};
                omat->ks = {mat->ks.x, mat->ks.y, mat->ks.z};
                omat->kr = {mat->kr.x, mat->kr.y, mat->kr.z};
                omat->kt = {mat->kt.x, mat->kt.y, mat->kt.z};
                omat->ns = (mat->rs) ? 2 / pow(mat->rs, 4.0f) - 2 : 1e6;
                omat->op = mat->op;
                omat->kd_txt = make_texture_info(mat->kd_txt, mat->kd_txt_info);
                omat->ks_txt = make_texture_info(mat->ks_txt, mat->ks_txt_info);
                omat->kr_txt = make_texture_info(mat->kr_txt, mat->kr_txt_info);
                omat->kt_txt = make_texture_info(mat->kt_txt, mat->kt_txt_info);
            } break;
            case material_type::metallic_roughness: {
                if (mat->rs == 1 && mat->ks.x == 0) {
                    omat->kd = mat->kd;
                    omat->ks = {0, 0, 0};
                    omat->ns = 1;
                } else {
                    auto kd = mat->kd * (1 - 0.04f) * (1 - mat->ks.x);
                    auto ks = mat->kd * mat->ks.x +
                              vec3f{0.04f, 0.04f, 0.04f} * (1 - mat->ks.x);
                    omat->kd = {kd.x, kd.y, kd.z};
                    omat->ks = {ks.x, ks.y, ks.z};
                    omat->ns = (mat->rs) ? 2 / pow(mat->rs, 4.0f) - 2 : 1e6;
                }
                omat->op = mat->op;
                if (mat->ks.x < 0.5f) {
                    omat->kd_txt =
                        make_texture_info(mat->kd_txt, mat->kd_txt_info);
                } else {
                    omat->ks_txt =
                        make_texture_info(mat->ks_txt, mat->ks_txt_info);
                }
            } break;
            case material_type::specular_glossiness: {
                omat->kd = {mat->kd.x, mat->kd.y, mat->kd.z};
                omat->ks = {mat->ks.x, mat->ks.y, mat->ks.z};
                omat->ns = (mat->rs) ? 2 / pow(1 - mat->rs, 4.0f) - 2 : 1e6;
                omat->op = mat->op;
                omat->kd_txt = make_texture_info(mat->kd_txt, mat->kd_txt_info);
                omat->ks_txt = make_texture_info(mat->ks_txt, mat->ks_txt_info);
            } break;
        }
        omat->bump_txt =
            make_texture_info(mat->bump_txt, mat->bump_txt_info, true);
        omat->disp_txt =
            make_texture_info(mat->disp_txt, mat->disp_txt_info, true);
        omat->norm_txt =
            make_texture_info(mat->norm_txt, mat->norm_txt_info, true);
        if (mat->op < 1 || mat->kt != zero3f) {
            omat->illum = 4;
        } else {
            omat->illum = 2;
        }
        obj->materials.push_back(omat);
    }

    // convert shapes
    for (auto sgr : scn->shapes) {
        auto oobj = new obj_object();
        oobj->name = sgr->name;
        for (auto shp : sgr->shapes) {
            auto offset = obj_vertex{(int)obj->pos.size(),
                (int)obj->texcoord.size(), (int)obj->norm.size(),
                (int)obj->color.size(), (int)obj->radius.size()};
            for (auto& v : shp->pos) obj->pos.push_back({v.x, v.y, v.z});
            for (auto& v : shp->norm) obj->norm.push_back({v.x, v.y, v.z});
            for (auto& v : shp->texcoord) obj->texcoord.push_back({v.x, v.y});
            for (auto& v : shp->color)
                obj->color.push_back({v.x, v.y, v.z, v.w});
            for (auto& v : shp->radius) obj->radius.push_back(v);
            auto group = new obj_group();
            group->groupname = shp->name;
            group->matname = (shp->mat) ? shp->mat->name : "";
            if (shp->subdivision)
                group->props["subdivision"].push_back(
                    std::to_string(shp->subdivision));
            if (shp->catmullclark) group->props["catmullclark"].push_back("1");
            for (auto point : shp->points) {
                group->elems.push_back({(uint32_t)group->verts.size(),
                    obj_element_type::point, 1});
                auto vert = obj_vertex{-1, -1, -1, -1, -1};
                if (!shp->pos.empty()) vert.pos = offset.pos + point;
                if (!shp->texcoord.empty())
                    vert.texcoord = offset.texcoord + point;
                if (!shp->norm.empty()) vert.norm = offset.norm + point;
                if (!shp->color.empty()) vert.color = offset.color + point;
                if (!shp->radius.empty()) vert.radius = offset.radius + point;
                group->verts.push_back(vert);
            }
            for (auto line : shp->lines) {
                group->elems.push_back(
                    {(uint32_t)group->verts.size(), obj_element_type::line, 2});
                for (auto vid : line) {
                    auto vert = obj_vertex{-1, -1, -1, -1, -1};
                    if (!shp->pos.empty()) vert.pos = offset.pos + vid;
                    if (!shp->texcoord.empty())
                        vert.texcoord = offset.texcoord + vid;
                    if (!shp->norm.empty()) vert.norm = offset.norm + vid;
                    if (!shp->color.empty()) vert.color = offset.color + vid;
                    if (!shp->radius.empty()) vert.radius = offset.radius + vid;
                    group->verts.push_back(vert);
                }
            }
            for (auto triangle : shp->triangles) {
                group->elems.push_back(
                    {(uint32_t)group->verts.size(), obj_element_type::face, 3});
                for (auto vid : triangle) {
                    auto vert = obj_vertex{-1, -1, -1, -1, -1};
                    if (!shp->pos.empty()) vert.pos = offset.pos + vid;
                    if (!shp->texcoord.empty())
                        vert.texcoord = offset.texcoord + vid;
                    if (!shp->norm.empty()) vert.norm = offset.norm + vid;
                    if (!shp->color.empty()) vert.color = offset.color + vid;
                    if (!shp->radius.empty()) vert.radius = offset.radius + vid;
                    group->verts.push_back(vert);
                }
            }
            for (auto quad : shp->quads) {
                group->elems.push_back(
                    {(uint32_t)group->verts.size(), obj_element_type::face,
                        (uint16_t)((quad.z == quad.w) ? 3 : 4)});
                if (group->elems.back().size == 3) {
                    for (auto vid : {quad.x, quad.y, quad.z}) {
                        auto vert = obj_vertex{-1, -1, -1, -1, -1};
                        if (!shp->pos.empty()) vert.pos = offset.pos + vid;
                        if (!shp->texcoord.empty())
                            vert.texcoord = offset.texcoord + vid;
                        if (!shp->norm.empty()) vert.norm = offset.norm + vid;
                        if (!shp->color.empty())
                            vert.color = offset.color + vid;
                        if (!shp->radius.empty())
                            vert.radius = offset.radius + vid;
                        group->verts.push_back(vert);
                    }
                } else {
                    for (auto vid : quad) {
                        auto vert = obj_vertex{-1, -1, -1, -1, -1};
                        if (!shp->pos.empty()) vert.pos = offset.pos + vid;
                        if (!shp->texcoord.empty())
                            vert.texcoord = offset.texcoord + vid;
                        if (!shp->norm.empty()) vert.norm = offset.norm + vid;
                        if (!shp->color.empty())
                            vert.color = offset.color + vid;
                        if (!shp->radius.empty())
                            vert.radius = offset.radius + vid;
                        group->verts.push_back(vert);
                    }
                }
            }
            for (auto fid = 0; fid < shp->quads_pos.size(); fid++) {
                group->elems.push_back(
                    {(uint32_t)group->verts.size(), obj_element_type::face, 4});
                auto last_vid = -1;
                for (auto i = 0; i < 4; i++) {
                    if (last_vid == shp->quads_pos[fid][i]) continue;
                    auto vert = obj_vertex{-1, -1, -1, -1, -1};
                    if (!shp->pos.empty() && !shp->quads_pos.empty())
                        vert.pos = offset.pos + shp->quads_pos[fid][i];
                    if (!shp->texcoord.empty() && !shp->quads_texcoord.empty())
                        vert.texcoord =
                            offset.texcoord + shp->quads_texcoord[fid][i];
                    if (!shp->norm.empty() && !shp->quads_norm.empty())
                        vert.norm = offset.norm + shp->quads_norm[fid][i];
                    group->verts.push_back(vert);
                    last_vid = shp->quads_pos[fid][i];
                }
            }
            for (auto bezier : shp->beziers) {
                group->elems.push_back({(uint32_t)group->verts.size(),
                    obj_element_type::bezier, 4});
                for (auto vid : bezier) {
                    auto vert = obj_vertex{-1, -1, -1, -1, -1};
                    if (!shp->pos.empty()) vert.pos = offset.pos + vid;
                    if (!shp->texcoord.empty())
                        vert.texcoord = offset.texcoord + vid;
                    if (!shp->norm.empty()) vert.norm = offset.norm + vid;
                    if (!shp->color.empty()) vert.color = offset.color + vid;
                    if (!shp->radius.empty()) vert.radius = offset.radius + vid;
                    group->verts.push_back(vert);
                }
            }
            oobj->groups.push_back(group);
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
        auto omat = new obj_material();
        omat->name = env->name + "_mat";
        omat->ke = env->ke;
        omat->ke_txt = make_texture_info(env->ke_txt, env->ke_txt_info);
        oenv->name = env->name;
        oenv->matname = omat->name;
        oenv->frame = env->frame;
        obj->materials.push_back(omat);
        obj->environments.push_back(oenv);
    }

    // convert hierarchy
    if (!obj->nodes.empty()) {
        for (auto nde : scn->nodes) {
            auto onde = new obj_node();
            onde->name = nde->name;
            if (nde->cam) onde->camname = nde->cam->name;
            if (nde->ist) onde->objname = nde->ist->name;
            if (nde->env) onde->envname = nde->env->name;
            onde->frame = nde->frame;
            onde->translation = nde->translation;
            onde->rotation = nde->rotation;
            onde->scaling = nde->scaling;
            obj->nodes.push_back(onde);
        }

        // parent
        for (auto idx = 0; idx < scn->nodes.size(); idx++) {
            auto nde = scn->nodes.at(idx);
            if (!nde->parent) continue;
            auto onde = obj->nodes.at(idx);
            onde->parent = nde->parent->name;
        }
    } else {
        for (auto ist : scn->instances) {
            auto onde = new obj_node();
            onde->name = ist->name;
            onde->objname = (ist->shp) ? ist->shp->name : "<undefined>";
            onde->frame = ist->frame;
            obj->nodes.push_back(onde);
        }
    }

    return obj;
}

// Save an obj scene
void save_obj_scene(
    const std::string& filename, const scene* scn, const save_options& opts) {
    auto oscn = std::unique_ptr<obj_scene>(scene_to_obj(scn));
    save_obj(filename, oscn.get(), opts.save_textures, opts.skip_missing,
        opts.obj_flip_texcoord, opts.obj_flip_tr);
}

#if YGL_GLTF

// Flattens a gltf file into a flattened asset.
scene* gltf_to_scene(const glTF* gltf, const load_options& opts) {
    // clear asset
    auto scn = new scene();

    // convert images
    for (auto gtxt : gltf->images) {
        auto txt = new texture();
        txt->name = gtxt->name;
        txt->path = (startswith(gtxt->uri, "data:")) ? std::string("inlines") :
                                                       gtxt->uri;
        if (!gtxt->data.datab.empty()) {
            txt->ldr = make_image(gtxt->data.width, gtxt->data.height,
                gtxt->data.ncomp, gtxt->data.datab.data(), vec4b{0, 0, 0, 255});
        } else if (!gtxt->data.dataf.empty()) {
            txt->hdr = make_image(gtxt->data.width, gtxt->data.height,
                gtxt->data.ncomp, gtxt->data.dataf.data(), vec4f{0, 0, 0, 1});
        }
        scn->textures.push_back(txt);
    }

    // add a texture
    auto get_texture = [gltf, scn](glTFTextureInfo* ginfo) -> texture* {
        if (!ginfo) return nullptr;
        auto gtxt = gltf->get(ginfo->index);
        if (!gtxt || !gtxt->source) return nullptr;
        return scn->textures.at((int)gtxt->source);
    };

    // add a texture
    auto make_texture_info = [gltf, scn](glTFTextureInfo* ginfo,
                                 bool normal = false, bool occlusion = false) {
        auto info = (texture_info*)nullptr;
        if (!ginfo) return info;
        auto gtxt = gltf->get(ginfo->index);
        if (!gtxt || !gtxt->source) return info;
        auto txt = scn->textures.at((int)gtxt->source);
        if (!txt) return info;
        auto gsmp = gltf->get(gtxt->sampler);
        info = new texture_info();
        if (gsmp) {
            info->linear = gsmp->magFilter != glTFSamplerMagFilter::Nearest;
            info->mipmap = gsmp->minFilter != glTFSamplerMinFilter::Linear &&
                           gsmp->minFilter != glTFSamplerMinFilter::Nearest;
            info->wrap_s = gsmp->wrapS != glTFSamplerWrapS::ClampToEdge;
            info->wrap_t = gsmp->wrapT != glTFSamplerWrapT::ClampToEdge;
        }
        if (normal) {
            auto ninfo = (glTFMaterialNormalTextureInfo*)ginfo;
            info->scale = ninfo->scale;
        }
        if (occlusion) {
            auto ninfo = (glTFMaterialOcclusionTextureInfo*)ginfo;
            info->scale = ninfo->strength;
        }
        return info;
    };

    // convert materials
    for (auto gmat : gltf->materials) {
        auto mat = new material();
        mat->name = gmat->name;
        mat->ke = gmat->emissiveFactor;
        mat->ke_txt = get_texture(gmat->emissiveTexture);
        mat->ke_txt_info = make_texture_info(gmat->emissiveTexture);
        if (gmat->pbrMetallicRoughness) {
            mat->type = material_type::metallic_roughness;
            auto gmr = gmat->pbrMetallicRoughness;
            mat->kd = {gmr->baseColorFactor[0], gmr->baseColorFactor[1],
                gmr->baseColorFactor[2]};
            mat->op = gmr->baseColorFactor[3];
            mat->ks = {
                gmr->metallicFactor, gmr->metallicFactor, gmr->metallicFactor};
            mat->rs = gmr->roughnessFactor;
            mat->kd_txt = get_texture(gmr->baseColorTexture);
            mat->kd_txt_info = make_texture_info(gmr->baseColorTexture);
            mat->ks_txt = get_texture(gmr->metallicRoughnessTexture);
            mat->ks_txt_info = make_texture_info(gmr->metallicRoughnessTexture);
        }
        if (gmat->pbrSpecularGlossiness) {
            mat->type = material_type::specular_glossiness;
            auto gsg = gmat->pbrSpecularGlossiness;
            mat->kd = {gsg->diffuseFactor[0], gsg->diffuseFactor[1],
                gsg->diffuseFactor[2]};
            mat->op = gsg->diffuseFactor[3];
            mat->ks = gsg->specularFactor;
            mat->rs = gsg->glossinessFactor;
            mat->kd_txt = get_texture(gsg->diffuseTexture);
            mat->kd_txt_info = make_texture_info(gsg->diffuseTexture);
            mat->ks_txt = get_texture(gsg->specularGlossinessTexture);
            mat->ks_txt_info =
                make_texture_info(gsg->specularGlossinessTexture);
        }
        mat->norm_txt = get_texture(gmat->normalTexture);
        mat->norm_txt_info =
            make_texture_info(gmat->normalTexture, true, false);
        mat->occ_txt = get_texture(gmat->occlusionTexture);
        mat->occ_txt_info =
            make_texture_info(gmat->occlusionTexture, false, true);
        mat->double_sided = gmat->doubleSided;
        scn->materials.push_back(mat);
    }

    // convert meshes
    for (auto gmesh : gltf->meshes) {
        auto sgr = new shape_group();
        sgr->name = gmesh->name;
        // primitives
        auto gid = 0;
        for (auto gprim : gmesh->primitives) {
            auto shp = new shape();
            shp->name = gmesh->name + "_" + std::to_string(gid++);
            if (gprim->material) {
                shp->mat = scn->materials[(int)gprim->material];
            }
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
            sgr->shapes.push_back(shp);
        }
        scn->shapes.push_back(sgr);
    }

    // convert cameras
    auto cameras = std::vector<camera>();
    for (auto gcam : gltf->cameras) {
        cameras.push_back({});
        auto cam = &cameras.back();
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
    }

    // convert nodes
    for (auto gnde : gltf->nodes) {
        auto nde = new node();
        nde->name = gnde->name;
        if (gnde->camera) {
            auto cam = new camera(cameras[(int)gnde->camera]);
            nde->cam = cam;
            scn->cameras.push_back(cam);
        }
        if (gnde->mesh) {
            auto ist = new instance();
            ist->name = gnde->name;
            ist->shp = scn->shapes[(int)gnde->mesh];
            nde->ist = ist;
            scn->instances.push_back(ist);
#if 0
            for (auto shp : node->msh->shapes) {
                if (node->morph_weights.size() < shp->morph_targets.size()) {
                    node->morph_weights.resize(shp->morph_targets.size());
                }
            }
#endif
        }
        nde->translation = gnde->translation;
        nde->rotation = gnde->rotation;
        nde->scaling = gnde->scale;
        nde->frame = mat_to_frame(gnde->matrix);
#if 0
        nde->weights = gnde->weights;
#endif
        scn->nodes.push_back(nde);
    }

    // set up parent pointers
    for (auto nid = 0; nid < gltf->nodes.size(); nid++) {
        auto gnde = gltf->nodes[nid];
        auto nde = scn->nodes[nid];
        for (auto cid : gnde->children) scn->nodes[(int)cid]->parent = nde;
    }

    // keyframe type conversion
    static auto keyframe_types =
        std::unordered_map<glTFAnimationSamplerInterpolation, keyframe_type>{
            {glTFAnimationSamplerInterpolation::NotSet, keyframe_type::linear},
            {glTFAnimationSamplerInterpolation::Linear, keyframe_type::linear},
            {glTFAnimationSamplerInterpolation::Step, keyframe_type::step},
            {glTFAnimationSamplerInterpolation::CubicSpline,
                keyframe_type::bezier},
            {glTFAnimationSamplerInterpolation::CatmullRomSpline,
                keyframe_type::catmull_rom},
        };

    // convert animations
    for (auto ganm : gltf->animations) {
        auto agr = new animation_group();
        agr->name = ganm->name;
        auto sampler_map = std::unordered_map<vec2i, animation*>();
        for (auto gchannel : ganm->channels) {
            if (sampler_map.find({(int)gchannel->sampler,
                    (int)gchannel->target->path}) == sampler_map.end()) {
                auto gsampler = ganm->get(gchannel->sampler);
                auto anm = new animation();
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
                            anm->rotation.push_back(
                                (quat4f)output_view.getv4f(i));
                    } break;
                    case glTFAnimationChannelTargetPath::Scale: {
                        anm->scaling.reserve(output_view.size());
                        for (auto i = 0; i < output_view.size(); i++)
                            anm->scaling.push_back(output_view.getv3f(i));
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
                sampler_map[{
                    (int)gchannel->sampler, (int)gchannel->target->path}] = anm;
                agr->animations.push_back(anm);
            }
            agr->targets.push_back({sampler_map.at({(int)gchannel->sampler,
                                        (int)gchannel->target->path}),
                scn->nodes[(int)gchannel->target->node]});
        }
        scn->animations.push_back(agr);
    }

    // compute transforms
    update_transforms(scn, 0);

    // remove hierarchy if necessary
    if (!opts.preserve_hierarchy) {
        for (auto nde : scn->nodes) delete nde;
        scn->nodes.clear();
        for (auto anm : scn->animations) delete anm;
        scn->animations.clear();
    }

    return scn;
}

// Load an gltf scene
scene* load_gltf_scene(const std::string& filename, const load_options& opts) {
    auto gscn = std::unique_ptr<glTF>(
        load_gltf(filename, true, opts.load_textures, opts.skip_missing));
    auto scn = std::unique_ptr<scene>(gltf_to_scene(gscn.get(), opts));
    if (!scn) {
        throw std::runtime_error("could not convert gltf scene");
        return nullptr;
    }
    return scn.release();
}

// Unflattnes gltf
glTF* scene_to_gltf(
    const scene* scn, const std::string& buffer_uri, bool separate_buffers) {
    auto gltf = std::unique_ptr<glTF>(new glTF());

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
        if (!txt->hdr.empty()) {
            gimg->data.width = txt->hdr.width();
            gimg->data.height = txt->hdr.height();
            gimg->data.ncomp = 4;
            gimg->data.dataf.assign((float*)data(txt->hdr),
                (float*)data(txt->hdr) +
                    txt->hdr.width() * txt->hdr.height() * 4);
        }
        if (!txt->ldr.empty()) {
            gimg->data.width = txt->ldr.width();
            gimg->data.height = txt->ldr.height();
            gimg->data.ncomp = 4;
            gimg->data.datab.assign((uint8_t*)data(txt->ldr),
                (uint8_t*)data(txt->ldr) +
                    txt->ldr.width() * txt->ldr.height() * 4);
        }
        gltf->images.push_back(gimg);
    }

    // index of an object
    auto index = [](const auto& vec, auto* val) -> int {
        auto pos = find(vec.begin(), vec.end(), val);
        if (pos == vec.end()) return -1;
        return (int)(pos - vec.begin());
    };

    // add a texture and sampler
    auto add_texture_info = [&gltf, &index, scn](const texture* txt,
                                const texture_info* info, bool norm = false,
                                bool occ = false) {
        if (!txt) return (glTFTextureInfo*)nullptr;
        auto gtxt = new glTFTexture();
        gtxt->name = txt->name;
        gtxt->source = glTFid<glTFImage>(index(scn->textures, txt));

        // check if it is default
        auto is_default = !info || (info->wrap_s && info->wrap_t &&
                                       info->linear && info->mipmap);
        if (!is_default) {
            auto gsmp = new glTFSampler();
            gsmp->wrapS = (!info || info->wrap_s) ?
                              glTFSamplerWrapS::Repeat :
                              glTFSamplerWrapS::ClampToEdge;
            gsmp->wrapT = (!info || info->wrap_t) ?
                              glTFSamplerWrapT::Repeat :
                              glTFSamplerWrapT::ClampToEdge;
            gsmp->minFilter = (!info || info->mipmap) ?
                                  glTFSamplerMinFilter::LinearMipmapLinear :
                                  glTFSamplerMinFilter::Nearest;
            gsmp->magFilter = (!info || info->linear) ?
                                  glTFSamplerMagFilter::Linear :
                                  glTFSamplerMagFilter::Nearest;
            gtxt->sampler = glTFid<glTFSampler>((int)gltf->samplers.size());
            gltf->samplers.push_back(gsmp);
        }
        gltf->textures.push_back(gtxt);
        if (norm) {
            auto ginfo = new glTFMaterialNormalTextureInfo();
            ginfo->index = glTFid<glTFTexture>{(int)gltf->textures.size() - 1};
            ginfo->scale = (info) ? info->scale : 1;
            return (glTFTextureInfo*)ginfo;
        } else if (occ) {
            auto ginfo = new glTFMaterialOcclusionTextureInfo();
            ginfo->index = glTFid<glTFTexture>{(int)gltf->textures.size() - 1};
            ginfo->strength = (info) ? info->scale : 1;
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
        gmat->emissiveTexture = add_texture_info(mat->ke_txt, mat->ke_txt_info);
        switch (mat->type) {
            case material_type::specular_roughness: {
                gmat->pbrSpecularGlossiness =
                    new glTFMaterialPbrSpecularGlossiness();
                auto gsg = gmat->pbrSpecularGlossiness;
                gsg->diffuseFactor = {
                    mat->kd[0], mat->kd[1], mat->kd[2], mat->op};
                gsg->specularFactor = mat->ks;
                gsg->glossinessFactor = 1 - mat->rs;
                gsg->diffuseTexture =
                    add_texture_info(mat->kd_txt, mat->kd_txt_info);
                gsg->specularGlossinessTexture =
                    add_texture_info(mat->ks_txt, mat->ks_txt_info);
            } break;
            case material_type::metallic_roughness: {
                gmat->pbrMetallicRoughness =
                    new glTFMaterialPbrMetallicRoughness();
                auto gmr = gmat->pbrMetallicRoughness;
                gmr->baseColorFactor = {
                    mat->kd.x, mat->kd.y, mat->kd.z, mat->op};
                gmr->metallicFactor = mat->ks.x;
                gmr->roughnessFactor = mat->rs;
                gmr->baseColorTexture =
                    add_texture_info(mat->kd_txt, mat->kd_txt_info);
                gmr->metallicRoughnessTexture =
                    add_texture_info(mat->ks_txt, mat->ks_txt_info);
            } break;
            case material_type::specular_glossiness: {
                gmat->pbrSpecularGlossiness =
                    new glTFMaterialPbrSpecularGlossiness();
                auto gsg = gmat->pbrSpecularGlossiness;
                gsg->diffuseFactor = {
                    mat->kd[0], mat->kd[1], mat->kd[2], mat->op};
                gsg->specularFactor = mat->ks;
                gsg->glossinessFactor = mat->rs;
                gsg->diffuseTexture =
                    add_texture_info(mat->kd_txt, mat->kd_txt_info);
                gsg->specularGlossinessTexture =
                    add_texture_info(mat->ks_txt, mat->ks_txt_info);
            } break;
        }
        gmat->normalTexture = (glTFMaterialNormalTextureInfo*)add_texture_info(
            mat->norm_txt, mat->norm_txt_info, true, false);
        gmat->occlusionTexture =
            (glTFMaterialOcclusionTextureInfo*)add_texture_info(
                mat->occ_txt, mat->occ_txt_info, false, true);
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
            switch (type) {
                case glTFAccessorType::Scalar: {
                    auto bbox = make_bbox(count, (vec1f*)data);
                    accessor->min = {bbox.min.x};
                    accessor->max = {bbox.max.x};
                } break;
                case glTFAccessorType::Vec2: {
                    auto bbox = make_bbox(count, (vec2f*)data);
                    accessor->min = {bbox.min.x, bbox.min.y};
                    accessor->max = {bbox.max.x, bbox.max.y};
                } break;
                case glTFAccessorType::Vec3: {
                    auto bbox = make_bbox(count, (vec3f*)data);
                    accessor->min = {bbox.min.x, bbox.min.y, bbox.min.z};
                    accessor->max = {bbox.max.x, bbox.max.y, bbox.max.z};
                } break;
                case glTFAccessorType::Vec4: {
                    auto bbox = make_bbox(count, (vec4f*)data);
                    accessor->min = {
                        bbox.min.x, bbox.min.y, bbox.min.z, bbox.min.w};
                    accessor->max = {
                        bbox.max.x, bbox.max.y, bbox.max.z, bbox.max.w};
                } break;
                default: break;
            }
        }
        return glTFid<glTFAccessor>((int)gltf->accessors.size() - 1);
    };

    // convert meshes
    for (auto sgr : scn->shapes) {
        auto gmesh = new glTFMesh();
        gmesh->name = sgr->name;
        auto gbuffer = add_opt_buffer(sgr->path);
        for (auto shp : sgr->shapes) {
            auto gprim = new glTFMeshPrimitive();
            gprim->material =
                glTFid<glTFMaterial>(index(scn->materials, shp->mat));
            if (!shp->pos.empty())
                gprim->attributes["POSITION"] = add_accessor(gbuffer,
                    shp->name + "_pos", glTFAccessorType::Vec3,
                    glTFAccessorComponentType::Float, (int)shp->pos.size(),
                    sizeof(vec3f), shp->pos.data(), true);
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
                    glTFAccessorComponentType::Float,
                    (int)shp->texcoord1.size(), sizeof(vec2f),
                    shp->texcoord1.data(), false);
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
                    glTFAccessorComponentType::UnsignedInt,
                    (int)shp->points.size(), sizeof(int),
                    (int*)shp->points.data(), false);
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
                    (int)triangles.size() * 3, sizeof(int),
                    (int*)triangles.data(), false);
                gprim->mode = glTFMeshPrimitiveMode::Triangles;
            } else if (!shp->quads_pos.empty()) {
                throw std::runtime_error("face varying not supported in glTF");
            } else {
                throw std::runtime_error("empty mesh");
            }
            gmesh->primitives.push_back(gprim);
        }
        gltf->meshes.push_back(gmesh);
    }

    // hierarchy
    if (scn->nodes.empty()) {
        // instances
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
        // nodes
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
            gnode->scale = nde->scaling;
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
        std::map<keyframe_type, glTFAnimationSamplerInterpolation>{
            {keyframe_type::step, glTFAnimationSamplerInterpolation::Step},
            {keyframe_type::linear, glTFAnimationSamplerInterpolation::Linear},
            {keyframe_type::bezier,
                glTFAnimationSamplerInterpolation::CubicSpline},
            {keyframe_type::catmull_rom,
                glTFAnimationSamplerInterpolation::CatmullRomSpline},
        };

    // animation
    for (auto agr : scn->animations) {
        auto ganm = new glTFAnimation();
        ganm->name = agr->name;
        auto gbuffer = add_opt_buffer(agr->path);
        auto count = 0;
        auto paths =
            std::unordered_map<animation*, glTFAnimationChannelTargetPath>();
        for (auto anm : agr->animations) {
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
            } else if (!anm->scaling.empty()) {
                gsmp->output = add_accessor(gbuffer, aid + "_scale",
                    glTFAccessorType::Vec3, glTFAccessorComponentType::Float,
                    (int)anm->scaling.size(), sizeof(vec3f),
                    anm->scaling.data(), false);
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
            ganm->samplers.push_back(gsmp);
            paths[anm] = path;
        }
        for (auto target : agr->targets) {
            auto kfr = target.first;
            auto node = target.second;
            auto gchan = new glTFAnimationChannel();
            gchan->sampler =
                glTFid<glTFAnimationSampler>{index(agr->animations, kfr)};
            gchan->target = new glTFAnimationChannelTarget();
            gchan->target->node = glTFid<glTFNode>{index(scn->nodes, node)};
            gchan->target->path = paths.at(kfr);
            ganm->channels.push_back(gchan);
        }

        gltf->animations.push_back(ganm);
    }

    // done
    return gltf.release();
}

// Save a gltf scene
void save_gltf_scene(
    const std::string& filename, const scene* scn, const save_options& opts) {
    auto buffer_uri = path_basename(filename) + ".bin";
    auto gscn = std::unique_ptr<glTF>(
        scene_to_gltf(scn, buffer_uri, opts.gltf_separate_buffers));
    save_gltf(filename, gscn.get(), true, opts.save_textures);
}

#endif

#if YGL_SVG

// Converts an svg scene
scene* svg_to_scene(const svg_scene* sscn, const load_options& opts) {
    auto scn = new scene();
    auto sid = 0;
    for (auto sshp : sscn->shapes) {
        // TODO: recursive shape
        auto sgr = new shape_group();
        sgr->name = "shape" + std::to_string(sid++);
        auto shp = new shape();
        shp->name = sgr->name;
        for (auto spth : sshp->paths) {
            auto o = (int)shp->pos.size();
            for (auto p : spth->pos) shp->pos.push_back({p.x, p.y, 0});
            for (auto vid = 1; vid < spth->pos.size(); vid += 3)
                shp->beziers.push_back(
                    {o + vid - 1, o + vid + 0, o + vid + 1, o + vid + 2});
        }
        sgr->shapes.push_back(shp);
        scn->shapes.push_back(sgr);
    }
    auto miny = flt_max, maxy = -flt_max;
    for (auto sgr : scn->shapes) {
        for (auto shp : sgr->shapes) {
            for (auto& p : shp->pos) {
                miny = min(miny, p.y);
                maxy = max(maxy, p.y);
            }
        }
    }
    auto mdly = (maxy + miny) / 2;
    for (auto sgr : scn->shapes) {
        for (auto shp : sgr->shapes) {
            for (auto& p : shp->pos) p.y = -(p.y - mdly) + mdly;
        }
    }
    return scn;
}

// Load an svg scene
scene* load_svg_scene(const std::string& filename, const load_options& opts) {
    auto sscn = std::unique_ptr<svg_scene>(load_svg(filename));
    auto scn = std::unique_ptr<scene>(svg_to_scene(sscn.get(), opts));
    return scn.release();
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
    return ks +
           (vec3f{1, 1, 1} - ks) * pow(clamp(1.0f - cosw, 0.0f, 1.0f), 5.0f);
}

// Schlick approximation of Fresnel term weighted by roughness.
// This is a hack, but works better than not doing it.
vec3f fresnel_schlick(const vec3f& ks, float cosw, float rs) {
    auto fks = fresnel_schlick(ks, cosw);
    return lerp(ks, fks, rs);
}

// Evaluates the GGX distribution and geometric term
float eval_ggx(float rs, float ndh, float ndi, float ndo) {
    // evaluate GGX
    auto alpha2 = rs * rs;
    auto di = (ndh * ndh) * (alpha2 - 1) + 1;
    auto d = alpha2 / (pif * di * di);
#ifndef YGL_GGX_SMITH
    auto lambda_o = (-1 + sqrt(1 + alpha2 * (1 - ndo * ndo) / (ndo * ndo))) / 2;
    auto lambda_i = (-1 + sqrt(1 + alpha2 * (1 - ndi * ndi) / (ndi * ndi))) / 2;
    auto g = 1 / (1 + lambda_o + lambda_i);
#else
    auto go = (2 * ndo) / (ndo + sqrt(alpha2 + (1 - alpha2) * ndo * ndo));
    auto gi = (2 * ndi) / (ndi + sqrt(alpha2 + (1 - alpha2) * ndi * ndi));
    auto g = go * gi;
#endif
    return d * g;
}

// Evaluates the GGX pdf
float sample_ggx_pdf(float rs, float ndh) {
    auto cos2 = ndh * ndh;
    auto tan2 = (1 - cos2) / cos2;
    auto alpha2 = rs * rs;
    auto d = alpha2 / (pif * cos2 * cos2 * (alpha2 + tan2) * (alpha2 + tan2));
    return d;
}

// Sample the GGX distribution
vec3f sample_ggx(float rs, const vec2f& rn) {
    auto tan2 = rs * rs * rn.y / (1 - rn.y);
    auto rz = sqrt(1 / (tan2 + 1)), rr = sqrt(1 - rz * rz),
         rphi = 2 * pif * rn.x;
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
float sample_next1f(trace_pixel& pxl, trace_rng_type type, int nsamples) {
    switch (type) {
        case trace_rng_type::uniform: {
            return clamp(next_rand1f(pxl.rng), 0.0f, 1 - flt_eps);
        } break;
        case trace_rng_type::stratified: {
            auto p = hash_uint64_32((uint64_t)pxl.i | (uint64_t)pxl.j << 16 |
                                    (uint64_t)pxl.dimension << 32);
            auto s = cmjs_permute(pxl.sample, nsamples, p);
            pxl.dimension += 1;
            return clamp(
                (s + next_rand1f(pxl.rng)) / nsamples, 0.0f, 1 - flt_eps);
        } break;
        default: {
            assert(false);
            return 0;
        }
    }
}

// Generates a 2-dimensional sample.
vec2f sample_next2f(trace_pixel& pxl, trace_rng_type type, int nsamples) {
    switch (type) {
        case trace_rng_type::uniform: {
            return {next_rand1f(pxl.rng), next_rand1f(pxl.rng)};
        } break;
        case trace_rng_type::stratified: {
            auto p = hash_uint64_32((uint64_t)pxl.i | (uint64_t)pxl.j << 16 |
                                    (uint64_t)pxl.dimension << 32);
            auto s = cmjs_permute(pxl.sample, nsamples, p);
            auto nsamples2 = (int)round(sqrt(nsamples));
            pxl.dimension += 2;
            return {clamp((s % nsamples2 + next_rand1f(pxl.rng)) / nsamples2,
                        0.0f, 1 - flt_eps),
                clamp((s / nsamples2 + next_rand1f(pxl.rng)) / nsamples2, 0.0f,
                    1 - flt_eps)};
        } break;
        default: {
            assert(false);
            return {0, 0};
        }
    }
}

// Surface point with geometry and material data. Supports point on
// envmap too. This is the key data manipulated in the path tracer.
struct trace_point {
    const shape* shp = nullptr;        // shape
    const environment* env = nullptr;  // environment
    vec3f pos = zero3f;                // pos
    vec3f norm = {0, 0, 1};            // norm
    vec2f texcoord = zero2f;           // texcoord
    vec3f ke = zero3f;                 // emission
    vec3f kd = {0, 0, 0};              // diffuse
    vec3f ks = {0, 0, 0};              // specular
    float rs = 0;                      // specular roughness
    vec3f kt = {0, 0, 0};              // transmission (thin glass)
    float op = 1.0f;                   // opacity
    bool has_brdf() const { return shp && kd + ks + kt != zero3f; }
    vec3f rho() const { return kd + ks + kt; }
    vec3f brdf_weights() const {
        auto w = vec3f{max_element_value(kd), max_element_value(ks),
            max_element_value(kt)};
        auto sw = w.x + w.y + w.z;
        if (!sw) return zero3f;
        return w / sw;
    }
};

// Evaluates emission.
vec3f eval_emission(const trace_point& pt, const vec3f& wo) {
    if (pt.shp && !pt.shp->triangles.empty() && dot(pt.norm, wo) <= 0)
        return zero3f;
    return pt.ke;
}

// Evaluates the BRDF scaled by the cosine of the incoming direction.
// - ggx from [Heitz 2014] and [Walter 2007] and [Lagarde 2014]
// "Understanding the Masking-Shadowing Function in Microfacet-Based
// BRDFs" http://jcgt.org/published/0003/02/03/
// - "Microfacet Models for Refraction through Rough Surfaces" EGSR 07
// https://www.cs.cornell.edu/~srm/publications/EGSR07-btdf.pdf
vec3f eval_ggx_brdfcos(const trace_point& pt, const vec3f& wo, const vec3f& wi,
    bool delta = false) {
    auto brdfcos = zero3f;
    auto wh = normalize(wo + wi);

    auto ndo = dot(pt.norm, wo), ndi = dot(pt.norm, wi),
         ndh = clamp(dot(wh, pt.norm), -1.0f, 1.0f);

    if (ndi > 0 && ndo > 0) {
        brdfcos += pt.kd * ndi / pif;
        if (ndh > 0 && pt.rs) {
            auto dg = eval_ggx(pt.rs, ndh, ndi, ndo);
            auto odh = clamp(dot(wo, wh), 0.0f, 1.0f);
            auto ks = fresnel_schlick(pt.ks, odh, pt.rs);
            brdfcos += ks * ndi * dg / (4 * ndi * ndo);
        }
        if (!pt.rs && delta) {
            auto ks = fresnel_schlick(pt.ks, ndo, pt.rs);
            brdfcos += ks;
        }
    }
    if (wo == -wi && delta) brdfcos += pt.kt;

    assert(isfinite(brdfcos.x) && isfinite(brdfcos.y) && isfinite(brdfcos.z));
    return brdfcos;
}

// Evaluates the BRDF scaled by the cosine of the incoming direction.
// - uses Kajiya-Kay for hair
vec3f eval_kajiyakay_brdfcos(const trace_point& pt, const vec3f& wo,
    const vec3f& wi, bool delta = false) {
    auto brdfcos = zero3f;
    auto wh = normalize(wo + wi);

    auto ndo = dot(pt.norm, wo), ndi = dot(pt.norm, wi),
         ndh = clamp(dot(pt.norm, wh), 0.0f, 1.0f);
    auto so = sqrt(clamp(1 - ndo * ndo, 0.0f, 1.0f)),
         si = sqrt(clamp(1 - ndi * ndi, 0.0f, 1.0f)),
         sh = sqrt(clamp(1 - ndh * ndh, 0.0f, 1.0f));

    if (si > 0 && so > 0) {
        brdfcos += pt.kd * si / pif;
        if (sh > 0 && pt.rs) {
            auto ns = 2 / (pt.rs * pt.rs) - 2;
            auto d = (ns + 2) * pow(sh, ns) / (2 + pif);
            brdfcos += pt.ks * si * d / (4.0f * si * so);
        }
    }
    if (wo == -wi && delta) brdfcos += pt.kt;

    assert(isfinite(brdfcos.x) && isfinite(brdfcos.y) && isfinite(brdfcos.z));
    return brdfcos;
}

// Evaluates the BRDF scaled by the cosine of the incoming direction.
// - uses a hack for points
vec3f eval_point_brdfcos(const trace_point& pt, const vec3f& wo,
    const vec3f& wi, bool delta = false) {
    auto brdfcos = zero3f;

    auto ido = dot(wo, wi);
    brdfcos += pt.kd * (2 * ido + 1) / (2 * pif);
    if (wo == -wi && delta) brdfcos += pt.kt;

    assert(isfinite(brdfcos.x) && isfinite(brdfcos.y) && isfinite(brdfcos.z));
    return brdfcos;
}

// Evaluates the BRDF scaled by the cosine of the incoming direction.
vec3f eval_brdfcos(const trace_point& pt, const vec3f& wo, const vec3f& wi,
    bool delta = false) {
    if (!pt.has_brdf()) return zero3f;
    if (!pt.shp->triangles.empty())
        return eval_ggx_brdfcos(pt, wo, wi, delta);
    else if (!pt.shp->lines.empty())
        return eval_kajiyakay_brdfcos(pt, wo, wi, delta);
    else if (!pt.shp->points.empty())
        return eval_point_brdfcos(pt, wo, wi, delta);
    else
        return zero3f;
}

// Compute the weight for sampling the BRDF
float weight_ggx_brdfcos(const trace_point& pt, const vec3f& wo,
    const vec3f& wi, bool delta = false) {
    auto weights = pt.brdf_weights();
    auto wh = normalize(wi + wo);
    auto ndo = dot(pt.norm, wo), ndi = dot(pt.norm, wi), ndh = dot(pt.norm, wh);

    auto pdf = 0.0f;
    if (ndo > 0 && ndi > 0) {
        pdf += weights.x * ndi / pif;
        if (ndh > 0 && pt.rs) {
            auto d = sample_ggx_pdf(pt.rs, ndh);
            auto hdo = dot(wo, wh);
            pdf += weights.y * d / (4 * hdo);
        }
        if (!pt.rs && delta) pdf += weights.y;
    }
    if (wi == -wo && delta) pdf += weights.z;

    assert(isfinite(pdf));
    if (!pdf) return 0;
    return 1 / pdf;
}

// Compute the weight for sampling the BRDF
float weight_kajiyakay_brdfcos(const trace_point& pt, const vec3f& wo,
    const vec3f& wi, bool delta = false) {
    auto weights = pt.brdf_weights();

    auto pdf = 0.0f;
    pdf += (weights.x + weights.y) / (4 * pif);
    if (wi == -wo && delta) pdf += weights.z;

    assert(isfinite(pdf));
    if (!pdf) return 0;
    return 1 / pdf;
}

// Compute the weight for sampling the BRDF
float weight_point_brdfcos(const trace_point& pt, const vec3f& wo,
    const vec3f& wi, bool delta = false) {
    auto weights = pt.brdf_weights();

    auto pdf = 0.0f;
    pdf += (weights.x + weights.y) / (4 * pif);
    if (wi == -wo && delta) pdf += weights.z;

    assert(isfinite(pdf));
    if (!pdf) return 0;
    return 1 / pdf;
}

// Compute the weight for sampling the BRDF
float weight_brdfcos(const trace_point& pt, const vec3f& wo, const vec3f& wi,
    bool delta = false) {
    if (!pt.has_brdf()) return 0;
    if (!pt.shp->triangles.empty())
        return weight_ggx_brdfcos(pt, wo, wi, delta);
    else if (!pt.shp->lines.empty())
        return weight_kajiyakay_brdfcos(pt, wo, wi, delta);
    else if (!pt.shp->points.empty())
        return weight_point_brdfcos(pt, wo, wi, delta);
    else
        return 0;
}

// Picks a direction based on the BRDF
std::tuple<vec3f, bool> sample_ggx_brdfcos(
    const trace_point& pt, const vec3f& wo, float rnl, const vec2f& rn) {
    auto weights = pt.brdf_weights();
    auto ndo = dot(pt.norm, wo);
    if (ndo <= 0) return {zero3f, false};

    // sample according to diffuse
    if (rnl < weights.x) {
        auto fp = make_frame_fromz(pt.pos, pt.norm);
        auto rz = sqrtf(rn.y), rr = sqrtf(1 - rz * rz), rphi = 2 * pif * rn.x;
        auto wi_local = vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
        return {transform_direction(fp, wi_local), false};
    }
    // sample according to specular GGX
    else if (rnl < weights.x + weights.y && pt.rs) {
        auto fp = make_frame_fromz(pt.pos, pt.norm);
        auto wh_local = sample_ggx(pt.rs, rn);
        auto wh = transform_direction(fp, wh_local);
        return {normalize(wh * 2.0f * dot(wo, wh) - wo), false};
    }
    // sample according to specular mirror
    else if (rnl < weights.x + weights.y && !pt.rs) {
        return {normalize(pt.norm * 2.0f * dot(wo, pt.norm) - wo), true};
    }
    // transmission hack
    else if (rnl < weights.x + weights.y + weights.z) {
        return {-wo, true};
    } else
        assert(false);

    return {zero3f, false};
}

// Picks a direction based on the BRDF
std::tuple<vec3f, bool> sample_kajiyakay_brdfcos(
    const trace_point& pt, const vec3f& wo, float rnl, const vec2f& rn) {
    auto weights = pt.brdf_weights();
    // diffuse and specular: samnple a uniform spherical direction
    if (rnl < weights.x + weights.y) {
        auto fp = make_frame_fromz(pt.pos, pt.norm);
        auto rz = 2 * rn.y - 1, rr = sqrtf(1 - rz * rz), rphi = 2 * pif * rn.x;
        auto wi_local = vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
        return {transform_direction(fp, wi_local), false};
    }
    // transmission hack
    else if (rnl < weights.x + weights.y + weights.z) {
        return {-wo, true};
    } else
        assert(false);
    return {zero3f, false};
}

// Picks a direction based on the BRDF
std::tuple<vec3f, bool> sample_point_brdfcos(
    const trace_point& pt, const vec3f& wo, float rnl, const vec2f& rn) {
    auto weights = pt.brdf_weights();
    // diffuse and specular: samnple a uniform spherical direction
    if (rnl < weights.x + weights.y) {
        auto fp = make_frame_fromz(pt.pos, pt.norm);
        auto rz = 2 * rn.y - 1, rr = sqrtf(1 - rz * rz), rphi = 2 * pif * rn.x;
        auto wi_local = vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
        return {transform_direction(fp, wi_local), false};
    }
    // transmission hack
    else if (rnl < weights.x + weights.y + weights.z) {
        // continue ray direction
        return {-wo, true};
    } else
        assert(false);
    return {zero3f, false};
}

// Picks a direction based on the BRDF
std::tuple<vec3f, bool> sample_brdfcos(
    const trace_point& pt, const vec3f& wo, float rnl, const vec2f& rn) {
    if (!pt.has_brdf()) return {zero3f, false};
    if (!pt.shp->triangles.empty())
        return sample_ggx_brdfcos(pt, wo, rnl, rn);
    else if (!pt.shp->lines.empty())
        return sample_kajiyakay_brdfcos(pt, wo, rnl, rn);
    else if (!pt.shp->points.empty())
        return sample_point_brdfcos(pt, wo, rnl, rn);
    else
        return {zero3f, false};
}

// Create a point for an environment map. Resolves material with
// textures.
trace_point eval_point(const environment* env, const vec3f& wo) {
    auto pt = trace_point();
    pt.env = env;
    pt.pos = wo * flt_max;
    pt.norm = -wo;
    pt.ke = env->ke;
    if (env->ke_txt) {
        auto w = transform_direction_inverse(env->frame, -wo);
        auto theta = acos(clamp(w.y, -1.0f, 1.0f));
        auto phi = atan2(w.z, w.x);
        auto texcoord = vec2f{0.5f + phi / (2 * pif), theta / pif};
        auto txt = eval_texture(env->ke_txt, env->ke_txt_info, texcoord);
        pt.ke *= {txt.x, txt.y, txt.z};
    }
    return pt;
}

// Create a point for a shape. Resolves geometry and material with
// textures.
trace_point eval_point(
    const instance* ist, int sid, int eid, const vec2f& euv, const vec3f& wo) {
    // default material
    static auto def_material = (material*)nullptr;
    if (!def_material) {
        def_material = new material();
        def_material->kd = {0.2f, 0.2f, 0.2f};
        def_material->rs = 1;
    }

    // point
    auto pt = trace_point();
    pt.shp = ist->shp->shapes.at(sid);
    pt.pos = eval_pos(pt.shp, eid, euv);
    pt.norm = eval_norm(pt.shp, eid, euv);
    pt.texcoord = eval_texcoord(pt.shp, eid, euv);
    // shortcuts
    auto mat = (pt.shp->mat) ? pt.shp->mat : def_material;

    // handle normal map
    if (mat->norm_txt) {
        auto tangsp = eval_tangsp(pt.shp, eid, euv);
        auto txt = eval_texture(
                       mat->norm_txt, mat->norm_txt_info, pt.texcoord, false) *
                       2.0f -
                   vec4f{1};
        auto ntxt = normalize(vec3f{txt.x, -txt.y, txt.z});
        auto frame = make_frame_fromzx(
            {0, 0, 0}, pt.norm, {tangsp.x, tangsp.y, tangsp.z});
        frame.y *= tangsp.w;
        pt.norm = transform_direction(frame, ntxt);
    }

    // move to world coordinates
    pt.pos = transform_point(ist->frame, pt.pos);
    pt.norm = transform_direction(ist->frame, pt.norm);

    // correct for double sided
    if (mat->double_sided && dot(pt.norm, wo) < 0) pt.norm = -pt.norm;

    // initialized material values
    auto kx = vec3f{1, 1, 1};
    pt.op = 1;
    if (!pt.shp->color.empty()) {
        auto col = eval_color(pt.shp, eid, euv);
        kx *= {col.x, col.y, col.z};
        pt.op *= col.w;
    }

    // handle occlusion
    if (mat->occ_txt) {
        auto txt = eval_texture(mat->occ_txt, mat->occ_txt_info, pt.texcoord);
        kx *= {txt.x, txt.y, txt.z};
    }

    // sample emission
    pt.ke = mat->ke * kx;
    if (mat->ke_txt) {
        auto txt = eval_texture(mat->ke_txt, mat->ke_txt_info, pt.texcoord);
        pt.ke *= {txt.x, txt.y, txt.z};
    }

    // sample reflectance
    switch (mat->type) {
        case material_type::specular_roughness: {
            pt.kd = mat->kd * kx;
            if (mat->kd_txt) {
                auto txt =
                    eval_texture(mat->kd_txt, mat->kd_txt_info, pt.texcoord);
                pt.kd *= {txt.x, txt.y, txt.z};
                pt.op *= txt.w;
            }
            pt.ks = mat->ks * kx;
            pt.rs = mat->rs;
            if (mat->ks_txt) {
                auto txt =
                    eval_texture(mat->ks_txt, mat->ks_txt_info, pt.texcoord);
                pt.ks *= {txt.x, txt.y, txt.z};
            }
            pt.kt = mat->kt * kx;
            if (mat->kt_txt) {
                auto txt =
                    eval_texture(mat->kt_txt, mat->kt_txt_info, pt.texcoord);
                pt.kt *= {txt.x, txt.y, txt.z};
            }
        } break;
        case material_type::metallic_roughness: {
            auto kb = mat->kd * kx;
            if (mat->kd_txt) {
                auto txt =
                    eval_texture(mat->kd_txt, mat->kd_txt_info, pt.texcoord);
                kb *= {txt.x, txt.y, txt.z};
                pt.op *= txt.w;
            }
            auto km = mat->ks.x;
            pt.rs = mat->rs;
            if (mat->ks_txt) {
                auto txt =
                    eval_texture(mat->ks_txt, mat->ks_txt_info, pt.texcoord);
                km *= txt.y;
                pt.rs *= txt.z;
            }
            pt.kd = kb * (1 - km);
            pt.ks = kb * km + vec3f{0.04f} * (1 - km);
        } break;
        case material_type::specular_glossiness: {
            pt.kd = mat->kd * kx;
            if (mat->kd_txt) {
                auto txt =
                    eval_texture(mat->kd_txt, mat->kd_txt_info, pt.texcoord);
                pt.kd *= {txt.x, txt.y, txt.z};
                pt.op *= txt.w;
            }
            pt.ks = mat->ks * kx;
            pt.rs = mat->rs;
            if (mat->ks_txt) {
                auto txt =
                    eval_texture(mat->ks_txt, mat->ks_txt_info, pt.texcoord);
                pt.ks *= {txt.x, txt.y, txt.z};
                pt.rs *= txt.w;
            }
            pt.rs = 1 - pt.rs;  // glossiness -> roughnes
            pt.kt = mat->kt * kx;
            if (mat->kt_txt) {
                auto txt =
                    eval_texture(mat->kt_txt, mat->kt_txt_info, pt.texcoord);
                pt.kt *= {txt.x, txt.y, txt.z};
            }
        } break;
    }

    // set up final values
    pt.ke *= pt.op;
    pt.kd *= pt.op;
    if (pt.ks != zero3f && pt.rs < 0.9999f) {
        pt.ks *= pt.op;
        pt.rs = pt.rs * pt.rs;
    } else {
        pt.ks = zero3f;
        pt.rs = 0;
    }
    if (pt.kt == zero3f) pt.kt = vec3f{1 - pt.op};

    // done
    return pt;
}

// Sample weight for a light point.
float weight_light(
    const trace_lights& lights, const trace_point& lpt, const trace_point& pt) {
    if (lpt.shp) {
        auto dist = length(lpt.pos - pt.pos);
        auto area = lights.shape_areas.at(lpt.shp);
        if (!lpt.shp->triangles.empty()) {
            return area * abs(dot(lpt.norm, normalize(lpt.pos - pt.pos))) /
                   (dist * dist);
        } else if (!lpt.shp->lines.empty()) {
            // TODO: fixme
            return 0;
        } else if (!lpt.shp->points.empty()) {
            return area / (dist * dist);
        }
    }
    if (lpt.env) { return 4 * pif; }
    return 0;
}

// Sample weight for a light point.
float weight_lights(
    const trace_lights& lights, const trace_point& lpt, const trace_point& pt) {
    return lights.lights.size() * weight_light(lights, lpt, pt);
}

// Picks a point on a light.
trace_point sample_light(const trace_lights& lights, const trace_light& lgt,
    const trace_point& pt, float rel, const vec2f& ruv) {
    if (lgt.ist) {
        auto shp = lgt.ist->shp->shapes.at(0);
        auto& cdf = lights.shape_cdfs.at(shp);
        auto eid = 0;
        auto euv = zero2f;
        if (!shp->triangles.empty()) {
            std::tie(eid, euv) = sample_triangles(cdf, rel, ruv);
        } else if (!shp->lines.empty()) {
            std::tie(eid, (float&)euv) = sample_lines(cdf, rel, ruv.x);
        } else if (!shp->lines.empty()) {
            eid = sample_points(cdf, rel);
        }
        return eval_point(lgt.ist, 0, eid, euv, zero3f);
    }
    if (lgt.env) {
        auto z = -1 + 2 * ruv.y;
        auto rr = sqrt(clamp(1 - z * z, 0.0f, 1.0f));
        auto phi = 2 * pif * ruv.x;
        auto wo = vec3f{cos(phi) * rr, z, sin(phi) * rr};
        return eval_point(lgt.env, wo);
    }
    return {};
}

// Picks a point on a light.
trace_point sample_lights(const trace_lights& lights, const trace_point& pt,
    float rnl, float rne, const vec2f& ruv) {
    auto& lgt = lights.lights.at(
        clamp((int)(rnl * lights.lights.size()), 0, (int)lights.lights.size()));
    return sample_light(lights, lgt, pt, rne, ruv);
}

// Intersects a ray with the scn and return the point (or env
// point).
trace_point intersect_scene(
    const scene* scn, const bvh_tree* bvh, const ray3f& ray) {
    auto iid = 0, sid = 0, eid = 0;
    auto euv = zero2f;
    auto ray_t = 0.0f;
    if (intersect_bvh(bvh, ray, false, ray_t, iid, sid, eid, euv)) {
        return eval_point(scn->instances[iid], sid, eid, euv, -ray.d);
    } else if (!scn->environments.empty()) {
        return eval_point(scn->environments[0], -ray.d);
    } else {
        return {};
    }
}

// Test occlusion
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
            cpt = intersect_scene(scn, bvh, ray);
            if (!cpt.shp) break;
            weight *= cpt.kt;
            if (weight == zero3f) break;
        }
        return weight;
    }
}

// Mis weight
float weight_mis(float w0, float w1) {
    if (!w0 || !w1) return 1;
    return (1 / w0) / (1 / w0 + 1 / w1);
}

// Recursive path tracing.
vec3f trace_path(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt_, const vec3f& wo_,
    trace_pixel& pxl, const trace_params& params) {
    auto pt = pt_;
    auto wo = wo_;

    // emission
    auto l = eval_emission(pt, wo);
    if (!pt.has_brdf() || lights.empty()) return l;

    // trace path
    auto weight = vec3f{1, 1, 1};
    auto emission = false;
    for (auto bounce = 0; bounce < params.max_depth; bounce++) {
        // emission
        if (emission) l += weight * eval_emission(pt, wo);

        // direct – light
        auto rll = sample_next1f(pxl, params.rng, params.nsamples);
        auto rle = sample_next1f(pxl, params.rng, params.nsamples);
        auto rluv = sample_next2f(pxl, params.rng, params.nsamples);
        auto& lgt = lights.lights[(int)(rll * lights.lights.size())];
        auto lpt = sample_light(lights, lgt, pt, rle, rluv);
        auto lw = weight_light(lights, lpt, pt) * (float)lights.size();
        auto lwi = normalize(lpt.pos - pt.pos);
        auto lke = eval_emission(lpt, -lwi);
        auto lbc = eval_brdfcos(pt, wo, lwi);
        auto lld = lke * lbc * lw;
        if (lld != zero3f) {
            l += weight * lld * eval_transmission(scn, bvh, pt, lpt, params) *
                 weight_mis(lw, weight_brdfcos(pt, wo, lwi));
        }

        // direct – brdf
        auto rbl = sample_next1f(pxl, params.rng, params.nsamples);
        auto rbuv = sample_next2f(pxl, params.rng, params.nsamples);
        auto bwi = zero3f;
        auto bdelta = false;
        std::tie(bwi, bdelta) = sample_brdfcos(pt, wo, rbl, rbuv);
        auto bpt = intersect_scene(scn, bvh, make_ray(pt.pos, bwi));
        auto bw = weight_brdfcos(pt, wo, bwi, bdelta);
        auto bke = eval_emission(bpt, -bwi);
        auto bbc = eval_brdfcos(pt, wo, bwi, bdelta);
        auto bld = bke * bbc * bw;
        if (bld != zero3f) {
            l += weight * bld * weight_mis(bw, weight_light(lights, bpt, pt));
        }

        // skip recursion if path ends
        if (bounce == params.max_depth - 1) break;
        if (!bpt.has_brdf()) break;

        // continue path
        weight *= eval_brdfcos(pt, wo, bwi) * weight_brdfcos(pt, wo, bwi);
        if (weight == zero3f) break;

        // roussian roulette
        if (bounce > 2) {
            auto rrprob = 1.0f - min(max_element_value(pt.rho()), 0.95f);
            if (sample_next1f(pxl, params.rng, params.nsamples) < rrprob) break;
            weight *= 1 / (1 - rrprob);
        }

        // continue path
        pt = bpt;
        wo = -bwi;
        emission = false;
    }

    return l;
}

// Recursive path tracing.
vec3f trace_path_nomis(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt_, const vec3f& wo_,
    trace_pixel& pxl, const trace_params& params) {
    // emission
    auto pt = pt_;
    auto wo = wo_;

    auto l = eval_emission(pt, wo);
    if (!pt.has_brdf() || lights.empty()) return l;

    // trace path
    auto weight = vec3f{1, 1, 1};
    auto emission = false;
    for (auto bounce = 0; bounce < params.max_depth; bounce++) {
        // emission
        if (emission) l += weight * eval_emission(pt, wo);

        // direct
        auto rll = sample_next1f(pxl, params.rng, params.nsamples);
        auto rle = sample_next1f(pxl, params.rng, params.nsamples);
        auto rluv = sample_next2f(pxl, params.rng, params.nsamples);
        auto& lgt = lights.lights[(int)(rll * lights.lights.size())];
        auto lpt = sample_light(lights, lgt, pt, rle, rluv);
        auto lwi = normalize(lpt.pos - pt.pos);
        auto ld = eval_emission(lpt, -lwi) * eval_brdfcos(pt, wo, lwi) *
                  weight_light(lights, lpt, pt) * (float)lights.size();
        if (ld != zero3f) {
            l += weight * ld * eval_transmission(scn, bvh, pt, lpt, params);
        }

        // skip recursion if path ends
        if (bounce == params.max_depth - 1) break;

        // roussian roulette
        if (bounce > 2) {
            auto rrprob = 1.0f - min(max_element_value(pt.rho()), 0.95f);
            if (sample_next1f(pxl, params.rng, params.nsamples) < rrprob) break;
            weight *= 1 / (1 - rrprob);
        }

        // continue path
        auto rbl = sample_next1f(pxl, params.rng, params.nsamples);
        auto rbuv = sample_next2f(pxl, params.rng, params.nsamples);
        auto bwi = zero3f;
        auto bdelta = false;
        std::tie(bwi, bdelta) = sample_brdfcos(pt, wo, rbl, rbuv);
        weight *= eval_brdfcos(pt, wo, bwi, bdelta) *
                  weight_brdfcos(pt, wo, bwi, bdelta);
        if (weight == zero3f) break;

        auto bpt = intersect_scene(scn, bvh, make_ray(pt.pos, bwi));
        emission = false;
        if (!bpt.has_brdf()) break;

        // continue path
        pt = bpt;
        wo = -bwi;
    }

    return l;
}

// Recursive path tracing.
vec3f trace_path_hack(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt_, const vec3f& wo_,
    trace_pixel& pxl, const trace_params& params) {
    auto pt = pt_;
    auto wo = wo_;

    // emission
    auto l = eval_emission(pt, wo);
    if (!pt.has_brdf() || lights.empty()) return l;

    // trace path
    auto weight = vec3f{1, 1, 1};
    for (auto bounce = 0; bounce < params.max_depth; bounce++) {
        // direct
        auto rll = sample_next1f(pxl, params.rng, params.nsamples);
        auto rle = sample_next1f(pxl, params.rng, params.nsamples);
        auto rluv = sample_next2f(pxl, params.rng, params.nsamples);
        auto& lgt = lights.lights[(int)(rll * lights.lights.size())];
        auto lpt = sample_light(lights, lgt, pt, rle, rluv);
        auto lwi = normalize(lpt.pos - pt.pos);
        auto ld = eval_emission(lpt, -lwi) * eval_brdfcos(pt, wo, -lwi) *
                  weight_light(lights, lpt, pt) * (float)lights.size();
        if (ld != zero3f) {
            l += weight * ld * eval_transmission(scn, bvh, pt, lpt, params);
        }

        // skip recursion if path ends
        if (bounce == params.max_depth - 1) break;

        // roussian roulette
        if (bounce > 2) {
            auto rrprob = 1.0f - min(max_element_value(pt.rho()), 0.95f);
            if (sample_next1f(pxl, params.rng, params.nsamples) < rrprob) break;
            weight *= 1 / (1 - rrprob);
        }

        // continue path
        auto bwi = zero3f;
        auto bdelta = false;
        std::tie(bwi, bdelta) = sample_brdfcos(pt, wo,
            sample_next1f(pxl, params.rng, params.nsamples),
            sample_next2f(pxl, params.rng, params.nsamples));
        weight *= eval_brdfcos(pt, wo, bwi, bdelta) *
                  weight_brdfcos(pt, wo, bwi, bdelta);
        if (weight == zero3f) break;

        auto bpt = intersect_scene(scn, bvh, make_ray(pt.pos, bwi));
        if (!bpt.has_brdf()) break;

        // continue path
        pt = bpt;
        wo = -bwi;
    }

    return l;
}

// Direct illumination.
vec3f trace_direct(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt, const vec3f& wo,
    int bounce, trace_pixel& pxl, const trace_params& params) {
    // emission
    auto l = eval_emission(pt, wo);
    if (!pt.has_brdf()) return l;

    // ambient
    l += params.ambient * pt.rho();

    // direct
    for (auto& lgt : lights.lights) {
        auto rle = sample_next1f(pxl, params.rng, params.nsamples);
        auto rluv = sample_next2f(pxl, params.rng, params.nsamples);
        auto lpt = sample_light(lights, lgt, pt, rle, rluv);
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
        auto rpt = intersect_scene(scn, bvh, make_ray(pt.pos, wi));
        l += pt.ks *
             trace_direct(scn, bvh, lights, rpt, -wi, bounce + 1, pxl, params);
    }

    // opacity
    if (pt.kt != zero3f) {
        auto opt = intersect_scene(scn, bvh, make_ray(pt.pos, -wo));
        l += pt.kt *
             trace_direct(scn, bvh, lights, opt, wo, bounce + 1, pxl, params);
    }

    // done
    return l;
}

// Direct illumination.
vec3f trace_direct(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt, const vec3f& wo,
    trace_pixel& pxl, const trace_params& params) {
    return trace_direct(scn, bvh, lights, pt, wo, 0, pxl, params);
}

// Eyelight for quick previewing.
vec3f trace_eyelight(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt, const vec3f& wo,
    int bounce, trace_pixel& pxl, const trace_params& params) {
    // emission
    auto l = eval_emission(pt, wo);
    if (!pt.has_brdf()) return l;

    // brdf*light
    l += eval_brdfcos(pt, wo, wo) * pif;

    // opacity
    if (bounce >= params.max_depth) return l;
    if (pt.kt != zero3f) {
        auto opt = intersect_scene(scn, bvh, make_ray(pt.pos, -wo));
        l += pt.kt *
             trace_eyelight(scn, bvh, lights, opt, wo, bounce + 1, pxl, params);
    }

    // done
    return l;
}

// Eyelight for quick previewing.
vec3f trace_eyelight(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt, const vec3f& wo,
    trace_pixel& pxl, const trace_params& params) {
    return trace_eyelight(scn, bvh, lights, pt, wo, 0, pxl, params);
}

// Debug previewing.
vec3f trace_debug_normal(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt, const vec3f& wo,
    trace_pixel& pxl, const trace_params& params) {
    return pt.norm * 0.5f + vec3f{0.5f, 0.5f, 0.5f};
}

// Debug previewing.
vec3f trace_debug_albedo(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt, const vec3f& wo,
    trace_pixel& pxl, const trace_params& params) {
    return pt.rho();
}

// Debug previewing.
vec3f trace_debug_texcoord(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt, const vec3f& wo,
    trace_pixel& pxl, const trace_params& params) {
    return {pt.texcoord.x, pt.texcoord.y, 0};
}

// Trace shader function
using trace_shader = vec3f (*)(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt, const vec3f& wo,
    trace_pixel& pxl, const trace_params& params);

// Trace filter function
using trace_filter = float (*)(float);

// map to convert trace samplers
static auto trace_shaders = std::unordered_map<trace_shader_type, trace_shader>{
    {trace_shader_type::eyelight, trace_eyelight},
    {trace_shader_type::direct, trace_direct},
    {trace_shader_type::pathtrace, trace_path},
    {trace_shader_type::pathtrace_nomis, trace_path_nomis},
    {trace_shader_type::debug_albedo, trace_debug_albedo},
    {trace_shader_type::debug_normal, trace_debug_normal},
    {trace_shader_type::debug_texcoord, trace_debug_texcoord},
};

// map to convert trace filters
static auto trace_filters = std::unordered_map<trace_filter_type, trace_filter>{
    {trace_filter_type::box, (trace_filter) nullptr},
    {trace_filter_type::triangle, filter_triangle},
    {trace_filter_type::cubic, filter_cubic},
    {trace_filter_type::catmull_rom, filter_catmullrom},
    {trace_filter_type::mitchell, filter_mitchell},
};

// map to convert trace filters
static auto trace_filter_sizes = std::unordered_map<trace_filter_type, int>{
    {trace_filter_type::box, 0},
    {trace_filter_type::triangle, 1},
    {trace_filter_type::cubic, 2},
    {trace_filter_type::catmull_rom, 2},
    {trace_filter_type::mitchell, 2},
};

// Trace a single sample
void trace_sample(const scene* scn, const camera* cam, const bvh_tree* bvh,
    const trace_lights& lights, trace_pixel& pxl, trace_shader shader,
    const trace_params& params) {
    pxl.sample += 1;
    pxl.dimension = 0;
    auto crn = sample_next2f(pxl, params.rng, params.nsamples);
    auto lrn = sample_next2f(pxl, params.rng, params.nsamples);
    auto uv = vec2f{(pxl.i + crn.x) / (cam->aspect * params.resolution),
        1 - (pxl.j + crn.y) / params.resolution};
    auto ray = eval_camera_ray(cam, uv, lrn);
    auto pt = intersect_scene(scn, bvh, ray);
    if (!pt.shp && params.envmap_invisible) return;
    auto l = shader(scn, bvh, lights, pt, -ray.d, pxl, params);
    if (!isfinite(l.x) || !isfinite(l.y) || !isfinite(l.z)) {
        log_error("NaN detected");
        return;
    }
    if (params.pixel_clamp > 0) l = clamplen(l, params.pixel_clamp);
    pxl.col += l;
    pxl.alpha += 1;
}

// Trace the next nsamples.
void trace_samples(const scene* scn, const camera* cam, const bvh_tree* bvh,
    const trace_lights& lights, image4f& img, image<trace_pixel>& pixels,
    int nsamples, const trace_params& params) {
    auto shader = trace_shaders.at(params.shader);
    if (params.parallel) {
        auto nthreads = std::thread::hardware_concurrency();
        auto threads = std::vector<std::thread>();
        for (auto tid = 0; tid < std::thread::hardware_concurrency(); tid++) {
            threads.push_back(std::thread([=, &img, &pixels, &params]() {
                for (auto j = tid; j < img.height(); j += nthreads) {
                    for (auto i = 0; i < img.width(); i++) {
                        auto& pxl = pixels.at(i, j);
                        for (auto s = 0; s < nsamples; s++)
                            trace_sample(
                                scn, cam, bvh, lights, pxl, shader, params);
                        img.at(i, j) =
                            vec4f{pxl.col.x, pxl.col.y, pxl.col.z, pxl.alpha};
                        img.at(i, j) /= pxl.sample;
                    }
                }
            }));
        }
        for (auto& t : threads) t.join();
        threads.clear();
    } else {
        auto shader = trace_shaders.at(params.shader);
        for (auto j = 0; j < img.height(); j++) {
            for (auto i = 0; i < img.width(); i++) {
                auto& pxl = pixels.at(i, j);
                for (auto s = 0; s < params.nsamples; s++)
                    trace_sample(scn, cam, bvh, lights, pxl, shader, params);
                img.at(i, j) =
                    vec4f{pxl.col.x, pxl.col.y, pxl.col.z, pxl.alpha};
                img.at(i, j) /= pxl.sample;
            }
        }
    }
}

// Trace a filtered sample of samples
void trace_sample_filtered(const scene* scn, const camera* cam,
    const bvh_tree* bvh, const trace_lights& lights, image4f& img,
    trace_pixel& pxl, trace_shader shader, trace_filter filter, int filter_size,
    std::mutex& image_mutex, const trace_params& params) {
    pxl.sample += 1;
    pxl.dimension = 0;
    auto crn = sample_next2f(pxl, params.rng, params.nsamples);
    auto lrn = sample_next2f(pxl, params.rng, params.nsamples);
    auto uv = vec2f{(pxl.i + crn.x) / (cam->aspect * params.resolution),
        1 - (pxl.j + crn.y) / params.resolution};
    auto ray = eval_camera_ray(cam, uv, lrn);
    auto pt = intersect_scene(scn, bvh, ray);
    if (!pt.shp && params.envmap_invisible) return;
    auto l = shader(scn, bvh, lights, pt, -ray.d, pxl, params);
    if (!isfinite(l.x) || !isfinite(l.y) || !isfinite(l.z)) {
        log_error("NaN detected");
        return;
    }
    if (params.pixel_clamp > 0) l = clamplen(l, params.pixel_clamp);
    if (params.filter == trace_filter_type::box) {
        pxl.col += l;
        pxl.alpha += 1;
        pxl.weight += 1;
    } else {
        std::lock_guard<std::mutex> lock(image_mutex);
        for (auto fj = max(0, pxl.j - filter_size);
             fj <= min(img.height() - 1, pxl.j + filter_size); fj++) {
            for (auto fi = max(0, pxl.i - filter_size);
                 fi <= min(img.width() - 1, pxl.i + filter_size); fi++) {
                auto w = filter((fi - pxl.i) - uv.x + 0.5f) *
                         filter((fj - pxl.j) - uv.y + 0.5f);
                pxl.col += l * w;
                pxl.alpha += w;
                pxl.weight += w;
            }
        }
    }
}

// Trace the next nsamples.
void trace_samples_filtered(const scene* scn, const camera* cam,
    const bvh_tree* bvh, const trace_lights& lights, image4f& img,
    image<trace_pixel>& pixels, int nsamples, const trace_params& params) {
    auto shader = trace_shaders.at(params.shader);
    auto filter = trace_filters.at(params.filter);
    auto filter_size = trace_filter_sizes.at(params.filter);
    std::mutex image_mutex;
    if (params.parallel) {
        auto nthreads = std::thread::hardware_concurrency();
        auto threads = std::vector<std::thread>();
        for (auto tid = 0; tid < std::thread::hardware_concurrency(); tid++) {
            threads.push_back(
                std::thread([=, &pixels, &params, &image_mutex, &img]() {
                    for (auto j = tid; j < img.height(); j += nthreads) {
                        for (auto i = 0; i < img.width(); i++) {
                            auto& pxl = pixels.at(i, j);
                            for (auto s = 0; s < nsamples; s++) {
                                trace_sample_filtered(scn, cam, bvh, lights,
                                    img, pxl, shader, filter, filter_size,
                                    image_mutex, params);
                            }
                        }
                    }
                }));
        }
        for (auto& t : threads) t.join();
        threads.clear();
    } else {
        for (auto j = 0; j < img.height(); j++) {
            for (auto i = 0; i < img.width(); i++) {
                auto& pxl = pixels.at(i, j);
                for (auto s = 0; s < params.nsamples; s++) {
                    trace_sample_filtered(scn, cam, bvh, lights, img, pxl,
                        shader, filter, filter_size, image_mutex, params);
                }
            }
        }
    }
    for (auto j = 0; j < img.height(); j++) {
        for (auto i = 0; i < img.width(); i++) {
            auto& pxl = pixels.at(i, j);
            img.at(i, j) = {pxl.col.x, pxl.col.y, pxl.col.z, pxl.alpha};
            img.at(i, j) /= pxl.weight;
        }
    }
}

// Starts an anyncrhounous renderer.
void trace_async_start(const scene* scn, const camera* cam, const bvh_tree* bvh,
    const trace_lights& lights, image4f& img, image<trace_pixel>& pixels,
    std::vector<std::thread>& threads, bool& stop_flag,
    const trace_params& params) {
    pixels = make_trace_pixels(img, params);
    auto nthreads = std::thread::hardware_concurrency();
    for (auto tid = 0; tid < std::thread::hardware_concurrency(); tid++) {
        threads.push_back(std::thread([=, &img, &pixels, &stop_flag]() {
            auto shader = trace_shaders.at(params.shader);
            for (auto s = 0; s < params.nsamples; s++) {
                for (auto j = tid; j < img.height(); j += nthreads) {
                    for (auto i = 0; i < img.width(); i++) {
                        if (stop_flag) return;
                        auto& pxl = pixels.at(i, j);
                        trace_sample(
                            scn, cam, bvh, lights, pxl, shader, params);
                        img.at(i, j) = {
                            pxl.col.x, pxl.col.y, pxl.col.z, pxl.alpha};
                        img.at(i, j) /= pxl.sample;
                    }
                }
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

// Initialize trace lights
trace_lights make_trace_lights(const scene* scn) {
    auto lights = trace_lights();
    for (auto ist : scn->instances) {
        auto shp = ist->shp->shapes.at(0);
        if (!shp->mat) continue;
        if (shp->mat->ke == zero3f) continue;
        auto lgt = trace_light();
        lgt.ist = ist;
        lights.lights.push_back(lgt);
        if (!contains(lights.shape_cdfs, shp)) {
            if (!shp->points.empty()) {
                lights.shape_cdfs[shp] = sample_points_cdf(shp->points.size());
            } else if (!shp->lines.empty()) {
                lights.shape_cdfs[shp] = sample_lines_cdf(shp->lines, shp->pos);
            } else if (!shp->triangles.empty()) {
                lights.shape_cdfs[shp] =
                    sample_triangles_cdf(shp->triangles, shp->pos);
            }
            lights.shape_areas[shp] = lights.shape_cdfs[shp].back();
        }
    }

    for (auto env : scn->environments) {
        if (env->ke == zero3f) continue;
        auto lgt = trace_light();
        lgt.env = env;
        lights.lights.push_back(lgt);
    }

    return lights;
}

// Initialize a rendering state
image<trace_pixel> make_trace_pixels(
    const image4f& img, const trace_params& params) {
    auto pixels = image<trace_pixel>(img.width(), img.height());
    for (auto j = 0; j < img.height(); j++) {
        for (auto i = 0; i < img.width(); i++) {
            auto pxl = trace_pixel();
            pxl.i = i;
            pxl.j = j;
            pxl.rng = init_rng(params.seed, (j * img.width() + i) * 2 + 1);
            pixels.at(i, j) = pxl;
        }
    }
    return pixels;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR WAVEFRONT OBJ
// -----------------------------------------------------------------------------
namespace ygl {

// cleanup
obj_object::~obj_object() {
    for (auto v : groups) delete v;
}

// cleanup
obj_scene::~obj_scene() {
    for (auto v : objects)
        if (v) delete v;
    for (auto v : materials)
        if (v) delete v;
    for (auto v : textures)
        if (v) delete v;
    for (auto v : cameras)
        if (v) delete v;
    for (auto v : environments)
        if (v) delete v;
    for (auto v : nodes)
        if (v) delete v;
}

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

// parse value
inline void obj_parse(char*& s, float& val) {
    obj_skipws(s);
    obj_parse_base(s, val);
}

// parse value
inline void obj_parse(char*& s, bool& val) {
    auto i = 0;
    obj_parse(s, i);
    val = i;
}

// parse value
inline void obj_parse(char*& s, std::string& val) {
    obj_skipws(s);
    obj_parse_base(s, val);
}

// parse value
inline void obj_parse(char*& s, char* val) {
    obj_skipws(s);
    obj_parse_base(s, val);
}

// parse value
template <typename T, int N>
inline void obj_parse(char*& s, vec<T, N>& val) {
    for (auto i = 0; i < N; i++) obj_parse(s, (&val.x)[i]);
}

// parse value
template <typename T, int N>
inline void obj_parse(char*& s, quat<T, N>& val) {
    for (auto i = 0; i < N; i++) obj_parse(s, (&val.x)[i]);
}

// parse value
template <typename T, int N>
inline void obj_parse(char*& s, frame<T, N>& val) {
    for (auto i = 0; i < N * (N + 1); i++) obj_parse(s, (&val.x.x)[i]);
}

// parse value
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
        } else if (obj_streq(cmd, "Kt") || obj_streq(cmd, "Tf")) {
            auto ntok = 0;
            obj_skipws(ss);
            while (*ss && ntok < 3) {
                obj_parse(ss, mat->kt[ntok++]);
                obj_skipws(ss);
            }
            if (ntok < 3) mat->kt = {mat->kt.x, mat->kt.x, mat->kt.x};
        } else if (obj_streq(cmd, "Tr")) {
            auto ntok = 0;
            obj_skipws(ss);
            while (*ss && ntok < 3) {
                obj_parse(ss, mat->kt[ntok++]);
                obj_skipws(ss);
            }
            if (ntok < 3) {
                materials.back()->op = (flip_tr) ? 1 - mat->kt.x : mat->kt.x;
                mat->kt = {0, 0, 0};
            }
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
    obj_scene* asset, const std::string& dirname, bool skip_missing) {
    for (auto txt : asset->textures) {
        auto filename = dirname + txt->path;
        for (auto& c : filename)
            if (c == '\\') c = '/';
#if YGL_IMAGEIO
        if (is_hdr_filename(filename)) {
            txt->dataf =
                load_imagef(filename, txt->width, txt->height, txt->ncomp);
        } else {
            txt->datab =
                load_image(filename, txt->width, txt->height, txt->ncomp);
        }
#endif
        if (txt->datab.empty() && txt->dataf.empty()) {
            if (skip_missing) continue;
            throw std::runtime_error("cannot laod image " + filename);
        }
    }
}

// Loads an OBJ
obj_scene* load_obj(const std::string& filename, bool load_txt,
    bool skip_missing, bool flip_texcoord, bool flip_tr) {
    // clear obj
    auto asset = std::unique_ptr<obj_scene>(new obj_scene());

    // open file
    auto fs = fopen(filename.c_str(), "rt");
    if (!fs) throw std::runtime_error("cannot open filename " + filename);

    // initializing obj
    asset->objects.push_back(new obj_object());
    asset->objects.back()->groups.push_back(new obj_group());

    // current parsing value
    auto matname = std::string();
    auto mtllibs = std::vector<std::string>();
    auto object = asset->objects.back();
    auto group = object->groups.back();
    auto elems = std::vector<obj_vertex>();

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
            asset->pos.push_back(zero3f);
            obj_parse(ss, asset->pos.back());
        } else if (obj_streq(cmd, "vn")) {
            vert_size.norm += 1;
            asset->norm.push_back(zero3f);
            obj_parse(ss, asset->norm.back());
        } else if (obj_streq(cmd, "vt")) {
            vert_size.texcoord += 1;
            asset->texcoord.push_back(zero2f);
            obj_parse(ss, asset->texcoord.back());
            if (flip_texcoord)
                asset->texcoord.back().y = 1 - asset->texcoord.back().y;
        } else if (obj_streq(cmd, "vc")) {
            vert_size.color += 1;
            asset->color.push_back(vec4f{0, 0, 0, 1});
            obj_parse(ss, asset->color.back());
        } else if (obj_streq(cmd, "vr")) {
            vert_size.radius += 1;
            asset->radius.push_back(0);
            obj_parse(ss, asset->radius.back());
        } else if (obj_streq(cmd, "f") || obj_streq(cmd, "l") ||
                   obj_streq(cmd, "p") || obj_streq(cmd, "b")) {
            group->elems.push_back(
                {(uint32_t)group->verts.size(), elem_type_map.at(cmd), 0});
            obj_skipws(ss);
            while (*ss) {
                auto vert = obj_vertex();
                obj_parse(ss, vert, vert_size);
                obj_skipws(ss);
                group->verts.push_back(vert);
                group->elems.back().size += 1;
            }
        } else if (obj_streq(cmd, "o")) {
            asset->objects.push_back(new obj_object());
            object = asset->objects.back();
            obj_parse(ss, object->name);
            object->groups.push_back(new obj_group());
            group = object->groups.back();
            group->matname = matname;
        } else if (obj_streq(cmd, "usemtl")) {
            obj_parse(ss, matname);
            object->groups.push_back(new obj_group());
            group = object->groups.back();
            group->matname = matname;
        } else if (obj_streq(cmd, "g")) {
            object->groups.push_back(new obj_group());
            group = object->groups.back();
            obj_parse(ss, group->groupname);
            group->matname = matname;
        } else if (obj_streq(cmd, "s")) {
            auto name = std::string();
            obj_parse(ss, name);
            auto smoothing = (name == "on");
            if (group->smoothing != smoothing) {
                auto gname = group->groupname;
                object->groups.push_back(new obj_group());
                group = object->groups.back();
                group->matname = matname;
                group->groupname = gname;
                group->smoothing = smoothing;
            }
        } else if (obj_streq(cmd, "gp")) {
            auto name = std::string();
            obj_parse(ss, name);
            obj_skipws(ss);
            while (*ss) {
                auto tok = std::string();
                obj_parse(ss, tok);
                obj_skipws(ss);
                group->props[name].push_back(tok);
            }
        } else if (obj_streq(cmd, "op")) {
            auto name = std::string();
            obj_parse(ss, name);
            obj_skipws(ss);
            while (*ss) {
                auto tok = std::string();
                obj_parse(ss, tok);
                obj_skipws(ss);
                object->props[name].push_back(tok);
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
            asset->cameras.push_back(cam);
        } else if (obj_streq(cmd, "e")) {
            auto env = new obj_environment();
            obj_parse(ss, env->name);
            obj_parse(ss, env->matname);
            obj_parse(ss, env->frame);
            asset->environments.push_back(env);
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
            obj_parse(ss, nde->scaling);
            if (nde->parent == "\"\"") nde->parent = "";
            if (nde->camname == "\"\"") nde->camname = "";
            if (nde->objname == "\"\"") nde->objname = "";
            if (nde->envname == "\"\"") nde->envname = "";
            asset->nodes.push_back(nde);
        } else {
            // unused
        }
    }

    // cleanup unused
    for (auto& o : asset->objects) {
        for (auto& g : o->groups) {
            if (g->verts.empty()) {
                delete g;
                g = nullptr;
            }
        }
        auto end = std::remove_if(o->groups.begin(), o->groups.end(),
            [](const obj_group* x) { return !x; });
        o->groups.erase(end, o->groups.end());
        if (o->groups.empty()) {
            delete o;
            o = nullptr;
        }
    }
    auto end = std::remove_if(asset->objects.begin(), asset->objects.end(),
        [](const obj_object* x) { return !x; });
    asset->objects.erase(end, asset->objects.end());

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
        asset->materials.insert(
            asset->materials.end(), materials.begin(), materials.end());
        for (auto& txt : textures) {
            if (texture_set.find(txt) != texture_set.end()) continue;
            asset->textures.push_back(new obj_texture());
            asset->textures.back()->path = txt;
            texture_set.insert(txt);
        }
    }

    // load textures
    if (load_txt) load_textures(asset.get(), dirname, skip_missing);

    // close file
    fclose(fs);

    // done
    return asset.release();
}

// Dumps a value
inline void obj_dump(char*& s, char* val) {
    while (*val) *s++ = *val++;
}
// Dumps a value
inline void obj_dump(char*& s, const char* val) {
    while (*val) *s++ = *val++;
}
// Dumps a value
inline void obj_dump(char*& s, const std::string& val) {
    auto val_ = val.c_str();
    while (*val_) *s++ = *val_++;
}

#if YGL_FASTOBJ
// Dumps a value
inline void obj_dump(char*& s, int val) {
    static auto digits = "0123456789";
    if (val < 0) {
        *s++ = '-';
        val = -val;
    }
    char buf[64];
    buf[sizeof(buf) - 1] = 0;
    auto ss = buf + sizeof(buf) - 2;
    while (val >= 10) {
        *--ss = digits[val % 10];
        val /= 10;
    }
    *--ss = digits[val];
    while (*ss) *s++ = *ss++;
}
#else
// Dumps a value
inline void obj_dump(char*& s, int val) { s += sprintf(s, "%d", val); }
#endif

// Dumps a value
inline void obj_dump(char*& s, float val) { s += sprintf(s, "%g", val); }
// Dumps a value
template <typename T, int N>
inline void obj_dump(char*& s, const vec<T, N>& val) {
    for (auto i = 0; i < N; i++) {
        if (i) *s++ = ' ';
        obj_dump(s, (&val.x)[i]);
    }
}
// Dumps a value
template <typename T, int N>
inline void obj_dump(char*& s, const quat<T, N>& val) {
    for (auto i = 0; i < N; i++) {
        if (i) *s++ = ' ';
        obj_dump(s, (&val.x)[i]);
    }
}
// Dumps a value
template <typename T, int N>
inline void obj_dump(char*& s, const frame<T, N>& val) {
    for (auto i = 0; i < N * (N + 1); i++) {
        if (i) *s++ = ' ';
        obj_dump(s, (&val.x.x)[i]);
    }
}

// Dumps a value
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

// Dumps a value
inline void obj_dump(char*& s, const std::vector<std::string>& vals) {
    for (auto i = 0; i < vals.size(); i++) {
        if (i) *s++ = ' ';
        obj_dump(s, vals[i]);
    }
}

// Dumps a value
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
        if (mat->kt != zero3f) obj_dump_line(fs, "  Tf", mat->kt);
        if (mat->ns != 0.0f) obj_dump_line(fs, "  Ns", mat->ns);
        if (mat->op != 1.0f) obj_dump_line(fs, "  d", mat->op);
        if (mat->ior != 0.0f) obj_dump_line(fs, "  Ni", mat->ior);
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
    const obj_scene* asset, const std::string& dirname, bool skip_missing) {
    for (auto txt : asset->textures) {
        if (txt->datab.empty() && txt->dataf.empty()) continue;
        auto filename = dirname + txt->path;
        for (auto& c : filename)
            if (c == '\\') c = '/';
        auto ok = false;
#if YGL_IMAGEIO
        if (!txt->datab.empty()) {
            ok = save_image(filename, txt->width, txt->height, txt->ncomp,
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
void save_obj(const std::string& filename, const obj_scene* asset,
    bool save_txt, bool skip_missing, bool flip_texcoord, bool flip_tr) {
    // open file
    auto fs = fopen(filename.c_str(), "wt");
    if (!fs) throw std::runtime_error("cannot open filename " + filename);

    // linkup to mtl
    auto dirname = path_dirname(filename);
    auto basename = filename.substr(dirname.length());
    basename = basename.substr(0, basename.length() - 4);
    if (!asset->materials.empty()) {
        obj_dump_line(fs, "mtllib", basename + ".mtl");
    }

    // save cameras
    for (auto cam : asset->cameras) {
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
    for (auto env : asset->environments) {
        obj_dump_sp(fs, "e");
        obj_dump_sp(fs, env->name);
        obj_dump_sp(fs, env->matname);
        obj_dump_sp(fs, env->frame);
        obj_dump_sp(fs, "\n");
    }

    // save nodes
    for (auto nde : asset->nodes) {
        obj_dump_sp(fs, "n");
        obj_dump_sp(fs, nde->name);
        obj_dump_sp(fs, (nde->parent.empty()) ? "\"\""s : nde->parent);
        obj_dump_sp(fs, (nde->camname.empty()) ? "\"\""s : nde->camname);
        obj_dump_sp(fs, (nde->objname.empty()) ? "\"\""s : nde->objname);
        obj_dump_sp(fs, (nde->envname.empty()) ? "\"\""s : nde->envname);
        obj_dump_sp(fs, nde->frame);
        obj_dump_sp(fs, nde->translation);
        obj_dump_sp(fs, nde->rotation);
        obj_dump_sp(fs, nde->scaling);
        obj_dump_sp(fs, "\n");
    }

    // save all vertex data
    for (auto& v : asset->pos) obj_dump_line(fs, "v", v);
    if (flip_texcoord) {
        for (auto& v : asset->texcoord)
            obj_dump_line(fs, "vt", vec2f{v.x, 1 - v.y});
    } else {
        for (auto& v : asset->texcoord) obj_dump_line(fs, "vt", v);
    }
    for (auto& v : asset->norm) obj_dump_line(fs, "vn", v);
    for (auto& v : asset->color) obj_dump_line(fs, "vc", v);
    for (auto& v : asset->radius) obj_dump_line(fs, "vr", v);

    // save element data
    static auto elem_labels = std::unordered_map<obj_element_type, std::string>{
        {obj_element_type::point, "p"}, {obj_element_type::line, "l"},
        {obj_element_type::face, "f"}, {obj_element_type::bezier, "b"}};
    for (auto object : asset->objects) {
        obj_dump_line(fs, "o", object->name);
        for (auto& kv : object->props) {
            obj_dump_line(fs, "op", join({kv.first}, kv.second));
        }
        for (auto group : object->groups) {
            if (!group->matname.empty())
                obj_dump_line(fs, "usemtl", group->matname);
            if (!group->groupname.empty())
                obj_dump_line(fs, "g", group->groupname);
            if (!group->smoothing) obj_dump_line(fs, "s", "off");
            for (auto& kv : group->props) {
                obj_dump_line(fs, "gp", join({kv.first}, kv.second));
            }
            for (auto elem : group->elems) {
                obj_dump_line(fs, elem_labels.at(elem.type).c_str(),
                    group->verts.data() + elem.start, elem.size);
            }
        }
    }

    fclose(fs);

    // save materials
    if (!asset->materials.empty())
        save_mtl(dirname + basename + ".mtl", asset->materials, flip_tr);

    // save textures
    if (save_txt) save_textures(asset, dirname, skip_missing);
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
template <typename T, int N>
void serialize(vec<T, N>& vals, json& js, bool reading) {
    serialize((std::array<T, N>&)vals, js, reading);
}

// Parse support function.
template <typename T, int N>
void serialize(quat<T, N>& vals, json& js, bool reading) {
    serialize((std::array<T, N>&)vals, js, reading);
}

// Parse support function.
template <typename T, int N>
void serialize(mat<T, N>& vals, json& js, bool reading) {
    serialize((std::array<T, N * N>&)vals, js, reading);
}

// Parse support function.
template <typename T, int N>
void serialize(frame<T, N>& vals, json& js, bool reading) {
    serialize((std::array<T, N*(N + 1)>&)vals, js, reading);
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

// #codegen begin gltf-func
// ----------------------------------------------------

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
            {"CATMULLROMSPLINE",
                glTFAnimationSamplerInterpolation::CatmullRomSpline},
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
void load_buffers(glTF* gltf, const std::string& dirname, bool skip_missing) {
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
                buffer->data =
                    load_binary(path_convert_eparator(dirname + buffer->uri));
                if (buffer->data.empty()) {
                    if (skip_missing) continue;
                    throw std::runtime_error(
                        "could not load binary file " +
                        path_convert_eparator(dirname + buffer->uri));
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
void load_images(glTF* gltf, const std::string& dirname, bool skip_missing) {
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
                image->data.datab = load_image_from_memory(filename, data,
                    data_size, image->data.width, image->data.height,
                    image->data.ncomp);
            }
        } else {
            filename = path_convert_eparator(dirname + image->uri);
            if (is_hdr_filename(filename)) {
                image->data.dataf = load_imagef(filename, image->data.width,
                    image->data.height, image->data.ncomp);
            } else {
                image->data.datab = load_image(filename, image->data.width,
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
    auto gltf = std::unique_ptr<glTF>(new glTF());

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
    auto gltf_ = gltf.get();
    try {
        serialize(gltf_, js, true);
    } catch (const std::exception& e) {
        throw std::runtime_error("error parsing gltf " + std::string(e.what()));
    }

    // load external resources
    auto dirname = path_dirname(filename);
    if (load_bin) load_buffers(gltf.get(), dirname, skip_missing);
    if (load_image) load_images(gltf.get(), dirname, skip_missing);

    // done
    return gltf.release();
}

// Save buffer data.
void save_buffers(
    const glTF* gltf, const std::string& dirname, bool skip_missing) {
    for (auto buffer : gltf->buffers) {
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
                ok = save_image(filename, image->data.width, image->data.height,
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
    auto gltf = std::unique_ptr<glTF>(new glTF());

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
    auto gltf_ = gltf.get();
    try {
        serialize(gltf_, js, true);
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
    if (load_bin) load_buffers(gltf.get(), dirname, skip_missing);
    if (load_image) load_images(gltf.get(), dirname, skip_missing);

    // close
    fclose(f);

    // done
    return gltf.release();
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

#if YGL_SVG

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR SVG
// -----------------------------------------------------------------------------
namespace ygl {

// Load SVG
svg_scene* load_svg(const std::string& filename) {
    auto svg = nsvgParseFromFile(filename.c_str(), "mm", 96);
    if (!svg) throw std::runtime_error("cannot load SVG");
    auto scn = new svg_scene();
    for (auto shape = svg->shapes; shape != nullptr; shape = shape->next) {
        auto shp = new svg_shape();
        scn->shapes.push_back(shp);
        for (auto path = shape->paths; path != nullptr; path = path->next) {
            auto pth = new svg_path();
            shp->paths.push_back(pth);
            pth->pos.resize(path->npts);
            for (int i = 0; i < path->npts; i += 1) {
                pth->pos[i] = {path->pts[i * 2 + 0], path->pts[i * 2 + 1]};
            }
        }
    }
    nsvgDelete(svg);
    return scn;
}

// Save SVG
void save_svg(const std::string& filename, const std::vector<svg_path>& paths) {
    throw std::runtime_error("not implemented yet");
}

}  // namespace ygl

#endif

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR SHAPE EXAMPLES
// -----------------------------------------------------------------------------
namespace ygl {

// Make a sphere. This is not watertight.
void make_uvsphere(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, int tesselation,
    bool flipped) {
    if (!flipped) {
        make_quads(quads, pos, norm, texcoord, pow2(tesselation + 2),
            pow2(tesselation + 1),
            [](auto uv, auto& pos, auto& norm, auto& texcoord) {
                auto a = vec2f{2 * pif * uv.x, pif * (1 - uv.y)};
                pos = {cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
                norm = {cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
                texcoord = uv;
            });
    } else {
        make_quads(quads, pos, norm, texcoord, pow2(tesselation + 2),
            pow2(tesselation + 1),
            [](auto uv, auto& pos, auto& norm, auto& texcoord) {
                auto a = vec2f{2 * pif * uv.x, pif * uv.y};
                pos = {cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
                norm = {-cos(a.x) * sin(a.y), -sin(a.x) * sin(a.y), -cos(a.y)};
                texcoord = {uv.x, 1 - uv.y};
            });
    }
}

// Make a geodesic sphere.
void make_geodesicsphere(
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

// Make a sphere. This is not watertight.
void make_uvhemisphere(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, int tesselation,
    bool flipped) {
    if (!flipped) {
        make_quads(quads, pos, norm, texcoord, pow2(tesselation + 2),
            pow2(tesselation),
            [](auto uv, auto& pos, auto& norm, auto& texcoord) {
                auto a = vec2f{2 * pif * uv.x, pif * 0.5f * (1 - uv.y)};
                pos = {cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
                norm = {cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
                texcoord = uv;
            });
    } else {
        make_quads(quads, pos, norm, texcoord, pow2(tesselation + 2),
            pow2(tesselation),
            [](auto uv, auto& pos, auto& norm, auto& texcoord) {
                auto a = vec2f{2 * pif * uv.x, pif * (0.5f + 0.5f * uv.y)};
                pos = {cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
                norm = {-cos(a.x) * sin(a.y), -sin(a.x) * sin(a.y), -cos(a.y)};
                texcoord = {uv.x, 1 - uv.y};
            });
    }
}

// Make a quad.
void make_uvquad(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, int tesselation) {
    make_quads(quads, pos, norm, texcoord, pow2(tesselation), pow2(tesselation),
        [](auto uv, auto& pos, auto& norm, auto& texcoord) {
            pos = vec3f{(-1 + uv.x * 2), (-1 + uv.y * 2), 0};
            norm = vec3f{0, 0, 1};
            texcoord = uv;
        });
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
    make_quads(quads_pos, pos, pow2(tesselation + 2), pow2(tesselation + 1),
        [](auto uv) {
            auto a = vec2f{2 * pif * uv.x, pif * (1 - uv.y)};
            return vec3f{cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
        },
        true, false, true, true);
    make_quads(quads_norm, norm, pow2(tesselation + 2), pow2(tesselation + 1),
        [](auto uv) {
            auto a = vec2f{2 * pif * uv.x, pif * (1 - uv.y)};
            return vec3f{cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
        },
        true, false, true, true);
    make_quads(quads_texcoord, texcoord, pow2(tesselation + 2),
        pow2(tesselation + 1), [](auto uv) { return uv; });
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

// Make a cube with uv. This is not watertight.
void make_uvcube(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, int tesselation) {
    frame3f frames[6] = {frame3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 1}},
        frame3f{{-1, 0, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, -1}},
        frame3f{{-1, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 0}},
        frame3f{{1, 0, 0}, {0, 0, 1}, {0, -1, 0}, {0, -1, 0}},
        frame3f{{0, 1, 0}, {0, 0, 1}, {1, 0, 0}, {1, 0, 0}},
        frame3f{{0, -1, 0}, {0, 0, 1}, {-1, 0, 0}, {-1, 0, 0}}};
    std::vector<vec3f> quad_pos, quad_norm;
    std::vector<vec2f> quad_texcoord;
    std::vector<vec4i> quad_quads;
    make_uvquad(quad_quads, quad_pos, quad_norm, quad_texcoord, tesselation);
    for (auto i = 0; i < 6; i++) {
        pos.insert(pos.end(), quad_pos.begin(), quad_pos.end());
        norm.insert(norm.end(), quad_norm.begin(), quad_norm.end());
        texcoord.insert(
            texcoord.end(), quad_texcoord.begin(), quad_texcoord.end());
        quads.insert(quads.end(), quad_quads.begin(), quad_quads.end());
    }
    auto quad_verts = quad_pos.size();
    for (auto i = 0; i < 6; i++) {
        for (auto j = quad_verts * i; j < quad_verts * (i + 1); j++)
            pos[j] = transform_point(frames[i], pos[j]);
        for (auto j = quad_verts * i; j < quad_verts * (i + 1); j++)
            norm[j] = transform_direction(frames[i], norm[j]);
    }
    auto quad_faces = quad_quads.size();
    for (auto i = 0; i < 6; i++) {
        for (auto j = quad_faces * i; j < quad_faces * (i + 1); j++) {
            quads[j].x += quad_verts * i;
            quads[j].y += quad_verts * i;
            quads[j].z += quad_verts * i;
            quads[j].w += quad_verts * i;
        }
    }
}

// Make a sphere from a cube. This is not watertight.
void make_uvspherecube(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, int tesselation) {
    make_uvcube(quads, pos, norm, texcoord, tesselation);
    for (auto i = 0; i < pos.size(); i++) {
        pos[i] = normalize(pos[i]);
        norm[i] = normalize(pos[i]);
    }
}

// Make a cube than stretch it towards a sphere. This is not watertight.
void make_uvspherizedcube(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, int tesselation,
    float radius) {
    make_uvcube(quads, pos, norm, texcoord, tesselation);
    for (auto i = 0; i < pos.size(); i++) {
        norm[i] = normalize(pos[i]);
        pos[i] *= 1 - radius;
        pos[i] += norm[i] * radius;
    }
    compute_normals(quads, pos, norm);
}

// Make a flipped sphere. This is not watertight.
void make_uvflipcapsphere(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, int tesselation,
    float z, bool flipped) {
    make_uvsphere(quads, pos, norm, texcoord, tesselation, flipped);
    for (auto i = 0; i < pos.size(); i++) {
        if (pos[i].z > z) {
            pos[i].z = 2 * z - pos[i].z;
            norm[i].x = -norm[i].x;
            norm[i].y = -norm[i].y;
        } else if (pos[i].z < -z) {
            pos[i].z = -2 * z - pos[i].z;
            norm[i].x = -norm[i].x;
            norm[i].y = -norm[i].y;
        }
    }
}

// Make a cutout sphere. This is not watertight.
void make_uvcutsphere(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, int tesselation,
    float z, bool flipped) {
    if (!flipped) {
        make_quads(quads, pos, norm, texcoord, pow2(tesselation + 2),
            pow2(tesselation + 1),
            [=](auto uv, auto& pos, auto& norm, auto& texcoord) {
                auto p = 1 - acos(z) / pif;
                auto a = vec2f{2 * pif * uv.x, pif * (1 - p * uv.y)};
                pos = {cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
                norm = {cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
                texcoord = uv;
            });
    } else {
        make_quads(quads, pos, norm, texcoord, pow2(tesselation + 2),
            pow2(tesselation + 1),
            [=](auto uv, auto& pos, auto& norm, auto& texcoord) {
                auto p = 1 - acos(z) / pif;
                auto a = vec2f{2 * pif * uv.x, pif * ((1 - p) + p * uv.y)};
                pos = {cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
                norm = {-cos(a.x) * sin(a.y), -sin(a.x) * sin(a.y), -cos(a.y)};
                texcoord = {uv.x, 1 - uv.y};
            });
    }
}

// Make a seashell. This is not watertight. Returns quads, pos, norm,
// texcoord.
void make_uvseashell(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, int tesselation,
    const make_seashell_params& params) {
    auto R = params.spiral_revolutions;
    auto D = -1.0f;
    auto a = params.spiral_angle;
    auto b = params.enlarging_angle;
    auto A = params.spiral_aperture;
    auto e = params.ellipse_axis;
    auto O = params.curve_rotation;  // (psi, Omega, mu)
    auto W = params.nodule_length;
    auto N = params.nodules_num;
    auto P = params.nodule_pos;
    auto L = params.nodule_height;

    auto cot_a = 1 / tan(a);

    make_quads(quads, pos, norm, texcoord, pow2(tesselation + 2),
        pow2(tesselation + 1 + (int)round(R)),
        [=](auto uv, auto& pos, auto& norm, auto& texcoord) {
            auto s = uv.x * 2 * pif;
            auto t = uv.y * 2 * pif * R - pif * R;
            auto re = 1 / sqrt(pow(cos(s) / e.x, 2) + pow(sin(s) / e.y, 2));
            if (L && W.x && W.y && t > 0) {
                auto l = (2 * pif / N) *
                         ((N * t) / (2 * pif) - floor((N * t) / (2 * pif)));
                auto rn = L * exp(-pow((2 * (s - P)) / W.x, 2) -
                                  pow((2 * l) / W.y, 2));
                re += rn;
            }
            pos.x = (A * sin(b) * cos(t) + cos(s + O.x) * cos(t + O.y) * re -
                        sin(O.z) * sin(t + O.y) * re) *
                    D * exp(t * cot_a);
            pos.y = (A * sin(b) * sin(t) + cos(s + O.x) * sin(t + O.y) * re +
                        sin(O.z) * sin(s + O.x) * cos(t + O.y) * re) *
                    exp(t * cot_a);
            pos.z =
                (-A * cos(b) + cos(O.z) * sin(s + O.x) * re) * exp(t * cot_a);
            norm = vec3f{0, 0, 1};
            texcoord = vec2f{uv.x, uv.y * R};
        });

    norm = std::vector<vec3f>(pos.size());
    compute_normals(quads, pos, norm);
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
    std::vector<float>& radius, int num, int tesselation,
    const std::vector<vec3i>& striangles, const std::vector<vec4i>& squads,
    const std::vector<vec3f>& spos, const std::vector<vec3f>& snorm,
    const std::vector<vec2f>& stexcoord, const make_hair_params& params) {
    std::vector<vec3f> bpos;
    std::vector<vec3f> bnorm;
    std::vector<vec2f> btexcoord;
    std::tie(bpos, bnorm, btexcoord) = sample_triangles_points(
        join(striangles, convert_quads_to_triangles(squads)), spos, snorm,
        stexcoord, num, params.seed);

    auto rng = init_rng(params.seed, 3);
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

    auto usteps = pow2(tesselation);
    make_lines_uv(lines, texcoord, num, usteps);
    pos = std::vector<vec3f>(texcoord.size());
    norm = std::vector<vec3f>(texcoord.size());
    radius = std::vector<float>(texcoord.size());
    for (auto i = 0; i < texcoord.size(); i++) {
        auto u = texcoord[i].x;
        auto bidx = i / (usteps + 1);
        pos[i] = bpos[bidx] + bnorm[bidx] * u * blen[bidx];
        norm[i] = bnorm[bidx];
        radius[i] = lerp(params.radius.x, params.radius.y, u);
        if (params.clump.x > 0) {
            pos[i] = lerp(pos[i], pos[i + (cidx[bidx] - bidx) * (usteps + 1)],
                u * params.clump.x);
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
// IMPLEMENTATION FOR EXAMPLE SCENES
// -----------------------------------------------------------------------------
namespace ygl {

// makes the cornell box scene
// http://graphics.cs.williams.edu/data
// http://www.graphics.cornell.edu/online/box/data.html
scene* make_cornell_box_scene() {
    auto make_camera = [](std::string name, vec3f from, vec3f to, float yfov,
                           float aperture, float aspect = 16.0f / 9.0f) {
        auto cam = new camera();
        cam->name = name;
        cam->frame = lookat_frame(from, to, {0, 1, 0});
        cam->aperture = aperture;
        cam->focus = length(from - to);
        cam->yfov = yfov * pif / 180;
        cam->aspect = aspect;
        return cam;
    };

    auto make_instance = [](std::string name, shape_group* shp, vec3f pos,
                             vec3f rot = {0, 0, 0}) {
        auto ist = new instance();
        ist->name = name;
        ist->shp = shp;
        ist->frame = translation_frame(pos) *
                     rotation_frame(vec3f{0, 0, 1}, rot[2] * pif / 180) *
                     rotation_frame(vec3f{0, 1, 0}, rot[1] * pif / 180) *
                     rotation_frame(vec3f{1, 0, 0}, rot[0] * pif / 180);
        return ist;
    };

    auto make_quad = [](std::string name, material* mat, float scale = 1) {
        auto sgr = new shape_group();
        sgr->name = name;
        auto shp = new shape();
        shp->mat = mat;
        shp->name = name;
        make_uvquad(shp->quads, shp->pos, shp->norm, shp->texcoord, 0);
        for (auto& p : shp->pos) p *= scale;
        sgr->shapes.push_back(shp);
        return sgr;
    };

    auto make_box = [](std::string name, material* mat, vec3f scale) {
        auto sgr = new shape_group();
        sgr->name = name;
        auto shp = new shape();
        shp->mat = mat;
        shp->name = name;
        make_uvcube(shp->quads, shp->pos, shp->norm, shp->texcoord, 0);
        for (auto& p : shp->pos) p *= scale;
        sgr->shapes.push_back(shp);
        return sgr;
    };

    auto make_material = [](std::string name, vec3f kd, vec3f ke = {0, 0, 0}) {
        auto mat = new material();
        mat->type = material_type::specular_roughness;
        mat->name = name;
        mat->ke = ke;
        mat->kd = kd;
        mat->ks = zero3f;
        mat->rs = 1;
        return mat;
    };

    auto scn = new scene();
    scn->cameras.push_back(
        make_camera("cb_cam", {0, 1, 5.15f}, {0, 1, 0}, 27, 0, 1));
    scn->materials.push_back(make_material("cb_white", {0.725f, 0.71f, 0.68f}));
    scn->materials.push_back(make_material("cb_red", {0.63f, 0.065f, 0.05f}));
    scn->materials.push_back(make_material("cb_green", {0.14f, 0.45f, 0.091f}));
    scn->materials.push_back(make_material("cb_light", zero3f, {17, 12, 4}));
    scn->shapes.push_back(make_quad("cb_floor", scn->materials[0]));
    scn->shapes.push_back(make_quad("cb_ceiling", scn->materials[0]));
    scn->shapes.push_back(make_quad("cb_back", scn->materials[0]));
    scn->shapes.push_back(make_quad("cb_left", scn->materials[2]));
    scn->shapes.push_back(make_quad("cb_right", scn->materials[1]));
    scn->shapes.push_back(
        make_box("cb_tallbox", scn->materials[0], {0.3f, 0.6f, 0.3f}));
    scn->shapes.push_back(
        make_box("cb_shortbox", scn->materials[0], {0.3f, 0.3f, 0.3f}));
    scn->shapes.push_back(make_quad("cb_light", scn->materials[3], 0.25f));
    scn->instances.push_back(
        make_instance("cb_floor", scn->shapes[0], {0, 0, 0}, {-90, 0, 0}));
    scn->instances.push_back(
        make_instance("cb_ceiling", scn->shapes[1], {0, 2, 0}, {90, 0, 0}));
    scn->instances.push_back(
        make_instance("cb_back", scn->shapes[2], {0, 1, -1}));
    scn->instances.push_back(
        make_instance("cb_left", scn->shapes[3], {+1, 1, 0}, {0, -90, 0}));
    scn->instances.push_back(
        make_instance("cb_right", scn->shapes[4], {-1, 1, 0}, {0, 90, 0}));
    scn->instances.push_back(make_instance(
        "cb_tallbox", scn->shapes[5], {-0.33f, 0.6f, -0.29f}, {0, 15, 0}));
    scn->instances.push_back(make_instance(
        "cb_shortbox", scn->shapes[6], {0.33f, 0.3f, 0.33f}, {0, -15, 0}));
    scn->instances.push_back(
        make_instance("cb_light", scn->shapes[7], {0, 1.999f, 0}, {90, 0, 0}));
    return scn;
}

// Make standard shape.
void make_uvhollowcutsphere(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, int tesselation,
    float radius) {
    std::vector<vec3f> mpos, mnorm;
    std::vector<vec2f> mtexcoord;
    std::vector<vec4i> mquads;
    make_uvcutsphere(mquads, mpos, mnorm, mtexcoord, tesselation, radius);
    for (auto& uv : mtexcoord) uv.y *= radius;
    merge_quads(quads, pos, norm, texcoord, mquads, mpos, mnorm, mtexcoord);

    make_uvcutsphere(mquads, mpos, mnorm, mtexcoord, tesselation, radius, true);
    for (auto& p : mpos) p *= radius;
    merge_quads(quads, pos, norm, texcoord, mquads, mpos, mnorm, mtexcoord);

    // dpdu = [- s r s0 s1, s r c0 s1, 0] === [- s0, c0, 0]
    // dpdv = [s c0 s1, s s0 s1, s c1] === [c0 s1, s0 s1, c1]
    // n = [c0 c1, - s0 c1, s1]
    make_quads(mquads, mpos, mnorm, mtexcoord, pow2(tesselation + 2),
        pow2(tesselation + 1),
        [=](auto uv, auto& pos, auto& norm, auto& texcoord) {
            auto a = vec2f{2 * pif * uv[0], pif * (1 - radius)};
            auto r = (1 - uv[1]) + uv[1] * radius;
            pos = {r * cos(a[0]) * sin(a[1]), r * sin(a[0]) * sin(a[1]),
                r * cos(a[1])};
            norm = {-cos(a[0]) * cos(a[1]), -sin(a[0]) * cos(a[1]), sin(a[1])};
            texcoord = vec2f{uv[0], radius + (1 - radius) * uv[1]};
        });
    merge_quads(quads, pos, norm, texcoord, mquads, mpos, mnorm, mtexcoord);
}

// Make standard shape.
void make_uvhollowcutsphere1(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, int tesselation,
    float radius) {
    std::vector<vec3f> mpos, mnorm;
    std::vector<vec2f> mtexcoord;
    std::vector<vec4i> mquads;
    std::vector<vec2i> _aux1;
    std::vector<vec3i> _aux2;

    make_uvcutsphere(mquads, mpos, mnorm, mtexcoord, tesselation, radius);
    for (auto& uv : mtexcoord) uv.y *= radius;
    for (auto i = (pow2(tesselation + 2) + 1) * pow2(tesselation + 1);
         i < mnorm.size(); i++)
        mnorm[i] = normalize(mnorm[i] + vec3f{0, 0, 1});
    merge_quads(quads, pos, norm, texcoord, mquads, mpos, mnorm, mtexcoord);

    make_uvcutsphere(
        mquads, mpos, mnorm, mtexcoord, tesselation, radius * 1.05f, true);
    for (auto& p : mpos) p *= 0.8f;
    merge_quads(quads, pos, norm, texcoord, mquads, mpos, mnorm, mtexcoord);

    make_quads(mquads, mpos, mnorm, mtexcoord, pow2(tesselation + 2),
        pow2(tesselation + 1) / 4,
        [=](auto uv, auto& pos, auto& norm, auto& texcoord) {
            auto p = 1 - acos(radius) / pif;
            auto v = p + uv[1] * (1 - p);
            auto a = vec2f{2 * pif * uv[0], pif * (1 - v)};
            pos = vec3f{cos(a[0]) * sin(a[1]), sin(a[0]) * sin(a[1]),
                (2 * radius - cos(a[1]))};
            norm = vec3f{
                -cos(a[0]) * sin(a[1]), -sin(a[0]) * sin(a[1]), cos(a[1])};
            texcoord = vec2f{uv[0], radius + (1 - radius) * uv[1]};
        });

    for (auto i = 0; i < (pow2(tesselation + 2) + 1); i++)
        mnorm[i] = normalize(mnorm[i] + vec3f{0, 0, 1});
    merge_quads(quads, pos, norm, texcoord, mquads, mpos, mnorm, mtexcoord);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// EXAMPLE SCENES
// -----------------------------------------------------------------------------
namespace ygl {

// cleanup
proc_shape::~proc_shape() {
    if (hair_params) delete hair_params;
}

// cleanup
proc_scene::~proc_scene() {
    for (auto v : cameras) delete v;
    for (auto v : textures) delete v;
    for (auto v : materials) delete v;
    for (auto v : shapes) delete v;
    for (auto v : instances) delete v;
    for (auto v : environments) delete v;
    for (auto v : nodes) delete v;
    for (auto v : animations) delete v;
}

// Makes/updates a test texture
void update_proc_elem(
    const scene* scn, texture* txt, const proc_texture* ptxt) {
    if (ptxt->name == "") throw std::runtime_error("cannot use empty name");

    txt->name = ptxt->name;
    txt->path = "";
    txt->ldr = {};
    txt->hdr = {};

    switch (ptxt->type) {
        case proc_texture_type::none: break;
        case proc_texture_type::grid: {
            txt->ldr = make_grid_image(ptxt->resolution, ptxt->resolution);
        } break;
        case proc_texture_type::checker: {
            txt->ldr = make_checker_image(ptxt->resolution, ptxt->resolution);
        } break;
        case proc_texture_type::colored: {
            txt->ldr = make_uvgrid_image(ptxt->resolution, ptxt->resolution);
        } break;
        case proc_texture_type::rcolored: {
            txt->ldr = make_recuvgrid_image(ptxt->resolution, ptxt->resolution);
        } break;
        case proc_texture_type::bump: {
            txt->ldr = make_bumpdimple_image(
                ptxt->resolution, ptxt->resolution, ptxt->tile_size);
        } break;
        case proc_texture_type::uv: {
            txt->ldr = make_uv_image(ptxt->resolution, ptxt->resolution);
        } break;
        case proc_texture_type::gamma: {
            txt->ldr = make_gammaramp_image(ptxt->resolution, ptxt->resolution);
        } break;
        case proc_texture_type::noise: {
            txt->ldr = make_noise_image(
                ptxt->resolution, ptxt->resolution, ptxt->noise_scale);
        } break;
        case proc_texture_type::ridge: {
            txt->ldr = make_ridge_image(
                ptxt->resolution, ptxt->resolution, ptxt->noise_scale);
        } break;
        case proc_texture_type::fbm: {
            txt->ldr = make_fbm_image(
                ptxt->resolution, ptxt->resolution, ptxt->noise_scale);
        } break;
        case proc_texture_type::turbulence: {
            txt->ldr = make_turbulence_image(
                ptxt->resolution, ptxt->resolution, ptxt->noise_scale);
        } break;
        case proc_texture_type::gammaf: {
            txt->hdr =
                make_gammaramp_imagef(ptxt->resolution, ptxt->resolution);
        } break;
        case proc_texture_type::sky: {
            txt->hdr = make_sunsky_image(ptxt->resolution, ptxt->sky_sunangle);
        } break;
        default: throw std::runtime_error("should not have gotten here");
    }

    if (ptxt->bump_to_normal) {
        txt->ldr = bump_to_normal_map(txt->ldr, ptxt->bump_scale);
    }

    if (!txt->ldr.empty()) txt->path = ptxt->name + ".png";
    if (!txt->hdr.empty()) txt->path = ptxt->name + ".hdr";
}

// Makes/updates a test material
void update_proc_elem(
    const scene* scn, material* mat, const proc_material* pmat) {
    if (pmat->name == "") throw std::runtime_error("cannot use empty name");

    auto txt = find_named_elem(scn->textures, pmat->texture);
    auto norm = find_named_elem(scn->textures, pmat->normal);

    mat->name = pmat->name;
    mat->type = material_type::specular_roughness;
    mat->ke = zero3f;
    mat->kd = zero3f;
    mat->rs = 1;
    mat->kr = zero3f;
    mat->kt = zero3f;
    mat->ke_txt = nullptr;
    mat->kd_txt = nullptr;
    mat->ks_txt = nullptr;
    mat->kr_txt = nullptr;
    mat->kt_txt = nullptr;

    switch (pmat->type) {
        case proc_material_type::none: break;
        case proc_material_type::emission: {
            mat->ke = pmat->emission * pmat->color;
            mat->ke_txt = txt;
        } break;
        case proc_material_type::matte: {
            mat->kd = pmat->color;
            mat->kd_txt = txt;
        } break;
        case proc_material_type::plastic: {
            mat->kd = pmat->color;
            mat->ks = {0.04f, 0.04f, 0.04f};
            mat->rs = pmat->roughness;
            mat->kd_txt = txt;
        } break;
        case proc_material_type::metal: {
            mat->ks = pmat->color;
            mat->rs = pmat->roughness;
            mat->ks_txt = txt;
        } break;
        case proc_material_type::transparent: {
            mat->kd = pmat->color;
            mat->op = pmat->opacity;
            mat->kd_txt = txt;
        } break;
        default: throw std::runtime_error("should not have gotten here");
    }

    mat->norm_txt = norm;
}

// Makes/updates a test shape
void update_proc_elem(
    const scene* scn, shape_group* sgr, const proc_shape* pshp) {
    if (pshp->name == "") throw std::runtime_error("cannot use empty name");
    auto nshapes = 1;
    if (pshp->type == proc_shape_type::matball ||
        pshp->type == proc_shape_type::hairball)
        nshapes = 2;
    if (sgr->shapes.size() != nshapes) {
        for (auto v : sgr->shapes) delete v;
        sgr->shapes.clear();
        for (auto sid = 0; sid < nshapes; sid++)
            sgr->shapes.push_back(new shape());
    }

    sgr->name = pshp->name;
    auto sid = 0;
    for (auto shp : sgr->shapes) {
        shp->name =
            pshp->name + ((sid > 0) ? std::to_string(sid++) : std::string(""));
        shp->mat = find_named_elem(scn->materials, pshp->material);
        shp->pos = {};
        shp->norm = {};
        shp->texcoord = {};
        shp->texcoord1 = {};
        shp->color = {};
        shp->radius = {};
        shp->tangsp = {};
        shp->points = {};
        shp->lines = {};
        shp->triangles = {};
        shp->quads = {};
        shp->quads_pos = {};
        shp->quads_norm = {};
        shp->quads_texcoord = {};
    }

    auto shp = sgr->shapes.front();
    switch (pshp->type) {
        case proc_shape_type::floor: {
            make_uvquad(shp->quads, shp->pos, shp->norm, shp->texcoord,
                (pshp->tesselation < 0) ? 5 : pshp->tesselation);
            for (auto& p : shp->pos) p = {-p.x, p.z, p.y};
            for (auto& n : shp->norm) n = {n.x, n.z, n.y};
            for (auto& p : shp->pos) p *= 20;
            for (auto& uv : shp->texcoord) uv *= 20;
        } break;
        case proc_shape_type::quad: {
            make_uvquad(shp->quads, shp->pos, shp->norm, shp->texcoord,
                (pshp->tesselation < 0) ? 0 : pshp->tesselation);
        } break;
        case proc_shape_type::cube: {
            make_uvcube(shp->quads, shp->pos, shp->norm, shp->texcoord,
                (pshp->tesselation < 0) ? 0 : pshp->tesselation);
        } break;
        case proc_shape_type::sphere: {
            make_uvsphere(shp->quads, shp->pos, shp->norm, shp->texcoord,
                (pshp->tesselation < 0) ? 5 : pshp->tesselation);
        } break;
        case proc_shape_type::spherecube: {
            make_uvspherecube(shp->quads, shp->pos, shp->norm, shp->texcoord,
                (pshp->tesselation < 0) ? 4 : pshp->tesselation);
        } break;
        case proc_shape_type::spherizedcube: {
            make_uvspherizedcube(shp->quads, shp->pos, shp->norm, shp->texcoord,
                (pshp->tesselation < 0) ? 4 : pshp->tesselation, 0.75f);
        } break;
        case proc_shape_type::geosphere: {
            make_geodesicsphere(shp->triangles, shp->pos,
                (pshp->tesselation < 0) ? 5 : pshp->tesselation);
            shp->norm = shp->pos;
        } break;
        case proc_shape_type::flipcapsphere: {
            make_uvflipcapsphere(shp->quads, shp->pos, shp->norm, shp->texcoord,
                (pshp->tesselation < 0) ? 5 : pshp->tesselation, 0.75f);
        } break;
        case proc_shape_type::suzanne: {
            make_suzanne(shp->quads, shp->pos,
                (pshp->tesselation < 0) ? 0 : pshp->tesselation);
        } break;
        case proc_shape_type::cubep: {
            make_cube(shp->quads, shp->pos,
                (pshp->tesselation < 0) ? 0 : pshp->tesselation);
        } break;
        case proc_shape_type::fvcube: {
            make_fvcube(shp->quads_pos, shp->pos, shp->quads_norm, shp->norm,
                shp->quads_texcoord, shp->texcoord,
                (pshp->tesselation < 0) ? 0 : pshp->tesselation);
        } break;
        case proc_shape_type::fvsphere: {
            make_uvsphere(shp->quads, shp->pos, shp->norm, shp->texcoord,
                (pshp->tesselation < 0) ? 5 : pshp->tesselation);
        } break;
        case proc_shape_type::matball: {
            make_uvflipcapsphere(shp->quads, shp->pos, shp->norm, shp->texcoord,
                (pshp->tesselation < 0) ? 5 : pshp->tesselation, 0.75f);
            auto shp1 = sgr->shapes.at(1);
            make_uvsphere(shp1->quads, shp1->pos, shp1->norm, shp1->texcoord,
                (pshp->tesselation < 0) ? 5 : pshp->tesselation);
            for (auto& p : shp1->pos) p *= 0.8f;
            shp1->mat = find_named_elem(scn->materials, pshp->interior);
        } break;
        case proc_shape_type::point: {
            shp->points.push_back(0);
            shp->pos.push_back({0, 0, 0});
            shp->norm.push_back({0, 0, 1});
            shp->radius.push_back(0.001f);
        } break;
        case proc_shape_type::pointscube: {
            auto npoints = (pshp->num < 0) ? 64 * 64 * 16 : pshp->num;
            auto radius = (pshp->radius < 0) ? 0.0025f : pshp->radius;
            make_points_uv(shp->points, shp->texcoord, npoints);
            shp->pos.reserve(shp->texcoord.size());
            shp->norm.resize(shp->texcoord.size(), {0, 0, 1});
            shp->radius.resize(shp->texcoord.size(), radius);
            auto rn = init_rng(0);
            for (auto i = 0; i < shp->texcoord.size(); i++) {
                shp->pos.push_back(vec3f{-1 + 2 * next_rand1f(rn),
                    -1 + 2 * next_rand1f(rn), -1 + 2 * next_rand1f(rn)});
            }
        } break;
        case proc_shape_type::hairball: {
            auto shp1 = sgr->shapes.at(1);
            make_uvspherecube(
                shp1->quads, shp1->pos, shp1->norm, shp1->texcoord, 5);
            shp1->mat = find_named_elem(scn->materials, pshp->interior);
            auto nhairs = (pshp->num < 0) ? 65536 : pshp->num;
            // auto radius = (pshp->radius < 0) ? vec2f{0.001f, 0.0001f}
            // :
            //                                   vec2f{pshp->radius,
            //                                   0.0001f};
            make_hair(shp->lines, shp->pos, shp->norm, shp->texcoord,
                shp->radius, nhairs, 2, {}, shp1->quads, shp1->pos, shp1->norm,
                shp1->texcoord,
                (pshp->hair_params) ? *pshp->hair_params : make_hair_params());
        } break;
        case proc_shape_type::beziercircle: {
            make_bezier_circle(shp->beziers, shp->pos);
            shp->subdivision = 2;
        } break;
        default: throw std::runtime_error("should not have gotten here");
    }

    if (pshp->scale != 1) {
        for (auto& p : shp->pos) p *= pshp->scale;
    }

    for (auto i = 0; i < pshp->subdivision; i++) {
        subdivide_shape_once(shp, true);
    }

    if (pshp->faceted) facet_shape(shp);
}

// Makes/updates a test shape.
void update_proc_elem(
    const scene* scn, instance* ist, const proc_instance* pist) {
    if (pist->name == "") throw std::runtime_error("cannot use empty name");
    ist->name = pist->name;
    ist->frame = pist->frame;
    if (pist->rotation != zero3f) {
        ist->frame =
            pist->frame *
            rotation_frame(vec3f{0, 0, 1}, pist->rotation.z * pif / 180) *
            rotation_frame(vec3f{0, 1, 0}, pist->rotation.y * pif / 180) *
            rotation_frame(vec3f{1, 0, 0}, pist->rotation.x * pif / 180);
    }
    ist->shp = find_named_elem(scn->shapes, pist->shape);
}

// Makes/updates a test shape
void update_proc_elem(const scene* scn, camera* cam, const proc_camera* pcam) {
    if (pcam->name == "") throw std::runtime_error("cannot use empty name");
    cam->name = pcam->name;
    cam->frame = lookat_frame(pcam->from, pcam->to, vec3f{0, 1, 0});
    cam->yfov = pcam->yfov;
    cam->aspect = pcam->aspect;
    cam->near = 0.01f;
    cam->far = 10000;
    cam->aperture = 0;
    cam->focus = length(pcam->from - pcam->to);
}

// Makes/updates a test shape
void update_proc_elem(
    const scene* scn, environment* env, const proc_environment* penv) {
    if (penv->name == "") throw std::runtime_error("cannot use empty name");
    env->name = penv->name;
    env->frame = identity_frame3f;
    if (penv->rotation) {
        env->frame = rotation_frame(vec3f{0, 1, 0}, penv->rotation);
    }
    env->ke = penv->emission * penv->color;
    env->ke_txt = find_named_elem(scn->textures, penv->texture);
}

// Makes/updates a test shape
void update_proc_elem(const scene* scn, node* nde, const proc_node* pnde) {
    if (pnde->name == "") throw std::runtime_error("cannot use empty name");
    nde->name = pnde->name;
    nde->frame = pnde->frame;
    nde->translation = pnde->translation;
    nde->rotation = pnde->rotation;
    nde->scaling = pnde->scaling;
    nde->cam = find_named_elem(scn->cameras, pnde->camera);
    nde->ist = find_named_elem(scn->instances, pnde->instance);
    nde->env = find_named_elem(scn->environments, pnde->environment);
}

// Makes/updates a test animation
void update_proc_elem(
    const scene* scn, animation_group* agr, const proc_animation* panm) {
    if (panm->name == "") throw std::runtime_error("cannot use empty name");
    if (agr->animations.size() != 1) {
        for (auto v : agr->animations) delete v;
        agr->animations.clear();
        agr->animations.push_back(new animation());
    }
    agr->name = panm->name;
    auto anm = agr->animations.front();
    anm->name = panm->name;
    anm->type = (!panm->bezier) ? keyframe_type::linear : keyframe_type::bezier;
    anm->times = panm->times;
    for (auto& v : anm->times) v *= panm->speed;
    anm->translation = panm->translation;
    anm->rotation = panm->rotation;
    anm->scaling = panm->scaling;
    for (auto& v : anm->translation) v *= panm->scale;
    for (auto& v : anm->scaling) v *= panm->scale;
    agr->targets.clear();
    for (auto& nde : panm->nodes)
        agr->targets.push_back({anm, find_named_elem(scn->nodes, nde)});
}

// Update test elements
template <typename T, typename T1>
void update_proc_scene_elem(scene* scn, std::vector<T*>& elems,
    const std::vector<T1*>& telems, const std::unordered_set<void*>& refresh) {
    auto emap = std::unordered_map<std::string, T*>();
    for (auto elem : elems) emap[elem->name] = elem;
    for (auto telem : telems) {
        if (!contains(emap, telem->name)) {
            elems.push_back(new T());
            update_proc_elem(scn, elems.back(), telem);
        } else {
            auto elem = emap.at(telem->name);
            if (contains(refresh, elem)) update_proc_elem(scn, elem, telem);
        }
    }
}

// Makes/updates a test scene
void update_proc_elems(scene* scn, const proc_scene* pscn,
    const std::unordered_set<void*>& refresh) {
    update_proc_scene_elem(scn, scn->cameras, pscn->cameras, refresh);
    update_proc_scene_elem(scn, scn->textures, pscn->textures, refresh);
    update_proc_scene_elem(scn, scn->materials, pscn->materials, refresh);
    update_proc_scene_elem(scn, scn->shapes, pscn->shapes, refresh);
    update_proc_scene_elem(scn, scn->instances, pscn->instances, refresh);
    update_proc_scene_elem(scn, scn->environments, pscn->environments, refresh);
    update_proc_scene_elem(scn, scn->nodes, pscn->nodes, refresh);
    update_proc_scene_elem(scn, scn->animations, pscn->animations, refresh);
}

// remove duplicate elems
template <typename T>
void remove_duplicate_elems(std::vector<T*>& telems) {
    auto names = std::unordered_set<std::string>();
    for (auto elem : telems) names.insert(elem->name);
    if (names.size() == telems.size()) return;
    names.clear();
    auto ntelems = std::vector<T*>();
    auto tdel = std::vector<T*>();
    for (auto elem : telems) {
        if (contains(names, elem->name)) {
            tdel.push_back(elem);
        } else {
            ntelems.push_back(elem);
            names.insert(elem->name);
        }
    }
    for (auto e : tdel) delete e;
    telems = ntelems;
}

// remove duplicates
void remove_duplicates(proc_scene* scn) {
    remove_duplicate_elems(scn->cameras);
    remove_duplicate_elems(scn->textures);
    remove_duplicate_elems(scn->materials);
    remove_duplicate_elems(scn->shapes);
    remove_duplicate_elems(scn->instances);
    remove_duplicate_elems(scn->environments);
}

std::unordered_map<std::string, proc_texture*>& proc_texture_presets() {
    static auto presets = std::unordered_map<std::string, proc_texture*>();
    if (!presets.empty()) return presets;

    auto make_test_texture = [](const std::string& name,
                                 proc_texture_type type) {
        auto params = new proc_texture();
        params->name = name;
        params->type = type;
        return params;
    };

    presets["grid"] = make_test_texture("grid", proc_texture_type::grid);
    presets["checker"] =
        make_test_texture("checker", proc_texture_type::checker);
    presets["colored"] =
        make_test_texture("colored", proc_texture_type::colored);
    presets["rcolored"] =
        make_test_texture("rcolored", proc_texture_type::rcolored);
    presets["bump"] = make_test_texture("bump", proc_texture_type::bump);
    presets["bump"]->tile_size = 32;
    presets["tgrid"] = make_test_texture("tgrid", proc_texture_type::bump);
    presets["tgrid"]->tile_size = 32;
    presets["uv"] = make_test_texture("uv", proc_texture_type::uv);
    presets["gamma"] = make_test_texture("gamma", proc_texture_type::gamma);
    presets["gridn"] = make_test_texture("gridn", proc_texture_type::grid);
    presets["gridn"]->bump_to_normal = true;
    presets["gridn"]->bump_to_normal = true;
    presets["gridn"]->bump_scale = 4;
    presets["tgridn"] = make_test_texture("tgridn", proc_texture_type::grid);
    presets["tgridn"]->tile_size = 32;
    presets["tgridn"]->bump_to_normal = true;
    presets["tgridn"]->bump_to_normal = true;
    presets["tgridn"]->bump_scale = 4;
    presets["bumpn"] = make_test_texture("bumpn", proc_texture_type::bump);
    presets["bumpn"]->tile_size = 32;
    presets["bumpn"]->bump_to_normal = true;
    presets["bumpn"]->bump_scale = 4;
    presets["noise"] = make_test_texture("noise", proc_texture_type::noise);
    presets["ridge"] = make_test_texture("ridge", proc_texture_type::ridge);
    presets["fbm"] = make_test_texture("fbm", proc_texture_type::fbm);
    presets["turbulence"] =
        make_test_texture("turbulence", proc_texture_type::turbulence);

    presets["gammaf"] = make_test_texture("gammaf", proc_texture_type::gammaf);
    presets["sky1"] = make_test_texture("sky1", proc_texture_type::sky);
    presets["sky1"]->sky_sunangle = pif / 4;
    presets["sky2"] = make_test_texture("sky2", proc_texture_type::sky);
    presets["sky2"]->sky_sunangle = pif / 2;

    return presets;
}

std::unordered_map<std::string, proc_material*>& proc_material_presets() {
    static auto presets = std::unordered_map<std::string, proc_material*>();
    if (!presets.empty()) return presets;

    auto make_test_material = [](const std::string& name,
                                  proc_material_type type, const vec3f& color,
                                  float roughness = 1) {
        auto params = new proc_material();
        params->name = name;
        params->type = type;
        params->color = color;
        params->roughness = roughness;
        return params;
    };
    auto make_test_materialt =
        [](const std::string& name, proc_material_type type,
            const std::string& txt, float roughness = 1) {
            auto params = new proc_material();
            params->name = name;
            params->type = type;
            params->color = {1, 1, 1};
            params->roughness = roughness;
            params->texture = txt;
            return params;
        };

    auto emission = proc_material_type::emission;
    auto matte = proc_material_type::matte;
    auto plastic = proc_material_type::plastic;
    auto metal = proc_material_type::metal;
    auto transparent = proc_material_type::transparent;

    auto gray = vec3f{0.2f, 0.2f, 0.2f};
    auto lgray = vec3f{0.5f, 0.5f, 0.5f};
    auto red = vec3f{0.5f, 0.2f, 0.2f};
    auto green = vec3f{0.2f, 0.5f, 0.2f};
    auto blue = vec3f{0.2f, 0.2f, 0.5f};
    auto white = vec3f{1, 1, 1};

    auto gold = vec3f{0.66f, 0.45f, 0.34f};

    auto rough = 0.25f;
    auto sharp = 0.05f;

    auto params = std::vector<proc_material>();

    presets["matte_floor"] = make_test_materialt("matte_floor", matte, "grid");

    presets["matte_gray"] = make_test_material("matte_gray", matte, gray);
    presets["matte_red"] = make_test_material("matte_red", matte, red);
    presets["matte_green"] = make_test_material("matte_green", matte, green);
    presets["matte_blue"] = make_test_material("matte_blue", matte, blue);
    presets["matte_grid"] = make_test_materialt("matte_grid", matte, "grid");
    presets["matte_colored"] =
        make_test_materialt("matte_colored", matte, "colored");
    presets["matte_uv"] = make_test_materialt("matte_uv", matte, "uv");

    presets["plastic_red"] =
        make_test_material("plastic_red", plastic, red, rough);
    presets["plastic_green"] =
        make_test_material("plastic_green", plastic, green, rough);
    presets["plastic_blue"] =
        make_test_material("plastic_blue", plastic, blue, sharp);
    presets["plastic_colored"] =
        make_test_materialt("plastic_colored", plastic, "colored", rough);
    presets["plastic_blue_bumped"] =
        make_test_material("plastic_blue_bumped", plastic, blue, sharp);
    presets["plastic_blue_bumped"]->normal = "bumpn";
    presets["plastic_colored_bumped"] = make_test_materialt(
        "plastic_colored_bumped", plastic, "colored", rough);
    presets["plastic_colored_bumped"]->normal = "bumpn";

    presets["silver_sharp"] =
        make_test_material("silver_sharp", metal, lgray, sharp);
    presets["silver_rough"] =
        make_test_material("silver_rough", metal, lgray, rough);
    presets["gold_sharp"] =
        make_test_material("gold_sharp", metal, gold, sharp);
    presets["gold_rough"] =
        make_test_material("gold_rough", metal, gold, rough);

    presets["transparent_red"] =
        make_test_material("transparent_red", transparent, red);
    presets["transparent_red"]->opacity = 0.9f;
    presets["transparent_green"] =
        make_test_material("transparent_green", transparent, green);
    presets["transparent_green"]->opacity = 0.5f;
    presets["transparent_blue"] =
        make_test_material("transparent_blue", transparent, blue);
    presets["transparent_blue"]->opacity = 0.2f;

    presets["pointlight"] = make_test_material("pointlight", emission, white);
    presets["pointlight"]->emission = 80;
    presets["arealight"] = make_test_material("arealight", emission, white);
    presets["arealight"]->emission = 80;

    return presets;
}

std::unordered_map<std::string, proc_shape*>& proc_shape_presets() {
    static auto presets = std::unordered_map<std::string, proc_shape*>();
    if (!presets.empty()) return presets;

    auto make_test_shape = [](const std::string& name, proc_shape_type type,
                               int tesselation = -1, int subdivision = 0,
                               bool faceted = false) {
        auto params = new proc_shape();
        params->name = name;
        params->type = type;
        params->tesselation = tesselation;
        params->subdivision = subdivision;
        params->faceted = faceted;
        return params;
    };

    presets["floor"] = make_test_shape("floor", proc_shape_type::floor);
    presets["quad"] = make_test_shape("quad", proc_shape_type::quad);
    presets["cube"] = make_test_shape("cube", proc_shape_type::cube);
    presets["sphere"] = make_test_shape("sphere", proc_shape_type::sphere);
    presets["spherecube"] =
        make_test_shape("spherecube", proc_shape_type::spherecube);
    presets["spherizedcube"] =
        make_test_shape("spherizedcube", proc_shape_type::spherizedcube);
    presets["flipcapsphere"] =
        make_test_shape("flipcapsphere", proc_shape_type::flipcapsphere);
    presets["geosphere"] =
        make_test_shape("geosphere", proc_shape_type::geosphere, 5);
    presets["geospheref"] =
        make_test_shape("geospheref", proc_shape_type::geosphere, 5, 0, true);
    presets["geospherel"] =
        make_test_shape("geospherel", proc_shape_type::geosphere, 4, 0, true);
    presets["cubep"] = make_test_shape("cubep", proc_shape_type::cubep);
    presets["cubes"] = make_test_shape("cubes", proc_shape_type::cubep, 0, 4);
    presets["suzanne"] = make_test_shape("suzanne", proc_shape_type::suzanne);
    presets["suzannes"] =
        make_test_shape("suzannes", proc_shape_type::suzanne, 0, 2);
    presets["cubefv"] = make_test_shape("cubefv", proc_shape_type::fvcube);
    presets["cubefvs"] =
        make_test_shape("cubefvs", proc_shape_type::fvcube, 0, 4);
    presets["spherefv"] =
        make_test_shape("spherefv", proc_shape_type::fvsphere);
    presets["matball"] = make_test_shape("matball", proc_shape_type::matball);
    presets["matballi"] = make_test_shape("matballi", proc_shape_type::sphere);
    presets["matballi"]->scale = 0.8f;
    presets["pointscube"] =
        make_test_shape("pointscube", proc_shape_type::pointscube);
    presets["hairball1"] =
        make_test_shape("hairball1", proc_shape_type::hairball);
    presets["hairball1"]->hair_params = new make_hair_params();
    presets["hairball1"]->hair_params->radius = {0.001f, 0.0001f};
    presets["hairball1"]->hair_params->length = {0.1f, 0.1f};
    presets["hairball1"]->hair_params->noise = {0.5f, 8};
    presets["hairball2"] =
        make_test_shape("hairball2", proc_shape_type::hairball);
    presets["hairball2"]->hair_params = new make_hair_params();
    presets["hairball2"]->hair_params->radius = {0.001f, 0.0001f};
    presets["hairball2"]->hair_params->length = {0.1f, 0.1f};
    presets["hairball2"]->hair_params->clump = {0.5f, 128};
    presets["hairball3"] =
        make_test_shape("hairball3", proc_shape_type::hairball);
    presets["hairball3"]->hair_params = new make_hair_params();
    presets["hairball3"]->hair_params->radius = {0.001f, 0.0001f};
    presets["hairball3"]->hair_params->length = {0.1f, 0.1f};
    presets["hairballi"] =
        make_test_shape("hairballi", proc_shape_type::sphere);
    presets["hairballi"]->scale = 0.8f;
    presets["beziercircle"] =
        make_test_shape("beziercircle", proc_shape_type::beziercircle);
    presets["point"] = make_test_shape("point", proc_shape_type::point);

    return presets;
}

std::unordered_map<std::string, proc_environment*>& proc_environment_presets() {
    static auto presets = std::unordered_map<std::string, proc_environment*>();
    if (!presets.empty()) return presets;

    auto make_test_environment = [](const std::string& name,
                                     const std::string& texture) {
        auto params = new proc_environment();
        params->name = name;
        params->color = {1, 1, 1};
        params->texture = texture;
        return params;
    };

    presets["const"] = make_test_environment("const", "");
    presets["sky1"] = make_test_environment("sky1", "sky1");
    presets["sky2"] = make_test_environment("sky2", "sky2");

    return presets;
}

std::unordered_map<std::string, proc_animation*>& proc_animation_presets() {
    static auto presets = std::unordered_map<std::string, proc_animation*>();
    if (!presets.empty()) return presets;

    auto make_test_animation = [](const std::string& name, bool bezier,
                                   const std::vector<float>& times,
                                   const std::vector<vec3f>& translation,
                                   const std::vector<quat4f>& rotation,
                                   const std::vector<vec3f>& scaling) {
        auto params = new proc_animation();
        params->name = name;
        params->speed = 1;
        params->scale = 1;
        params->bezier = bezier;
        params->times = times;
        params->translation = translation;
        params->rotation = rotation;
        params->scaling = scaling;
        return params;
    };

    presets["bounce"] = make_test_animation(
        "bounce", false, {0, 1, 2}, {{0, 0, 0}, {0, 1, 0}, {0, 0, 0}}, {}, {});
    presets["scale"] = make_test_animation("scale", false, {0, 1, 2}, {}, {},
        {{1, 1, 1}, {0.1f, 0.1f, 0.1f}, {1, 1, 1}});
    presets["rotation"] = make_test_animation("rotation", false, {0, 1, 2}, {},
        {rotation_quat<float>({0, 1, 0}, 0),
            rotation_quat<float>({0, 1, 0}, pif),
            rotation_quat<float>({0, 1, 0}, 0)},
        {});

    return presets;
}

std::unordered_map<std::string, proc_camera*>& proc_camera_presets() {
    static auto presets = std::unordered_map<std::string, proc_camera*>();
    if (!presets.empty()) return presets;

    auto make_test_camera = [](const std::string& name, const vec3f& from,
                                const vec3f& to, float yfov, float aspect) {
        auto params = new proc_camera();
        params->name = name;
        params->from = from;
        params->to = to;
        params->yfov = yfov;
        params->aspect = aspect;
        return params;
    };

    presets["cam1"] =
        make_test_camera("cam1", {0, 4, 10}, {0, 1, 0}, 15 * pif / 180, 1);
    presets["cam2"] = make_test_camera(
        "cam2", {0, 4, 10}, {0, 1, 0}, 15 * pif / 180, 16.0f / 9.0f);
    presets["cam3"] = make_test_camera(
        "cam3", {0, 6, 24}, {0, 1, 0}, 7.5f * pif / 180, 2.35f / 1.0f);

    return presets;
}

std::unordered_map<std::string, proc_scene*>& proc_scene_presets() {
    static auto presets = std::unordered_map<std::string, proc_scene*>();
    if (!presets.empty()) return presets;

    auto make_test_scene = [](const std::string& name) {
        auto params = new proc_scene();
        params->name = name;
        return params;
    };
    auto make_test_instance = [](const std::string& name,
                                  const std::string& shape,
                                  const vec3f& pos = {0, 0, 0}) {
        auto params = new proc_instance();
        params->name = name;
        params->shape = shape;
        params->frame.o = pos;
        return params;
    };
    auto make_texture_preset = [](const std::string& name) {
        return new proc_texture(*proc_texture_presets().at(name));
    };
    auto make_shape_preset = [](const std::string& name) {
        auto npshp = new proc_shape(*proc_shape_presets().at(name));
        if (npshp->hair_params)
            npshp->hair_params = new make_hair_params(*npshp->hair_params);
        return npshp;
    };
    auto make_environment_preset = [](const std::string& name) {
        return new proc_environment(*proc_environment_presets().at(name));
    };
    auto make_camera_preset = [](const std::string& name) {
        return new proc_camera(*proc_camera_presets().at(name));
    };
    auto make_material_preset = [](const std::string& name) {
        return new proc_material(*proc_material_presets().at(name));
    };
    auto make_animation_preset = [](const std::string& name) {
        return new proc_animation(*proc_animation_presets().at(name));
    };

    // textures
    presets["textures"] = make_test_scene("textures");
    for (auto& txt_kv : proc_texture_presets())
        presets["textures"]->textures.push_back(
            make_texture_preset(txt_kv.first));

    // shapes
    presets["shapes"] = make_test_scene("shapes");
    for (auto& shp_kv : proc_shape_presets())
        presets["shapes"]->shapes.push_back(make_shape_preset(shp_kv.first));

    // envmap
    presets["environments"] = make_test_scene("envmaps");
    for (auto& env_kv : proc_environment_presets())
        presets["environments"]->environments.push_back(
            make_environment_preset(env_kv.first));

    // simple scenes shared functions
    auto make_simple_scene = [&](const std::string& name,
                                 const std::vector<std::string>& shapes,
                                 const std::vector<std::string>& mats,
                                 const std::string& lights, bool nodes = false,
                                 const std::vector<std::string>& animations =
                                     {}) {
        auto pos =
            std::vector<vec3f>{{-2.50f, 1, 0}, {0, 1, 0}, {+2.50f, 1, 0}};
        auto params = make_test_scene(name);
        params->cameras.push_back(make_camera_preset("cam3"));
        params->materials.push_back(make_material_preset("matte_floor"));
        if (params->materials.back()->texture != "")
            params->textures.push_back(
                make_texture_preset(params->materials.back()->texture));
        params->shapes.push_back(make_shape_preset("floor"));
        params->shapes.back()->material = params->materials.back()->name;
        params->instances.push_back(
            make_test_instance("floor", "floor", {0, 0, 0}));
        for (auto i = 0; i < shapes.size(); i++) {
            auto name = "obj" + std::to_string(i + 1);
            params->materials.push_back(make_material_preset(mats[i]));
            params->shapes.push_back(make_shape_preset(shapes[i]));
            params->shapes.back()->name = name;
            params->shapes.back()->material = mats[i];
            if (params->shapes.back()->type == proc_shape_type::hairball ||
                params->shapes.back()->type == proc_shape_type::matball) {
                params->materials.push_back(make_material_preset("matte_gray"));
                params->shapes.back()->interior = "matte_gray";
            }
            params->instances.push_back(make_test_instance(name, name, pos[i]));
            if (!animations.empty()) {
                params->animations.push_back(
                    make_animation_preset(animations[i]));
                params->animations.back()->nodes.push_back(name);
            }
        }
        if (lights == "pointlights" || lights == "arealights" ||
            lights == "arealights1") {
            auto emission = 120;
            auto shp = "point";
            auto mat = "pointlight";
            auto pos = std::vector<vec3f>{{-2, 10, 8}, {+2, 10, 8}};
            auto scale = 1.0f;
            if (lights == "arealights") {
                emission = 8;
                shp = "quad";
                mat = "arealight";
                pos = {{0, 16, 0}, {0, 16, 16}};
                scale = 8;
            }
            if (lights == "arealights1") {
                emission = 80;
                shp = "quad";
                mat = "arealight";
                pos = {{-4, 5, 8}, {+4, 5, 8}};
                scale = 2;
            }
            for (auto i = 0; i < 2; i++) {
                auto name = "light" + std::to_string(i + 1);
                params->materials.push_back(make_material_preset(mat));
                params->materials.back()->name = name;
                params->materials.back()->emission = emission;
                params->shapes.push_back(make_shape_preset(shp));
                params->shapes.back()->name = name;
                params->shapes.back()->material = name;
                params->shapes.back()->scale = scale;
                params->instances.push_back(
                    make_test_instance(name, name, pos[i]));
                if (lights == "arealights" || lights == "arealights1")
                    params->instances.back()->frame =
                        lookat_frame(pos[i], {0, 1, 0}, {0, 0, 1}, true);
            }
        }
        if (lights == "envlights") {
            params->environments.push_back(make_environment_preset("sky1"));
            params->environments.back()->frame =
                lookat_frame(pos[1], pos[1] + vec3f{0, 0, 1}, {0, 1, 0}, true);
            params->textures.push_back(make_texture_preset("sky1"));
        }
        if (!animations.empty() || nodes) {
            for (auto cam : params->cameras) {
                auto nde = new proc_node();
                nde->name = cam->name;
                nde->frame = lookat_frame(cam->from, cam->to, vec3f{0, 1, 0});
                nde->camera = cam->name;
                params->nodes.push_back(nde);
            }
            for (auto ist : params->instances) {
                auto nde = new proc_node();
                nde->name = ist->name;
                nde->frame = ist->frame;
                nde->instance = ist->name;
                params->nodes.push_back(nde);
            }
            for (auto env : params->environments) {
                auto nde = new proc_node();
                nde->name = env->name;
                nde->frame = env->frame;
                nde->environment = env->name;
                params->nodes.push_back(nde);
            }
        }
        return params;
    };

    // plane only
    presets["plane_pl"] = make_simple_scene("plane_pl", {}, {}, "pointlights");
    presets["plane_al"] = make_simple_scene("plane_al", {}, {}, "arealights");

    // basic shapes
    presets["basic_pl"] = make_simple_scene("basic_pl",
        {"flipcapsphere", "spherecube", "spherizedcube"},
        {"plastic_red", "plastic_green", "plastic_blue"}, "pointlights");
    presets["basic_al"] = make_simple_scene("basic_al",
        {"flipcapsphere", "spherecube", "spherizedcube"},
        {"plastic_red", "plastic_green", "plastic_blue"}, "arealights");

    // simple shapes
    presets["simple_pl"] = make_simple_scene("simple_pl",
        {"flipcapsphere", "spherecube", "spherizedcube"},
        {"plastic_colored", "plastic_colored", "plastic_colored"},
        "pointlights");
    presets["simple_al"] = make_simple_scene("simple_al",
        {"flipcapsphere", "spherecube", "spherizedcube"},
        {"plastic_colored", "plastic_colored", "plastic_colored"},
        "arealights");
    presets["simple_el"] = make_simple_scene("simple_el",
        {"flipcapsphere", "spherecube", "spherizedcube"},
        {"plastic_colored", "plastic_colored", "plastic_colored"}, "envlights");

    // transparent shapes
    presets["transparent_al"] =
        make_simple_scene("transparent_al", {"quad", "quad", "quad"},
            {"transparent_red", "transparent_green", "transparent_blue"},
            "arealights");

    // points shapes
    presets["points_al"] = make_simple_scene("transparent_al",
        {"pointscube", "pointscube", "pointscube"},
        {"matte_gray", "matte_gray", "matte_gray"}, "arealights");

    // lines shapes
    presets["lines_al"] =
        make_simple_scene("lines_al", {"hairball1", "hairball2", "hairball3"},
            {"matte_gray", "matte_gray", "matte_gray"}, "arealights");

    // subdiv shapes
    presets["subdiv_al"] =
        make_simple_scene("subdiv_al", {"cubes", "suzannes", "suzannes"},
            {"plastic_red", "plastic_green", "plastic_blue"}, "arealights");

    // plastics shapes
    presets["plastics_al"] =
        make_simple_scene("plastics_al", {"matball", "matball", "matball"},
            {"matte_green", "plastic_green", "plastic_colored"}, "arealights",
            true);
    presets["plastics_el"] =
        make_simple_scene("plastics_el", {"matball", "matball", "matball"},
            {"matte_green", "plastic_green", "plastic_colored"}, "envlights");

    // metals shapes
    presets["metals_al"] =
        make_simple_scene("metals_al", {"matball", "matball", "matball"},
            {"gold_rough", "gold_sharp", "silver_sharp"}, "arealights");
    presets["metals_el"] =
        make_simple_scene("metals_el", {"matball", "matball", "matball"},
            {"gold_rough", "gold_sharp", "silver_sharp"}, "envlights");

    // tesselation shapes
    presets["tesselation_pl"] = make_simple_scene("tesselation_pl",
        {"geospherel", "geospheref", "geosphere"},
        {"matte_gray", "matte_gray", "matte_gray"}, "pointlights");

    // textureuv shapes
    presets["textureuv_pl"] = make_simple_scene("textureuv_pl",
        {"flipcapsphere", "flipcapsphere", "flipcapsphere"},
        {"matte_green", "matte_colored", "matte_uv"}, "pointlights");

    // normalmap shapes
    presets["normalmap_pl"] = make_simple_scene("normalmap_pl",
        {"flipcapsphere", "flipcapsphere", "flipcapsphere"},
        {"plastic_blue", "plastic_blue_bumped", "plastic_colored_bumped"},
        "pointlights");

    // animated shapes
    presets["animated_pl"] = make_simple_scene("animated_pl",
        {"flipcapsphere", "spherecube", "spherizedcube"},
        {"plastic_colored", "plastic_colored", "plastic_colored"},
        "pointlights", true, {"bounce", "scale", "rotation"});

    // instances shared functions
    auto make_random_scene = [&](const std::string& name, const vec2i& num,
                                 const bbox2f& bbox, uint32_t seed = 13) {
        auto rscale = 0.9f * 0.25f *
                      min((bbox.max.x - bbox.min.x) / num.x,
                          (bbox.max.x - bbox.min.x) / num.y);

        auto params = make_test_scene(name);
        params->cameras.push_back(proc_camera_presets().at("cam3"));
        params->materials.push_back(proc_material_presets().at("matte_floor"));
        params->shapes.push_back(make_shape_preset("floor"));
        params->instances.push_back(
            make_test_instance("floor", "floor", {0, 0, 0}));
        auto shapes = std::vector<std::string>();
        for (auto mat : {"plastic_red", "plastic_green", "plastic_blue"})
            params->materials.push_back(make_material_preset(mat));
        for (auto shp : {"sphere", "flipcapsphere", "cube"}) {
            for (auto mat : {"plastic_red", "plastic_green", "plastic_blue"}) {
                params->shapes.push_back(make_shape_preset(shp));
                params->shapes.back()->name += "_"s + mat;
                params->shapes.back()->material = mat;
                params->shapes.back()->scale *= rscale;
                shapes.push_back(params->shapes.back()->name);
            }
        }

        auto rng = init_rng(seed, 7);
        auto count = 0;
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
                params->instances.push_back(
                    make_test_instance("instance" + std::to_string(count++),
                        shapes[next_rand1i(rng, (int)shapes.size())], pos));
            }
        }

        auto pos = std::vector<vec3f>{{-2, 10, 8}, {+2, 10, 8}};
        for (auto i = 0; i < 2; i++) {
            auto name = "light" + std::to_string(i + 1);
            params->materials.push_back(make_material_preset("pointlight"));
            params->materials.back()->name = name;
            params->materials.back()->emission = 80;
            params->shapes.push_back(make_shape_preset("point"));
            params->shapes.back()->name = name;
            params->shapes.back()->material = name;
            params->instances.push_back(make_test_instance(name, name, pos[i]));
        }

        return params;
    };

    // instances
    presets["instances_pl"] =
        make_random_scene("instances_pl", {10, 10}, {{-3, -3}, {3, 3}});
    presets["instancel_pl"] =
        make_random_scene("instancel_pl", {100, 100}, {{-3, -3}, {3, 3}});

#if 0
        else if (otype == "normdisp") {
            auto mat = std::vector<material*>{
                add_test_material(scn, test_material_type::plastic_bumped),
                add_test_material(scn, test_material_type::plastic_bumped)};
            add_test_instance(scn, "obj01",
                add_uvspherecube(scn, "base_obj01", mat[0], 4), {-1.25f, 1, 0});
            add_test_instance(scn, "obj03",
                add_uvspherecube(scn, "subdiv_02_obj02", mat[1], 4), {1.25f, 1, 0});
            }
#endif

    // add missing textures
    for (auto kv : presets) {
        auto preset = kv.second;
        auto used = std::unordered_set<std::string>();
        for (auto mat : preset->materials) used.insert(mat->texture);
        for (auto mat : preset->materials) used.insert(mat->normal);
        for (auto env : preset->environments) used.insert(env->texture);
        used.erase("");
        for (auto txt : preset->textures) used.erase(txt->name);
        for (auto txt : used)
            preset->textures.push_back(make_texture_preset(txt));
    }

    // remove duplicates
    for (auto& kv : presets) remove_duplicates(kv.second);

    return presets;
}

// Parses a test_camera object
void serialize(proc_camera& val, json& js, bool reading) {
    static auto def = proc_camera();
    serialize_obj(js, reading);
    serialize_attr(val.name, js, "name", reading);
    serialize_attr(val.from, js, "from", reading);
    serialize_attr(val.to, js, "to", reading);
    serialize_attr(val.yfov, js, "yfov", reading);
    serialize_attr(val.aspect, js, "aspect", reading);
}

// Parses a test_camera object
void serialize(proc_texture_type& val, json& js, bool reading) {
    serialize(val, js, reading, enum_names(val));
}

// Parses a test_camera object
void serialize(proc_texture& val, json& js, bool reading) {
    static auto def = proc_texture();
    serialize_obj(js, reading);
    serialize_attr(val.name, js, "name", reading);
    serialize_attr(val.type, js, "type", reading);
    serialize_attr(val.resolution, js, "resolution", reading);
    serialize_attr(
        val.tile_size, js, "tile_size", reading, false, def.tile_size);
    serialize_attr(
        val.noise_scale, js, "noise_scale", reading, false, def.noise_scale);
    serialize_attr(
        val.sky_sunangle, js, "sky_sunangle", reading, false, def.sky_sunangle);
    serialize_attr(val.bump_to_normal, js, "bump_to_normal", reading, false,
        def.bump_to_normal);
    serialize_attr(
        val.bump_scale, js, "bump_scale", reading, false, def.bump_scale);
}

// Parses a test_camera object
void serialize(proc_material_type& val, json& js, bool reading) {
    serialize(val, js, reading, enum_names(val));
}

// Parses a test_camera object
void serialize(proc_material& val, json& js, bool reading) {
    static auto def = proc_material();
    serialize_obj(js, reading);
    serialize_attr(val.name, js, "name", reading);
    serialize_attr(val.type, js, "type", reading);
    serialize_attr(val.emission, js, "emission", reading, false, def.emission);
    serialize_attr(val.color, js, "color", reading, false, def.color);
    serialize_attr(val.opacity, js, "opacity", reading, false, def.opacity);
    serialize_attr(
        val.roughness, js, "roughness", reading, false, def.roughness);
    serialize_attr(val.texture, js, "texture", reading, false, def.texture);
    serialize_attr(val.normal, js, "normal", reading, false, def.normal);
}

// Parses a test_camera object
void serialize(proc_shape_type& val, json& js, bool reading) {
    serialize(val, js, reading, enum_names(val));
}

// Parses a test_camera object
void serialize(proc_shape& val, json& js, bool reading) {
    static auto def = proc_shape();
    serialize_obj(js, reading);
    serialize_attr(val.name, js, "name", reading);
    serialize_attr(val.type, js, "type", reading);
    serialize_attr(val.material, js, "material", reading, false, def.material);
    serialize_attr(
        val.tesselation, js, "tesselation", reading, false, def.tesselation);
    serialize_attr(
        val.subdivision, js, "subdivision", reading, false, def.subdivision);
    serialize_attr(val.scale, js, "scale", reading, false, def.scale);
    serialize_attr(val.radius, js, "radius", reading, false, def.radius);
    serialize_attr(val.faceted, js, "faceted", reading, false, def.faceted);
    serialize_attr(val.num, js, "num", reading, false, def.num);
    // TODO: hair parameters
}

// Parses a test_camera object
void serialize(proc_instance& val, json& js, bool reading) {
    static auto def = proc_instance();
    serialize_obj(js, reading);
    serialize_attr(val.name, js, "name", reading);
    serialize_attr(val.shape, js, "shape", reading);
    serialize_attr(val.frame, js, "frame", reading, false, def.frame);
    serialize_attr(val.rotation, js, "rotation", reading, false, def.rotation);
}

// Parses a test_camera object
void serialize(proc_environment& val, json& js, bool reading) {
    static auto def = proc_environment();
    serialize_obj(js, reading);
    serialize_attr(val.name, js, "name", reading);
    serialize_attr(val.emission, js, "emission", reading, false, def.emission);
    serialize_attr(val.color, js, "color", reading, false, def.color);
    serialize_attr(val.texture, js, "texture", reading, false, def.texture);
    serialize_attr(val.frame, js, "frame", reading, false, def.frame);
    serialize_attr(val.rotation, js, "rotation", reading, false, def.rotation);
}

// Parses a test_camera object
void serialize(proc_scene& val, json& js, bool reading) {
    static auto def = proc_scene();
    serialize_obj(js, reading);
    serialize_attr(val.name, js, "name", reading, false, def.name);
    serialize_attr(val.cameras, js, "cameras", reading, false, def.cameras);
    serialize_attr(val.textures, js, "textures", reading, false, def.textures);
    serialize_attr(
        val.materials, js, "materials", reading, false, def.materials);
    serialize_attr(val.shapes, js, "shapes", reading, false, def.shapes);
    serialize_attr(
        val.environments, js, "environments", reading, false, def.environments);
}

// Load test scene
proc_scene* load_proc_scene(const std::string& filename) {
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

    // clear data
    auto scn = new proc_scene();
    try {
        serialize(scn, js, true);
    } catch (const std::exception& e) {
        throw std::runtime_error(
            "error parsing test scene " + std::string(e.what()));
    }

    // done
    return scn;
}

// Save test scene
void save_proc_scene(const std::string& filename, const proc_scene* scn) {
    auto js = json();
    serialize((proc_scene&)*scn, js, false);
    save_text(filename, js.dump(2));
}

}  // namespace ygl

#if YGL_OPENGL

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR OPENGL
// -----------------------------------------------------------------------------
namespace ygl {
// Checks for GL error and then prints
bool gl_check_error(bool print) {
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
void gl_clear_buffers(const vec4f& background) {
    assert(gl_check_error());
    glClearColor(background[0], background[1], background[2], background[3]);
    glClearDepth(1);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    assert(gl_check_error());
}

// Enable/disable depth test
void gl_enable_depth_test(bool enabled) {
    assert(gl_check_error());
    if (enabled) {
        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LEQUAL);
    } else {
        glDisable(GL_DEPTH_TEST);
    }
    assert(gl_check_error());
}

// Enable/disable culling
void gl_enable_culling(bool enabled, bool front, bool back) {
    assert(gl_check_error());
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
    assert(gl_check_error());
}

// Enable/disable wireframe
void gl_enable_wireframe(bool enabled) {
    assert(gl_check_error());
    if (enabled)
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    else
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
    assert(gl_check_error());
}

// Enable/disable blending
void gl_enable_blending(bool enabled) {
    assert(gl_check_error());
    if (enabled) {
        glEnable(GL_BLEND);
    } else {
        glDisable(GL_BLEND);
    }
    assert(gl_check_error());
}

// Set blending to over operator
void gl_set_blend_over() {
    assert(gl_check_error());
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    assert(gl_check_error());
}

// Line width
void gl_line_width(float w) {
    assert(gl_check_error());
    glLineWidth(min(max(w, 0.0f), 1.0f));
    assert(gl_check_error());
}

// Set viewport
void gl_set_viewport(const vec4i& v) {
    assert(gl_check_error());
    glViewport(v.x, v.y, v.z, v.w);
    assert(gl_check_error());
}

// Set viewport
void gl_set_viewport(const vec2i& v) {
    assert(gl_check_error());
    glViewport(0, 0, v.x, v.y);
    assert(gl_check_error());
}

// This is a public API. See above for documentation.
void gl_read_imagef(float* pixels, int w, int h, int nc) {
    assert(gl_check_error());
    int formats[4] = {GL_RED, GL_RG, GL_RGB, GL_RGBA};
    glReadPixels(0, 0, w, h, formats[nc - 1], GL_FLOAT, pixels);
    assert(gl_check_error());
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF OPENGL TEXTURE FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Implementation of update_texture.
void _update_texture(gl_texture& txt, int w, int h, int nc, const void* pixels,
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
    assert(gl_check_error());
    if (w * h) {
        int formats_ub[4] = {GL_RED, GL_RG, GL_RGB, GL_RGBA};
        int formats_sub[4] = {GL_RED, GL_RG, GL_SRGB, GL_SRGB_ALPHA};
        int formats_f[4] = {GL_R32F, GL_RG32F, GL_RGB32F, GL_RGBA32F};
        int* formats =
            (as_float) ? formats_f : ((as_srgb) ? formats_sub : formats_ub);
        assert(gl_check_error());
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
    assert(gl_check_error());
}

// Binds a texture to a texture unit
void bind_texture(const gl_texture& txt, uint unit) {
    glActiveTexture(GL_TEXTURE0 + unit);
    glBindTexture(GL_TEXTURE_2D, txt.tid);
}

// Unbinds
void unbind_texture(const gl_texture& txt, uint unit) {
    glActiveTexture(GL_TEXTURE0 + unit);
    glBindTexture(GL_TEXTURE_2D, 0);
}

// Destroys the texture tid.
void clear_texture(gl_texture& txt) {
    assert(gl_check_error());
    glDeleteTextures(1, &txt.tid);
    txt.tid = 0;
    assert(gl_check_error());
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF OPENGL VERTEX ARRAY BUFFER
// -----------------------------------------------------------------------------
namespace ygl {

// Updates the buffer with new data.
void _update_vertex_buffer(gl_vertex_buffer& buf, int n, int nc,
    const void* values, bool as_float, bool dynamic) {
    auto resize =
        !buf.bid || n * nc != buf.num * buf.ncomp || as_float != buf.as_float;
    buf.num = n;
    buf.ncomp = nc;
    buf.as_float = as_float;
    assert(gl_check_error());
    if (n) {
        if (!buf.bid) glGenBuffers(1, &buf.bid);
        glBindBuffer(GL_ARRAY_BUFFER, buf.bid);
        if (resize) {
            glBufferData(GL_ARRAY_BUFFER,
                buf.num * buf.ncomp *
                    ((as_float) ? sizeof(float) : sizeof(int)),
                values, (dynamic) ? GL_DYNAMIC_DRAW : GL_STATIC_DRAW);
        } else {
            glBufferSubData(GL_ARRAY_BUFFER, 0,
                buf.num * buf.ncomp *
                    ((as_float) ? sizeof(float) : sizeof(int)),
                values);
        }
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    } else {
        if (buf.bid) {
            glBindBuffer(GL_ARRAY_BUFFER, buf.bid);
            glDeleteBuffers(1, &buf.bid);
            buf.bid = 0;
            glBindBuffer(GL_ARRAY_BUFFER, 0);
        }
    }
    assert(gl_check_error());
}

// Bind the buffer at a particular attribute location
void bind_vertex_buffer(const gl_vertex_buffer& buf, uint vattr) {
    glEnableVertexAttribArray(vattr);
    glBindBuffer(GL_ARRAY_BUFFER, buf.bid);
    glVertexAttribPointer(vattr, buf.ncomp, GL_FLOAT, false, 0, 0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

// Unbind the buffer
void unbind_vertex_buffer(const gl_vertex_buffer& buf, uint vattr) {
    glDisableVertexAttribArray(vattr);
}

// Unbind the buffer
void unbind_vertex_buffer(uint vattr) { glDisableVertexAttribArray(vattr); }

// Destroys the buffer
void clear_vertex_buffer(gl_vertex_buffer& buf) {
    assert(gl_check_error());
    glDeleteBuffers(1, &buf.bid);
    buf.bid = 0;
    assert(gl_check_error());
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF OPENGL VERTEX ELEMENTS BUFFER
// -----------------------------------------------------------------------------
namespace ygl {

// Updates the buffer bid with new data.
void _update_element_buffer(
    gl_element_buffer& buf, int n, int nc, const int* values, bool dynamic) {
    auto resize = !buf.bid || n * nc != buf.num * buf.ncomp;
    buf.num = n;
    buf.ncomp = nc;
    assert(gl_check_error());
    if (n) {
        if (!buf.bid) glGenBuffers(1, &buf.bid);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf.bid);
        if (resize) {
            glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                buf.num * buf.ncomp * sizeof(int), values,
                (dynamic) ? GL_DYNAMIC_DRAW : GL_STATIC_DRAW);
        } else {
            glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0,
                buf.num * buf.ncomp * sizeof(int), values);
        }
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    } else {
        if (buf.bid) {
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf.bid);
            glDeleteBuffers(1, &buf.bid);
            buf.bid = 0;
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        }
    }
    assert(gl_check_error());
}

// Draws elements.
void draw_elems(const gl_element_buffer& buf) {
    if (!buf.bid) return;
    assert(gl_check_error());
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
    assert(gl_check_error());
}

// Destroys the buffer
void clear_element_buffer(gl_element_buffer& buf) {
    assert(gl_check_error());
    glDeleteBuffers(1, &buf.bid);
    buf.bid = 0;
    assert(gl_check_error());
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF OPENGL PROGRAM FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Creates and OpenGL program from vertex and fragment code. Returns the
// program id. Optionally return vertex and fragment shader ids. A VAO
// is created.
gl_program make_program(
    const std::string& vertex, const std::string& fragment) {
    auto prog = gl_program();

    assert(gl_check_error());
    glGenVertexArrays(1, &prog.vao);
    glBindVertexArray(prog.vao);
    assert(gl_check_error());

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
    assert(gl_check_error());

    glBindVertexArray(0);
    assert(gl_check_error());

    return prog;
}

// Destroys the program pid and optionally the sahders vid and fid.
void clear_program(gl_program& prog) {
    assert(gl_check_error());
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
    assert(gl_check_error());
}

// Get uniform location (simple GL wrapper that avoids GL includes)
int get_program_uniform_location(
    const gl_program& prog, const std::string& name) {
    assert(gl_check_error());
    return glGetUniformLocation(prog.pid, name.c_str());
    assert(gl_check_error());
}

// Get uniform location (simple GL wrapper that avoids GL includes)
int get_program_attrib_location(
    const gl_program& prog, const std::string& name) {
    assert(gl_check_error());
    return glGetAttribLocation(prog.pid, name.c_str());
    assert(gl_check_error());
}

// Get the names of all uniforms
std::vector<std::pair<std::string, int>> get_program_uniforms_names(
    const gl_program& prog) {
    auto num = 0;
    assert(gl_check_error());
    glGetProgramiv(prog.pid, GL_ACTIVE_UNIFORMS, &num);
    assert(gl_check_error());
    auto names = std::vector<std::pair<std::string, int>>();
    for (auto i = 0; i < num; i++) {
        char name[4096];
        auto size = 0, length = 0;
        GLenum type;
        glGetActiveUniform(prog.pid, i, 4096, &length, &size, &type, name);
        if (length > 3 && name[length - 1] == ']' && name[length - 2] == '0' &&
            name[length - 3] == '[')
            name[length - 3] = 0;
        auto loc = glGetUniformLocation(prog.pid, name);
        if (loc < 0) continue;
        names.push_back({name, loc});
        assert(gl_check_error());
    }
    return names;
}

// Get the names of all attributes
std::vector<std::pair<std::string, int>> get_program_attributes_names(
    const gl_program& prog) {
    auto num = 0;
    assert(gl_check_error());
    glGetProgramiv(prog.pid, GL_ACTIVE_ATTRIBUTES, &num);
    assert(gl_check_error());
    auto names = std::vector<std::pair<std::string, int>>();
    for (auto i = 0; i < num; i++) {
        char name[4096];
        auto size = 0;
        GLenum type;
        glGetActiveAttrib(prog.pid, i, 4096, nullptr, &size, &type, name);
        auto loc = glGetAttribLocation(prog.pid, name);
        if (loc < 0) continue;
        names.push_back({name, loc});
        assert(gl_check_error());
    }
    return names;
}

// Set uniform integer values val for program pid and variable loc.
// The values have nc number of components (1-4) and count elements
// (for arrays).
bool set_program_uniform(
    gl_program& prog, int pos, const int* val, int ncomp, int count) {
    assert(ncomp >= 1 && ncomp <= 4);
    assert(gl_check_error());
    if (pos < 0) return false;
    switch (ncomp) {
        case 1: glUniform1iv(pos, count, val); break;
        case 2: glUniform2iv(pos, count, val); break;
        case 3: glUniform3iv(pos, count, val); break;
        case 4: glUniform4iv(pos, count, val); break;
        default: assert(false);
    }
    assert(gl_check_error());
    return true;
}

// Set uniform float values val for program pid and variable var.
// The values have nc number of components (1-4) and count elements
// (for arrays).
bool set_program_uniform(
    gl_program& prog, int pos, const float* val, int ncomp, int count) {
    assert((ncomp >= 1 && ncomp <= 4) || (ncomp == 16) || (ncomp == 12));
    assert(gl_check_error());
    if (pos < 0) return false;
    switch (ncomp) {
        case 1: glUniform1fv(pos, count, val); break;
        case 2: glUniform2fv(pos, count, val); break;
        case 3: glUniform3fv(pos, count, val); break;
        case 4: glUniform4fv(pos, count, val); break;
        case 12: glUniformMatrix4x3fv(pos, count, false, val); break;
        case 16: glUniformMatrix4fv(pos, count, false, val); break;
        default: assert(false); return 0;
    }
    assert(gl_check_error());
    return true;
}

// Set uniform texture id tid and unit tunit for program pid and
// variable var.
bool set_program_uniform_texture(
    gl_program& prog, int pos, const gl_texture_info& tinfo, uint tunit) {
    static const auto wrap_mode_map =
        std::map<gl_texture_wrap, uint>{{gl_texture_wrap::repeat, GL_REPEAT},
            {gl_texture_wrap::clamp, GL_CLAMP_TO_EDGE},
            {gl_texture_wrap::mirror, GL_MIRRORED_REPEAT}};
    static const auto filter_mode_map = std::map<gl_texture_filter, uint>{
        {gl_texture_filter::nearest, GL_NEAREST},
        {gl_texture_filter::linear, GL_LINEAR},
        {gl_texture_filter::nearest_mipmap_nearest, GL_NEAREST_MIPMAP_NEAREST},
        {gl_texture_filter::linear_mipmap_nearest, GL_LINEAR_MIPMAP_NEAREST},
        {gl_texture_filter::nearest_mipmap_linear, GL_NEAREST_MIPMAP_LINEAR},
        {gl_texture_filter::linear_mipmap_linear, GL_LINEAR_MIPMAP_LINEAR}};

    assert(gl_check_error());
    if (pos < 0) return false;
    if (is_texture_valid(tinfo.txt)) {
        bind_texture(tinfo.txt, tunit);
        if (tinfo.wrap_s != gl_texture_wrap::not_set)
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S,
                wrap_mode_map.at(tinfo.wrap_s));
        if (tinfo.wrap_t != gl_texture_wrap::not_set)
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T,
                wrap_mode_map.at(tinfo.wrap_t));
        if (tinfo.filter_min != gl_texture_filter::not_set)
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER,
                filter_mode_map.at(tinfo.filter_min));
        if (tinfo.filter_mag != gl_texture_filter::not_set)
            glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER,
                filter_mode_map.at(tinfo.filter_mag));
        glUniform1i(pos, tunit);
    } else {
        unbind_texture(tinfo.txt, tunit);
        glUniform1i(pos, tunit);
    }
    assert(gl_check_error());
    return true;
}

// Sets a constant value for a vertex attribute for program pid and
// variable var. The attribute has nc components.
bool set_program_vertattr(
    gl_program& prog, int pos, const float* value, int nc) {
    assert(nc >= 1 && nc <= 4);
    assert(gl_check_error());
    if (pos < 0) return false;
    glDisableVertexAttribArray(pos);
    switch (nc) {
        case 1: glVertexAttrib1fv(pos, value); break;
        case 2: glVertexAttrib2fv(pos, value); break;
        case 3: glVertexAttrib3fv(pos, value); break;
        case 4: glVertexAttrib4fv(pos, value); break;
        default: assert(false); break;
    }
    assert(gl_check_error());
    return true;
}

// Sets a constant value for a vertex attribute for program pid and
// variable var. The attribute has nc components.
bool set_program_vertattr(gl_program& prog, int pos, const int* value, int nc) {
    assert(nc >= 1 && nc <= 4);
    assert(gl_check_error());
    if (pos < 0) return false;
    glDisableVertexAttribArray(pos);
    switch (nc) {
        case 1: glVertexAttribI1iv(pos, value); break;
        case 2: glVertexAttribI2iv(pos, value); break;
        case 3: glVertexAttribI3iv(pos, value); break;
        case 4: glVertexAttribI4iv(pos, value); break;
        default: assert(false); break;
    }
    assert(gl_check_error());
    return true;
}

// Sets a vartex attribute for program pid and variable var to the
// buffer bid. The attribute has nc components and per-vertex values
// values.
bool set_program_vertattr(
    gl_program& prog, const std::string& var, const gl_vertex_buffer& buf) {
    assert(gl_check_error());
    int pos = glGetAttribLocation(prog.pid, var.c_str());
    if (pos < 0) return false;
    if (is_vertex_buffer_valid(buf)) {
        glEnableVertexAttribArray(pos);
        bind_vertex_buffer(buf, pos);
    } else {
        glDisableVertexAttribArray(pos);
        unbind_vertex_buffer(buf, pos);
    }
    assert(gl_check_error());
    return true;
}

// Sets a vartex attribute for program pid and variable var. The
// attribute has nc components and either buffer bid or a single value
// def (if bid is zero). Convenience wrapper to above functions.
bool set_program_vertattr(gl_program& prog, int pos,
    const gl_vertex_buffer& buf, int nc, const float* def) {
    assert(nc >= 1 && nc <= 4);
    assert(gl_check_error());
    if (pos < 0) return false;
    if (is_vertex_buffer_valid(buf)) {
        assert(gl_check_error());
        glEnableVertexAttribArray(pos);
        assert(gl_check_error());
        bind_vertex_buffer(buf, pos);
        assert(gl_check_error());
    } else {
        glDisableVertexAttribArray(pos);
        unbind_vertex_buffer(buf, pos);
        if (def) {
            switch (nc) {
                case 1: glVertexAttrib1fv(pos, def); break;
                case 2: glVertexAttrib2fv(pos, def); break;
                case 3: glVertexAttrib3fv(pos, def); break;
                case 4: glVertexAttrib4fv(pos, def); break;
                default: assert(false); break;
            }
        }
    }
    assert(gl_check_error());
    return true;
}

// Binds a program
void bind_program(const gl_program& prog) {
    assert(gl_check_error());
    if (!prog.pid) return;
    glBindVertexArray(prog.vao);
    glUseProgram(prog.pid);
    assert(gl_check_error());
}

// Unbind a program
void unbind_program(const gl_program& prog) {
    assert(gl_check_error());
    glUseProgram(0);
    glBindVertexArray(0);
    assert(gl_check_error());
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SCENE SHADER FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Clear gl state
void clear_shape(gl_shape& shp) {
    clear_vertex_buffer(shp.pos);
    clear_vertex_buffer(shp.norm);
    clear_vertex_buffer(shp.texcoord);
    clear_vertex_buffer(shp.color);
    clear_vertex_buffer(shp.tangsp);
    clear_element_buffer(shp.points);
    clear_element_buffer(shp.lines);
    clear_element_buffer(shp.triangles);
    clear_element_buffer(shp.quads);
    clear_element_buffer(shp.edges);
}

// Init shading
void update_gl_texture(
    std::unordered_map<texture*, gl_texture>& gtextures, const texture* txt) {
    if (gtextures.find((texture*)txt) == gtextures.end())
        gtextures[(texture*)txt] = gl_texture();
    auto& gtxt = gtextures.at((texture*)txt);
    if (!txt->hdr.empty()) {
        update_texture(gtxt, txt->hdr, true, true, true);
    } else if (!txt->ldr.empty()) {
        update_texture(gtxt, txt->ldr, true, true, true);
    } else
        assert(false);
}

// Update shading
std::unordered_map<texture*, gl_texture> make_gl_textures(const scene* scn) {
    auto gtextures = std::unordered_map<texture*, gl_texture>();
    for (auto txt : scn->textures) update_gl_texture(gtextures, txt);
    return gtextures;
}

// Update shading
void update_gl_shape(
    std::unordered_map<shape*, gl_shape>& gshapes, const shape* shp) {
    if (gshapes.find((shape*)shp) == gshapes.end())
        gshapes[(shape*)shp] = gl_shape();
    auto& gshp = gshapes.at((shape*)shp);
    if (!shp->quads_pos.empty()) {
        auto pos = std::vector<vec3f>();
        auto norm = std::vector<vec3f>();
        auto texcoord = std::vector<vec2f>();
        auto quads = std::vector<vec4i>();
        std::tie(quads, pos, norm, texcoord) =
            convert_face_varying(shp->quads_pos, shp->quads_norm,
                shp->quads_texcoord, shp->pos, shp->norm, shp->texcoord);
        update_vertex_buffer(gshp.pos, pos);
        update_vertex_buffer(gshp.norm, norm);
        update_vertex_buffer(gshp.texcoord, texcoord);
        update_element_buffer(gshp.quads, convert_quads_to_triangles(quads));
        update_element_buffer(gshp.edges, get_edges({}, {}, shp->quads));
        update_vertex_buffer(gshp.color, std::vector<vec4f>{});
        update_vertex_buffer(gshp.tangsp, std::vector<vec4f>{});
    } else {
        update_vertex_buffer(gshp.pos, shp->pos);
        update_vertex_buffer(gshp.norm, shp->norm);
        update_vertex_buffer(gshp.texcoord, shp->texcoord);
        update_vertex_buffer(gshp.color, shp->color);
        update_vertex_buffer(gshp.tangsp, shp->tangsp);
        update_element_buffer(gshp.points, shp->points);
        update_element_buffer(gshp.lines, shp->lines);
        update_element_buffer(gshp.triangles, shp->triangles);
        update_element_buffer(
            gshp.quads, convert_quads_to_triangles(shp->quads));
        update_element_buffer(
            gshp.beziers, convert_bezier_to_lines(shp->beziers));
        update_element_buffer(
            gshp.edges, get_edges({}, shp->triangles, shp->quads));
    }
}

// Init shading
std::unordered_map<shape*, gl_shape> make_gl_shapes(const scene* scn) {
    auto gshapes = std::unordered_map<shape*, gl_shape>();
    for (auto sgr : scn->shapes)
        for (auto shp : sgr->shapes) update_gl_shape(gshapes, shp);
    return gshapes;
}

// Clear scene textures on the GPU.
void clear_gl_textures(std::unordered_map<texture*, gl_texture>& gtextures) {
    for (auto& kv : gtextures) clear_texture(kv.second);
    gtextures.clear();
}

// Clear scene shapes on the GPU.
void clear_gl_shapes(std::unordered_map<shape*, gl_shape>& gshapes) {
    for (auto& kv : gshapes) clear_shape(kv.second);
    gshapes.clear();
}

// Initialize gl lights
gl_lights make_gl_lights(const scene* scn) {
    auto lights = gl_lights();
    for (auto ist : scn->instances) {
        if (!ist->shp) continue;
        for (auto shp : ist->shp->shapes) {
            if (!shp->mat) continue;
            if (shp->mat->ke == zero3f) continue;
            if (lights.pos.size() >= 16) break;
            if (!shp->points.empty()) {
                for (auto p : shp->points) {
                    if (lights.pos.size() >= 16) break;
                    lights.pos.push_back(
                        transform_point(ist->frame, shp->pos[p]));
                    lights.ke.push_back(shp->mat->ke);
                    lights.type.push_back(gl_light_type::point);
                }
            } else {
                auto bbox = make_bbox(shp->pos.size(), shp->pos.data());
                auto pos = bbox_center(bbox);
                auto area = 0.0f;
                for (auto l : shp->lines)
                    area += line_length(shp->pos[l.x], shp->pos[l.y]);
                for (auto t : shp->triangles)
                    area += triangle_area(
                        shp->pos[t.x], shp->pos[t.y], shp->pos[t.z]);
                for (auto t : shp->quads)
                    area += quad_area(shp->pos[t.x], shp->pos[t.y],
                        shp->pos[t.z], shp->pos[t.w]);
                auto ke = shp->mat->ke * area;
                if (lights.pos.size() < 16) {
                    lights.pos.push_back(transform_point(ist->frame, pos));
                    lights.ke.push_back(ke);
                    lights.type.push_back(gl_light_type::point);
                }
            }
        }
    }
    return lights;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR OPRNGL IMAGE SHADER FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Initialize the program. Call with true only after the GL is
// initialized.
gl_stdimage_program make_stdimage_program() {
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
            bool filmic;
            float exposure;
            float gamma;
        };
        uniform Tonemap tonemap;

        vec3 eval_filmic(vec3 x) {
            float a = 2.51f;
            float b = 0.03f;
            float c = 2.43f;
            float d = 0.59f;
            float e = 0.14f;
            return clamp((x*(a*x+b))/(x*(c*x+d)+e),0,1);
        }

        vec3 eval_tonemap(vec3 c) {
            // final color correction
            c = c*pow(2,tonemap.exposure);
            if(tonemap.filmic) {
                c = eval_filmic(c);
            } else {
                c = pow(c,vec3(1/tonemap.gamma));
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

    auto prog = gl_stdimage_program();
    prog.prog =
        make_program(_header + _vert, _header + _frag_tonemap + _frag_main);

    prog.vbo = make_vertex_buffer(
        std::vector<vec2f>{{0, 0}, {0, 1}, {1, 1}, {1, 0}}, false);
    prog.ebo =
        make_element_buffer(std::vector<vec3i>{{0, 1, 2}, {0, 2, 3}}, false);
    return prog;
}

// Draws the stdimage program.
void draw_image(gl_stdimage_program& prog, const gl_texture& txt,
    const vec2i& win_size, const vec2f& offset, float zoom, float exposure,
    float gamma, bool filmic) {
    assert(is_texture_valid(txt));

    bind_program(prog.prog);

    gl_enable_blending(true);
    gl_set_blend_over();

    bind_texture(txt, 0);
    set_program_uniform(prog.prog, "zoom", zoom);
    set_program_uniform(
        prog.prog, "win_size", vec2f{(float)win_size.x, (float)win_size.y});
    set_program_uniform(prog.prog, "offset", offset);
    set_program_uniform(prog.prog, "tonemap.filmic", filmic);
    set_program_uniform(prog.prog, "tonemap.exposure", exposure);
    set_program_uniform(prog.prog, "tonemap.gamma", gamma);
    set_program_uniform_texture(prog.prog, "img", txt, 0);

    set_program_vertattr(prog.prog, "vert_texcoord", prog.vbo, vec2f{0, 0});
    draw_elems(prog.ebo);

    unbind_program(prog.prog);

    gl_enable_blending(false);

    assert(gl_check_error());
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR OPRNGL STANDARD SURFACE SHADER FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Initialize a standard shader. Call with true only after the gl has
// been initialized
gl_stdsurface_program make_stdsurface_program() {
#ifndef _WIN32
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverlength-strings"
#endif
    std::string _vert_header =
        R"(
        #version 330

        )";

    std::string _vert_skinning =
        R"(
        #define SKIN_NONE 0
        #define SKIN_STD 1
        #define SKIN_GLTF 2

        uniform int skin_type = 0;
        uniform mat4 skin_xforms[32];
        layout(location = 6) in vec4 vert_skin_weights;            // vertex skinning weights
        layout(location = 7) in ivec4 vert_skin_joints;            // vertex skinning joints (in mesh coordinate frame)

        vec3 transform_point(mat4 m, vec3 p) {
            vec4 p4 = m * vec4(p,1);
            return p4.xyz / p4.w;
        }

        vec3 transform_normal(mat4 m, vec3 p) {
            vec4 p4 = m * vec4(p,0);
            return p4.xyz;
        }

        void apply_skin(inout vec3 pos, inout vec3 norm) {
            if(skin_type == 0) {
                return;
            } else if(skin_type == SKIN_STD) {
                vec4 w = vert_skin_weights;
                ivec4 j = ivec4(vert_skin_joints);
                pos = transform_point( skin_xforms[j.x], pos ) * w.x +
                transform_point( skin_xforms[j.y], pos ) * w.y +
                transform_point( skin_xforms[j.z], pos ) * w.z +
                transform_point( skin_xforms[j.w], pos ) * w.w;
                norm = normalize(
                                 transform_normal( skin_xforms[j.x], norm ) * w.x +
                                 transform_normal( skin_xforms[j.y], norm ) * w.y +
                                 transform_normal( skin_xforms[j.z], norm ) * w.z +
                                 transform_normal( skin_xforms[j.w], norm ) * w.w);
            } else if(skin_type == SKIN_GLTF) {
                vec4 w = vert_skin_weights;
                ivec4 j = ivec4(vert_skin_joints);
                mat4 xf = skin_xforms[j.x] * w.x + skin_xforms[j.y] * w.y +
                skin_xforms[j.z] * w.z + skin_xforms[j.w] * w.w;
                pos = transform_point(xf, pos);
                norm = normalize(transform_normal(xf, norm));
            }
        }
        )";

    std::string _vert_main =
        R"(
        layout(location = 0) in vec3 vert_pos;            // vertex position (in mesh coordinate frame)
        layout(location = 1) in vec3 vert_norm;           // vertex normal (in mesh coordinate frame)
        layout(location = 2) in vec2 vert_texcoord;       // vertex texcoords
        layout(location = 3) in vec4 vert_color;          // vertex color
        layout(location = 4) in vec4 vert_tangsp;         // vertex tangent space

        uniform mat4 shape_xform;           // shape transform
        uniform float shape_normal_offset;           // shape normal offset

        struct Camera {
            mat4 xform;          // camera xform
            mat4 xform_inv;      // inverse of the camera frame (as a matrix)
            mat4 proj;           // camera projection
        };
        uniform Camera camera;      // camera data

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

            // skinning
            apply_skin(pos, norm);

            // copy other vertex properties
            texcoord = vert_texcoord;
            color = vert_color;

            // clip
            gl_Position = camera.proj * camera.xform_inv * vec4(pos,1);
        }
        )";

    std::string _frag_header =
        R"(
        #version 330

        float pi = 3.14159265;

        )";

    std::string _frag_tonemap =
        R"(
        struct Tonemap {
            bool filmic;       // tonemap type (TM_...)
            float exposure; // image exposure
            float gamma;    // image gamma
        };
        uniform Tonemap tonemap;

        vec3 eval_filmic(vec3 x) {
            float a = 2.51f;
            float b = 0.03f;
            float c = 2.43f;
            float d = 0.59f;
            float e = 0.14f;
            return clamp((x*(a*x+b))/(x*(c*x+d)+e),0,1);
        }

        vec3 eval_tonemap(vec3 c) {
            // final color correction
            c = c*pow(2,tonemap.exposure);
            if(tonemap.filmic) {
                c = eval_filmic(c);
            } else {
                c = pow(c,vec3(1/tonemap.gamma));
            }
            return c;
        }

        )";

    std::string _frag_lighting =
        R"(
        struct Lighting {
            bool eyelight;        // eyelight shading
            vec3 amb;             // ambient light
            int lnum;              // number of lights
            int ltype[16];         // light type (0 -> point, 1 -> directional)
            vec3 lpos[16];         // light positions
            vec3 lke[16];          // light intensities
        };
        uniform Lighting lighting;

        void eval_light(int lid, vec3 pos, out vec3 cl, out vec3 wi) {
            cl = vec3(0,0,0);
            wi = vec3(0,0,0);
            if(lighting.ltype[lid] == 0) {
                // compute point light color at pos
                cl = lighting.lke[lid] / pow(length(lighting.lpos[lid]-pos),2);
                // compute light direction at pos
                wi = normalize(lighting.lpos[lid]-pos);
            }
            else if(lighting.ltype[lid] == 1) {
                // compute light color
                cl = lighting.lke[lid];
                // compute light direction
                wi = normalize(lighting.lpos[lid]);
            }
        }

        )";

    std::string _frag_brdf =
        R"(
        struct Brdf {
            int type;
            vec3 ke;
            vec3 kd;
            vec3 ks;
            float rs;
            float op;
            bool cutout;
        };

        vec3 brdfcos(Brdf brdf, vec3 n, vec3 wi, vec3 wo) {
            if(brdf.type == 0) return vec3(0);
            vec3 wh = normalize(wi+wo);
            float ns = 2/(brdf.rs*brdf.rs)-2;
            float ndi = dot(wi,n), ndo = dot(wo,n), ndh = dot(wh,n);
            if(brdf.type == 1) {
                return ((1+dot(wo,wi))/2) * brdf.kd/pi;
            } else if(brdf.type == 2) {
                float si = sqrt(1-ndi*ndi);
                float so = sqrt(1-ndo*ndo);
                float sh = sqrt(1-ndh*ndh);
                if(si <= 0) return vec3(0);
                vec3 diff = si * brdf.kd / pi;
                if(sh<=0) return diff;
                float d = ((2+ns)/(2*pi)) * pow(si,ns);
                vec3 spec = si * brdf.ks * d / (4*si*so);
                return diff+spec;
            } else if(brdf.type == 3 || brdf.type == 4) {
                if(ndi<=0 || ndo <=0) return vec3(0);
                vec3 diff = ndi * brdf.kd / pi;
                if(ndh<=0) return diff;
                if(brdf.type == 4) {
                    float d = ((2+ns)/(2*pi)) * pow(ndh,ns);
                    vec3 spec = ndi * brdf.ks * d / (4*ndi*ndo);
                    return diff+spec;
                } else {
                    float cos2 = ndh * ndh;
                    float tan2 = (1 - cos2) / cos2;
                    float alpha2 = brdf.rs * brdf.rs;
                    float d = alpha2 / (pi * cos2 * cos2 * (alpha2 + tan2) * (alpha2 + tan2));
                    float lambda_o = (-1 + sqrt(1 + (1 - ndo * ndo) / (ndo * ndo))) / 2;
                    float lambda_i = (-1 + sqrt(1 + (1 - ndi * ndi) / (ndi * ndi))) / 2;
                    float g = 1 / (1 + lambda_o + lambda_i);
                    vec3 spec = ndi * brdf.ks * d * g / (4*ndi*ndo);
                    return diff+spec;
                }
            }
        }

        )";

    std::string _frag_material =
        R"(
        struct Material {
            int type;         // material type
            int etype;         // element type
            vec3 ke;           // material ke
            vec3 kd;           // material kd
            vec3 ks;           // material ks
            float rs;          // material rs
            float op;          // material op

            bool txt_ke_on;    // material ke texture on
            sampler2D txt_ke;  // material ke texture
            bool txt_kd_on;    // material kd texture on
            sampler2D txt_kd;  // material kd texture
            bool txt_ks_on;    // material ks texture on
            sampler2D txt_ks;  // material ks texture
            bool txt_rs_on;    // material rs texture on
            sampler2D txt_rs;  // material rs texture

            bool txt_norm_on;    // material norm texture on
            sampler2D txt_norm;  // material norm texture
            sampler2D txt_norm_scale;  // material norm scale

            bool txt_occ_on;    // material occ texture on
            sampler2D txt_occ;  // material occ texture
            sampler2D txt_occ_scale;  // material occ scale

            bool use_phong;       // material use phong
            bool double_sided;    // material double sided
            bool alpha_cutout;    // material alpha cutout
        };
        uniform Material material;

        void eval_material(vec2 texcoord, vec4 color, out int type, out vec3 ke,
                           out vec3 kd, out vec3 ks, out float rs, out float op, out bool cutout) {
            if(material.type == 0) {
                type = 0;
                ke = material.ke;
                kd = vec3(0,0,0);
                ks = vec3(0,0,0);
                op = 1;
            }

            ke = color.xyz * material.ke;
            kd = color.xyz * material.kd;
            ks = color.xyz * material.ks;
            rs = material.rs;
            op = color.w * material.op;

            vec3 ke_txt = (material.txt_ke_on) ? texture(material.txt_ke,texcoord).xyz : vec3(1);
            vec4 kd_txt = (material.txt_kd_on) ? texture(material.txt_kd,texcoord) : vec4(1);
            vec4 ks_txt = (material.txt_ks_on) ? texture(material.txt_ks,texcoord) : vec4(1);
            float rs_txt = (material.txt_rs_on) ? texture(material.txt_rs,texcoord).x : 1;

            // scale common values
            ke *= ke_txt;

            // get material color from textures and adjust values
            if(material.type == 0) {
                type = 0;
            } else if(material.type == 1) {
                type = material.etype;
                kd *= kd_txt.xyz;
                ks *= ks_txt.xyz;
                rs *= rs_txt;
                rs = rs*rs;
            } else if(material.type == 2) {
                type = material.etype;
                vec3 kb = kd * kd_txt.xyz;
                float km = ks.x * ks_txt.z;
                kd = kb * (1 - km);
                ks = kb * km + vec3(0.04) * (1 - km);
                rs *= ks_txt.y;
                rs = rs*rs;
                op *= kd_txt.w;
            } else if(material.type == 3) {
                type = material.etype;
                kd *= kd_txt.xyz;
                ks *= ks_txt.xyz;
                rs *= ks_txt.w;
                rs = (1 - rs) * (1 - rs);
                op *= kd_txt.w;
            }

            cutout = material.alpha_cutout && op == 0;
        }

        vec3 apply_normal_map(vec2 texcoord, vec3 norm, vec4 tangsp) {
            if(!material.txt_norm_on) return norm;
            vec3 tangu = normalize(tangsp.xyz);
            vec3 tangv = normalize(cross(tangu, norm));
            if(tangsp.w < 0) tangv = -tangv;
            vec3 txt = 2 * pow(texture(material.txt_norm,texcoord).xyz, vec3(1/2.2)) - 1;
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

        struct Camera {
            mat4 xform;          // camera xform
            mat4 xform_inv;      // inverse of the camera frame (as a matrix)
            mat4 proj;           // camera projection
        };
        uniform Camera camera;      // camera data

        uniform vec4 highlight;   // highlighted color

        out vec4 frag_color;        // eyelight shading

        // main
        void main() {
            // view vector
            vec3 wo = normalize( (camera.xform*vec4(0,0,0,1)).xyz - pos );

            // re-normalize normals
            vec3 n = normalize(norm);

            // apply normal map
            n = apply_normal_map(texcoord, n, tangsp);

            // use faceforward to ensure the normals points toward us
            if(material.double_sided) n = faceforward(n,-wo,n);

            // get material color from textures
            Brdf brdf;
            eval_material(texcoord, color, brdf.type, brdf.ke, brdf.kd, brdf.ks, brdf.rs, brdf.op, brdf.cutout);

            // exit if needed
            if(brdf.cutout) discard;

            // check const color
            if(brdf.type == 0) {
                frag_color = vec4(brdf.ke,brdf.op);
                return;
            }

            // emission
            vec3 c = brdf.ke;

            // check early exit
            if(brdf.kd != vec3(0,0,0) || brdf.ks != vec3(0,0,0)) {
                // eyelight shading
                if(lighting.eyelight) {
                    vec3 wi = wo;
                    c += pi * brdfcos(brdf,n,wi,wo);
                } else {
                    // accumulate ambient
                    c += lighting.amb * brdf.kd;
                    // foreach light
                    for(int lid = 0; lid < lighting.lnum; lid ++) {
                        vec3 cl = vec3(0,0,0); vec3 wi = vec3(0,0,0);
                        eval_light(lid, pos, cl, wi);
                        c += cl * brdfcos(brdf,n,wi,wo);
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
            frag_color = vec4(c,brdf.op);
        }
        )";
#ifndef _WIN32
#pragma GCC diagnostic pop
#endif

    assert(gl_check_error());
    auto prog = gl_stdsurface_program();
    prog.prog = make_program(_vert_header + _vert_skinning + _vert_main,
        _frag_header + _frag_tonemap + _frag_lighting + _frag_brdf +
            _frag_material + _frag_main);
    assert(gl_check_error());
    return prog;
}

// Starts a frame by setting exposure/gamma values, camera transforms
// and projection. Sets also whether to use full shading or a quick
// eyelight preview.
void begin_stdsurface_frame(gl_stdsurface_program& prog, bool shade_eyelight,
    float tonemap_exposure, float tonemap_gamma, bool tonemap_filmic,
    const mat4f& camera_xform, const mat4f& camera_xform_inv,
    const mat4f& camera_proj) {
    static auto eyelight_id =
        get_program_uniform_location(prog.prog, "lighting.eyelight");
    static auto exposure_id =
        get_program_uniform_location(prog.prog, "tonemap.exposure");
    static auto gamma_id =
        get_program_uniform_location(prog.prog, "tonemap.gamma");
    static auto filmic_id =
        get_program_uniform_location(prog.prog, "tonemap.filmic");
    static auto xform_id =
        get_program_uniform_location(prog.prog, "camera.xform");
    static auto xform_inv_id =
        get_program_uniform_location(prog.prog, "camera.xform_inv");
    static auto proj_id =
        get_program_uniform_location(prog.prog, "camera.proj");
    assert(gl_check_error());
    bind_program(prog.prog);
    set_program_uniform(prog.prog, eyelight_id, shade_eyelight);
    set_program_uniform(prog.prog, exposure_id, tonemap_exposure);
    set_program_uniform(prog.prog, gamma_id, tonemap_gamma);
    set_program_uniform(prog.prog, filmic_id, tonemap_filmic);
    set_program_uniform(prog.prog, xform_id, camera_xform);
    set_program_uniform(prog.prog, xform_inv_id, camera_xform_inv);
    set_program_uniform(prog.prog, proj_id, camera_proj);
    assert(gl_check_error());
}

// Ends a frame.
void end_stdsurface_frame(gl_stdsurface_program& prog) {
    assert(gl_check_error());
    unbind_program(prog.prog);
    //    glBindVertexArray(0);
    //    glUseProgram(0);
    assert(gl_check_error());
}

// Set num lights with position pos, color ke, type ltype. Also set the
// ambient illumination amb.
void set_stdsurface_lights(
    gl_stdsurface_program& prog, const vec3f& amb, const gl_lights& lights) {
    static auto amb_id =
        get_program_uniform_location(prog.prog, "lighting.amb");
    static auto lnum_id =
        get_program_uniform_location(prog.prog, "lighting.lnum");
    static auto lpos_id =
        get_program_uniform_location(prog.prog, "lighting.lpos");
    static auto lke_id =
        get_program_uniform_location(prog.prog, "lighting.lke");
    static auto ltype_id =
        get_program_uniform_location(prog.prog, "lighting.ltype");
    assert(gl_check_error());
    set_program_uniform(prog.prog, amb_id, amb);
    set_program_uniform(prog.prog, lnum_id, (int)lights.pos.size());
    set_program_uniform(
        prog.prog, lpos_id, lights.pos.data(), (int)lights.pos.size());
    set_program_uniform(
        prog.prog, lke_id, lights.ke.data(), (int)lights.pos.size());
    set_program_uniform(
        prog.prog, ltype_id, (int*)lights.type.data(), (int)lights.pos.size());
    assert(gl_check_error());
}

// Begins drawing a shape with transform xform.
void begin_stdsurface_shape(
    gl_stdsurface_program& prog, const mat4f& xform, float normal_offset) {
    static auto xform_id =
        get_program_uniform_location(prog.prog, "shape_xform");
    static auto normal_offset_id =
        get_program_uniform_location(prog.prog, "shape_normal_offset");
    assert(gl_check_error());
    set_program_uniform(prog.prog, xform_id, xform);
    set_program_uniform(prog.prog, normal_offset_id, normal_offset);
    assert(gl_check_error());
}

// End shade drawing.
void end_stdsurface_shape(gl_stdsurface_program& prog) {
    assert(gl_check_error());
    for (int i = 0; i < 16; i++) unbind_vertex_buffer(i);
    assert(gl_check_error());
}

// Sets normal offset.
void set_stdsurface_normaloffset(
    gl_stdsurface_program& prog, float normal_offset) {
    static auto normal_offset_id =
        get_program_uniform_location(prog.prog, "shape_normal_offset");
    assert(gl_check_error());
    set_program_uniform(prog.prog, normal_offset_id, normal_offset);
    assert(gl_check_error());
}

// Set the object as highlighted.
void set_stdsurface_highlight(
    gl_stdsurface_program& prog, const vec4f& highlight) {
    static auto highlight_id =
        get_program_uniform_location(prog.prog, "highlight");
    set_program_uniform(prog.prog, highlight_id, highlight);
}

// Set material values with emission ke, diffuse kd, specular ks and
// specular roughness rs, opacity op. Indicates textures ids with the
// correspoinding XXX_txt variables. Sets also normal and occlusion
// maps. Works for points/lines/triangles (diffuse for points,
// Kajiya-Kay for lines, GGX/Phong for triangles).
// Material type matches the scene material type.
void set_stdsurface_material(gl_stdsurface_program& prog, material_type type,
    gl_elem_type etype, const vec3f& ke, const vec3f& kd, const vec3f& ks,
    float rs, float op, const gl_texture_info& ke_txt,
    const gl_texture_info& kd_txt, const gl_texture_info& ks_txt,
    const gl_texture_info& rs_txt, const gl_texture_info& norm_txt,
    const gl_texture_info& occ_txt, bool use_phong, bool double_sided,
    bool alpha_cutout) {
    static auto mtype_id =
        get_program_uniform_location(prog.prog, "material.type");
    static auto etype_id =
        get_program_uniform_location(prog.prog, "material.etype");
    static auto ke_id = get_program_uniform_location(prog.prog, "material.ke");
    static auto kd_id = get_program_uniform_location(prog.prog, "material.kd");
    static auto ks_id = get_program_uniform_location(prog.prog, "material.ks");
    static auto rs_id = get_program_uniform_location(prog.prog, "material.rs");
    static auto op_id = get_program_uniform_location(prog.prog, "material.op");
    static auto ke_txt_id =
        get_program_uniform_location(prog.prog, "material.txt_ke");
    static auto ke_txt_on_id =
        get_program_uniform_location(prog.prog, "material.txt_ke_on");
    static auto kd_txt_id =
        get_program_uniform_location(prog.prog, "material.txt_kd");
    static auto kd_txt_on_id =
        get_program_uniform_location(prog.prog, "material.txt_kd_on");
    static auto ks_txt_id =
        get_program_uniform_location(prog.prog, "material.txt_ks");
    static auto ks_txt_on_id =
        get_program_uniform_location(prog.prog, "material.txt_ks_on");
    static auto rs_txt_id =
        get_program_uniform_location(prog.prog, "material.txt_rs");
    static auto rs_txt_on_id =
        get_program_uniform_location(prog.prog, "material.txt_rs_on");
    static auto norm_txt_id =
        get_program_uniform_location(prog.prog, "material.txt_norm");
    static auto norm_txt_on_id =
        get_program_uniform_location(prog.prog, "material.txt_norm_on");
    static auto occ_txt_id =
        get_program_uniform_location(prog.prog, "material.txt_occ");
    static auto occ_txt_on_id =
        get_program_uniform_location(prog.prog, "material.txt_occ_on");
    static auto norm_scale_id =
        get_program_uniform_location(prog.prog, "material.norm_scale");
    static auto occ_scale_id =
        get_program_uniform_location(prog.prog, "material.occ_scale");
    static auto use_phong_id =
        get_program_uniform_location(prog.prog, "material.use_phong");
    static auto double_sided_id =
        get_program_uniform_location(prog.prog, "material.double_sided");
    static auto alpha_cutout_id =
        get_program_uniform_location(prog.prog, "material.alpha_cutout");

    static auto mtypes = std::unordered_map<material_type, int>{
        {material_type::specular_roughness, 1},
        {material_type::metallic_roughness, 2},
        {material_type::specular_glossiness, 3}};

    assert(gl_check_error());
    set_program_uniform(prog.prog, mtype_id, mtypes.at(type));
    set_program_uniform(prog.prog, etype_id, (int)etype);
    set_program_uniform(prog.prog, ke_id, ke);
    set_program_uniform(prog.prog, kd_id, kd);
    set_program_uniform(prog.prog, ks_id, ks);
    set_program_uniform(prog.prog, rs_id, rs);
    set_program_uniform(prog.prog, op_id, op);
    set_program_uniform_texture(prog.prog, ke_txt_id, ke_txt_on_id, ke_txt, 0);
    set_program_uniform_texture(prog.prog, kd_txt_id, kd_txt_on_id, kd_txt, 1);
    set_program_uniform_texture(prog.prog, ks_txt_id, ks_txt_on_id, ks_txt, 2);
    set_program_uniform_texture(prog.prog, rs_txt_id, rs_txt_on_id, rs_txt, 3);
    set_program_uniform_texture(
        prog.prog, norm_txt_id, norm_txt_on_id, norm_txt, 4);
    set_program_uniform_texture(
        prog.prog, occ_txt_id, occ_txt_on_id, occ_txt, 5);
    set_program_uniform(prog.prog, norm_scale_id, norm_txt.scale);
    set_program_uniform(prog.prog, occ_scale_id, occ_txt.scale);
    set_program_uniform(prog.prog, use_phong_id, use_phong);
    set_program_uniform(prog.prog, double_sided_id, double_sided);
    set_program_uniform(prog.prog, alpha_cutout_id, alpha_cutout);
    assert(gl_check_error());
}

// Set constant material values with emission ke.
void set_stdsurface_constmaterial(
    gl_stdsurface_program& prog, const vec3f& ke, float op) {
    static auto mtype_id =
        get_program_uniform_location(prog.prog, "material.type");
    static auto etype_id =
        get_program_uniform_location(prog.prog, "material.etype");
    static auto ke_id = get_program_uniform_location(prog.prog, "material.ke");
    static auto op_id = get_program_uniform_location(prog.prog, "material.op");

    assert(gl_check_error());
    set_program_uniform(prog.prog, mtype_id, 0);
    set_program_uniform(prog.prog, etype_id, 0);
    set_program_uniform(prog.prog, ke_id, ke);
    set_program_uniform(prog.prog, op_id, op);
    assert(gl_check_error());
}

// Set vertex data with buffers for position pos, normals norm, texture
// coordinates texcoord, per-vertex color color and tangent space
// tangsp.
void set_stdsurface_vert(gl_stdsurface_program& prog,
    const gl_vertex_buffer& pos, const gl_vertex_buffer& norm,
    const gl_vertex_buffer& texcoord, const gl_vertex_buffer& color,
    const gl_vertex_buffer& tangsp) {
    static auto pos_id = get_program_attrib_location(prog.prog, "vert_pos");
    static auto norm_id = get_program_attrib_location(prog.prog, "vert_norm");
    static auto texcoord_id =
        get_program_attrib_location(prog.prog, "vert_texcoord");
    static auto color_id = get_program_attrib_location(prog.prog, "vert_color");
    static auto tangsp_id =
        get_program_attrib_location(prog.prog, "vert_tangsp");
    assert(gl_check_error());
    set_program_vertattr(prog.prog, pos_id, pos, zero3f);
    set_program_vertattr(prog.prog, norm_id, norm, zero3f);
    set_program_vertattr(prog.prog, texcoord_id, texcoord, zero2f);
    set_program_vertattr(prog.prog, color_id, color, vec4f{1, 1, 1, 1});
    set_program_vertattr(prog.prog, tangsp_id, tangsp, zero4f);
    assert(gl_check_error());
}

// Set vertex data with buffers for skinning.
void set_stdsurface_vert_skinning(gl_stdsurface_program& prog,
    const gl_vertex_buffer& weights, const gl_vertex_buffer& joints,
    int nxforms, const mat4f* xforms) {
    static auto type_id = get_program_uniform_location(prog.prog, "skin_type");
    static auto xforms_id =
        get_program_uniform_location(prog.prog, "skin_xforms");
    static auto weights_id =
        get_program_attrib_location(prog.prog, "vert_skin_weights");
    static auto joints_id =
        get_program_attrib_location(prog.prog, "vert_skin_joints");
    int type = 1;
    set_program_uniform(prog.prog, type_id, type);
    set_program_uniform(prog.prog, xforms_id, xforms, min(nxforms, 32));
    set_program_vertattr(prog.prog, weights_id, weights, zero4f);
    set_program_vertattr(prog.prog, joints_id, joints, zero4f);
}

// Set vertex data with buffers for skinning.
void set_stdsurface_vert_gltf_skinning(gl_stdsurface_program& prog,
    const gl_vertex_buffer& weights, const gl_vertex_buffer& joints,
    int nxforms, const mat4f* xforms) {
    static auto type_id = get_program_uniform_location(prog.prog, "skin_type");
    static auto xforms_id =
        get_program_uniform_location(prog.prog, "skin_xforms");
    static auto weights_id =
        get_program_attrib_location(prog.prog, "vert_skin_weights");
    static auto joints_id =
        get_program_attrib_location(prog.prog, "vert_skin_joints");
    int type = 2;
    set_program_uniform(prog.prog, type_id, type);
    set_program_uniform(prog.prog, xforms_id, xforms, min(nxforms, 32));
    set_program_vertattr(prog.prog, weights_id, weights, zero4f);
    set_program_vertattr(prog.prog, joints_id, joints, zero4f);
}

// Disables vertex skinning.
void set_stdsurface_vert_skinning_off(gl_stdsurface_program& prog) {
    static auto type_id = get_program_uniform_location(prog.prog, "skin_type");
    // static auto xforms_id = get_program_uniform_location(prog.prog,
    // "skin_xforms");
    static auto weights_id =
        get_program_attrib_location(prog.prog, "vert_skin_weights");
    static auto joints_id =
        get_program_attrib_location(prog.prog, "vert_skin_joints");
    int type = 0;
    set_program_uniform(prog.prog, type_id, type);
    set_program_vertattr(prog.prog, weights_id, {}, zero4f);
    set_program_vertattr(prog.prog, joints_id, {}, zero4f);
}

// Draw a shape
void draw_stdsurface_shape(const shape* shp, const mat4f& xform,
    bool highlighted, gl_stdsurface_program& prog,
    std::unordered_map<shape*, gl_shape>& shapes,
    std::unordered_map<texture*, gl_texture>& textures,
    const gl_stdsurface_params& params) {
    static auto default_material = material();
    default_material.kd = {0.2f, 0.2f, 0.2f};

    begin_stdsurface_shape(prog, xform);

    auto etype = gl_elem_type::triangle;
    if (!shp->lines.empty()) etype = gl_elem_type::line;
    if (!shp->points.empty()) etype = gl_elem_type::point;

    auto txt = [&textures](texture* txt) -> gl_texture_info {
        if (!txt) return {};
        return textures.at(txt);
    };

    auto mat = (shp->mat) ? shp->mat : &default_material;
    set_stdsurface_material(prog, mat->type, etype, mat->ke, mat->kd, mat->ks,
        mat->rs, mat->op, txt(mat->ke_txt), txt(mat->kd_txt), txt(mat->ks_txt),
        txt(mat->rs_txt), txt(mat->norm_txt), txt(mat->occ_txt), false,
        mat->double_sided, params.cutout);

    auto& gshp = shapes.at((shape*)shp);
    set_stdsurface_vert(
        prog, gshp.pos, gshp.norm, gshp.texcoord, gshp.color, gshp.tangsp);

    draw_elems(gshp.points);
    draw_elems(gshp.lines);
    draw_elems(gshp.triangles);
    draw_elems(gshp.quads);
    draw_elems(gshp.beziers);

    if ((params.edges && !params.wireframe) || highlighted) {
        gl_enable_culling(false);
        set_stdsurface_constmaterial(prog,
            (highlighted) ? params.highlight_color : params.edge_color,
            (highlighted) ? 1 : mat->op);
        set_stdsurface_normaloffset(prog, params.edge_offset);
        draw_elems(gshp.edges);
        gl_enable_culling(params.cull_backface);
    }

    if (highlighted && false) {
        set_stdsurface_constmaterial(prog, params.highlight_color, 1);
        set_stdsurface_normaloffset(prog, params.edge_offset);
        draw_elems(gshp.points);
        draw_elems(gshp.lines);
        draw_elems(gshp.triangles);
        draw_elems(gshp.quads);
        draw_elems(gshp.beziers);
    }

    end_stdsurface_shape(prog);
}

// Display a scene
void draw_stdsurface_scene(const scene* scn, const camera* cam,
    gl_stdsurface_program& prog, std::unordered_map<shape*, gl_shape>& shapes,
    std::unordered_map<texture*, gl_texture>& textures, const gl_lights& lights,
    const vec2i& viewport_size, void* highlighted,
    const gl_stdsurface_params& params, const tonemap_params& tmparams) {
    // begin frame
    gl_enable_depth_test(true);
    gl_enable_culling(params.cull_backface);
    gl_enable_wireframe(params.wireframe);
    gl_set_viewport(viewport_size);

    auto camera_xform = frame_to_mat(cam->frame);
    auto camera_view = frame_to_mat(inverse(cam->frame));
    auto camera_proj = perspective_mat(cam->yfov,
        (float)viewport_size.x / (float)viewport_size.y, cam->near, cam->far);

    begin_stdsurface_frame(prog, params.eyelight, tmparams.exposure,
        tmparams.gamma, tmparams.filmic, camera_xform, camera_view,
        camera_proj);

    if (!params.eyelight) {
        set_stdsurface_lights(prog, params.ambient, lights);
    }

    if (!scn->instances.empty()) {
        for (auto ist : scn->instances) {
            for (auto shp : ist->shp->shapes) {
                draw_stdsurface_shape(shp, frame_to_mat(ist->frame),
                    ist == highlighted || ist->shp == highlighted ||
                        shp == highlighted,
                    prog, shapes, textures, params);
            }
        }
    } else {
        for (auto sgr : scn->shapes) {
            for (auto shp : sgr->shapes) {
                draw_stdsurface_shape(shp, identity_mat4f,
                    sgr == highlighted || shp == highlighted, prog, shapes,
                    textures, params);
            }
        }
    }

    end_stdsurface_frame(prog);
    gl_enable_wireframe(false);
}

// Support
void _glfw_error_cb(int error, const char* description) {
    printf("GLFW error: %s\n", description);
}

// Support
void _glfw_text_cb(GLFWwindow* gwin, unsigned key) {
    auto win = (gl_window*)glfwGetWindowUserPointer(gwin);
    if (win->widget_enabled) { ImGui_ImplGlfwGL3_CharCallback(win->gwin, key); }
    if (win->text_cb) win->text_cb(win, key);
}

// Support
void _glfw_key_cb(
    GLFWwindow* gwin, int key, int scancode, int action, int mods) {
    auto win = (gl_window*)glfwGetWindowUserPointer(gwin);
    if (win->widget_enabled) {
        ImGui_ImplGlfwGL3_KeyCallback(win->gwin, key, scancode, action, mods);
    }
}

// Support
void _glfw_mouse_cb(GLFWwindow* gwin, int button, int action, int mods) {
    auto win = (gl_window*)glfwGetWindowUserPointer(gwin);
    if (win->widget_enabled) {
        ImGui_ImplGlfwGL3_MouseButtonCallback(win->gwin, button, action, mods);
    }
    if (win->mouse_cb) win->mouse_cb(win, button, action == GLFW_PRESS, mods);
}

// Support
void _glfw_scroll_cb(GLFWwindow* gwin, double xoffset, double yoffset) {
    auto win = (gl_window*)glfwGetWindowUserPointer(gwin);
    if (win->widget_enabled) {
        ImGui_ImplGlfwGL3_ScrollCallback(win->gwin, xoffset, yoffset);
    }
}

// Support
void _glfw_refresh_cb(GLFWwindow* gwin) {
    auto win = (gl_window*)glfwGetWindowUserPointer(gwin);
    if (win->refresh_cb) win->refresh_cb(win);
}

// Initialize gl_window
gl_window* make_window(
    int width, int height, const std::string& title, void* user_pointer) {
    auto win = new gl_window();
    // gl_window
    win->user_pointer = user_pointer;

    // gl_window
    if (!glfwInit()) throw std::runtime_error("cannot open gl_window");

    // profile creation
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
#if __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    win->gwin = glfwCreateWindow(width, height, title.c_str(), 0, 0);
    glfwMakeContextCurrent(win->gwin);
    glfwSetWindowUserPointer(win->gwin, win);

    glfwSetErrorCallback(_glfw_error_cb);

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

// Set gl_window callbacks
void set_window_callbacks(gl_window* win, gl_text_callback text_cb,
    gl_mouse_callback mouse_cb, gl_refresh_callback refresh_cb) {
    win->text_cb = text_cb;
    win->mouse_cb = mouse_cb;
    win->refresh_cb = refresh_cb;
    if (win->text_cb) glfwSetCharCallback(win->gwin, _glfw_text_cb);
}

// Clear gl_window
void clear_window(gl_window* win) {
    if (win->gwin) {
        glfwDestroyWindow(win->gwin);
        glfwTerminate();
        win->gwin = nullptr;
    }
    if (win->widget_enabled) {
        ImGui_ImplGlfwGL3_Shutdown();
        win->widget_enabled = false;
    }
}

// Set gl_window title
void set_window_title(gl_window* win, const std::string& title) {
    glfwSetWindowTitle(win->gwin, title.c_str());
}

// Wait events
void wait_events(gl_window* win) { glfwWaitEvents(); }

// Poll events
void poll_events(gl_window* win) { glfwPollEvents(); }

// Swap buffers
void swap_buffers(gl_window* win) { glfwSwapBuffers(win->gwin); }

// Should close
bool should_close(gl_window* win) { return glfwWindowShouldClose(win->gwin); }

// Mouse button
int get_mouse_button(gl_window* win) {
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
vec2i get_mouse_pos(gl_window* win) {
    double x, y;
    glfwGetCursorPos(win->gwin, &x, &y);
    return {(int)x, (int)y};
}

// Mouse position
vec2f get_mouse_posf(gl_window* win) {
    double x, y;
    glfwGetCursorPos(win->gwin, &x, &y);
    return {(float)x, (float)y};
}

// Window size
vec2i get_window_size(gl_window* win) {
    auto ret = vec2i{0, 0};
    glfwGetWindowSize(win->gwin, &ret.x, &ret.y);
    return ret;
}

// Check if a key is pressed (not all keys are supported)
bool get_key(gl_window* win, int key) {
    key = std::toupper(key);
    return glfwGetKey(win->gwin, key) == GLFW_PRESS;
}

// Framebuffer size
vec2i get_framebuffer_size(gl_window* win) {
    auto ret = vec2i{0, 0};
    glfwGetFramebufferSize(win->gwin, &ret.x, &ret.y);
    return ret;
}

// Read pixels
image4b take_screenshot4b(gl_window* win, bool flipy, bool back) {
    auto wh = get_framebuffer_size(win);
    auto img = image4b(wh.x, wh.y);
    glReadBuffer((back) ? GL_BACK : GL_FRONT);
    glReadPixels(
        0, 0, img.width(), img.height(), GL_RGBA, GL_UNSIGNED_BYTE, data(img));
    if (flipy) {
        for (int j = 0; j < img.height() / 2; j++) {
            for (auto i = 0; i < img.width(); i++) {
                std::swap(img.at(i, j), img.at(i, img.height() - 1 - j));
            }
        }
    }
    return img;
}

// Handles camera navigation
bool handle_camera_navigation(
    gl_window* win, camera* cam, bool navigation_fps) {
    static auto mouse_last = zero2f;
    auto mouse_pos = get_mouse_posf(win);
    auto mouse_button = get_mouse_button(win);

    // updated
    auto updated = false;

    // handle mouse and keyboard for navigation
    if (mouse_button && !get_widget_active(win)) {
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
    if (!get_widget_active(win) && navigation_fps) {
        auto transl = zero3f;
        if (get_key(win, 'a')) transl.x -= 1;
        if (get_key(win, 'd')) transl.x += 1;
        if (get_key(win, 's')) transl.z += 1;
        if (get_key(win, 'w')) transl.z -= 1;
        if (get_key(win, 'e')) transl.y += 1;
        if (get_key(win, 'q')) transl.y -= 1;
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

// Initialize widgets
void init_widgets(gl_window* win, bool light_style, bool alt_font) {
    ImGui_ImplGlfwGL3_Init(win->gwin, false);
    ImGui::GetStyle().WindowRounding = 0;
    ImGui::GetIO().IniFilename = nullptr;
    auto size = get_window_size(win);
    ImGui::SetNextWindowPos({(float)size[0] - 320, 0});
    ImGui::SetNextWindowSize({(float)320, (float)size[1]});
    if (light_style) ImGui::StyleColorsLight();
    if (alt_font) {
        ImGuiIO& io = ImGui::GetIO();
        io.Fonts->AddFontFromMemoryCompressedTTF(
            imgui_extrafont_compressed_data(),
            imgui_extrafont_compressed_size(), 16);
    } else {
        ImGuiIO& io = ImGui::GetIO();
        io.Fonts->AddFontDefault();
    }
    win->widget_enabled = true;
}

// Begin draw widget
bool begin_widgets(gl_window* win, const std::string& title) {
    static bool first_time = true;
    ImGui_ImplGlfwGL3_NewFrame();
    if (first_time) {
        auto size = get_window_size(win);
        ImGui::SetNextWindowPos({(float)size[0] - 320, 0});
        ImGui::SetNextWindowSize({(float)320, (float)size[1]});
        ImGui::SetNextWindowCollapsed(true);
        first_time = false;
    }
    ImGui::Begin(title.c_str(), nullptr);
    // ImGui::ShowTestWindow();
    // ImGui::ShowStyleEditor();
    return true;
}

// End draw widget
void end_widgets(gl_window* win) {
    ImGui::End();
    ImGui::Render();
}

// Whether widget are active
bool get_widget_active(gl_window* win) {
    if (!win->widget_enabled) return false;
    auto io = &ImGui::GetIO();
    return io->WantTextInput || io->WantCaptureMouse || io->WantCaptureKeyboard;
}

// Horizontal separator
void draw_separator_widget(gl_window* win) { ImGui::Separator(); }

// Indent widget
void draw_indent_widget_begin(gl_window* win) { ImGui::Indent(); }

// Indent widget
void draw_indent_widget_end(gl_window* win) { ImGui::Unindent(); }

// Continue line with next widget
void draw_continue_widget(gl_window* win) { ImGui::SameLine(); }

// Label widget
void draw_label_widget(
    gl_window* win, const std::string& lbl, const std::string& msg) {
    ImGui::LabelText(lbl.c_str(), "%s", msg.c_str());
}

// Value widget
bool draw_text_widget(
    gl_window* win, const std::string& lbl, std::string& str) {
    char buf[4096];
    if (str.length() >= 4096) throw std::runtime_error("bad memory");
    memcpy(buf, str.c_str(), str.size());
    buf[str.size()] = 0;
    auto ret = ImGui::InputText(lbl.c_str(), buf, 4096);
    str = buf;
    return ret;
}

// Value widget
bool draw_slider_widget(
    gl_window* win, const std::string& lbl, int& val, int min, int max) {
    return ImGui::SliderInt(lbl.c_str(), &val, min, max);
}
// Value widget
bool draw_slider_widget(
    gl_window* win, const std::string& lbl, vec2i& val, int min, int max) {
    return ImGui::SliderInt2(lbl.c_str(), &val.x, min, max);
}
// Value widget
bool draw_slider_widget(
    gl_window* win, const std::string& lbl, vec3i& val, int min, int max) {
    return ImGui::SliderInt3(lbl.c_str(), &val.x, min, max);
}
// Value widget
bool draw_slider_widget(
    gl_window* win, const std::string& lbl, vec4i& val, int min, int max) {
    return ImGui::SliderInt4(lbl.c_str(), &val.x, min, max);
}

// Value widget
bool draw_slider_widget(
    gl_window* win, const std::string& lbl, float& val, float min, float max) {
    return ImGui::SliderFloat(lbl.c_str(), &val, min, max);
}
// Value widget
bool draw_slider_widget(
    gl_window* win, const std::string& lbl, vec2f& val, float min, float max) {
    return ImGui::SliderFloat2(lbl.c_str(), &val.x, min, max);
}
// Value widget
bool draw_slider_widget(
    gl_window* win, const std::string& lbl, vec3f& val, float min, float max) {
    return ImGui::SliderFloat3(lbl.c_str(), &val.x, min, max);
}
// Value widget
bool draw_slider_widget(
    gl_window* win, const std::string& lbl, vec4f& val, float min, float max) {
    return ImGui::SliderFloat4(lbl.c_str(), &val.x, min, max);
}

// Slider widget.
bool draw_slider_widget(gl_window* win, const std::string& lbl,
    mat<float, 4>& val, float min, float max) {
    auto modx = draw_slider_widget(win, lbl + ".x", val.x, min, max);
    auto mody = draw_slider_widget(win, lbl + ".y", val.y, min, max);
    auto modz = draw_slider_widget(win, lbl + ".z", val.z, min, max);
    auto modw = draw_slider_widget(win, lbl + ".w", val.w, min, max);
    return modx || mody || modz || modw;
}
// Slider widget.
bool draw_slider_widget(gl_window* win, const std::string& lbl,
    frame<float, 3>& val, float min, float max) {
    auto modx = draw_slider_widget(win, lbl + ".x", val.x, -1, 1);
    auto mody = draw_slider_widget(win, lbl + ".y", val.y, -1, 1);
    auto modz = draw_slider_widget(win, lbl + ".z", val.z, -1, 1);
    auto modo = draw_slider_widget(win, lbl + ".o", val.o, min, max);
    // TODO: orthonormalize
    return modx || mody || modz || modo;
}
// Slider widget
bool draw_slider_widget(
    gl_window* win, const std::string& lbl, quat4f& val, float min, float max) {
    auto mod = ImGui::SliderFloat4(lbl.c_str(), &val.x, min, max);
    if (mod) val = normalize(val);
    return mod;
}

// Color widget
bool draw_color_widget(gl_window* win, const std::string& lbl, vec4f& val) {
    return ImGui::ColorEdit4(lbl.c_str(), (float*)&val.x);
}
// Color widget
bool draw_color_widget(gl_window* win, const std::string& lbl, vec3f& val) {
    return ImGui::ColorEdit3(lbl.c_str(), (float*)&val.x);
}
// Color widget
bool draw_color_widget(gl_window* win, const std::string& lbl, vec4b& val) {
    auto valf = ImGui::ColorConvertU32ToFloat4(*(uint32_t*)&val);
    if (ImGui::ColorEdit4(lbl.c_str(), &valf.x)) {
        auto valb = ImGui::ColorConvertFloat4ToU32(valf);
        *(uint32_t*)&val = valb;
        return true;
    }
    return false;
}

// Bool widget
bool draw_checkbox_widget(gl_window* win, const std::string& lbl, bool& val) {
    return ImGui::Checkbox(lbl.c_str(), &val);
}

// Combo widget
bool draw_combo_widget_begin(
    gl_window* win, const std::string& lbl, const std::string& label) {
    return ImGui::BeginCombo(lbl.c_str(), label.c_str());
}

// Combo widget
bool draw_combo_widget_item(
    gl_window* win, const std::string& label, int idx, bool selected) {
    ImGui::PushID((void*)(intptr_t)idx);
    auto clicked = ImGui::Selectable(label.c_str(), selected);
    if (selected) ImGui::SetItemDefaultFocus();
    ImGui::PopID();
    return clicked;
}

// Combo widget
void draw_combo_widget_end(gl_window* win) { ImGui::EndCombo(); }

// Button widget
bool draw_button_widget(gl_window* win, const std::string& lbl) {
    return ImGui::Button(lbl.c_str());
}

// Collapsible header
bool draw_header_widget(gl_window* win, const std::string& lbl) {
    return ImGui::CollapsingHeader(lbl.c_str());
}

// Start tree node
bool draw_tree_widget_begin(gl_window* win, const std::string& lbl) {
    return ImGui::TreeNode(lbl.c_str());
}

// Collapsible header
void draw_tree_widget_end(gl_window* win) { ImGui::TreePop(); }

// Start selectable tree node
bool draw_tree_widget_begin(
    gl_window* win, const std::string& lbl, void*& selection, void* content) {
    ImGuiTreeNodeFlags node_flags =
        ImGuiTreeNodeFlags_OpenOnArrow | ImGuiTreeNodeFlags_OpenOnDoubleClick;
    if (selection == content) node_flags |= ImGuiTreeNodeFlags_Selected;
    auto open = ImGui::TreeNodeEx(content, node_flags, "%s", lbl.c_str());
    if (ImGui::IsItemClicked()) selection = content;
    return open;
}

// Start selectable tree node
bool draw_tree_widget_begin(gl_window* win, const std::string& lbl,
    void*& selection, void* content, const vec4f& col) {
    ImGui::PushStyleColor(ImGuiCol_Text, {col.x, col.y, col.z, col.w});
    auto ret = draw_tree_widget_begin(win, lbl, selection, content);
    ImGui::PopStyleColor();
    return ret;
}

// End selectable tree node
void draw_tree_widget_end(gl_window* win, void* content) { ImGui::TreePop(); }

// Selectable tree leaf node
void draw_tree_widget_leaf(
    gl_window* win, const std::string& lbl, void*& selection, void* content) {
    ImGuiTreeNodeFlags node_flags =
        ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen;
    if (selection == content) node_flags |= ImGuiTreeNodeFlags_Selected;
    ImGui::TreeNodeEx(content, node_flags, "%s", lbl.c_str());
    if (ImGui::IsItemClicked()) selection = content;
}

// Selectable tree leaf node
void draw_tree_widget_leaf(gl_window* win, const std::string& lbl,
    void*& selection, void* content, const vec4f& col) {
    ImGui::PushStyleColor(ImGuiCol_Text, {col.x, col.y, col.z, col.w});
    draw_tree_widget_leaf(win, lbl, selection, content);
    ImGui::PopStyleColor();
}

// Image widget
void draw_image_widget(
    gl_window* win, int tid, const vec2i& size, const vec2i& imsize) {
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
void draw_image_widget(gl_window* win, gl_texture& txt, const vec2i& size) {
    draw_image_widget(win, get_texture_id(txt), size, {txt.width, txt.height});
}

// Scroll region
void draw_scroll_widget_begin(
    gl_window* win, const std::string& lbl, int height, bool border) {
    ImGui::BeginChild(lbl.c_str(), ImVec2(0, height), border);
}
// Scroll region
void draw_scroll_widget_end(gl_window* win) { ImGui::EndChild(); }
// Scroll region
void draw_scroll_widget_here(gl_window* win) { ImGui::SetScrollHere(); }

// Group ids
void draw_groupid_widget_begin(gl_window* win, int gid) { ImGui::PushID(gid); }
// Group ids
void draw_groupid_widget_begin(gl_window* win, void* gid) {
    ImGui::PushID(gid);
}
// Group ids
void draw_groupid_widget_begin(gl_window* win, const char* gid) {
    ImGui::PushID(gid);
}
// Group ids
void draw_groupid_widget_end(gl_window* win) { ImGui::PopID(); }

// Text color
void draw_tree_widget_color_begin(gl_window* win, const vec4f& color) {
    ImGui::PushStyleColor(
        ImGuiCol_Text, {color[0], color[1], color[2], color[3]});
}

// Text color
void draw_tree_widget_color_end(gl_window* win) { ImGui::PopStyleColor(); }

// Image inspection widgets.
void draw_imageinspect_widgets(gl_window* win, const std::string& lbl,
    const image4f& hdr, const image4b& ldr, const vec2f& mouse_pos,
    const gl_stdimage_params& params) {
    auto xy = (mouse_pos - params.offset) / params.zoom;
    auto i = (int)round(xy.x), j = (int)round(xy.y);
    auto v4f = zero4f;
    auto v4b = zero4b;
    if (!hdr.empty()) {
        auto w = hdr.width(), h = hdr.height();
        if (i >= 0 && i < w && j >= 0 && j < h) {
            v4f = hdr.at(i, j);
            v4b = linear_to_srgb(hdr.at(i, j));
        }
    } else if (!ldr.empty()) {
        auto w = ldr.width(), h = ldr.height();
        if (i >= 0 && i < w && j >= 0 && j < h) {
            v4f = srgb_to_linear(ldr.at(i, j));
            v4b = ldr.at(i, j);
        }
    }
    char buf[1024];
    sprintf(buf, "%5d %5d", i, j);
    draw_label_widget(win, lbl + "mouse pos", buf);
    sprintf(buf, "%4.4g %4.4g %4.4g %4.4g", v4f.x, v4f.y, v4f.z, v4f.w);
    draw_label_widget(win, lbl + "hdr val", buf);
    sprintf(
        buf, "%3d %3d %3d %3d", (int)v4b.x, (int)v4b.y, (int)v4b.z, (int)v4b.w);
    draw_label_widget(win, lbl + "ldr val", buf);
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR SIMPLE SCENE UI
// -----------------------------------------------------------------------------
namespace ygl {

// Implementation of draw camera
struct draw_camera_visitor {
    gl_window* win = nullptr;
    int edited = 0;

    // constructor
    draw_camera_visitor(gl_window* win_) : win(win_) {}

    template <typename T>
    void operator()(T& val, const visit_var& var) {
        edited += draw_value_widget(win, var.name, val, var.min, var.max,
            var.type == visit_var_type::color);
    }

    template <typename T>
    void operator()(
        T* val, const visit_var& var = visit_var{"", visit_var_type::object}) {
        if (!val) return;
        // TODO: fix this separator
        draw_separator_widget(win);
        draw_groupid_widget_begin(win, val);
        draw_separator_widget(win);
        visit(*val, *this);
        draw_groupid_widget_end(win);
    }
};

// Implementation of camera selection
bool draw_camera_widgets(gl_window* win, const std::string& lbl, camera* cam) {
    if (!cam) return false;
    if (draw_header_widget(win, lbl)) {
        draw_groupid_widget_begin(win, cam);
        auto visitor = draw_camera_visitor{win};
        visitor(cam);
        draw_groupid_widget_end(win);
        return visitor.edited;
    } else {
        return false;
    }
}

// get typed selection
template <typename T>
inline T*& get_typed_selection(scene_selection& sel) {
    return nullptr;
}
template <>
inline camera*& get_typed_selection<camera>(scene_selection& sel) {
    return sel.cam;
}
template <>
inline shape_group*& get_typed_selection<shape_group>(scene_selection& sel) {
    return sel.sgr;
}
template <>
inline shape*& get_typed_selection<shape>(scene_selection& sel) {
    return sel.shp;
}
template <>
inline material*& get_typed_selection<material>(scene_selection& sel) {
    return sel.mat;
}
template <>
inline texture*& get_typed_selection<texture>(scene_selection& sel) {
    return sel.txt;
}
template <>
inline instance*& get_typed_selection<instance>(scene_selection& sel) {
    return sel.ist;
}
template <>
inline environment*& get_typed_selection<environment>(scene_selection& sel) {
    return sel.env;
}
template <>
inline node*& get_typed_selection<node>(scene_selection& sel) {
    return sel.nde;
}
template <>
inline animation_group*& get_typed_selection<animation_group>(
    scene_selection& sel) {
    return sel.agr;
}
template <>
inline animation*& get_typed_selection<animation>(scene_selection& sel) {
    return sel.anm;
}

// Implementation of draw tree
struct draw_tree_visitor {
    gl_window* win = nullptr;
    scene_selection& sel;

    // constructor
    draw_tree_visitor(gl_window* win_, scene_selection& sel_)
        : win(win_), sel(sel_) {}

    // generic callback
    template <typename T>
    void operator()(T& val, const visit_var&) {}
    // callback for texture_info
    void operator()(texture_info* val, const visit_var&) {}
    // callback for array
    template <typename T>
    void operator()(std::vector<T*>& val, const visit_var& var) {
        if (draw_tree_widget_begin(win, var.name)) {
            for (auto v : val) (*this)(v, var);
            draw_tree_widget_end(win);
        }
    }
    // callback for pointer
    template <typename T>
    void operator()(T* val, const visit_var& var) {
        if (!val) return;
        auto lbl = val->name;
        if (!var.name.empty()) lbl = var.name + ": " + val->name;
        void* selection = get_typed_selection<T>(sel);
        auto open = draw_tree_widget_begin(win, lbl, selection, val);
        if (selection == val) {
            clear_selection(sel);
            get_typed_selection<T>(sel) = val;
        }
        if (open) {
            visit(val, *this);
            draw_tree_widget_end(win);
        }
    }
    // callback for pointer
    void operator()(scene* val, const visit_var& var) { visit(val, *this); }
};

// Implementation of draw elements
struct draw_elem_visitor {
    gl_window* win = nullptr;
    scene* scn = nullptr;
    scene_selection& sel;
    const std::unordered_map<texture*, gl_texture>* gl_txt;
    int edited = 0;

    // constructor
    draw_elem_visitor(gl_window* win_, scene* scn_, scene_selection& sel_,
        const std::unordered_map<texture*, gl_texture>* gl_txt_)
        : win(win_), scn(scn_), sel(sel_), gl_txt(gl_txt_) {}

    template <typename T>
    void operator()(T& val, const visit_var& var) {
        edited += draw_value_widget(win, var.name, val, var.min, var.max,
            var.type == visit_var_type::color);
    }

    void operator()(texture_info* val, const visit_var& var) {}
    void operator()(make_hair_params* val, const visit_var& var) {}

    template <typename T>
    void operator()(image<T>& val, const visit_var& var) {
        if (empty(val)) return;
        auto size = format("{} x {}", val.width(), val.height());
        draw_label_widget(win, var.name, size);
    }

    template <typename T>
    void operator()(std::vector<T*>& val, const visit_var& var) {
        if (var.type != visit_var_type::reference) {
            for (auto idx = 0; idx < val.size(); idx++) visit(val[idx], *this);
        }
    }
    template <typename T>
    void operator()(std::vector<T>& val, const visit_var& var) {
        if (val.empty()) return;
        draw_label_widget(win, var.name, (int)val.size());
    }
    template <typename T1, typename T2>
    void operator()(std::vector<std::pair<T1*, T2*>>& val, const visit_var&) {
        // TODO
    }

    template <typename T>
    void operator()(
        T* val, const visit_var& var = visit_var{"", visit_var_type::object}) {
        if (!val) return;
        if (var.type == visit_var_type::reference) {
            edited += draw_combo_widget(win, var.name, val, elems(val));
        } else {
            if (!val) return;
            // TODO: fix this separator
            draw_separator_widget(win);
            draw_groupid_widget_begin(win, val);
            draw_separator_widget(win);
            visit(*val, *this);
            preview(val);
            draw_groupid_widget_end(win);
        }
    }

    template <typename T, typename TT>
    void operator()(T* elem, std::vector<TT*>& telems) {
        if (!elem) return;
        auto telem = (TT*)nullptr;
        for (auto te : telems)
            if (te->name == elem->name) telem = te;
        if (telem) {
            auto last = edited;
            operator()(telem);
            if (last != edited) update_proc_elem(scn, elem, telem);
        }
        operator()(elem);
    }

    template <typename T>
    const std::vector<T*>& elems(T* aux) {
        static auto v = std::vector<T*>();
        return v;
    }
    const std::vector<texture*>& elems(texture* aux) { return scn->textures; }
    const std::vector<shape_group*>& elems(shape_group* aux) {
        return scn->shapes;
    }
    const std::vector<camera*> elems(camera* aux) { return scn->cameras; }
    const std::vector<material*> elems(material* aux) { return scn->materials; }
    const std::vector<instance*> elems(instance* aux) { return scn->instances; }
    const std::vector<environment*> elems(environment* aux) {
        return scn->environments;
    }
    const std::vector<node*> elems(node* aux) { return scn->nodes; }
    const std::vector<animation_group*> elems(animation_group* aux) {
        return scn->animations;
    }

    template <typename T>
    void preview(T* val) {}
    void preview(texture* txt) {
        if (!gl_txt) return;
        if (!contains(*gl_txt, txt)) return;
        draw_image_widget(win, (gl_texture&)gl_txt->at(txt), {128, 128});
    }
};

template <typename T, typename T1>
inline bool draw_add_elem_widgets(gl_window* win, scene* scn,
    const std::string& lbl, std::vector<T*>& elems,
    std::vector<T1*>& proc_elems, scene_selection& sel,
    std::vector<ygl::scene_selection>& update_list) {
    static auto count = 0;
    if (draw_button_widget(win, "add " + lbl)) {
        auto name = lbl + "_" + std::to_string(count++);
        elems.push_back(new T());
        elems.back()->name = name;
        proc_elems.push_back(new T1());
        proc_elems.back()->name = name;
        clear_selection(sel);
        get_typed_selection<T>(sel) = elems.back();
        update_list.push_back(sel);
        update_proc_elem(scn, elems.back(), proc_elems.back());
        return true;
    }
    draw_continue_widget(win);
    return false;
}

bool draw_scene_widgets(gl_window* win, const std::string& lbl, scene* scn,
    scene_selection& sel, std::vector<ygl::scene_selection>& update_list,
    const std::unordered_map<texture*, gl_texture>& gl_txt,
    proc_scene* test_scn) {
    static auto test_scn_def = proc_scene();

    if (!scn) return false;
    if (!lbl.empty() && !draw_header_widget(win, lbl)) return false;
    draw_groupid_widget_begin(win, scn);
    // draw_scroll_widget_begin(win, "model", 240, false);
    auto tree_visitor = draw_tree_visitor{win, sel};
    tree_visitor(scn, visit_var{"", visit_var_type::object});
    // draw_scroll_widget_end(win);

    auto update_len = update_list.size();
    if (test_scn) {
        draw_add_elem_widgets(
            win, scn, "cam", scn->cameras, test_scn->cameras, sel, update_list);
        draw_add_elem_widgets(win, scn, "txt", scn->textures,
            test_scn->textures, sel, update_list);
        draw_add_elem_widgets(win, scn, "mat", scn->materials,
            test_scn->materials, sel, update_list);
        auto shp_add = draw_add_elem_widgets(
            win, scn, "shp", scn->shapes, test_scn->shapes, sel, update_list);
        if (shp_add) {
            scn->instances.push_back(new instance());
            scn->instances.back()->name = scn->shapes.back()->name;
            scn->instances.back()->shp = scn->shapes.back();
            test_scn->instances.push_back(new proc_instance());
            test_scn->instances.back()->name = scn->instances.back()->name;
            test_scn->instances.back()->shape =
                scn->instances.back()->shp->name;
            update_list.push_back({});
            update_list.back().ist = scn->instances.back();
        }
        draw_add_elem_widgets(win, scn, "ist", scn->instances,
            test_scn->instances, sel, update_list);
        draw_add_elem_widgets(win, scn, "env", scn->environments,
            test_scn->environments, sel, update_list);
        draw_add_elem_widgets(win, scn, "anim", scn->animations,
            test_scn->animations, sel, update_list);
    }

    auto test_scn_res = (test_scn) ? test_scn : &test_scn_def;
    auto elem_visitor = draw_elem_visitor{win, scn, sel, &gl_txt};
    elem_visitor(sel.cam, test_scn_res->cameras);
    elem_visitor(sel.shp);
    elem_visitor(sel.sgr, test_scn_res->shapes);
    elem_visitor(sel.sgr);
    elem_visitor(sel.ist, test_scn_res->instances);
    elem_visitor(sel.txt, test_scn_res->textures);
    elem_visitor(sel.mat, test_scn_res->materials);
    elem_visitor(sel.env, test_scn_res->environments);
    elem_visitor(sel.nde, test_scn_res->nodes);
    elem_visitor(sel.anm);
    elem_visitor(sel.agr, test_scn_res->animations);
    if (elem_visitor.edited) update_list.push_back(sel);

    draw_groupid_widget_end(win);
    return update_list.size() != update_len;
}

}  // namespace ygl

#endif
