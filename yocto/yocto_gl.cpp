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

//
// # Todo
//
// ## Next
//
// - remove elem cdf from shape
//    - put in lights for trace
//    - make lights a structure instead of a vector
//    - make light not a pointer
// - move shape CDF in lights
// - move pixel sampling to trace_block and not in sampler
// - remove stdprogram state
//
// ## General
//
// - remove python operators
//
// ## Infrastructure
//
// - remove default environment
// - meshes with multiple shapes
// - uniform serialization
//    - consider simpler serialization code based on input flag
//    - consider json archive model (use this to define to_son/from_from)
// - lights in scene for viewers
//    - do this only if necessary
// - check rotation and decompoaition of rotations
//    - see euclideanspace.com
//
// ## BVH
//
// - consider merging axis with internal
// - simplify build node
//
// ## Trace
//
// - handle missing environment
// - envmap sampling
// - sampler simplification
//     https://lemire.me/blog/2017/09/18/visiting-all-values-in-an-array-exactly-once-in-random-order/
//     https://lemire.me/blog/2017/09/26/benchmarking-algorithms-to-visit-all-values-in-an-array-in-random-order/
// - look at simpler denoiser
// - yitrace: check editing
// - yitrace: consider update
//
// ## Shade
//
// - start in edit mode
// - show edit scene
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
void camera_turntable(vec3f& from, vec3f& to, vec3f& up, const vec3f& rotate,
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
    if (rotate.x || rotate.y) {
        auto phi = atan2(frame.z.z, frame.z.x) + rotate.x;
        auto theta = acos(frame.z.y) + rotate.y;
        theta = clamp(theta, 0.001f, pif - 0.001f);
        auto new_z =
            vec3f{sin(theta) * cos(phi), cos(theta), sin(theta) * sin(phi)};
        auto new_center = frame.o - frame.z * focus;
        auto new_o = new_center + new_z * focus;
        frame = lookat_frame3(new_o, new_center, {0, 1, 0});
        focus = length(new_o - new_center);
    }

    // pan if necessary
    if (dolly) {
        auto c = frame.o - frame.z * focus;
        focus = max(focus + dolly, 0.001f);
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

    frame.rot() = rotation_mat3(vec3f{1, 0, 0}, rotate.y) * frame.rot() *
                  rotation_mat3(vec3f{0, 1, 0}, rotate.x);
    frame.pos() += transl.x * x + transl.y * y + transl.z * z;
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

// Compute per-vertex normals/tangents for lines, triangles and quads with
// positions pos. Weighted indicated whether the normals/tangents are
// weighted by line length.
vector<vec3f> compute_normals(const vector<vec2i>& lines,
    const vector<vec3i>& triangles, const vector<vec4i>& quads,
    const vector<vec3f>& pos, bool weighted) {
    auto norm = vector<vec3f>(pos.size(), zero3f);
    for (auto& l : lines) {
        auto n = pos[l.y] - pos[l.x];
        if (!weighted) n = normalize(n);
        for (auto vid : l) norm[vid] += n;
    }
    for (auto& t : triangles) {
        auto n = cross(pos[t.y] - pos[t.x], pos[t.z] - pos[t.x]);
        if (!weighted) n = normalize(n);
        for (auto vid : t) norm[vid] += n;
    }
    for (auto& q : quads) {
        auto n = cross(pos[q.y] - pos[q.x], pos[q.w] - pos[q.x]) +
                 cross(pos[q.w] - pos[q.z], pos[q.x] - pos[q.z]);
        if (!weighted) n = normalize(n);
        for (auto vid : q) norm[vid] += n;
    }
    for (auto& n : norm) n = normalize(n);
    return norm;
}

// Compute per-vertex tangent frame for triangle meshes.
// Tangent space is defined by a four component vector.
// The first three components are the tangent with respect to the U texcoord.
// The fourth component is the sign of the tangent wrt the V texcoord.
// Tangent frame is useful in normal mapping.
vector<vec4f> compute_tangent_frames(const vector<vec3i>& triangles,
    const vector<vec3f>& pos, const vector<vec3f>& norm,
    const vector<vec2f>& texcoord, bool weighted) {
    auto tangu = vector<vec3f>(pos.size(), zero3f);
    auto tangv = vector<vec3f>(pos.size(), zero3f);
    for (auto& t : triangles) {
        auto tutv = triangle_tangents_fromuv(pos[t.x], pos[t.y], pos[t.z],
            texcoord[t.x], texcoord[t.y], texcoord[t.z]);
        if (!weighted) tutv = {normalize(tutv.first), normalize(tutv.second)};
        for (auto vid : t) tangu[vid] += tutv.first;
        for (auto vid : t) tangv[vid] += tutv.second;
    }
    for (auto& t : tangu) t = normalize(t);
    for (auto& t : tangv) t = normalize(t);
    auto tangsp = vector<vec4f>(pos.size(), zero4f);
    for (auto i = 0; i < pos.size(); i++) {
        tangu[i] = orthonormalize(tangu[i], norm[i]);
        auto s = (dot(cross(norm[i], tangu[i]), tangv[i]) < 0) ? -1.0f : 1.0f;
        tangsp[i] = {tangu[i].x, tangu[i].y, tangu[i].z, s};
    }
    return tangsp;
}

// Apply skinning
void compute_skinning(const vector<vec3f>& pos, const vector<vec3f>& norm,
    const vector<vec4f>& weights, const vector<vec4i>& joints,
    const vector<mat4f>& xforms, vector<vec3f>& skinned_pos,
    vector<vec3f>& skinned_norm) {
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
void compute_skinning(const vector<vec3f>& pos, const vector<vec3f>& norm,
    const vector<vec4f>& weights, const vector<vec4i>& joints,
    const vector<frame3f>& xforms, vector<vec3f>& skinned_pos,
    vector<vec3f>& skinned_norm) {
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
void compute_matrix_skinning(const vector<vec3f>& pos,
    const vector<vec3f>& norm, const vector<vec4f>& weights,
    const vector<vec4i>& joints, const vector<mat4f>& xforms,
    vector<vec3f>& skinned_pos, vector<vec3f>& skinned_norm) {
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
vector<vec2i> get_edges(const vector<vec2i>& lines,
    const vector<vec3i>& triangles, const vector<vec4i>& quads) {
    auto edges = vector<vec2i>();
    auto eset = unordered_set<vec2i>();
    for (auto e : lines) {
        e = {min(e.x, e.y), max(e.x, e.y)};
        if (!eset.insert(e).second) continue;
        eset.insert({e.y, e.x});
        edges += e;
    }
    for (auto& t : triangles) {
        for (auto e : {vec2i{t.x, t.y}, vec2i{t.y, t.z}, vec2i{t.z, t.x}}) {
            e = {min(e.x, e.y), max(e.x, e.y)};
            if (!eset.insert(e).second) continue;
            eset.insert({e.y, e.x});
            edges += e;
        }
    }
    for (auto& q : quads) {
        for (auto e : {vec2i{q.x, q.y}, vec2i{q.y, q.z}, vec2i{q.z, q.w},
                 vec2i{q.w, q.x}}) {
            if (e.x == e.y) continue;
            e = {min(e.x, e.y), max(e.x, e.y)};
            if (!eset.insert(e).second) continue;
            eset.insert({e.y, e.x});
            edges += e;
        }
    }

    return edges;
}

// Create an array of boundary edges. Lines are always considered boundaries.
vector<vec2i> get_boundary_edges(const vector<vec2i>& lines,
    const vector<vec3i>& triangles, const vector<vec4i>& quads) {
    auto ecount = unordered_map<vec2i, int>();

    // lines are added manually later
    for (auto l : lines) { ecount.insert({l, 2}); }
    for (auto t : triangles) {
        for (auto e : {vec2i{t.x, t.y}, vec2i{t.y, t.z}, vec2i{t.z, t.x}}) {
            e = {min(e.x, e.y), max(e.x, e.y)};
            auto ins = ecount.insert({e, 1});
            if (!ins.second) ins.first->second += 1;
        }
    }
    for (auto q : quads) {
        for (auto e : {vec2i{q.x, q.y}, vec2i{q.y, q.z}, vec2i{q.z, q.w},
                 vec2i{q.w, q.x}}) {
            if (e.x == e.y) continue;
            e = {min(e.x, e.y), max(e.x, e.y)};
            auto ins = ecount.insert({e, 1});
            if (!ins.second) ins.first->second += 1;
        }
    }

    auto boundary = lines;
    for (auto ec : ecount) {
        if (ec.second > 1) continue;
        boundary += ec.first;
    }

    return boundary;
}

// Get a list of all unique vertices.
vector<int> get_verts(const vector<vec2i>& lines,
    const vector<vec3i>& triangles, const vector<vec4i>& quads) {
    auto verts = vector<int>();
    auto vset = unordered_set<int>();
    for (auto l : lines)
        for (auto vid : l)
            if (vset.insert(vid).second) verts += vid;
    for (auto t : triangles)
        for (auto vid : t)
            if (vset.insert(vid).second) verts += vid;
    for (auto q : quads)
        for (auto vid : q)
            if (vset.insert(vid).second) verts += vid;
    return verts;
}

// Create an array of boundary vertices. Lines are always considered
// boundaries.
vector<int> get_boundary_verts(const vector<vec2i>& lines,
    const vector<vec3i>& triangles, const vector<vec4i>& quads) {
    return get_verts(get_boundary_edges(lines, triangles, quads), {}, {});
}

// Convert quads to triangles
vector<vec3i> convert_quads_to_triangles(const vector<vec4i>& quads) {
    auto triangles = vector<vec3i>();
    triangles.reserve(quads.size() * 2);
    for (auto& q : quads) {
        triangles += {q.x, q.y, q.w};
        if (q.z != q.w) triangles += {q.z, q.w, q.y};
    }
    return triangles;
}

// Convert quads to triangles with a diamond-like topology.
// Quads have to be consecutive one row after another.
vector<vec3i> convert_quads_to_triangles(
    const vector<vec4i>& quads, int row_length) {
    auto triangles = vector<vec3i>();
    triangles.reserve(quads.size() * 2);
    for (auto& q : quads) {
        triangles += {q.x, q.y, q.w};
        if (q.z != q.w) triangles += {q.z, q.w, q.y};
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
vector<vec2i> convert_bezier_to_lines(const vector<vec4i>& beziers) {
    auto lines = vector<vec2i>();
    lines.reserve(beziers.size() * 3);
    for (auto& b : beziers) {
        lines += {b.x, b.y};
        lines += {b.y, b.z};
        lines += {b.z, b.w};
    }
    return lines;
}

// Convert face varying data to single primitives. Returns the quads indices
// and filled vectors for pos, norm and texcoord.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>
convert_face_varying(const vector<vec4i>& quads_pos,
    const vector<vec4i>& quads_norm, const vector<vec4i>& quads_texcoord,
    const vector<vec3f>& pos, const vector<vec3f>& norm,
    const vector<vec2f>& texcoord) {
    // make faces unique
    unordered_map<vec3i, int> vert_map;
    auto quads = vector<vec4i>(quads_pos.size());
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
    auto qpos = vector<vec3f>();
    if (!pos.empty()) {
        qpos.resize(vert_map.size());
        for (auto& kv : vert_map) { qpos[kv.second] = pos[kv.first.x]; }
    }
    auto qnorm = vector<vec3f>();
    if (!norm.empty()) {
        qnorm.resize(vert_map.size());
        for (auto& kv : vert_map) { qnorm[kv.second] = norm[kv.first.y]; }
    }
    auto qtexcoord = vector<vec2f>();
    if (!texcoord.empty()) {
        qtexcoord.resize(vert_map.size());
        for (auto& kv : vert_map) {
            qtexcoord[kv.second] = texcoord[kv.first.z];
        }
    }

    // done
    return {quads, qpos, qnorm, qtexcoord};
}

// wrapper for implementation below
inline float _subdivide_normalize(float x) { return x; }
inline vec2f _subdivide_normalize(const vec2f& x) { return normalize(x); }
inline vec3f _subdivide_normalize(const vec3f& x) { return normalize(x); }
inline vec4f _subdivide_normalize(const vec4f& x) { return normalize(x); }

// Tesselate lines, triangles and quads by spolitting edges.
// Returns the tesselated elements and dictionaries for vertex calculations.
tuple<vector<vec2i>, vector<vec3i>, vector<vec4i>, vector<vec2i>, vector<vec4i>>
subdivide_elems_linear(const vector<vec2i>& lines,
    const vector<vec3i>& triangles, const vector<vec4i>& quads, int nverts) {
    if (!nverts) return {};
    auto emap = unordered_map<vec2i, int>();
    auto edges = vector<vec2i>();
    for (auto e : lines) {
        if (contains(emap, e)) continue;
        emap[{e.x, e.y}] = nverts + (int)edges.size();
        emap[{e.y, e.x}] = nverts + (int)edges.size();
        edges += e;
    }
    for (auto& t : triangles) {
        for (auto e : {vec2i{t.x, t.y}, vec2i{t.y, t.z}, vec2i{t.z, t.x}}) {
            if (contains(emap, e)) continue;
            emap[{e.x, e.y}] = nverts + (int)edges.size();
            emap[{e.y, e.x}] = nverts + (int)edges.size();
            edges += e;
        }
    }
    for (auto& q : quads) {
        for (auto e : {vec2i{q.x, q.y}, vec2i{q.y, q.z}, vec2i{q.z, q.w},
                 vec2i{q.w, q.x}}) {
            if (e.x == e.y) continue;
            if (contains(emap, e)) continue;
            emap[{e.x, e.y}] = nverts + (int)edges.size();
            emap[{e.y, e.x}] = nverts + (int)edges.size();
            edges += e;
        }
    }

    auto tlines = vector<vec2i>();
    tlines.reserve(lines.size() * 2);
    for (auto& l : lines) {
        tlines += {l.x, emap.at(l)};
        tlines += {emap.at(l), l.y};
    }

    auto ttriangles = vector<vec3i>();
    ttriangles.reserve(triangles.size() * 4);
    for (auto& t : triangles) {
        ttriangles.push_back({t.x, emap.at({t.x, t.y}), emap.at({t.z, t.x})});
        ttriangles.push_back({t.y, emap.at({t.y, t.z}), emap.at({t.x, t.y})});
        ttriangles.push_back({t.z, emap.at({t.z, t.x}), emap.at({t.y, t.z})});
        ttriangles.push_back(
            {emap.at({t.x, t.y}), emap.at({t.y, t.z}), emap.at({t.z, t.x})});
    }

    auto tquads = vector<vec4i>();
    tquads.reserve(quads.size() * 4);
    for (auto fkv : enumerate(quads)) {
        auto f = fkv.second;
        auto fvert = nverts + (int)edges.size() + fkv.first;
        if (f.z != f.w) {
            tquads += {f.x, emap.at({f.x, f.y}), fvert, emap.at({f.w, f.x})};
            tquads += {f.y, emap.at({f.y, f.z}), fvert, emap.at({f.x, f.y})};
            tquads += {f.z, emap.at({f.z, f.w}), fvert, emap.at({f.y, f.z})};
            tquads += {f.w, emap.at({f.w, f.x}), fvert, emap.at({f.z, f.w})};
        } else {
            tquads += {f.x, emap.at({f.x, f.y}), fvert, emap.at({f.z, f.x})};
            tquads += {f.y, emap.at({f.y, f.z}), fvert, emap.at({f.x, f.y})};
            tquads += {f.z, emap.at({f.z, f.x}), fvert, emap.at({f.y, f.z})};
        }
    }
    tquads.shrink_to_fit();

    return {tlines, ttriangles, tquads, edges, quads};
}

// Subdivide vertex properties given the maps
template <typename T>
vector<T> subdivide_vert_linear(const vector<T>& vert,
    const vector<vec2i>& edges, const vector<vec4i>& faces, bool normalized) {
    if (vert.empty()) return {};

    auto tvert = vector<T>();
    tvert.reserve(vert.size() + edges.size() + faces.size());

    tvert += vert;
    for (auto e : edges) tvert += (vert[e.x] + vert[e.y]) / 2;
    for (auto f : faces) {
        if (f.z != f.w)
            tvert += (vert[f.x] + vert[f.y] + vert[f.z] + vert[f.w]) / 4;
        else
            tvert += (vert[f.x] + vert[f.y] + vert[f.z]) / 3;
    }

    if (normalized) {
        for (auto& n : tvert) n = _subdivide_normalize(n);
    }

    return tvert;
}

// template instantiations
template vector<float> subdivide_vert_linear(const vector<float>& vert,
    const vector<vec2i>& edges, const vector<vec4i>& faces, bool normalized);
template vector<vec2f> subdivide_vert_linear(const vector<vec2f>& vert,
    const vector<vec2i>& edges, const vector<vec4i>& faces, bool normalized);
template vector<vec3f> subdivide_vert_linear(const vector<vec3f>& vert,
    const vector<vec2i>& edges, const vector<vec4i>& faces, bool normalized);
template vector<vec4f> subdivide_vert_linear(const vector<vec4f>& vert,
    const vector<vec2i>& edges, const vector<vec4i>& faces, bool normalized);

// Performs the smoothing step of Catmull-Clark. Start with a tesselate quad
// mesh obtained with subdivide_elems_linear() and subdivide_vert_linear(). To
// handle open meshes with boundary, get the boundary from make_boundary_edge()
// and pass it as crease_lines. To fix the boundary entirely, just get the
// boundary vertices and pass it as creases.
template <typename T>
vector<T> subdivide_vert_catmullclark(const vector<vec4i>& quads,
    const vector<T>& vert, const vector<vec2i>& crease_tlines,
    const vector<int>& crease_tpoints, bool normalized) {
    if (quads.empty() || vert.empty()) return vert;

    // define vertex valence ---------------------------
    auto val = vector<int>(vert.size(), 2);
    for (auto e : crease_tlines)
        for (auto vid : e) val[vid] = 1;
    for (auto vid : crease_tpoints) val[vid] = 0;

    // averaging pass ----------------------------------
    auto tvert = vector<T>(vert.size(), T());
    auto count = vector<int>(vert.size(), 0);
    for (auto p : crease_tpoints) {
        auto c = vert[p];
        if (val[p] == 0) tvert[p] += c;
        if (val[p] == 0) count[p] += 1;
    }
    for (auto e : crease_tlines) {
        auto c = (vert[e.x] + vert[e.y]) / 2.0f;
        for (auto vid : e) {
            if (val[vid] == 1) tvert[vid] += c;
            if (val[vid] == 1) count[vid] += 1;
        }
    }
    for (auto& f : quads) {
        auto c = (vert[f.x] + vert[f.y] + vert[f.z] + vert[f.w]) / 4.0f;
        for (auto vid : f) {
            if (val[vid] == 2) tvert[vid] += c;
            if (val[vid] == 2) count[vid] += 1;
        }
    }
    for (auto i = 0; i < vert.size(); i++) { tvert[i] /= (float)count[i]; }

    // correction pass ----------------------------------
    // p = p + (avg_p - p) * (4/avg_count)
    for (auto i = 0; i < vert.size(); i++) {
        if (val[i] != 2) continue;
        tvert[i] = vert[i] + (tvert[i] - vert[i]) * (4.0f / count[i]);
    }

    if (normalized) {
        for (auto& v : tvert) v = _subdivide_normalize(v);
    }

    return tvert;
}

// explicit instantiations
template vector<float> subdivide_vert_catmullclark(const vector<vec4i>& quads,
    const vector<float>& vert, const vector<vec2i>& crease_tlines,
    const vector<int>& crease_tpoints, bool normalized);
template vector<vec2f> subdivide_vert_catmullclark(const vector<vec4i>& quads,
    const vector<vec2f>& vert, const vector<vec2i>& crease_tlines,
    const vector<int>& crease_tpoints, bool normalized);
template vector<vec3f> subdivide_vert_catmullclark(const vector<vec4i>& quads,
    const vector<vec3f>& vert, const vector<vec2i>& crease_tlines,
    const vector<int>& crease_tpoints, bool normalized);
template vector<vec4f> subdivide_vert_catmullclark(const vector<vec4i>& quads,
    const vector<vec4f>& vert, const vector<vec2i>& crease_tlines,
    const vector<int>& crease_tpoints, bool normalized);

// Subdivide bezier recursive by splitting each segment into two in the middle.
// Returns the tesselated elements and dictionaries for vertex calculations.
tuple<vector<vec4i>, vector<int>, vector<vec4i>> subdivide_bezier_recursive(
    const vector<vec4i>& beziers, int nverts) {
    if (!nverts) return {};
    auto vmap = unordered_map<int, int>();
    auto verts = vector<int>();
    for (auto& b : beziers) {
        if (!contains(vmap, b.x)) {
            vmap[b.x] = verts.size();
            verts += b.x;
        }
        if (!contains(vmap, b.w)) {
            vmap[b.w] = verts.size();
            verts += b.w;
        }
    }
    auto tbeziers = vector<vec4i>();
    tbeziers.reserve(beziers.size() * 2);
    for (auto b_kv : enumerate(beziers)) {
        auto b = b_kv.second;
        auto bo = (int)verts.size() + b_kv.first * 5;
        tbeziers += {vmap.at(b.x), bo + 0, bo + 1, bo + 2};
        tbeziers += {bo + 2, bo + 3, bo + 4, vmap.at(b.w)};
    }
    return {tbeziers, verts, beziers};
}

// Subdivide vertex properties given the maps
template <typename T>
vector<T> subdivide_vert_bezier(const vector<T>& vert, const vector<int>& verts,
    const vector<vec4i>& segments, bool normalized) {
    if (vert.empty()) return {};

    auto tvert = vector<T>();
    tvert.reserve(verts.size() + segments.size() * 5);

    for (auto v : verts) tvert += vert[v];
    for (auto s : segments) {
        tvert += vert[s.x] * (1.f / 2) + vert[s.y] * (1.f / 2);
        tvert += vert[s.x] * (1.f / 4) + vert[s.y] * (1.f / 2) +
                 vert[s.z] * (1.f / 4);
        tvert += vert[s.x] * (1.f / 8) + vert[s.y] * (3.f / 8) +
                 vert[s.z] * (3.f / 8) + vert[s.w] * (1.f / 8);
        tvert += vert[s.y] * (1.f / 4) + vert[s.z] * (1.f / 2) +
                 vert[s.w] * (1.f / 4);
        tvert += vert[s.z] * (1.f / 2) + vert[s.w] * (1.f / 2);
    }

    if (normalized) {
        for (auto& n : tvert) n = _subdivide_normalize(n);
    }

    return tvert;
}

// template instantiations
template vector<float> subdivide_vert_bezier(const vector<float>& vert,
    const vector<int>& verts, const vector<vec4i>& segments, bool normalized);
template vector<vec2f> subdivide_vert_bezier(const vector<vec2f>& vert,
    const vector<int>& verts, const vector<vec4i>& segments, bool normalized);
template vector<vec3f> subdivide_vert_bezier(const vector<vec3f>& vert,
    const vector<int>& verts, const vector<vec4i>& segments, bool normalized);
template vector<vec4f> subdivide_vert_bezier(const vector<vec4f>& vert,
    const vector<int>& verts, const vector<vec4i>& segments, bool normalized);

// Generate a rectangular grid of usteps x vsteps uv values for parametric
// surface generation.
tuple<vector<vec4i>, vector<vec2f>> make_uvquads(
    int usteps, int vsteps, bool uwrap, bool vwrap, bool vpole0, bool vpole1) {
    auto uvert = (uwrap) ? usteps : usteps + 1;
    auto vvert = (vwrap) ? vsteps : vsteps + 1;
    auto vid = [=](int i, int j) {
        if (uwrap) i = i % usteps;
        if (vwrap) j = j % vsteps;
        return j * uvert + i;
    };

    auto uv = vector<vec2f>(uvert * vvert);
    for (auto j = 0; j < vvert; j++) {
        for (auto i = 0; i < uvert; i++) {
            uv[vid(i, j)] = {i / (float)usteps, j / (float)vsteps};
        }
    }

    auto quads = vector<vec4i>(usteps * vsteps);
    for (auto j = 0; j < vsteps; j++) {
        for (auto i = 0; i < usteps; i++) {
            quads[j * usteps + i] = {
                vid(i, j), vid(i + 1, j), vid(i + 1, j + 1), vid(i, j + 1)};
        }
    }

    if (vpole0) {
        if (vwrap) throw runtime_error("cannot have a pole with wrapping");
        uv = vector<vec2f>(uv.begin() + uvert, uv.end());
        uv.insert(uv.begin(), {0, 0});
        for (auto& q : quads) {
            for (auto& vid : q) { vid = (vid < usteps) ? 0 : vid - uvert + 1; }
            if (q.x == 0 && q.y == 0) q = {q.z, q.w, q.x, q.y};
        }
    }

    if (vpole1) {
        if (vwrap) throw runtime_error("cannot have a pole with wrapping");
        auto pid = (int)uv.size() - uvert;
        uv = vector<vec2f>(uv.begin(), uv.end() - uvert);
        uv.insert(uv.end(), {0, 1});
        for (auto& q : quads) {
            for (auto& vid : q) { vid = (vid < pid) ? vid : pid; }
        }
    }

    return {quads, uv};
}

// Generate parametric num lines of usteps segments.
tuple<vector<vec2i>, vector<vec2f>> make_uvlines(int num, int usteps) {
    auto vid = [usteps](int i, int j) { return j * (usteps + 1) + i; };
    auto uv = vector<vec2f>((usteps + 1) * num);
    for (auto j = 0; j < num; j++) {
        for (auto i = 0; i <= usteps; i++) {
            uv[vid(i, j)] = {i / (float)usteps, j / (float)num};
        }
    }

    auto lines = vector<vec2i>(usteps * num);
    for (int j = 0; j < num; j++) {
        for (int i = 0; i < usteps; i++) {
            lines[j * usteps + i] = {vid(i, j), vid(i + 1, j)};
        }
    }

    return {lines, uv};
}

// Generate a parametric point set. Mostly here for completeness.
tuple<vector<int>, vector<vec2f>> make_uvpoints(int num) {
    auto uv = vector<vec2f>(num);
    for (auto i = 0; i < num; i++) { uv[i] = {i / (float)num, 0}; }

    auto points = vector<int>(num);
    for (auto i = 0; i < num; i++) points[i] = i;

    return {points, uv};
}

// Merge elements between shapes. The elements are merged by increasing the
// array size of the second array by the number of vertices of the first.
// Vertex data can then be concatenated successfully.
tuple<vector<vec2i>, vector<vec3i>, vector<vec4i>> merge_elems(int nverts,
    const vector<vec2i>& lines1, const vector<vec3i>& triangles1,
    const vector<vec4i>& quads1, const vector<vec2i>& lines2,
    const vector<vec3i>& triangles2, const vector<vec4i>& quads2) {
    auto lines = lines1 + lines2;
    auto triangles = triangles1 + triangles2;
    auto quads = quads1 + quads2;
    for (auto i = lines1.size(); i < lines.size(); i++)
        lines[i] += {nverts, nverts};
    for (auto i = triangles1.size(); i < triangles.size(); i++)
        triangles[i] += {nverts, nverts, nverts};
    for (auto i = quads1.size(); i < quads.size(); i++)
        quads[i] += {nverts, nverts, nverts, nverts};
    return {lines, triangles, quads};
}

// Unshare shape data by duplicating all vertex data for each element,
// giving a faceted look. Note that faceted tangents are not computed.
tuple<vector<vec2i>, vector<vec3i>, vector<vec4i>, vector<int>> facet_elems(
    const vector<vec2i>& lines, const vector<vec3i>& triangles,
    const vector<vec4i>& quads) {
    auto verts = vector<int>();
    auto nlines = vector<vec2i>();
    for (auto l : lines) {
        nlines.push_back({(int)verts.size(), (int)verts.size() + 1});
        for (auto v : l) verts += v;
    }

    auto ntriangles = vector<vec3i>();
    for (auto t : triangles) {
        ntriangles.push_back(
            {(int)verts.size(), (int)verts.size() + 1, (int)verts.size() + 2});
        for (auto v : t) verts += v;
    }

    auto nquads = vector<vec4i>();
    for (auto q : quads) {
        if (q.z != q.w) {
            nquads.push_back({(int)verts.size(), (int)verts.size() + 1,
                (int)verts.size() + 2, (int)verts.size() + 3});
            for (auto v : q) verts += v;
        } else {
            nquads.push_back({(int)verts.size(), (int)verts.size() + 1,
                (int)verts.size() + 2, (int)verts.size() + 2});
            for (auto v : q.xyz()) verts += v;
        }
    }

    return {nlines, ntriangles, nquads, verts};
}

// Unshare vertices for faceting
template <typename T>
vector<T> facet_vert(const vector<T>& vert, const vector<int>& vmap) {
    if (vert.empty()) return vert;
    auto tvert = vector<T>(vmap.size());
    for (auto vkv : enumerate(vmap)) tvert[vkv.first] = vert[vkv.second];
    return tvert;
}

// instantations
template vector<float> facet_vert<float>(
    const vector<float>& vert, const vector<int>& vmap);
template vector<vec2f> facet_vert<vec2f>(
    const vector<vec2f>& vert, const vector<int>& vmap);
template vector<vec3f> facet_vert<vec3f>(
    const vector<vec3f>& vert, const vector<int>& vmap);
template vector<vec4f> facet_vert<vec4f>(
    const vector<vec4f>& vert, const vector<int>& vmap);

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION OF SHAPE SAMPLING
// -----------------------------------------------------------------------------
namespace ygl {

// Pick a point
int sample_points(int npoints, float re) {
    return clamp(0, npoints - 1, (int)(re * npoints));
}

// Compute a distribution for sampling points uniformly
vector<float> sample_points_cdf(int npoints) {
    auto cdf = vector<float>(npoints);
    for (auto i = 0; i < npoints; i++) cdf[i] = i + 1;
    return cdf;
}

// Pick a point
int sample_points(const vector<float>& cdf, float re) {
    re = clamp(re * cdf.back(), 0.0f, cdf.back() - 0.00001f);
    return (int)(std::upper_bound(cdf.begin(), cdf.end(), re) - cdf.begin());
}

// Compute a distribution for sampling lines uniformly
vector<float> sample_lines_cdf(
    const vector<vec2i>& lines, const vector<vec3f>& pos) {
    auto cdf = vector<float>(lines.size());
    for (auto i = 0; i < lines.size(); i++)
        cdf[i] = length(pos[lines[i].x] - pos[lines[i].y]);
    for (auto i = 1; i < lines.size(); i++) cdf[i] += cdf[i - 1];
    return cdf;
}

// Pick a point on lines
pair<int, vec2f> sample_lines(const vector<float>& cdf, float re, float ruv) {
    re = clamp(re * cdf.back(), 0.0f, cdf.back() - 0.00001f);
    auto eid =
        (int)(std::upper_bound(cdf.begin(), cdf.end(), re) - cdf.begin());
    return {eid, {1 - ruv, ruv}};
}

// Compute a distribution for sampling triangle meshes uniformly
vector<float> sample_triangles_cdf(
    const vector<vec3i>& triangles, const vector<vec3f>& pos) {
    auto cdf = vector<float>(triangles.size());
    for (auto i = 0; i < triangles.size(); i++)
        cdf[i] = triangle_area(
            pos[triangles[i].x], pos[triangles[i].y], pos[triangles[i].z]);
    for (auto i = 1; i < triangles.size(); i++) cdf[i] += cdf[i - 1];
    return cdf;
}

// Pick a point on a triangle mesh
pair<int, vec3f> sample_triangles(
    const vector<float>& cdf, float re, const vec2f& ruv) {
    re = clamp(re * cdf.back(), 0.0f, cdf.back() - 0.00001f);
    auto eid =
        (int)(std::upper_bound(cdf.begin(), cdf.end(), re) - cdf.begin());
    return {
        eid, {sqrt(ruv.x) * (1 - ruv.y), 1 - sqrt(ruv.x), ruv.y * sqrt(ruv.x)}};
}

// Compute a distribution for sampling quad meshes uniformly
vector<float> sample_quads_cdf(
    const vector<vec4i>& quads, const vector<vec3f>& pos) {
    auto cdf = vector<float>(quads.size());
    for (auto i = 0; i < quads.size(); i++)
        cdf[i] = quad_area(
            pos[quads[i].x], pos[quads[i].y], pos[quads[i].z], pos[quads[i].w]);
    for (auto i = 1; i < quads.size(); i++) cdf[i] += cdf[i - 1];
    return cdf;
}

// Pick a point on a quad mesh
pair<int, vec4f> sample_quads(
    const vector<float>& cdf, float re, const vec2f& ruv) {
    if (ruv.x < 0.5f) {
        auto eid = 0;
        auto euv = zero3f;
        std::tie(eid, euv) = sample_triangles(cdf, re, {ruv.x * 2, ruv.y});
        return {eid, {euv.x, euv.y, 0, euv.z}};
    } else {
        auto eid = 0;
        auto euv = zero3f;
        std::tie(eid, euv) =
            sample_triangles(cdf, re, {(ruv.x - 0.5f) * 2, ruv.y});
        return {eid, {0, euv.z, euv.x, euv.y}};
    }
}

// Samples a set of points over a triangle mesh uniformly. The rng function
// takes the point index and returns vec3f numbers uniform directibuted in
// [0,1]^3. unorm and texcoord are optional.
tuple<vector<vec3f>, vector<vec3f>, vector<vec2f>> sample_triangles_points(
    const vector<vec3i>& triangles, const vector<vec3f>& pos,
    const vector<vec3f>& norm, const vector<vec2f>& texcoord, int npoints,
    uint64_t seed) {
    auto sampled_pos = vector<vec3f>(npoints);
    auto sampled_norm = vector<vec3f>(norm.empty() ? 0 : npoints);
    auto sampled_texcoord = vector<vec2f>(texcoord.empty() ? 0 : npoints);
    auto cdf = sample_triangles_cdf(triangles, pos);
    auto rng = init_rng(seed);
    for (auto i = 0; i < npoints; i++) {
        auto eid = 0;
        auto euv = zero3f;
        std::tie(eid, euv) = sample_triangles(
            cdf, next_rand1f(rng), {next_rand1f(rng), next_rand1f(rng)});
        auto t = triangles[eid];
        sampled_pos[i] = pos[t.x] * euv.x + pos[t.y] * euv.y + pos[t.z] * euv.z;
        if (!sampled_norm.empty())
            sampled_norm[i] = normalize(
                norm[t.x] * euv.x + norm[t.y] * euv.y + norm[t.z] * euv.z);
        if (!sampled_texcoord.empty())
            sampled_texcoord[i] = texcoord[t.x] * euv.x +
                                  texcoord[t.y] * euv.y + texcoord[t.z] * euv.z;
    }

    return {sampled_pos, sampled_norm, sampled_texcoord};
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR IMAGEIO
// -----------------------------------------------------------------------------
namespace ygl {

// check hdr extensions
bool is_hdr_filename(const string& filename) {
    auto ext = path_extension(filename);
    return ext == ".hdr" || ext == ".exr";
}

// Loads an ldr image.
image4b load_image4b(const string& filename) {
    auto w = 0, h = 0, c = 0;
    auto pixels = unique_ptr<byte>(stbi_load(filename.c_str(), &w, &h, &c, 4));
    if (!pixels) return {};
    return make_image(w, h, (vec4b*)pixels.get());
}

// Loads an hdr image.
image4f load_image4f(const string& filename) {
    auto ext = path_extension(filename);
    auto w = 0, h = 0, c = 0;
    auto pixels = unique_ptr<float>(nullptr);
    if (ext == ".exr") {
        auto pixels_ = (float*)nullptr;
        if (!LoadEXR(&pixels_, &w, &h, filename.c_str(), nullptr))
            pixels = unique_ptr<float>(pixels_);
    } else {
        pixels = unique_ptr<float>(stbi_loadf(filename.c_str(), &w, &h, &c, 4));
    }
    if (!pixels) return {};
    return make_image(w, h, (vec4f*)pixels.get());
}

// Saves an ldr image.
bool save_image4b(const string& filename, const image4b& img) {
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
bool save_image4f(const string& filename, const image4f& img) {
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
vector<float> load_imagef(
    const string& filename, int& width, int& height, int& ncomp) {
    auto pixels = stbi_loadf(filename.c_str(), &width, &height, &ncomp, 0);
    if (!pixels) return {};
    auto ret = vector<float>(pixels, pixels + width * height * ncomp);
    free(pixels);
    return ret;
}

// Loads an image
vector<byte> load_image(
    const string& filename, int& width, int& height, int& ncomp) {
    auto pixels = stbi_load(filename.c_str(), &width, &height, &ncomp, 0);
    if (!pixels) return {};
    auto ret = vector<byte>(pixels, pixels + width * height * ncomp);
    free(pixels);
    return ret;
}

// Loads an image from memory.
vector<float> load_imagef_from_memory(const string& filename, const byte* data,
    int length, int& width, int& height, int& ncomp) {
    auto pixels =
        stbi_loadf_from_memory(data, length, &width, &height, &ncomp, 0);
    if (!pixels) return {};
    auto ret = vector<float>(pixels, pixels + width * height * ncomp);
    free(pixels);
    return ret;
}

// Loads an image from memory.
vector<byte> load_image_from_memory(const string& filename, const byte* data,
    int length, int& width, int& height, int& ncomp) {
    auto pixels =
        stbi_load_from_memory(data, length, &width, &height, &ncomp, 0);
    if (!pixels) return {};
    auto ret = vector<byte>(pixels, pixels + width * height * ncomp);
    free(pixels);
    return ret;
}

// Saves an image
bool save_imagef(const string& filename, int width, int height, int ncomp,
    const float* hdr) {
    if (path_extension(filename) == ".hdr") {
        return stbi_write_hdr(filename.c_str(), width, height, ncomp, hdr);
    } else {
        return false;
    }
}

// Saves an image
bool save_image(
    const string& filename, int width, int height, int ncomp, const byte* ldr) {
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
bool save_image(const string& filename, const image4f& hdr, float exposure,
    float gamma, bool filmic) {
    if (is_hdr_filename(filename)) {
        return save_image4f(filename, hdr);
    } else {
        auto ldr = tonemap_image(hdr, exposure, gamma, filmic);
        return save_image4b(filename, ldr);
    }
}

// Resize image.
void resize_image(const image4f& img, image4f& res_img, resize_filter filter,
    resize_edge edge, bool premultiplied_alpha) {
    static const auto filter_map = map<resize_filter, stbir_filter>{
        {resize_filter::def, STBIR_FILTER_DEFAULT},
        {resize_filter::box, STBIR_FILTER_BOX},
        {resize_filter::triangle, STBIR_FILTER_TRIANGLE},
        {resize_filter::cubic_spline, STBIR_FILTER_CUBICBSPLINE},
        {resize_filter::catmull_rom, STBIR_FILTER_CATMULLROM},
        {resize_filter::mitchell, STBIR_FILTER_MITCHELL}};

    static const auto edge_map =
        map<resize_edge, stbir_edge>{{resize_edge::def, STBIR_EDGE_CLAMP},
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
    static const auto filter_map = map<resize_filter, stbir_filter>{
        {resize_filter::def, STBIR_FILTER_DEFAULT},
        {resize_filter::box, STBIR_FILTER_BOX},
        {resize_filter::triangle, STBIR_FILTER_TRIANGLE},
        {resize_filter::cubic_spline, STBIR_FILTER_CUBICBSPLINE},
        {resize_filter::catmull_rom, STBIR_FILTER_CATMULLROM},
        {resize_filter::mitchell, STBIR_FILTER_MITCHELL}};

    static const auto edge_map =
        map<resize_edge, stbir_edge>{{resize_edge::def, STBIR_EDGE_CLAMP},
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
image4b tonemap_image(
    const image4f& hdr, float exposure, float gamma, bool filmic) {
    auto ldr = image4b(hdr.width(), hdr.height());
    auto scale = pow(2.0f, exposure);
    for (auto j = 0; j < hdr.height(); j++) {
        for (auto i = 0; i < hdr.width(); i++) {
            auto h = hdr[{i, j}];
            h.xyz() *= scale;
            if (filmic) {
                h.xyz() = {tonemap_filmic(h.x), tonemap_filmic(h.y),
                    tonemap_filmic(h.z)};
            } else {
                h.xyz() = {pow(h.x, 1 / gamma), pow(h.y, 1 / gamma),
                    pow(h.z, 1 / gamma)};
            }
            ldr[{i, j}] = float_to_byte(h);
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
    for (auto i : range(3)) {
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
            img[{i, j}] = {sky(theta, gamma) + sun(theta, gamma), 1};
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
            img[{i, j}] = float_to_byte({g, g, g, 1});
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
            img[{i, j}] = float_to_byte({g, g, g, 1});
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
            img[{i, j}] = float_to_byte({g, g, g, 1});
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
            img[{i, j}] = float_to_byte({g, g, g, 1});
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
    auto r = r0 * (1 - s) + r1 * s;
    if (dot(p01, p01) > r * r) return false;

    // intersection occurred: set params and exit
    ray_t = t;
    euv = {1 - s, s};

    return true;
}

// Intersect a ray with a triangle
bool intersect_triangle(const ray3f& ray, const vec3f& v0, const vec3f& v1,
    const vec3f& v2, float& ray_t, vec3f& euv) {
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
    euv = {1 - u - v, u, v};

    return true;
}

// Intersect a ray with a quad.
bool intersect_quad(const ray3f& ray, const vec3f& v0, const vec3f& v1,
    const vec3f& v2, const vec3f& v3, float& ray_t, vec4f& euv) {
    auto hit = false;
    auto tray = ray;
    if (intersect_triangle(tray, v0, v1, v3, ray_t, (vec3f&)euv)) {
        euv = {euv.x, euv.y, 0, euv.z};
        tray.tmax = ray_t;
        hit = true;
    }
    if (intersect_triangle(tray, v2, v3, v1, ray_t, (vec3f&)euv)) {
        euv = {0, 1 - euv.y, euv.y + euv.z - 1, 1 - euv.z};
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
    auto tuv = zero3f;
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
vec2f closestuv_line(const vec3f& pos, const vec3f& v0, const vec3f& v1) {
    auto ab = v1 - v0;
    auto d = dot(ab, ab);
    // Project c onto ab, computing parameterized position d(t) = a + t*(b –
    // a)
    auto u = dot(pos - v0, ab) / d;
    u = clamp(u, (float)0, (float)1);
    return {1 - u, u};
}

// TODO: documentation
bool overlap_line(const vec3f& pos, float dist_max, const vec3f& v0,
    const vec3f& v1, float r0, float r1, float& dist, vec2f& euv) {
    auto uv = closestuv_line(pos, v0, v1);
    // Compute projected position from the clamped t d = a + t * ab;
    auto p = lerp(v0, v1, uv.y);
    auto r = lerp(r0, r1, uv.y);
    auto d2 = dot(pos - p, pos - p);
    // check distance
    if (d2 > (dist_max + r) * (dist_max + r)) return false;
    // done
    dist = sqrt(d2);
    euv = uv;
    return true;
}

// TODO: documentation
// this is a complicated test -> I probably prefer to use a sequence of test
// (triangle body, and 3 edges)
vec3f closestuv_triangle(
    const vec3f& pos, const vec3f& v0, const vec3f& v1, const vec3f& v2) {
    auto ab = v1 - v0;
    auto ac = v2 - v0;
    auto ap = pos - v0;

    auto d1 = dot(ab, ap);
    auto d2 = dot(ac, ap);

    // corner and edge cases
    if (d1 <= 0 && d2 <= 0) return {1, 0, 0};

    auto bp = pos - v1;
    auto d3 = dot(ab, bp);
    auto d4 = dot(ac, bp);
    if (d3 >= 0 && d4 <= d3) return {0, 1, 0};

    auto vc = d1 * d4 - d3 * d2;
    if ((vc <= 0) && (d1 >= 0) && (d3 <= 0))
        return {1 - d1 / (d1 - d3), d1 / (d1 - d3), 0};

    auto cp = pos - v2;
    auto d5 = dot(ab, cp);
    auto d6 = dot(ac, cp);
    if (d6 >= 0 && d5 <= d6) return {0, 0, 1};

    auto vb = d5 * d2 - d1 * d6;
    if ((vb <= 0) && (d2 >= 0) && (d6 <= 0))
        return {1 - d2 / (d2 - d6), 0, d2 / (d2 - d6)};

    auto va = d3 * d6 - d5 * d4;
    if ((va <= 0) && (d4 - d3 >= 0) && (d5 - d6 >= 0)) {
        auto w = (d4 - d3) / ((d4 - d3) + (d5 - d6));
        return {0, 1 - w, w};
    }

    // face case
    auto denom = 1 / (va + vb + vc);
    auto v = vb * denom;
    auto w = vc * denom;
    return {1 - v - w, v, w};
}

// TODO: documentation
bool overlap_triangle(const vec3f& pos, float dist_max, const vec3f& v0,
    const vec3f& v1, const vec3f& v2, float r0, float r1, float r2, float& dist,
    vec3f& euv) {
    auto uv = closestuv_triangle(pos, v0, v1, v2);
    auto p = v0 * uv.x + v1 * uv.y + v2 * uv.z;
    auto r = r0 * uv.x + r1 * uv.y + r2 * uv.z;
    auto dd = dot(p - pos, p - pos);
    if (dd > (dist_max + r) * (dist_max + r)) return false;
    dist = sqrt(dd);
    euv = uv;
    return true;
}

// TODO: documentation
bool overlap_quad(const vec3f& pos, float dist_max, const vec3f& v0,
    const vec3f& v1, const vec3f& v2, const vec3f& v3, float r0, float r1,
    float r2, float r3, float& dist, vec4f& euv) {
    auto hit = false;
    if (overlap_triangle(
            pos, dist_max, v0, v1, v3, r0, r1, r3, dist, (vec3f&)euv)) {
        euv = {euv.x, euv.y, 0, euv.z};
        dist_max = dist;
        hit = true;
    }
    if (overlap_triangle(
            pos, dist_max, v2, v3, v1, r2, r3, r1, dist, (vec3f&)euv)) {
        // dist_max = dist;
        euv = {0, 1 - euv.y, euv.y + euv.z - 1, 1 - euv.z};
        hit = true;
    }
    return hit;
}

// TODO: documentation
bool overlap_tetrahedron(const vec3f& pos, const vec3f& v0, const vec3f& v1,
    const vec3f& v2, const vec3f& v3, vec4f& euv) {
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
    // check interior
    if (overlap_tetrahedron(pos, v0, v1, v2, v3, euv)) {
        dist = 0;
        return true;
    }

    // check faces
    auto hit = false;
    auto tuv = zero3f;
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

// Struct that pack a bounding box, its associate primitive index, and other
// data for faster hierarchy build.
// This is internal only and should not be used externally.
struct bvh_bound_prim {
    bbox3f bbox;   // bounding box
    vec3f center;  // bounding box center (for faster sort)
    int pid;       // primitive id
};

// Comparison function for each axis
struct bvh_bound_prim_comp {
    int axis;
    float middle;

    bvh_bound_prim_comp(int a, float m = 0) : axis(a), middle(m) {}

    bool operator()(const bvh_bound_prim& a, const bvh_bound_prim& b) const {
        return a.center[axis] < b.center[axis];
    }

    bool operator()(const bvh_bound_prim& a) const {
        return a.center[axis] < middle;
    }
};

// number of primitives to avoid splitting on
const int bvh_minprims = 4;

// Initializes the BVH node node that contains the primitives sorted_prims
// from start to end, by either splitting it into two other nodes,
// or initializing it as a leaf. When splitting, the heuristic heuristic is
// used and nodes added sequentially in the preallocated nodes array and
// the number of nodes nnodes is updated.
void make_bvh_node(bvh_tree* bvh, int nodeid, bvh_bound_prim* sorted_prims,
    int start, int end, bool equalsize) {
    // get node
    auto node = &bvh->nodes.at(nodeid);
    // compute node bounds
    node->bbox = invalid_bbox3f;
    for (auto i = start; i < end; i++) node->bbox += sorted_prims[i].bbox;

    // decide whether to create a leaf
    if (end - start <= bvh_minprims) {
        // makes a leaf node
        node->type = bvh->type;
        node->start = start;
        node->count = end - start;
    } else {
        // choose the split axis and position
        // init to default values
        auto axis = 0;
        auto mid = (start + end) / 2;

        // compute primintive bounds and size
        auto centroid_bbox = invalid_bbox3f;
        for (auto i = start; i < end; i++)
            centroid_bbox += sorted_prims[i].center;
        auto centroid_size = bbox_diagonal(centroid_bbox);

        // check if it is not possible to split
        if (centroid_size == zero3f) {
            // we failed to split for some reasons
            node->type = bvh->type;
            node->start = start;
            node->count = end - start;
        } else {
            // split along largest
            auto largest_axis = max_element(centroid_size);

            // check heuristic
            if (equalsize) {
                // split the space in the middle along the largest axis
                axis = largest_axis;
                mid = (int)(std::partition(sorted_prims + start,
                                sorted_prims + end,
                                bvh_bound_prim_comp(largest_axis,
                                    bbox_center(centroid_bbox)[largest_axis])) -
                            sorted_prims);
            } else {
                // balanced tree split: find the largest axis of the bounding
                // box and split along this one right in the middle
                axis = largest_axis;
                mid = (start + end) / 2;
                std::nth_element(sorted_prims + start, sorted_prims + mid,
                    sorted_prims + end, bvh_bound_prim_comp(largest_axis));
            }

            // check correctness
            assert(axis >= 0 && mid > 0);
            assert(mid > start && mid < end);

            // makes an internal node
            node->type = bvh_node_type::internal;
            // perform the splits by preallocating the child nodes and recurring
            node->axis = axis;
            node->start = (int)bvh->nodes.size();
            node->count = 2;
            bvh->nodes.emplace_back();
            bvh->nodes.emplace_back();
            // build child nodes
            make_bvh_node(
                bvh, node->start, sorted_prims, start, mid, equalsize);
            make_bvh_node(
                bvh, node->start + 1, sorted_prims, mid, end, equalsize);
        }
    }
}

// Build a BVH from a set of primitives.
bvh_tree* make_bvh(const vector<int>& points, const vector<vec2i>& lines,
    const vector<vec3i>& triangles, const vector<vec4i>& quads,
    const vector<vec3f>& pos, const vector<float>& radius, float def_radius,
    bool equalsize) {
    // allocate the bvh
    auto bvh = new bvh_tree();

    // set values
    bvh->points = points;
    bvh->lines = lines;
    bvh->triangles = triangles;
    bvh->quads = quads;
    bvh->pos = pos;
    bvh->radius =
        (radius.empty()) ? vector<float>(pos.size(), def_radius) : radius;

    // get the number of primitives and the primitive type
    auto bboxes = vector<bbox3f>();
    if (!bvh->points.empty()) {
        bboxes.reserve(points.size());
        for (auto& p : bvh->points) {
            bboxes.push_back(point_bbox(bvh->pos[p], bvh->radius[p]));
        }
        bvh->type = bvh_node_type::point;
    } else if (!bvh->lines.empty()) {
        bboxes.reserve(bvh->lines.size());
        for (auto& l : bvh->lines) {
            bboxes.push_back(line_bbox(bvh->pos[l.x], bvh->pos[l.y],
                bvh->radius[l.x], bvh->radius[l.y]));
        }
        bvh->type = bvh_node_type::line;
    } else if (!bvh->triangles.empty()) {
        bboxes.reserve(bvh->triangles.size());
        for (auto& t : bvh->triangles) {
            bboxes.push_back(
                triangle_bbox(bvh->pos[t.x], bvh->pos[t.y], bvh->pos[t.z]));
        }
        bvh->type = bvh_node_type::triangle;
    } else if (!bvh->lines.empty()) {
        bboxes.reserve(quads.size());
        for (auto& q : bvh->quads) {
            bboxes.push_back(quad_bbox(
                bvh->pos[q.x], bvh->pos[q.y], bvh->pos[q.z], bvh->pos[q.w]));
        }
        bvh->type = bvh_node_type::quad;
    } else if (!bvh->lines.empty()) {
        bboxes.reserve(lines.size());
        for (auto i = 0; i < bvh->pos.size(); i++) {
            bboxes.push_back(point_bbox(bvh->pos[i], bvh->radius[i]));
        }
        bvh->type = bvh_node_type::vertex;
    }

    // create buonded primitived for sorting
    auto bound_prims = vector<bvh_bound_prim>(bboxes.size());
    for (auto i = 0; i < bboxes.size(); i++) {
        bound_prims[i].pid = i;
        bound_prims[i].bbox = bboxes[i];
        bound_prims[i].center = bbox_center(bboxes[i]);
    }

    // clear bvh
    bvh->nodes.clear();
    bvh->sorted_prim.clear();

    // allocate nodes (over-allocate now then shrink)
    bvh->nodes.reserve(bound_prims.size() * 2);

    // start recursive splitting
    bvh->nodes.emplace_back();
    make_bvh_node(
        bvh, 0, bound_prims.data(), 0, (int)bound_prims.size(), equalsize);

    // shrink back
    bvh->nodes.shrink_to_fit();

    // init sorted element arrays
    // for shared memory, stored pointer to the external data
    // store the sorted primitive order for BVH walk
    bvh->sorted_prim.resize(bound_prims.size());
    for (int i = 0; i < bound_prims.size(); i++) {
        bvh->sorted_prim[i] = bound_prims[i].pid;
    }

    // done
    return bvh;
}

// Build a BVH from a set of shape instances.
bvh_tree* make_bvh(const vector<frame3f>& frames,
    const vector<frame3f>& frames_inv, const vector<int>& ist_bvhs,
    const vector<bvh_tree*>& shape_bvhs, bool own_shape_bvhs, bool equal_size) {
    // allocate the bvh
    auto bvh = new bvh_tree();

    // set values
    bvh->ist_frames = frames;
    bvh->ist_frames_inv = frames_inv;
    bvh->ist_bvhs = ist_bvhs;
    bvh->shape_bvhs = shape_bvhs;
    bvh->own_shape_bvhs = own_shape_bvhs;

    // get the number of primitives and the primitive type
    auto bboxes = vector<bbox3f>();
    bboxes.reserve(bvh->ist_bvhs.size());
    for (auto idx = 0; idx < bvh->ist_bvhs.size(); idx++) {
        auto sbvh = bvh->shape_bvhs[bvh->ist_bvhs[idx]];
        bboxes.push_back(
            transform_bbox(bvh->ist_frames[idx], sbvh->nodes[0].bbox));
    }
    bvh->type = bvh_node_type::instance;

    // create buonded primitived for sorting
    auto bound_prims = vector<bvh_bound_prim>(bboxes.size());
    for (auto i = 0; i < bboxes.size(); i++) {
        bound_prims[i].pid = i;
        bound_prims[i].bbox = bboxes[i];
        bound_prims[i].center = bbox_center(bboxes[i]);
    }

    // clear bvh
    bvh->nodes.clear();
    bvh->sorted_prim.clear();

    // allocate nodes (over-allocate now then shrink)
    bvh->nodes.reserve(bound_prims.size() * 2);

    // start recursive splitting
    bvh->nodes.emplace_back();
    make_bvh_node(
        bvh, 0, bound_prims.data(), 0, (int)bound_prims.size(), equal_size);

    // shrink back
    bvh->nodes.shrink_to_fit();

    // init sorted element arrays
    // for shared memory, stored pointer to the external data
    // store the sorted primitive order for BVH walk
    bvh->sorted_prim.resize(bound_prims.size());
    for (int i = 0; i < bound_prims.size(); i++) {
        bvh->sorted_prim[i] = bound_prims[i].pid;
    }

    // done
    return bvh;
}

// Recursively recomputes the node bounds for a shape bvh
void refit_bvh(bvh_tree* bvh, int nodeid) {
    // refit
    auto node = &bvh->nodes[nodeid];
    node->bbox = invalid_bbox3f;
    switch (node->type) {
        case bvh_node_type::internal: {
            for (auto i = 0; i < node->count; i++) {
                auto idx = node->start + i;
                refit_bvh(bvh, idx);
                node->bbox += bvh->nodes[idx].bbox;
            }
        } break;
        case bvh_node_type::point: {
            for (auto i = 0; i < node->count; i++) {
                auto idx = bvh->sorted_prim[node->start + i];
                auto& p = bvh->points[idx];
                node->bbox += point_bbox(bvh->pos[p], bvh->radius[p]);
            }
        } break;
        case bvh_node_type::line: {
            for (auto i = 0; i < node->count; i++) {
                auto idx = bvh->sorted_prim[node->start + i];
                auto& l = bvh->lines[idx];
                node->bbox += line_bbox(bvh->pos[l.x], bvh->pos[l.y],
                    bvh->radius[l.x], bvh->radius[l.y]);
            }
        } break;
        case bvh_node_type::triangle: {
            for (auto i = 0; i < node->count; i++) {
                auto idx = bvh->sorted_prim[node->start + i];
                auto& t = bvh->triangles[idx];
                node->bbox +=
                    triangle_bbox(bvh->pos[t.x], bvh->pos[t.y], bvh->pos[t.z]);
            }
        } break;
        case bvh_node_type::quad: {
            for (auto i = 0; i < node->count; i++) {
                auto idx = bvh->sorted_prim[node->start + i];
                auto& q = bvh->quads[idx];
                node->bbox += quad_bbox(
                    bvh->pos[q.x], bvh->pos[q.y], bvh->pos[q.z], bvh->pos[q.w]);
            }
        } break;
        case bvh_node_type::vertex: {
            for (auto i = 0; i < node->count; i++) {
                auto idx = bvh->sorted_prim[node->start + i];
                node->bbox += point_bbox(bvh->pos[idx], bvh->radius[idx]);
            }
        } break;
        case bvh_node_type::instance: {
            for (auto i = 0; i < node->count; i++) {
                auto idx = bvh->sorted_prim[node->start + i];
                auto sbvh = bvh->shape_bvhs[bvh->ist_bvhs[idx]];
                node->bbox +=
                    transform_bbox(bvh->ist_frames[idx], sbvh->nodes[0].bbox);
            }
        } break;
    }
}

// Recursively recomputes the node bounds for a shape bvh
void refit_bvh(bvh_tree* bvh, const vector<vec3f>& pos,
    const vector<float>& radius, float def_radius) {
    // set values
    bvh->pos = pos;
    bvh->radius =
        (radius.empty()) ? vector<float>(pos.size(), def_radius) : radius;

    // refit
    refit_bvh(bvh, 0);
}

// Recursively recomputes the node bounds for a scene bvh
void refit_bvh(bvh_tree* bvh, const vector<frame3f>& frames,
    const vector<frame3f>& frames_inv) {
    // set values
    bvh->ist_frames = frames;
    bvh->ist_frames_inv = frames_inv;

    // refit
    refit_bvh(bvh, 0);
}

// Intersect ray with a bvh.
bool intersect_bvh(const bvh_tree* bvh, const ray3f& ray_, bool early_exit,
    float& ray_t, int& iid, int& eid, vec4f& ew) {
    // node stack
    int node_stack[64];
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
    auto ray_reverse = array<bool, 4>{
        {(bool)ray_dsign.x, (bool)ray_dsign.y, (bool)ray_dsign.z, false}};

    // walking stack
    while (node_cur) {
        // grab node
        auto node = bvh->nodes[node_stack[--node_cur]];

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
                    for (auto i = 0; i < node.count; i++) {
                        auto idx = node.start + i;
                        node_stack[node_cur++] = idx;
                        assert(node_cur < 64);
                    }
                } else {
                    for (auto i = node.count - 1; i >= 0; i--) {
                        auto idx = node.start + i;
                        node_stack[node_cur++] = idx;
                        assert(node_cur < 64);
                    }
                }
            } break;
            case bvh_node_type::point: {
                for (auto i = 0; i < node.count; i++) {
                    auto idx = bvh->sorted_prim[node.start + i];
                    auto& p = bvh->points[idx];
                    if (intersect_point(
                            ray, bvh->pos[p], bvh->radius[p], ray_t)) {
                        hit = true;
                        ray.tmax = ray_t;
                        eid = idx;
                        ew = {1, 0, 0, 0};
                    }
                }
            } break;
            case bvh_node_type::line: {
                for (auto i = 0; i < node.count; i++) {
                    auto idx = bvh->sorted_prim[node.start + i];
                    auto& l = bvh->lines[idx];
                    if (intersect_line(ray, bvh->pos[l.x], bvh->pos[l.y],
                            bvh->radius[l.x], bvh->radius[l.y], ray_t,
                            (vec2f&)ew)) {
                        hit = true;
                        ray.tmax = ray_t;
                        eid = idx;
                        ew = {ew.x, ew.y, 0, 0};
                    }
                }
            } break;
            case bvh_node_type::triangle: {
                for (auto i = 0; i < node.count; i++) {
                    auto idx = bvh->sorted_prim[node.start + i];
                    auto& t = bvh->triangles[idx];
                    if (intersect_triangle(ray, bvh->pos[t.x], bvh->pos[t.y],
                            bvh->pos[t.z], ray_t, (vec3f&)ew)) {
                        hit = true;
                        ray.tmax = ray_t;
                        eid = idx;
                        ew = {ew.x, ew.y, ew.z, 0};
                    }
                }
            } break;
            case bvh_node_type::quad: {
                for (auto i = 0; i < node.count; i++) {
                    auto idx = bvh->sorted_prim[node.start + i];
                    auto& q = bvh->quads[idx];
                    if (intersect_quad(ray, bvh->pos[q.x], bvh->pos[q.y],
                            bvh->pos[q.z], bvh->pos[q.w], ray_t, ew)) {
                        hit = true;
                        ray.tmax = ray_t;
                        eid = idx;
                        ew = {ew.x, ew.y, ew.z, ew.w};
                    }
                }
            } break;
            case bvh_node_type::vertex: {
                for (auto i = 0; i < node.count; i++) {
                    auto idx = bvh->sorted_prim[node.start + i];
                    if (intersect_point(
                            ray, bvh->pos[idx], bvh->radius[idx], ray_t)) {
                        hit = true;
                        ray.tmax = ray_t;
                        eid = idx;
                        ew = {1, 0, 0, 0};
                    }
                }
            } break;
            case bvh_node_type::instance: {
                for (auto i = 0; i < node.count; i++) {
                    auto idx = bvh->sorted_prim[node.start + i];
                    auto sbvh = bvh->shape_bvhs[bvh->ist_bvhs[idx]];
                    if (intersect_bvh(sbvh,
                            transform_ray(bvh->ist_frames_inv[idx], ray),
                            early_exit, ray_t, iid, eid, ew)) {
                        hit = true;
                        ray.tmax = ray_t;
                        iid = idx;
                    }
                }
            } break;
        }

        // check for early exit
        if (early_exit && hit) return true;
    }

    return hit;
}

// Finds the closest element with a bvh.
bool overlap_bvh(const bvh_tree* bvh, const vec3f& pos, float max_dist,
    bool early_exit, float& dist, int& iid, int& eid, vec4f& ew) {
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
                for (auto idx = node.start; idx < node.start + node.count;
                     idx++) {
                    node_stack[node_cur++] = idx;
                    assert(node_cur < 64);
                }
            } break;
            case bvh_node_type::point: {
                for (auto i = 0; i < node.count; i++) {
                    auto idx = bvh->sorted_prim[node.start + i];
                    auto& p = bvh->points[idx];
                    if (overlap_point(
                            pos, max_dist, bvh->pos[p], bvh->radius[p], dist)) {
                        hit = true;
                        max_dist = dist;
                        eid = idx;
                        ew = {1, 0, 0, 0};
                    }
                }
            } break;
            case bvh_node_type::line: {
                for (auto i = 0; i < node.count; i++) {
                    auto idx = bvh->sorted_prim[node.start + i];
                    auto& l = bvh->lines[idx];
                    if (overlap_line(pos, max_dist, bvh->pos[l.x],
                            bvh->pos[l.y], bvh->radius[l.x], bvh->radius[l.y],
                            dist, (vec2f&)ew)) {
                        hit = true;
                        max_dist = dist;
                        eid = idx;
                        ew = {ew.x, ew.y, 0, 0};
                    }
                }
            } break;
            case bvh_node_type::triangle: {
                for (auto i = 0; i < node.count; i++) {
                    auto idx = bvh->sorted_prim[node.start + i];
                    auto& t = bvh->triangles[idx];
                    if (overlap_triangle(pos, max_dist, bvh->pos[t.x],
                            bvh->pos[t.y], bvh->pos[t.z], bvh->radius[t.x],
                            bvh->radius[t.y], bvh->radius[t.z], dist,
                            (vec3f&)ew)) {
                        hit = true;
                        max_dist = dist;
                        eid = idx;
                        ew = {ew.x, ew.y, ew.z, 0};
                    }
                }
            } break;
            case bvh_node_type::quad: {
                for (auto i = 0; i < node.count; i++) {
                    auto idx = bvh->sorted_prim[node.start + i];
                    auto& q = bvh->quads[idx];
                    if (overlap_quad(pos, max_dist, bvh->pos[q.x],
                            bvh->pos[q.y], bvh->pos[q.z], bvh->pos[q.w],
                            bvh->radius[q.x], bvh->radius[q.y],
                            bvh->radius[q.z], bvh->radius[q.w], dist, ew)) {
                        hit = true;
                        max_dist = dist;
                        eid = idx;
                        ew = {ew.x, ew.y, ew.z, ew.w};
                    }
                }
            } break;
            case bvh_node_type::vertex: {
                for (auto i = 0; i < node.count; i++) {
                    auto idx = bvh->sorted_prim[node.start + i];
                    if (overlap_point(pos, max_dist, bvh->pos[idx],
                            bvh->radius[idx], dist)) {
                        hit = true;
                        max_dist = dist;
                        eid = idx;
                        ew = {1, 0, 0, 0};
                    }
                }
            } break;
            case bvh_node_type::instance: {
                for (auto i = 0; i < node.count; i++) {
                    auto idx = bvh->sorted_prim[node.start + i];
                    auto sbvh = bvh->shape_bvhs[bvh->ist_bvhs[idx]];
                    if (overlap_bvh(sbvh,
                            transform_point(bvh->ist_frames_inv[idx], pos),
                            max_dist, early_exit, dist, iid, eid, ew)) {
                        hit = true;
                        max_dist = dist;
                        iid = idx;
                    }
                }
            } break;
        }
    }

    return hit;
}

// Intersect ray with a bvh (convenience wrapper).
intersection_point intersect_bvh(
    const bvh_tree* bvh, const ray3f& ray, bool early_exit) {
    auto isec = intersection_point();
    if (!intersect_bvh(
            bvh, ray, early_exit, isec.dist, isec.iid, isec.eid, isec.euv))
        return {};
    return isec;
}

// Finds the closest element with a bvh (convenience wrapper).
intersection_point overlap_bvh(
    const bvh_tree* bvh, const vec3f& pos, float max_dist, bool early_exit) {
    auto isec = intersection_point();
    if (!overlap_bvh(bvh, pos, max_dist, early_exit, isec.dist, isec.iid,
            isec.eid, isec.euv))
        return {};
    return isec;
}

#if 0
    // Finds the overlap between BVH leaf nodes.
    template <typename OverlapElem>
    void overlap_bvh_elems(const bvh_tree* bvh1, const bvh_tree* bvh2,
                           bool skip_duplicates, bool skip_self, vector<vec2i>& overlaps,
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
animation::~animation() {
    for (auto v : keyframes) delete v;
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
T eval_barycentric(
    const shape* shp, const vector<T>& vals, int eid, const vec4f& euv) {
    if (vals.empty()) return T();
    if (!shp->triangles.empty()) {
        return eval_barycentric_triangle(
            vals, shp->triangles[eid], vec3f{euv.x, euv.y, euv.z});
    } else if (!shp->lines.empty()) {
        return eval_barycentric_line(
            vals, shp->lines[eid], vec2f{euv.x, euv.y});
    } else if (!shp->points.empty()) {
        return eval_barycentric_point(vals, shp->points[eid], euv.x);
    } else if (!shp->quads.empty()) {
        return eval_barycentric_quad(vals, shp->quads[eid], euv);
    } else {
        return vals[eid];  // points
    }
}

// Shape position interpolated using barycentric coordinates
vec3f eval_pos(const shape* shp, int eid, const vec4f& euv) {
    return eval_barycentric(shp, shp->pos, eid, euv);
}

// Shape normal interpolated using barycentric coordinates
vec3f eval_norm(const shape* shp, int eid, const vec4f& euv) {
    return normalize(eval_barycentric(shp, shp->norm, eid, euv));
}

// Shape texcoord interpolated using barycentric coordinates
vec2f eval_texcoord(const shape* shp, int eid, const vec4f& euv) {
    return eval_barycentric(shp, shp->texcoord, eid, euv);
}

// Shape texcoord interpolated using barycentric coordinates
vec4f eval_color(const shape* shp, int eid, const vec4f& euv) {
    return eval_barycentric(shp, shp->color, eid, euv);
}

// Shape tangent space interpolated using barycentric coordinates
vec4f eval_tangsp(const shape* shp, int eid, const vec4f& euv) {
    return eval_barycentric(shp, shp->tangsp, eid, euv);
}

// Instance position interpolated using barycentric coordinates
vec3f eval_pos(const instance* ist, int eid, const vec4f& euv) {
    return transform_point(
        ist->frame, eval_barycentric(ist->shp, ist->shp->pos, eid, euv));
}

// Instance normal interpolated using barycentric coordinates
vec3f eval_norm(const instance* ist, int eid, const vec4f& euv) {
    return transform_direction(ist->frame,
        normalize(eval_barycentric(ist->shp, ist->shp->norm, eid, euv)));
}

// Evaluate a texture
vec4f eval_texture(const texture_info& info, const vec2f& texcoord, bool srgb,
    const vec4f& def) {
    if (!info.txt) return def;

    // get texture
    auto txt = info.txt;
    assert(txt->hdr || txt->ldr);

    auto lookup = [&def, &txt, &srgb](int i, int j) {
        if (txt->ldr)
            return (srgb) ? srgb_to_linear(txt->ldr[{i, j}]) :
                            byte_to_float(txt->ldr[{i, j}]);
        else if (txt->hdr)
            return txt->hdr[{i, j}];
        else
            return def;
    };

    // get image width/height
    auto w = txt->width(), h = txt->height();

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

// Subdivides shape elements.
void subdivide_shape_once(shape* shp, bool subdiv) {
    if (!shp->lines.empty() || !shp->triangles.empty() || !shp->quads.empty()) {
        vector<vec2i> edges;
        vector<vec4i> faces;
        tie(shp->lines, shp->triangles, shp->quads, edges, faces) =
            subdivide_elems_linear(
                shp->lines, shp->triangles, shp->quads, (int)shp->pos.size());
        shp->pos = subdivide_vert_linear(shp->pos, edges, faces);
        shp->norm = subdivide_vert_linear(shp->norm, edges, faces);
        shp->texcoord = subdivide_vert_linear(shp->texcoord, edges, faces);
        shp->color = subdivide_vert_linear(shp->color, edges, faces);
        shp->radius = subdivide_vert_linear(shp->radius, edges, faces);
        if (subdiv && !shp->quads.empty()) {
            auto boundary = get_boundary_edges({}, {}, shp->quads);
            shp->pos =
                subdivide_vert_catmullclark(shp->quads, shp->pos, boundary, {});
            shp->norm = subdivide_vert_catmullclark(
                shp->quads, shp->norm, boundary, {});
            shp->texcoord = subdivide_vert_catmullclark(
                shp->quads, shp->texcoord, boundary, {});
            shp->color = subdivide_vert_catmullclark(
                shp->quads, shp->color, boundary, {});
            shp->radius = subdivide_vert_catmullclark(
                shp->quads, shp->radius, boundary, {});
            shp->norm = compute_normals({}, {}, shp->quads, shp->pos);
        }
    } else if (!shp->quads_pos.empty()) {
        vector<vec2i> _lines;
        vector<vec3i> _triangles;
        vector<vec2i> edges;
        vector<vec4i> faces;
        tie(_lines, _triangles, shp->quads_pos, edges, faces) =
            subdivide_elems_linear({}, {}, shp->quads_pos, shp->pos.size());
        shp->pos = subdivide_vert_linear(shp->pos, edges, faces);
        tie(_lines, _triangles, shp->quads_norm, edges, faces) =
            subdivide_elems_linear({}, {}, shp->quads_norm, shp->norm.size());
        shp->norm = subdivide_vert_linear(shp->norm, edges, faces);
        tie(_lines, _triangles, shp->quads_texcoord, edges, faces) =
            subdivide_elems_linear(
                {}, {}, shp->quads_texcoord, shp->texcoord.size());
        shp->texcoord = subdivide_vert_linear(shp->texcoord, edges, faces);
        if (subdiv) {
            shp->pos = subdivide_vert_catmullclark(shp->quads_pos, shp->pos,
                get_boundary_edges({}, {}, shp->quads_pos), {});
            shp->norm = subdivide_vert_catmullclark(shp->quads_norm, shp->norm,
                get_boundary_edges({}, {}, shp->quads_norm), {});
            shp->texcoord =
                subdivide_vert_catmullclark(shp->quads_texcoord, shp->texcoord,
                    {}, get_boundary_verts({}, {}, shp->quads_texcoord));
        }
    } else if (!shp->beziers.empty()) {
        vector<int> verts;
        vector<vec4i> segments;
        tie(shp->beziers, verts, segments) =
            subdivide_bezier_recursive(shp->beziers, (int)shp->pos.size());
        shp->pos = subdivide_vert_bezier(shp->pos, verts, segments);
        shp->norm = subdivide_vert_bezier(shp->norm, verts, segments);
        shp->texcoord = subdivide_vert_bezier(shp->texcoord, verts, segments);
        shp->color = subdivide_vert_bezier(shp->color, verts, segments);
        shp->radius = subdivide_vert_bezier(shp->radius, verts, segments);
    }
}

// Facet a shape. Supports only non-facevarying shapes
void facet_shape(shape* shp, bool recompute_normals) {
    if (!shp->lines.empty() || !shp->triangles.empty() || !shp->quads.empty()) {
        vector<int> verts;
        tie(shp->lines, shp->triangles, shp->quads, verts) =
            facet_elems(shp->lines, shp->triangles, shp->quads);
        shp->pos = facet_vert(shp->pos, verts);
        shp->norm = facet_vert(shp->norm, verts);
        shp->texcoord = facet_vert(shp->texcoord, verts);
        shp->color = facet_vert(shp->color, verts);
        shp->radius = facet_vert(shp->radius, verts);
        if (recompute_normals) {
            shp->norm = compute_normals(
                shp->lines, shp->triangles, shp->quads, shp->pos);
        }
    }
}

// Tesselate a shape into basic primitives
void tesselate_shape(shape* shp, bool subdivide,
    bool facevarying_to_sharedvertex, bool quads_to_triangles,
    bool bezier_to_lines) {
    if (subdivide && shp->subdivision_level) {
        for (auto l = 0; l < shp->subdivision_level; l++) {
            subdivide_shape_once(shp, shp->subdivision_catmullclark);
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
    for (auto shp : scn->shapes) {
        tesselate_shape(shp, subdivide, facevarying_to_sharedvertex,
            quads_to_triangles, bezier_to_lines);
    }
}

// Update animation transforms
void update_transforms(animation* anm, float time) {
    auto interpolate = [](keyframe_type type, const vector<float>& times,
                           const auto& vals, float time) {
        switch (type) {
            case keyframe_type::step:
                return eval_keyframed_step(times, vals, time);
            case keyframe_type::linear:
                return eval_keyframed_linear(times, vals, time);
            case keyframe_type::catmull_rom: return vals.at(0);
            case keyframe_type::bezier:
                return eval_keyframed_bezier(times, vals, time);
            default: throw runtime_error("should not have been here");
        }
        return vals.at(0);
    };

    for (auto kfr : anm->keyframes) {
        if (!kfr->translation.empty()) {
            auto val =
                interpolate(kfr->type, kfr->times, kfr->translation, time);
            for (auto target : anm->targets)
                if (target.first == kfr) target.second->translation = val;
        }
        if (!kfr->rotation.empty()) {
            auto val = interpolate(kfr->type, kfr->times, kfr->rotation, time);
            for (auto target : anm->targets)
                if (target.first == kfr) target.second->rotation = val;
        }
        if (!kfr->scaling.empty()) {
            auto val = interpolate(kfr->type, kfr->times, kfr->scaling, time);
            for (auto target : anm->targets)
                if (target.first == kfr) target.second->scaling = val;
        }
    }
}

// Update node transforms
void update_transforms(node* nde, const frame3f& parent = identity_frame3f) {
    auto frame = parent * nde->frame * translation_frame3(nde->translation) *
                 rotation_frame3(nde->rotation) * scaling_frame3(nde->scaling);
    for (auto ist : nde->ists) ist->frame = frame;
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
    for (auto anm : scn->animations) {
        for (auto kfr : anm->keyframes) {
            range.x = min(range.x, kfr->times.front());
            range.y = max(range.y, kfr->times.back());
        }
    }
    return range;
}

// Add missing values and elements
void add_elements(scene* scn, const add_elements_options& opts) {
    if (opts.smooth_normals) {
        for (auto shp : scn->shapes) {
            if (!shp->norm.empty()) continue;
            shp->norm.resize(shp->pos.size(), {0, 0, 1});
            if (!shp->lines.empty() || !shp->triangles.empty() ||
                !shp->quads.empty()) {
                shp->norm = compute_normals(
                    shp->lines, shp->triangles, shp->quads, shp->pos);
            }
            if (!shp->quads_pos.empty()) {
                if (!shp->quads_norm.empty())
                    throw runtime_error("bad normals");
                shp->quads_norm = shp->quads_pos;
                shp->norm = compute_normals({}, {}, shp->quads_pos, shp->pos);
            }
        }
    }

    if (opts.tangent_space) {
        for (auto shp : scn->shapes) {
            if (!shp->tangsp.empty() || shp->texcoord.empty() || !shp->mat ||
                !(shp->mat->norm_txt.txt || shp->mat->bump_txt.txt))
                continue;
            if (!shp->triangles.empty()) {
                shp->tangsp = compute_tangent_frames(
                    shp->triangles, shp->pos, shp->norm, shp->texcoord);
            } else if (!shp->quads.empty()) {
                auto triangles = convert_quads_to_triangles(shp->quads);
                shp->tangsp = compute_tangent_frames(
                    triangles, shp->pos, shp->norm, shp->texcoord);
            }
        }
    }

    if (opts.pointline_radius > 0) {
        for (auto shp : scn->shapes) {
            if ((shp->points.empty() && shp->lines.empty()) ||
                !shp->radius.empty())
                continue;
            shp->radius.resize(shp->pos.size(), opts.pointline_radius);
        }
    }

    if (opts.texture_data) {
        for (auto txt : scn->textures) {
            if (!txt->hdr && !txt->ldr) {
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
    }

    // default environment
    if (opts.default_environment && scn->environments.empty()) {
        auto env = new environment();
        env->name = "default_environment";
        scn->environments.push_back(env);
    }
}

// Make a view camera either copying a given one or building a default one.
camera* make_view_camera(const scene* scn, int camera_id) {
    if (scn->cameras.empty()) {
        auto bbox = compute_bounds(scn);
        auto bbox_center = (bbox.max + bbox.min) / 2.0f;
        auto bbox_size = bbox.max - bbox.min;
        auto bbox_msize = max(bbox_size[0], max(bbox_size[1], bbox_size[2]));
        // set up camera
        auto cam = new camera();
        cam->name = "default_camera";
        auto camera_dir = vec3f{1, 0.4f, 1};
        auto from = camera_dir * bbox_msize + bbox_center;
        auto to = bbox_center;
        auto up = vec3f{0, 1, 0};
        cam->frame = lookat_frame3(from, to, up);
        cam->ortho = false;
        cam->aspect = 16.0f / 9.0f;
        cam->yfov = 2 * atanf(0.5f);
        cam->aperture = 0;
        cam->focus = length(to - from);
        return cam;
    } else {
        camera_id = clamp(camera_id, 0, (int)scn->cameras.size());
        return new camera(*scn->cameras[camera_id]);
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

// Computes a shape bounding box (quick computation that ignores radius)
bbox3f compute_bounds(const shape* shp) {
    auto bbox = invalid_bbox3f;
    for (auto p : shp->pos) bbox += vec3f(p);
    return bbox;
}

// Updates the instance bounding box
bbox3f compute_bounds(const instance* ist) {
    return transform_bbox(ist->frame, compute_bounds(ist->shp));
}

// Updates the scene and scene's instances bounding boxes
bbox3f compute_bounds(const scene* scn) {
    auto shape_bboxes = unordered_map<shape*, bbox3f>();
    for (auto shp : scn->shapes) shape_bboxes[shp] = compute_bounds(shp);
    auto bbox = invalid_bbox3f;
    if (!scn->instances.empty()) {
        for (auto ist : scn->instances) {
            bbox += transform_bbox(ist->frame, shape_bboxes.at(ist->shp));
        }
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
        auto nshp = new shape(*ist->shp);
        for (auto& p : nshp->pos) p = transform_point(ist->frame, p);
        for (auto& n : nshp->norm) n = transform_direction(ist->frame, n);
        scn->shapes.push_back(nshp);
    }
    for (auto e : shapes) delete e;
    for (auto e : instances) delete e;
}

// Build a shape BVH
bvh_tree* make_bvh(const shape* shp, float def_radius, bool equalsize) {
    return make_bvh(shp->points, shp->lines, shp->triangles, shp->quads,
        shp->pos, shp->radius, def_radius, equalsize);
}

// Build a scene BVH
bvh_tree* make_bvh(const scene* scn, float def_radius, bool equalsize) {
    // do shapes
    auto shape_bvhs = vector<bvh_tree*>(scn->shapes.size());
    auto smap = unordered_map<shape*, int>();
    for (auto sid = 0; sid < scn->shapes.size(); sid++) {
        auto shp = scn->shapes[sid];
        smap[shp] = sid;
        shape_bvhs[sid] = make_bvh(shp, def_radius, equalsize);
    }

    // tree bvh
    auto ist_frames = vector<frame3f>(scn->instances.size());
    auto ist_frames_inv = vector<frame3f>(scn->instances.size());
    auto ist_bvh = vector<int>(scn->instances.size());
    for (auto iid = 0; iid < scn->instances.size(); iid++) {
        auto ist = scn->instances[iid];
        ist_frames[iid] = ist->frame;
        ist_frames_inv[iid] = inverse(ist->frame);
        ist_bvh[iid] = smap.at(ist->shp);
    }
    return make_bvh(
        ist_frames, ist_frames_inv, ist_bvh, shape_bvhs, true, equalsize);
}

// Refits a scene BVH
void refit_bvh(bvh_tree* bvh, const shape* shp, float def_radius) {
    refit_bvh(bvh, shp->pos, shp->radius, def_radius);
}

// Refits a scene BVH
void refit_bvh(
    bvh_tree* bvh, const scene* scn, bool do_shapes, float def_radius) {
    if (do_shapes) {
        for (auto sid = 0; sid < scn->shapes.size(); sid++) {
            refit_bvh(get_shape_bvhs(bvh).at(sid), scn->shapes[sid]->pos,
                scn->shapes[sid]->radius, def_radius);
        }
    }
    auto ist_frames = vector<frame3f>();
    auto ist_frames_inv = vector<frame3f>();
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
    for (auto shp : scn->shapes) {
        nverts += shp->pos.size();
        nnorms += shp->norm.size();
        ntexcoords += shp->texcoord.size();
        npoints += shp->points.size();
        nlines += shp->lines.size();
        ntriangles += shp->triangles.size();
        nquads += shp->quads.size();
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

    struct obj_obj_vertex_hash {
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
    auto tmap = unordered_map<string, texture*>{{"", nullptr}};
    for (auto otxt : obj->textures) {
        auto txt = new texture();
        txt->name = otxt->path;
        txt->path = otxt->path;
        if (!otxt->datab.empty()) {
            txt->ldr = image4b(otxt->width, otxt->height);
            for (auto j = 0; j < otxt->height; j++) {
                for (auto i = 0; i < otxt->width; i++) {
                    auto v = otxt->datab.data() +
                             (otxt->width * j + i) * otxt->ncomp;
                    switch (otxt->ncomp) {
                        case 1:
                            txt->ldr.at(i, j) = {v[0], v[0], v[0], 255};
                            break;
                        case 2: txt->ldr.at(i, j) = {v[0], v[1], 0, 255}; break;
                        case 3:
                            txt->ldr.at(i, j) = {v[0], v[1], v[2], 255};
                            break;
                        case 4:
                            txt->ldr.at(i, j) = {v[0], v[1], v[2], v[3]};
                            break;
                        default: assert(false); break;
                    }
                }
            }
        } else if (!otxt->dataf.empty()) {
            txt->hdr = image4f(otxt->width, otxt->height);
            for (auto j = 0; j < otxt->height; j++) {
                for (auto i = 0; i < otxt->width; i++) {
                    auto v = otxt->dataf.data() +
                             (otxt->width * j + i) * otxt->ncomp;
                    switch (otxt->ncomp) {
                        case 1:
                            txt->hdr.at(i, j) = {v[0], v[0], v[0], 1};
                            break;
                        case 2: txt->hdr.at(i, j) = {v[0], v[1], 0, 1}; break;
                        case 3:
                            txt->hdr.at(i, j) = {v[0], v[1], v[2], 1};
                            break;
                        case 4:
                            txt->hdr.at(i, j) = {v[0], v[1], v[2], v[3]};
                            break;
                        default: assert(false); break;
                    }
                }
            }
        }
        scn->textures.push_back(txt);
        tmap[txt->path] = txt;
    }

    auto add_texture = [&tmap](const obj_texture_info& oinfo) {
        auto info = texture_info();
        if (oinfo.path == "") return info;
        info.txt = tmap.at(oinfo.path);
        info.wrap_s = !oinfo.clamp;
        info.wrap_t = !oinfo.clamp;
        info.scale = oinfo.scale;
        return info;
    };

    // convert materials and build textures
    auto mmap = unordered_map<string, material*>{{"", nullptr}};
    for (auto omat : obj->materials) {
        auto mat = new material();
        mat->name = omat->name;
        mat->type = material_type::specular_roughness;
        mat->ke = {omat->ke.x, omat->ke.y, omat->ke.z};
        mat->kd = {omat->kd.x, omat->kd.y, omat->kd.z};
        mat->ks = {omat->ks.x, omat->ks.y, omat->ks.z};
        mat->kr = {omat->kr.x, omat->kr.y, omat->kr.z};
        mat->kt = {omat->kt.x, omat->kt.y, omat->kt.z};
        mat->rs = (omat->ns >= 1e6f) ? 0 : pow(2 / (omat->ns + 2), 1 / 4.0f);
        mat->op = omat->op;
        mat->ke_txt = add_texture(omat->ke_txt);
        mat->kd_txt = add_texture(omat->kd_txt);
        mat->ks_txt = add_texture(omat->ks_txt);
        mat->kr_txt = add_texture(omat->kr_txt);
        mat->kt_txt = add_texture(omat->kt_txt);
        mat->rs_txt = add_texture(omat->ns_txt);
        mat->norm_txt = add_texture(omat->norm_txt);
        mat->bump_txt = add_texture(omat->bump_txt);
        mat->disp_txt = add_texture(omat->disp_txt);
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
    auto omap = unordered_map<string, vector<shape*>>{{"", {}}};
    for (auto omsh : obj->objects) {
        omap[omsh->name] = {};
        for (auto& oshp : omsh->groups) {
            if (oshp.verts.empty()) continue;
            if (oshp.elems.empty()) continue;

            auto shp = new shape();
            shp->name = omsh->name + oshp.groupname;
            shp->mat = mmap[oshp.matname];
            shp->subdivision_level = oshp.subdivision_level;
            shp->subdivision_catmullclark = oshp.subdivision_catmullclark;

            // check to see if this shuold be face-varying or flat quads
            auto as_facevarying = false, as_quads = false;
            if (opts.preserve_quads || opts.preserve_facevarying) {
                auto m = 10000, M = -1;
                for (auto& elem : oshp.elems) {
                    if (elem.type != obj_element_type::face) {
                        m = 2;
                        break;
                    } else {
                        m = min(m, (int)elem.size);
                        M = max(M, (int)elem.size);
                    }
                }
                if (m >= 3 && M == 4) as_quads = opts.preserve_quads;
                if (m >= 3 && M <= 4)
                    as_facevarying = opts.preserve_facevarying;
            }

            // in case of facevarying, check if there is really a need for it
            if (as_facevarying) {
                auto need_facevarying = false;
                for (auto& elem : oshp.elems) {
                    for (auto i = elem.start; i < elem.start + elem.size; i++) {
                        auto& v = oshp.verts[i];
                        if (v.norm >= 0 && v.pos != v.norm)
                            need_facevarying = true;
                        if (v.texcoord >= 0 && v.pos != v.texcoord)
                            need_facevarying = true;
                        if (v.color >= 0) as_facevarying = false;
                        if (v.radius >= 0) as_facevarying = false;
                    }
                    if (!as_facevarying) break;
                }
                as_facevarying = need_facevarying;
            }

            if (!as_facevarying) {
                // insert all vertices
                unordered_map<obj_vertex, int, obj_obj_vertex_hash> vert_map;
                vector<int> vert_ids;
                for (auto& vert : oshp.verts) {
                    if (vert_map.find(vert) == vert_map.end()) {
                        auto s = (int)vert_map.size();
                        vert_map[vert] = s;
                    }
                    vert_ids.push_back(vert_map.at(vert));
                }

                // convert elements
                for (auto& elem : oshp.elems) {
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
                            if (as_quads) {
                                shp->quads.push_back({vert_ids[elem.start + 0],
                                    vert_ids[elem.start + 1],
                                    vert_ids[elem.start + 2],
                                    vert_ids[elem.start +
                                             ((elem.size == 3) ? 2 : 3)]});
                            } else if (elem.size == 3) {
                                shp->triangles.push_back(
                                    {vert_ids[elem.start + 0],
                                        vert_ids[elem.start + 1],
                                        vert_ids[elem.start + 2]});
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
                                throw runtime_error("bad obj bezier");
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
                auto v = oshp.verts[0];
                if (v.pos >= 0) shp->pos.resize(vert_map.size());
                if (v.texcoord >= 0) shp->texcoord.resize(vert_map.size());
                if (v.norm >= 0) shp->norm.resize(vert_map.size());
                if (v.color >= 0) shp->color.resize(vert_map.size());
                if (v.radius >= 0) shp->radius.resize(vert_map.size());
                for (auto& kv : vert_map) {
                    if (v.pos >= 0 && kv.first.pos >= 0) {
                        auto v = obj->pos[kv.first.pos];
                        shp->pos[kv.second] = {v.x, v.y, v.z};
                    }
                    if (v.texcoord >= 0 && kv.first.texcoord >= 0) {
                        auto v = obj->texcoord[kv.first.texcoord];
                        shp->texcoord[kv.second] = {v.x, v.y};
                    }
                    if (v.norm >= 0 && kv.first.norm >= 0) {
                        auto v = obj->norm[kv.first.norm];
                        shp->norm[kv.second] = {v.x, v.y, v.z};
                    }
                    if (v.color >= 0 && kv.first.color >= 0) {
                        auto v = obj->color[kv.first.color];
                        shp->color[kv.second] = {v.x, v.y, v.z, v.w};
                    }
                    if (v.radius >= 0 && kv.first.radius >= 0) {
                        shp->radius[kv.second] = obj->radius[kv.first.radius];
                    }
                }

                // fix smoothing
                if (!oshp.smoothing && opts.obj_facet_non_smooth) {
                    auto faceted = new shape();
                    faceted->name = shp->name;
                    auto pidx = vector<int>();
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
                unordered_map<int, int> pos_map, norm_map, texcoord_map;
                vector<int> pos_ids, norm_ids, texcoord_ids;
                for (auto& vert : oshp.verts) {
                    if (vert.pos >= 0) {
                        if (pos_map.find(vert.pos) == pos_map.end()) {
                            auto s = (int)pos_map.size();
                            pos_map[vert.pos] = s;
                        }
                        pos_ids.push_back(pos_map.at(vert.pos));
                    } else {
                        if (!pos_ids.empty())
                            throw runtime_error("malformed obj");
                    }
                    if (vert.norm >= 0) {
                        if (norm_map.find(vert.norm) == norm_map.end()) {
                            auto s = (int)norm_map.size();
                            norm_map[vert.norm] = s;
                        }
                        norm_ids.push_back(norm_map.at(vert.norm));
                    } else {
                        if (!norm_ids.empty())
                            throw runtime_error("malformed obj");
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
                            throw runtime_error("malformed obj");
                    }
                }

                // convert elements
                for (auto& elem : oshp.elems) {
                    if (elem.type != obj_element_type::face)
                        throw runtime_error("malformed obj");
                    if (elem.size < 3 || elem.size > 4)
                        throw runtime_error("malformed obj");
                    if (!pos_ids.empty()) {
                        shp->quads_pos.push_back({pos_ids[elem.start + 0],
                            pos_ids[elem.start + 1], pos_ids[elem.start + 2],
                            pos_ids[elem.start + ((elem.size == 3) ? 2 : 3)]});
                    }
                    if (!texcoord_ids.empty()) {
                        shp->quads_texcoord.push_back(
                            {texcoord_ids[elem.start + 0],
                                texcoord_ids[elem.start + 1],
                                texcoord_ids[elem.start + 2],
                                texcoord_ids[elem.start +
                                             ((elem.size == 3) ? 2 : 3)]});
                    }
                    if (!norm_ids.empty()) {
                        shp->quads_norm.push_back({norm_ids[elem.start + 0],
                            norm_ids[elem.start + 1], norm_ids[elem.start + 2],
                            norm_ids[elem.start + ((elem.size == 3) ? 2 : 3)]});
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
                if (!oshp.smoothing && opts.obj_facet_non_smooth) {}
            }
            scn->shapes.push_back(shp);
            omap[omsh->name].push_back(shp);
        }
    }

    // convert cameras
    auto cmap = unordered_map<string, camera*>{{"", nullptr}};
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
    unordered_set<material*> env_mat;
    auto emap = unordered_map<string, environment*>{{"", nullptr}};
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
    for (auto shp : scn->shapes) env_mat.erase(shp->mat);
    for (auto mat : env_mat) {
        auto end =
            std::remove(scn->materials.begin(), scn->materials.end(), mat);
        scn->materials.erase(end, scn->materials.end());
        delete mat;
    }

    if (obj->nodes.empty()) {
        // instance cameras
        // instance environments
    } else {
        // convert nodes
        for (auto onde : obj->nodes) {
            auto nde = new node();
            nde->name = onde->name;
            nde->cam = cmap.at(onde->camname);
            if (!onde->objname.empty()) {
                for (auto shp : omap.at(onde->objname)) {
                    auto ist = new instance();
                    ist->name = onde->name;
                    ist->shp = shp;
                    nde->ists.push_back(ist);
                    scn->instances.push_back(ist);
                }
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
scene* load_obj_scene(const string& filename, const load_options& opts) {
    auto oscn = unique_ptr<obj_scene>(load_obj(filename, opts.load_textures,
        opts.skip_missing, opts.obj_flip_texcoord, opts.obj_flip_tr));
    auto scn = unique_ptr<scene>(obj_to_scene(oscn.get(), opts));
    return scn.release();
}

// Save an scene
obj_scene* scene_to_obj(const scene* scn) {
    auto obj = new obj_scene();

    auto add_texture = [](const texture_info& info, bool bump = false) {
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
        if (txt->hdr) {
            otxt->width = txt->hdr.width();
            otxt->height = txt->hdr.height();
            otxt->ncomp = 4;
            otxt->dataf.assign((float*)data(txt->hdr),
                (float*)data(txt->hdr) +
                    txt->hdr.width() * txt->hdr.height() * 4);
        }
        if (txt->ldr) {
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
        omat->ke_txt = add_texture(mat->ke_txt);
        switch (mat->type) {
            case material_type::specular_roughness: {
                omat->kd = {mat->kd.x, mat->kd.y, mat->kd.z};
                omat->ks = {mat->ks.x, mat->ks.y, mat->ks.z};
                omat->kr = {mat->kr.x, mat->kr.y, mat->kr.z};
                omat->kt = {mat->kt.x, mat->kt.y, mat->kt.z};
                omat->ns = (mat->rs) ? 2 / pow(mat->rs, 4.0f) - 2 : 1e6;
                omat->op = mat->op;
                omat->kd_txt = add_texture(mat->kd_txt);
                omat->ks_txt = add_texture(mat->ks_txt);
                omat->kr_txt = add_texture(mat->kr_txt);
                omat->kt_txt = add_texture(mat->kt_txt);
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
                    omat->kd_txt = add_texture(mat->kd_txt);
                } else {
                    omat->ks_txt = add_texture(mat->ks_txt);
                }
            } break;
            case material_type::specular_glossiness: {
                omat->kd = {mat->kd.x, mat->kd.y, mat->kd.z};
                omat->ks = {mat->ks.x, mat->ks.y, mat->ks.z};
                omat->ns = (mat->rs) ? 2 / pow(1 - mat->rs, 4.0f) - 2 : 1e6;
                omat->op = mat->op;
                omat->kd_txt = add_texture(mat->kd_txt);
                omat->ks_txt = add_texture(mat->ks_txt);
            } break;
        }
        omat->bump_txt = add_texture(mat->bump_txt, true);
        omat->disp_txt = add_texture(mat->disp_txt, true);
        omat->norm_txt = add_texture(mat->norm_txt, true);
        if (mat->op < 1 || mat->kt != zero3f) {
            omat->illum = 4;
        } else {
            omat->illum = 2;
        }
        obj->materials.push_back(omat);
    }

    // convert shapes
    auto shapes = unordered_map<shape*, unique_ptr<obj_group>>();
    for (auto shp : scn->shapes) {
        auto offset = obj_vertex{(int)obj->pos.size(),
            (int)obj->texcoord.size(), (int)obj->norm.size(),
            (int)obj->color.size(), (int)obj->radius.size()};
        for (auto& v : shp->pos) obj->pos.push_back({v.x, v.y, v.z});
        for (auto& v : shp->norm) obj->norm.push_back({v.x, v.y, v.z});
        for (auto& v : shp->texcoord) obj->texcoord.push_back({v.x, v.y});
        for (auto& v : shp->color) obj->color.push_back({v.x, v.y, v.z, v.w});
        for (auto& v : shp->radius) obj->radius.push_back(v);
        auto group = new obj_group();
        group->groupname = shp->name;
        group->matname = (shp->mat) ? shp->mat->name : "";
        group->subdivision_level = shp->subdivision_level;
        group->subdivision_catmullclark = shp->subdivision_catmullclark;
        for (auto point : shp->points) {
            group->elems.push_back(
                {(uint32_t)group->verts.size(), obj_element_type::point, 1});
            auto vert = obj_vertex{-1, -1, -1, -1, -1};
            if (!shp->pos.empty()) vert.pos = offset.pos + point;
            if (!shp->texcoord.empty()) vert.texcoord = offset.texcoord + point;
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
                for (auto vid : quad.xyz()) {
                    auto vert = obj_vertex{-1, -1, -1, -1, -1};
                    if (!shp->pos.empty()) vert.pos = offset.pos + vid;
                    if (!shp->texcoord.empty())
                        vert.texcoord = offset.texcoord + vid;
                    if (!shp->norm.empty()) vert.norm = offset.norm + vid;
                    if (!shp->color.empty()) vert.color = offset.color + vid;
                    if (!shp->radius.empty()) vert.radius = offset.radius + vid;
                    group->verts.push_back(vert);
                }
            } else {
                for (auto vid : quad) {
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
            group->elems.push_back(
                {(uint32_t)group->verts.size(), obj_element_type::bezier, 4});
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
        shapes[shp] = unique_ptr<obj_group>(group);
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
        omat->ke_txt = add_texture(env->ke_txt);
        oenv->name = env->name;
        oenv->matname = omat->name;
        oenv->frame = env->frame;
        obj->materials.push_back(omat);
        obj->environments.push_back(oenv);
    }

    // convert hierarchy
    if (obj->nodes.empty()) {
        // meshes
        for (auto shp : scn->shapes) {
            auto omsh = new obj_object();
            omsh->name = shp->name;
            omsh->groups.push_back(*shapes.at(shp));
            obj->objects.push_back(omsh);
        }

        // instances
        for (auto ist : scn->instances) {
            auto onde = new obj_node();
            onde->name = ist->name;
            onde->objname = (ist->shp) ? ist->shp->name : "<undefined>";
            onde->frame = ist->frame;
            obj->nodes.emplace_back(onde);
        }
    } else {
        // meshes
        auto meshes = map<vector<instance*>, obj_object*>();
        for (auto nde : scn->nodes) {
            if (nde->ists.empty()) continue;
            if (contains(meshes, nde->ists)) continue;
            auto omsh = new obj_object();
            omsh->name = nde->ists.front()->name;
            for (auto ist : nde->ists)
                omsh->groups.push_back(*shapes.at(ist->shp));
            meshes[nde->ists] = omsh;
            obj->objects.push_back(omsh);
        }

        // nodes
        for (auto nde : scn->nodes) {
            auto onde = new obj_node();
            onde->name = nde->name;
            if (nde->cam) onde->camname = nde->cam->name;
            if (!nde->ists.empty()) onde->objname = meshes.at(nde->ists)->name;
            if (nde->env) onde->envname = nde->env->name;
            onde->frame = to_mat(nde->frame);
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
    }

    return obj;
}

// Save an obj scene
void save_obj_scene(
    const string& filename, const scene* scn, const save_options& opts) {
    auto oscn = unique_ptr<obj_scene>(scene_to_obj(scn));
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
        txt->path =
            (startswith(gtxt->uri, "data:")) ? string("inlines") : gtxt->uri;
        if (!gtxt->data.datab.empty()) {
            txt->ldr = image4b(gtxt->data.width, gtxt->data.height);
            for (auto j = 0; j < gtxt->data.height; j++) {
                for (auto i = 0; i < gtxt->data.width; i++) {
                    auto v = gtxt->data.datab.data() +
                             (gtxt->data.width * j + i) * gtxt->data.ncomp;
                    switch (gtxt->data.ncomp) {
                        case 1:
                            txt->ldr.at(i, j) = {v[0], v[0], v[0], 255};
                            break;
                        case 2: txt->ldr.at(i, j) = {v[0], v[1], 0, 255}; break;
                        case 3:
                            txt->ldr.at(i, j) = {v[0], v[1], v[2], 255};
                            break;
                        case 4:
                            txt->ldr.at(i, j) = {v[0], v[1], v[2], v[3]};
                            break;
                        default: assert(false); break;
                    }
                }
            }
        } else if (!gtxt->data.dataf.empty()) {
            txt->hdr = image4f(gtxt->data.width, gtxt->data.height);
            for (auto j = 0; j < gtxt->data.height; j++) {
                for (auto i = 0; i < gtxt->data.width; i++) {
                    auto v = gtxt->data.dataf.data() +
                             (gtxt->data.width * j + i) * gtxt->data.ncomp;
                    switch (gtxt->data.ncomp) {
                        case 1:
                            txt->hdr.at(i, j) = {v[0], v[0], v[0], 1};
                            break;
                        case 2: txt->hdr.at(i, j) = {v[0], v[1], 0, 1}; break;
                        case 3:
                            txt->hdr.at(i, j) = {v[0], v[1], v[2], 1};
                            break;
                        case 4:
                            txt->hdr.at(i, j) = {v[0], v[1], v[2], v[3]};
                            break;
                        default: assert(false); break;
                    }
                }
            }
        }
        scn->textures.push_back(txt);
    }

    // add a texture
    auto add_texture = [gltf, scn](glTFTextureInfo* ginfo, bool normal = false,
                           bool occlusion = false) {
        auto info = texture_info();
        if (!ginfo) return info;
        auto gtxt = gltf->get(ginfo->index);
        if (!gtxt || !gtxt->source) return info;
        auto txt = scn->textures.at((int)gtxt->source);
        if (!txt) return info;
        info.txt = scn->textures.at((int)gtxt->source);
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
        mat->ke_txt = add_texture(gmat->emissiveTexture);
        if (gmat->pbrMetallicRoughness) {
            mat->type = material_type::metallic_roughness;
            auto gmr = gmat->pbrMetallicRoughness;
            mat->kd = {gmr->baseColorFactor[0], gmr->baseColorFactor[1],
                gmr->baseColorFactor[2]};
            mat->op = gmr->baseColorFactor[3];
            mat->ks = {
                gmr->metallicFactor, gmr->metallicFactor, gmr->metallicFactor};
            mat->rs = gmr->roughnessFactor;
            mat->kd_txt = add_texture(gmr->baseColorTexture);
            mat->ks_txt = add_texture(gmr->metallicRoughnessTexture);
        }
        if (gmat->pbrSpecularGlossiness) {
            mat->type = material_type::specular_glossiness;
            auto gsg = gmat->pbrSpecularGlossiness;
            mat->kd = {gsg->diffuseFactor[0], gsg->diffuseFactor[1],
                gsg->diffuseFactor[2]};
            mat->op = gsg->diffuseFactor[3];
            mat->ks = gsg->specularFactor;
            mat->rs = gsg->glossinessFactor;
            mat->kd_txt = add_texture(gsg->diffuseTexture);
            mat->ks_txt = add_texture(gsg->specularGlossinessTexture);
        }
        mat->norm_txt = add_texture(gmat->normalTexture, true, false);
        mat->occ_txt = add_texture(gmat->occlusionTexture, false, true);
        mat->double_sided = gmat->doubleSided;
        scn->materials.push_back(mat);
    }

    // convert meshes
    auto meshes = vector<vector<shape*>>();
    for (auto gmesh : gltf->meshes) {
        meshes.push_back({});
        // primitives
        for (auto gprim : gmesh->primitives) {
            auto shp = new shape();
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
            scn->shapes.push_back(shp);
            meshes.back().push_back(shp);
        }
    }

    // convert cameras
    auto cameras = vector<camera>();
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
            for (auto shp : meshes[(int)gnde->mesh]) {
                auto ist = new instance();
                ist->name = gnde->name;
                ist->shp = shp;
                nde->ists.push_back(ist);
                scn->instances.push_back(ist);
            }
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
        nde->frame = to_frame(gnde->matrix);
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
        unordered_map<glTFAnimationSamplerInterpolation, keyframe_type>{
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
        auto anm = new animation();
        anm->name = ganm->name;
        auto sampler_map = unordered_map<vec2i, keyframe*>();
        for (auto gchannel : ganm->channels) {
            if (sampler_map.find({(int)gchannel->sampler,
                    (int)gchannel->target->path}) == sampler_map.end()) {
                auto gsampler = ganm->get(gchannel->sampler);
                auto kfr = new keyframe();
                auto input_view =
                    accessor_view(gltf, gltf->get(gsampler->input));
                kfr->times.resize(input_view.size());
                for (auto i = 0; i < input_view.size(); i++)
                    kfr->times[i] = input_view.get(i);
                kfr->type = keyframe_types.at(gsampler->interpolation);
                auto output_view =
                    accessor_view(gltf, gltf->get(gsampler->output));
                switch (gchannel->target->path) {
                    case glTFAnimationChannelTargetPath::Translation: {
                        kfr->translation.reserve(output_view.size());
                        for (auto i = 0; i < output_view.size(); i++)
                            kfr->translation.push_back(output_view.getv3f(i));
                    } break;
                    case glTFAnimationChannelTargetPath::Rotation: {
                        kfr->rotation.reserve(output_view.size());
                        for (auto i = 0; i < output_view.size(); i++)
                            kfr->rotation.push_back(
                                (quat4f)output_view.getv4f(i));
                    } break;
                    case glTFAnimationChannelTargetPath::Scale: {
                        kfr->scaling.reserve(output_view.size());
                        for (auto i = 0; i < output_view.size(); i++)
                            kfr->scaling.push_back(output_view.getv3f(i));
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
                            kfr->weights.resize(values.size() / ncomp);
                            for (auto i = 0; i < kfr->weights.size(); i++) {
                                kfr->weights[i].resize(ncomp);
                                for (auto j = 0; j < ncomp; j++)
                                    kfr->weights[i][j] = values[i * ncomp + j];
                            }
                        }
                    } break;
                    default: {
                        throw runtime_error("should not have gotten here");
                    }
                }
                sampler_map[{
                    (int)gchannel->sampler, (int)gchannel->target->path}] = kfr;
                anm->keyframes.push_back(kfr);
            }
            anm->targets.push_back({sampler_map.at({(int)gchannel->sampler,
                                        (int)gchannel->target->path}),
                scn->nodes[(int)gchannel->target->node]});
        }
        scn->animations.push_back(anm);
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
scene* load_gltf_scene(const string& filename, const load_options& opts) {
    auto gscn = unique_ptr<glTF>(
        load_gltf(filename, true, opts.load_textures, opts.skip_missing));
    auto scn = unique_ptr<scene>(gltf_to_scene(gscn.get(), opts));
    if (!scn) {
        throw runtime_error("could not convert gltf scene");
        return nullptr;
    }
    return scn.release();
}

// Unflattnes gltf
glTF* scene_to_gltf(
    const scene* scn, const string& buffer_uri, bool separate_buffers) {
    auto gltf = unique_ptr<glTF>(new glTF());

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
        if (txt->hdr) {
            gimg->data.width = txt->hdr.width();
            gimg->data.height = txt->hdr.height();
            gimg->data.ncomp = 4;
            gimg->data.dataf.assign((float*)data(txt->hdr),
                (float*)data(txt->hdr) +
                    txt->hdr.width() * txt->hdr.height() * 4);
        }
        if (txt->ldr) {
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
    auto add_texture = [&gltf, &index, scn](const texture_info& info,
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
        gmat->emissiveTexture = add_texture(mat->ke_txt);
        switch (mat->type) {
            case material_type::specular_roughness: {
                gmat->pbrSpecularGlossiness =
                    new glTFMaterialPbrSpecularGlossiness();
                auto gsg = gmat->pbrSpecularGlossiness;
                gsg->diffuseFactor = {
                    mat->kd[0], mat->kd[1], mat->kd[2], mat->op};
                gsg->specularFactor = mat->ks;
                gsg->glossinessFactor = 1 - mat->rs;
                gsg->diffuseTexture = add_texture(mat->kd_txt);
                gsg->specularGlossinessTexture = add_texture(mat->ks_txt);
            } break;
            case material_type::metallic_roughness: {
                gmat->pbrMetallicRoughness =
                    new glTFMaterialPbrMetallicRoughness();
                auto gmr = gmat->pbrMetallicRoughness;
                gmr->baseColorFactor = {
                    mat->kd.x, mat->kd.y, mat->kd.z, mat->op};
                gmr->metallicFactor = mat->ks.x;
                gmr->roughnessFactor = mat->rs;
                gmr->baseColorTexture = add_texture(mat->kd_txt);
                gmr->metallicRoughnessTexture = add_texture(mat->ks_txt);
            } break;
            case material_type::specular_glossiness: {
                gmat->pbrSpecularGlossiness =
                    new glTFMaterialPbrSpecularGlossiness();
                auto gsg = gmat->pbrSpecularGlossiness;
                gsg->diffuseFactor = {
                    mat->kd[0], mat->kd[1], mat->kd[2], mat->op};
                gsg->specularFactor = mat->ks;
                gsg->glossinessFactor = mat->rs;
                gsg->diffuseTexture = add_texture(mat->kd_txt);
                gsg->specularGlossinessTexture = add_texture(mat->ks_txt);
            } break;
        }
        gmat->normalTexture = (glTFMaterialNormalTextureInfo*)add_texture(
            mat->norm_txt, true, false);
        gmat->occlusionTexture = (glTFMaterialOcclusionTextureInfo*)add_texture(
            mat->occ_txt, false, true);
        gmat->doubleSided = mat->double_sided;
        gltf->materials.push_back(gmat);
    }

    // add buffer
    auto add_buffer = [&gltf](const string& buffer_uri) {
        auto gbuffer = new glTFBuffer();
        gltf->buffers.push_back(gbuffer);
        gbuffer->uri = buffer_uri;
        return gbuffer;
    };

    // init buffers
    auto gbuffer_global = add_buffer(buffer_uri);

    // add an optional buffer
    auto add_opt_buffer = [&gbuffer_global, buffer_uri, &add_buffer,
                              separate_buffers](const string& uri) {
        if (separate_buffers && uri != "") {
            return add_buffer(uri);
        } else {
            if (!gbuffer_global) gbuffer_global = add_buffer(buffer_uri);
            return gbuffer_global;
        }
    };

    // attribute handling
    auto add_accessor = [&gltf, &index](glTFBuffer* gbuffer, const string& name,
                            glTFAccessorType type,
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
    auto shapes = unordered_map<shape*, unique_ptr<glTFMeshPrimitive>>();
    for (auto shp : scn->shapes) {
        auto gbuffer = add_opt_buffer(shp->path);
        auto gprim = new glTFMeshPrimitive();
        gprim->material = glTFid<glTFMaterial>(index(scn->materials, shp->mat));
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
            throw runtime_error("face varying not supported in glTF");
        } else {
            throw runtime_error("empty mesh");
        }
        shapes[shp] = unique_ptr<glTFMeshPrimitive>(gprim);
    }

    // hierarchy
    if (scn->nodes.empty()) {
        // meshes
        for (auto shp : scn->shapes) {
            auto gmsh = new glTFMesh();
            gmsh->name = shp->name;
            gmsh->primitives.push_back(new glTFMeshPrimitive(*shapes.at(shp)));
            gltf->meshes.push_back(gmsh);
        }

        // instances
        for (auto ist : scn->instances) {
            auto gnode = new glTFNode();
            gnode->name = ist->name;
            gnode->mesh = glTFid<glTFMesh>(index(scn->shapes, ist->shp));
            gnode->matrix = to_mat(ist->frame);
            gltf->nodes.push_back(gnode);
        }

        // cameras
        for (auto cam : scn->cameras) {
            auto gnode = new glTFNode();
            gnode->name = cam->name;
            gnode->camera = glTFid<glTFCamera>(index(scn->cameras, cam));
            gnode->matrix = to_mat(cam->frame);
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
        // meshes
        auto meshes = map<vector<instance*>, int>();
        for (auto nde : scn->nodes) {
            if (nde->ists.empty()) continue;
            if (contains(meshes, nde->ists)) continue;
            auto gmsh = new glTFMesh();
            gmsh->name = nde->ists.front()->name;
            for (auto ist : nde->ists)
                gmsh->primitives.push_back(
                    new glTFMeshPrimitive(*shapes.at(ist->shp)));
            meshes[nde->ists] = (int)gltf->meshes.size();
            gltf->meshes.push_back(gmsh);
        }

        // nodes
        for (auto nde : scn->nodes) {
            auto gnode = new glTFNode();
            gnode->name = nde->name;
            if (nde->cam) {
                gnode->camera =
                    glTFid<glTFCamera>(index(scn->cameras, nde->cam));
            }
            if (!nde->ists.empty()) {
                gnode->mesh = glTFid<glTFMesh>(meshes.at(nde->ists));
            }
            gnode->matrix = to_mat(nde->frame);
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
        auto is_root = vector<bool>(gltf->nodes.size(), true);
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
    for (auto anm : scn->animations) {
        auto ganm = new glTFAnimation();
        ganm->name = anm->name;
        auto gbuffer = add_opt_buffer(anm->path);
        auto count = 0;
        auto paths = unordered_map<keyframe*, glTFAnimationChannelTargetPath>();
        for (auto kfr : anm->keyframes) {
            auto aid = ganm->name + "_" + std::to_string(count++);
            auto gsmp = new glTFAnimationSampler();
            gsmp->input =
                add_accessor(gbuffer, aid + "_time", glTFAccessorType::Scalar,
                    glTFAccessorComponentType::Float, (int)kfr->times.size(),
                    sizeof(float), kfr->times.data(), false);
            auto path = glTFAnimationChannelTargetPath::NotSet;
            if (!kfr->translation.empty()) {
                gsmp->output = add_accessor(gbuffer, aid + "_translation",
                    glTFAccessorType::Vec3, glTFAccessorComponentType::Float,
                    (int)kfr->translation.size(), sizeof(vec3f),
                    kfr->translation.data(), false);
                path = glTFAnimationChannelTargetPath::Translation;
            } else if (!kfr->rotation.empty()) {
                gsmp->output = add_accessor(gbuffer, aid + "_rotation",
                    glTFAccessorType::Vec4, glTFAccessorComponentType::Float,
                    (int)kfr->rotation.size(), sizeof(vec4f),
                    kfr->rotation.data(), false);
                path = glTFAnimationChannelTargetPath::Rotation;
            } else if (!kfr->scaling.empty()) {
                gsmp->output = add_accessor(gbuffer, aid + "_scale",
                    glTFAccessorType::Vec3, glTFAccessorComponentType::Float,
                    (int)kfr->scaling.size(), sizeof(vec3f),
                    kfr->scaling.data(), false);
                path = glTFAnimationChannelTargetPath::Scale;
            } else if (!kfr->weights.empty()) {
                auto values = std::vector<float>();
                values.reserve(kfr->weights.size() * kfr->weights[0].size());
                for (auto i = 0; i < kfr->weights.size(); i++) {
                    values.insert(values.end(), kfr->weights[i].begin(),
                        kfr->weights[i].end());
                }
                gsmp->output = add_accessor(gbuffer, aid + "_weights",
                    glTFAccessorType::Scalar, glTFAccessorComponentType::Float,
                    (int)values.size(), sizeof(float), values.data(), false);
                path = glTFAnimationChannelTargetPath::Weights;
            } else {
                throw runtime_error("should not have gotten here");
            }
            gsmp->interpolation = interpolation_map.at(kfr->type);
            ganm->samplers.push_back(gsmp);
            paths[kfr] = path;
        }
        for (auto target : anm->targets) {
            auto kfr = target.first;
            auto node = target.second;
            auto gchan = new glTFAnimationChannel();
            gchan->sampler =
                glTFid<glTFAnimationSampler>{index(anm->keyframes, kfr)};
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
    const string& filename, const scene* scn, const save_options& opts) {
    auto buffer_uri = path_basename(filename) + ".bin";
    auto gscn = unique_ptr<glTF>(
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
        auto shp = new shape();
        shp->name = "shape" + to_string(sid++);
        for (auto spth : sshp->paths) {
            auto o = (int)shp->pos.size();
            for (auto p : spth->pos) shp->pos += {p.x, p.y, 0};
            for (auto vid = 1; vid < spth->pos.size(); vid += 3)
                shp->beziers +=
                    {o + vid - 1, o + vid + 0, o + vid + 1, o + vid + 2};
        }
        scn->shapes += shp;
    }
    auto miny = flt_max, maxy = -flt_max;
    for (auto shp : scn->shapes) {
        for (auto& p : shp->pos) {
            miny = min(miny, p.y);
            maxy = max(maxy, p.y);
        }
    }
    auto mdly = (maxy + miny) / 2;
    for (auto shp : scn->shapes) {
        for (auto& p : shp->pos) p.y = -(p.y - mdly) + mdly;
    }
    return scn;
}

// Load an svg scene
scene* load_svg_scene(const string& filename, const load_options& opts) {
    auto sscn = unique_ptr<svg_scene>(load_svg(filename));
    auto scn = unique_ptr<scene>(svg_to_scene(sscn.get(), opts));
    return scn.release();
}

#endif

// Load a scene
scene* load_scene(const string& filename, const load_options& opts) {
    auto ext = path_extension(filename);
    if (ext == ".obj" || ext == ".OBJ") return load_obj_scene(filename, opts);
#if YGL_GLTF
    if (ext == ".gltf" || ext == ".GLTF")
        return load_gltf_scene(filename, opts);
#endif
#if YGL_SVG
    if (ext == ".svg" || ext == ".SVG") return load_svg_scene(filename, opts);
#endif
    throw runtime_error("unsupported extension " + ext);
    return nullptr;
}

// Save a scene
void save_scene(
    const string& filename, const scene* scn, const save_options& opts) {
    auto ext = path_extension(filename);
    if (ext == ".obj" || ext == ".OBJ")
        return save_obj_scene(filename, scn, opts);
#if YGL_GLTF
    if (ext == ".gltf" || ext == ".GLTF")
        return save_gltf_scene(filename, scn, opts);
#endif
    throw runtime_error("unsupported extension " + ext);
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
float trace_next1f(trace_pixel& pxl) {
    switch (pxl.rtype) {
        case trace_rng_type::uniform: {
            return clamp(next_rand1f(pxl.rng), 0.0f, 1 - flt_eps);
        } break;
        case trace_rng_type::stratified: {
            auto p = hash_uint64_32((uint64_t)pxl.i | (uint64_t)pxl.j << 16 |
                                    (uint64_t)pxl.dimension << 32);
            auto s = hash_permute(pxl.sample, pxl.nsamples, p);
            pxl.dimension += 1;
            return clamp(
                (s + next_rand1f(pxl.rng)) / pxl.nsamples, 0.0f, 1 - flt_eps);
        } break;
        default: {
            assert(false);
            return 0;
        }
    }
}

// Generates a 2-dimensional sample.
vec2f trace_next2f(trace_pixel& pxl) {
    switch (pxl.rtype) {
        case trace_rng_type::uniform: {
            return {next_rand1f(pxl.rng), next_rand1f(pxl.rng)};
        } break;
        case trace_rng_type::stratified: {
            auto p = hash_uint64_32((uint64_t)pxl.i | (uint64_t)pxl.j << 16 |
                                    (uint64_t)pxl.dimension << 32);
            auto s = hash_permute(pxl.sample, pxl.nsamples, p);
            auto nsamples2 = (int)round(sqrt(pxl.nsamples));
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

// Creates a 1-dimensional sample in [0,num-1]
inline int trace_next1i(trace_pixel& pxl, int num) {
    return clamp(int(trace_next1f(pxl) * num), 0, num - 1);
}

// Brdf type
enum struct trace_brdf_type {
    none = 0,
    microfacet = 1,
    kajiya_kay = 2,
    point = 3
};

// Brdf
struct trace_brdf {
    trace_brdf_type type = trace_brdf_type::none;  // type
    vec3f kd = {0, 0, 0};                          // diffuse
    vec3f ks = {0, 0, 0};                          // specular
    float rs = 0;                                  // specular roughness
    vec3f kt = {0, 0, 0};                          // transmission (thin glass)
    operator bool() const { return type != trace_brdf_type::none; }
    vec3f rho() const { return kd + ks + kt; }
};

// Emission type
enum struct trace_emission_type {
    none = 0,
    diffuse = 1,
    point = 2,
    line = 3,
    env = 4,
};

// Emission
struct trace_emission {
    trace_emission_type type = trace_emission_type::none;
    vec3f ke = zero3f;
    operator bool() const { return type != trace_emission_type::none; }
};

// Surface point with geometry and material data. Supports point on envmap too.
// This is the key data manipulated in the path tracer.
struct trace_point {
    const instance* ist = nullptr;     // instance
    const environment* env = nullptr;  // environment
    frame3f frame = identity_frame3f;  // local frame
    vec3f wo = zero3f;                 // outgoing direction
    trace_emission em = {};            // emission
    trace_brdf fr = {};                // brdf
};

// Evaluates emission.
inline vec3f trace_eval_emission(const trace_point& pt) {
    auto& em = pt.em;
    auto& wo = pt.wo;
    auto& wn = pt.frame.z;

    if (!em) return zero3f;
    auto ke = zero3f;
    switch (em.type) {
        case trace_emission_type::diffuse:
            ke += (dot(wn, wo) > 0) ? em.ke : zero3f;
            break;
        case trace_emission_type::point: ke += em.ke; break;
        case trace_emission_type::line: ke += em.ke; break;
        case trace_emission_type::env: ke += em.ke; break;
        default: assert(false); break;
    }
    return ke;
}

// Evaluates the BRDF scaled by the cosine of the incoming direction.
//
// Implementation notes:
// - ggx from [Heitz 2014] and [Walter 2007] and [Lagarde 2014]
// "Understanding the Masking-Shadowing Function in Microfacet-Based BRDFs"
// http://jcgt.org/published/0003/02/03/
// - "Microfacet Models for Refraction through Rough Surfaces" EGSR 07
// https://www.cs.cornell.edu/~srm/publications/EGSR07-btdf.pdf
// - uses Kajiya-Kay for hair
// - uses a hack for points
inline vec3f trace_eval_brdfcos(
    const trace_point& pt, const vec3f& wi, bool delta = false) {
    // grab variables
    auto& fr = pt.fr;
    auto& wn = pt.frame.z;
    auto& wo = pt.wo;

    // exit if not needed
    if (!fr) return zero3f;

    // accumulate brdfcos for each lobe
    auto brdfcos = zero3f;
    switch (fr.type) {
        // reflection terms
        case trace_brdf_type::microfacet: {
            // compute wh
            auto wh = normalize(wo + wi);

            // compute dot products
            auto ndo = dot(wn, wo), ndi = dot(wn, wi),
                 ndh = clamp(dot(wh, wn), (float)-1, (float)1);

            // diffuse term
            if (fr.kd != zero3f && ndi > 0 && ndo > 0) {
                brdfcos += fr.kd * ndi / pif;
            }

            // specular term (GGX)
            if (fr.ks != zero3f && ndi > 0 && ndo > 0 && ndh > 0 && fr.rs) {
                // microfacet term
                auto dg = eval_ggx(fr.rs, ndh, ndi, ndo);

                // handle fresnel
                auto odh = clamp(dot(wo, wh), 0.0f, 1.0f);
                auto ks = fresnel_schlick(fr.ks, odh, fr.rs);

                // sum up
                brdfcos += ks * ndi * dg / (4 * ndi * ndo);
            }

            // specular term (mirror)
            if (fr.ks != zero3f && ndi > 0 && ndo > 0 && !fr.rs && delta) {
                // handle fresnel
                auto ks = fresnel_schlick(fr.ks, ndo, fr.rs);

                // sum up
                brdfcos += ks;
            }

            // transmission hack
            if (fr.kt != zero3f && wo == -wi) brdfcos += fr.kt;
        } break;
        // hair (Kajiya-Kay)
        case trace_brdf_type::kajiya_kay: {
            // compute wh
            auto wh = normalize(wo + wi);

            // compute dot products
            auto ndo = dot(wn, wo), ndi = dot(wn, wi),
                 ndh = clamp(dot(wh, wn), (float)0, (float)1);

            // take sines
            auto so = sqrt(clamp(1 - ndo * ndo, (float)0, (float)1)),
                 si = sqrt(clamp(1 - ndi * ndi, (float)0, (float)1)),
                 sh = sqrt(clamp(1 - ndh * ndh, (float)0, (float)1));

            // diffuse term (Kajiya-Kay)
            if (fr.kd != zero3f && si > 0 && so > 0) {
                brdfcos += fr.kd * si / pif;
            }

            // specular term (Kajiya-Kay)
            if (fr.ks != zero3f && si > 0 && so > 0 && sh > 0) {
                auto ns = 2 / (fr.rs * fr.rs) - 2;
                auto d = (ns + 2) * pow(sh, ns) / (2 + pif);
                brdfcos += fr.ks * si * d / (4.0f * si * so);
            }

            // transmission hack
            if (fr.kt != zero3f && wo == -wi) brdfcos += fr.kt;
        } break;
        // points
        case trace_brdf_type::point: {
            // diffuse term
            auto ido = dot(wo, wi);
            brdfcos += fr.kd * (2 * ido + 1) / (2 * pif);

            // transmission hack
            if (fr.kt != zero3f && wo == -wi) brdfcos += fr.kt;
        } break;
        default: assert(false); break;
    }

    // check
    assert(isfinite(brdfcos.x) && isfinite(brdfcos.y) && isfinite(brdfcos.z));

    // done
    return brdfcos;
}

// Compute the weight for sampling the BRDF
inline float trace_weight_brdfcos(
    const trace_point& pt, const vec3f& wi, bool delta = false) {
    // grab variables
    auto& fr = pt.fr;
    auto& wn = pt.frame.z;
    auto& wo = pt.wo;

    // skip if no component
    if (!fr) return 0;

    // probability of each lobe
    auto kdw = max_element_value(fr.kd), ksw = max_element_value(fr.ks),
         ktw = max_element_value(fr.kt);
    auto kaw = kdw + ksw + ktw;
    kdw /= kaw;
    ksw /= kaw;
    ktw /= kaw;

    // accumulate the probability over all lobes
    auto pdf = 0.0f;
    // sample the lobe
    switch (fr.type) {
        // reflection term
        case trace_brdf_type::microfacet: {
            // compute wh
            auto wh = normalize(wi + wo);

            // compute dot products
            auto ndo = dot(wn, wo), ndi = dot(wn, wi), ndh = dot(wn, wh);

            // diffuse term (hemipherical cosine probability)
            if (kdw && ndo > 0 && ndi > 0) { pdf += kdw * ndi / pif; }

            // specular term (GGX)
            if (ksw && ndo > 0 && ndi > 0 && ndh > 0 && fr.rs) {
                // probability proportional to d adjusted by wh projection
                auto d = sample_ggx_pdf(fr.rs, ndh);
                auto hdo = dot(wo, wh);
                pdf += ksw * d / (4 * hdo);
            }

            // specular term (mirror)
            if (ksw && ndo > 0 && ndi > 0 && !fr.rs && delta) {
                // probability proportional to d adjusted by wh projection
                pdf += ksw;
            }

            // transmission hack
            if (ktw && wi == -wo) pdf += ktw;

            // check
            assert(isfinite(pdf));
        } break;
        // hair (Kajiya-Kay)
        case trace_brdf_type::kajiya_kay: {
            // diffuse and specular
            pdf += (kdw + ksw) * 4 * pif;
            // transmission hack
            if (wi == -wo) pdf += ktw;
        } break;
        // point
        case trace_brdf_type::point: {
            // diffuse and specular
            pdf += (kdw + ksw) * 4 * pif;
            // transmission hack
            if (wi == -wo) pdf += ktw;
        } break;
        default: assert(false); break;
    }

    // check for missed pdf
    if (!pdf) return 0;

    // check
    assert(isfinite(pdf));

    // done
    return 1 / pdf;
}

// reflected vector
inline vec3f reflect(const vec3f& w, const vec3f& n) {
    return -w + 2 * dot(n, w) * n;
}

// refracted vector
inline vec3f refract(const vec3f& w, const vec3f& n, float eta) {
    // auto k = 1.0 - eta * eta * (1.0 - dot(n, w) * dot(n, w));
    auto k = 1 - eta * eta * max(0.0f, 1 - dot(n, w) * dot(n, w));
    if (k < 0) return zero3f;  // tir
    return -w * eta + (eta * dot(n, w) - sqrt(k)) * n;
}

// Picks a direction based on the BRDF
inline tuple<vec3f, bool> trace_sample_brdfcos(
    const trace_point& pt, float rnl, const vec2f& rn) {
    // grab variables
    auto& fr = pt.fr;
    auto& wn = pt.frame.z;
    auto& fp = pt.frame;
    auto& wo = pt.wo;

    // skip if no component
    if (!fr) return {zero3f, false};

    // probability of each lobe
    auto kdw = max_element_value(fr.kd), ksw = max_element_value(fr.ks),
         ktw = max_element_value(fr.kt);
    auto kaw = kdw + ksw + ktw;
    kdw /= kaw;
    ksw /= kaw;
    ktw /= kaw;

    // sample selected lobe
    switch (fr.type) {
        // reflection term
        case trace_brdf_type::microfacet: {
            // compute cosine
            auto ndo = dot(wn, wo);

            // check to make sure we are above the surface
            if (ndo <= 0) return {zero3f, false};

            // sample according to diffuse
            if (rnl < kdw) {
                // sample wi with hemispherical cosine distribution
                auto rz = sqrtf(rn.y), rr = sqrtf(1 - rz * rz),
                     rphi = 2 * pif * rn.x;
                // set to wi
                auto wi_local = vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
                return {transform_direction(fp, wi_local), false};
            }
            // sample according to specular GGX
            else if (rnl < kdw + ksw && fr.rs) {
                // sample wh with ggx distribution
                auto wh_local = sample_ggx(fr.rs, rn);
                auto wh = transform_direction(fp, wh_local);
                // compute wi
                return {normalize(wh * 2.0f * dot(wo, wh) - wo), false};
            }
            // sample according to specular mirror
            else if (rnl < kdw + ksw && !fr.rs) {
                // compute wi
                return {normalize(wn * 2.0f * dot(wo, wn) - wo), true};
            }
            // transmission hack
            else if (rnl < kdw + ksw + ktw) {
                // continue ray direction
                return {-wo, true};
            } else
                assert(false);
        } break;
        // hair (Kajiya-Kay)
        case trace_brdf_type::kajiya_kay: {
            // diffuse and specular
            if (rnl < kdw + ksw) {
                // sample wi with uniform spherical distribution
                auto rz = 2 * rn.y - 1, rr = sqrtf(1 - rz * rz),
                     rphi = 2 * pif * rn.x;
                auto wi_local = vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
                return {transform_direction(fp, wi_local), false};
            }
            // transmission hack
            else if (rnl < kdw + ksw + ktw) {
                // continue ray direction
                return {-wo, true};
            } else
                assert(false);
        } break;
        // diffuse term point
        case trace_brdf_type::point: {
            // diffuse and specular
            if (rnl < kdw + ksw) {
                // sample wi with uniform spherical distribution
                auto rz = 2 * rn.y - 1, rr = sqrtf(1 - rz * rz),
                     rphi = 2 * pif * rn.x;
                auto wi_local = vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
                return {transform_direction(fp, wi_local), false};
            }
            // transmission hack
            else if (rnl < kdw + ksw + ktw) {
                // continue ray direction
                return {-wo, true};
            } else
                assert(false);
        } break;
        default: assert(false); break;
    }

    // done
    return {zero3f, false};
}

// Create a point for an environment map. Resolves material with textures.
inline trace_point trace_eval_point(const environment* env, const vec3f& wo) {
    // set shape data
    auto pt = trace_point();

    // env
    pt.env = env;

    // direction
    pt.wo = wo;

    // maerial
    auto ke = env->ke;
    if (env->ke_txt) {
        auto w = transform_direction_inverse(env->frame, -wo);
        auto theta = acos(clamp(w.y, (float)-1, (float)1));
        auto phi = atan2(w.z, w.x);
        auto texcoord = vec2f{0.5f + phi / (2 * pif), theta / pif};
        ke *= eval_texture(env->ke_txt, texcoord).xyz();
    }

    // create emission lobe
    if (ke != zero3f) {
        pt.em.type = trace_emission_type::env;
        pt.em.ke = ke;
    }

    // done
    return pt;
}

// Create a point for a shape. Resolves geometry and material with textures.
inline trace_point trace_eval_point(
    const instance* ist, int eid, const vec4f& euv, const vec3f& wo) {
    // default material
    static auto def_material = (material*)nullptr;
    if (!def_material) {
        def_material = new material();
        def_material->kd = {0.2f, 0.2f, 0.2f};
        def_material->rs = 1;
    }

    // set shape data
    auto pt = trace_point();

    // instance
    pt.ist = ist;

    // direction
    pt.wo = wo;

    // shortcuts
    auto shp = ist->shp;
    auto mat = ist->shp->mat;
    if (!mat) mat = def_material;

    // compute points and weights
    auto pos = eval_pos(ist->shp, eid, euv);
    auto norm = eval_norm(ist->shp, eid, euv);
    auto texcoord = eval_texcoord(ist->shp, eid, euv);
    auto color = eval_color(ist->shp, eid, euv);

    // handle normal map
    if (mat->norm_txt) {
        auto tangsp = eval_tangsp(ist->shp, eid, euv);
        auto txt = eval_texture(mat->norm_txt, texcoord, false).xyz() * 2.0f -
                   vec3f{1, 1, 1};
        auto ntxt = normalize(vec3f{txt.x, -txt.y, txt.z});
        auto frame =
            make_frame_fromzx({0, 0, 0}, norm, {tangsp.x, tangsp.y, tangsp.z});
        frame.y *= tangsp.w;
        norm = transform_direction(frame, ntxt);
    }

    // correct for double sided
    if (mat->double_sided && dot(norm, wo) < 0) norm = -norm;

    // creating frame
    pt.frame = make_frame_fromz(transform_point(ist->frame, pos),
        transform_direction(ist->frame, norm));

    // handle color
    auto kx_scale = vec4f{1, 1, 1, 1};
    if (!shp->color.empty()) kx_scale *= color;

    // handle occlusion
    if (mat->occ_txt)
        kx_scale.xyz() *= eval_texture(mat->occ_txt, texcoord).xyz();

    // sample emission
    auto ke = mat->ke * kx_scale.xyz();

    // sample reflectance
    auto kd = zero4f, ks = zero4f, kt = zero4f;
    switch (mat->type) {
        case material_type::specular_roughness: {
            kd = vec4f{mat->kd, mat->op} * kx_scale *
                 eval_texture(mat->kd_txt, texcoord);
            ks = vec4f{mat->ks, mat->rs} * vec4f{kx_scale.xyz(), 1} *
                 eval_texture(mat->ks_txt, texcoord);
            kt = vec4f{mat->kt, mat->rs} * vec4f{kx_scale.xyz(), 1} *
                 eval_texture(mat->kt_txt, texcoord);
        } break;
        case material_type::metallic_roughness: {
            auto kb = vec4f{mat->kd, mat->op} * kx_scale *
                      eval_texture(mat->kd_txt, texcoord);
            auto km = vec2f{mat->ks.x, mat->rs};
            if (mat->ks_txt) {
                auto ks_txt = eval_texture(mat->ks_txt, texcoord);
                km.x *= ks_txt.y;
                km.y *= ks_txt.z;
            }
            kd = vec4f{kb.xyz() * (1 - km.x), kb.w};
            ks =
                vec4f{kb.xyz() * km.x + vec3f{0.04f, 0.04f, 0.04f} * (1 - km.x),
                    km.y};
        } break;
        case material_type::specular_glossiness: {
            kd = vec4f{mat->kd, mat->op} * kx_scale *
                 eval_texture(mat->kd_txt, texcoord);
            ks = vec4f{mat->ks, mat->rs} * vec4f{kx_scale.xyz(), 1} *
                 eval_texture(mat->ks_txt, texcoord);
            ks.w = 1 - ks.w;  // glossiness -> roughness
        } break;
    }

    // set up final values
    pt.em.ke = ke * kd.w;
    pt.fr.kd = kd.xyz() * kd.w;
    pt.fr.ks =
        (ks.xyz() != zero3f && ks.w < 0.9999f) ? ks.xyz() * kd.w : zero3f;
    pt.fr.rs = (ks.xyz() != zero3f && ks.w < 0.9999f) ? ks.w * ks.w : 0;
    pt.fr.kt = {1 - kd.w, 1 - kd.w, 1 - kd.w};
    if (kt.xyz() != zero3f) pt.fr.kt *= kt.xyz();

    // setup brdf and emission
    if (!shp->points.empty()) {
        if (kd.xyz() != zero3f || ks.xyz() != zero3f || kt.xyz() != zero3f)
            pt.fr.type = trace_brdf_type::point;
        if (ke != zero3f) pt.em.type = trace_emission_type::point;
    } else if (!shp->lines.empty()) {
        if (kd.xyz() != zero3f || ks.xyz() != zero3f || kt.xyz() != zero3f)
            pt.fr.type = trace_brdf_type::kajiya_kay;
        if (ke != zero3f) pt.em.type = trace_emission_type::line;
    } else if (!shp->triangles.empty()) {
        if (kd.xyz() != zero3f || ks.xyz() != zero3f || kt.xyz() != zero3f)
            pt.fr.type = trace_brdf_type::microfacet;
        if (ke != zero3f) pt.em.type = trace_emission_type::diffuse;
    }

    // done
    return pt;
}

// Sample weight for a light point.
inline float trace_weight_light(const trace_lights& lights,
                                const trace_point& lpt, const trace_point& pt) {
    if (!lpt.em) return 0;
    // support only one lobe for now
    switch (lpt.em.type) {
        case trace_emission_type::env: {
            return 4 * pif;
        } break;
        case trace_emission_type::point: {
            auto d = length(lpt.frame.o - pt.frame.o);
            return lights.shape_areas.at(lpt.ist->shp) / (d * d);
        } break;
        case trace_emission_type::line: {
            assert(false);
            return 0;
        } break;
        case trace_emission_type::diffuse: {
            auto d = length(lpt.frame.o - pt.frame.o);
            return lights.shape_areas.at(lpt.ist->shp) *
                   abs(dot(lpt.frame.z, lpt.wo)) / (d * d);
        } break;
        default: {
            assert(false);
            return 0;
        } break;
    }
}

// Picks a point on a light.
inline trace_point trace_sample_light(const trace_lights& lights,
    const trace_light& lgt, const trace_point& pt, float rne, const vec2f& rn) {
    if (lgt.ist) {
        auto eid = 0;
        auto euv = zero4f;
        if (!lgt.ist->shp->points.empty()) {
            eid = sample_points(lights.shape_cdfs.at(lgt.ist->shp), rne);
            euv = {1, 0, 0, 0};
        } else if (!lgt.ist->shp->lines.empty()) {
            std::tie(eid, (vec2f&)euv) = sample_lines(lights.shape_cdfs.at(lgt.ist->shp), rne, rn.x);
        } else if (!lgt.ist->shp->triangles.empty()) {
            std::tie(eid, (vec3f&)euv) =
                sample_triangles(lights.shape_cdfs.at(lgt.ist->shp), rne, rn);
        } else if (!lgt.ist->shp->quads.empty()) {
            std::tie(eid, (vec4f&)euv) = sample_quads(lights.shape_cdfs.at(lgt.ist->shp), rne, rn);
        } else {
            assert(false);
        }
        auto lpt = trace_eval_point(lgt.ist, eid, euv, zero3f);
        lpt.wo = normalize(pt.frame.o - lpt.frame.o);
        return lpt;
    } else if (lgt.env) {
        auto z = -1 + 2 * rn.y;
        auto rr = sqrt(clamp(1 - z * z, (float)0, (float)1));
        auto phi = 2 * pif * rn.x;
        auto wo = vec3f{cos(phi) * rr, z, sin(phi) * rr};
        auto lpt = trace_eval_point(lgt.env, wo);
        return lpt;
    } else {
        assert(false);
        return {};
    }
}

// Intersects a ray with the scn and return the point (or env
// point).
inline trace_point trace_intersect_scene(
    const scene* scn, const bvh_tree* bvh, const ray3f& ray) {
    auto iid = 0, eid = 0;
    auto euv = zero4f;
    auto ray_t = 0.0f;
    if (intersect_bvh(bvh, ray, false, ray_t, iid, eid, euv)) {
        return trace_eval_point(scn->instances[iid], eid, euv, -ray.d);
    } else if (!scn->environments.empty()) {
        return trace_eval_point(scn->environments[0], -ray.d);
    } else {
        return {};
    }
}

// Test occlusion
inline vec3f trace_eval_transmission(const scene* scn, const bvh_tree* bvh,
    const trace_point& pt, const trace_point& lpt, const trace_params& params) {
    if (params.shadow_notransmission) {
        auto shadow_ray = (lpt.env) ? make_ray(pt.frame.o, -lpt.wo) :
                                      make_segment(pt.frame.o, lpt.frame.o);
        // auto shadow_ray = ray3f{pt.frame.o, -lpt.wo, 0.01f, flt_max};
        return (intersect_bvh(bvh, shadow_ray, true)) ? zero3f : vec3f{1, 1, 1};
    } else {
        auto cpt = pt;
        auto weight = vec3f{1, 1, 1};
        for (auto bounce = 0; bounce < params.max_depth; bounce++) {
            auto ray = (lpt.env) ? make_ray(cpt.frame.o, -lpt.wo) :
                                   make_segment(cpt.frame.o, lpt.frame.o);
            cpt = trace_intersect_scene(scn, bvh, ray);
            if (!cpt.ist) break;
            weight *= cpt.fr.kt;
            if (weight == zero3f) break;
        }
        return weight;
    }
}

// Mis weight
inline float trace_weight_mis(float w0, float w1) {
    if (!w0 || !w1) return 1;
    return (1 / w0) / (1 / w0 + 1 / w1);
}

// Recursive path tracing.
inline vec3f trace_path(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const ray3f& ray, trace_pixel& pxl,
    const trace_params& params, bool& hit) {
    // intersection
    auto pt = trace_intersect_scene(scn, bvh, ray);
    hit = pt.ist;

    // emission
    auto l = trace_eval_emission(pt);
    if (!pt.fr || lights.empty()) return l;

    // trace path
    auto weight = vec3f{1, 1, 1};
    auto emission = false;
    for (auto bounce = 0; bounce < params.max_depth; bounce++) {
        // emission
        if (emission) l += weight * trace_eval_emission(pt);

        // direct – light
        auto lgt = lights.lights[trace_next1i(pxl, (int)lights.size())];
        auto lpt =
            trace_sample_light(lights, lgt, pt, trace_next1f(pxl), trace_next2f(pxl));
        auto lw = trace_weight_light(lights, lpt, pt) * (float)lights.size();
        auto lke = trace_eval_emission(lpt);
        auto lbc = trace_eval_brdfcos(pt, -lpt.wo);
        auto lld = lke * lbc * lw;
        if (lld != zero3f) {
            l += weight * lld *
                 trace_eval_transmission(scn, bvh, pt, lpt, params) *
                 trace_weight_mis(lw, trace_weight_brdfcos(pt, -lpt.wo));
        }

        // direct – brdf
        auto bwi = zero3f;
        auto bdelta = false;
        std::tie(bwi, bdelta) =
            trace_sample_brdfcos(pt, trace_next1f(pxl), trace_next2f(pxl));
        auto bpt = trace_intersect_scene(scn, bvh, make_ray(pt.frame.o, bwi));
        auto bw = trace_weight_brdfcos(pt, -bpt.wo, bdelta);
        auto bke = trace_eval_emission(bpt);
        auto bbc = trace_eval_brdfcos(pt, -bpt.wo, bdelta);
        auto bld = bke * bbc * bw;
        if (bld != zero3f) {
            l += weight * bld *
                 trace_weight_mis(bw, trace_weight_light(lights, bpt, pt));
        }

        // skip recursion if path ends
        if (bounce == params.max_depth - 1) break;
        if (!bpt.fr) break;

        // continue path
        weight *=
            trace_eval_brdfcos(pt, -bpt.wo) * trace_weight_brdfcos(pt, -bpt.wo);
        if (weight == zero3f) break;

        // roussian roulette
        if (bounce > 2) {
            auto rrprob = 1.0f - min(max_element_value(pt.fr.rho()), 0.95f);
            if (trace_next1f(pxl) < rrprob) break;
            weight *= 1 / (1 - rrprob);
        }

        // continue path
        pt = bpt;
        emission = false;
    }

    return l;
}

// Recursive path tracing.
inline vec3f trace_path_nomis(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const ray3f& ray, trace_pixel& pxl,
    const trace_params& params, bool& hit) {
    // intersection
    auto pt = trace_intersect_scene(scn, bvh, ray);
    hit = pt.ist;

    // emission
    auto l = trace_eval_emission(pt);
    if (!pt.fr) return l;

    // trace path
    auto weight = vec3f{1, 1, 1};
    auto emission = false;
    for (auto bounce = 0; bounce < params.max_depth; bounce++) {
        // emission
        if (emission) l += weight * trace_eval_emission(pt);

        // direct
        auto& lgt = lights.lights[trace_next1i(pxl, (int)lights.size())];
        auto lpt =
            trace_sample_light(lights, lgt, pt, trace_next1f(pxl), trace_next2f(pxl));
        auto ld = trace_eval_emission(lpt) * trace_eval_brdfcos(pt, -lpt.wo) *
                  trace_weight_light(lights, lpt, pt) * (float)lights.size();
        if (ld != zero3f) {
            l += weight * ld *
                 trace_eval_transmission(scn, bvh, pt, lpt, params);
        }

        // skip recursion if path ends
        if (bounce == params.max_depth - 1) break;

        // roussian roulette
        if (bounce > 2) {
            auto rrprob = 1.0f - min(max_element_value(pt.fr.rho()), 0.95f);
            if (trace_next1f(pxl) < rrprob) break;
            weight *= 1 / (1 - rrprob);
        }

        // continue path
        auto bwi = zero3f;
        auto bdelta = false;
        std::tie(bwi, bdelta) =
            trace_sample_brdfcos(pt, trace_next1f(pxl), trace_next2f(pxl));
        weight *= trace_eval_brdfcos(pt, bwi, bdelta) *
                  trace_weight_brdfcos(pt, bwi, bdelta);
        if (weight == zero3f) break;

        auto bpt = trace_intersect_scene(scn, bvh, make_ray(pt.frame.o, bwi));
        emission = false;
        if (!bpt.fr) break;

        // continue path
        pt = bpt;
    }

    return l;
}

// Recursive path tracing.
inline vec3f trace_path_hack(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const ray3f& ray, trace_pixel& pxl,
    const trace_params& params, bool& hit) {
    // intersection
    auto pt = trace_intersect_scene(scn, bvh, ray);
    hit = pt.ist;

    // emission
    auto l = trace_eval_emission(pt);
    if (!pt.fr || lights.empty()) return l;

    // trace path
    auto weight = vec3f{1, 1, 1};
    for (auto bounce = 0; bounce < params.max_depth; bounce++) {
        // direct
        auto& lgt = lights.lights[trace_next1i(pxl, (int)lights.size())];
        auto lpt =
            trace_sample_light(lights, lgt, pt, trace_next1f(pxl), trace_next2f(pxl));
        auto ld = trace_eval_emission(lpt) * trace_eval_brdfcos(pt, -lpt.wo) *
                  trace_weight_light(lights, lpt, pt) * (float)lights.size();
        if (ld != zero3f) {
            l += weight * ld *
                 trace_eval_transmission(scn, bvh, pt, lpt, params);
        }

        // skip recursion if path ends
        if (bounce == params.max_depth - 1) break;

        // roussian roulette
        if (bounce > 2) {
            auto rrprob = 1.0f - min(max_element_value(pt.fr.rho()), 0.95f);
            if (trace_next1f(pxl) < rrprob) break;
            weight *= 1 / (1 - rrprob);
        }

        // continue path
        auto bwi = zero3f;
        auto bdelta = false;
        std::tie(bwi, bdelta) =
            trace_sample_brdfcos(pt, trace_next1f(pxl), trace_next2f(pxl));
        weight *= trace_eval_brdfcos(pt, bwi, bdelta) *
                  trace_weight_brdfcos(pt, bwi, bdelta);
        if (weight == zero3f) break;

        auto bpt = trace_intersect_scene(scn, bvh, make_ray(pt.frame.o, bwi));
        if (!bpt.fr) break;

        // continue path
        pt = bpt;
    }

    return l;
}

// Direct illumination.
inline vec3f trace_direct(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const ray3f& ray, int bounce, trace_pixel& pxl,
    const trace_params& params, bool& hit) {
    // intersection
    auto pt = trace_intersect_scene(scn, bvh, ray);
    if (!bounce) hit = pt.ist;

    // emission
    auto l = trace_eval_emission(pt);
    if (!pt.fr || lights.empty()) return l;

    // ambient
    l += params.ambient * pt.fr.rho();

    // direct
    for (auto& lgt : lights.lights) {
        auto lpt =
            trace_sample_light(lights, lgt, pt, trace_next1f(pxl), trace_next2f(pxl));
        auto ld = trace_eval_emission(lpt) * trace_eval_brdfcos(pt, -lpt.wo) *
                  trace_weight_light(lights, lpt, pt);
        if (ld == zero3f) continue;
        l += ld * trace_eval_transmission(scn, bvh, pt, lpt, params);
    }

    // exit if needed
    if (bounce >= params.max_depth) return l;

    // reflection
    if (pt.fr.ks != zero3f && !pt.fr.rs) {
        auto wi = reflect(pt.wo, pt.frame.z);
        auto ray = make_ray(pt.frame.o, wi);
        l += pt.fr.ks *
             trace_direct(scn, bvh, lights, ray, bounce + 1, pxl, params, hit);
    }

    // opacity
    if (pt.fr.kt != zero3f) {
        auto ray = make_ray(pt.frame.o, -pt.wo);
        l += pt.fr.kt *
             trace_direct(scn, bvh, lights, ray, bounce + 1, pxl, params, hit);
    }

    // done
    return l;
}

// Direct illumination.
inline vec3f trace_direct(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const ray3f& ray, trace_pixel& pxl,
    const trace_params& params, bool& hit) {
    return trace_direct(scn, bvh, lights, ray, 0, pxl, params, hit);
}

// Eyelight for quick previewing.
inline vec3f trace_eyelight(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const ray3f& ray, int bounce, trace_pixel& pxl,
    const trace_params& params, bool& hit) {
    // intersection
    auto pt = trace_intersect_scene(scn, bvh, ray);
    if (!bounce) hit = pt.ist;

    // emission
    auto l = trace_eval_emission(pt);
    if (!pt.fr) return l;

    // brdf*light
    l += trace_eval_brdfcos(pt, pt.wo) * pif;

    // opacity
    if (bounce >= params.max_depth) return l;
    if (pt.fr.kt != zero3f) {
        auto ray = make_ray(pt.frame.o, -pt.wo);
        l += pt.fr.kt * trace_eyelight(scn, bvh, lights, ray, bounce + 1, pxl,
                            params, hit);
    }

    // done
    return l;
}

// Eyelight for quick previewing.
inline vec3f trace_eyelight(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const ray3f& ray, trace_pixel& pxl,
    const trace_params& params, bool& hit) {
    return trace_eyelight(scn, bvh, lights, ray, 0, pxl, params, hit);
}

// Debug previewing.
inline vec3f trace_debug_normal(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const ray3f& ray, trace_pixel& pxl,
    const trace_params& params, bool& hit) {
    // intersection
    auto isec = intersect_bvh(bvh, ray, false);
    hit = (bool)isec;
    if (!hit) return {0, 0, 0};

    // texcoord
    auto norm = eval_norm(scn->instances[isec.iid], isec.eid, isec.euv);
    return norm * 0.5f + vec3f{0.5f, 0.5f, 0.5f};
}

// Debug previewing.
inline vec3f trace_debug_albedo(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const ray3f& ray, trace_pixel& pxl,
    const trace_params& params, bool& hit) {
    // intersection
    auto pt = trace_intersect_scene(scn, bvh, ray);
    hit = pt.ist;

    return pt.fr.rho();
}

// Debug previewing.
inline vec3f trace_debug_texcoord(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const ray3f& ray, trace_pixel& pxl,
    const trace_params& params, bool& hit) {
    // intersection
    auto isec = intersect_bvh(bvh, ray, false);
    hit = (bool)isec;
    if (!hit) return {0, 0, 0};

    // texcoord
    auto texcoord =
        eval_texcoord(scn->instances[isec.iid]->shp, isec.eid, isec.euv);
    return {texcoord.x, texcoord.y, 0};
}

/// Trace shader function
using trace_shader = vec3f (*)(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const ray3f& ray, trace_pixel& pxl,
    const trace_params& params, bool& hit);

/// Trace filter function
using trace_filter = float (*)(float);

// map to convert trace samplers
static auto trace_shaders = unordered_map<trace_shader_type, trace_shader>{
    {trace_shader_type::eyelight, trace_eyelight},
    {trace_shader_type::direct, trace_direct},
    {trace_shader_type::pathtrace, trace_path},
    {trace_shader_type::pathtrace_nomis, trace_path_nomis},
    {trace_shader_type::debug_albedo, trace_debug_albedo},
    {trace_shader_type::debug_normal, trace_debug_normal},
    {trace_shader_type::debug_texcoord, trace_debug_texcoord},
};

// map to convert trace filters
static auto trace_filters = unordered_map<trace_filter_type, trace_filter>{
    {trace_filter_type::box, (trace_filter) nullptr},
    {trace_filter_type::triangle, filter_triangle},
    {trace_filter_type::cubic, filter_cubic},
    {trace_filter_type::catmull_rom, filter_catmullrom},
    {trace_filter_type::mitchell, filter_mitchell},
};

// map to convert trace filters
static auto trace_filter_sizes = unordered_map<trace_filter_type, int>{
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
    auto crn = trace_next2f(pxl);
    auto lrn = trace_next2f(pxl);
    auto uv = vec2f{
        (pxl.i + crn.x) / params.width, 1 - (pxl.j + crn.y) / params.height};
    auto ray = eval_camera_ray(cam, uv, lrn);
    bool hit = false;
    auto l = shader(scn, bvh, lights, ray, pxl, params, hit);
    if (!hit && params.envmap_invisible) return;
    if (!isfinite(l.x) || !isfinite(l.y) || !isfinite(l.z)) {
        log_error("NaN detected");
        return;
    }
    if (params.pixel_clamp > 0) l = clamplen(l, params.pixel_clamp);
    pxl.acc += {l, 1};
}

// Trace the next nsamples.
void trace_samples(const scene* scn, const camera* cam, const bvh_tree* bvh,
    const trace_lights& lights, image4f& img, image<trace_pixel>& pixels,
    int nsamples, const trace_params& params) {
    auto shader = trace_shaders.at(params.stype);
    if (params.parallel) {
        auto nthreads = std::thread::hardware_concurrency();
        auto threads = vector<std::thread>();
        for (auto tid = 0; tid < std::thread::hardware_concurrency(); tid++) {
            threads.push_back(std::thread([=, &img, &pixels, &params]() {
                for (auto j = tid; j < params.height; j += nthreads) {
                    for (auto i = 0; i < params.width; i++) {
                        auto& pxl = pixels.at(i, j);
                        for (auto s = 0; s < nsamples; s++)
                            trace_sample(
                                scn, cam, bvh, lights, pxl, shader, params);
                        img.at(i, j) = pxl.acc / pxl.sample;
                    }
                }
            }));
        }
        for (auto& t : threads) t.join();
        threads.clear();
    } else {
        auto shader = trace_shaders.at(params.stype);
        for (auto j = 0; j < params.height; j++) {
            for (auto i = 0; i < params.width; i++) {
                auto& pxl = pixels.at(i, j);
                for (auto s = 0; s < params.nsamples; s++)
                    trace_sample(scn, cam, bvh, lights, pxl, shader, params);
                img.at(i, j) = pxl.acc / pxl.sample;
            }
        }
    }
}

// Trace a filtered sample of samples
void trace_sample_filtered(const scene* scn, const camera* cam,
    const bvh_tree* bvh, const trace_lights& lights, trace_pixel& pxl,
    trace_shader shader, trace_filter filter, int filter_size,
    std::mutex& image_mutex, const trace_params& params) {
    pxl.sample += 1;
    pxl.dimension = 0;
    auto crn = trace_next2f(pxl);
    auto lrn = trace_next2f(pxl);
    auto uv = vec2f{
        (pxl.i + crn.x) / params.width, 1 - (pxl.j + crn.y) / params.height};
    auto ray = eval_camera_ray(cam, uv, lrn);
    auto hit = false;
    auto l = shader(scn, bvh, lights, ray, pxl, params, hit);
    if (!hit && params.envmap_invisible) return;
    if (!isfinite(l.x) || !isfinite(l.y) || !isfinite(l.z)) {
        log_error("NaN detected");
        return;
    }
    if (params.pixel_clamp > 0) l = clamplen(l, params.pixel_clamp);
    if (params.ftype == trace_filter_type::box) {
        pxl.acc += {l, 1};
        pxl.weight += 1;
    } else {
        std::lock_guard<std::mutex> lock(image_mutex);
        for (auto fj = max(0, pxl.j - filter_size);
             fj <= min(params.height - 1, pxl.j + filter_size); fj++) {
            for (auto fi = max(0, pxl.i - filter_size);
                 fi <= min(params.width - 1, pxl.i + filter_size); fi++) {
                auto w = filter((fi - pxl.i) - uv.x + 0.5f) *
                         filter((fj - pxl.j) - uv.y + 0.5f);
                pxl.acc += {l * w, w};
                pxl.weight += w;
            }
        }
    }
}

// Trace the next nsamples.
void trace_samples_filtered(const scene* scn, const camera* cam,
    const bvh_tree* bvh, const trace_lights& lights, image4f& img,
    image<trace_pixel>& pixels, int nsamples, const trace_params& params) {
    auto shader = trace_shaders.at(params.stype);
    auto filter = trace_filters.at(params.ftype);
    auto filter_size = trace_filter_sizes.at(params.ftype);
    std::mutex image_mutex;
    if (params.parallel) {
        auto nthreads = std::thread::hardware_concurrency();
        auto threads = vector<std::thread>();
        for (auto tid = 0; tid < std::thread::hardware_concurrency(); tid++) {
            threads.push_back(
                std::thread([=, &pixels, &params, &image_mutex]() {
                    for (auto j = tid; j < params.height; j += nthreads) {
                        for (auto i = 0; i < params.width; i++) {
                            auto& pxl = pixels.at(i, j);
                            for (auto s = 0; s < nsamples; s++) {
                                trace_sample_filtered(scn, cam, bvh, lights,
                                    pxl, shader, filter, filter_size,
                                    image_mutex, params);
                            }
                        }
                    }
                }));
        }
        for (auto& t : threads) t.join();
        threads.clear();
    } else {
        for (auto j = 0; j < params.height; j++) {
            for (auto i = 0; i < params.width; i++) {
                auto& pxl = pixels.at(i, j);
                for (auto s = 0; s < params.nsamples; s++) {
                    trace_sample_filtered(scn, cam, bvh, lights, pxl, shader,
                        filter, filter_size, image_mutex, params);
                }
            }
        }
    }
    for (auto j = 0; j < params.height; j++) {
        for (auto i = 0; i < params.width; i++) {
            img.at(i, j) = pixels.at(i, j).acc / pixels.at(i, j).weight;
        }
    }
}

// Starts an anyncrhounous renderer.
void trace_async_start(const scene* scn, const camera* cam, const bvh_tree* bvh,
    const trace_lights& lights, image4f& img, image<trace_pixel>& pixels,
    vector<std::thread>& threads, bool& stop_flag, const trace_params& params) {
    pixels = make_trace_pixels(params);
    auto nthreads = std::thread::hardware_concurrency();
    for (auto tid = 0; tid < std::thread::hardware_concurrency(); tid++) {
        threads.push_back(std::thread([=, &img, &pixels, &stop_flag]() {
            auto shader = trace_shaders.at(params.stype);
            for (auto s = 0; s < params.nsamples; s++) {
                for (auto j = tid; j < params.height; j += nthreads) {
                    for (auto i = 0; i < params.width; i++) {
                        if (stop_flag) return;
                        auto& pxl = pixels.at(i, j);
                        trace_sample(
                            scn, cam, bvh, lights, pxl, shader, params);
                        img.at(i, j) = pxl.acc / pxl.sample;
                    }
                }
            }
        }));
    }
}

// Stop the asynchronous renderer.
void trace_async_stop(vector<std::thread>& threads, bool& stop_flag) {
    stop_flag = true;
    for (auto& t : threads) t.join();
    stop_flag = false;
    threads.clear();
}

// Initialize trace lights
trace_lights make_trace_lights(const scene* scn) {
    auto lights = trace_lights();
    for (auto ist : scn->instances) {
        if (!ist->shp->mat) continue;
        if (ist->shp->mat->ke == zero3f) continue;
        lights.lights.push_back({ist, nullptr});
        if (!contains(lights.shape_cdfs, ist->shp)) {
            if (!ist->shp->points.empty()) {
                lights.shape_cdfs[ist->shp] = sample_points_cdf(ist->shp->points.size());
            } else if (!ist->shp->lines.empty()) {
                lights.shape_cdfs[ist->shp] = sample_lines_cdf(ist->shp->lines, ist->shp->pos);
            } else if (!ist->shp->triangles.empty()) {
                lights.shape_cdfs[ist->shp] = sample_triangles_cdf(ist->shp->triangles, ist->shp->pos);
            }
            lights.shape_areas[ist->shp] = lights.shape_cdfs[ist->shp].back();
        }
    }

    for (auto env : scn->environments) {
        if (env->ke == zero3f) continue;
        lights.lights.push_back({nullptr, env});
    }

    return lights;
}

// Initialize a rendering state
image<trace_pixel> make_trace_pixels(const trace_params& params) {
    auto pixels = image<trace_pixel>(params.width, params.height);
    for (auto j = 0; j < params.height; j++) {
        for (auto i = 0; i < params.width; i++) {
            pixels.at(i, j) = {zero4f,
                init_rng(params.seed, (j * params.width + i) * 2 + 1), i, j, 0,
                0, 0, params.nsamples, params.rtype};
        }
    }
    return pixels;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR WAVEFRONT OBJ
// -----------------------------------------------------------------------------
namespace ygl {

// Parse a value
template <typename T>
inline void obj_parse_val(stringstream& ss, T& v) {
    ss >> v;
}

// Parse a value
inline void obj_parse_val(stringstream& ss, string& v) {
    ss >> v;
    if (v.length() == 0) return;
    if ((v.front() == '"' && v.back() == '"') ||
        (v.front() == '\'' && v.back() == '\'')) {
        v = v.substr(1, v.length() - 2);
        return;
    }
}

// Parse a value
template <typename T, int N>
inline void obj_parse_val(stringstream& ss, vec<T, N>& v) {
    for (auto i = 0; i < N; i++) obj_parse_val(ss, v[i]);
}

// Parse a value
template <typename T, int N>
inline void obj_parse_val(stringstream& ss, quat<T, N>& v) {
    for (auto i = 0; i < N; i++) obj_parse_val(ss, v[i]);
}

// Parse a value
template <typename T>
inline void obj_parse_val(stringstream& ss, frame<T, 3>& v) {
    obj_parse_val(ss, v.x);
    obj_parse_val(ss, v.y);
    obj_parse_val(ss, v.z);
    obj_parse_val(ss, v.o);
}

// Parse texture options and name
inline void parse_texture(stringstream& ss, obj_texture_info& info,
    vector<string>& textures, unordered_set<string>& texture_set,
    bool bump = false) {
    // get tokens
    auto tokens = vector<string>();
    while (ss) {
        auto s = string();
        ss >> s;
        if (!s.empty()) tokens.push_back(s);
    }

    // exit if no tokens
    if (tokens.empty()) return;

    // texture name
    info.path = tokens.back();
    for (auto& c : info.path)
        if (c == '\\') c = '/';

    // texture options
    auto last = string();
    for (auto& tok : tokens) {
        if (tok == tokens.back()) break;
        if (tok[0] == '-') {
            last = tok;
            info.unknown_props[last] = {};
        } else {
            info.unknown_props[last].push_back(tok);
        }
    }

    // clamp
    if (info.unknown_props.find("-clamp") != info.unknown_props.end() &&
        info.unknown_props.at("-clamp").size() > 0) {
        auto& clamp_vec = info.unknown_props.at("-clamp");
        auto clamp_str = (clamp_vec.empty()) ? "" : clamp_vec.front();
        info.clamp = clamp_str == "on" || clamp_str == "1";
        info.unknown_props.erase("-clamp");
    }

    if (info.unknown_props.find("-bm") != info.unknown_props.end() &&
        info.unknown_props.at("-bm").size() > 0) {
        auto& bm_vec = info.unknown_props.at("-bm");
        auto bm_str = (bm_vec.empty()) ? "" : bm_vec.front();
        info.scale = std::atof(bm_str.c_str());
        info.unknown_props.erase("-bm");
    }

    // insert texture
    if (!info.path.empty() &&
        texture_set.find(info.path) == texture_set.end()) {
        textures.push_back(info.path);
        texture_set.insert(info.path);
    }
}

// Load MTL
inline vector<obj_material*> load_mtl(
    const string& filename, bool flip_tr, vector<string>& textures) {
    // clear materials
    auto materials = vector<obj_material*>();

    // clear textures
    textures.clear();
    auto texture_set = unordered_set<string>();

    // open file
    auto fs = fstream(filename, ios_base::in);
    if (!fs) throw runtime_error("cannot open filename " + filename);
    // fs.exceptions(ios_base::failbit);

    // add a material preemptively to avoid crashes
    materials.push_back(new obj_material());

    // read the file line by line
    string line;
    auto linenum = 0;
    while (getline(fs, line)) {
        // prepare to parse
        linenum += 1;
        auto ss = stringstream(line);
        auto cmd = string();
        ss >> cmd;

        // skip empty and comments
        if (cmd.empty() || cmd[0] == '#') continue;

        // possible token values
        if (cmd == "newmtl") {
            materials.push_back(new obj_material());
            obj_parse_val(ss, materials.back()->name);
        } else if (cmd == "illum") {
            obj_parse_val(ss, materials.back()->illum);
        } else if (cmd == "Ke") {
            obj_parse_val(ss, materials.back()->ke);
        } else if (cmd == "Ka") {
            obj_parse_val(ss, materials.back()->ka);
        } else if (cmd == "Kd") {
            obj_parse_val(ss, materials.back()->kd);
        } else if (cmd == "Ks") {
            obj_parse_val(ss, materials.back()->ks);
        } else if (cmd == "Kr") {
            obj_parse_val(ss, materials.back()->kr);
        } else if (cmd == "Kt" || cmd == "Tf") {
            auto vals = zero3f;
            auto ntok = 0;
            while (ss) obj_parse_val(ss, vals[ntok++]);
            if (ntok >= 3)
                materials.back()->kt = vals;
            else
                materials.back()->kt = {vals.x, vals.x, vals.x};
        } else if (cmd == "Tr") {
            auto vals = zero3f;
            auto ntok = 0;
            while (ss) obj_parse_val(ss, vals[ntok++]);
            if (ntok >= 3)
                materials.back()->kt = vals;
            else
                materials.back()->op = (flip_tr) ? 1 - vals.x : vals.x;
        } else if (cmd == "Ns") {
            obj_parse_val(ss, materials.back()->ns);
        } else if (cmd == "d") {
            obj_parse_val(ss, materials.back()->op);
        } else if (cmd == "Ni") {
            obj_parse_val(ss, materials.back()->ior);
        } else if (cmd == "map_Ke") {
            parse_texture(ss, materials.back()->ke_txt, textures, texture_set);
        } else if (cmd == "map_Ka") {
            parse_texture(ss, materials.back()->ka_txt, textures, texture_set);
        } else if (cmd == "map_Kd") {
            parse_texture(ss, materials.back()->kd_txt, textures, texture_set);
        } else if (cmd == "map_Ks") {
            parse_texture(ss, materials.back()->ks_txt, textures, texture_set);
        } else if (cmd == "map_Kr") {
            parse_texture(ss, materials.back()->kr_txt, textures, texture_set);
        } else if (cmd == "map_Tr") {
            parse_texture(ss, materials.back()->kt_txt, textures, texture_set);
        } else if (cmd == "map_Ns") {
            parse_texture(ss, materials.back()->ns_txt, textures, texture_set);
        } else if (cmd == "map_d") {
            parse_texture(ss, materials.back()->op_txt, textures, texture_set);
        } else if (cmd == "map_Ni") {
            parse_texture(ss, materials.back()->ior_txt, textures, texture_set);
        } else if (cmd == "map_bump" || cmd == "bump") {
            parse_texture(
                ss, materials.back()->bump_txt, textures, texture_set);
        } else if (cmd == "map_disp" || cmd == "disp") {
            parse_texture(
                ss, materials.back()->disp_txt, textures, texture_set);
        } else if (cmd == "map_norm" || cmd == "norm") {
            parse_texture(
                ss, materials.back()->norm_txt, textures, texture_set);
        } else {
            // copy into strings
            while (ss) {
                materials.back()->unknown_props[cmd].push_back({});
                obj_parse_val(ss, materials.back()->unknown_props[cmd].back());
            }
        }
    }

    // remove first fake material
    materials.erase(materials.begin());

    // done
    return materials;
}

// Loads textures for an scene.
inline void load_textures(
    obj_scene* asset, const string& dirname, bool skip_missing) {
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
            throw runtime_error("cannot laod image " + filename);
        }
    }
}

// Parses an OBJ vertex list. Handles negative values.
inline void obj_parse_vertlist(
    stringstream& ss, vector<obj_vertex>& elems, const obj_vertex& vert_size) {
    elems.clear();
    while (true) {
        auto tok = string();
        obj_parse_val(ss, tok);
        if (tok.empty()) break;
        auto toks = split(tok, "/");
        if (toks.empty()) break;
        auto v = obj_vertex{-1, -1, -1, -1, -1};
        for (auto i = 0; i < min(5, (int)toks.size()); i++) {
            if (toks[i] == "") continue;
            ((int*)&v)[i] = (int)atoi(toks[i].c_str());
            ((int*)&v)[i] = (((int*)&v)[i] < 0) ?
                                ((int*)&vert_size)[i] + ((int*)&v)[i] :
                                ((int*)&v)[i] - 1;
        }
        elems.push_back(v);
    }
}

// Loads an OBJ
inline obj_scene* load_obj(const string& filename, bool load_txt,
    bool skip_missing, bool flip_texcoord, bool flip_tr) {
    // clear obj
    auto asset = unique_ptr<obj_scene>(new obj_scene());

    // open file
    auto fs = fstream(filename, ios_base::in);
    if (!fs) throw runtime_error("cannot open filename " + filename);
    // fs.exceptions(ios_base::failbit);

    // initializing obj
    asset->objects.push_back(new obj_object());
    asset->objects.back()->groups.push_back({});

    // allocate buffers to avoid re-allocing
    auto cur_elems = vector<obj_vertex>();
    auto cur_matname = string();
    auto cur_mtllibs = vector<string>();

    // keep track of array lengths
    auto vert_size = obj_vertex{0, 0, 0, 0, 0};

    // read the file line by line
    string line;
    auto linenum = 0;
    while (getline(fs, line)) {
        // prepare to parse
        linenum += 1;
        auto ss = stringstream(line);
        auto cmd = string();
        ss >> cmd;

        // skip empty and comments
        if (cmd.empty() || cmd[0] == '#') continue;

        // possible token values
        if (cmd == "v") {
            vert_size.pos += 1;
            asset->pos.push_back({});
            obj_parse_val(ss, asset->pos.back());
        } else if (cmd == "vn") {
            vert_size.norm += 1;
            asset->norm.push_back({});
            obj_parse_val(ss, asset->norm.back());
        } else if (cmd == "vt") {
            vert_size.texcoord += 1;
            asset->texcoord.push_back({});
            obj_parse_val(ss, asset->texcoord.back());
            if (flip_texcoord)
                asset->texcoord.back().y = 1 - asset->texcoord.back().y;
        } else if (cmd == "vc") {
            vert_size.color += 1;
            asset->color.push_back({});
            obj_parse_val(ss, asset->color.back());
        } else if (cmd == "vr") {
            vert_size.radius += 1;
            asset->radius.push_back({});
            obj_parse_val(ss, asset->radius.back());
        } else if (cmd == "f") {
            obj_parse_vertlist(ss, cur_elems, vert_size);
            auto& g = asset->objects.back()->groups.back();
            g.elems.push_back({(uint32_t)g.verts.size(), obj_element_type::face,
                (uint16_t)cur_elems.size()});
            g.verts.insert(g.verts.end(), cur_elems.begin(), cur_elems.end());
        } else if (cmd == "l") {
            obj_parse_vertlist(ss, cur_elems, vert_size);
            auto& g = asset->objects.back()->groups.back();
            g.elems.push_back({(uint32_t)g.verts.size(), obj_element_type::line,
                (uint16_t)cur_elems.size()});
            g.verts.insert(g.verts.end(), cur_elems.begin(), cur_elems.end());
        } else if (cmd == "p") {
            obj_parse_vertlist(ss, cur_elems, vert_size);
            auto& g = asset->objects.back()->groups.back();
            g.elems.push_back({(uint32_t)g.verts.size(),
                obj_element_type::point, (uint16_t)cur_elems.size()});
            g.verts.insert(g.verts.end(), cur_elems.begin(), cur_elems.end());
        } else if (cmd == "b") {
            obj_parse_vertlist(ss, cur_elems, vert_size);
            auto& g = asset->objects.back()->groups.back();
            g.elems.push_back({(uint32_t)g.verts.size(),
                obj_element_type::bezier, (uint16_t)cur_elems.size()});
            g.verts.insert(g.verts.end(), cur_elems.begin(), cur_elems.end());
        } else if (cmd == "t") {
            obj_parse_vertlist(ss, cur_elems, vert_size);
            auto& g = asset->objects.back()->groups.back();
            g.elems.push_back({(uint32_t)g.verts.size(),
                obj_element_type::tetra, (uint16_t)cur_elems.size()});
            g.verts.insert(g.verts.end(), cur_elems.begin(), cur_elems.end());
        } else if (cmd == "o") {
            auto name = string();
            obj_parse_val(ss, name);
            asset->objects.push_back(new obj_object{name, {}});
            asset->objects.back()->groups.push_back({});
            asset->objects.back()->groups.back().matname = cur_matname;
        } else if (cmd == "usemtl") {
            auto name = string();
            obj_parse_val(ss, name);
            cur_matname = name;
            asset->objects.back()->groups.push_back({});
            asset->objects.back()->groups.back().matname = cur_matname;
        } else if (cmd == "g") {
            auto name = string();
            obj_parse_val(ss, name);
            asset->objects.back()->groups.push_back({});
            asset->objects.back()->groups.back().matname = cur_matname;
            asset->objects.back()->groups.back().groupname = name;
        } else if (cmd == "s") {
            auto name = string();
            obj_parse_val(ss, name);
            auto smoothing = (name == "on");
            if (asset->objects.back()->groups.empty()) {
                asset->objects.back()->groups.push_back({});
                asset->objects.back()->groups.back().matname = cur_matname;
                asset->objects.back()->groups.back().smoothing = smoothing;
            } else if (asset->objects.back()->groups.back().smoothing !=
                       smoothing) {
                auto gname = asset->objects.back()->groups.back().groupname;
                asset->objects.back()->groups.push_back({});
                asset->objects.back()->groups.back().matname = cur_matname;
                asset->objects.back()->groups.back().groupname = gname;
                asset->objects.back()->groups.back().smoothing = smoothing;
            }
        } else if (cmd == "sl") {
            auto subdiv = zero2i;
            obj_parse_val(ss, subdiv);
            if (asset->objects.back()->groups.empty()) {
                asset->objects.back()->groups.push_back({});
                asset->objects.back()->groups.back().matname = cur_matname;
            }
            asset->objects.back()->groups.back().subdivision_level = subdiv.x;
            asset->objects.back()->groups.back().subdivision_catmullclark =
                (bool)subdiv.y;
        } else if (cmd == "mtllib") {
            auto name = string();
            obj_parse_val(ss, name);
            if (name != string("")) {
                auto found = false;
                for (auto lib : cur_mtllibs) {
                    if (lib == name) {
                        found = true;
                        break;
                    }
                }
                if (!found) cur_mtllibs.push_back(name);
            }
        } else if (cmd == "c") {
            auto cam = new obj_camera();
            obj_parse_val(ss, cam->name);
            obj_parse_val(ss, cam->ortho);
            obj_parse_val(ss, cam->yfov);
            obj_parse_val(ss, cam->aspect);
            obj_parse_val(ss, cam->aperture);
            obj_parse_val(ss, cam->focus);
            obj_parse_val(ss, cam->frame);
            asset->cameras.push_back(cam);
        } else if (cmd == "e") {
            auto env = new obj_environment();
            obj_parse_val(ss, env->name);
            obj_parse_val(ss, env->matname);
            obj_parse_val(ss, env->frame);
            asset->environments.push_back(env);
        } else if (cmd == "n") {
            auto nde = new obj_node();
            obj_parse_val(ss, nde->name);
            obj_parse_val(ss, nde->parent);
            obj_parse_val(ss, nde->camname);
            obj_parse_val(ss, nde->objname);
            obj_parse_val(ss, nde->envname);
            obj_parse_val(ss, nde->frame);
            obj_parse_val(ss, nde->translation);
            obj_parse_val(ss, nde->rotation);
            obj_parse_val(ss, nde->scaling);
            asset->nodes.push_back(nde);
        } else {
            // unused
        }
    }

    // cleanup unused
    for (auto o : asset->objects) {
        auto end = std::remove_if(o->groups.begin(), o->groups.end(),
            [](const obj_group& x) { return x.verts.empty(); });
        o->groups.erase(end, o->groups.end());
    }
    // TODO: possible memory leak
    auto end = std::remove_if(asset->objects.begin(), asset->objects.end(),
        [](const obj_object* x) { return x->groups.empty(); });
    asset->objects.erase(end, asset->objects.end());

    // parse materials
    auto dirname = path_dirname(filename);
    unordered_set<string> texture_set;
    for (auto mtllib : cur_mtllibs) {
        auto mtlname = dirname + mtllib;
        vector<string> textures;
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

    // done
    return asset.release();
}

// write to stream
template <typename T>
inline void obj_dump_val(fstream& fs, const T& v) {
    fs << v;
}

// write to stream
template <typename T, int N>
inline void obj_dump_val(fstream& fs, const vec<T, N>& v) {
    for (auto i = 0; i < N; i++) {
        if (i) fs << ' ';
        obj_dump_val(fs, v[i]);
    }
}

// write to stream
template <typename T, int N>
inline void obj_dump_val(fstream& fs, const quat<T, N>& v) {
    for (auto i = 0; i < N; i++) {
        if (i) fs << ' ';
        obj_dump_val(fs, v[i]);
    }
}

// write to stream
template <typename T>
inline void obj_dump_val(fstream& fs, const frame<T, 3>& v) {
    obj_dump_val(fs, v.x);
    fs << ' ';
    obj_dump_val(fs, v.y);
    fs << ' ';
    obj_dump_val(fs, v.z);
    fs << ' ';
    obj_dump_val(fs, v.o);
}

// write to stream
inline void obj_dump_val(fstream& fs, const obj_texture_info& v) {
    for (auto&& kv : v.unknown_props) {
        obj_dump_val(fs, kv.first + " ");
        for (auto&& vv : kv.second) obj_dump_val(fs, vv + " ");
    }
    if (v.clamp) obj_dump_val(fs, "-clamp on ");
    obj_dump_val(fs, v.path);
}

// write to stream
template <typename T>
inline void obj_dump_named_val(fstream& fs, const string& name, const T& v) {
    obj_dump_val(fs, name);
    fs << ' ';
    obj_dump_val(fs, v);
    fs << '\n';
}

// write to stream
template <typename T>
inline void obj_dump_opt_val(
    fstream& fs, const string& name, const T& v, const T& def = {}) {
    if (v == def) return;
    obj_dump_named_val(fs, name, v);
}

// write an OBJ vertex triplet using only the indices that are active
inline void obj_dump_objverts(
    fstream& fs, const char* str, int nv, const obj_vertex* verts) {
    obj_dump_val(fs, str);
    for (auto v = 0; v < nv; v++) {
        auto& vert = verts[v];
        auto vert_ptr = &vert.pos;
        auto nto_write = 0;
        for (auto i = 0; i < 5; i++) {
            if (vert_ptr[i] >= 0) nto_write = i + 1;
        }
        for (auto i = 0; i < nto_write; i++) {
            if (vert_ptr[i] >= 0) {
                obj_dump_val(fs, ((i == 0) ? ' ' : '/'));
                obj_dump_val(fs, vert_ptr[i] + 1);
            } else {
                obj_dump_val(fs, '/');
            }
        }
    }
    obj_dump_val(fs, "\n");
}

// Save an MTL file
inline void save_mtl(const string& filename,
    const vector<obj_material*>& materials, bool flip_tr) {
    // open file
    auto fs = fstream(filename, ios_base::out);
    if (!fs) throw runtime_error("cannot open filename " + filename);
    fs.exceptions(ios_base::failbit);

    // for each material, dump all the values
    for (auto mat : materials) {
        obj_dump_named_val(fs, "newmtl", mat->name);
        obj_dump_named_val(fs, "  illum", mat->illum);
        obj_dump_opt_val(fs, "  Ke", mat->ke);
        obj_dump_opt_val(fs, "  Ka", mat->ka);
        obj_dump_opt_val(fs, "  Kd", mat->kd);
        obj_dump_opt_val(fs, "  Ks", mat->ks);
        obj_dump_opt_val(fs, "  Kr", mat->kr);
        obj_dump_opt_val(fs, "  Tf", mat->kt);
        obj_dump_opt_val(fs, "  Ns", mat->ns, 0.0f);
        obj_dump_opt_val(fs, "  d", mat->op, 1.0f);
        obj_dump_opt_val(fs, "  Ni", mat->ior, 1.0f);
        obj_dump_opt_val(fs, "  map_Ke", mat->ke_txt);
        obj_dump_opt_val(fs, "  map_Ka", mat->ka_txt);
        obj_dump_opt_val(fs, "  map_Kd", mat->kd_txt);
        obj_dump_opt_val(fs, "  map_Ks", mat->ks_txt);
        obj_dump_opt_val(fs, "  map_Kr", mat->kr_txt);
        obj_dump_opt_val(fs, "  map_Kt", mat->kt_txt);
        obj_dump_opt_val(fs, "  map_Ns", mat->ns_txt);
        obj_dump_opt_val(fs, "  map_d", mat->op_txt);
        obj_dump_opt_val(fs, "  map_Ni", mat->ior_txt);
        obj_dump_opt_val(fs, "  map_bump", mat->bump_txt);
        obj_dump_opt_val(fs, "  map_disp", mat->disp_txt);
        obj_dump_opt_val(fs, "  map_norm", mat->norm_txt);
        for (auto&& kv : mat->unknown_props) {
            obj_dump_val(fs, kv.first);
            for (auto&& v : kv.second) {
                obj_dump_val(fs, " ");
                obj_dump_val(fs, v);
            }
            obj_dump_val(fs, "\n");
        }
        obj_dump_val(fs, "\n");
    }
}

// Loads textures for an scene.
inline void save_textures(
    const obj_scene* asset, const string& dirname, bool skip_missing) {
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
            throw runtime_error("cannot save image " + filename);
        }
    }
}

// Save an OBJ
inline void save_obj(const string& filename, const obj_scene* asset,
    bool save_txt, bool skip_missing, bool flip_texcoord, bool flip_tr) {
    // open file
    auto fs = fstream(filename, ios_base::out);
    if (!fs) throw runtime_error("cannot open filename " + filename);
    fs.exceptions(ios_base::failbit);

    // linkup to mtl
    auto dirname = path_dirname(filename);
    auto basename = filename.substr(dirname.length());
    basename = basename.substr(0, basename.length() - 4);
    if (!asset->materials.empty()) {
        obj_dump_named_val(fs, "mtllib", basename + ".mtl");
    }

    // save cameras
    for (auto cam : asset->cameras) {
        obj_dump_val(fs, "c ");
        obj_dump_val(fs, cam->name);
        obj_dump_val(fs, " ");
        obj_dump_val(fs, cam->ortho);
        obj_dump_val(fs, " ");
        obj_dump_val(fs, cam->yfov);
        obj_dump_val(fs, " ");
        obj_dump_val(fs, cam->aspect);
        obj_dump_val(fs, " ");
        obj_dump_val(fs, cam->aperture);
        obj_dump_val(fs, " ");
        obj_dump_val(fs, cam->focus);
        obj_dump_val(fs, " ");
        obj_dump_val(fs, cam->frame);
        obj_dump_val(fs, "\n");
    }

    // save envs
    for (auto env : asset->environments) {
        obj_dump_val(fs, "e ");
        obj_dump_val(fs, env->name);
        obj_dump_val(fs, " ");
        obj_dump_val(fs, env->matname);
        obj_dump_val(fs, " ");
        obj_dump_val(fs, env->frame);
        obj_dump_val(fs, "\n");
    }

    // save nodes
    for (auto nde : asset->nodes) {
        obj_dump_val(fs, "n ");
        obj_dump_val(fs, nde->name);
        obj_dump_val(fs, " ");
        obj_dump_val(fs, (nde->parent.empty()) ? string("\"\"") : nde->parent);
        obj_dump_val(fs, " ");
        obj_dump_val(
            fs, (nde->camname.empty()) ? string("\"\"") : nde->camname);
        obj_dump_val(fs, " ");
        obj_dump_val(
            fs, (nde->objname.empty()) ? string("\"\"") : nde->objname);
        obj_dump_val(fs, " ");
        obj_dump_val(
            fs, (nde->envname.empty()) ? string("\"\"") : nde->envname);
        obj_dump_val(fs, " ");
        obj_dump_val(fs, nde->frame);
        obj_dump_val(fs, " ");
        obj_dump_val(fs, nde->translation);
        obj_dump_val(fs, " ");
        obj_dump_val(fs, nde->rotation);
        obj_dump_val(fs, " ");
        obj_dump_val(fs, nde->scaling);
        obj_dump_val(fs, "\n");
    }

    // save all vertex data
    for (auto& v : asset->pos) obj_dump_named_val(fs, "v", v);
    if (flip_texcoord) {
        for (auto& v : asset->texcoord)
            obj_dump_named_val(fs, "vt", vec2f{v.x, 1 - v.y});
    } else {
        for (auto& v : asset->texcoord) obj_dump_named_val(fs, "vt", v);
    }
    for (auto& v : asset->norm) obj_dump_named_val(fs, "vn", v);
    for (auto& v : asset->color) obj_dump_named_val(fs, "vc", v);
    for (auto& v : asset->radius) obj_dump_named_val(fs, "vr", v);

    // save element data
    static auto elem_labels =
        unordered_map<obj_element_type, string>{{obj_element_type::point, "p"},
            {obj_element_type::line, "l"}, {obj_element_type::face, "f"},
            {obj_element_type::bezier, "b"}, {obj_element_type::tetra, "t"}};
    for (auto object : asset->objects) {
        obj_dump_named_val(fs, "o", object->name);
        for (auto& group : object->groups) {
            obj_dump_opt_val(fs, "usemtl", group.matname);
            obj_dump_opt_val(fs, "g", group.groupname);
            if (!group.smoothing) obj_dump_named_val(fs, "s", "off");
            if (group.subdivision_level) {
                auto sl = vec2i{group.subdivision_level,
                    (group.subdivision_catmullclark) ? 1 : 0};
                obj_dump_named_val(fs, "sl", sl);
            }
            for (auto elem : group.elems) {
                auto lbl = "";
                switch (elem.type) {
                    case obj_element_type::point: lbl = "p"; break;
                    case obj_element_type::line: lbl = "l"; break;
                    case obj_element_type::face: lbl = "f"; break;
                    case obj_element_type::bezier: lbl = "b"; break;
                    case obj_element_type::tetra: lbl = "t"; break;
                    default: throw runtime_error("should not have gotten here");
                }
                obj_dump_objverts(
                    fs, lbl, elem.size, group.verts.data() + elem.start);
            }
        }
    }

    // save materials
    if (!asset->materials.empty())
        save_mtl(dirname + basename + ".mtl", asset->materials, flip_tr);

    // save textures
    if (save_txt) save_textures(asset, dirname, skip_missing);
}

// A hash function for vecs
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

// Flattens an scene
inline obj_mesh* get_mesh(
    const obj_scene* model, const obj_object& oshape, bool facet_non_smooth) {
    // convert meshes
    auto msh = new obj_mesh();
    msh->name = oshape.name;
    for (auto& group : oshape.groups) {
        if (group.verts.empty()) continue;
        if (group.elems.empty()) continue;
        msh->shapes.emplace_back();
        auto prim = &msh->shapes.back();
        prim->name = group.groupname;
        prim->matname = group.matname;

        // insert all vertices
        unordered_map<obj_vertex, int, obj_vertex_hash> vert_map;
        vector<int> vert_ids;
        // vert_map.clear();
        // vert_ids.clear();
        for (auto& vert : group.verts) {
            if (vert_map.find(vert) == vert_map.end()) {
                // split in two to avoid undefined behaviour
                auto size = (int)vert_map.size();
                vert_map[vert] = size;
            }
            vert_ids.push_back(vert_map.at(vert));
        }

        // convert elements
        for (auto& elem : group.elems) {
            switch (elem.type) {
                case obj_element_type::point: {
                    for (auto i = elem.start; i < elem.start + elem.size; i++) {
                        prim->points.push_back(vert_ids[i]);
                    }
                } break;
                case obj_element_type::line: {
                    for (auto i = elem.start; i < elem.start + elem.size - 1;
                         i++) {
                        prim->lines.push_back({vert_ids[i], vert_ids[i + 1]});
                    }
                } break;
                case obj_element_type::face: {
                    for (auto i = elem.start + 2; i < elem.start + elem.size;
                         i++) {
                        prim->triangles.push_back({vert_ids[elem.start],
                            vert_ids[i - 1], vert_ids[i]});
                    }
                } break;
                case obj_element_type::bezier: {
                    if ((elem.size - 1) % 3)
                        throw runtime_error("bad obj bezier");
                    for (auto i = elem.start + 1; i < elem.start + elem.size;
                         i += 3) {
                        prim->bezier.push_back({vert_ids[i - 1], vert_ids[i],
                            vert_ids[i + 1], vert_ids[i + 2]});
                    }
                } break;
                case obj_element_type::tetra: {
                    for (auto i = elem.start; i < elem.start + elem.size;
                         i += 4) {
                        if (i + 3 >= vert_ids.size()) continue;
                        prim->tetras.push_back({vert_ids[i], vert_ids[i + 1],
                            vert_ids[i + 2], vert_ids[i + 3]});
                    }
                } break;
                default: { assert(false); }
            }
        }

        // check for errors
        // copy vertex data
        auto v = group.verts[0];
        if (v.pos >= 0) prim->pos.resize(vert_map.size());
        if (v.texcoord >= 0) prim->texcoord.resize(vert_map.size());
        if (v.norm >= 0) prim->norm.resize(vert_map.size());
        if (v.color >= 0) prim->color.resize(vert_map.size());
        if (v.radius >= 0) prim->radius.resize(vert_map.size());
        for (auto& kv : vert_map) {
            if (v.pos >= 0 && kv.first.pos >= 0) {
                prim->pos[kv.second] = model->pos[kv.first.pos];
            }
            if (v.texcoord >= 0 && kv.first.texcoord >= 0) {
                prim->texcoord[kv.second] = model->texcoord[kv.first.texcoord];
            }
            if (v.norm >= 0 && kv.first.norm >= 0) {
                prim->norm[kv.second] = model->norm[kv.first.norm];
            }
            if (v.color >= 0 && kv.first.color >= 0) {
                prim->color[kv.second] = model->color[kv.first.color];
            }
            if (v.radius >= 0 && kv.first.radius >= 0) {
                prim->radius[kv.second] = model->radius[kv.first.radius];
            }
        }

        // fix smoothing
        if (!group.smoothing && facet_non_smooth) {
            auto faceted_ = obj_shape();
            auto faceted = &faceted_;
            faceted->name = prim->name;
            faceted->matname = prim->matname;
            auto pidx = vector<int>();
            for (auto point : prim->points) {
                faceted->points.push_back((int)pidx.size());
                pidx.push_back(point);
            }
            for (auto line : prim->lines) {
                faceted->lines.push_back(
                    {(int)pidx.size() + 0, (int)pidx.size() + 1});
                pidx.push_back(line.x);
                pidx.push_back(line.y);
            }
            for (auto triangle : prim->triangles) {
                faceted->triangles.push_back({(int)pidx.size() + 0,
                    (int)pidx.size() + 1, (int)pidx.size() + 2});
                pidx.push_back(triangle.x);
                pidx.push_back(triangle.y);
                pidx.push_back(triangle.z);
            }
            for (auto tetra : prim->tetras) {
                faceted->tetras.push_back(
                    {(int)pidx.size() + 0, (int)pidx.size() + 1,
                        (int)pidx.size() + 2, (int)pidx.size() + 3});
                pidx.push_back(tetra.x);
                pidx.push_back(tetra.y);
                pidx.push_back(tetra.z);
                pidx.push_back(tetra.w);
            }
            for (auto idx : pidx) {
                if (!prim->pos.empty()) faceted->pos.push_back(prim->pos[idx]);
                if (!prim->norm.empty())
                    faceted->norm.push_back(prim->norm[idx]);
                if (!prim->texcoord.empty())
                    faceted->texcoord.push_back(prim->texcoord[idx]);
                if (!prim->color.empty())
                    faceted->color.push_back(prim->color[idx]);
                if (!prim->radius.empty())
                    faceted->radius.push_back(prim->radius[idx]);
            }
            *prim = *faceted;
        }
    }

    // done
    return msh;
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
inline void serialize_from_json(int& val, const json& js) {
    if (!js.is_number_integer()) throw runtime_error("integer expected");
    val = js;
}

// Parse float function.
inline void serialize_from_json(float& val, const json& js) {
    if (!js.is_number()) throw runtime_error("number expected");
    val = js;
}

// Parse bool function.
inline void serialize_from_json(bool& val, const json& js) {
    if (!js.is_boolean()) throw runtime_error("bool expected");
    val = js;
}

// Parse std::string function.
inline void serialize_from_json(string& val, const json& js) {
    if (!js.is_string()) throw runtime_error("string expected");
    val = js;
}

// Parse json function.
inline void serialize_from_json(json& val, const json& js) { val = js; }

// Parse support function.
template <typename T>
inline void serialize_from_json(T*& val, const json& js) {
    if (js.is_null()) {
        val = nullptr;
        return;
    }
    if (!js.is_object()) throw runtime_error("object expected");
    if (!val) val = new T();
    serialize_from_json(*val, js);
}

// Parse support function.
template <typename T>
inline void serialize_from_json(vector<T>& vals, const json& js) {
    if (!js.is_array()) throw runtime_error("array expected");
    vals.resize(js.size());
    for (auto i = 0; i < js.size(); i++) {
        // this is contrived to support for vector<bool>
        auto v = T();
        serialize_from_json(v, js[i]);
        vals[i] = v;
    }
}

// Parse support function.
template <typename T>
inline void serialize_from_json(map<string, T>& vals, const json& js) {
    if (!js.is_object()) throw runtime_error("object expected");
    for (auto kv = js.begin(); kv != js.end(); ++kv) {
        serialize_from_json(vals[kv.key()], kv.value());
    }
}

// Parse support function.
template <typename T, size_t N>
inline void serialize_from_json(array<T, N>& vals, const json& js) {
    if (!js.is_array()) throw runtime_error("array expected");
    if (N != js.size()) throw runtime_error("wrong array size");
    for (auto i = 0; i < N; i++) serialize_from_json(vals[i], js.at(i));
}

// Parse support function.
template <typename T, typename T1>
inline void serialize_from_json(
    T& val, const json& js, const vector<pair<T1, T>>& table) {
    auto v = T1();
    serialize_from_json(v, js);
    auto found = false;
    for (auto& kv : table) {
        if (kv.first == v) {
            val = kv.second;
            found = true;
            break;
        }
    }
    if (!found) throw runtime_error("bad enum value");
}

// Parse support function.
inline void serialize_from_json(vec2f& vals, const json& js) {
    serialize_from_json((array<float, 2>&)vals, js);
}

// Parse support function.
inline void serialize_from_json(vec3f& vals, const json& js) {
    serialize_from_json((array<float, 3>&)vals, js);
}

// Parse support function.
inline void serialize_from_json(vec4f& vals, const json& js) {
    serialize_from_json((array<float, 4>&)vals, js);
}

// Parse support function.
inline void serialize_from_json(quat4f& vals, const json& js) {
    serialize_from_json((array<float, 4>&)vals, js);
}

// Parse support function.
inline void serialize_from_json(mat4f& vals, const json& js) {
    serialize_from_json((array<float, 16>&)vals, js);
}

// Parse support function.
inline void serialize_from_json(frame3f& vals, const json& js) {
    serialize_from_json((array<float, 12>&)vals, js);
}

// Parse support function.
inline void serialize_from_json_obj(const json& js) {
    if (!js.is_object()) throw runtime_error("object expected");
}

// Parse support function.
template <typename T>
inline void serialize_from_json_attr(T& val, const json& js, const char* name,
    bool required = true, const T& def = {}) {
    if (required) {
        if (!js.count(name)) throw runtime_error("missing value");
        serialize_from_json(val, js.at(name));
    } else {
        if (!js.count(name))
            val = def;
        else
            serialize_from_json(val, js.at(name));
    }
}

// Converts int to json.
inline void serialize_to_json(int val, json& js) { js = val; }

// Converts float to json.
inline void serialize_to_json(float val, json& js) { js = val; }

// Converts bool to json.
inline void serialize_to_json(bool val, json& js) { js = val; }

// Converts string to json.
inline void serialize_to_json(const string& val, json& js) { js = val; }

// Converts json to json.
inline void serialize_to_json(const json& val, json& js) { js = val; }

// Dump support function.
template <typename T>
inline void serialize_to_json(const T* val, json& js) {
    if (!val) {
        js = nullptr;
        return;
    }
    if (!js.is_object()) js = json::object();
    serialize_to_json(*val, js);
}

// Dump support function.
template <typename T, size_t N>
inline void serialize_to_json(const array<T, N>& vals, json& js) {
    js = json::array();
    for (auto i = 0; i < N; i++) serialize_to_json(vals[i], js[i]);
}

// Dump support function.
template <typename T>
inline void serialize_to_json(const vector<T>& vals, json& js) {
    js = json::array();
    for (auto i = 0; i < vals.size(); i++) serialize_to_json(vals[i], js[i]);
}

// Dump support function.
template <typename T>
inline void serialize_to_json(const map<string, T>& vals, json& js) {
    js = json::object();
    for (auto& kv : vals) serialize_to_json(kv.second, js[kv.first]);
}

// Dump support function.
template <typename T, typename T1>
inline void serialize_to_json(
    const T& val, json& js, const vector<pair<T1, T>>& table) {
    auto found = false;
    auto v = T1();
    for (auto& kv : table) {
        if (kv.second == val) {
            v = kv.first;
            found = true;
            break;
        }
    }
    if (!found) throw runtime_error("invalid value");
    serialize_to_json(v, js);
}

// Dump support function.
inline void serialize_to_json(const vec2f& vals, json& js) {
    serialize_to_json((const array<float, 2>&)vals, js);
}

// Dump support function.
inline void serialize_to_json(const vec3f& vals, json& js) {
    serialize_to_json((const array<float, 3>&)vals, js);
}

// Dump support function.
inline void serialize_to_json(const vec4f& vals, json& js) {
    serialize_to_json((const array<float, 4>&)vals, js);
}

// Dump support function.
inline void serialize_to_json(const quat4f& vals, json& js) {
    serialize_to_json((const array<float, 4>&)vals, js);
}

// Dump support function.
inline void serialize_to_json(const mat4f& vals, json& js) {
    serialize_to_json((const array<float, 16>&)vals, js);
}

// Dump support function.
inline void serialize_to_json(const frame3f& vals, json& js) {
    serialize_to_json((const array<float, 12>&)vals, js);
}

// Dump support function.
inline void serialize_to_json_obj(json& js) {
    if (!js.is_object()) js = json::object();
}

// Dump support function.
template <typename T>
inline void serialize_to_json_attr(const T& val, json& js, const char* name,
    bool required = true, const T& def = {}) {
    if (required || val != def) serialize_to_json(val, js[name]);
}

// Dump support function.
template <typename T>
inline void serialize_to_json_attr(const vector<T>& val, json& js,
    const char* name, bool required = true, const vector<T>& def = {}) {
    if (required || !val.empty()) serialize_to_json(val, js[name]);
}

// #codegen begin func ---------------------------------------------------------

// Parse id function.
template <typename T>
inline void serialize_from_json(glTFid<T>& val, const json& js) {
    if (!js.is_number_integer()) throw runtime_error("int expected");
    val = glTFid<T>((int)js);
}

// Parses a glTFProperty object
inline void serialize_from_json(glTFProperty& val, const json& js) {
    if (!js.is_object()) throw runtime_error("object expected");
#if YGL_GLTFJSON
    if (js.count("extensions"))
        serialize_from_json(val.extensions, js.at("extensions"));
    if (js.count("extras")) serialize_from_json(val.extras, js.at("extras"));
#endif
}

// Parses a glTFChildOfRootProperty object
inline void serialize_from_json(glTFChildOfRootProperty& val, const json& js) {
    static auto def = glTFChildOfRootProperty();
    serialize_from_json_obj(js);
    serialize_from_json((glTFProperty&)val, js);
    serialize_from_json_attr(val.name, js, "name", false, def.name);
}
// Parse a glTFAccessorSparseIndicesComponentType enum
inline void serialize_from_json(
    glTFAccessorSparseIndicesComponentType& val, const json& js) {
    static vector<pair<int, glTFAccessorSparseIndicesComponentType>> table = {
        {5121, glTFAccessorSparseIndicesComponentType::UnsignedByte},
        {5123, glTFAccessorSparseIndicesComponentType::UnsignedShort},
        {5125, glTFAccessorSparseIndicesComponentType::UnsignedInt},
    };
    serialize_from_json(val, js, table);
}

// Parses a glTFAccessorSparseIndices object
inline void serialize_from_json(
    glTFAccessorSparseIndices& val, const json& js) {
    static auto def = glTFAccessorSparseIndices();
    serialize_from_json_obj(js);
    serialize_from_json((glTFProperty&)val, js);
    serialize_from_json_attr(
        val.bufferView, js, "bufferView", true, def.bufferView);
    serialize_from_json_attr(
        val.byteOffset, js, "byteOffset", false, def.byteOffset);
    serialize_from_json_attr(
        val.componentType, js, "componentType", true, def.componentType);
}

// Parses a glTFAccessorSparseValues object
inline void serialize_from_json(glTFAccessorSparseValues& val, const json& js) {
    static auto def = glTFAccessorSparseValues();
    serialize_from_json_obj(js);
    serialize_from_json((glTFProperty&)val, js);
    serialize_from_json_attr(
        val.bufferView, js, "bufferView", true, def.bufferView);
    serialize_from_json_attr(
        val.byteOffset, js, "byteOffset", false, def.byteOffset);
}

// Parses a glTFAccessorSparse object
inline void serialize_from_json(glTFAccessorSparse& val, const json& js) {
    static auto def = glTFAccessorSparse();
    serialize_from_json_obj(js);
    serialize_from_json((glTFProperty&)val, js);
    serialize_from_json_attr(val.count, js, "count", true, def.count);
    serialize_from_json_attr(val.indices, js, "indices", true, def.indices);
    serialize_from_json_attr(val.values, js, "values", true, def.values);
}
// Parse a glTFAccessorComponentType enum
inline void serialize_from_json(
    glTFAccessorComponentType& val, const json& js) {
    static vector<pair<int, glTFAccessorComponentType>> table = {
        {5120, glTFAccessorComponentType::Byte},
        {5121, glTFAccessorComponentType::UnsignedByte},
        {5122, glTFAccessorComponentType::Short},
        {5123, glTFAccessorComponentType::UnsignedShort},
        {5125, glTFAccessorComponentType::UnsignedInt},
        {5126, glTFAccessorComponentType::Float},
    };
    serialize_from_json(val, js, table);
}

// Parse a glTFAccessorType enum
inline void serialize_from_json(glTFAccessorType& val, const json& js) {
    static vector<pair<string, glTFAccessorType>> table = {
        {"SCALAR", glTFAccessorType::Scalar},
        {"VEC2", glTFAccessorType::Vec2},
        {"VEC3", glTFAccessorType::Vec3},
        {"VEC4", glTFAccessorType::Vec4},
        {"MAT2", glTFAccessorType::Mat2},
        {"MAT3", glTFAccessorType::Mat3},
        {"MAT4", glTFAccessorType::Mat4},
    };
    serialize_from_json(val, js, table);
}

// Parses a glTFAccessor object
inline void serialize_from_json(glTFAccessor& val, const json& js) {
    static auto def = glTFAccessor();
    serialize_from_json_obj(js);
    serialize_from_json((glTFChildOfRootProperty&)val, js);
    serialize_from_json_attr(
        val.bufferView, js, "bufferView", false, def.bufferView);
    serialize_from_json_attr(
        val.byteOffset, js, "byteOffset", false, def.byteOffset);
    serialize_from_json_attr(
        val.componentType, js, "componentType", true, def.componentType);
    serialize_from_json_attr(
        val.normalized, js, "normalized", false, def.normalized);
    serialize_from_json_attr(val.count, js, "count", true, def.count);
    serialize_from_json_attr(val.type, js, "type", true, def.type);
    serialize_from_json_attr(val.max, js, "max", false, def.max);
    serialize_from_json_attr(val.min, js, "min", false, def.min);
    serialize_from_json_attr(val.sparse, js, "sparse", false, def.sparse);
}
// Parse a glTFAnimationChannelTargetPath enum
inline void serialize_from_json(
    glTFAnimationChannelTargetPath& val, const json& js) {
    static vector<pair<string, glTFAnimationChannelTargetPath>> table = {
        {"translation", glTFAnimationChannelTargetPath::Translation},
        {"rotation", glTFAnimationChannelTargetPath::Rotation},
        {"scale", glTFAnimationChannelTargetPath::Scale},
        {"weights", glTFAnimationChannelTargetPath::Weights},
    };
    serialize_from_json(val, js, table);
}

// Parses a glTFAnimationChannelTarget object
inline void serialize_from_json(
    glTFAnimationChannelTarget& val, const json& js) {
    static auto def = glTFAnimationChannelTarget();
    serialize_from_json_obj(js);
    serialize_from_json((glTFProperty&)val, js);
    serialize_from_json_attr(val.node, js, "node", true, def.node);
    serialize_from_json_attr(val.path, js, "path", true, def.path);
}

// Parses a glTFAnimationChannel object
inline void serialize_from_json(glTFAnimationChannel& val, const json& js) {
    static auto def = glTFAnimationChannel();
    serialize_from_json_obj(js);
    serialize_from_json((glTFProperty&)val, js);
    serialize_from_json_attr(val.sampler, js, "sampler", true, def.sampler);
    serialize_from_json_attr(val.target, js, "target", true, def.target);
}
// Parse a glTFAnimationSamplerInterpolation enum
inline void serialize_from_json(
    glTFAnimationSamplerInterpolation& val, const json& js) {
    static vector<pair<string, glTFAnimationSamplerInterpolation>> table = {
        {"LINEAR", glTFAnimationSamplerInterpolation::Linear},
        {"STEP", glTFAnimationSamplerInterpolation::Step},
        {"CATMULLROMSPLINE",
            glTFAnimationSamplerInterpolation::CatmullRomSpline},
        {"CUBICSPLINE", glTFAnimationSamplerInterpolation::CubicSpline},
    };
    serialize_from_json(val, js, table);
}

// Parses a glTFAnimationSampler object
inline void serialize_from_json(glTFAnimationSampler& val, const json& js) {
    static auto def = glTFAnimationSampler();
    serialize_from_json_obj(js);
    serialize_from_json((glTFProperty&)val, js);
    serialize_from_json_attr(val.input, js, "input", true, def.input);
    serialize_from_json_attr(
        val.interpolation, js, "interpolation", false, def.interpolation);
    serialize_from_json_attr(val.output, js, "output", true, def.output);
}

// Parses a glTFAnimation object
inline void serialize_from_json(glTFAnimation& val, const json& js) {
    static auto def = glTFAnimation();
    serialize_from_json_obj(js);
    serialize_from_json((glTFChildOfRootProperty&)val, js);
    serialize_from_json_attr(val.channels, js, "channels", true, def.channels);
    serialize_from_json_attr(val.samplers, js, "samplers", true, def.samplers);
}

// Parses a glTFAsset object
inline void serialize_from_json(glTFAsset& val, const json& js) {
    static auto def = glTFAsset();
    serialize_from_json_obj(js);
    serialize_from_json((glTFProperty&)val, js);
    serialize_from_json_attr(
        val.copyright, js, "copyright", false, def.copyright);
    serialize_from_json_attr(
        val.generator, js, "generator", false, def.generator);
    serialize_from_json_attr(val.version, js, "version", true, def.version);
    serialize_from_json_attr(
        val.minVersion, js, "minVersion", false, def.minVersion);
}

// Parses a glTFBuffer object
inline void serialize_from_json(glTFBuffer& val, const json& js) {
    static auto def = glTFBuffer();
    serialize_from_json_obj(js);
    serialize_from_json((glTFChildOfRootProperty&)val, js);
    serialize_from_json_attr(val.uri, js, "uri", false, def.uri);
    serialize_from_json_attr(
        val.byteLength, js, "byteLength", true, def.byteLength);
}
// Parse a glTFBufferViewTarget enum
inline void serialize_from_json(glTFBufferViewTarget& val, const json& js) {
    static vector<pair<int, glTFBufferViewTarget>> table = {
        {34962, glTFBufferViewTarget::ArrayBuffer},
        {34963, glTFBufferViewTarget::ElementArrayBuffer},
    };
    serialize_from_json(val, js, table);
}

// Parses a glTFBufferView object
inline void serialize_from_json(glTFBufferView& val, const json& js) {
    static auto def = glTFBufferView();
    serialize_from_json_obj(js);
    serialize_from_json((glTFChildOfRootProperty&)val, js);
    serialize_from_json_attr(val.buffer, js, "buffer", true, def.buffer);
    serialize_from_json_attr(
        val.byteOffset, js, "byteOffset", false, def.byteOffset);
    serialize_from_json_attr(
        val.byteLength, js, "byteLength", true, def.byteLength);
    serialize_from_json_attr(
        val.byteStride, js, "byteStride", false, def.byteStride);
    serialize_from_json_attr(val.target, js, "target", false, def.target);
}

// Parses a glTFCameraOrthographic object
inline void serialize_from_json(glTFCameraOrthographic& val, const json& js) {
    static auto def = glTFCameraOrthographic();
    serialize_from_json_obj(js);
    serialize_from_json((glTFProperty&)val, js);
    serialize_from_json_attr(val.xmag, js, "xmag", true, def.xmag);
    serialize_from_json_attr(val.ymag, js, "ymag", true, def.ymag);
    serialize_from_json_attr(val.zfar, js, "zfar", true, def.zfar);
    serialize_from_json_attr(val.znear, js, "znear", true, def.znear);
}

// Parses a glTFCameraPerspective object
inline void serialize_from_json(glTFCameraPerspective& val, const json& js) {
    static auto def = glTFCameraPerspective();
    serialize_from_json_obj(js);
    serialize_from_json((glTFProperty&)val, js);
    serialize_from_json_attr(
        val.aspectRatio, js, "aspectRatio", false, def.aspectRatio);
    serialize_from_json_attr(val.yfov, js, "yfov", true, def.yfov);
    serialize_from_json_attr(val.zfar, js, "zfar", false, def.zfar);
    serialize_from_json_attr(val.znear, js, "znear", true, def.znear);
}
// Parse a glTFCameraType enum
inline void serialize_from_json(glTFCameraType& val, const json& js) {
    static vector<pair<string, glTFCameraType>> table = {
        {"perspective", glTFCameraType::Perspective},
        {"orthographic", glTFCameraType::Orthographic},
    };
    serialize_from_json(val, js, table);
}

// Parses a glTFCamera object
inline void serialize_from_json(glTFCamera& val, const json& js) {
    static auto def = glTFCamera();
    serialize_from_json_obj(js);
    serialize_from_json((glTFChildOfRootProperty&)val, js);
    serialize_from_json_attr(
        val.orthographic, js, "orthographic", false, def.orthographic);
    serialize_from_json_attr(
        val.perspective, js, "perspective", false, def.perspective);
    serialize_from_json_attr(val.type, js, "type", true, def.type);
}
// Parse a glTFImageMimeType enum
inline void serialize_from_json(glTFImageMimeType& val, const json& js) {
    static vector<pair<string, glTFImageMimeType>> table = {
        {"image/jpeg", glTFImageMimeType::ImageJpeg},
        {"image/png", glTFImageMimeType::ImagePng},
    };
    serialize_from_json(val, js, table);
}

// Parses a glTFImage object
inline void serialize_from_json(glTFImage& val, const json& js) {
    static auto def = glTFImage();
    serialize_from_json_obj(js);
    serialize_from_json((glTFChildOfRootProperty&)val, js);
    serialize_from_json_attr(val.uri, js, "uri", false, def.uri);
    serialize_from_json_attr(val.mimeType, js, "mimeType", false, def.mimeType);
    serialize_from_json_attr(
        val.bufferView, js, "bufferView", false, def.bufferView);
}

// Parses a glTFTextureInfo object
inline void serialize_from_json(glTFTextureInfo& val, const json& js) {
    static auto def = glTFTextureInfo();
    serialize_from_json_obj(js);
    serialize_from_json((glTFProperty&)val, js);
    serialize_from_json_attr(val.index, js, "index", true, def.index);
    serialize_from_json_attr(val.texCoord, js, "texCoord", false, def.texCoord);
}

// Parses a glTFTexture object
inline void serialize_from_json(glTFTexture& val, const json& js) {
    static auto def = glTFTexture();
    serialize_from_json_obj(js);
    serialize_from_json((glTFChildOfRootProperty&)val, js);
    serialize_from_json_attr(val.sampler, js, "sampler", false, def.sampler);
    serialize_from_json_attr(val.source, js, "source", false, def.source);
}

// Parses a glTFMaterialNormalTextureInfo object
inline void serialize_from_json(
    glTFMaterialNormalTextureInfo& val, const json& js) {
    static auto def = glTFMaterialNormalTextureInfo();
    serialize_from_json_obj(js);
    serialize_from_json((glTFTextureInfo&)val, js);
    serialize_from_json_attr(val.scale, js, "scale", false, def.scale);
}

// Parses a glTFMaterialOcclusionTextureInfo object
inline void serialize_from_json(
    glTFMaterialOcclusionTextureInfo& val, const json& js) {
    static auto def = glTFMaterialOcclusionTextureInfo();
    serialize_from_json_obj(js);
    serialize_from_json((glTFTextureInfo&)val, js);
    serialize_from_json_attr(val.strength, js, "strength", false, def.strength);
}

// Parses a glTFMaterialPbrMetallicRoughness object
inline void serialize_from_json(
    glTFMaterialPbrMetallicRoughness& val, const json& js) {
    static auto def = glTFMaterialPbrMetallicRoughness();
    serialize_from_json_obj(js);
    serialize_from_json((glTFProperty&)val, js);
    serialize_from_json_attr(
        val.baseColorFactor, js, "baseColorFactor", false, def.baseColorFactor);
    serialize_from_json_attr(val.baseColorTexture, js, "baseColorTexture",
        false, def.baseColorTexture);
    serialize_from_json_attr(
        val.metallicFactor, js, "metallicFactor", false, def.metallicFactor);
    serialize_from_json_attr(
        val.roughnessFactor, js, "roughnessFactor", false, def.roughnessFactor);
    serialize_from_json_attr(val.metallicRoughnessTexture, js,
        "metallicRoughnessTexture", false, def.metallicRoughnessTexture);
}

// Parses a glTFMaterialPbrSpecularGlossiness object
inline void serialize_from_json(
    glTFMaterialPbrSpecularGlossiness& val, const json& js) {
    static auto def = glTFMaterialPbrSpecularGlossiness();
    serialize_from_json_obj(js);
    serialize_from_json((glTFProperty&)val, js);
    serialize_from_json_attr(
        val.diffuseFactor, js, "diffuseFactor", false, def.diffuseFactor);
    serialize_from_json_attr(
        val.diffuseTexture, js, "diffuseTexture", false, def.diffuseTexture);
    serialize_from_json_attr(
        val.specularFactor, js, "specularFactor", false, def.specularFactor);
    serialize_from_json_attr(val.glossinessFactor, js, "glossinessFactor",
        false, def.glossinessFactor);
    serialize_from_json_attr(val.specularGlossinessTexture, js,
        "specularGlossinessTexture", false, def.specularGlossinessTexture);
}
// Parse a glTFMaterialAlphaMode enum
inline void serialize_from_json(glTFMaterialAlphaMode& val, const json& js) {
    static vector<pair<string, glTFMaterialAlphaMode>> table = {
        {"OPAQUE", glTFMaterialAlphaMode::Opaque},
        {"MASK", glTFMaterialAlphaMode::Mask},
        {"BLEND", glTFMaterialAlphaMode::Blend},
    };
    serialize_from_json(val, js, table);
}

// Parses a glTFMaterial object
inline void serialize_from_json(glTFMaterial& val, const json& js) {
    static auto def = glTFMaterial();
    serialize_from_json_obj(js);
    serialize_from_json((glTFChildOfRootProperty&)val, js);
    serialize_from_json_attr(val.pbrMetallicRoughness, js,
        "pbrMetallicRoughness", false, def.pbrMetallicRoughness);
    serialize_from_json_attr(
        val.normalTexture, js, "normalTexture", false, def.normalTexture);
    serialize_from_json_attr(val.occlusionTexture, js, "occlusionTexture",
        false, def.occlusionTexture);
    serialize_from_json_attr(
        val.emissiveTexture, js, "emissiveTexture", false, def.emissiveTexture);
    serialize_from_json_attr(
        val.emissiveFactor, js, "emissiveFactor", false, def.emissiveFactor);
    serialize_from_json_attr(
        val.alphaMode, js, "alphaMode", false, def.alphaMode);
    serialize_from_json_attr(
        val.alphaCutoff, js, "alphaCutoff", false, def.alphaCutoff);
    serialize_from_json_attr(
        val.doubleSided, js, "doubleSided", false, def.doubleSided);
    if (js.count("extensions")) {
        auto& js_ext = js["extensions"];
        serialize_from_json_attr(val.pbrSpecularGlossiness, js_ext,
            "KHR_materials_pbrSpecularGlossiness", false,
            def.pbrSpecularGlossiness);
    }
}
// Parse a glTFMeshPrimitiveMode enum
inline void serialize_from_json(glTFMeshPrimitiveMode& val, const json& js) {
    static vector<pair<int, glTFMeshPrimitiveMode>> table = {
        {0, glTFMeshPrimitiveMode::Points},
        {1, glTFMeshPrimitiveMode::Lines},
        {2, glTFMeshPrimitiveMode::LineLoop},
        {3, glTFMeshPrimitiveMode::LineStrip},
        {4, glTFMeshPrimitiveMode::Triangles},
        {5, glTFMeshPrimitiveMode::TriangleStrip},
        {6, glTFMeshPrimitiveMode::TriangleFan},
    };
    serialize_from_json(val, js, table);
}

// Parses a glTFMeshPrimitive object
inline void serialize_from_json(glTFMeshPrimitive& val, const json& js) {
    static auto def = glTFMeshPrimitive();
    serialize_from_json_obj(js);
    serialize_from_json((glTFProperty&)val, js);
    serialize_from_json_attr(
        val.attributes, js, "attributes", true, def.attributes);
    serialize_from_json_attr(val.indices, js, "indices", false, def.indices);
    serialize_from_json_attr(val.material, js, "material", false, def.material);
    serialize_from_json_attr(val.mode, js, "mode", false, def.mode);
    serialize_from_json_attr(val.targets, js, "targets", false, def.targets);
}

// Parses a glTFMesh object
inline void serialize_from_json(glTFMesh& val, const json& js) {
    static auto def = glTFMesh();
    serialize_from_json_obj(js);
    serialize_from_json((glTFChildOfRootProperty&)val, js);
    serialize_from_json_attr(
        val.primitives, js, "primitives", true, def.primitives);
    serialize_from_json_attr(val.weights, js, "weights", false, def.weights);
}

// Parses a glTFNode object
inline void serialize_from_json(glTFNode& val, const json& js) {
    static auto def = glTFNode();
    serialize_from_json_obj(js);
    serialize_from_json((glTFChildOfRootProperty&)val, js);
    serialize_from_json_attr(val.camera, js, "camera", false, def.camera);
    serialize_from_json_attr(val.children, js, "children", false, def.children);
    serialize_from_json_attr(val.skin, js, "skin", false, def.skin);
    serialize_from_json_attr(val.matrix, js, "matrix", false, def.matrix);
    serialize_from_json_attr(val.mesh, js, "mesh", false, def.mesh);
    serialize_from_json_attr(val.rotation, js, "rotation", false, def.rotation);
    serialize_from_json_attr(val.scale, js, "scale", false, def.scale);
    serialize_from_json_attr(
        val.translation, js, "translation", false, def.translation);
    serialize_from_json_attr(val.weights, js, "weights", false, def.weights);
}
// Parse a glTFSamplerMagFilter enum
inline void serialize_from_json(glTFSamplerMagFilter& val, const json& js) {
    static vector<pair<int, glTFSamplerMagFilter>> table = {
        {9728, glTFSamplerMagFilter::Nearest},
        {9729, glTFSamplerMagFilter::Linear},
    };
    serialize_from_json(val, js, table);
}

// Parse a glTFSamplerMinFilter enum
inline void serialize_from_json(glTFSamplerMinFilter& val, const json& js) {
    static vector<pair<int, glTFSamplerMinFilter>> table = {
        {9728, glTFSamplerMinFilter::Nearest},
        {9729, glTFSamplerMinFilter::Linear},
        {9984, glTFSamplerMinFilter::NearestMipmapNearest},
        {9985, glTFSamplerMinFilter::LinearMipmapNearest},
        {9986, glTFSamplerMinFilter::NearestMipmapLinear},
        {9987, glTFSamplerMinFilter::LinearMipmapLinear},
    };
    serialize_from_json(val, js, table);
}

// Parse a glTFSamplerWrapS enum
inline void serialize_from_json(glTFSamplerWrapS& val, const json& js) {
    static vector<pair<int, glTFSamplerWrapS>> table = {
        {33071, glTFSamplerWrapS::ClampToEdge},
        {33648, glTFSamplerWrapS::MirroredRepeat},
        {10497, glTFSamplerWrapS::Repeat},
    };
    serialize_from_json(val, js, table);
}

// Parse a glTFSamplerWrapT enum
inline void serialize_from_json(glTFSamplerWrapT& val, const json& js) {
    static vector<pair<int, glTFSamplerWrapT>> table = {
        {33071, glTFSamplerWrapT::ClampToEdge},
        {33648, glTFSamplerWrapT::MirroredRepeat},
        {10497, glTFSamplerWrapT::Repeat},
    };
    serialize_from_json(val, js, table);
}

// Parses a glTFSampler object
inline void serialize_from_json(glTFSampler& val, const json& js) {
    static auto def = glTFSampler();
    serialize_from_json_obj(js);
    serialize_from_json((glTFChildOfRootProperty&)val, js);
    serialize_from_json_attr(
        val.magFilter, js, "magFilter", false, def.magFilter);
    serialize_from_json_attr(
        val.minFilter, js, "minFilter", false, def.minFilter);
    serialize_from_json_attr(val.wrapS, js, "wrapS", false, def.wrapS);
    serialize_from_json_attr(val.wrapT, js, "wrapT", false, def.wrapT);
}

// Parses a glTFScene object
inline void serialize_from_json(glTFScene& val, const json& js) {
    static auto def = glTFScene();
    serialize_from_json_obj(js);
    serialize_from_json((glTFChildOfRootProperty&)val, js);
    serialize_from_json_attr(val.nodes, js, "nodes", false, def.nodes);
}

// Parses a glTFSkin object
inline void serialize_from_json(glTFSkin& val, const json& js) {
    static auto def = glTFSkin();
    serialize_from_json_obj(js);
    serialize_from_json((glTFChildOfRootProperty&)val, js);
    serialize_from_json_attr(val.inverseBindMatrices, js, "inverseBindMatrices",
        false, def.inverseBindMatrices);
    serialize_from_json_attr(val.skeleton, js, "skeleton", false, def.skeleton);
    serialize_from_json_attr(val.joints, js, "joints", true, def.joints);
}

// Parses a glTF object
inline void serialize_from_json(glTF& val, const json& js) {
    static auto def = glTF();
    serialize_from_json_obj(js);
    serialize_from_json((glTFProperty&)val, js);
    serialize_from_json_attr(
        val.extensionsUsed, js, "extensionsUsed", false, def.extensionsUsed);
    serialize_from_json_attr(val.extensionsRequired, js, "extensionsRequired",
        false, def.extensionsRequired);
    serialize_from_json_attr(
        val.accessors, js, "accessors", false, def.accessors);
    serialize_from_json_attr(
        val.animations, js, "animations", false, def.animations);
    serialize_from_json_attr(val.asset, js, "asset", true, def.asset);
    serialize_from_json_attr(val.buffers, js, "buffers", false, def.buffers);
    serialize_from_json_attr(
        val.bufferViews, js, "bufferViews", false, def.bufferViews);
    serialize_from_json_attr(val.cameras, js, "cameras", false, def.cameras);
    serialize_from_json_attr(val.images, js, "images", false, def.images);
    serialize_from_json_attr(
        val.materials, js, "materials", false, def.materials);
    serialize_from_json_attr(val.meshes, js, "meshes", false, def.meshes);
    serialize_from_json_attr(val.nodes, js, "nodes", false, def.nodes);
    serialize_from_json_attr(val.samplers, js, "samplers", false, def.samplers);
    serialize_from_json_attr(val.scene, js, "scene", false, def.scene);
    serialize_from_json_attr(val.scenes, js, "scenes", false, def.scenes);
    serialize_from_json_attr(val.skins, js, "skins", false, def.skins);
    serialize_from_json_attr(val.textures, js, "textures", false, def.textures);
}

// Converts glTFid to json.
template <typename T>
inline void serialize_to_json(const glTFid<T>& val, json& js) {
    js = (int)val;
}

// Check for default value
template <typename T>
inline bool operator==(const glTFid<T>& a, const glTFid<T>& b) {
    return (int)a == (int)b;
}

// Check for default value
template <typename T>
inline bool operator!=(const glTFid<T>& a, const glTFid<T>& b) {
    return (int)a != (int)b;
}

// Converts a glTFProperty object to JSON
inline void serialize_to_json(const glTFProperty& val, json& js) {
    if (!js.is_object()) js = json::object();
#if YGL_GLTFJSON
    if (!val.extensions.empty())
        serialize_to_json(val.extensions, js["extensions"]);
    if (!val.extras.is_null()) dump_attr(val.extras, "extras", js);
#endif
}

// Converts a glTFChildOfRootProperty object to JSON
inline void serialize_to_json(const glTFChildOfRootProperty& val, json& js) {
    static auto def = glTFChildOfRootProperty();
    serialize_to_json_obj(js);
    serialize_to_json((const glTFProperty&)val, js);
    serialize_to_json_attr(val.name, js, "name", false, def.name);
}
// Converts a glTFAccessorSparseIndicesComponentType enum to JSON
inline void serialize_to_json(
    const glTFAccessorSparseIndicesComponentType& val, json& js) {
    static vector<pair<int, glTFAccessorSparseIndicesComponentType>> table = {
        {5121, glTFAccessorSparseIndicesComponentType::UnsignedByte},
        {5123, glTFAccessorSparseIndicesComponentType::UnsignedShort},
        {5125, glTFAccessorSparseIndicesComponentType::UnsignedInt},
    };
    serialize_to_json(val, js, table);
}

// Converts a glTFAccessorSparseIndices object to JSON
inline void serialize_to_json(const glTFAccessorSparseIndices& val, json& js) {
    static auto def = glTFAccessorSparseIndices();
    serialize_to_json_obj(js);
    serialize_to_json((const glTFProperty&)val, js);
    serialize_to_json_attr(
        val.bufferView, js, "bufferView", true, def.bufferView);
    serialize_to_json_attr(
        val.byteOffset, js, "byteOffset", false, def.byteOffset);
    serialize_to_json_attr(
        val.componentType, js, "componentType", true, def.componentType);
}

// Converts a glTFAccessorSparseValues object to JSON
inline void serialize_to_json(const glTFAccessorSparseValues& val, json& js) {
    static auto def = glTFAccessorSparseValues();
    serialize_to_json_obj(js);
    serialize_to_json((const glTFProperty&)val, js);
    serialize_to_json_attr(
        val.bufferView, js, "bufferView", true, def.bufferView);
    serialize_to_json_attr(
        val.byteOffset, js, "byteOffset", false, def.byteOffset);
}

// Converts a glTFAccessorSparse object to JSON
inline void serialize_to_json(const glTFAccessorSparse& val, json& js) {
    static auto def = glTFAccessorSparse();
    serialize_to_json_obj(js);
    serialize_to_json((const glTFProperty&)val, js);
    serialize_to_json_attr(val.count, js, "count", true, def.count);
    serialize_to_json_attr(val.indices, js, "indices", true, def.indices);
    serialize_to_json_attr(val.values, js, "values", true, def.values);
}
// Converts a glTFAccessorComponentType enum to JSON
inline void serialize_to_json(const glTFAccessorComponentType& val, json& js) {
    static vector<pair<int, glTFAccessorComponentType>> table = {
        {5120, glTFAccessorComponentType::Byte},
        {5121, glTFAccessorComponentType::UnsignedByte},
        {5122, glTFAccessorComponentType::Short},
        {5123, glTFAccessorComponentType::UnsignedShort},
        {5125, glTFAccessorComponentType::UnsignedInt},
        {5126, glTFAccessorComponentType::Float},
    };
    serialize_to_json(val, js, table);
}

// Converts a glTFAccessorType enum to JSON
inline void serialize_to_json(const glTFAccessorType& val, json& js) {
    static vector<pair<string, glTFAccessorType>> table = {
        {"SCALAR", glTFAccessorType::Scalar},
        {"VEC2", glTFAccessorType::Vec2},
        {"VEC3", glTFAccessorType::Vec3},
        {"VEC4", glTFAccessorType::Vec4},
        {"MAT2", glTFAccessorType::Mat2},
        {"MAT3", glTFAccessorType::Mat3},
        {"MAT4", glTFAccessorType::Mat4},
    };
    serialize_to_json(val, js, table);
}

// Converts a glTFAccessor object to JSON
inline void serialize_to_json(const glTFAccessor& val, json& js) {
    static auto def = glTFAccessor();
    serialize_to_json_obj(js);
    serialize_to_json((const glTFChildOfRootProperty&)val, js);
    serialize_to_json_attr(
        val.bufferView, js, "bufferView", false, def.bufferView);
    serialize_to_json_attr(
        val.byteOffset, js, "byteOffset", false, def.byteOffset);
    serialize_to_json_attr(
        val.componentType, js, "componentType", true, def.componentType);
    serialize_to_json_attr(
        val.normalized, js, "normalized", false, def.normalized);
    serialize_to_json_attr(val.count, js, "count", true, def.count);
    serialize_to_json_attr(val.type, js, "type", true, def.type);
    serialize_to_json_attr(val.max, js, "max", false, def.max);
    serialize_to_json_attr(val.min, js, "min", false, def.min);
    serialize_to_json_attr(val.sparse, js, "sparse", false, def.sparse);
}
// Converts a glTFAnimationChannelTargetPath enum to JSON
inline void serialize_to_json(
    const glTFAnimationChannelTargetPath& val, json& js) {
    static vector<pair<string, glTFAnimationChannelTargetPath>> table = {
        {"translation", glTFAnimationChannelTargetPath::Translation},
        {"rotation", glTFAnimationChannelTargetPath::Rotation},
        {"scale", glTFAnimationChannelTargetPath::Scale},
        {"weights", glTFAnimationChannelTargetPath::Weights},
    };
    serialize_to_json(val, js, table);
}

// Converts a glTFAnimationChannelTarget object to JSON
inline void serialize_to_json(const glTFAnimationChannelTarget& val, json& js) {
    static auto def = glTFAnimationChannelTarget();
    serialize_to_json_obj(js);
    serialize_to_json((const glTFProperty&)val, js);
    serialize_to_json_attr(val.node, js, "node", true, def.node);
    serialize_to_json_attr(val.path, js, "path", true, def.path);
}

// Converts a glTFAnimationChannel object to JSON
inline void serialize_to_json(const glTFAnimationChannel& val, json& js) {
    static auto def = glTFAnimationChannel();
    serialize_to_json_obj(js);
    serialize_to_json((const glTFProperty&)val, js);
    serialize_to_json_attr(val.sampler, js, "sampler", true, def.sampler);
    serialize_to_json_attr(val.target, js, "target", true, def.target);
}
// Converts a glTFAnimationSamplerInterpolation enum to JSON
inline void serialize_to_json(
    const glTFAnimationSamplerInterpolation& val, json& js) {
    static vector<pair<string, glTFAnimationSamplerInterpolation>> table = {
        {"LINEAR", glTFAnimationSamplerInterpolation::Linear},
        {"STEP", glTFAnimationSamplerInterpolation::Step},
        {"CATMULLROMSPLINE",
            glTFAnimationSamplerInterpolation::CatmullRomSpline},
        {"CUBICSPLINE", glTFAnimationSamplerInterpolation::CubicSpline},
    };
    serialize_to_json(val, js, table);
}

// Converts a glTFAnimationSampler object to JSON
inline void serialize_to_json(const glTFAnimationSampler& val, json& js) {
    static auto def = glTFAnimationSampler();
    serialize_to_json_obj(js);
    serialize_to_json((const glTFProperty&)val, js);
    serialize_to_json_attr(val.input, js, "input", true, def.input);
    serialize_to_json_attr(
        val.interpolation, js, "interpolation", false, def.interpolation);
    serialize_to_json_attr(val.output, js, "output", true, def.output);
}

// Converts a glTFAnimation object to JSON
inline void serialize_to_json(const glTFAnimation& val, json& js) {
    static auto def = glTFAnimation();
    serialize_to_json_obj(js);
    serialize_to_json((const glTFChildOfRootProperty&)val, js);
    serialize_to_json_attr(val.channels, js, "channels", true, def.channels);
    serialize_to_json_attr(val.samplers, js, "samplers", true, def.samplers);
}

// Converts a glTFAsset object to JSON
inline void serialize_to_json(const glTFAsset& val, json& js) {
    static auto def = glTFAsset();
    serialize_to_json_obj(js);
    serialize_to_json((const glTFProperty&)val, js);
    serialize_to_json_attr(
        val.copyright, js, "copyright", false, def.copyright);
    serialize_to_json_attr(
        val.generator, js, "generator", false, def.generator);
    serialize_to_json_attr(val.version, js, "version", true, def.version);
    serialize_to_json_attr(
        val.minVersion, js, "minVersion", false, def.minVersion);
}

// Converts a glTFBuffer object to JSON
inline void serialize_to_json(const glTFBuffer& val, json& js) {
    static auto def = glTFBuffer();
    serialize_to_json_obj(js);
    serialize_to_json((const glTFChildOfRootProperty&)val, js);
    serialize_to_json_attr(val.uri, js, "uri", false, def.uri);
    serialize_to_json_attr(
        val.byteLength, js, "byteLength", true, def.byteLength);
}
// Converts a glTFBufferViewTarget enum to JSON
inline void serialize_to_json(const glTFBufferViewTarget& val, json& js) {
    static vector<pair<int, glTFBufferViewTarget>> table = {
        {34962, glTFBufferViewTarget::ArrayBuffer},
        {34963, glTFBufferViewTarget::ElementArrayBuffer},
    };
    serialize_to_json(val, js, table);
}

// Converts a glTFBufferView object to JSON
inline void serialize_to_json(const glTFBufferView& val, json& js) {
    static auto def = glTFBufferView();
    serialize_to_json_obj(js);
    serialize_to_json((const glTFChildOfRootProperty&)val, js);
    serialize_to_json_attr(val.buffer, js, "buffer", true, def.buffer);
    serialize_to_json_attr(
        val.byteOffset, js, "byteOffset", false, def.byteOffset);
    serialize_to_json_attr(
        val.byteLength, js, "byteLength", true, def.byteLength);
    serialize_to_json_attr(
        val.byteStride, js, "byteStride", false, def.byteStride);
    serialize_to_json_attr(val.target, js, "target", false, def.target);
}

// Converts a glTFCameraOrthographic object to JSON
inline void serialize_to_json(const glTFCameraOrthographic& val, json& js) {
    static auto def = glTFCameraOrthographic();
    serialize_to_json_obj(js);
    serialize_to_json((const glTFProperty&)val, js);
    serialize_to_json_attr(val.xmag, js, "xmag", true, def.xmag);
    serialize_to_json_attr(val.ymag, js, "ymag", true, def.ymag);
    serialize_to_json_attr(val.zfar, js, "zfar", true, def.zfar);
    serialize_to_json_attr(val.znear, js, "znear", true, def.znear);
}

// Converts a glTFCameraPerspective object to JSON
inline void serialize_to_json(const glTFCameraPerspective& val, json& js) {
    static auto def = glTFCameraPerspective();
    serialize_to_json_obj(js);
    serialize_to_json((const glTFProperty&)val, js);
    serialize_to_json_attr(
        val.aspectRatio, js, "aspectRatio", false, def.aspectRatio);
    serialize_to_json_attr(val.yfov, js, "yfov", true, def.yfov);
    serialize_to_json_attr(val.zfar, js, "zfar", false, def.zfar);
    serialize_to_json_attr(val.znear, js, "znear", true, def.znear);
}
// Converts a glTFCameraType enum to JSON
inline void serialize_to_json(const glTFCameraType& val, json& js) {
    static vector<pair<string, glTFCameraType>> table = {
        {"perspective", glTFCameraType::Perspective},
        {"orthographic", glTFCameraType::Orthographic},
    };
    serialize_to_json(val, js, table);
}

// Converts a glTFCamera object to JSON
inline void serialize_to_json(const glTFCamera& val, json& js) {
    static auto def = glTFCamera();
    serialize_to_json_obj(js);
    serialize_to_json((const glTFChildOfRootProperty&)val, js);
    serialize_to_json_attr(
        val.orthographic, js, "orthographic", false, def.orthographic);
    serialize_to_json_attr(
        val.perspective, js, "perspective", false, def.perspective);
    serialize_to_json_attr(val.type, js, "type", true, def.type);
}
// Converts a glTFImageMimeType enum to JSON
inline void serialize_to_json(const glTFImageMimeType& val, json& js) {
    static vector<pair<string, glTFImageMimeType>> table = {
        {"image/jpeg", glTFImageMimeType::ImageJpeg},
        {"image/png", glTFImageMimeType::ImagePng},
    };
    serialize_to_json(val, js, table);
}

// Converts a glTFImage object to JSON
inline void serialize_to_json(const glTFImage& val, json& js) {
    static auto def = glTFImage();
    serialize_to_json_obj(js);
    serialize_to_json((const glTFChildOfRootProperty&)val, js);
    serialize_to_json_attr(val.uri, js, "uri", false, def.uri);
    serialize_to_json_attr(val.mimeType, js, "mimeType", false, def.mimeType);
    serialize_to_json_attr(
        val.bufferView, js, "bufferView", false, def.bufferView);
}

// Converts a glTFTextureInfo object to JSON
inline void serialize_to_json(const glTFTextureInfo& val, json& js) {
    static auto def = glTFTextureInfo();
    serialize_to_json_obj(js);
    serialize_to_json((const glTFProperty&)val, js);
    serialize_to_json_attr(val.index, js, "index", true, def.index);
    serialize_to_json_attr(val.texCoord, js, "texCoord", false, def.texCoord);
}

// Converts a glTFTexture object to JSON
inline void serialize_to_json(const glTFTexture& val, json& js) {
    static auto def = glTFTexture();
    serialize_to_json_obj(js);
    serialize_to_json((const glTFChildOfRootProperty&)val, js);
    serialize_to_json_attr(val.sampler, js, "sampler", false, def.sampler);
    serialize_to_json_attr(val.source, js, "source", false, def.source);
}

// Converts a glTFMaterialNormalTextureInfo object to JSON
inline void serialize_to_json(
    const glTFMaterialNormalTextureInfo& val, json& js) {
    static auto def = glTFMaterialNormalTextureInfo();
    serialize_to_json_obj(js);
    serialize_to_json((const glTFTextureInfo&)val, js);
    serialize_to_json_attr(val.scale, js, "scale", false, def.scale);
}

// Converts a glTFMaterialOcclusionTextureInfo object to JSON
inline void serialize_to_json(
    const glTFMaterialOcclusionTextureInfo& val, json& js) {
    static auto def = glTFMaterialOcclusionTextureInfo();
    serialize_to_json_obj(js);
    serialize_to_json((const glTFTextureInfo&)val, js);
    serialize_to_json_attr(val.strength, js, "strength", false, def.strength);
}

// Converts a glTFMaterialPbrMetallicRoughness object to JSON
inline void serialize_to_json(
    const glTFMaterialPbrMetallicRoughness& val, json& js) {
    static auto def = glTFMaterialPbrMetallicRoughness();
    serialize_to_json_obj(js);
    serialize_to_json((const glTFProperty&)val, js);
    serialize_to_json_attr(
        val.baseColorFactor, js, "baseColorFactor", false, def.baseColorFactor);
    serialize_to_json_attr(val.baseColorTexture, js, "baseColorTexture", false,
        def.baseColorTexture);
    serialize_to_json_attr(
        val.metallicFactor, js, "metallicFactor", false, def.metallicFactor);
    serialize_to_json_attr(
        val.roughnessFactor, js, "roughnessFactor", false, def.roughnessFactor);
    serialize_to_json_attr(val.metallicRoughnessTexture, js,
        "metallicRoughnessTexture", false, def.metallicRoughnessTexture);
}

// Converts a glTFMaterialPbrSpecularGlossiness object to JSON
inline void serialize_to_json(
    const glTFMaterialPbrSpecularGlossiness& val, json& js) {
    static auto def = glTFMaterialPbrSpecularGlossiness();
    serialize_to_json_obj(js);
    serialize_to_json((const glTFProperty&)val, js);
    serialize_to_json_attr(
        val.diffuseFactor, js, "diffuseFactor", false, def.diffuseFactor);
    serialize_to_json_attr(
        val.diffuseTexture, js, "diffuseTexture", false, def.diffuseTexture);
    serialize_to_json_attr(
        val.specularFactor, js, "specularFactor", false, def.specularFactor);
    serialize_to_json_attr(val.glossinessFactor, js, "glossinessFactor", false,
        def.glossinessFactor);
    serialize_to_json_attr(val.specularGlossinessTexture, js,
        "specularGlossinessTexture", false, def.specularGlossinessTexture);
}
// Converts a glTFMaterialAlphaMode enum to JSON
inline void serialize_to_json(const glTFMaterialAlphaMode& val, json& js) {
    static vector<pair<string, glTFMaterialAlphaMode>> table = {
        {"OPAQUE", glTFMaterialAlphaMode::Opaque},
        {"MASK", glTFMaterialAlphaMode::Mask},
        {"BLEND", glTFMaterialAlphaMode::Blend},
    };
    serialize_to_json(val, js, table);
}

// Converts a glTFMaterial object to JSON
inline void serialize_to_json(const glTFMaterial& val, json& js) {
    static auto def = glTFMaterial();
    serialize_to_json_obj(js);
    serialize_to_json((const glTFChildOfRootProperty&)val, js);
    serialize_to_json_attr(val.pbrMetallicRoughness, js, "pbrMetallicRoughness",
        false, def.pbrMetallicRoughness);
    serialize_to_json_attr(
        val.normalTexture, js, "normalTexture", false, def.normalTexture);
    serialize_to_json_attr(val.occlusionTexture, js, "occlusionTexture", false,
        def.occlusionTexture);
    serialize_to_json_attr(
        val.emissiveTexture, js, "emissiveTexture", false, def.emissiveTexture);
    serialize_to_json_attr(
        val.emissiveFactor, js, "emissiveFactor", false, def.emissiveFactor);
    serialize_to_json_attr(
        val.alphaMode, js, "alphaMode", false, def.alphaMode);
    serialize_to_json_attr(
        val.alphaCutoff, js, "alphaCutoff", false, def.alphaCutoff);
    serialize_to_json_attr(
        val.doubleSided, js, "doubleSided", false, def.doubleSided);

    if (val.pbrSpecularGlossiness != nullptr) {
        auto& js_ext = js["extensions"];
        serialize_to_json_attr(val.pbrSpecularGlossiness, js_ext,
            "KHR_materials_pbrSpecularGlossiness", false,
            def.pbrSpecularGlossiness);
    }
}
// Converts a glTFMeshPrimitiveMode enum to JSON
inline void serialize_to_json(const glTFMeshPrimitiveMode& val, json& js) {
    static vector<pair<int, glTFMeshPrimitiveMode>> table = {
        {0, glTFMeshPrimitiveMode::Points},
        {1, glTFMeshPrimitiveMode::Lines},
        {2, glTFMeshPrimitiveMode::LineLoop},
        {3, glTFMeshPrimitiveMode::LineStrip},
        {4, glTFMeshPrimitiveMode::Triangles},
        {5, glTFMeshPrimitiveMode::TriangleStrip},
        {6, glTFMeshPrimitiveMode::TriangleFan},
    };
    serialize_to_json(val, js, table);
}

// Converts a glTFMeshPrimitive object to JSON
inline void serialize_to_json(const glTFMeshPrimitive& val, json& js) {
    static auto def = glTFMeshPrimitive();
    serialize_to_json_obj(js);
    serialize_to_json((const glTFProperty&)val, js);
    serialize_to_json_attr(
        val.attributes, js, "attributes", true, def.attributes);
    serialize_to_json_attr(val.indices, js, "indices", false, def.indices);
    serialize_to_json_attr(val.material, js, "material", false, def.material);
    serialize_to_json_attr(val.mode, js, "mode", false, def.mode);
    serialize_to_json_attr(val.targets, js, "targets", false, def.targets);
}

// Converts a glTFMesh object to JSON
inline void serialize_to_json(const glTFMesh& val, json& js) {
    static auto def = glTFMesh();
    serialize_to_json_obj(js);
    serialize_to_json((const glTFChildOfRootProperty&)val, js);
    serialize_to_json_attr(
        val.primitives, js, "primitives", true, def.primitives);
    serialize_to_json_attr(val.weights, js, "weights", false, def.weights);
}

// Converts a glTFNode object to JSON
inline void serialize_to_json(const glTFNode& val, json& js) {
    static auto def = glTFNode();
    serialize_to_json_obj(js);
    serialize_to_json((const glTFChildOfRootProperty&)val, js);
    serialize_to_json_attr(val.camera, js, "camera", false, def.camera);
    serialize_to_json_attr(val.children, js, "children", false, def.children);
    serialize_to_json_attr(val.skin, js, "skin", false, def.skin);
    serialize_to_json_attr(val.matrix, js, "matrix", false, def.matrix);
    serialize_to_json_attr(val.mesh, js, "mesh", false, def.mesh);
    serialize_to_json_attr(val.rotation, js, "rotation", false, def.rotation);
    serialize_to_json_attr(val.scale, js, "scale", false, def.scale);
    serialize_to_json_attr(
        val.translation, js, "translation", false, def.translation);
    serialize_to_json_attr(val.weights, js, "weights", false, def.weights);
}
// Converts a glTFSamplerMagFilter enum to JSON
inline void serialize_to_json(const glTFSamplerMagFilter& val, json& js) {
    static vector<pair<int, glTFSamplerMagFilter>> table = {
        {9728, glTFSamplerMagFilter::Nearest},
        {9729, glTFSamplerMagFilter::Linear},
    };
    serialize_to_json(val, js, table);
}

// Converts a glTFSamplerMinFilter enum to JSON
inline void serialize_to_json(const glTFSamplerMinFilter& val, json& js) {
    static vector<pair<int, glTFSamplerMinFilter>> table = {
        {9728, glTFSamplerMinFilter::Nearest},
        {9729, glTFSamplerMinFilter::Linear},
        {9984, glTFSamplerMinFilter::NearestMipmapNearest},
        {9985, glTFSamplerMinFilter::LinearMipmapNearest},
        {9986, glTFSamplerMinFilter::NearestMipmapLinear},
        {9987, glTFSamplerMinFilter::LinearMipmapLinear},
    };
    serialize_to_json(val, js, table);
}

// Converts a glTFSamplerWrapS enum to JSON
inline void serialize_to_json(const glTFSamplerWrapS& val, json& js) {
    static vector<pair<int, glTFSamplerWrapS>> table = {
        {33071, glTFSamplerWrapS::ClampToEdge},
        {33648, glTFSamplerWrapS::MirroredRepeat},
        {10497, glTFSamplerWrapS::Repeat},
    };
    serialize_to_json(val, js, table);
}

// Converts a glTFSamplerWrapT enum to JSON
inline void serialize_to_json(const glTFSamplerWrapT& val, json& js) {
    static vector<pair<int, glTFSamplerWrapT>> table = {
        {33071, glTFSamplerWrapT::ClampToEdge},
        {33648, glTFSamplerWrapT::MirroredRepeat},
        {10497, glTFSamplerWrapT::Repeat},
    };
    serialize_to_json(val, js, table);
}

// Converts a glTFSampler object to JSON
inline void serialize_to_json(const glTFSampler& val, json& js) {
    static auto def = glTFSampler();
    serialize_to_json_obj(js);
    serialize_to_json((const glTFChildOfRootProperty&)val, js);
    serialize_to_json_attr(
        val.magFilter, js, "magFilter", false, def.magFilter);
    serialize_to_json_attr(
        val.minFilter, js, "minFilter", false, def.minFilter);
    serialize_to_json_attr(val.wrapS, js, "wrapS", false, def.wrapS);
    serialize_to_json_attr(val.wrapT, js, "wrapT", false, def.wrapT);
}

// Converts a glTFScene object to JSON
inline void serialize_to_json(const glTFScene& val, json& js) {
    static auto def = glTFScene();
    serialize_to_json_obj(js);
    serialize_to_json((const glTFChildOfRootProperty&)val, js);
    serialize_to_json_attr(val.nodes, js, "nodes", false, def.nodes);
}

// Converts a glTFSkin object to JSON
inline void serialize_to_json(const glTFSkin& val, json& js) {
    static auto def = glTFSkin();
    serialize_to_json_obj(js);
    serialize_to_json((const glTFChildOfRootProperty&)val, js);
    serialize_to_json_attr(val.inverseBindMatrices, js, "inverseBindMatrices",
        false, def.inverseBindMatrices);
    serialize_to_json_attr(val.skeleton, js, "skeleton", false, def.skeleton);
    serialize_to_json_attr(val.joints, js, "joints", true, def.joints);
}

// Converts a glTF object to JSON
inline void serialize_to_json(const glTF& val, json& js) {
    static auto def = glTF();
    serialize_to_json_obj(js);
    serialize_to_json((const glTFProperty&)val, js);
    serialize_to_json_attr(
        val.extensionsUsed, js, "extensionsUsed", false, def.extensionsUsed);
    serialize_to_json_attr(val.extensionsRequired, js, "extensionsRequired",
        false, def.extensionsRequired);
    serialize_to_json_attr(
        val.accessors, js, "accessors", false, def.accessors);
    serialize_to_json_attr(
        val.animations, js, "animations", false, def.animations);
    serialize_to_json_attr(val.asset, js, "asset", true, def.asset);
    serialize_to_json_attr(val.buffers, js, "buffers", false, def.buffers);
    serialize_to_json_attr(
        val.bufferViews, js, "bufferViews", false, def.bufferViews);
    serialize_to_json_attr(val.cameras, js, "cameras", false, def.cameras);
    serialize_to_json_attr(val.images, js, "images", false, def.images);
    serialize_to_json_attr(
        val.materials, js, "materials", false, def.materials);
    serialize_to_json_attr(val.meshes, js, "meshes", false, def.meshes);
    serialize_to_json_attr(val.nodes, js, "nodes", false, def.nodes);
    serialize_to_json_attr(val.samplers, js, "samplers", false, def.samplers);
    serialize_to_json_attr(val.scene, js, "scene", false, def.scene);
    serialize_to_json_attr(val.scenes, js, "scenes", false, def.scenes);
    serialize_to_json_attr(val.skins, js, "skins", false, def.skins);
    serialize_to_json_attr(val.textures, js, "textures", false, def.textures);
}
// #codegen end func

// Encode in base64
inline string base64_encode(
    unsigned char const* bytes_to_encode, unsigned int in_len) {
    static const string base64_chars =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        "abcdefghijklmnopqrstuvwxyz"
        "0123456789+/";

    string ret;
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
inline string base64_decode(string const& encoded_string) {
    static const string base64_chars =
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
    string ret;

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
void load_buffers(glTF* gltf, const string& dirname, bool skip_missing) {
    for (auto buffer : gltf->buffers) {
        if (buffer->uri == "") continue;
        try {
            if (startswith(buffer->uri, "data:")) {
                // assume it is base64 and find ','
                auto pos = buffer->uri.find(',');
                if (pos == buffer->uri.npos) {
                    if (skip_missing) continue;
                    throw runtime_error("could not decode base64 data");
                }
                // decode
                auto data = base64_decode(buffer->uri.substr(pos + 1));
                buffer->data =
                    vector<unsigned char>((unsigned char*)data.c_str(),
                        (unsigned char*)data.c_str() + data.length());
            } else {
                buffer->data =
                    load_binary(path_convert_eparator(dirname + buffer->uri));
                if (buffer->data.empty()) {
                    if (skip_missing) continue;
                    throw runtime_error(
                        "could not load binary file " +
                        path_convert_eparator(dirname + buffer->uri));
                }
            }
            if (buffer->byteLength != buffer->data.size()) {
                throw runtime_error("mismatched buffer size");
            }
        } catch (const std::exception&) {
            if (skip_missing) continue;
            throw;
        }
    }
}

// Loads images.
void load_images(glTF* gltf, const string& dirname, bool skip_missing) {
    for (auto image : gltf->images) {
        image->data = image_data();
        auto filename = string();
#if YGL_IMAGEIO
        if (image->bufferView || startswith(image->uri, "data:")) {
            auto buffer = string();
            auto data = (unsigned char*)nullptr;
            auto data_size = 0;
            if (image->bufferView) {
                auto view = gltf->get(image->bufferView);
                auto buffer = gltf->get(view->buffer);
                if (!view || !buffer || view->byteStride) {
                    if (skip_missing) continue;
                    throw runtime_error("invalid image buffer view");
                }
                if (image->mimeType == glTFImageMimeType::ImagePng)
                    filename = "internal_data.png";
                else if (image->mimeType == glTFImageMimeType::ImageJpeg)
                    filename = "internal_data.jpg";
                else {
                    if (skip_missing) continue;
                    throw runtime_error("unsupported image format");
                }
                data = buffer->data.data() + view->byteOffset;
                data_size = view->byteLength;
            } else {
                // assume it is base64 and find ','
                auto pos = image->uri.find(',');
                if (pos == image->uri.npos) {
                    if (skip_missing) continue;
                    throw runtime_error("could not decode base64 data");
                }
                auto header = image->uri.substr(0, pos);
                for (auto format : {"png", "jpg", "jpeg", "tga", "ppm", "hdr"})
                    if (header.find(format) != header.npos)
                        filename = string("fake.") + format;
                if (is_hdr_filename(filename)) {
                    if (skip_missing) continue;
                    throw runtime_error("unsupported embedded image format " +
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
            throw runtime_error("cannot load image " + filename);
        }
    }
}

// Loads a gltf.
glTF* load_gltf(
    const string& filename, bool load_bin, bool load_image, bool skip_missing) {
    // clear data
    auto gltf = unique_ptr<glTF>(new glTF());

    // load json
    std::ifstream stream(filename.c_str());
    if (!stream) throw runtime_error("could not load json " + filename);
    auto js = json();
    try {
        stream >> js;
    } catch (const exception& e) {
        throw runtime_error(
            string("could not load json with error ") + e.what());
    }

    // parse json
    auto gltf_ = gltf.get();
    try {
        serialize_from_json(gltf_, js);
    } catch (const exception& e) {
        throw runtime_error("error parsing gltf " + string(e.what()));
    }

    // load external resources
    auto dirname = path_dirname(filename);
    if (load_bin) load_buffers(gltf.get(), dirname, skip_missing);
    if (load_image) load_images(gltf.get(), dirname, skip_missing);

    // done
    return gltf.release();
}

// Save buffer data.
void save_buffers(const glTF* gltf, const string& dirname, bool skip_missing) {
    for (auto buffer : gltf->buffers) {
        try {
            if (startswith(buffer->uri, "data:")) {
                throw runtime_error("saving of embedded data not supported");
            }
            save_binary(dirname + buffer->uri, buffer->data);
        } catch (const std::exception&) {
            if (skip_missing) continue;
            throw;
        }
    }
}

// Save images.
void save_images(const glTF* gltf, const string& dirname, bool skip_missing) {
    for (auto image : gltf->images) {
        try {
            if (startswith(image->uri, "data:")) {
                throw runtime_error("saving of embedded data not supported");
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
            if (!ok) { throw runtime_error("cannot save image " + filename); }
        } catch (const std::exception&) {
            if (skip_missing) continue;
            throw;
        }
    }
}

// Saves a gltf.
void save_gltf(
    const string& filename, const glTF* gltf, bool save_bin, bool save_image) {
    // dumps json
    auto js = json();
    serialize_to_json(gltf, js);

    // save json
    save_text(filename, js.dump(2));

    // save external resources
    auto dirname = path_dirname(filename);
    if (save_bin) save_buffers(gltf, dirname, false);
    if (save_image) save_images(gltf, dirname, false);
}

// reading shortcut
template <typename T>
inline void gltf_fread(FILE* f, T* v, int count) {
    if (fread(v, sizeof(T), count, f) != count)
        throw runtime_error("could not read binary file");
}

// writing shortcut
template <typename T>
inline void gltf_fwrite(FILE* f, const T* v, int count) {
    if (fwrite(v, sizeof(T), count, f) != count)
        runtime_error("could not write binary file");
}

// Loads a binary gltf.
glTF* load_binary_gltf(
    const string& filename, bool load_bin, bool load_image, bool skip_missing) {
    // clear data
    auto gltf = unique_ptr<glTF>(new glTF());

    // opens binary file
    auto f = fopen(filename.c_str(), "rb");
    if (!f) throw runtime_error("could not load binary file " + filename);

    // read magic
    uint32_t magic;
    gltf_fread(f, &magic, 1);
    if (magic != 0x46546C67) throw runtime_error("corrupted glb format");

    // read version
    uint32_t version;
    gltf_fread(f, &version, 1);
    if (version != 1 && version != 2)
        throw runtime_error("unsupported glb version");

    // read length
    uint32_t length;
    gltf_fread(f, &length, 1);

    // data
    auto json_bytes = vector<char>();
    auto buffer_bytes = vector<unsigned char>();
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
            throw runtime_error("corrupt binary format");
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
            throw runtime_error("corrupt binary format");

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
    } catch (const exception& e) {
        throw runtime_error(
            string("could not load json with error ") + e.what());
    }

    // parse json
    auto gltf_ = gltf.get();
    try {
        serialize_from_json(gltf_, js);
    } catch (const exception& e) {
        throw runtime_error("cannot parse gltf json " + string(e.what()));
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
void save_binary_gltf(
    const string& filename, const glTF* gltf, bool save_bin, bool save_image) {
    // opens binary file
    auto f = fopen(filename.c_str(), "wb");
    if (!f) throw runtime_error("could not write binary file");

    // dumps json
    auto js = json();
    serialize_to_json(gltf, js);

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
    if (!_valid) throw runtime_error("corrupted glTF accessor view");
}

float accessor_view::get(int idx, int c) const {
    auto i = min(max(c, 0), ncomp() - 1);
    auto valb = _data + _stride * idx + i * _ctype_size(_ctype);
    // use double for integer conversion to attempt to maintain precision
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
                throw runtime_error("bad enum value");
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
                throw runtime_error("bad enum value");
                break;
        }
    }
    return 0;
}

int accessor_view::geti(int idx, int c) const {
    auto i = min(max(c, 0), ncomp() - 1);
    auto valb = _data + _stride * idx + i * _ctype_size(_ctype);
    // use double for integer conversion to attempt to maintain precision
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
            throw runtime_error("bad enum value");
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
inline svg_scene* load_svg(const string& filename) {
    auto svg = nsvgParseFromFile(filename.c_str(), "mm", 96);
    if (!svg) throw runtime_error("cannot load SVG");
    auto scn = new svg_scene();
    for (auto shape = svg->shapes; shape != nullptr; shape = shape->next) {
        auto shp = new svg_shape();
        scn->shapes += shp;
        for (auto path = shape->paths; path != nullptr; path = path->next) {
            auto pth = new svg_path();
            shp->paths += pth;
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
inline void save_svg(const string& filename, const vector<svg_path>& paths) {
    throw runtime_error("not implemented yet");
}

}  // namespace ygl

#endif

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR SHAPE EXAMPLES
// -----------------------------------------------------------------------------
namespace ygl {

// Make a sphere. This is not watertight.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>> make_uvsphere(
    int tesselation, bool flipped) {
    auto quads = vector<vec4i>();
    auto texcoord = vector<vec2f>();
    tie(quads, texcoord) =
        make_uvquads(pow2(tesselation + 2), pow2(tesselation + 1));
    auto pos = vector<vec3f>(texcoord.size());
    auto norm = vector<vec3f>(texcoord.size());
    if (!flipped) {
        for (auto i = 0; i < texcoord.size(); i++) {
            auto uv = texcoord[i];
            auto a = vec2f{2 * pif * uv.x, pif * (1 - uv.y)};
            pos[i] = {cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
            norm[i] = {cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
        }
    } else {
        for (auto i = 0; i < texcoord.size(); i++) {
            auto uv = texcoord[i];
            auto a = vec2f{2 * pif * uv.x, pif * uv.y};
            pos[i] = {cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
            norm[i] = {-cos(a.x) * sin(a.y), -sin(a.x) * sin(a.y), -cos(a.y)};
            texcoord[i] = {uv.x, 1 - uv.y};
        }
    }
    return {quads, pos, norm, texcoord};
}

// Make a geodesic sphere.
tuple<vector<vec3i>, vector<vec3f>> make_geodesicsphere(int tesselation) {
    // https://stackoverflow.com/questions/17705621/algorithm-for-a-geodesic-sphere
    const float X = 0.525731112119133606f;
    const float Z = 0.850650808352039932f;
    auto pos = vector<vec3f>{{-X, 0.0, Z}, {X, 0.0, Z}, {-X, 0.0, -Z},
        {X, 0.0, -Z}, {0.0, Z, X}, {0.0, Z, -X}, {0.0, -Z, X}, {0.0, -Z, -X},
        {Z, X, 0.0}, {-Z, X, 0.0}, {Z, -X, 0.0}, {-Z, -X, 0.0}};
    auto triangles = vector<vec3i>{{0, 1, 4}, {0, 4, 9}, {9, 4, 5}, {4, 8, 5},
        {4, 1, 8}, {8, 1, 10}, {8, 10, 3}, {5, 8, 3}, {5, 3, 2}, {2, 3, 7},
        {7, 3, 10}, {7, 10, 6}, {7, 6, 11}, {11, 6, 0}, {0, 6, 1}, {6, 10, 1},
        {9, 11, 0}, {9, 2, 11}, {9, 5, 2}, {7, 11, 2}};
    for (auto l = 0; l < tesselation - 2; l++) {
        vector<vec2i> _lines;
        vector<vec4i> _quads;
        vector<vec2i> edges;
        vector<vec4i> faces;
        tie(_lines, triangles, _quads, edges, faces) =
            subdivide_elems_linear({}, triangles, {}, (int)pos.size());
        pos = subdivide_vert_linear(pos, edges, faces);
    }
    for (auto& p : pos) p = normalize(p);
    return {triangles, pos};
}

// Make a sphere. This is not watertight.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>
make_uvhemisphere(int tesselation, bool flipped) {
    auto quads = vector<vec4i>();
    auto texcoord = vector<vec2f>();
    tie(quads, texcoord) =
        make_uvquads(pow2(tesselation + 2), pow2(tesselation));
    auto pos = vector<vec3f>(texcoord.size());
    auto norm = vector<vec3f>(texcoord.size());
    if (!flipped) {
        for (auto i = 0; i < texcoord.size(); i++) {
            auto uv = texcoord[i];
            auto a = vec2f{2 * pif * uv.x, pif * 0.5f * (1 - uv.y)};
            pos[i] = {cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
            norm[i] = {cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
        }
    } else {
        for (auto i = 0; i < texcoord.size(); i++) {
            auto uv = texcoord[i];
            auto a = vec2f{2 * pif * uv.x, pif * (0.5f + 0.5f * uv.y)};
            pos[i] = {cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
            norm[i] = {-cos(a.x) * sin(a.y), -sin(a.x) * sin(a.y), -cos(a.y)};
            texcoord[i] = {uv.x, 1 - uv.y};
        }
    }
    return {quads, pos, norm, texcoord};
}

// Make a quad.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>> make_uvquad(
    int tesselation) {
    auto quads = vector<vec4i>();
    auto texcoord = vector<vec2f>();
    tie(quads, texcoord) = make_uvquads(pow2(tesselation), pow2(tesselation));
    auto pos = vector<vec3f>(texcoord.size());
    auto norm = vector<vec3f>(texcoord.size());
    for (auto i = 0; i < texcoord.size(); i++) {
        auto uv = texcoord[i];
        pos[i] = {(-1 + uv.x * 2), (-1 + uv.y * 2), 0};
        norm[i] = {0, 0, 1};
    }
    return {quads, pos, norm, texcoord};
}

// Make a cube with unique vertices. This is watertight but has no
// texture coordinates or normals.
tuple<vector<vec4i>, vector<vec3f>> make_cube(int tesselation) {
    static auto pos = vector<vec3f>{{-1, -1, -1}, {-1, +1, -1}, {+1, +1, -1},
        {+1, -1, -1}, {-1, -1, +1}, {-1, +1, +1}, {+1, +1, +1}, {+1, -1, +1}};
    static auto quads = vector<vec4i>{{0, 1, 2, 3}, {7, 6, 5, 4}, {4, 5, 1, 0},
        {6, 7, 3, 2}, {2, 1, 5, 6}, {0, 3, 7, 4}};

    auto tpos = pos;
    auto tquads = quads;
    for (auto l = 0; l < tesselation; l++) {
        auto lines = vector<vec2i>();
        auto triangles = vector<vec3i>();
        auto edges = vector<vec2i>();
        auto faces = vector<vec4i>();
        tie(lines, triangles, tquads, edges, faces) =
            subdivide_elems_linear({}, {}, tquads, (int)tpos.size());
        tpos = subdivide_vert_linear(tpos, edges, faces);
    }

    return {tquads, tpos};
}

// Make a facevarying cube with unique vertices but different texture
// coordinates.
tuple<vector<vec4i>, vector<vec3f>, vector<vec4i>, vector<vec3f>, vector<vec4i>,
    vector<vec2f>>
make_fvcube(int tesselation) {
    static auto pos = vector<vec3f>{{-1, -1, -1}, {-1, +1, -1}, {+1, +1, -1},
        {+1, -1, -1}, {-1, -1, +1}, {-1, +1, +1}, {+1, +1, +1}, {+1, -1, +1}};
    static auto quads_pos = vector<vec4i>{{0, 1, 2, 3}, {7, 6, 5, 4},
        {4, 5, 1, 0}, {6, 7, 3, 2}, {2, 1, 5, 6}, {0, 3, 7, 4}};
    static auto norm = vector<vec3f>{{0, 0, -1}, {0, 0, -1}, {0, 0, -1},
        {0, 0, -1}, {0, 0, +1}, {0, 0, +1}, {0, 0, +1}, {0, 0, +1}, {-1, 0, 0},
        {-1, 0, 0}, {-1, 0, 0}, {-1, 0, 0}, {+1, 0, 0}, {+1, 0, 0}, {+1, 0, 0},
        {+1, 0, 0}, {0, +1, 0}, {0, +1, 0}, {0, +1, 0}, {0, +1, 0}, {0, -1, 0},
        {0, -1, 0}, {0, -1, 0}, {0, -1, 0}};
    static auto quads_norm = vector<vec4i>{{0, 1, 2, 3}, {4, 5, 6, 7},
        {8, 9, 10, 11}, {12, 13, 14, 15}, {16, 17, 18, 19}, {20, 21, 22, 23}};
    static auto texcoord = vector<vec2f>{{0, 0}, {1, 0}, {1, 1}, {0, 1}, {0, 0},
        {1, 0}, {1, 1}, {0, 1}, {0, 0}, {1, 0}, {1, 1}, {0, 1}, {0, 0}, {1, 0},
        {1, 1}, {0, 1}, {0, 0}, {1, 0}, {1, 1}, {0, 1}, {0, 0}, {1, 0}, {1, 1},
        {0, 1}};
    static auto quads_texcoord = vector<vec4i>{{0, 1, 2, 3}, {4, 5, 6, 7},
        {8, 9, 10, 11}, {12, 13, 14, 15}, {16, 17, 18, 19}, {20, 21, 22, 23}};

    if (!tesselation)
        return {quads_pos, pos, quads_norm, norm, quads_texcoord, texcoord};

    auto tpos = pos, tnorm = norm;
    auto ttexcoord = texcoord;
    auto tquads_pos = quads_pos, tquads_norm = quads_norm,
         tquads_texcoord = quads_texcoord;
    for (auto l = 0; l < tesselation; l++) {
        auto lines = vector<vec2i>();
        auto triangles = vector<vec3i>();
        auto edges = vector<vec2i>();
        auto faces = vector<vec4i>();
        tie(lines, triangles, tquads_pos, edges, faces) =
            subdivide_elems_linear({}, {}, tquads_pos, (int)tpos.size());
        tpos = subdivide_vert_linear(tpos, edges, faces);
        tie(lines, triangles, tquads_norm, edges, faces) =
            subdivide_elems_linear({}, {}, tquads_norm, (int)tnorm.size());
        tnorm = subdivide_vert_linear(tnorm, edges, faces);
        tie(lines, triangles, tquads_texcoord, edges, faces) =
            subdivide_elems_linear(
                {}, {}, tquads_texcoord, (int)ttexcoord.size());
        ttexcoord = subdivide_vert_linear(ttexcoord, edges, faces);
    }

    return {tquads_pos, tpos, tquads_norm, tnorm, tquads_texcoord, ttexcoord};
}

// Make a facevarying sphere with unique vertices but different texture
// coordinates.
tuple<vector<vec4i>, vector<vec3f>, vector<vec4i>, vector<vec3f>, vector<vec4i>,
    vector<vec2f>>
make_fvsphere(int tesselation) {
    auto usteps = pow2(tesselation + 2), vsteps = pow2(tesselation + 1);
    auto qpos = vector<vec4i>(), qtexcoord = vector<vec4i>();
    auto uvpos = vector<vec2f>(), texcoord = vector<vec2f>();
    tie(qpos, uvpos) = make_uvquads(usteps, vsteps, true, false, true, true);
    tie(qtexcoord, texcoord) = make_uvquads(usteps, vsteps);
    auto pos = vector<vec3f>(uvpos.size());
    for (auto i = 0; i < uvpos.size(); i++) {
        auto uv = uvpos[i];
        auto a = vec2f{2 * pif * uv.x, pif * (1 - uv.y)};
        pos[i] = {cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
    }
    return {qpos, pos, qpos, pos, qtexcoord, texcoord};
}

// Make a suzanne monkey model for testing. Note that some quads are
// degenerate.
tuple<vector<vec4i>, vector<vec3f>> make_suzanne(int tesselation) {
    static auto suzanne_pos = vector<vec3f>{{0.4375, 0.1640625, 0.765625},
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
    static auto suzanne_triangles = vector<vec3i>{{60, 64, 48}, {49, 65, 61},
        {62, 64, 60}, {61, 65, 63}, {60, 58, 62}, {63, 59, 61}, {60, 56, 58},
        {59, 57, 61}, {60, 54, 56}, {57, 55, 61}, {60, 52, 54}, {55, 53, 61},
        {60, 50, 52}, {53, 51, 61}, {60, 48, 50}, {51, 49, 61}, {224, 228, 226},
        {227, 229, 225}, {72, 283, 73}, {73, 284, 72}, {341, 347, 383},
        {384, 348, 342}, {299, 345, 343}, {344, 346, 300}, {323, 379, 351},
        {352, 380, 324}, {441, 443, 445}, {446, 444, 442}, {463, 491, 465},
        {466, 492, 464}, {495, 497, 499}, {500, 498, 496}};
    static auto suzanne_quads = vector<vec4i>{{46, 0, 2, 44}, {3, 1, 47, 45},
        {44, 2, 4, 42}, {5, 3, 45, 43}, {2, 8, 6, 4}, {7, 9, 3, 5},
        {0, 10, 8, 2}, {9, 11, 1, 3}, {10, 12, 14, 8}, {15, 13, 11, 9},
        {8, 14, 16, 6}, {17, 15, 9, 7}, {14, 20, 18, 16}, {19, 21, 15, 17},
        {12, 22, 20, 14}, {21, 23, 13, 15}, {22, 24, 26, 20}, {27, 25, 23, 21},
        {20, 26, 28, 18}, {29, 27, 21, 19}, {26, 32, 30, 28}, {31, 33, 27, 29},
        {24, 34, 32, 26}, {33, 35, 25, 27}, {34, 36, 38, 32}, {39, 37, 35, 33},
        {32, 38, 40, 30}, {41, 39, 33, 31}, {38, 44, 42, 40}, {43, 45, 39, 41},
        {36, 46, 44, 38}, {45, 47, 37, 39}, {46, 36, 50, 48}, {51, 37, 47, 49},
        {36, 34, 52, 50}, {53, 35, 37, 51}, {34, 24, 54, 52}, {55, 25, 35, 53},
        {24, 22, 56, 54}, {57, 23, 25, 55}, {22, 12, 58, 56}, {59, 13, 23, 57},
        {12, 10, 62, 58}, {63, 11, 13, 59}, {10, 0, 64, 62}, {65, 1, 11, 63},
        {0, 46, 48, 64}, {49, 47, 1, 65}, {88, 173, 175, 90},
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

    auto tpos = suzanne_pos;
    auto tquads = suzanne_quads;
    tquads.reserve(suzanne_quads.size() + suzanne_triangles.size());
    for (auto& t : suzanne_triangles) {
        tquads.push_back({t.x, t.y, t.z, t.z});
    }

    if (!tesselation) return {tquads, tpos};

    for (auto l = 0; l < tesselation; l++) {
        auto lines = vector<vec2i>();
        auto triangles = vector<vec3i>();
        auto edges = vector<vec2i>();
        auto faces = vector<vec4i>();
        tie(lines, triangles, tquads, edges, faces) =
            subdivide_elems_linear({}, {}, tquads, (int)tpos.size());
        tpos = subdivide_vert_linear(tpos, edges, faces);
    }

    return {tquads, tpos};
}

// Make a cube with uv. This is not watertight.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>> make_uvcube(
    int tesselation) {
    frame3f frames[6] = {frame3f{{1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {0, 0, 1}},
        frame3f{{-1, 0, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, -1}},
        frame3f{{-1, 0, 0}, {0, 0, 1}, {0, 1, 0}, {0, 1, 0}},
        frame3f{{1, 0, 0}, {0, 0, 1}, {0, -1, 0}, {0, -1, 0}},
        frame3f{{0, 1, 0}, {0, 0, 1}, {1, 0, 0}, {1, 0, 0}},
        frame3f{{0, -1, 0}, {0, 0, 1}, {-1, 0, 0}, {-1, 0, 0}}};
    vector<vec3f> quad_pos, quad_norm;
    vector<vec2f> quad_texcoord;
    vector<vec4i> quad_quads;
    tie(quad_quads, quad_pos, quad_norm, quad_texcoord) =
        make_uvquad(tesselation);
    vector<vec3f> pos, norm;
    vector<vec2f> texcoord;
    vector<vec4i> quads;
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
    return {quads, pos, norm, texcoord};
}

// Make a sphere from a cube. This is not watertight.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>
make_uvspherecube(int tesselation) {
    vector<vec3f> pos, norm;
    vector<vec2f> texcoord;
    vector<vec4i> quads;
    tie(quads, pos, norm, texcoord) = make_uvcube(tesselation);
    for (auto i = 0; i < pos.size(); i++) {
        pos[i] = normalize(pos[i]);
        norm[i] = normalize(pos[i]);
    }
    return {quads, pos, norm, texcoord};
}

// Make a cube than stretch it towards a sphere. This is not watertight.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>
make_uvspherizedcube(int tesselation, float radius) {
    vector<vec3f> pos, norm;
    vector<vec2f> texcoord;
    vector<vec4i> quads;
    tie(quads, pos, norm, texcoord) = make_uvcube(tesselation);
    for (auto i = 0; i < pos.size(); i++) {
        norm[i] = normalize(pos[i]);
        pos[i] *= 1 - radius;
        pos[i] += norm[i] * radius;
    }
    norm = compute_normals({}, {}, quads, pos, true);
    return {quads, pos, norm, texcoord};
}

// Make a flipped sphere. This is not watertight.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>
make_uvflipcapsphere(int tesselation, float z, bool flipped) {
    vector<vec3f> pos, norm;
    vector<vec2f> texcoord;
    vector<vec4i> quads;
    tie(quads, pos, norm, texcoord) = make_uvsphere(tesselation, flipped);
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
    return {quads, pos, norm, texcoord};
}

// Make a cutout sphere. This is not watertight.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>
make_uvcutsphere(int tesselation, float z, bool flipped) {
    auto quads = vector<vec4i>();
    auto texcoord = vector<vec2f>();
    tie(quads, texcoord) =
        make_uvquads(pow2(tesselation + 2), pow2(tesselation + 1));
    auto pos = vector<vec3f>(texcoord.size());
    auto norm = vector<vec3f>(texcoord.size());
    if (!flipped) {
        for (auto i = 0; i < texcoord.size(); i++) {
            auto uv = texcoord[i];
            auto p = 1 - acos(z) / pif;
            auto a = vec2f{2 * pif * uv.x, pif * (1 - p * uv.y)};
            pos[i] = {cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
            norm[i] = {cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
        }
    } else {
        for (auto i = 0; i < texcoord.size(); i++) {
            auto uv = texcoord[i];
            auto p = 1 - acos(z) / pif;
            auto a = vec2f{2 * pif * uv.x, pif * ((1 - p) + p * uv.y)};
            pos[i] = {cos(a.x) * sin(a.y), sin(a.x) * sin(a.y), cos(a.y)};
            norm[i] = {-cos(a.x) * sin(a.y), -sin(a.x) * sin(a.y), -cos(a.y)};
            texcoord[i] = {uv.x, (1 - uv.y)};
        }
    }
    return {quads, pos, norm, texcoord};
}

// Make a seashell. This is not watertight. Returns quads, pos, norm, texcoord.
tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>
make_uvseashell(int tesselation, const make_seashell_params& params) {
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

    auto quads = vector<vec4i>();
    auto texcoord = vector<vec2f>();
    tie(quads, texcoord) = make_uvquads(
        pow2(tesselation + 2), pow2(tesselation + 1 + (int)round(R)));
    auto pos = vector<vec3f>(texcoord.size());
    for (auto i = 0; i < texcoord.size(); i++) {
        auto uv = texcoord[i];
        auto s = uv.x * 2 * pif;
        auto t = uv.y * 2 * pif * R - pif * R;
        auto re = 1 / sqrt(pow(cos(s) / e.x, 2) + pow(sin(s) / e.y, 2));
        if (L && W.x && W.y && t > 0) {
            auto l = (2 * pif / N) *
                     ((N * t) / (2 * pif) - floor((N * t) / (2 * pif)));
            auto rn =
                L * exp(-pow((2 * (s - P)) / W.x, 2) - pow((2 * l) / W.y, 2));
            re += rn;
        }
        pos[i].x = (A * sin(b) * cos(t) + cos(s + O.x) * cos(t + O.y) * re -
                       sin(O.z) * sin(t + O.y) * re) *
                   D * exp(t * cot_a);
        pos[i].y = (A * sin(b) * sin(t) + cos(s + O.x) * sin(t + O.y) * re +
                       sin(O.z) * sin(s + O.x) * cos(t + O.y) * re) *
                   exp(t * cot_a);
        pos[i].z =
            (-A * cos(b) + cos(O.z) * sin(s + O.x) * re) * exp(t * cot_a);
        texcoord[i] = {uv.x, uv.y * R};
    }
    auto norm = compute_normals({}, {}, quads, pos);
    return {quads, pos, norm, texcoord};
}

// Make a bezier circle. Returns bezier, pos.
tuple<vector<vec4i>, vector<vec3f>> make_bezier_circle() {
    // constant from http://spencermortensen.com/articles/bezier-circle/
    static auto c = 0.551915024494f;
    static auto pos = vector<vec3f>{{1, 0, 0}, {1, c, 0}, {c, 1, 0}, {0, 1, 0},
        {-c, 1, 0}, {-1, c, 0}, {-1, 0, 0}, {-1, -c, 0}, {-c, -1, 0},
        {0, -1, 0}, {c, -1, 0}, {1, -c, 0}};
    static auto bezier =
        vector<vec4i>{{0, 1, 2, 3}, {3, 4, 5, 6}, {6, 7, 8, 9}, {9, 10, 11, 0}};
    return {bezier, pos};
}

// Make a hair ball around a shape
tuple<vector<vec2i>, vector<vec3f>, vector<vec3f>, vector<vec2f>, vector<float>>
make_hair(int num, int tesselation, const vector<vec3i>& striangles,
    const vector<vec4i>& squads, const vector<vec3f>& spos,
    const vector<vec3f>& snorm, const vector<vec2f>& stexcoord,
    const make_hair_params& params) {
    vector<vec3f> bpos;
    vector<vec3f> bnorm;
    vector<vec2f> btexcoord;
    tie(bpos, bnorm, btexcoord) =
        sample_triangles_points(striangles + convert_quads_to_triangles(squads),
            spos, snorm, stexcoord, num, params.seed);

    auto rng = init_rng(params.seed, 3);
    auto blen = vector<float>(bpos.size());
    for (auto& l : blen)
        l = lerp(params.length.x, params.length.y, next_rand1f(rng));

    auto cidx = vector<int>();
    if (params.clump.x > 0) {
        for (auto bidx : range(bpos.size())) {
            cidx += 0;
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
    auto lines = vector<vec2i>();
    auto texcoord = vector<vec2f>();
    tie(lines, texcoord) = make_uvlines(num, usteps);
    auto pos = vector<vec3f>(texcoord.size());
    auto norm = vector<vec3f>(texcoord.size());
    auto radius = vector<float>(texcoord.size());
    for (auto i : range(texcoord.size())) {
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
        norm = compute_normals(lines, {}, {}, pos);

    return {lines, pos, norm, texcoord, radius};
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
    auto make_camera = [](string name, vec3f from, vec3f to, float yfov,
                           float aperture, float aspect = 16.0f / 9.0f) {
        auto cam = new camera();
        cam->name = name;
        cam->frame = lookat_frame3(from, to, {0, 1, 0});
        cam->aperture = aperture;
        cam->focus = length(from - to);
        cam->yfov = yfov * pif / 180;
        cam->aspect = aspect;
        return cam;
    };

    auto make_instance = [](string name, shape* shp, vec3f pos,
                             vec3f rot = {0, 0, 0}) {
        auto ist = new instance();
        ist->name = name;
        ist->shp = shp;
        ist->frame = {rotation_mat3(vec3f{0, 0, 1}, rot[2] * pif / 180) *
                          rotation_mat3(vec3f{0, 1, 0}, rot[1] * pif / 180) *
                          rotation_mat3(vec3f{1, 0, 0}, rot[0] * pif / 180),
            pos};
        return ist;
    };

    auto make_quad = [](string name, material* mat, float scale = 1) {
        auto shp = new shape();
        shp->mat = mat;
        shp->name = name;
        tie(shp->quads, shp->pos, shp->norm, shp->texcoord) = make_uvquad(0);
        for (auto& p : shp->pos) p *= scale;
        return shp;
    };

    auto make_box = [](string name, material* mat, vec3f scale) {
        auto shp = new shape();
        shp->mat = mat;
        shp->name = name;
        tie(shp->quads, shp->pos, shp->norm, shp->texcoord) = make_uvcube(0);
        for (auto& p : shp->pos) p *= scale;
        return shp;
    };

    auto make_material = [](string name, vec3f kd, vec3f ke = {0, 0, 0}) {
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
    scn->cameras += make_camera("cb_cam", {0, 1, 5.15f}, {0, 1, 0}, 27, 0, 1);
    scn->materials += make_material("cb_white", {0.725f, 0.71f, 0.68f});
    scn->materials += make_material("cb_red", {0.63f, 0.065f, 0.05f});
    scn->materials += make_material("cb_green", {0.14f, 0.45f, 0.091f});
    scn->materials += make_material("cb_light", zero3f, {17, 12, 4});
    scn->shapes += make_quad("cb_floor", scn->materials[0]);
    scn->shapes += make_quad("cb_ceiling", scn->materials[0]);
    scn->shapes += make_quad("cb_back", scn->materials[0]);
    scn->shapes += make_quad("cb_left", scn->materials[2]);
    scn->shapes += make_quad("cb_right", scn->materials[1]);
    scn->shapes +=
        make_box("cb_tallbox", scn->materials[0], {0.3f, 0.6f, 0.3f});
    scn->shapes +=
        make_box("cb_shortbox", scn->materials[0], {0.3f, 0.3f, 0.3f});
    scn->shapes += make_quad("cb_light", scn->materials[3], 0.25f);
    scn->instances +=
        make_instance("cb_floor", scn->shapes[0], {0, 0, 0}, {-90, 0, 0});
    scn->instances +=
        make_instance("cb_ceiling", scn->shapes[1], {0, 2, 0}, {90, 0, 0});
    scn->instances += make_instance("cb_back", scn->shapes[2], {0, 1, -1});
    scn->instances +=
        make_instance("cb_left", scn->shapes[3], {+1, 1, 0}, {0, -90, 0}),
        scn->instances +=
        make_instance("cb_right", scn->shapes[4], {-1, 1, 0}, {0, 90, 0}),
        scn->instances += make_instance(
            "cb_tallbox", scn->shapes[5], {-0.33f, 0.6f, -0.29f}, {0, 15, 0}),
        scn->instances += make_instance(
            "cb_shortbox", scn->shapes[6], {0.33f, 0.3f, 0.33f}, {0, -15, 0}),
        scn->instances +=
        make_instance("cb_light", scn->shapes[7], {0, 1.999f, 0}, {90, 0, 0});
    return scn;
}

//
// Make standard shape. Public API described above.
//
inline tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>
make_uvhollowcutsphere(int tesselation, float radius) {
    auto quads = vector<vec4i>();
    auto pos = vector<vec3f>();
    auto norm = vector<vec3f>();
    auto texcoord = vector<vec2f>();

    vector<vec3f> mpos, mnorm;
    vector<vec2f> mtexcoord;
    vector<vec4i> mquads;
    vector<vec2i> _aux1;
    vector<vec3i> _aux2;

    tie(mquads, mpos, mnorm, mtexcoord) = make_uvcutsphere(tesselation, radius);
    for (auto& uv : mtexcoord) uv.y *= radius;
    tie(_aux1, _aux2, quads) =
        merge_elems((int)pos.size(), {}, {}, quads, {}, {}, mquads);
    pos += mpos;
    norm += mnorm;
    texcoord += mtexcoord;

    tie(mquads, mpos, mnorm, mtexcoord) =
        make_uvcutsphere(tesselation, radius, true);
    for (auto& p : mpos) p *= radius;
    tie(_aux1, _aux2, quads) =
        merge_elems((int)pos.size(), {}, {}, quads, {}, {}, mquads);
    pos += mpos;
    norm += mnorm;
    texcoord += mtexcoord;

    // dpdu = [- s r s0 s1, s r c0 s1, 0] === [- s0, c0, 0]
    // dpdv = [s c0 s1, s s0 s1, s c1] === [c0 s1, s0 s1, c1]
    // n = [c0 c1, - s0 c1, s1]
    tie(mquads, mtexcoord) =
        make_uvquads(pow2(tesselation + 2), pow2(tesselation + 1));
    mpos.resize(mtexcoord.size());
    mnorm.resize(mtexcoord.size());
    for (auto i = 0; i < mtexcoord.size(); i++) {
        auto uv = mtexcoord[i];
        auto a = vec2f{2 * pif * uv[0], pif * (1 - radius)};
        auto r = (1 - uv[1]) + uv[1] * radius;
        mpos[i] = {r * cos(a[0]) * sin(a[1]), r * sin(a[0]) * sin(a[1]),
            r * cos(a[1])};
        mnorm[i] = {-cos(a[0]) * cos(a[1]), -sin(a[0]) * cos(a[1]), sin(a[1])};
        mtexcoord[i] = {uv[0], radius + (1 - radius) * uv[1]};
    }
    tie(_aux1, _aux2, quads) =
        merge_elems((int)pos.size(), {}, {}, quads, {}, {}, mquads);
    pos += mpos;
    norm += mnorm;
    texcoord += mtexcoord;
    return {quads, pos, norm, texcoord};
}

//
// Make standard shape. Public API described above.
//
inline tuple<vector<vec4i>, vector<vec3f>, vector<vec3f>, vector<vec2f>>
make_uvhollowcutsphere1(int tesselation, float radius) {
    auto quads = vector<vec4i>();
    auto pos = vector<vec3f>();
    auto norm = vector<vec3f>();
    auto texcoord = vector<vec2f>();

    vector<vec3f> mpos, mnorm;
    vector<vec2f> mtexcoord;
    vector<vec4i> mquads;
    vector<vec2i> _aux1;
    vector<vec3i> _aux2;

    tie(mquads, mpos, mnorm, mtexcoord) = make_uvcutsphere(tesselation, radius);
    for (auto& uv : mtexcoord) uv.y *= radius;
    for (auto i = (pow2(tesselation + 2) + 1) * pow2(tesselation + 1);
         i < mnorm.size(); i++)
        mnorm[i] = normalize(mnorm[i] + vec3f{0, 0, 1});
    tie(_aux1, _aux2, quads) =
        merge_elems((int)pos.size(), {}, {}, quads, {}, {}, mquads);
    pos += mpos;
    norm += mnorm;
    texcoord += mtexcoord;

    tie(mquads, mpos, mnorm, mtexcoord) =
        make_uvcutsphere(tesselation, radius * 1.05f, true);
    for (auto& p : mpos) p *= 0.8f;
    tie(_aux1, _aux2, quads) =
        merge_elems((int)pos.size(), {}, {}, quads, {}, {}, mquads);
    pos += mpos;
    norm += mnorm;
    texcoord += mtexcoord;

    tie(mquads, mtexcoord) =
        make_uvquads(pow2(tesselation + 2), pow2(tesselation + 1) / 4);
    mpos.resize(mtexcoord.size());
    mnorm.resize(mtexcoord.size());
    for (auto i = 0; i < mtexcoord.size(); i++) {
        auto uv = mtexcoord[i];
        auto p = 1 - acos(radius) / pif;
        auto v = p + uv[1] * (1 - p);
        auto a = vec2f{2 * pif * uv[0], pif * (1 - v)};
        mpos[i] = {cos(a[0]) * sin(a[1]), sin(a[0]) * sin(a[1]),
            (2 * radius - cos(a[1]))};
        mnorm[i] = {-cos(a[0]) * sin(a[1]), -sin(a[0]) * sin(a[1]), cos(a[1])};
        mtexcoord[i] = {uv[0], radius + (1 - radius) * uv[1]};
    }
    for (auto i = 0; i < (pow2(tesselation + 2) + 1); i++)
        mnorm[i] = normalize(mnorm[i] + vec3f{0, 0, 1});
    tie(_aux1, _aux2, quads) =
        merge_elems((int)pos.size(), {}, {}, quads, {}, {}, mquads);
    pos += mpos;
    norm += mnorm;
    texcoord += mtexcoord;
    return {quads, pos, norm, texcoord};
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// EXAMPLE SCENES
// -----------------------------------------------------------------------------
namespace ygl {

// Makes/updates a test texture
void update_test_texture(
    const scene* scn, texture* txt, const test_texture_params& ttxt) {
    if (ttxt.name == "") throw runtime_error("cannot use empty name");

    txt->name = ttxt.name;
    txt->path = "";
    txt->ldr = {};
    txt->hdr = {};

    switch (ttxt.type) {
        case test_texture_type::none: break;
        case test_texture_type::grid: {
            txt->ldr = make_grid_image(ttxt.resolution, ttxt.resolution);
        } break;
        case test_texture_type::checker: {
            txt->ldr = make_checker_image(ttxt.resolution, ttxt.resolution);
        } break;
        case test_texture_type::colored: {
            txt->ldr = make_uvgrid_image(ttxt.resolution, ttxt.resolution);
        } break;
        case test_texture_type::rcolored: {
            txt->ldr = make_recuvgrid_image(ttxt.resolution, ttxt.resolution);
        } break;
        case test_texture_type::bump: {
            txt->ldr = make_bumpdimple_image(
                ttxt.resolution, ttxt.resolution, ttxt.tile_size);
        } break;
        case test_texture_type::uv: {
            txt->ldr = make_uv_image(ttxt.resolution, ttxt.resolution);
        } break;
        case test_texture_type::gamma: {
            txt->ldr = make_gammaramp_image(ttxt.resolution, ttxt.resolution);
        } break;
        case test_texture_type::noise: {
            txt->ldr = make_noise_image(
                ttxt.resolution, ttxt.resolution, ttxt.noise_scale);
        } break;
        case test_texture_type::ridge: {
            txt->ldr = make_ridge_image(
                ttxt.resolution, ttxt.resolution, ttxt.noise_scale);
        } break;
        case test_texture_type::fbm: {
            txt->ldr = make_fbm_image(
                ttxt.resolution, ttxt.resolution, ttxt.noise_scale);
        } break;
        case test_texture_type::turbulence: {
            txt->ldr = make_turbulence_image(
                ttxt.resolution, ttxt.resolution, ttxt.noise_scale);
        } break;
        case test_texture_type::gammaf: {
            txt->hdr = make_gammaramp_imagef(ttxt.resolution, ttxt.resolution);
        } break;
        case test_texture_type::sky: {
            txt->hdr = make_sunsky_image(ttxt.resolution, ttxt.sky_sunangle);
        } break;
        default: throw runtime_error("should not have gotten here");
    }

    if (ttxt.bump_to_normal) {
        txt->ldr = bump_to_normal_map(txt->ldr, ttxt.bump_scale);
    }

    if (txt->ldr) txt->path = ttxt.name + ".png";
    if (txt->hdr) txt->path = ttxt.name + ".sky";
}

// Makes/updates a test material
void update_test_material(
    const scene* scn, material* mat, const test_material_params& tmat) {
    if (tmat.name == "") throw runtime_error("cannot use empty name");

    auto txt = find_named_elem(scn->textures, tmat.texture);
    auto norm = find_named_elem(scn->textures, tmat.normal);

    mat->name = tmat.name;
    mat->type = material_type::specular_roughness;
    mat->ke = zero3f;
    mat->kd = zero3f;
    mat->rs = 1;
    mat->kr = zero3f;
    mat->kt = zero3f;
    mat->ke_txt.txt = nullptr;
    mat->kd_txt.txt = nullptr;
    mat->ks_txt.txt = nullptr;
    mat->kr_txt.txt = nullptr;
    mat->kt_txt.txt = nullptr;

    switch (tmat.type) {
        case test_material_type::none: break;
        case test_material_type::emission: {
            mat->ke = tmat.emission * tmat.color;
            mat->ke_txt.txt = txt;
        } break;
        case test_material_type::matte: {
            mat->kd = tmat.color;
            mat->kd_txt.txt = txt;
        } break;
        case test_material_type::plastic: {
            mat->kd = tmat.color;
            mat->ks = {0.04f, 0.04f, 0.04f};
            mat->rs = tmat.roughness;
            mat->kd_txt.txt = txt;
        } break;
        case test_material_type::metal: {
            mat->ks = tmat.color;
            mat->rs = tmat.roughness;
            mat->ks_txt.txt = txt;
        } break;
        case test_material_type::transparent: {
            mat->kd = tmat.color;
            mat->op = tmat.opacity;
            mat->kd_txt.txt = txt;
        } break;
        default: throw runtime_error("should not have gotten here");
    }

    mat->norm_txt.txt = norm;
}

// Makes/updates a test shape
void update_test_shape(
    const scene* scn, shape* shp, const test_shape_params& tshp) {
    if (tshp.name == "") throw runtime_error("cannot use empty name");
    auto mat = (material*)nullptr;
    if (scn && tshp.material != "") {
        for (auto elem : scn->materials)
            if (elem->name == tshp.material) mat = elem;
    }

    shp->name = tshp.name;
    shp->mat = mat;
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

    switch (tshp.type) {
        case test_shape_type::floor: {
            tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
                make_uvquad((tshp.tesselation < 0) ? 5 : tshp.tesselation);
            for (auto& p : shp->pos) p = {-p.x, p.z, p.y};
            for (auto& n : shp->norm) n = {n.x, n.z, n.y};
            for (auto& p : shp->pos) p *= 20;
            for (auto& uv : shp->texcoord) uv *= 20;
        } break;
        case test_shape_type::quad: {
            tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
                make_uvquad((tshp.tesselation < 0) ? 0 : tshp.tesselation);
        } break;
        case test_shape_type::cube: {
            tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
                make_uvcube((tshp.tesselation < 0) ? 0 : tshp.tesselation);
        } break;
        case test_shape_type::sphere: {
            tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
                make_uvsphere((tshp.tesselation < 0) ? 5 : tshp.tesselation);
        } break;
        case test_shape_type::spherecube: {
            tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
                make_uvspherecube(
                    (tshp.tesselation < 0) ? 4 : tshp.tesselation);
        } break;
        case test_shape_type::spherizedcube: {
            tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
                make_uvspherizedcube(
                    (tshp.tesselation < 0) ? 4 : tshp.tesselation, 0.75f);
        } break;
        case test_shape_type::geosphere: {
            tie(shp->triangles, shp->pos) = make_geodesicsphere(
                (tshp.tesselation < 0) ? 5 : tshp.tesselation);
            shp->norm = shp->pos;
        } break;
        case test_shape_type::flipcapsphere: {
            tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
                make_uvflipcapsphere(
                    (tshp.tesselation < 0) ? 5 : tshp.tesselation, 0.75f);
        } break;
        case test_shape_type::suzanne: {
            tie(shp->quads, shp->pos) =
                make_suzanne((tshp.tesselation < 0) ? 0 : tshp.tesselation);
        } break;
        case test_shape_type::cubep: {
            tie(shp->quads, shp->pos) =
                make_cube((tshp.tesselation < 0) ? 0 : tshp.tesselation);
        } break;
        case test_shape_type::fvcube: {
            tie(shp->quads_pos, shp->pos, shp->quads_norm, shp->norm,
                shp->quads_texcoord, shp->texcoord) =
                make_fvcube((tshp.tesselation < 0) ? 0 : tshp.tesselation);
        } break;
        case test_shape_type::fvsphere: {
            tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
                make_uvsphere((tshp.tesselation < 0) ? 5 : tshp.tesselation);
        } break;
        case test_shape_type::matball: {
            tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
                make_uvflipcapsphere(
                    (tshp.tesselation < 0) ? 5 : tshp.tesselation, 0.75f);
        } break;
        case test_shape_type::point: {
            shp->points.push_back(0);
            shp->pos.push_back({0, 0, 0});
            shp->norm.push_back({0, 0, 1});
            shp->radius.push_back(0.001f);
        } break;
        case test_shape_type::pointscube: {
            auto npoints = (tshp.num < 0) ? 64 * 64 * 16 : tshp.num;
            auto radius = (tshp.radius < 0) ? 0.0025f : tshp.radius;
            tie(shp->points, shp->texcoord) = make_uvpoints(npoints);
            shp->pos.reserve(shp->texcoord.size());
            shp->norm.resize(shp->texcoord.size(), {0, 0, 1});
            shp->radius.resize(shp->texcoord.size(), radius);
            auto rn = init_rng(0);
            for (auto i = 0; i < shp->texcoord.size(); i++) {
                shp->pos += vec3f{-1 + 2 * next_rand1f(rn),
                    -1 + 2 * next_rand1f(rn), -1 + 2 * next_rand1f(rn)};
            }
        } break;
        case test_shape_type::hairball: {
            auto nhairs = (tshp.num < 0) ? 65536 : tshp.num;
            auto radius = (tshp.radius < 0) ? vec2f{0.001f, 0.0001f} :
                                              vec2f{tshp.radius, 0.0001f};
            tie(shp->quads, shp->pos, shp->norm, shp->texcoord) =
                make_uvspherecube(5);
            tie(shp->lines, shp->pos, shp->norm, shp->texcoord, shp->radius) =
                make_hair(nhairs, 2, {}, shp->quads, shp->pos, shp->norm,
                    shp->texcoord, tshp.hair_params);
            shp->quads.clear();
        } break;
        case test_shape_type::beziercircle: {
            tie(shp->beziers, shp->pos) = make_bezier_circle();
            shp->subdivision_level = 2;
        } break;
        default: throw runtime_error("should not have gotten here");
    }

    if (tshp.scale != 1) {
        for (auto& p : shp->pos) p *= tshp.scale;
    }

    for (auto i = 0; i < tshp.subdivision; i++) {
        subdivide_shape_once(shp, true);
    }

    if (tshp.faceted) facet_shape(shp);
}

// Makes/updates a test shape.
void update_test_instance(
    const scene* scn, instance* ist, const test_instance_params& tist) {
    if (tist.name == "") throw runtime_error("cannot use empty name");
    auto shp = (shape*)nullptr;
    if (scn && tist.shape != "") {
        for (auto elem : scn->shapes)
            if (elem->name == tist.shape) shp = elem;
    }

    ist->name = tist.name;
    ist->frame = tist.frame;
    if (tist.rotation != zero3f) {
        auto rot = rotation_mat3(vec3f{0, 0, 1}, tist.rotation.z * pif / 180) *
                   rotation_mat3(vec3f{0, 1, 0}, tist.rotation.y * pif / 180) *
                   rotation_mat3(vec3f{1, 0, 0}, tist.rotation.x * pif / 180);
        ist->frame.rot() = ist->frame.rot() * rot;
    }
    ist->shp = find_named_elem(scn->shapes, tist.shape);
}

// Makes/updates a test shape
void update_test_camera(
    const scene* scn, camera* cam, const test_camera_params& tcam) {
    if (tcam.name == "") throw runtime_error("cannot use empty name");
    cam->name = tcam.name;
    cam->frame = lookat_frame3(tcam.from, tcam.to, vec3f{0, 1, 0});
    cam->yfov = tcam.yfov;
    cam->aspect = tcam.aspect;
    cam->near = 0.01f;
    cam->far = 10000;
    cam->aperture = 0;
    cam->focus = length(tcam.from - tcam.to);
}

// Makes/updates a test shape
void update_test_environment(
    const scene* scn, environment* env, const test_environment_params& tenv) {
    if (tenv.name == "") throw runtime_error("cannot use empty name");
    env->name = tenv.name;
    env->frame = identity_frame3f;
    if (tenv.rotation) {
        env->frame = rotation_frame3({0, 1, 0}, tenv.rotation);
    }
    env->ke = tenv.emission * tenv.color;
    env->ke_txt.txt = find_named_elem(scn->textures, tenv.texture);
}

// Makes/updates a test shape
void update_test_node(
    const scene* scn, node* nde, const test_node_params& tnde) {
    if (tnde.name == "") throw runtime_error("cannot use empty name");
    nde->name = tnde.name;
    nde->frame = tnde.frame;
    nde->translation = tnde.translation;
    nde->rotation = tnde.rotation;
    nde->scaling = tnde.scaling;
    nde->cam = find_named_elem(scn->cameras, tnde.camera);
    auto ist = find_named_elem(scn->instances, tnde.instance);
    nde->ists = (ist) ? vector<instance*>{ist} : vector<instance*>{};
    nde->env = find_named_elem(scn->environments, tnde.environment);
}

// Makes/updates a test animation
void update_test_animation(
    const scene* scn, animation* anm, const test_animation_params& tanm) {
    if (tanm.name == "") throw runtime_error("cannot use empty name");
    if (anm->keyframes.size() != 1) {
        for (auto v : anm->keyframes) delete v;
        anm->keyframes.clear();
        anm->keyframes.push_back(new keyframe());
    }
    anm->name = tanm.name;
    auto kfr = anm->keyframes.front();
    kfr->name = tanm.name;
    kfr->type = (!tanm.bezier) ? keyframe_type::linear : keyframe_type::bezier;
    kfr->times = tanm.times;
    for (auto& v : kfr->times) v *= tanm.speed;
    kfr->translation = tanm.translation;
    kfr->rotation = tanm.rotation;
    kfr->scaling = tanm.scaling;
    for (auto& v : kfr->translation) v *= tanm.scale;
    for (auto& v : kfr->scaling) v *= tanm.scale;
    anm->targets.clear();
    for (auto& nde : tanm.nodes)
        anm->targets.push_back({kfr, find_named_elem(scn->nodes, nde)});
}

// Update test elements
template <typename T, typename T1>
inline void update_test_scene_elem(scene* scn, vector<T*>& elems,
    const vector<T1>& telems, void (*update)(const scene*, T*, const T1&),
    const unordered_set<void*>& refresh) {
    auto emap = unordered_map<string, T*>();
    for (auto elem : elems) emap[elem->name] = elem;
    for (auto& telem : telems) {
        if (!contains(emap, telem.name)) {
            elems += new T();
            update(scn, elems.back(), telem);
        } else {
            auto elem = emap.at(telem.name);
            if (contains(refresh, elem)) update(scn, elem, telem);
        }
    }
}

// Makes/updates a test scene
void update_test_scene(scene* scn, const test_scene_params& params,
    const unordered_set<void*>& refresh) {
    update_test_scene_elem(
        scn, scn->cameras, params.cameras, update_test_camera, refresh);
    update_test_scene_elem(
        scn, scn->textures, params.textures, update_test_texture, refresh);
    update_test_scene_elem(
        scn, scn->materials, params.materials, update_test_material, refresh);
    update_test_scene_elem(
        scn, scn->shapes, params.shapes, update_test_shape, refresh);
    update_test_scene_elem(
        scn, scn->instances, params.instances, update_test_instance, refresh);
    update_test_scene_elem(scn, scn->environments, params.environments,
        update_test_environment, refresh);
    update_test_scene_elem(
        scn, scn->nodes, params.nodes, update_test_node, refresh);
    update_test_scene_elem(scn, scn->animations, params.animations,
        update_test_animation, refresh);
}

// remove duplicate elems
template <typename T>
inline void remove_duplicate_elems(vector<T>& telems) {
    auto names = unordered_set<string>();
    for (auto& elem : telems) names.insert(elem.name);
    if (names.size() == telems.size()) return;
    names.clear();
    auto ntelems = vector<T>();
    ntelems.reserve(names.size());
    for (auto& elem : telems) {
        if (contains(names, elem.name)) continue;
        ntelems += elem;
        names.insert(elem.name);
    }
    telems = ntelems;
}

// remove duplicates
void remove_duplicates(test_scene_params& scn) {
    remove_duplicate_elems(scn.cameras);
    remove_duplicate_elems(scn.textures);
    remove_duplicate_elems(scn.materials);
    remove_duplicate_elems(scn.shapes);
    remove_duplicate_elems(scn.instances);
    remove_duplicate_elems(scn.environments);
}

unordered_map<string, test_texture_params>& test_texture_presets() {
    static auto presets = unordered_map<string, test_texture_params>();
    if (!presets.empty()) return presets;

    auto make_test_texture = [](const string& name, test_texture_type type) {
        auto params = test_texture_params();
        params.name = name;
        params.type = type;
        return params;
    };

    presets["grid"] = make_test_texture("grid", test_texture_type::grid);
    presets["checker"] =
        make_test_texture("checker", test_texture_type::checker);
    presets["colored"] =
        make_test_texture("colored", test_texture_type::colored);
    presets["rcolored"] =
        make_test_texture("rcolored", test_texture_type::rcolored);
    presets["bump"] = make_test_texture("bump", test_texture_type::bump);
    presets["bump"].tile_size = 32;
    presets["tgrid"] = make_test_texture("tgrid", test_texture_type::bump);
    presets["tgrid"].tile_size = 32;
    presets["uv"] = make_test_texture("uv", test_texture_type::uv);
    presets["gamma"] = make_test_texture("gamma", test_texture_type::gamma);
    presets["gridn"] = make_test_texture("gridn", test_texture_type::grid);
    presets["gridn"].bump_to_normal = true;
    presets["gridn"].bump_to_normal = true;
    presets["gridn"].bump_scale = 4;
    presets["tgridn"] = make_test_texture("tgridn", test_texture_type::grid);
    presets["tgridn"].tile_size = 32;
    presets["tgridn"].bump_to_normal = true;
    presets["tgridn"].bump_to_normal = true;
    presets["tgridn"].bump_scale = 4;
    presets["bumpn"] = make_test_texture("bumpn", test_texture_type::bump);
    presets["bumpn"].tile_size = 32;
    presets["bumpn"].bump_to_normal = true;
    presets["bumpn"].bump_scale = 4;
    presets["noise"] = make_test_texture("noise", test_texture_type::noise);
    presets["ridge"] = make_test_texture("ridge", test_texture_type::ridge);
    presets["fbm"] = make_test_texture("fbm", test_texture_type::fbm);
    presets["turbulence"] =
        make_test_texture("turbulence", test_texture_type::turbulence);

    presets["gammaf"] = make_test_texture("gammaf", test_texture_type::gammaf);
    presets["sky1"] = make_test_texture("sky1", test_texture_type::sky);
    presets["sky1"].sky_sunangle = pif / 4;
    presets["sky2"] = make_test_texture("sky2", test_texture_type::sky);
    presets["sky2"].sky_sunangle = pif / 2;

    return presets;
}

unordered_map<string, test_material_params>& test_material_presets() {
    static auto presets = unordered_map<string, test_material_params>();
    if (!presets.empty()) return presets;

    auto make_test_material = [](const string& name, test_material_type type,
                                  const vec3f& color, float roughness = 1) {
        auto params = test_material_params();
        params.name = name;
        params.type = type;
        params.color = color;
        params.roughness = roughness;
        return params;
    };
    auto make_test_materialt = [](const string& name, test_material_type type,
                                   const string& txt, float roughness = 1) {
        auto params = test_material_params();
        params.name = name;
        params.type = type;
        params.color = {1, 1, 1};
        params.roughness = roughness;
        params.texture = txt;
        return params;
    };

    auto emission = test_material_type::emission;
    auto matte = test_material_type::matte;
    auto plastic = test_material_type::plastic;
    auto metal = test_material_type::metal;
    auto transparent = test_material_type::transparent;

    auto gray = vec3f{0.2f, 0.2f, 0.2f};
    auto lgray = vec3f{0.5f, 0.5f, 0.5f};
    auto red = vec3f{0.5f, 0.2f, 0.2f};
    auto green = vec3f{0.2f, 0.5f, 0.2f};
    auto blue = vec3f{0.2f, 0.2f, 0.5f};
    auto white = vec3f{1, 1, 1};

    auto gold = vec3f{0.66f, 0.45f, 0.34f};

    auto rough = 0.25f;
    auto sharp = 0.05f;

    auto params = vector<test_material_params>();

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
    presets["plastic_blue_bumped"].normal = "bumpn";
    presets["plastic_colored_bumped"] = make_test_materialt(
        "plastic_colored_bumped", plastic, "colored", rough);
    presets["plastic_colored_bumped"].normal = "bumpn";

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
    presets["transparent_red"].opacity = 0.9f;
    presets["transparent_green"] =
        make_test_material("transparent_green", transparent, green);
    presets["transparent_green"].opacity = 0.5f;
    presets["transparent_blue"] =
        make_test_material("transparent_blue", transparent, blue);
    presets["transparent_blue"].opacity = 0.2f;

    presets["pointlight"] = make_test_material("pointlight", emission, white);
    presets["pointlight"].emission = 80;
    presets["arealight"] = make_test_material("arealight", emission, white);
    presets["arealight"].emission = 80;

    return presets;
}

unordered_map<string, test_shape_params>& test_shape_presets() {
    static auto presets = unordered_map<string, test_shape_params>();
    if (!presets.empty()) return presets;

    auto make_test_shape = [](const string& name, test_shape_type type,
                               int tesselation = -1, int subdivision = 0,
                               bool faceted = false) {
        auto params = test_shape_params();
        params.name = name;
        params.type = type;
        params.tesselation = tesselation;
        params.subdivision = subdivision;
        params.faceted = faceted;
        return params;
    };

    presets["floor"] = make_test_shape("floor", test_shape_type::floor);
    presets["quad"] = make_test_shape("quad", test_shape_type::quad);
    presets["cube"] = make_test_shape("cube", test_shape_type::cube);
    presets["sphere"] = make_test_shape("sphere", test_shape_type::sphere);
    presets["spherecube"] =
        make_test_shape("spherecube", test_shape_type::spherecube);
    presets["spherizedcube"] =
        make_test_shape("spherizedcube", test_shape_type::spherizedcube);
    presets["flipcapsphere"] =
        make_test_shape("flipcapsphere", test_shape_type::flipcapsphere);
    presets["geosphere"] =
        make_test_shape("geosphere", test_shape_type::geosphere, 5);
    presets["geospheref"] =
        make_test_shape("geospheref", test_shape_type::geosphere, 5, 0, true);
    presets["geospherel"] =
        make_test_shape("geospherel", test_shape_type::geosphere, 4, 0, true);
    presets["cubep"] = make_test_shape("cubep", test_shape_type::cubep);
    presets["cubes"] = make_test_shape("cubes", test_shape_type::cubep, 0, 4);
    presets["suzanne"] = make_test_shape("suzanne", test_shape_type::suzanne);
    presets["suzannes"] =
        make_test_shape("suzannes", test_shape_type::suzanne, 0, 2);
    presets["cubefv"] = make_test_shape("cubefv", test_shape_type::fvcube);
    presets["cubefvs"] =
        make_test_shape("cubefvs", test_shape_type::fvcube, 0, 4);
    presets["spherefv"] =
        make_test_shape("spherefv", test_shape_type::fvsphere);
    presets["matball"] = make_test_shape("matball", test_shape_type::matball);
    presets["matballi"] = make_test_shape("matballi", test_shape_type::sphere);
    presets["matballi"].scale = 0.8f;
    presets["pointscube"] =
        make_test_shape("pointscube", test_shape_type::pointscube);
    presets["hairball1"] =
        make_test_shape("hairball1", test_shape_type::hairball);
    presets["hairball1"].hair_params.radius = {0.001f, 0.0001f};
    presets["hairball1"].hair_params.length = {0.1f, 0.1f};
    presets["hairball1"].hair_params.noise = {0.5f, 8};
    presets["hairball2"] =
        make_test_shape("hairball2", test_shape_type::hairball);
    presets["hairball2"].hair_params.radius = {0.001f, 0.0001f};
    presets["hairball2"].hair_params.length = {0.1f, 0.1f};
    presets["hairball2"].hair_params.clump = {0.5f, 128};
    presets["hairball3"] =
        make_test_shape("hairball3", test_shape_type::hairball);
    presets["hairball3"].hair_params.radius = {0.001f, 0.0001f};
    presets["hairball3"].hair_params.length = {0.1f, 0.1f};
    presets["hairballi"] =
        make_test_shape("hairballi", test_shape_type::sphere);
    presets["hairballi"].scale = 0.8f;
    presets["beziercircle"] =
        make_test_shape("beziercircle", test_shape_type::beziercircle);
    presets["point"] = make_test_shape("point", test_shape_type::point);

    return presets;
}

unordered_map<string, test_environment_params>& test_environment_presets() {
    static auto presets = unordered_map<string, test_environment_params>();
    if (!presets.empty()) return presets;

    auto make_test_environment = [](const string& name, const string& texture) {
        auto params = test_environment_params();
        params.name = name;
        params.color = {1, 1, 1};
        params.texture = texture;
        return params;
    };

    presets["const"] = make_test_environment("const", "");
    presets["sky1"] = make_test_environment("sky1", "");
    presets["sky2"] = make_test_environment("sky2", "");

    return presets;
}

unordered_map<string, test_animation_params>& test_animation_presets() {
    static auto presets = unordered_map<string, test_animation_params>();
    if (!presets.empty()) return presets;

    auto make_test_animation =
        [](const string& name, bool bezier, const vector<float>& times,
            const vector<vec3f>& translation, const vector<quat4f>& rotation,
            const vector<vec3f>& scaling) {
            auto params = test_animation_params();
            params.name = name;
            params.speed = 1;
            params.scale = 1;
            params.bezier = bezier;
            params.times = times;
            params.translation = translation;
            params.rotation = rotation;
            params.scaling = scaling;
            return params;
        };

    presets["bounce"] = make_test_animation(
        "bounce", false, {0, 1, 2}, {{0, 0, 0}, {0, 1, 0}, {0, 0, 0}}, {}, {});
    presets["scale"] = make_test_animation("scale", false, {0, 1, 2}, {}, {},
        {{1, 1, 1}, {0.1f, 0.1f, 0.1f}, {1, 1, 1}});
    presets["rotation"] = make_test_animation("rotation", false, {0, 1, 2}, {},
        {{{0, 1, 0}, 0}, {{0, 1, 0}, pif}, {{0, 1, 0}, 0}}, {});

    return presets;
}

unordered_map<string, test_camera_params>& test_camera_presets() {
    static auto presets = unordered_map<string, test_camera_params>();
    if (!presets.empty()) return presets;

    auto make_test_camera = [](const string& name, const vec3f& from,
                                const vec3f& to, float yfov, float aspect) {
        auto params = test_camera_params();
        params.name = name;
        params.from = from;
        params.to = to;
        params.yfov = yfov;
        params.aspect = aspect;
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

unordered_map<string, test_scene_params>& test_scene_presets() {
    static auto presets = unordered_map<string, test_scene_params>();
    if (!presets.empty()) return presets;

    auto make_test_scene = [](const string& name) {
        auto params = test_scene_params();
        params.name = name;
        return params;
    };
    auto make_test_instance = [](const string& name, const string& shape,
                                  const vec3f& pos = {0, 0, 0}) {
        auto params = test_instance_params();
        params.name = name;
        params.shape = shape;
        params.frame.o = pos;
        return params;
    };

    // textures
    presets["textures"] = make_test_scene("textures");
    for (auto& txt_kv : test_texture_presets())
        presets["textures"].textures += txt_kv.second;

    // shapes
    presets["shapes"] = make_test_scene("shapes");
    for (auto& shp_kv : test_shape_presets())
        presets["shapes"].shapes += shp_kv.second;

    // envmap
    presets["environments"] = make_test_scene("envmaps");
    for (auto& env_kv : test_environment_presets())
        presets["environments"].environments += env_kv.second;

    // simple scenes shared functions
    auto make_simple_scene = [&](const string& name,
                                 const vector<string>& shapes,
                                 const vector<string>& mats,
                                 const string& lights, bool interior = false,
                                 bool nodes = false,
                                 const vector<string>& animations = {}) {
        auto pos = vector<vec3f>{{-2.50f, 1, 0}, {0, 1, 0}, {+2.50f, 1, 0}};
        auto params = make_test_scene(name);
        params.cameras += test_camera_presets().at("cam3");
        params.materials += test_material_presets().at("matte_floor");
        if (params.materials.back().texture != "")
            params.textures +=
                test_texture_presets().at(params.materials.back().texture);
        params.shapes += test_shape_presets().at("floor");
        params.shapes.back().material = params.materials.back().name;
        params.instances += make_test_instance("floor", "floor", {0, 0, 0});
        if (interior) {
            params.materials += test_material_presets().at("matte_gray");
            params.shapes += test_shape_presets().at("sphere");
            params.shapes.back().name = "interior";
            params.shapes.back().material = "matte_gray";
            params.shapes.back().scale = 0.8f;
        }
        for (auto i : range(shapes.size())) {
            auto name = "obj" + to_string(i + 1);
            params.materials += test_material_presets().at(mats[i]);
            params.shapes += test_shape_presets().at(shapes[i]);
            params.shapes.back().name = name;
            params.shapes.back().material = mats[i];
            params.instances += make_test_instance(name, name, pos[i]);
            if (interior) {
                params.instances +=
                    make_test_instance(name + "i", "interior", pos[i]);
            }
            if (!animations.empty()) {
                params.animations += test_animation_presets().at(animations[i]);
                params.animations.back().nodes += name;
            }
        }
        if (lights == "pointlights" || lights == "arealights" ||
            lights == "arealights1") {
            auto emission = 120;
            auto shp = "point";
            auto mat = "pointlight";
            auto pos = vector<vec3f>{{-2, 10, 8}, {+2, 10, 8}};
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
            for (auto i : range(2)) {
                auto name = "light" + to_string(i + 1);
                params.materials += test_material_presets().at(mat);
                params.materials.back().name = name;
                params.materials.back().emission = emission;
                params.shapes += test_shape_presets().at(shp);
                params.shapes.back().name = name;
                params.shapes.back().material = name;
                params.shapes.back().scale = scale;
                params.instances += make_test_instance(name, name, pos[i]);
                if (lights == "arealights" || lights == "arealights1")
                    params.instances.back().frame =
                        lookat_frame3(pos[i], {0, 1, 0}, {0, 0, 1}, true);
            }
        }
        if (lights == "envlights") {
            // env = add_test_environment(params,
            // test_environment_type::sky,
            //     lookat_frame3f({0, 1, 0}, {0, 1, 1}, {0, 1, 0}, true));
        }
        if (!animations.empty() || nodes) {
            for (auto& cam : params.cameras) {
                auto nde = test_node_params();
                nde.name = cam.name;
                nde.frame = lookat_frame3(cam.from, cam.to, vec3f{0, 1, 0});
                nde.camera = cam.name;
                params.nodes += nde;
            }
            for (auto& ist : params.instances) {
                auto nde = test_node_params();
                nde.name = ist.name;
                nde.frame = ist.frame;
                nde.instance = ist.name;
                params.nodes += nde;
            }
            for (auto& env : params.environments) {
                auto nde = test_node_params();
                nde.name = env.name;
                nde.frame = env.frame;
                nde.environment = env.name;
                params.nodes += nde;
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
            {"matte_gray", "matte_gray", "matte_gray"}, "arealights", true);

    // subdiv shapes
    presets["subdiv_al"] =
        make_simple_scene("subdiv_al", {"cubes", "suzannes", "suzannes"},
            {"plastic_red", "plastic_green", "plastic_blue"}, "arealights");

    // plastics shapes
    presets["plastics_al"] =
        make_simple_scene("plastics_al", {"matball", "matball", "matball"},
            {"matte_green", "plastic_green", "plastic_colored"}, "arealights",
            true);
    presets["plastics_el"] = make_simple_scene("plastics_el",
        {"matball", "matball", "matball"},
        {"matte_green", "plastic_green", "plastic_colored"}, "envlights", true);

    // metals shapes
    presets["metals_al"] =
        make_simple_scene("metals_al", {"matball", "matball", "matball"},
            {"gold_rough", "gold_sharp", "silver_sharp"}, "arealights", true);
    presets["metals_el"] =
        make_simple_scene("metals_el", {"matball", "matball", "matball"},
            {"gold_rough", "gold_sharp", "silver_sharp"}, "envlights", true);

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
        "pointlights", false, true, {"bounce", "scale", "rotation"});

    // instances shared functions
    auto make_random_scene = [&](const string& name, const vec2i& num,
                                 const bbox2f& bbox, uint32_t seed = 13) {
        auto rscale = 0.9f * 0.25f *
                      min((bbox.max.x - bbox.min.x) / num.x,
                          (bbox.max.x - bbox.min.x) / num.y);

        auto params = make_test_scene(name);
        params.cameras += test_camera_presets().at("cam3");
        params.materials += test_material_presets().at("matte_floor");
        params.shapes += test_shape_presets().at("floor");
        params.instances += make_test_instance("floor", "floor", {0, 0, 0});
        auto shapes = vector<string>();
        for (auto mat : {"plastic_red", "plastic_green", "plastic_blue"})
            params.materials += test_material_presets().at(mat);
        for (auto shp : {"sphere", "flipcapsphere", "cube"}) {
            for (auto mat : {"plastic_red", "plastic_green", "plastic_blue"}) {
                params.shapes += test_shape_presets().at(shp);
                params.shapes.back().name += "_"s + mat;
                params.shapes.back().material = mat;
                params.shapes.back().scale *= rscale;
                shapes += params.shapes.back().name;
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
                params.instances +=
                    make_test_instance("instance" + to_string(count++),
                        shapes[next_rand1i(rng, (int)shapes.size())], pos);
            }
        }

        auto pos = vector<vec3f>{{-2, 10, 8}, {+2, 10, 8}};
        for (auto i : range(2)) {
            auto name = "light" + to_string(i + 1);
            params.materials += test_material_presets().at("pointlight");
            params.materials.back().name = name;
            params.materials.back().emission = 80;
            params.shapes += test_shape_presets().at("point");
            params.shapes.back().name = name;
            params.shapes.back().material = name;
            params.instances += make_test_instance(name, name, pos[i]);
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
            auto mat = vector<material*>{
                add_test_material(scn, test_material_type::plastic_bumped),
                add_test_material(scn, test_material_type::plastic_bumped)};
            add_test_instance(scn, "obj01",
                add_uvspherecube(scn, "base_obj01", mat[0], 4), {-1.25f, 1, 0});
            add_test_instance(scn, "obj03",
                add_uvspherecube(scn, "subdiv_02_obj02", mat[1], 4), {1.25f, 1, 0});
            }
#endif

    // add missing textures
    for (auto& kv : presets) {
        auto& preset = kv.second;
        auto used = unordered_set<string>();
        for (auto& mat : preset.materials) used.insert(mat.texture);
        for (auto& mat : preset.materials) used.insert(mat.normal);
        for (auto& env : preset.environments) used.insert(env.texture);
        used.erase("");
        for (auto& txt : preset.textures) used.erase(txt.name);
        for (auto& txt : used)
            preset.textures += test_texture_presets().at(txt);
    }

    return presets;
}

// Parses a test_camera object
inline void serialize_from_json(test_camera_params& val, const json& js) {
    static auto def = test_camera_params();
    serialize_from_json_obj(js);
    serialize_from_json_attr(val.name, js, "name");
    serialize_from_json_attr(val.from, js, "from");
    serialize_from_json_attr(val.to, js, "to");
    serialize_from_json_attr(val.yfov, js, "yfov");
    serialize_from_json_attr(val.aspect, js, "aspect");
}

// Parses a test_camera object
inline void serialize_from_json(test_texture_type& val, const json& js) {
    serialize_from_json(val, js, test_texture_names());
}

// Parses a test_camera object
inline void serialize_from_json(test_texture_params& val, const json& js) {
    static auto def = test_texture_params();
    serialize_from_json_obj(js);
    serialize_from_json_attr(val.name, js, "name");
    serialize_from_json_attr(val.type, js, "type");
    serialize_from_json_attr(val.resolution, js, "resolution");
    serialize_from_json_attr(
        val.tile_size, js, "tile_size", false, def.tile_size);
    serialize_from_json_attr(
        val.noise_scale, js, "noise_scale", false, def.noise_scale);
    serialize_from_json_attr(
        val.sky_sunangle, js, "sky_sunangle", false, def.sky_sunangle);
    serialize_from_json_attr(
        val.bump_to_normal, js, "bump_to_normal", false, def.bump_to_normal);
    serialize_from_json_attr(
        val.bump_scale, js, "bump_scale", false, def.bump_scale);
}

// Parses a test_camera object
inline void serialize_from_json(test_material_type& val, const json& js) {
    serialize_from_json(val, js, test_material_names());
}

// Parses a test_camera object
inline void serialize_from_json(test_material_params& val, const json& js) {
    static auto def = test_material_params();
    serialize_from_json_obj(js);
    serialize_from_json_attr(val.name, js, "name");
    serialize_from_json_attr(val.type, js, "type");
    serialize_from_json_attr(val.emission, js, "emission", false, def.emission);
    serialize_from_json_attr(val.color, js, "color", false, def.color);
    serialize_from_json_attr(val.opacity, js, "opacity", false, def.opacity);
    serialize_from_json_attr(
        val.roughness, js, "roughness", false, def.roughness);
    serialize_from_json_attr(val.texture, js, "texture", false, def.texture);
    serialize_from_json_attr(val.normal, js, "normal", false, def.normal);
}

// Parses a test_camera object
inline void serialize_from_json(test_shape_type& val, const json& js) {
    serialize_from_json(val, js, test_shape_names());
}

// Parses a test_camera object
inline void serialize_from_json(test_shape_params& val, const json& js) {
    static auto def = test_shape_params();
    serialize_from_json_obj(js);
    serialize_from_json_attr(val.name, js, "name");
    serialize_from_json_attr(val.type, js, "type");
    serialize_from_json_attr(val.material, js, "material", false, def.material);
    serialize_from_json_attr(
        val.tesselation, js, "tesselation", false, def.tesselation);
    serialize_from_json_attr(
        val.subdivision, js, "subdivision", false, def.subdivision);
    serialize_from_json_attr(val.scale, js, "scale", false, def.scale);
    serialize_from_json_attr(val.radius, js, "radius", false, def.radius);
    serialize_from_json_attr(val.faceted, js, "faceted", false, def.faceted);
    serialize_from_json_attr(val.num, js, "num", false, def.num);
    // TODO: hair parameters
}

// Parses a test_camera object
inline void serialize_from_json(test_instance_params& val, const json& js) {
    static auto def = test_instance_params();
    serialize_from_json_obj(js);
    serialize_from_json_attr(val.name, js, "name");
    serialize_from_json_attr(val.shape, js, "shape");
    serialize_from_json_attr(val.frame, js, "frame", false, def.frame);
    serialize_from_json_attr(val.rotation, js, "rotation", false, def.rotation);
}

// Parses a test_camera object
inline void serialize_from_json(test_environment_params& val, const json& js) {
    static auto def = test_environment_params();
    serialize_from_json_obj(js);
    serialize_from_json_attr(val.name, js, "name");
    serialize_from_json_attr(val.emission, js, "emission", false, def.emission);
    serialize_from_json_attr(val.color, js, "color", false, def.color);
    serialize_from_json_attr(val.texture, js, "texture", false, def.texture);
    serialize_from_json_attr(val.frame, js, "frame", false, def.frame);
    serialize_from_json_attr(val.rotation, js, "rotation", false, def.rotation);
}

// Parses a test_camera object
inline void serialize_from_json(test_scene_params& val, const json& js) {
    static auto def = test_scene_params();
    serialize_from_json_obj(js);
    serialize_from_json_attr(val.name, js, "name", false, def.name);
    serialize_from_json_attr(val.cameras, js, "cameras", false, def.cameras);
    serialize_from_json_attr(val.textures, js, "textures", false, def.textures);
    serialize_from_json_attr(
        val.materials, js, "materials", false, def.materials);
    serialize_from_json_attr(val.shapes, js, "shapes", false, def.shapes);
    serialize_from_json_attr(
        val.environments, js, "environments", false, def.environments);
}

// Converts a test_camera object to JSON
inline void serialize_to_json(const test_camera_params& val, json& js) {
    static auto def = test_camera_params();
    serialize_to_json_obj(js);
    serialize_to_json_attr(val.name, js, "name");
    serialize_to_json_attr(val.from, js, "from");
    serialize_to_json_attr(val.to, js, "to");
    serialize_to_json_attr(val.yfov, js, "yfov");
    serialize_to_json_attr(val.aspect, js, "aspect");
}

// Converts a test_camera object to JSON
inline void serialize_to_json(const test_texture_type& val, json& js) {
    serialize_to_json(val, js, test_texture_names());
}

// Converts a test_camera object to JSON
inline void serialize_to_json(const test_texture_params& val, json& js) {
    static auto def = test_texture_params();
    serialize_to_json_obj(js);
    serialize_to_json_attr(val.name, js, "name");
    serialize_to_json_attr(val.type, js, "type");
    serialize_to_json_attr(val.resolution, js, "resolution");
    serialize_to_json_attr(
        val.tile_size, js, "tile_size", false, def.tile_size);
    serialize_to_json_attr(
        val.noise_scale, js, "noise_scale", false, def.noise_scale);
    serialize_to_json_attr(
        val.sky_sunangle, js, "sky_sunangle", false, def.sky_sunangle);
    serialize_to_json_attr(
        val.bump_to_normal, js, "bump_to_normal", false, def.bump_to_normal);
    serialize_to_json_attr(
        val.bump_scale, js, "bump_scale", false, def.bump_scale);
}

// Converts a test_camera object to JSON
inline void serialize_to_json(const test_material_type& val, json& js) {
    serialize_to_json(val, js, test_material_names());
}

// Converts a test_camera object to JSON
inline void serialize_to_json(const test_material_params& val, json& js) {
    static auto def = test_material_params();
    serialize_to_json_obj(js);
    serialize_to_json_attr(val.name, js, "name");
    serialize_to_json_attr(val.type, js, "type");
    serialize_to_json_attr(val.emission, js, "emission", false, def.emission);
    serialize_to_json_attr(val.color, js, "color", false, def.color);
    serialize_to_json_attr(val.opacity, js, "opacity", false, def.opacity);
    serialize_to_json_attr(
        val.roughness, js, "roughness", false, def.roughness);
    serialize_to_json_attr(val.texture, js, "texture", false, def.texture);
    serialize_to_json_attr(val.normal, js, "normal", false, def.normal);
}

// Parses a test_camera object
inline void serialize_to_json(const test_shape_type& val, json& js) {
    serialize_to_json(val, js, test_shape_names());
}

// Parses a test_camera object
inline void serialize_to_json(const test_shape_params& val, json& js) {
    static auto def = test_shape_params();
    serialize_to_json_obj(js);
    serialize_to_json_attr(val.name, js, "name");
    serialize_to_json_attr(val.type, js, "type");
    serialize_to_json_attr(val.material, js, "material", false, def.material);
    serialize_to_json_attr(
        val.tesselation, js, "tesselation", false, def.tesselation);
    serialize_to_json_attr(
        val.subdivision, js, "subdivision", false, def.subdivision);
    serialize_to_json_attr(val.scale, js, "scale", false, def.scale);
    serialize_to_json_attr(val.radius, js, "radius", false, def.radius);
    serialize_to_json_attr(val.faceted, js, "faceted", false, def.faceted);
    serialize_to_json_attr(val.num, js, "num", false, def.num);
    // TODO: hair parameters
}

// Parses a test_camera object
inline void serialize_to_json(const test_instance_params& val, json& js) {
    static auto def = test_instance_params();
    serialize_to_json_obj(js);
    serialize_to_json_attr(val.name, js, "name");
    serialize_to_json_attr(val.shape, js, "shape");
    serialize_to_json_attr(val.frame, js, "frame", false, def.frame);
    serialize_to_json_attr(val.rotation, js, "rotation", false, def.rotation);
}

// Parses a test_camera object
inline void serialize_to_json(const test_environment_params& val, json& js) {
    static auto def = test_environment_params();
    serialize_to_json_obj(js);
    serialize_to_json_attr(val.name, js, "name");
    serialize_to_json_attr(val.emission, js, "emission", false, def.emission);
    serialize_to_json_attr(val.color, js, "color", false, def.color);
    serialize_to_json_attr(val.texture, js, "texture", false, def.texture);
    serialize_to_json_attr(val.frame, js, "frame", false, def.frame);
    serialize_to_json_attr(val.rotation, js, "rotation", false, def.rotation);
}

// Parses a test_camera object
inline void serialize_to_json(const test_scene_params& val, json& js) {
    static auto def = test_scene_params();
    serialize_to_json_obj(js);
    serialize_to_json_attr(val.name, js, "name", false, ""s);
    serialize_to_json_attr(val.cameras, js, "cameras", false, def.cameras);
    serialize_to_json_attr(val.textures, js, "textures", false, def.textures);
    serialize_to_json_attr(
        val.materials, js, "materials", false, def.materials);
    serialize_to_json_attr(val.shapes, js, "shapes", false, def.shapes);
    serialize_to_json_attr(
        val.environments, js, "environments", false, def.environments);
}

// Load test scene
test_scene_params load_test_scene(const string& filename) {
    // load json
    std::ifstream stream(filename.c_str());
    if (!stream) throw runtime_error("could not load json " + filename);
    auto js = json();
    try {
        stream >> js;
    } catch (const exception& e) {
        throw runtime_error(
            string("could not load json with error ") + e.what());
    }

    // clear data
    auto scn = test_scene_params();
    try {
        serialize_from_json(scn, js);
    } catch (const exception& e) {
        throw runtime_error("error parsing test scene " + string(e.what()));
    }

    // done
    return scn;
}

// Save test scene
void save_test_scene(const string& filename, const test_scene_params& scn) {
    auto js = json();
    serialize_to_json(scn, js);
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

// Enable/disable edges. Attempts to avoid z-fighting but the method is not
// robust.
void gl_enable_edges(bool enabled, float tolerance) {
    assert(gl_check_error());
    if (enabled) {
        glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
        glDepthRange(0, tolerance);
    } else {
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
        glDepthRange(0, 1);
    }
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
// TEXTURE FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Implementation of update_texture.
void _update_texture(gl_texture& txt, int w, int h, int nc, const void* pixels,
    bool floats, bool linear, bool mipmap, bool as_float, bool as_srgb) {
    auto refresh = !txt._tid || txt._width != w || txt._height != h ||
                   txt._ncomp != nc || txt._float != as_float ||
                   txt._srgb != as_srgb || txt._mipmap != mipmap ||
                   txt._linear != linear;
    txt._width = w;
    txt._height = h;
    txt._ncomp = nc;
    txt._float = as_float;
    txt._srgb = as_srgb;
    txt._mipmap = mipmap;
    txt._linear = linear;
    assert(!as_srgb || !as_float);
    assert(gl_check_error());
    if (w * h) {
        int formats_ub[4] = {GL_RED, GL_RG, GL_RGB, GL_RGBA};
        int formats_sub[4] = {GL_RED, GL_RG, GL_SRGB, GL_SRGB_ALPHA};
        int formats_f[4] = {GL_R32F, GL_RG32F, GL_RGB32F, GL_RGBA32F};
        int* formats =
            (as_float) ? formats_f : ((as_srgb) ? formats_sub : formats_ub);
        assert(gl_check_error());
        if (!txt._tid) glGenTextures(1, &txt._tid);
        glBindTexture(GL_TEXTURE_2D, txt._tid);
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
        if (txt._tid) {
            glBindTexture(GL_TEXTURE_2D, txt._tid);
            glDeleteTextures(1, &txt._tid);
            txt._tid = 0;
            glBindTexture(GL_TEXTURE_2D, 0);
        }
    }
    assert(gl_check_error());
}

// Binds a texture to a texture unit
void bind_texture(const gl_texture& txt, uint unit) {
    glActiveTexture(GL_TEXTURE0 + unit);
    glBindTexture(GL_TEXTURE_2D, txt._tid);
}

// Unbinds
void unbind_texture(const gl_texture& txt, uint unit) {
    glActiveTexture(GL_TEXTURE0 + unit);
    glBindTexture(GL_TEXTURE_2D, 0);
}

// Destroys the texture tid.
void clear_texture(gl_texture& txt) {
    assert(gl_check_error());
    glDeleteTextures(1, &txt._tid);
    txt._tid = 0;
    assert(gl_check_error());
}

// -----------------------------------------------------------------------------
// VERTEX ARRAY BUFFER
// -----------------------------------------------------------------------------

// Updates the buffer with new data.
void _update_vertex_buffer(gl_vertex_buffer& buf, int n, int nc,
    const void* values, bool as_float, bool dynamic) {
    auto resize =
        !buf._bid || n * nc != buf._num * buf._ncomp || as_float != buf._float;
    buf._num = n;
    buf._ncomp = nc;
    buf._float = as_float;
    assert(gl_check_error());
    if (n) {
        if (!buf._bid) glGenBuffers(1, &buf._bid);
        glBindBuffer(GL_ARRAY_BUFFER, buf._bid);
        if (resize) {
            glBufferData(GL_ARRAY_BUFFER,
                buf._num * buf._ncomp *
                    ((as_float) ? sizeof(float) : sizeof(int)),
                values, (dynamic) ? GL_DYNAMIC_DRAW : GL_STATIC_DRAW);
        } else {
            glBufferSubData(GL_ARRAY_BUFFER, 0,
                buf._num * buf._ncomp *
                    ((as_float) ? sizeof(float) : sizeof(int)),
                values);
        }
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    } else {
        if (buf._bid) {
            glBindBuffer(GL_ARRAY_BUFFER, buf._bid);
            glDeleteBuffers(1, &buf._bid);
            buf._bid = 0;
            glBindBuffer(GL_ARRAY_BUFFER, 0);
        }
    }
    assert(gl_check_error());
}

// Bind the buffer at a particular attribute location
void bind_vertex_buffer(const gl_vertex_buffer& buf, uint vattr) {
    glEnableVertexAttribArray(vattr);
    glBindBuffer(GL_ARRAY_BUFFER, buf._bid);
    glVertexAttribPointer(vattr, buf._ncomp, GL_FLOAT, false, 0, 0);
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
    glDeleteBuffers(1, &buf._bid);
    buf._bid = 0;
    assert(gl_check_error());
}

// -----------------------------------------------------------------------------
// VERTEX ELEMENTS BUFFER
// -----------------------------------------------------------------------------

// Updates the buffer bid with new data.
void _update_element_buffer(
    gl_element_buffer& buf, int n, int nc, const int* values, bool dynamic) {
    auto resize = !buf._bid || n * nc != buf._num * buf._ncomp;
    buf._num = n;
    buf._ncomp = nc;
    assert(gl_check_error());
    if (n) {
        if (!buf._bid) glGenBuffers(1, &buf._bid);
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf._bid);
        if (resize) {
            glBufferData(GL_ELEMENT_ARRAY_BUFFER,
                buf._num * buf._ncomp * sizeof(int), values,
                (dynamic) ? GL_DYNAMIC_DRAW : GL_STATIC_DRAW);
        } else {
            glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0,
                buf._num * buf._ncomp * sizeof(int), values);
        }
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    } else {
        if (buf._bid) {
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf._bid);
            glDeleteBuffers(1, &buf._bid);
            buf._bid = 0;
            glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        }
    }
    assert(gl_check_error());
}

// Draws elements.
void draw_elems(const gl_element_buffer& buf) {
    if (!buf._bid) return;
    assert(gl_check_error());
    int mode = 0;
    switch (buf._ncomp) {
        case 1: mode = GL_POINTS; break;
        case 2: mode = GL_LINES; break;
        case 3: mode = GL_TRIANGLES; break;
        case 4: mode = GL_QUADS; break;
        default: assert(false);
    };
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, buf._bid);
    glDrawElements(mode, buf._ncomp * buf._num, GL_UNSIGNED_INT, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    assert(gl_check_error());
}

// Destroys the buffer
void clear_element_buffer(gl_element_buffer& buf) {
    assert(gl_check_error());
    glDeleteBuffers(1, &buf._bid);
    buf._bid = 0;
    assert(gl_check_error());
}

// -----------------------------------------------------------------------------
// PROGRAM FUNCTIONS
// -----------------------------------------------------------------------------
// Creates and OpenGL program from vertex and fragment code. Returns the
// program id. Optionally return vertex and fragment shader ids. A VAO is
// created.
gl_program make_program(const string& vertex, const string& fragment) {
    auto prog = gl_program();

    assert(gl_check_error());
    glGenVertexArrays(1, &prog._vao);
    glBindVertexArray(prog._vao);
    assert(gl_check_error());

    int errflags;
    char errbuf[10000];

    // create vertex
    prog._vid = glCreateShader(GL_VERTEX_SHADER);
    const char* vertex_str = vertex.c_str();
    glShaderSource(prog._vid, 1, &vertex_str, NULL);
    glCompileShader(prog._vid);
    glGetShaderiv(prog._vid, GL_COMPILE_STATUS, &errflags);
    if (!errflags) {
        glGetShaderInfoLog(prog._vid, 10000, 0, errbuf);
        throw runtime_error(string("shader not compiled\n\n") + errbuf);
    }
    assert(glGetError() == GL_NO_ERROR);

    // create fragment
    prog._fid = glCreateShader(GL_FRAGMENT_SHADER);
    const char* fragment_str = fragment.c_str();
    glShaderSource(prog._fid, 1, &fragment_str, NULL);
    glCompileShader(prog._fid);
    glGetShaderiv(prog._fid, GL_COMPILE_STATUS, &errflags);
    if (!errflags) {
        glGetShaderInfoLog(prog._fid, 10000, 0, errbuf);
        throw runtime_error(string("shader not compiled\n\n") + errbuf);
    }
    assert(glGetError() == GL_NO_ERROR);

    // create program
    prog._pid = glCreateProgram();
    glAttachShader(prog._pid, prog._vid);
    glAttachShader(prog._pid, prog._fid);
    glLinkProgram(prog._pid);
    glValidateProgram(prog._pid);
    glGetProgramiv(prog._pid, GL_LINK_STATUS, &errflags);
    if (!errflags) {
        glGetProgramInfoLog(prog._pid, 10000, 0, errbuf);
        throw runtime_error(string("program not linked\n\n") + errbuf);
    }
    glGetProgramiv(prog._pid, GL_VALIDATE_STATUS, &errflags);
    if (!errflags) {
        glGetProgramInfoLog(prog._pid, 10000, 0, errbuf);
        throw runtime_error(string("program not linked\n\n") + errbuf);
    }
    assert(gl_check_error());

    glBindVertexArray(0);
    assert(gl_check_error());

    return prog;
}

// Destroys the program pid and optionally the sahders vid and fid.
void clear_program(gl_program& prog) {
    assert(gl_check_error());
    glDetachShader(prog._pid, prog._vid);
    glDeleteShader(prog._vid);
    prog._vid = 0;
    glDetachShader(prog._pid, prog._fid);
    glDeleteShader(prog._fid);
    prog._fid = 0;
    glDeleteProgram(prog._pid);
    prog._pid = 0;
    glDeleteVertexArrays(1, &prog._vao);
    prog._vao = 0;
    assert(gl_check_error());
}

// Get uniform location (simple GL wrapper that avoids GL includes)
int get_program_uniform_location(const gl_program& prog, const string& name) {
    assert(gl_check_error());
    return glGetUniformLocation(prog._pid, name.c_str());
    assert(gl_check_error());
}

// Get uniform location (simple GL wrapper that avoids GL includes)
int get_program_attrib_location(const gl_program& prog, const string& name) {
    assert(gl_check_error());
    return glGetAttribLocation(prog._pid, name.c_str());
    assert(gl_check_error());
}

// Get the names of all uniforms
vector<pair<string, int>> get_program_uniforms_names(const gl_program& prog) {
    auto num = 0;
    assert(gl_check_error());
    glGetProgramiv(prog._pid, GL_ACTIVE_UNIFORMS, &num);
    assert(gl_check_error());
    auto names = vector<pair<string, int>>();
    for (auto i = 0; i < num; i++) {
        char name[4096];
        auto size = 0, length = 0;
        GLenum type;
        glGetActiveUniform(prog._pid, i, 4096, &length, &size, &type, name);
        if (length > 3 && name[length - 1] == ']' && name[length - 2] == '0' &&
            name[length - 3] == '[')
            name[length - 3] = 0;
        auto loc = glGetUniformLocation(prog._pid, name);
        if (loc < 0) continue;
        names.push_back({name, loc});
        assert(gl_check_error());
    }
    return names;
}

// Get the names of all attributes
vector<pair<string, int>> get_program_attributes_names(const gl_program& prog) {
    auto num = 0;
    assert(gl_check_error());
    glGetProgramiv(prog._pid, GL_ACTIVE_ATTRIBUTES, &num);
    assert(gl_check_error());
    auto names = vector<pair<string, int>>();
    for (auto i = 0; i < num; i++) {
        char name[4096];
        auto size = 0;
        GLenum type;
        glGetActiveAttrib(prog._pid, i, 4096, nullptr, &size, &type, name);
        auto loc = glGetAttribLocation(prog._pid, name);
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

// Set uniform texture id tid and unit tunit for program pid and variable
// var.
bool set_program_uniform_texture(
    gl_program& prog, int pos, const gl_texture_info& tinfo, uint tunit) {
    static const auto wrap_mode_map =
        map<gl_texture_wrap, uint>{{gl_texture_wrap::repeat, GL_REPEAT},
            {gl_texture_wrap::clamp, GL_CLAMP_TO_EDGE},
            {gl_texture_wrap::mirror, GL_MIRRORED_REPEAT}};
    static const auto filter_mode_map = map<gl_texture_filter, uint>{
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

// Sets a vartex attribute for program pid and variable var to the buffer
// bid. The attribute has nc components and per-vertex values values.
bool set_program_vertattr(
    gl_program& prog, const string& var, const gl_vertex_buffer& buf) {
    assert(gl_check_error());
    int pos = glGetAttribLocation(prog._pid, var.c_str());
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

// Sets a vartex attribute for program pid and variable var. The attribute
// has nc components and either buffer bid or a single value def
// (if bid is zero). Convenience wrapper to above functions.
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
    if (!prog._pid) return;
    glBindVertexArray(prog._vao);
    glUseProgram(prog._pid);
    assert(gl_check_error());
}

// Unbind a program
void unbind_program(const gl_program& prog) {
    assert(gl_check_error());
    glUseProgram(0);
    glBindVertexArray(0);
    assert(gl_check_error());
}

// Initialize the program. Call with true only after the GL is initialized.
gl_stdimage_program make_stdimage_program() {
    string _header =
        R"(
        #version 330

        float pi = 3.14159265;

        uniform vec2 offset;
        uniform vec2 win_size;
        uniform float zoom;

        uniform sampler2D img;

        )";

    string _vert =
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

    string _frag_tonemap =
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

    string _frag_main =
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
    prog._prog =
        make_program(_header + _vert, _header + _frag_tonemap + _frag_main);

    prog._vbo = make_vertex_buffer(
        vector<vec2f>{{0, 0}, {0, 1}, {1, 1}, {1, 0}}, false);
    prog._ebo = make_element_buffer(vector<vec3i>{{0, 1, 2}, {0, 2, 3}}, false);
    return prog;
}

// Initialize a standard shader. Call with true only after the gl has
// been initialized
gl_stdsurface_program make_stdsurface_program() {
#ifndef _WIN32
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverlength-strings"
#endif
    string _vert_header =
        R"(
        #version 330

        )";

    string _vert_skinning =
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

    string _vert_main =
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

    string _frag_header =
        R"(
        #version 330

        float pi = 3.14159265;

        )";

    string _frag_tonemap =
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

    string _frag_lighting =
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

    string _frag_brdf =
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

    string _frag_material =
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

    string _frag_main =
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
    prog._prog = make_program(_vert_header + _vert_skinning + _vert_main,
        _frag_header + _frag_tonemap + _frag_lighting + _frag_brdf +
            _frag_material + _frag_main);
    assert(gl_check_error());
    return prog;
}

// Initialize gl_stdsurface_program draw state
gl_stdsurface_state* make_stdsurface_state() {
    auto st = new gl_stdsurface_state();
    if (!is_program_valid(st->prog)) st->prog = make_stdsurface_program();
    st->txt[nullptr] = {};
    return st;
}

// Clear fl_stdsurface_program state
void clear_stdsurface_state(gl_stdsurface_state* st, bool clear_program) {
    for (auto& txt_kv : st->txt) {
        auto& txt = txt_kv.second;
        clear_texture(txt);
    }
    st->txt.clear();
    for (auto& vbo_kv : st->vbo) {
        auto& vbo = vbo_kv.second;
        clear_vertex_buffer(vbo.pos);
        clear_vertex_buffer(vbo.norm);
        clear_vertex_buffer(vbo.texcoord);
        clear_vertex_buffer(vbo.color);
        clear_vertex_buffer(vbo.tangsp);
        clear_element_buffer(vbo.points);
        clear_element_buffer(vbo.lines);
        clear_element_buffer(vbo.triangles);
        clear_element_buffer(vbo.quads);
        clear_element_buffer(vbo.edges);
    }
    st->vbo.clear();
    if (clear_program) st->vbo.clear();
}

// Initialize stdsurface lights
gl_stdsurface_lights make_stdsurface_lights(const scene* scn) {
    auto lights = gl_stdsurface_lights();
    for (auto ist : scn->instances) {
        if (!ist->shp) continue;
        if (!ist->shp->mat) continue;
        if (ist->shp->mat->ke == zero3f) continue;
        if (lights.pos.size() >= 16) break;
        if (!ist->shp->points.empty()) {
            for (auto p : ist->shp->points) {
                if (lights.pos.size() >= 16) break;
                lights.pos += transform_point(ist->frame, ist->shp->pos[p]);
                lights.ke += ist->shp->mat->ke;
                lights.ltype += gl_ltype::point;
            }
        } else {
            auto bbox = make_bbox(ist->shp->pos.size(), ist->shp->pos.data());
            auto pos = bbox_center(bbox);
            auto area = 0.0f;
            for (auto l : ist->shp->lines)
                area += line_length(ist->shp->pos[l.x], ist->shp->pos[l.y]);
            for (auto t : ist->shp->triangles)
                area += triangle_area(
                    ist->shp->pos[t.x], ist->shp->pos[t.y], ist->shp->pos[t.z]);
            for (auto t : ist->shp->quads)
                area += quad_area(ist->shp->pos[t.x], ist->shp->pos[t.y],
                    ist->shp->pos[t.z], ist->shp->pos[t.w]);
            auto ke = ist->shp->mat->ke * area;
            if (lights.pos.size() < 16) {
                lights.pos += transform_point(ist->frame, pos);
                lights.ke += ke;
                lights.ltype += gl_ltype::point;
            }
        }
    }
    return lights;
}

// Init shading
void update_stdsurface_state(gl_stdsurface_state* st, const scene* scn,
    const gl_stdsurface_params& params, const unordered_set<void*>& refresh) {
    // update textures -----------------------------------------------------
    for (auto txt : scn->textures) {
        if (st->txt.find(txt) == st->txt.end()) {
            st->txt[txt] = gl_texture();
        } else {
            if (refresh.find(txt) == refresh.end()) continue;
        }
        if (txt->hdr) {
            update_texture(st->txt[txt], txt->hdr, true, true, true);
        } else if (txt->ldr) {
            update_texture(st->txt[txt], txt->ldr, true, true, true);
        } else
            assert(false);
    }

    // refresh vbos --------------------------------------------------------
    for (auto shp : scn->shapes) {
        if (st->vbo.find(shp) == st->vbo.end()) {
            st->vbo[shp] = gl_stdsurface_vbo();
        } else {
            if (refresh.find(shp) == refresh.end()) continue;
        }
        if (!shp->quads_pos.empty()) {
            auto pos = vector<vec3f>();
            auto norm = vector<vec3f>();
            auto texcoord = vector<vec2f>();
            auto quads = vector<vec4i>();
            tie(quads, pos, norm, texcoord) =
                convert_face_varying(shp->quads_pos, shp->quads_norm,
                    shp->quads_texcoord, shp->pos, shp->norm, shp->texcoord);
            update_vertex_buffer(st->vbo[shp].pos, pos);
            update_vertex_buffer(st->vbo[shp].norm, norm);
            update_vertex_buffer(st->vbo[shp].texcoord, texcoord);
            update_element_buffer(
                st->vbo[shp].quads, convert_quads_to_triangles(quads));
            update_element_buffer(
                st->vbo[shp].edges, get_edges({}, {}, shp->quads));
            update_vertex_buffer(st->vbo[shp].color, vector<vec4f>{});
            update_vertex_buffer(st->vbo[shp].tangsp, vector<vec4f>{});
        } else {
            update_vertex_buffer(st->vbo[shp].pos, shp->pos);
            update_vertex_buffer(st->vbo[shp].norm, shp->norm);
            update_vertex_buffer(st->vbo[shp].texcoord, shp->texcoord);
            update_vertex_buffer(st->vbo[shp].color, shp->color);
            update_vertex_buffer(st->vbo[shp].tangsp, shp->tangsp);
            update_element_buffer(st->vbo[shp].points, shp->points);
            update_element_buffer(st->vbo[shp].lines, shp->lines);
            update_element_buffer(st->vbo[shp].triangles, shp->triangles);
            update_element_buffer(
                st->vbo[shp].quads, convert_quads_to_triangles(shp->quads));
            update_element_buffer(
                st->vbo[shp].beziers, convert_bezier_to_lines(shp->beziers));
            update_element_buffer(
                st->vbo[shp].edges, get_edges({}, shp->triangles, shp->quads));
        }
    }
}

// Draw a shape
inline void draw_stdsurface_shape(gl_stdsurface_state* st, const shape* shp,
    const mat4f& xform, bool highlighted, const gl_stdsurface_params& params) {
    static auto default_material = material();
    default_material.kd = {0.2f, 0.2f, 0.2f};

    begin_stdsurface_shape(st->prog, xform);

    auto etype = gl_etype::triangle;
    if (!shp->lines.empty()) etype = gl_etype::line;
    if (!shp->points.empty()) etype = gl_etype::point;

    auto txt = [&st](texture_info& info) -> gl_texture_info {
        if (!info.txt) return {};
        return st->txt.at(info.txt);
    };

    auto mat = (shp->mat) ? shp->mat : &default_material;
    set_stdsurface_material(st->prog, mat->type, etype, mat->ke, mat->kd,
        mat->ks, mat->rs, mat->op, txt(mat->ke_txt), txt(mat->kd_txt),
        txt(mat->ks_txt), txt(mat->rs_txt), txt(mat->norm_txt),
        txt(mat->occ_txt), false, mat->double_sided, params.cutout);

    auto& vbo = st->vbo.at((shape*)shp);
    set_stdsurface_vert(
        st->prog, vbo.pos, vbo.norm, vbo.texcoord, vbo.color, vbo.tangsp);

    draw_elems(vbo.points);
    draw_elems(vbo.lines);
    draw_elems(vbo.triangles);
    draw_elems(vbo.quads);
    draw_elems(vbo.beziers);

    if ((params.edges && !params.wireframe) || highlighted) {
        gl_enable_culling(false);
        set_stdsurface_constmaterial(st->prog,
            (highlighted) ? params.highlight_color : params.edge_color,
            (highlighted) ? 1 : mat->op);
        set_stdsurface_normaloffset(st->prog, params.edge_offset);
        draw_elems(vbo.edges);
        gl_enable_culling(params.cull_backface);
    }

    if (highlighted && false) {
        set_stdsurface_constmaterial(st->prog, params.highlight_color, 1);
        set_stdsurface_normaloffset(st->prog, params.edge_offset);
        draw_elems(vbo.points);
        draw_elems(vbo.lines);
        draw_elems(vbo.triangles);
        draw_elems(vbo.quads);
        draw_elems(vbo.beziers);
    }

    end_stdsurface_shape(st->prog);
}

// Display a scene
void draw_stdsurface_scene(gl_stdsurface_state* st, const scene* scn,
    const camera* view, const gl_stdsurface_lights& lights,
    const gl_stdsurface_params& params) {
    // begin frame
    gl_enable_depth_test(true);
    gl_enable_culling(params.cull_backface);
    gl_enable_wireframe(params.wireframe);
    gl_set_viewport({params.width, params.height});

    auto cam = (params.camera_id < 0) ? view : scn->cameras[params.camera_id];
    auto camera_xform = to_mat(cam->frame);
    auto camera_view = to_mat(inverse(cam->frame));
    auto camera_proj = perspective_mat4(cam->yfov,
        (float)params.width / (float)params.height, cam->near, cam->far);

    begin_stdsurface_frame(st->prog, params.camera_lights, params.exposure,
        params.gamma, params.filmic, camera_xform, camera_view, camera_proj);

    if (!params.camera_lights) {
        set_stdsurface_lights(st->prog, params.ambient, (int)lights.pos.size(),
            lights.pos.data(), lights.ke.data(), lights.ltype.data());
    }

    if (!scn->instances.empty()) {
        for (auto ist : scn->instances) {
            draw_stdsurface_shape(st, ist->shp, to_mat(ist->frame),
                (ist == params.highlighted || ist->shp == params.highlighted),
                params);
        }
    } else {
        for (auto shp : scn->shapes) {
            draw_stdsurface_shape(
                st, shp, identity_mat4f, shp == params.highlighted, params);
        }
    }

    end_stdsurface_frame(st->prog);
    gl_enable_wireframe(false);
}

// Support
inline void _glfw_error_cb(int error, const char* description) {
    printf("GLFW error: %s\n", description);
}

// Support
inline void _glfw_text_cb(GLFWwindow* gwin, unsigned key) {
    auto win = (gl_window*)glfwGetWindowUserPointer(gwin);
    if (win->_widget_enabled) {
        ImGui_ImplGlfwGL3_CharCallback(win->_gwin, key);
    }
    if (win->_text_cb) win->_text_cb(win, key);
}

// Support
inline void _glfw_key_cb(
    GLFWwindow* gwin, int key, int scancode, int action, int mods) {
    auto win = (gl_window*)glfwGetWindowUserPointer(gwin);
    if (win->_widget_enabled) {
        ImGui_ImplGlfwGL3_KeyCallback(win->_gwin, key, scancode, action, mods);
    }
}

// Support
inline void _glfw_mouse_cb(GLFWwindow* gwin, int button, int action, int mods) {
    auto win = (gl_window*)glfwGetWindowUserPointer(gwin);
    if (win->_widget_enabled) {
        ImGui_ImplGlfwGL3_MouseButtonCallback(win->_gwin, button, action, mods);
    }
    if (win->_mouse_cb) win->_mouse_cb(win, button, action == GLFW_PRESS, mods);
}

// Support
inline void _glfw_scroll_cb(GLFWwindow* gwin, double xoffset, double yoffset) {
    auto win = (gl_window*)glfwGetWindowUserPointer(gwin);
    if (win->_widget_enabled) {
        ImGui_ImplGlfwGL3_ScrollCallback(win->_gwin, xoffset, yoffset);
    }
}

// Support
inline void _glfw_refresh_cb(GLFWwindow* gwin) {
    auto win = (gl_window*)glfwGetWindowUserPointer(gwin);
    if (win->_refresh_cb) win->_refresh_cb(win);
}

// Initialize gl_window
gl_window* make_window(
    int width, int height, const string& title, void* user_pointer) {
    auto win = new gl_window();
    // gl_window
    win->_user_pointer = user_pointer;

    // gl_window
    if (!glfwInit()) throw runtime_error("cannot open gl_window");

    // profile creation
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
#if __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

    win->_gwin = glfwCreateWindow(width, height, title.c_str(), 0, 0);
    glfwMakeContextCurrent(win->_gwin);
    glfwSetWindowUserPointer(win->_gwin, win);

    glfwSetErrorCallback(_glfw_error_cb);

    glfwSetCharCallback(win->_gwin, _glfw_text_cb);
    glfwSetKeyCallback(win->_gwin, _glfw_key_cb);
    glfwSetMouseButtonCallback(win->_gwin, _glfw_mouse_cb);
    glfwSetScrollCallback(win->_gwin, _glfw_scroll_cb);

    glfwSetWindowRefreshCallback(win->_gwin, _glfw_refresh_cb);

// init gl extensions
#ifndef __APPLE__
    if (glewInit() != GLEW_OK) return nullptr;
#endif
    return win;
}

// Set gl_window callbacks
void set_window_callbacks(gl_window* win, gl_text_callback text_cb,
    gl_mouse_callback mouse_cb, gl_refresh_callback refresh_cb) {
    win->_text_cb = text_cb;
    win->_mouse_cb = mouse_cb;
    win->_refresh_cb = refresh_cb;
    if (win->_text_cb) glfwSetCharCallback(win->_gwin, _glfw_text_cb);
}

// Clear gl_window
void clear_window(gl_window* win) {
    if (win->_gwin) {
        glfwDestroyWindow(win->_gwin);
        glfwTerminate();
        win->_gwin = nullptr;
    }
    if (win->_widget_enabled) {
        ImGui_ImplGlfwGL3_Shutdown();
        win->_widget_enabled = false;
    }
}

// Set gl_window title
void set_window_title(gl_window* win, const string& title) {
    glfwSetWindowTitle(win->_gwin, title.c_str());
}

// Wait events
void wait_events(gl_window* win) { glfwWaitEvents(); }

// Poll events
void poll_events(gl_window* win) { glfwPollEvents(); }

// Swap buffers
void swap_buffers(gl_window* win) { glfwSwapBuffers(win->_gwin); }

// Should close
bool should_close(gl_window* win) { return glfwWindowShouldClose(win->_gwin); }

// Mouse button
int get_mouse_button(gl_window* win) {
    auto mouse1 =
        glfwGetMouseButton(win->_gwin, GLFW_MOUSE_BUTTON_1) == GLFW_PRESS;
    auto mouse2 =
        glfwGetMouseButton(win->_gwin, GLFW_MOUSE_BUTTON_2) == GLFW_PRESS;
    auto mouse3 =
        glfwGetMouseButton(win->_gwin, GLFW_MOUSE_BUTTON_3) == GLFW_PRESS;
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
    glfwGetCursorPos(win->_gwin, &x, &y);
    return {(int)x, (int)y};
}

// Mouse position
vec2f get_mouse_posf(gl_window* win) {
    double x, y;
    glfwGetCursorPos(win->_gwin, &x, &y);
    return {(float)x, (float)y};
}

// Window size
vec2i get_window_size(gl_window* win) {
    auto ret = vec2i{0, 0};
    glfwGetWindowSize(win->_gwin, &ret.x, &ret.y);
    return ret;
}

// Check if a key is pressed (not all keys are supported)
bool get_key(gl_window* win, int key) {
    key = std::toupper(key);
    return glfwGetKey(win->_gwin, key) == GLFW_PRESS;
}

// Framebuffer size
vec2i get_framebuffer_size(gl_window* win) {
    auto ret = vec2i{0, 0};
    glfwGetFramebufferSize(win->_gwin, &ret.x, &ret.y);
    return ret;
}

// Read pixels
vector<vec4b> get_screenshot(gl_window* win, vec2i& wh, bool flipy, bool back) {
    wh = get_framebuffer_size(win);
    auto pixels = vector<vec4b>(wh[0] * wh[1]);
    glReadBuffer((back) ? GL_BACK : GL_FRONT);
    glReadPixels(0, 0, wh[0], wh[1], GL_RGBA, GL_UNSIGNED_BYTE, pixels.data());
    if (flipy) {
        vector<vec4b> line(wh[0]);
        for (int j = 0; j < wh[1] / 2; j++) {
            memcpy(line.data(), pixels.data() + j * wh[0] * 4, wh[0] * 4);
            memcpy(pixels.data() + j * wh[0] * 4,
                pixels.data() + (wh[1] - 1 - j) * wh[0] * 4, wh[0] * 4);
            memcpy(pixels.data() + (wh[1] - 1 - j) * wh[0] * 4, line.data(),
                wh[0] * 4);
        }
    }
    return pixels;
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
    ImGui_ImplGlfwGL3_Init(win->_gwin, false);
    ImGui::GetStyle().WindowRounding = 0;
    ImGui::GetIO().IniFilename = nullptr;
    ImGui::SetNextWindowPos({0, 0});
    auto size = get_window_size(win);
    ImGui::SetNextWindowSize({(float)win->_widget_width, (float)size[1]});
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
    win->_widget_enabled = true;
}

// Begin draw widget
bool begin_widgets(gl_window* win, const string& title) {
    static bool first_time = true;
    ImGui_ImplGlfwGL3_NewFrame();
    // ImGui::SetNextWindowSize({(float)win->_widget_width, (float)size[1]});
    // ImGui::SetNextWindowPos({(float)(size[0] - win->_widget_width),
    // (float)0});
    // auto flags = ImGuiWindowFlags_NoMove | ImGuiWindowFlags_NoResize;
    // ImGui::Begin(title.c_str(), nullptr, flags);
    if (first_time) {
        ImGui::SetNextWindowPos({0, 0});
        auto size = get_window_size(win);
        ImGui::SetNextWindowSize({(float)win->_widget_width, (float)size[1]});
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
    if (!win->_widget_enabled) return false;
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
void draw_label_widget(gl_window* win, const string& lbl, const string& msg) {
    ImGui::LabelText(lbl.c_str(), "%s", msg.c_str());
}

// Value widget
bool draw_value_widget(gl_window* win, const string& lbl, string& str) {
    char buf[4096];
    if (str.length() >= 4096) throw runtime_error("bad memory");
    memcpy(buf, str.c_str(), str.size());
    buf[str.size()] = 0;
    auto ret = ImGui::InputText(lbl.c_str(), buf, 4096);
    str = buf;
    return ret;
}

// Value widget
bool draw_value_widget(gl_window* win, const string& lbl, int* val, int ncomp,
    int min, int max, int incr) {
    switch (ncomp) {
        case 1: return ImGui::SliderInt(lbl.c_str(), val, min, max);
        case 2: return ImGui::SliderInt2(lbl.c_str(), val, min, max);
        case 3: return ImGui::SliderInt3(lbl.c_str(), val, min, max);
        case 4: return ImGui::SliderInt4(lbl.c_str(), val, min, max);
        default: throw runtime_error("bad number of components"); return false;
    }
}

// Value widget
bool draw_value_widget(gl_window* win, const string& lbl, float* val, int ncomp,
    float min, float max, float incr) {
    switch (ncomp) {
        case 1: return ImGui::SliderFloat(lbl.c_str(), val, min, max);
        case 2: return ImGui::SliderFloat2(lbl.c_str(), val, min, max);
        case 3: return ImGui::SliderFloat3(lbl.c_str(), val, min, max);
        case 4: return ImGui::SliderFloat4(lbl.c_str(), val, min, max);
        default: throw runtime_error("bad number of components"); return false;
    }
}

// Color widget
bool draw_color_widget(gl_window* win, const string& lbl, vec4f& val) {
    return ImGui::ColorEdit4(lbl.c_str(), (float*)&val.x);
}

// Color widget
bool draw_color_widget(gl_window* win, const string& lbl, vec3f& val) {
    return ImGui::ColorEdit3(lbl.c_str(), (float*)&val.x);
}

// Color widget
bool draw_color_widget(gl_window* win, const string& lbl, vec4b& val) {
    auto valf = ImGui::ColorConvertU32ToFloat4(*(uint32_t*)&val);
    if (ImGui::ColorEdit4(lbl.c_str(), &valf.x)) {
        auto valb = ImGui::ColorConvertFloat4ToU32(valf);
        *(uint32_t*)&val = valb;
        return true;
    }
    return false;
}

// Support
inline bool _enum_widget_labels_ptr(void* data, int idx, const char** out) {
    auto labels = (vector<pair<string, int>>*)data;
    *out = labels->at(idx).first.c_str();
    return true;
}

// Support
inline bool _enum_widget_labels_int(void* data, int idx, const char** out) {
    auto labels = (vector<pair<string, int>>*)data;
    *out = labels->at(idx).first.c_str();
    return true;
}

// Support
inline bool _enum_widget_labels_str(void* data, int idx, const char** out) {
    auto labels = (vector<string>*)data;
    *out = labels->at(idx).c_str();
    return true;
}

// Value widget
bool draw_value_widget(gl_window* win, const string& lbl, string& val,
    const vector<string>& labels) {
    auto cur = -1;
    for (auto idx = 0; idx < labels.size(); idx++) {
        if (labels[idx] == val) cur = idx;
    }
    assert(cur >= 0);
    auto ok = ImGui::Combo(lbl.c_str(), &cur, _enum_widget_labels_str,
        (void*)&labels, (int)labels.size());
    val = labels[cur];
    return ok;
}

// Value widget
bool draw_value_widget(gl_window* win, const string& lbl, int& val,
    const vector<pair<string, int>>& labels) {
    auto cur = -1;
    for (auto idx = 0; idx < labels.size(); idx++) {
        if (labels[idx].second == val) cur = idx;
    }
    assert(cur >= 0);
    auto ok = ImGui::Combo(lbl.c_str(), &cur, _enum_widget_labels_int,
        (void*)&labels, (int)labels.size());
    val = labels[cur].second;
    return ok;
}

// Value widget
bool draw_value_widget(gl_window* win, const string& lbl, void*& val,
    const vector<pair<string, void*>>& labels) {
    auto cur = -1;
    for (auto idx = 0; idx < labels.size(); idx++) {
        if (labels[idx].second == val) cur = idx;
    }
    assert(cur >= 0);
    auto ok = ImGui::Combo(lbl.c_str(), &cur, _enum_widget_labels_int,
        (void*)&labels, (int)labels.size());
    val = labels[cur].second;
    return ok;
}

// Bool widget
bool draw_value_widget(gl_window* win, const string& lbl, bool& val) {
    return ImGui::Checkbox(lbl.c_str(), &val);
}

// Button widget
bool draw_button_widget(gl_window* win, const string& lbl) {
    return ImGui::Button(lbl.c_str());
}

// Collapsible header
bool draw_header_widget(gl_window* win, const string& lbl) {
    return ImGui::CollapsingHeader(lbl.c_str());
}

// Start tree node
bool draw_tree_widget_begin(gl_window* win, const string& lbl) {
    return ImGui::TreeNode(lbl.c_str());
}

// Collapsible header
void draw_tree_widget_end(gl_window* win) { ImGui::TreePop(); }

// Start selectable tree node
bool draw_tree_widget_begin(
    gl_window* win, const string& lbl, void*& selection, void* content) {
    ImGuiTreeNodeFlags node_flags =
        ImGuiTreeNodeFlags_OpenOnArrow | ImGuiTreeNodeFlags_OpenOnDoubleClick;
    if (selection == content) node_flags |= ImGuiTreeNodeFlags_Selected;
    auto open = ImGui::TreeNodeEx(content, node_flags, "%s", lbl.c_str());
    if (ImGui::IsItemClicked()) selection = content;
    return open;
}

// Start selectable tree node
bool draw_tree_widget_begin(gl_window* win, const string& lbl, void*& selection,
    void* content, const vec4f& col) {
    ImGui::PushStyleColor(ImGuiCol_Text, {col.x, col.y, col.z, col.w});
    auto ret = draw_tree_widget_begin(win, lbl, selection, content);
    ImGui::PopStyleColor();
    return ret;
}

// End selectable tree node
void draw_tree_widget_end(gl_window* win, void* content) { ImGui::TreePop(); }

// Selectable tree leaf node
void draw_tree_widget_leaf(
    gl_window* win, const string& lbl, void*& selection, void* content) {
    ImGuiTreeNodeFlags node_flags =
        ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen;
    if (selection == content) node_flags |= ImGuiTreeNodeFlags_Selected;
    ImGui::TreeNodeEx(content, node_flags, "%s", lbl.c_str());
    if (ImGui::IsItemClicked()) selection = content;
}

// Selectable tree leaf node
void draw_tree_widget_leaf(gl_window* win, const string& lbl, void*& selection,
    void* content, const vec4f& col) {
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

// Scroll region
void draw_scroll_widget_begin(
    gl_window* win, const string& lbl, int height, bool border) {
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
void draw_groupid_widget_end(gl_window* win) { ImGui::PopID(); }

// Text color
void draw_tree_widget_color_begin(gl_window* win, const vec4f& color) {
    ImGui::PushStyleColor(
        ImGuiCol_Text, {color[0], color[1], color[2], color[3]});
}

// Text color
void draw_tree_widget_color_end(gl_window* win) { ImGui::PopStyleColor(); }

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR SIMPLE SCENE UI
// -----------------------------------------------------------------------------
namespace ygl {

inline void draw_tree_widgets(
    gl_window* win, const string& lbl, camera* cam, void*& selection) {
    draw_tree_widget_leaf(win, lbl + cam->name, selection, cam);
}

inline void draw_tree_widgets(
    gl_window* win, const string& lbl, texture* txt, void*& selection) {
    draw_tree_widget_leaf(win, lbl + txt->path, selection, txt);
}

inline void draw_tree_widgets(
    gl_window* win, const string& lbl, texture_info* info, void*& selection) {
    draw_tree_widgets(win, lbl, info->txt, selection);
}

inline void draw_tree_widgets(
    gl_window* win, const string& lbl, material* mat, void*& selection) {
    if (draw_tree_widget_begin(win, lbl + mat->name, selection, mat)) {
        if (mat->ke_txt.txt)
            draw_tree_widgets(win, "ke: ", mat->ke_txt.txt, selection);
        if (mat->kd_txt.txt)
            draw_tree_widgets(win, "kd: ", mat->kd_txt.txt, selection);
        if (mat->ks_txt.txt)
            draw_tree_widgets(win, "ks: ", mat->ks_txt.txt, selection);
        if (mat->rs_txt.txt)
            draw_tree_widgets(win, "rs: ", mat->rs_txt.txt, selection);
        if (mat->norm_txt.txt)
            draw_tree_widgets(win, "norm: ", mat->norm_txt.txt, selection);
        if (mat->bump_txt.txt)
            draw_tree_widgets(win, "bump: ", mat->bump_txt.txt, selection);
        if (mat->disp_txt.txt)
            draw_tree_widgets(win, "disp: ", mat->disp_txt.txt, selection);
        draw_tree_widget_end(win);
    }
}

inline void draw_tree_widgets(
    gl_window* win, const string& lbl, shape* shp, void*& selection) {
    if (draw_tree_widget_begin(win, lbl + shp->name, selection, shp)) {
        if (shp->mat) draw_tree_widgets(win, "mat: ", shp->mat, selection);
        draw_tree_widget_end(win);
    }
}

inline void draw_tree_widgets(
    gl_window* win, const string& lbl, instance* ist, void*& selection) {
    if (draw_tree_widget_begin(win, lbl + ist->name, selection, ist)) {
        if (ist->shp) draw_tree_widgets(win, "shp: ", ist->shp, selection);
        draw_tree_widget_end(win);
    }
}

inline void draw_tree_widgets(
    gl_window* win, const string& lbl, environment* env, void*& selection) {
    if (draw_tree_widget_begin(win, lbl + env->name, selection, env)) {
        if (env->ke_txt.txt)
            draw_tree_widgets(win, "txt: ", env->ke_txt.txt, selection);
        draw_tree_widget_end(win);
    }
}

inline void draw_tree_widgets(
    gl_window* win, const string& lbl, node* nde, void*& selection) {
    if (draw_tree_widget_begin(win, lbl + nde->name, selection, nde)) {
        if (nde->cam) draw_tree_widgets(win, "cam: ", nde->cam, selection);
        for (auto idx = 0; idx < nde->ists.size(); idx++)
            draw_tree_widgets(
                win, "ist" + to_string(idx) + ": ", nde->ists[idx], selection);
        if (nde->env) draw_tree_widgets(win, "env: ", nde->env, selection);
        if (nde->parent)
            draw_tree_widgets(win, "parent", nde->parent, selection);
        for (auto idx = 0; idx < nde->children_.size(); idx++) {
            draw_tree_widgets(win, "", nde->children_[idx], selection);
        }
        draw_tree_widget_end(win);
    }
}

inline void draw_tree_widgets(
    gl_window* win, const string& lbl, animation* anm, void*& selection) {
    draw_tree_widget_leaf(win, lbl + anm->name, selection, anm);
}

template <typename T>
inline void draw_scene_tree_widgets(gl_window* win, const string& lbl,
    const vector<T*>& elems, void*& selection) {
    if (draw_tree_widget_begin(win, lbl)) {
        for (auto elem : elems) draw_tree_widgets(win, "", elem, selection);
        draw_tree_widget_end(win);
    }
}

inline void draw_tree_widgets(
    gl_window* win, const string& lbl, scene* scn, void*& selection) {
    draw_scene_tree_widgets(win, lbl + "cameras", scn->cameras, selection);
    draw_scene_tree_widgets(win, lbl + "textures", scn->textures, selection);
    draw_scene_tree_widgets(win, lbl + "materials", scn->materials, selection);
    draw_scene_tree_widgets(win, lbl + "shapes", scn->shapes, selection);
    draw_scene_tree_widgets(win, lbl + "instances", scn->instances, selection);
    draw_scene_tree_widgets(
        win, lbl + "environments", scn->environments, selection);
    draw_scene_tree_widgets(win, lbl + "nodes", scn->nodes, selection);
    draw_scene_tree_widgets(
        win, lbl + "animations", scn->animations, selection);
}

inline bool draw_elem_widgets(gl_window* win, scene* scn, texture* txt,
    void*& selection, const unordered_map<texture*, gl_texture>& gl_txt) {
    draw_separator_widget(win);
    draw_label_widget(win, "path", txt->path);
    auto size = format("{} x {} @ 4 {}", txt->width(), txt->height(),
        (txt->ldr) ? "byte" : "float");
    draw_label_widget(win, "size", size);
    if (contains(gl_txt, txt)) {
        draw_image_widget(win, get_texture_id(gl_txt.at(txt)), {128, 128},
            {txt->width(), txt->height()});
    }
    return false;
}

inline bool draw_elem_widgets(gl_window* win, scene* scn, material* mat,
    void*& selection, const unordered_map<texture*, gl_texture>& gl_txt) {
    auto txt_names = vector<pair<string, texture*>>{{"<none>", nullptr}};
    for (auto txt : scn->textures) txt_names.push_back({txt->path, txt});

    auto edited = vector<bool>();
    draw_separator_widget(win);
    draw_label_widget(win, "name", mat->name);
    edited += draw_value_widget(win, "type", mat->type, material_type_names());
    auto ke_l = max_element_value(mat->ke);
    auto ke_c = (ke_l) ? mat->ke / ke_l : zero3f;
    edited += draw_value_widget(win, "ke l", ke_l, 0, 100);
    edited += draw_color_widget(win, "ke", ke_c);
    mat->ke = ke_l * ke_c;
    edited += draw_color_widget(win, "kd", mat->kd);
    edited += draw_color_widget(win, "ks", mat->ks);
    edited += draw_color_widget(win, "kt", mat->kt);
    edited += draw_value_widget(win, "rs", mat->rs);

    auto txt_widget = [&txt_names](gl_window* win, const string& lbl,
                          texture_info& info) {
        auto edited = vector<bool>();
        edited += draw_value_widget(win, lbl, info.txt, txt_names);
        if (info.txt) {
            edited += draw_value_widget(win, lbl + " wrap_s", info.wrap_s);
            edited += draw_value_widget(win, lbl + " wrap_t", info.wrap_t);
            edited += draw_value_widget(win, lbl + " linear", info.linear);
            edited += draw_value_widget(win, lbl + " mipmap", info.mipmap);
        }
        return std::any_of(
            edited.begin(), edited.end(), [](auto x) { return x; });
    };
    edited += txt_widget(win, "ke_txt", mat->ke_txt);
    edited += txt_widget(win, "kd_txt", mat->kd_txt);
    edited += txt_widget(win, "ks_txt", mat->ks_txt);
    edited += txt_widget(win, "kt_txt", mat->kt_txt);
    edited += txt_widget(win, "norm_txt", mat->norm_txt);
    edited += txt_widget(win, "bump_txt", mat->bump_txt);
    edited += txt_widget(win, "disp_txt", mat->disp_txt);
    return std::any_of(edited.begin(), edited.end(), [](auto x) { return x; });
}

inline bool draw_elem_widgets(gl_window* win, scene* scn, shape* shp,
    void*& selection, const unordered_map<texture*, gl_texture>& gl_txt) {
    auto mat_names = vector<pair<string, material*>>{{"<none>", nullptr}};
    for (auto mat : scn->materials) mat_names.push_back({mat->name, mat});

    auto edited = vector<bool>();
    draw_separator_widget(win);
    draw_label_widget(win, "name", shp->name);
    edited += draw_value_widget(win, "material", shp->mat, mat_names);
    draw_label_widget(win, "pos", shp->pos, true);
    draw_label_widget(win, "norm", shp->norm, true);
    draw_label_widget(win, "texcoord", shp->texcoord, true);
    draw_label_widget(win, "color", shp->color, true);
    draw_label_widget(win, "tangsp", shp->tangsp, true);
    draw_label_widget(win, "radius", shp->radius, true);
    draw_label_widget(win, "triangles", shp->triangles, true);
    draw_label_widget(win, "quads", shp->quads, true);
    draw_label_widget(win, "quads_pos", shp->quads_pos, true);
    draw_label_widget(win, "quads_norm", shp->quads_norm, true);
    draw_label_widget(win, "quads_texcoord", shp->quads_texcoord, true);
    draw_label_widget(win, "lines", shp->lines, true);
    draw_label_widget(win, "points", shp->points, true);
    draw_label_widget(win, "beziers", shp->beziers, true);
    return std::any_of(edited.begin(), edited.end(), [](auto x) { return x; });
}

inline bool draw_elem_widgets(gl_window* win, scene* scn, camera* cam,
    void*& selection, const unordered_map<texture*, gl_texture>& gl_txt) {
    auto edited = vector<bool>();
    draw_separator_widget(win);
    draw_label_widget(win, "name", cam->name);
    edited += draw_value_widget(win, "frame", cam->frame, -10, 10);
    edited += draw_value_widget(win, "yfov", cam->yfov, 0.1, 4);
    edited += draw_value_widget(win, "aspect", cam->aspect, 0.1, 4);
    edited += draw_value_widget(win, "focus", cam->focus, 0.01, 10);
    edited += draw_value_widget(win, "aperture", cam->aperture);
    return std::any_of(edited.begin(), edited.end(), [](auto x) { return x; });
}

inline bool draw_elem_widgets(gl_window* win, scene* scn, instance* ist,
    void*& selection, const unordered_map<texture*, gl_texture>& gl_txt) {
    auto shp_names = vector<pair<string, shape*>>{{"<none>", nullptr}};
    for (auto shp : scn->shapes) shp_names.push_back({shp->name, shp});

    auto edited = vector<bool>();
    draw_separator_widget(win);
    draw_label_widget(win, "name", ist->name);
    edited += draw_value_widget(win, "frame", ist->frame, -10, 10);
    edited += draw_value_widget(win, "shape", ist->shp, shp_names);
    return std::any_of(edited.begin(), edited.end(), [](auto x) { return x; });
}

inline bool draw_elem_widgets(gl_window* win, scene* scn, environment* env,
    void*& selection, const unordered_map<texture*, gl_texture>& gl_txt) {
    auto txt_names = vector<pair<string, texture*>>{{"<none>", nullptr}};
    for (auto txt : scn->textures) txt_names.push_back({txt->path, txt});

    auto txt_widget = [&txt_names](gl_window* win, const string& lbl,
                          texture_info& info) {
        auto edited = vector<bool>();
        edited += draw_value_widget(win, lbl, info.txt, txt_names);
        if (info.txt) {
            edited += draw_value_widget(win, lbl + " wrap_s", info.wrap_s);
            edited += draw_value_widget(win, lbl + " wrap_t", info.wrap_t);
            edited += draw_value_widget(win, lbl + " linear", info.linear);
            edited += draw_value_widget(win, lbl + " mipmap", info.mipmap);
        }
        return std::any_of(
            edited.begin(), edited.end(), [](auto x) { return x; });
    };

    auto edited = vector<bool>();
    draw_separator_widget(win);
    draw_label_widget(win, "name", env->name);
    edited += draw_value_widget(win, "frame", env->frame, 0, 0);
    auto ke_l = max_element_value(env->ke);
    auto ke_c = (ke_l) ? env->ke / ke_l : zero3f;
    edited += draw_value_widget(win, "ke l", ke_l, 0, 100);
    edited += draw_color_widget(win, "ke", ke_c);
    env->ke = ke_l * ke_c;
    edited += txt_widget(win, "ke_txt", env->ke_txt);
    return std::any_of(edited.begin(), edited.end(), [](auto x) { return x; });
}

inline bool draw_elem_widgets(gl_window* win, scene* scn, node* nde,
    void*& selection, const unordered_map<texture*, gl_texture>& gl_txt) {
    auto cam_names = vector<pair<string, camera*>>{{"<none>", nullptr}};
    for (auto cam : scn->cameras) cam_names.push_back({cam->name, cam});
    auto ist_names = vector<pair<string, instance*>>{{"<none>", nullptr}};
    for (auto ist : scn->instances) ist_names.push_back({ist->name, ist});
    auto env_names = vector<pair<string, environment*>>{{"<none>", nullptr}};
    for (auto env : scn->environments) env_names.push_back({env->name, env});
    auto nde_names = vector<pair<string, node*>>{{"<none>", nullptr}};
    for (auto nde : scn->nodes) nde_names.push_back({nde->name, nde});

    auto edited = vector<bool>();
    draw_separator_widget(win);
    draw_label_widget(win, "name", nde->name);
    edited += draw_value_widget(win, "parent", nde->parent, nde_names);
    edited += draw_value_widget(win, "frame", nde->frame, -10, 10);
    edited += draw_value_widget(win, "camera", nde->cam, cam_names);
    for (auto idx = 0; idx < nde->ists.size(); idx++)
        edited += draw_value_widget(
            win, "instance " + to_string(idx), nde->ists[idx], ist_names);
    edited += draw_value_widget(win, "environment", nde->env, env_names);
    for (auto idx = 0; idx < nde->children_.size(); idx++)
        edited += draw_value_widget(
            win, "child " + to_string(idx), nde->children_[idx], nde_names);
    return std::any_of(edited.begin(), edited.end(), [](auto x) { return x; });
}

inline bool draw_elem_widgets(gl_window* win, scene* scn, animation* anm,
    void*& selection, const unordered_map<texture*, gl_texture>& gl_txt) {
    auto nde_names = vector<pair<string, node*>>{{"<none>", nullptr}};
    for (auto nde : scn->nodes) nde_names.push_back({nde->name, nde});

    auto edited = vector<bool>();
    draw_separator_widget(win);
    draw_label_widget(win, "name", anm->name);
    for (auto kfr_kv : enumerate(anm->keyframes)) {
        auto kfr = kfr_kv.second;
        auto ids = to_string(kfr_kv.first);
        edited += draw_value_widget(win, "name " + ids, kfr->name);
        edited += draw_value_widget(
            win, "type " + ids, kfr->type, keyframe_type_names());
        draw_label_widget(win, "times " + ids, kfr->times, true);
        draw_label_widget(win, "translation " + ids, kfr->translation, true);
        draw_label_widget(win, "rotation " + ids, kfr->rotation, true);
        draw_label_widget(win, "scale " + ids, kfr->scaling, true);
    }
    for (auto target_kv : enumerate(anm->targets)) {
        auto kfr = target_kv.second.first;
        auto nde = target_kv.second.second;
        auto ids = to_string(target_kv.first);
        draw_label_widget(win, "target " + ids, kfr->name + " -> " + nde->name);
    }
    return std::any_of(edited.begin(), edited.end(), [](auto x) { return x; });
}

inline bool draw_elem_widgets(gl_window* win, test_scene_params* scn,
    test_texture_params* txt, void*& selection,
    const unordered_map<texture*, gl_texture>& gl_txt) {
    draw_separator_widget(win);
    draw_label_widget(win, "name", txt->name);
    draw_value_widget(win, "type", txt->type, test_texture_names());
    draw_value_widget(win, "resolution", txt->resolution);
    draw_value_widget(win, "tile size", txt->tile_size);
    draw_value_widget(win, "noise size", txt->noise_scale);
    draw_value_widget(win, "sky sun angle", txt->sky_sunangle);
    draw_value_widget(win, "normal map", txt->bump_to_normal);
    draw_value_widget(win, "bump scale", txt->bump_scale);
#if 0
        if (contains(gl_txt, txt)) {
            draw_image_widget(win, get_texture_id(gl_txt.at(txt)), {128, 128},
                              {txt->width(), txt->height()});
        }
#endif
    return false;
}

inline bool draw_elem_widgets(gl_window* win, test_scene_params* scn,
    test_material_params* mat, void*& selection,
    const unordered_map<texture*, gl_texture>& gl_txt) {
    auto edited = vector<bool>();
    draw_separator_widget(win);
    edited += draw_value_widget(win, "name", mat->name);
    edited += draw_value_widget(win, "type", mat->type, test_material_names());
    edited += draw_value_widget(win, "emission", mat->emission, 0, 100);
    edited += draw_color_widget(win, "color", mat->color);
    edited += draw_value_widget(win, "roughness", mat->roughness);
    edited += draw_value_widget(win, "texture", mat->texture);
    edited += draw_value_widget(win, "normal", mat->normal);
    return std::any_of(edited.begin(), edited.end(), [](auto x) { return x; });
}

inline bool draw_elem_widgets(gl_window* win, test_scene_params* scn,
    test_shape_params* shp, void*& selection,
    const unordered_map<texture*, gl_texture>& gl_txt) {
    auto edited = vector<bool>();
    draw_separator_widget(win);
    edited += draw_value_widget(win, "name", shp->name);
    edited += draw_value_widget(win, "type", shp->type, test_shape_names());
    edited += draw_value_widget(win, "material", shp->material);
    edited += draw_value_widget(win, "tesselation", shp->tesselation, -1, 8);
    edited += draw_value_widget(win, "subdivision", shp->subdivision, -1, 8);
    edited += draw_value_widget(win, "scale", shp->scale, 0.01f, 10.0f);
    edited += draw_value_widget(win, "faceted", shp->faceted);
    edited += draw_value_widget(win, "num lines/points", shp->num, -1, 100000);
    edited += draw_value_widget(
        win, "radius lines/points", shp->radius, 0.0001f, 0.01f);
    return std::any_of(edited.begin(), edited.end(), [](auto x) { return x; });
}

inline bool draw_elem_widgets(gl_window* win, test_scene_params* scn,
    test_camera_params* cam, void*& selection,
    const unordered_map<texture*, gl_texture>& gl_txt) {
    auto edited = vector<bool>();
    draw_separator_widget(win);
    edited += draw_value_widget(win, "name", cam->name);
    edited += draw_value_widget(win, "from", cam->from, -10, 10);
    edited += draw_value_widget(win, "to", cam->to, -10, 10);
    edited += draw_value_widget(win, "yfov", cam->yfov, 0.1, 4);
    edited += draw_value_widget(win, "aspect", cam->aspect, 0.1, 4);
    return std::any_of(edited.begin(), edited.end(), [](auto x) { return x; });
}

inline bool draw_elem_widgets(gl_window* win, test_scene_params* scn,
    test_environment_params* env, void*& selection,
    const unordered_map<texture*, gl_texture>& gl_txt) {
    auto edited = vector<bool>();
    draw_separator_widget(win);
    edited += draw_value_widget(win, "name", env->name);
    edited += draw_value_widget(win, "rotation", env->rotation, -pif, pif);
    edited += draw_value_widget(win, "emission", env->emission);
    edited += draw_value_widget(win, "txt", env->texture);
    return std::any_of(edited.begin(), edited.end(), [](auto x) { return x; });
}

inline bool draw_elem_widgets(gl_window* win, test_scene_params* scn,
    test_node_params* nde, void*& selection,
    const unordered_map<texture*, gl_texture>& gl_txt) {
    auto edited = vector<bool>();
    draw_separator_widget(win);
    edited += draw_value_widget(win, "name", nde->name);
    edited += draw_value_widget(win, "frame", nde->frame);
    edited += draw_value_widget(win, "camera", nde->camera);
    edited += draw_value_widget(win, "parent", nde->parent);
    edited += draw_value_widget(win, "instance", nde->instance);
    edited += draw_value_widget(win, "environment", nde->environment);
    return std::any_of(edited.begin(), edited.end(), [](auto x) { return x; });
}

inline bool draw_elem_widgets(gl_window* win, test_scene_params* scn,
    test_animation_params* anm, void*& selection,
    const unordered_map<texture*, gl_texture>& gl_txt) {
    auto edited = vector<bool>();
    draw_separator_widget(win);
    edited += draw_value_widget(win, "name", anm->name);
    edited += draw_value_widget(win, "bezier", anm->bezier);
    edited += draw_value_widget(win, "speed", anm->speed);
    edited += draw_value_widget(win, "scale", anm->scale);
    for (auto idx = 0; idx < anm->times.size(); idx++)
        edited +=
            draw_value_widget(win, "time " + to_string(idx), anm->times[idx]);
    for (auto idx = 0; idx < anm->translation.size(); idx++)
        edited += draw_value_widget(
            win, "translation " + to_string(idx), anm->translation[idx]);
    for (auto idx = 0; idx < anm->rotation.size(); idx++)
        edited += draw_value_widget(
            win, "rotation " + to_string(idx), anm->rotation[idx]);
    for (auto idx = 0; idx < anm->scaling.size(); idx++)
        edited += draw_value_widget(
            win, "scaling " + to_string(idx), anm->scaling[idx]);
    for (auto idx = 0; idx < anm->nodes.size(); idx++)
        edited +=
            draw_value_widget(win, "node " + to_string(idx), anm->nodes[idx]);
    return std::any_of(edited.begin(), edited.end(), [](auto x) { return x; });
}

inline bool draw_elem_widgets(gl_window* win, test_scene_params* scn,
    test_instance_params* ist, void*& selection,
    const unordered_map<texture*, gl_texture>& gl_txt) {
    auto edited = vector<bool>();
    draw_separator_widget(win);
    edited += draw_value_widget(win, "name", ist->name);
    edited += draw_value_widget(win, "frame", ist->frame, -10, 10);
    edited += draw_value_widget(win, "shape", ist->shape);
    return std::any_of(edited.begin(), edited.end(), [](auto x) { return x; });
}

template <typename T, typename T1>
inline bool draw_selected_elem_widgets(gl_window* win, scene* scn,
    test_scene_params* test_scn, const vector<T*>& elems,
    vector<T1>& test_elems, void*& selection,
    void (*update_test_elem)(const scene*, T*, const T1&),
    const unordered_map<texture*, gl_texture>& gl_txt) {
    auto selected = (T*)nullptr;
    for (auto elem : elems)
        if (elem == selection) {
            selected = elem;
            break;
        }
    if (!selected) return false;
    auto test_selected = (T1*)nullptr;
    for (auto& test_elem : test_elems)
        if (test_elem.name == selected->name) test_selected = &test_elem;
    auto edited = vector<bool>();
    if (test_selected) {
        edited +=
            draw_elem_widgets(win, test_scn, test_selected, selection, gl_txt);
        if (edited.back()) update_test_elem(scn, selected, *test_selected);
    }
    edited += draw_elem_widgets(win, scn, selected, selection, gl_txt);
    return std::any_of(edited.begin(), edited.end(), [](auto x) { return x; });
}

template <typename T, typename T1>
inline bool draw_add_elem_widgets(gl_window* win, scene* scn, const string& lbl,
    vector<T*>& elems, vector<T1>& test_elems, void*& selection,
    void (*update_test_elem)(const scene*, T*, const T1&)) {
    static auto count = 0;
    if (draw_button_widget(win, "add " + lbl)) {
        auto name = lbl + "_" + to_string(count++);
        elems += new T();
        elems.back()->name = name;
        test_elems += T1();
        test_elems.back().name = name;
        selection = elems.back();
        update_test_elem(scn, elems.back(), test_elems.back());
        return true;
    }
    draw_continue_widget(win);
    return false;
}

bool draw_scene_widgets(gl_window* win, const string& lbl, scene* scn,
    void*& selection, const unordered_map<texture*, gl_texture>& gl_txt,
    test_scene_params* test_scn) {
    static auto test_scn_def = test_scene_params();

    if (draw_header_widget(win, lbl)) {
        // draw_scroll_widget_begin(win, "model", 240, false);
        draw_tree_widgets(win, "", scn, selection);
        // draw_scroll_widget_end(win);

        auto edited = vector<bool>();

        if (test_scn) {
            edited += draw_add_elem_widgets(win, scn, "cam", scn->cameras,
                test_scn->cameras, selection, update_test_camera);
            edited += draw_add_elem_widgets(win, scn, "txt", scn->textures,
                test_scn->textures, selection, update_test_texture);
            edited += draw_add_elem_widgets(win, scn, "mat", scn->materials,
                test_scn->materials, selection, update_test_material);
            edited += draw_add_elem_widgets(win, scn, "shp", scn->shapes,
                test_scn->shapes, selection, update_test_shape);
            if (edited.back()) {
                scn->instances += new instance();
                scn->instances.back()->name = scn->shapes.back()->name;
                scn->instances.back()->shp = scn->shapes.back();
                test_scn->instances += test_instance_params();
                test_scn->instances.back().name = scn->instances.back()->name;
                test_scn->instances.back().shape =
                    scn->instances.back()->shp->name;
            }
            edited += draw_add_elem_widgets(win, scn, "ist", scn->instances,
                test_scn->instances, selection, update_test_instance);
            edited += draw_add_elem_widgets(win, scn, "env", scn->environments,
                test_scn->environments, selection, update_test_environment);
            edited += draw_add_elem_widgets(win, scn, "anim", scn->animations,
                test_scn->animations, selection, update_test_animation);
        }

        auto test_scn_res = (test_scn) ? test_scn : &test_scn_def;
        edited += draw_selected_elem_widgets(win, scn, test_scn, scn->cameras,
            test_scn_res->cameras, selection, update_test_camera, gl_txt);
        edited += draw_selected_elem_widgets(win, scn, test_scn, scn->textures,
            test_scn_res->textures, selection, update_test_texture, gl_txt);
        edited += draw_selected_elem_widgets(win, scn, test_scn, scn->materials,
            test_scn_res->materials, selection, update_test_material, gl_txt);
        edited += draw_selected_elem_widgets(win, scn, test_scn, scn->shapes,
            test_scn_res->shapes, selection, update_test_shape, gl_txt);
        edited += draw_selected_elem_widgets(win, scn, test_scn, scn->instances,
            test_scn_res->instances, selection, update_test_instance, gl_txt);
        edited += draw_selected_elem_widgets(win, scn, test_scn,
            scn->environments, test_scn->environments, selection,
            update_test_environment, gl_txt);
        edited += draw_selected_elem_widgets(win, scn, test_scn, scn->nodes,
            test_scn->nodes, selection, update_test_node, gl_txt);
        edited +=
            draw_selected_elem_widgets(win, scn, test_scn, scn->animations,
                test_scn->animations, selection, update_test_animation, gl_txt);
        return std::any_of(
            edited.begin(), edited.end(), [](auto x) { return x; });
    } else
        return false;
}

#if 0
        inline void draw_edit_widgets(gl_window* win, scene* scn,
            void*&  selection, const yshade_state* state) {
            static auto shape_names =
                vector<pair<string, int>>{{"cube", 0}, {"sphere", 1}};
            static auto shape_type = 0;
            static char txt_filename[1024] = "grid.png";

            auto selected_ist = (instance*)nullptr;
            auto selected_shp = (shape*)nullptr;
            auto selected_cam = (camera*)nullptr;
            auto selected_txt = (texture*)nullptr;
            auto selected_mat = (material*)nullptr;

            if (*selection) {
                for (auto ptr : scn->instances)
                    if (ptr == *selection) selected_node = ptr;
                for (auto ptr : scn->scenes)
                    if (ptr == *selection) selected_scene = ptr;
                for (auto ptr : scn->meshes)
                    if (ptr == *selection) selected_mesh = ptr;
                for (auto ptr : scn->cameras)
                    if (ptr == *selection) selected_cam = ptr;
                for (auto ptr : scn->materials)
                    if (ptr == *selection) selected_mat = ptr;
                for (auto ptr : scn->textures)
                    if (ptr == *selection) selected_txt = ptr;
            }

            static auto auto_parent = true;
            draw_value_widget(win, "set parent from selection", &auto_parent);

            if (draw_button_widget(win, "add mesh")) {
                static auto count = 0;
                auto mesh = new ygltf::mesh();
                mesh->name = "<new mesh " + to_string(count++) + ">";
                auto shp = new ygltf::shape();
                mesh->shapes.push_back(shp);
                shp->name = "<new shape " + to_string(count - 1) + ">";
                switch (shape_type) {
                    case 0: {
                        make_uvcube(
                            1, 1, shp->triangles, shp->pos, shp->norm, shp->texcoord);
                    } break;
                    case 1: {
                        make_uvsphere(
                            1, 1, shp->triangles, shp->pos, shp->norm, shp->texcoord);
                    } break;
                }
                if (auto_parent && selected_node) selected_node->msh = mesh;
                gscn->meshes.push_back(mesh);
                *selection = mesh;
            }
            draw_value_widget(win, "shape type", &shape_type, shape_names);

            if (draw_button_widget(win, "add camera")) {
                static auto count = 0;
                auto cam = new ygltf::camera();
                cam->name = "<new camera " + to_string(count++) + ">";
                if (auto_parent && selected_node) selected_node->cam = cam;
                gscn->cameras.push_back(cam);
                *selection = cam;
            }

            if (draw_button_widget(win, "add node")) {
                static auto count = 0;
                auto node = new ygltf::node();
                node->name = "<new node " + to_string(count++) + ">";
                if (auto_parent && selected_node)
                    selected_node->children.push_back(node);
                gscn->nodes.push_back(node);
                *selection = node;
            }

            if (draw_button_widget(win, "add texture")) {
                static auto count = 0;
                auto txt = new ygltf::texture();
                txt->name = "<new texture " + to_string(count++) + ">";
                txt->path = txt_filename;
                auto scn = (app_state*)get_user_pointer(win,);
                auto dirname = yu::path::get_dirname(scn->filename);
                try {
                    if (yimg::is_hdr_filename(txt->path)) {
                        txt->hdr = yimg::load_image4f(dirname + txt->path);
                    } else {
                        txt->ldr = yimg::load_image4b(dirname + txt->path);
                    }
                } catch (...) { txt->ldr = image4b(1, 1, {255, 255, 255, 255}); }
                gscn->textures.push_back(txt);
                *selection = txt;
            }
            draw_text_widget(win, "texture", txt_filename, sizeof(txt_filename));

            if (draw_button_widget(win, "delete")) {
                if (selected_cam) {
                    for (auto node : gscn->nodes)
                        if (node->cam == selected_cam) node->cam = nullptr;
                    remove(gscn->cameras, selected_cam);
                    delete selected_cam;
                    *selection = nullptr;
                }
                if (selected_mesh) {
                    for (auto node : gscn->nodes)
                        if (node->msh == selected_mesh) node->msh = nullptr;
                    remove(gscn->meshes, selected_mesh);
                    delete selected_mesh;
                    *selection = nullptr;
                }
                if (selected_mat) {
                    for (auto mesh : gscn->meshes)
                        for (auto shp : mesh->shapes)
                            if (shp->mat == selected_mat) shp->mat = nullptr;
                    remove(gscn->materials, selected_mat);
                    delete selected_mat;
                    *selection = nullptr;
                }
                if (selected_txt) {
                    for (auto mat : gscn->materials) {
                        if (mat->emission_txt == selected_txt)
                            mat->emission_txt = nullptr;
                        if (mat->normal_txt == selected_txt) mat->normal_txt = nullptr;
                        if (mat->occlusion_txt == selected_txt)
                            mat->occlusion_txt = nullptr;
                        if (mat->metallic_roughness) {
                            if (mat->metallic_roughness->base_txt == selected_txt)
                                mat->metallic_roughness->base_txt = nullptr;
                            if (mat->metallic_roughness->metallic_txt == selected_txt)
                                mat->metallic_roughness->metallic_txt = nullptr;
                        }
                        if (mat->specular_glossiness) {
                            if (mat->specular_glossiness->diffuse_txt == selected_txt)
                                mat->specular_glossiness->diffuse_txt = nullptr;
                            if (mat->specular_glossiness->specular_txt == selected_txt)
                                mat->specular_glossiness->specular_txt = nullptr;
                        }
                    }
                    remove(gscn->textures, selected_txt);
                    delete selected_txt;
                    *selection = nullptr;
                }
            }
        }
#endif

}  // namespace ygl

#endif

// HACK to avoid compilation with MSVC2015 without dirtying code
#ifdef constexpr
#undef constexpr
#endif
