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
// IMPLEMENTATION FOR MONETACARLO SAMPLING FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {
// Make a one-dimensional distribution for sampling.
distribution1f make_distribution(const std::vector<float>& weights) {
    auto dist = distribution1f();
    dist.weights = weights;
    dist.cdf = weights;
    for (auto i = 1; i < weights.size(); i++) dist.cdf[i] += dist.cdf[i - 1];
    return dist;
}

// Sample a discrete distribution.
int sample_distribution_discrete(const distribution1f& dist, float r) {
    // todo: implement binary search better
    r = clamp(r * dist.cdf.back(), 0.0f, dist.cdf.back() - 0.00001f);
    for (auto i = 0; i < dist.cdf.size(); i++) {
        if (dist.cdf[i] > r) return i;
    }
    return (int)dist.cdf.size() - 1;
}

// Sample a discrete distribution.
float sample_distribution_discrete_pdf(const distribution1f& dist, int idx) {
    if (idx == 0) return dist.cdf.at(0);
    return 1 / (dist.cdf.at(idx) - dist.cdf.at(idx - 1));
}

// Get the total weight of a distribution.
float sample_distribution_weightsum(const distribution1f& dist) {
    return dist.cdf.back();
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// FILE LOADING AND SAVING
// -----------------------------------------------------------------------------
namespace ygl {

// Loads the contents of a binary file in an in-memory array.
std::vector<unsigned char> load_binary(const std::string& filename) {
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
}

// Loads the contents of a text file into a string.
std::string load_text(const std::string& filename) {
    auto buf = load_binary(filename);
#ifdef _WIN32
    auto nbuf = std::vector<unsigned char>();
    nbuf.reserve(buf.size());
    for (auto i = 0; i < buf.size(); i++) {
        if (buf[i] == '\r') continue;
        nbuf.push_back(buf[i]);
    }
    buf = nbuf;
#endif
    return std::string((char*)buf.data(), (char*)buf.data() + buf.size());
}

// Saves binary data to a file.
void save_binary(
    const std::string& filename, const std::vector<unsigned char>& data) {
    auto f = fopen(filename.c_str(), "wb");
    if (!f) throw std::runtime_error("cannot write file " + filename);
    auto num = fwrite(data.data(), 1, data.size(), f);
    if (num != data.size())
        throw std::runtime_error("cannot write file " + filename);
    fclose(f);
}

// Saves a string to a text file.
void save_text(const std::string& filename, const std::string& str) {
    auto f = fopen(filename.c_str(), "wb");
    if (!f) throw std::runtime_error("cannot write file " + filename);
    auto num = fwrite(str.c_str(), 1, str.size(), f);
    if (num != str.size())
        throw std::runtime_error("cannot write file " + filename);
    fclose(f);
}

// Load INI file. The implementation does not handle escaping.
std::unordered_map<std::string, std::unordered_map<std::string, std::string>>
load_ini(const std::string& filename) {
    auto txt = load_text(filename);
    auto lines = splitlines(txt, false);
    auto ret = std::unordered_map<std::string,
        std::unordered_map<std::string, std::string>>();
    auto cur_group = ""s;
    ret[""] = {};
    for (auto line : lines) {
        line = strip(line);
        if (line.empty()) continue;
        if (line.front() == ';') continue;
        if (line.front() == '#') continue;
        if (line.front() == '[') {
            if (line.back() != ']') throw std::runtime_error("bad INI format");
            cur_group = line.substr(1, line.length() - 2);
            ret[cur_group] = {};
        } else if (line.find('=') != line.npos) {
            auto var = line.substr(0, line.find('='));
            auto val = line.substr(line.find('=') + 1);
            ret[cur_group][var] = val;
        } else {
            throw std::runtime_error("bad INI format");
        }
    }
    return ret;
}

// Save INI file. The implementation does not handle escaping.
void save_ini(const std::string& filename,
    std::unordered_map<std::string,
        std::unordered_map<std::string, std::string>>& values) {
    auto txt = ""s;
    for (auto& gkv : values) {
        txt += "[" + gkv.first + "]\n\n";
        for (auto& vkv : gkv.second) txt += vkv.first + "=" + vkv.second + "\n";
        txt += "\n";
    }
    save_text(filename, txt);
}

}  // namespace ygl

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
void merge_lines(std::vector<vec2i>& lines, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord,
    const std::vector<vec2i>& lines1, const std::vector<vec3f>& pos1,
    const std::vector<vec3f>& norm1, const std::vector<vec2f>& texcoord1) {
    lines.reserve(lines.size() + lines1.size());
    auto nverts = (int)pos.size();
    for (auto& l : lines1) lines.push_back({l.x + nverts, l.y + nverts});
    append(pos, pos1);
    append(norm, norm1);
    append(texcoord, texcoord1);
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
    append(pos, pos1);
    append(norm, norm1);
    append(texcoord, texcoord1);
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
    append(pos, pos1);
    append(norm, norm1);
    append(texcoord, texcoord1);
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
    auto dst = make_triangle_distribution(triangles, pos);
    auto rng = make_rng(seed);
    for (auto i = 0; i < npoints; i++) {
        auto eid = 0;
        auto euv = zero2f;
        std::tie(eid, euv) = sample_triangles(
            dst, next_rand1f(rng), {next_rand1f(rng), next_rand1f(rng)});
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
        auto pc = vec3f{abs(pos[i].x), abs(pos[i].y), abs(pos[i].z)};
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
        auto a = vec2f{2 * pif * uv.x, pif * (1 - uv.y)};
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
        auto phi = 2 * pif * uv.x;
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
        auto phi = 2 * pif * uv.x;
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
        auto pc = vec2f{pcyl.y, abs(pcyl.z)};
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

// Make a seashell. This is not watertight. Returns quads, pos, norm,
// texcoord.
void make_uvseashell(std::vector<vec4i>& quads, std::vector<vec3f>& pos,
    std::vector<vec3f>& norm, std::vector<vec2f>& texcoord, const vec2i& steps,
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

    make_quad(quads, pos, norm, texcoord, steps, {1, 1}, {1, 1});
    for (auto i = 0; i < pos.size(); i++) {
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
        texcoord[i] = vec2f{uv.x * params.uvsize.x, uv.y * params.uvsize.y * R};
    }

    compute_normals(quads, pos, norm);
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
    std::tie(bpos, bnorm, btexcoord) = sample_triangles_points(
        join(striangles, convert_quads_to_triangles(squads)), spos, snorm,
        stexcoord, steps.y, params.seed);

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
            pos[i] = lerp(pos[i], pos[i + (cidx[bidx] - bidx) * (steps.x + 1)],
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
// IMPLEMENTATION FOR IMAGEIO
// -----------------------------------------------------------------------------
namespace ygl {

// Pfm load
float* load_pfm(const char* filename, int* w, int* h, int* nc, int req) {
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
    return make_image(w, h, (vec4b*)pixels.get());
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
    } else if (path_extension(filename) == ".pfm") {
        return save_pfm(
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
    auto ext = path_extension(filename);
    auto pixels = (float*)nullptr;
    if (ext == ".exr") {
        if (LoadEXR(&pixels, &width, &height, filename.c_str(), nullptr) < 0)
            return {};
        ncomp = 4;
    } else if(ext == ".pfm") {
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
    } else if (path_extension(filename) == ".pfm") {
        return save_pfm(filename.c_str(), width, height, ncomp, hdr);
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
    return (x * (6.2f * x + vec3f{0.5f})) /
           (x * (6.2f * x + vec3f{1.7f}) + vec3f{0.06f});
}

inline vec3f tonemap_filmic2(const vec3f& hdr) {
    // https://knarkowicz.wordpress.com/2016/01/06/aces-filmic-tone-mapping-curve/
    auto x = hdr;
    // x *= 0.6; // brings it back to ACES range
    x = (x * (2.51f * x + vec3f{0.03f})) /
        (x * (2.43f * x + vec3f{0.59f}) + vec3f{0.14f});
    return tonemap_gamma(x);
}

inline vec3f tonemap_filmic3(const vec3f& hdr) {
    // https://github.com/TheRealMJP/BakingLab/blob/master/BakingLab/ACES.hlsl

    // sRGB => XYZ => D65_2_D60 => AP1 => RRT_SAT
    static const mat3f ACESInputMat = transpose(mat3f(
        vec3f(0.59719, 0.35458, 0.04823), vec3f(0.07600, 0.90834, 0.01566),
        vec3f(0.02840, 0.13383, 0.83777)));

    // ODT_SAT => XYZ => D60_2_D65 => sRGB
    static const mat3f ACESOutputMat = transpose(mat3f(
        vec3f(1.60475, -0.53108, -0.07367), vec3f(-0.10208, 1.10813, -0.00605),
        vec3f(-0.00327, -0.07276, 1.07602)));

    auto x = hdr;
    x = 2 * x;  // matches standard range
    x = ACESInputMat * x;
    // Apply RRT and ODT
    vec3f a = x * (x + vec3f{0.0245786f}) - vec3f{0.000090537f};
    vec3f b = x * (0.983729f * x + vec3f{0.4329510f}) + vec3f{0.238081f};
    x = a / b;
    x = ACESOutputMat * x;
    return tonemap_gamma(x);
}

// Tone mapping HDR to LDR images.
image4b tonemap_image(const image4f& hdr, const tonemap_params& params) {
    auto ldr = image4b(hdr.width(), hdr.height());
    auto scale = pow(2.0f, params.exposure);
    for (auto j = 0; j < hdr.height(); j++) {
        for (auto i = 0; i < hdr.width(); i++) {
            auto h4 = hdr.at(i, j);
            auto h = vec3f{h4.x, h4.y, h4.z} * scale;
            auto a = h4.w;
            switch (params.type) {
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
        auto centroid_size = bbox_size(centroid_bbox);

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

            // check whether we were able to split
            if (mid > start && mid < end) {
                // makes an internal node
                node.type = bvh_node_type::internal;
                // perform the splits by preallocating children and recurring
                node.axis = axis;
                node.start = (int)nodes.size();
                node.count = 2;
                nodes.emplace_back();
                nodes.emplace_back();
                // build child nodes
                make_bvh_node(nodes, node.start, sorted_prims, start, mid,
                    bboxes, type, equal_size);
                make_bvh_node(nodes, node.start + 1, sorted_prims, mid, end,
                    bboxes, type, equal_size);
            }
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
            auto sbvh = bvh->shape_bvhs[ist.sid];
            bboxes.push_back(transform_bbox(ist.frame, sbvh->nodes[0].bbox));
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
                auto sbvh = bvh->shape_bvhs[ist.sid];
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
                    auto sbvh = bvh->shape_bvhs[ist.sid];
                    if (intersect_bvh(sbvh, transform_ray(ist.frame_inv, ray),
                            find_any, ray_t, iid, sid, eid, euv)) {
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
                    auto sbvh = bvh->shape_bvhs[ist.sid];
                    if (overlap_bvh(sbvh, transform_point(ist.frame_inv, pos),
                            max_dist, find_any, dist, iid, sid, eid, euv)) {
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

// Check if a shape has emission.
bool has_emission(const shape* shp) {
    for (auto grp : shp->groups) {
        auto mat = grp.mat;
        if (mat && mat->ke != zero3f) return true;
    }
    return false;
}
// Gets the material for a shape element.
material* get_material(const shape* shp, int eid) {
    if (shp->groups.empty()) return nullptr;
    if (shp->group_ids.empty()) return shp->groups.at(0).mat;
    return shp->groups.at(shp->group_ids[eid]).mat;
}
// Returns is the shape is simple.
bool is_shape_simple(const shape* shp, bool split_facevarying) {
    if (split_facevarying && !shp->quads_pos.empty()) return false;
    if (!shp->group_ids.empty() && shp->groups.size() > 1) return false;
    return true;
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
vec3f eval_pos(const shape* shp, int eid, const vec2f& euv, bool transformed) {
    auto pos = eval_elem(shp, shp->pos, eid, euv, {0, 0, 0});
    if (transformed) pos = transform_point(shp->frame, pos);
    return pos;
}
// Shape normal interpolated using barycentric coordinates
vec3f eval_norm(const shape* shp, int eid, const vec2f& euv, bool transformed) {
    auto norm = normalize(eval_elem(shp, shp->norm, eid, euv, {0, 0, 1}));
    if (transformed) norm = transform_direction(shp->frame, norm);
    return norm;
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
        vec3f{uv.x * pif * 2, uv.y * pif, environment_distance});
    if (transformed) pos = transform_point(env->frame, pos);
    return pos;
}
// Environment normal interpolated using uv parametrization.
vec3f eval_norm(const environment* env, const vec2f& uv, bool transformed) {
    auto norm = normalize(
        -sphericaly_to_cartesian(vec3f{uv.x * pif * 2, uv.y * pif, 1.0f}));
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
    return {sh.x / (2 * pif), sh.y / pif};
}

// Evaluate a texture
vec4f eval_texture(const texture_info& info, const vec2f& texcoord, bool srgb,
    const vec4f& def) {
    auto txt = info.txt;
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
        vec2f{(ij.x + puv.x) / (cam->aspect * res), 1 - (ij.y - puv.y) / res};
    return eval_camera_ray(cam, uv, luv);
}

// Generate a distribution for sampling a shape uniformly based on area/length.
distribution1f make_shape_distribution(const shape* shp) {
    switch (get_shape_type(shp)) {
        case shape_elem_type::none: {
            throw std::runtime_error("type not supported");
        } break;
        case shape_elem_type::triangles: {
            return make_triangle_distribution(shp->triangles, shp->pos);
        } break;
        case shape_elem_type::lines: {
            return make_line_distribution(shp->lines, shp->pos);
        } break;
        case shape_elem_type::points: {
            return make_point_distribution(shp->points.size());
        } break;
        case shape_elem_type::quads: {
            return make_quad_distribution(shp->quads, shp->pos);
        } break;
        case shape_elem_type::beziers: {
            throw std::runtime_error("type not supported");
        } break;
        case shape_elem_type::vertices: {
            return make_point_distribution(shp->pos.size());
        } break;
        case shape_elem_type::facevarying: {
            return make_quad_distribution(shp->quads_pos, shp->pos);
        } break;
    }
}

// Sample a shape based on a distribution.
std::pair<int, vec2f> sample_shape(
    const shape* shp, const distribution1f& dst, float re, const vec2f& ruv) {
    switch (get_shape_type(shp)) {
        case shape_elem_type::none: {
            throw std::runtime_error("type not supported");
        } break;
        case shape_elem_type::triangles: {
            return sample_triangles(dst, re, ruv);
        } break;
        case shape_elem_type::lines: {
            return {sample_lines(dst, re, ruv.x).first, ruv};
        } break;
        case shape_elem_type::points: {
            return {sample_points(dst, re), ruv};
        } break;
        case shape_elem_type::quads: {
            return sample_quads(dst, re, ruv);
        } break;
        case shape_elem_type::beziers: {
            throw std::runtime_error("type not supported");
        } break;
        case shape_elem_type::vertices: {
            return {sample_points(dst, re), ruv};
        } break;
        case shape_elem_type::facevarying: {
            return sample_quads(dst, re, ruv);
        } break;
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

// Facet a shape. Supports only non-facevarying shapes
void facet_shape(shape* shp, bool recompute_normals) {
    switch (get_shape_type(shp)) {
        case shape_elem_type::none: break;
        case shape_elem_type::points: break;
        case shape_elem_type::vertices: break;
        case shape_elem_type::triangles: {
            facet_triangles(shp->triangles, shp->pos, shp->norm, shp->texcoord,
                shp->color, shp->radius);
        } break;
        case shape_elem_type::lines: {
            facet_lines(shp->lines, shp->pos, shp->norm, shp->texcoord,
                shp->color, shp->radius);
        } break;
        case shape_elem_type::quads: {
            std::vector<int> verts;
            facet_quads(shp->quads, shp->pos, shp->norm, shp->texcoord,
                shp->color, shp->radius);
        } break;
        case shape_elem_type::beziers: {
            throw std::runtime_error("type not supported");
        } break;
        case shape_elem_type::facevarying: {
            throw std::runtime_error("type not supported");
        } break;
    }
    if (recompute_normals) compute_normals(shp);
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

int add_shape_vert(shape* shp, const vec<int, 5>& vert,
    std::unordered_map<vec<int, 5>, int>& vmap, const std::vector<vec3f>& pos,
    const std::vector<vec3f>& norm, const std::vector<vec2f>& texcoord,
    const std::vector<vec4f>& color, const std::vector<float>& radius) {
    if (contains(vmap, vert)) return vmap.at(vert);
    auto vid = (int)vmap.size();
    vmap[vert] = vid;
    if (!pos.empty() && vert[0] >= 0) shp->pos.push_back(pos[vert[0]]);
    if (!norm.empty() && vert[1] >= 0) shp->norm.push_back(norm[vert[1]]);
    if (!texcoord.empty() && vert[2] >= 0)
        shp->texcoord.push_back(texcoord[vert[2]]);
    if (!color.empty() && vert[3] >= 0) shp->color.push_back(color[vert[3]]);
    if (!radius.empty() && vert[4] >= 0) shp->radius.push_back(radius[vert[4]]);
    return vid;
}

int add_shape_vert(shape* shp, int vert, std::unordered_map<int, int>& vmap,
    const std::vector<vec3f>& pos, const std::vector<vec3f>& norm,
    const std::vector<vec2f>& texcoord, const std::vector<vec4f>& color,
    const std::vector<float>& radius) {
    if (contains(vmap, vert)) return vmap.at(vert);
    auto vid = (int)vmap.size();
    vmap[vert] = vid;
    if (!pos.empty()) shp->pos.push_back(pos[vert]);
    if (!norm.empty()) shp->norm.push_back(norm[vert]);
    if (!texcoord.empty()) shp->texcoord.push_back(texcoord[vert]);
    if (!color.empty()) shp->color.push_back(color[vert]);
    if (!radius.empty()) shp->radius.push_back(radius[vert]);
    return vid;
}

// Convert face-varying shapes to shapes.
std::vector<shape*> split_shape(const shape* shp, bool split_facevarying) {
    if (is_shape_simple(shp, split_facevarying)) return {};
    auto sshps = std::vector<shape*>();
    for (auto gid = 0; gid < shp->groups.size(); gid++) {
        auto grp = shp->groups[gid];
        auto sshp = new shape();
        sshp->name = (grp.name == "") ? shp->name : grp.name;
        sshp->groups.push_back(grp);
        auto vmap = std::unordered_map<int, int>();
        auto add_vert = [&shp, &sshp, &vmap](int vert) {
            return add_shape_vert(sshp, vert, vmap, shp->pos, shp->norm,
                shp->texcoord, shp->color, shp->radius);
        };
        auto add_fvvert = [&vmap](auto& svals, auto& vals, int vert) {
            if (vals.empty()) return -1;
            if (contains(vmap, vert)) return vmap.at(vert);
            auto vid = (int)vmap.size();
            vmap[vert] = vid;
            svals.push_back(vals[vert]);
            return vid;
        };
        for (auto p : shp->points) sshp->points.push_back(add_vert(p));
        for (auto l : shp->lines)
            sshp->lines.push_back({add_vert(l.x), add_vert(l.y)});
        for (auto t : shp->triangles)
            sshp->triangles.push_back(
                {add_vert(t.x), add_vert(t.y), add_vert(t.z)});
        for (auto q : shp->quads)
            sshp->quads.push_back(
                {add_vert(q.x), add_vert(q.y), add_vert(q.z), add_vert(q.w)});
        for (auto b : shp->beziers)
            sshp->beziers.push_back(
                {add_vert(b.x), add_vert(b.y), add_vert(b.z), add_vert(b.w)});
        if (split_facevarying) {
            auto fvmap = std::unordered_map<vec<int, 5>, int>();
            for (auto eid = 0; eid < shp->quads_pos.size(); eid++) {
                sshp->quads.push_back({});
                for (auto c = 0; c < 4; c++) {
                    auto vert = vec<int, 5>{shp->quads_pos[eid][c],
                        (shp->quads_norm.empty()) ? -1 :
                                                    shp->quads_norm[eid][c],
                        (shp->quads_texcoord.empty()) ?
                            -1 :
                            shp->quads_texcoord[eid][c],
                        -1, -1};
                    sshp->quads.back()[c] = add_shape_vert(sshp, vert, fvmap,
                        shp->pos, shp->norm, shp->texcoord, {}, {});
                }
            }
        } else {
            vmap.clear();
            for (auto q : shp->quads_pos) {
                sshp->quads_pos.push_back({add_fvvert(sshp->pos, shp->pos, q.x),
                    add_fvvert(sshp->pos, shp->pos, q.y),
                    add_fvvert(sshp->pos, shp->pos, q.z),
                    add_fvvert(sshp->pos, shp->pos, q.w)});
            }
            vmap.clear();
            for (auto q : shp->quads_norm) {
                sshp->quads_norm.push_back(
                    {add_fvvert(sshp->norm, shp->norm, q.x),
                        add_fvvert(sshp->norm, shp->norm, q.y),
                        add_fvvert(sshp->norm, shp->norm, q.z),
                        add_fvvert(sshp->norm, shp->norm, q.w)});
            }
            vmap.clear();
            for (auto q : shp->quads_texcoord) {
                sshp->quads_texcoord.push_back(
                    {add_fvvert(sshp->texcoord, shp->texcoord, q.x),
                        add_fvvert(sshp->texcoord, shp->texcoord, q.y),
                        add_fvvert(sshp->texcoord, shp->texcoord, q.z),
                        add_fvvert(sshp->texcoord, shp->texcoord, q.w)});
            }
        }
        sshps.push_back(sshp);
    }
    return sshps;
}

// Convert a list of shapes into one shape.
shape* merge_shapes(const std::vector<shape*>& shps, bool pad_vertex) {
    if (shps.size() < 2) return nullptr;
    auto gshp = new shape();
    gshp->name = shps.front()->name;
    gshp->frame = shps.front()->frame;
    auto pad_norm = false, pad_texcoord = false, pad_texcoord1 = false,
         pad_color = false, pad_radius = false, pad_tangsp = false;
    for (auto i = 1; i < shps.size(); i++) {
        pad_norm = shps[i]->norm.empty() != shps[i - 1]->norm.empty();
        pad_texcoord =
            shps[i]->texcoord.empty() != shps[i - 1]->texcoord.empty();
        pad_texcoord1 =
            shps[i]->texcoord1.empty() != shps[i - 1]->texcoord1.empty();
        pad_color = shps[i]->color.empty() != shps[i - 1]->color.empty();
        pad_radius = shps[i]->radius.empty() != shps[i - 1]->radius.empty();
        pad_tangsp = shps[i]->tangsp.empty() != shps[i - 1]->tangsp.empty();
        if (shps[i]->frame != shps[i - 1]->frame)
            log_warning("frames are not corrected during shape merge");
    }
    if (!pad_vertex && (pad_norm || pad_texcoord || pad_texcoord1 ||
                           pad_color || pad_radius || pad_tangsp))
        throw std::runtime_error("need padding");
    for (auto shp : shps) {
        if (shp->pos.empty()) continue;
        auto goffset = (uint16_t)gshp->groups.size();
        auto voffset = (int)gshp->pos.size();
        if (shp->groups.empty()) {
            gshp->groups.push_back({});
        } else {
            append(gshp->groups, shp->groups);
        }
        append(gshp->pos, shp->pos);
        if (shp->norm.empty() && pad_norm)
            append(gshp->norm, shp->pos.size(), zero3f);
        else
            append(gshp->norm, shp->norm);
        if (shp->texcoord.empty() && pad_texcoord)
            append(gshp->texcoord, shp->pos.size(), zero2f);
        else
            append(gshp->texcoord, shp->texcoord);
        if (shp->texcoord1.empty() && pad_texcoord1)
            append(gshp->texcoord1, shp->pos.size(), zero2f);
        else
            append(gshp->texcoord1, shp->texcoord1);
        if (shp->color.empty() && pad_color)
            append(gshp->color, shp->pos.size(), {1, 1, 1, 1});
        else
            append(gshp->color, shp->color);
        if (shp->radius.empty() && pad_radius)
            append(gshp->radius, shp->pos.size(), 0.0f);
        else
            append(gshp->radius, shp->radius);
        if (shp->tangsp.empty() && pad_tangsp)
            append(gshp->tangsp, shp->pos.size(), {0, 0, 1, 1});
        else
            append(gshp->tangsp, shp->tangsp);
        for (auto p : shp->points) gshp->points.push_back(p + voffset);
        for (auto l : shp->lines) gshp->lines.push_back(l + vec2i{voffset});
        for (auto t : shp->triangles)
            gshp->triangles.push_back(t + vec3i{voffset});
        for (auto q : shp->quads) gshp->quads.push_back(q + vec4i{voffset});
        for (auto b : shp->beziers) gshp->beziers.push_back(b + vec4i{voffset});
        if (shp->group_ids.empty()) {
            auto nelems = (int)max({shp->points.size(), shp->lines.size(),
                shp->triangles.size(), shp->quads.size(), shp->beziers.size()});
            for (auto i = 0; i < nelems; i++)
                gshp->group_ids.push_back(goffset);
        } else {
            for (auto gid : shp->group_ids)
                gshp->group_ids.push_back((uint16_t)(goffset + gid));
        }
    }
    return gshp;
}

// Convert a list of shapes into one shape. Pad shape values if needed.
void merge_into(shape* shp, const shape* shp1, bool pad_vertex) {
    auto mshp = merge_shapes({shp, (shape*)shp1}, pad_vertex);
    *shp = *mshp;
    delete mshp;
}

// Update animation transforms
void update_transforms(const animation_group* agr, float time) {
    auto interpolate = [](keyframe_type type, const std::vector<float>& times,
                           const auto& vals, float time) {
        switch (type) {
            case keyframe_type::step:
                return eval_keyframed_step(times, vals, time);
            case keyframe_type::linear:
                return eval_keyframed_linear(times, vals, time);
            case keyframe_type::bezier:
                return eval_keyframed_bezier(times, vals, time);
            default: throw std::runtime_error("should not have been here");
        }
        return vals.at(0);
    };

    for (auto& anm_ : agr->animations) {
        auto anm = &anm_;
        if (!anm->translation.empty()) {
            auto val =
                interpolate(anm->type, anm->times, anm->translation, time);
            for (auto target : anm->targets) target->translation = val;
        }
        if (!anm->rotation.empty()) {
            auto val = interpolate(anm->type, anm->times, anm->rotation, time);
            for (auto target : anm->targets) target->rotation = val;
        }
        if (!anm->scaling.empty()) {
            auto val = interpolate(anm->type, anm->times, anm->scaling, time);
            for (auto target : anm->targets) target->scaling = val;
        }
    }
}

// Update node transforms
void update_transforms(node* nde, const frame3f& parent = identity_frame3f) {
    nde->frame_ = parent * nde->local * translation_frame(nde->translation) *
                  rotation_frame(nde->rotation) * scaling_frame(nde->scaling);
    if (nde->shp) nde->shp->frame = nde->frame_;
    if (nde->cam) nde->cam->frame = nde->frame_;
    if (nde->env) nde->env->frame = nde->frame_;
    for (auto child : nde->children_) update_transforms(child, nde->frame_);
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
        for (auto& anm_ : agr->animations) {
            auto anm = &anm_;
            range.x = min(range.x, anm->times.front());
            range.y = max(range.y, anm->times.back());
        }
    }
    return range;
}

// Add missing values and elements
void add_elements(scene* scn, const add_elements_options& opts) {
    if (opts.smooth_normals) {
        for (auto shp : scn->shapes) {
            if (!shp->norm.empty()) continue;
            auto type = get_shape_type(shp);
            if (type == shape_elem_type::facevarying) {
                if (!shp->quads_norm.empty())
                    throw std::runtime_error("bad normals");
                shp->quads_norm = shp->quads_pos;
                compute_normals(shp->quads_pos, shp->pos, shp->norm);
            } else if (type == shape_elem_type::beziers) {
                shp->norm.resize(shp->pos.size(), {0, 0, 1});
            } else {
                compute_normals(shp);
            }
        }
    }

    if (opts.tangent_space) {
        for (auto shp : scn->shapes) {
            if (!shp->tangsp.empty() || shp->texcoord.empty()) continue;
            auto mat_needs_tangents = false;
            for (auto grp : shp->groups) {
                if (grp.mat && (grp.mat->norm_txt.txt || grp.mat->bump_txt.txt))
                    mat_needs_tangents = true;
            }
            if (!mat_needs_tangents) continue;
            auto type = get_shape_type(shp);
            if (type == shape_elem_type::triangles) {
                compute_tangent_frames(shp->triangles, shp->pos, shp->norm,
                    shp->texcoord, shp->tangsp);
            } else if (type == shape_elem_type::quads) {
                auto triangles = convert_quads_to_triangles(shp->quads);
                compute_tangent_frames(
                    triangles, shp->pos, shp->norm, shp->texcoord, shp->tangsp);
            } else {
                throw std::runtime_error("type not supported");
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

    if (opts.node_hierarchy && scn->nodes.empty()) {
        for (auto cam : scn->cameras) {
            auto nde = new node();
            nde->name = cam->name;
            nde->local = cam->frame;
            nde->cam = cam;
            scn->nodes.push_back(nde);
        }
        for (auto shp : scn->shapes) {
            auto nde = new node();
            nde->name = shp->name;
            nde->local = shp->frame;
            nde->shp = shp;
            scn->nodes.push_back(nde);
        }
        for (auto env : scn->environments) {
            auto nde = new node();
            nde->name = env->name;
            nde->local = env->frame;
            nde->env = env;
            scn->nodes.push_back(nde);
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

        auto eid = 0;
        for (auto env : scn->environments) {
            if (env->name.empty())
                env->name = "unnamed_environment_" + std::to_string(eid);
            eid++;
        }

        auto nid = 0;
        for (auto nde : scn->nodes) {
            if (nde->name.empty())
                nde->name = "unnamed_node_" + std::to_string(nid);
            nid++;
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
bbox3f compute_bounds(const shape* shp) {
    auto bbox = invalid_bbox3f;
    for (auto p : shp->pos) bbox += transform_point(shp->frame, p);
    return bbox;
}

// Updates the scene and scene's instances bounding boxes
bbox3f compute_bounds(const scene* scn, bool skip_emitting) {
    if (scn->nodes.empty()) {
        auto bbox = invalid_bbox3f;
        for (auto shp : scn->shapes) {
            if (skip_emitting && has_emission(shp)) continue;
            bbox += compute_bounds(shp);
        }
        return bbox;
    } else {
        auto shape_bboxes = std::unordered_map<shape*, bbox3f>();
        for (auto shp : scn->shapes) {
            if (skip_emitting && has_emission(shp)) continue;
            shape_bboxes[shp] += compute_bounds(shp);
        }
        auto bbox = invalid_bbox3f;
        for (auto nde : scn->nodes) {
            if (skip_emitting && has_emission(nde->shp)) continue;
            bbox += transform_bbox(nde->frame_ * inverse(nde->shp->frame),
                shape_bboxes.at(nde->shp));
        }
        return bbox;
    }
}

// Flatten scene instances into separate meshes.
void flatten_instances(scene* scn) {
    auto shape_count = std::unordered_map<shape*, int>();
    for (auto shp : scn->shapes) shape_count[shp] = 0;
    for (auto nde : scn->nodes) {
        if (!nde->shp) continue;
        shape_count[nde->shp] += 1;
        if (shape_count[nde->shp] < 2) continue;
        auto shp = new shape(*nde->shp);
        shp->name += "_" + std::to_string(shape_count[nde->shp]);
        scn->shapes.push_back(shp);
        nde->shp = shp;
    }
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
    for (auto sid = 0; sid < scn->shapes.size(); sid++) {
        auto shp = scn->shapes[sid];
        shape_bvhs.push_back(make_bvh(shp, def_radius, equalsize));
        smap[shp] = sid;
    }

    // tree bvh
    auto bists = std::vector<bvh_instance>();
    if (scn->nodes.empty()) {
        for (auto iid = 0; iid < scn->shapes.size(); iid++) {
            auto shp = scn->shapes[iid];
            auto bist = bvh_instance();
            bist.frame = shp->frame;
            bist.frame_inv = inverse(shp->frame);
            bist.iid = iid;
            bist.sid = smap.at(shp);
            bists.push_back(bist);
        }
    } else {
        for (auto iid = 0; iid < scn->nodes.size(); iid++) {
            auto nde = scn->nodes[iid];
            if (!nde->shp) continue;
            auto bist = bvh_instance();
            bist.frame = nde->frame_;
            bist.frame_inv = inverse(nde->frame_);
            bist.iid = iid;
            bist.sid = smap.at(nde->shp);
            bists.push_back(bist);
        }
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
            refit_bvh(
                get_shape_bvhs(bvh).at(sid), shp->pos, shp->radius, def_radius);
            sid++;
        }
    }
    auto ist_frames = std::vector<frame3f>();
    auto ist_frames_inv = std::vector<frame3f>();
    if (scn->nodes.empty()) {
        for (auto sgr : scn->shapes) {
            ist_frames.push_back(sgr->frame);
            ist_frames_inv.push_back(inverse(sgr->frame));
        }
    } else {
        for (auto nde : scn->nodes) {
            if (!nde->shp) continue;
            ist_frames.push_back(nde->frame_);
            ist_frames_inv.push_back(inverse(nde->frame_));
        }
    }
    refit_bvh(bvh, ist_frames, ist_frames_inv);
}

// Compute stats.
scene_stats compute_stats(const scene* scn) {
    auto stats = scene_stats();
    stats.num_cameras = scn->cameras.size();
    stats.num_shapes = scn->shapes.size();
    stats.num_materials = scn->materials.size();
    stats.num_textures = scn->textures.size();
    stats.num_environments = scn->environments.size();
    stats.num_nodes = scn->nodes.size();
    stats.num_animation_groups = scn->animations.size();

    for (auto nde : scn->nodes) {
        if (nde->shp) stats.num_instances += 1;
    }
    for (auto shp : scn->shapes) {
        stats.elem_points += shp->points.size();
        stats.elem_lines += shp->lines.size();
        stats.elem_triangles += shp->triangles.size();
        stats.elem_quads += shp->quads.size();
        stats.vert_pos += shp->pos.size();
        stats.vert_norm += shp->norm.size();
        stats.vert_texcoord += shp->texcoord.size();
        stats.vert_color += shp->color.size();
        stats.vert_radius += shp->radius.size();
        stats.vert_tangsp += shp->tangsp.size();
    }
    stats.memory_elems =
        stats.elem_points * sizeof(int) + stats.elem_lines * sizeof(vec2i) +
        stats.elem_triangles * sizeof(vec3i) + stats.elem_quads * sizeof(vec4i);
    stats.memory_verts =
        stats.vert_pos * sizeof(vec3f) + stats.vert_norm * sizeof(vec3f) +
        stats.vert_texcoord * sizeof(vec3f) + stats.vert_color * sizeof(vec4f) +
        stats.vert_tangsp * sizeof(vec4f) + stats.vert_radius * sizeof(float);

    for (auto agr : scn->animations) {
        stats.num_animations += agr->animations.size();
    }

    for (auto txt : scn->textures) {
        stats.texel_ldrs = txt->ldr.width() * txt->ldr.height();
        stats.texel_hdrs = txt->hdr.width() * txt->hdr.height();
    }
    stats.memory_ldrs = stats.texel_ldrs * sizeof(vec4b);
    stats.memory_hdrs = stats.texel_hdrs * sizeof(vec4f);

    stats.bbox_scn = compute_bounds(scn);
    stats.bbox_nolights = compute_bounds(scn, true);

    return stats;
}

// Stream out
std::ostream& operator<<(std::ostream& os, const scene_stats& stats) {
    os << "num_cameras: " << stats.num_cameras << "\n";
    os << "num_shape_groups: " << stats.num_shape_groups << "\n";
    os << "num_shapes: " << stats.num_shapes << "\n";
    os << "num_instances: " << stats.num_instances << "\n";
    os << "num_materials: " << stats.num_materials << "\n";
    os << "num_textures: " << stats.num_textures << "\n";
    os << "num_environments: " << stats.num_environments << "\n";
    os << "num_nodes: " << stats.num_nodes << "\n";
    os << "num_animation_groups: " << stats.num_animation_groups << "\n";
    os << "num_animations: " << stats.num_animations << "\n";
    os << "elem_points: " << stats.elem_points << "\n";
    os << "elem_lines: " << stats.elem_lines << "\n";
    os << "elem_triangles: " << stats.elem_triangles << "\n";
    os << "elem_quads: " << stats.elem_quads << "\n";
    os << "vert_pos: " << stats.vert_pos << "\n";
    os << "vert_norm: " << stats.vert_norm << "\n";
    os << "vert_texcoord: " << stats.vert_texcoord << "\n";
    os << "vert_color: " << stats.vert_color << "\n";
    os << "vert_radius: " << stats.vert_radius << "\n";
    os << "vert_tangsp: " << stats.vert_tangsp << "\n";
    os << "texel_ldrs: " << stats.texel_ldrs << "\n";
    os << "texel_hdrs: " << stats.texel_hdrs << "\n";
    os << "memory_ldrs: " << stats.memory_ldrs << "\n";
    os << "memory_hdrs: " << stats.memory_hdrs << "\n";
    os << "memory_elems: " << stats.memory_elems << "\n";
    os << "memory_verts: " << stats.memory_verts << "\n";
    os << "bbox_scn: " << stats.bbox_scn << "\n";
    os << "bbox_nolights: " << stats.bbox_nolights << "\n";

    os << "bbox_min   : " << stats.bbox_scn.min << "\n";
    os << "bbox_max   : " << stats.bbox_scn.max << "\n";
    os << "bbox_size  : " << bbox_size(stats.bbox_scn) << "\n";
    os << "bbox_center: " << bbox_center(stats.bbox_scn) << "\n";

    return os;
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
        mat->rs = (omat->ns >= 1e6f) ? 0 : pow(2 / (omat->ns + 2), 1 / 4.0f);
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
    auto omap = std::unordered_map<std::string, shape*>{{"", nullptr}};
    for (auto omsh : obj->objects) {
        auto shp = new shape();
        shp->name = omsh->name;
        shp->frame = omsh->frame;
        if (omsh->verts.empty()) continue;
        if (omsh->elems.empty()) continue;

        for (auto ogrp : omsh->groups)
            shp->groups.push_back(
                {ogrp.name, mmap[ogrp.matname], ogrp.faceted});
        if (omsh->props.find("subdivision") != omsh->props.end()) {
            shp->subdivision =
                atoi(omsh->props.at("subdivision").at(0).c_str());
        }
        if (omsh->props.find("catmullclark") != omsh->props.end()) {
            shp->catmullclark =
                (bool)atoi(omsh->props.at("catmullclark").at(1).c_str());
        }

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
                    for (auto i = elem.start; i < elem.start + elem.size; i++) {
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
                switch (elem.type) {
                    case obj_element_type::point: {
                        for (auto i = elem.start; i < elem.start + elem.size;
                             i++) {
                            shp->points.push_back(vert_ids[i]);
                            shp->group_ids.push_back(elem.groupid);
                        }
                    } break;
                    case obj_element_type::line: {
                        for (auto i = elem.start;
                             i < elem.start + elem.size - 1; i++) {
                            shp->lines.push_back(
                                {vert_ids[i], vert_ids[i + 1]});
                            shp->group_ids.push_back(elem.groupid);
                        }
                    } break;
                    case obj_element_type::face: {
                        if (as_quads && elem.size == 4) {
                            shp->quads.push_back({vert_ids[elem.start + 0],
                                vert_ids[elem.start + 1],
                                vert_ids[elem.start + 2],
                                vert_ids[elem.start + 3]});
                            shp->group_ids.push_back(elem.groupid);
                        } else if (as_quads && elem.size != 4) {
                            for (auto i = elem.start + 2;
                                 i < elem.start + elem.size; i++) {
                                shp->quads.push_back({vert_ids[elem.start],
                                    vert_ids[i - 1], vert_ids[i], vert_ids[i]});
                                shp->group_ids.push_back(elem.groupid);
                            }
                        } else {
                            for (auto i = elem.start + 2;
                                 i < elem.start + elem.size; i++) {
                                shp->triangles.push_back({vert_ids[elem.start],
                                    vert_ids[i - 1], vert_ids[i]});
                                shp->group_ids.push_back(elem.groupid);
                            }
                        }
                    } break;
                    case obj_element_type::bezier: {
                        if ((elem.size - 1) % 3)
                            throw std::runtime_error("bad obj bezier");
                        for (auto i = elem.start + 1;
                             i < elem.start + elem.size; i += 3) {
                            shp->beziers.push_back({vert_ids[i - 1],
                                vert_ids[i], vert_ids[i + 1], vert_ids[i + 2]});
                            shp->group_ids.push_back(elem.groupid);
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
                if (elem.size == 4) {
                    if (!pos_ids.empty()) {
                        shp->quads_pos.push_back({pos_ids[elem.start + 0],
                            pos_ids[elem.start + 1], pos_ids[elem.start + 2],
                            pos_ids[elem.start + 3]});
                        shp->group_ids.push_back(elem.groupid);
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
                            norm_ids[elem.start + 1], norm_ids[elem.start + 2],
                            norm_ids[elem.start + 3]});
                    }
                } else {
                    if (!pos_ids.empty()) {
                        for (auto i = elem.start + 2;
                             i < elem.start + elem.size; i++) {
                            shp->quads_pos.push_back({pos_ids[elem.start],
                                pos_ids[i - 1], pos_ids[i], pos_ids[i]});
                            shp->group_ids.push_back(elem.groupid);
                        }
                    }
                    if (!texcoord_ids.empty()) {
                        for (auto i = elem.start + 2;
                             i < elem.start + elem.size; i++) {
                            shp->quads_texcoord.push_back(
                                {texcoord_ids[elem.start], texcoord_ids[i - 1],
                                    texcoord_ids[i], texcoord_ids[i]});
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
        omap[omsh->name] = shp;
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
    for (auto shp : scn->shapes)
        for (auto grp : shp->groups) env_mat.erase(grp.mat);
    for (auto mat : env_mat) {
        auto end =
            std::remove(scn->materials.begin(), scn->materials.end(), mat);
        scn->materials.erase(end, scn->materials.end());
    }

    // adjust positions and normals
    for (auto shp : scn->shapes) {
        auto inv_frame = inverse(shp->frame);
        for (auto& p : shp->pos) p = transform_point(inv_frame, p);
        for (auto& n : shp->norm) n = transform_direction(inv_frame, n);
    }

    // convert nodes
    if (!obj->nodes.empty()) {
        for (auto onde : obj->nodes) {
            auto nde = new node();
            nde->name = onde->name;
            nde->cam = cmap.at(onde->camname);
            nde->shp = omap.at(onde->objname);
            nde->env = emap.at(onde->envname);
            nde->translation = onde->translation;
            nde->rotation = onde->rotation;
            nde->scaling = onde->scaling;
            nde->local = onde->local;
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
obj_scene* scene_to_obj(const scene* scn) {
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
        omat->ke_txt = make_texture_info(mat->ke_txt);
        switch (mat->type) {
            case material_type::specular_roughness: {
                omat->kd = {mat->kd.x, mat->kd.y, mat->kd.z};
                omat->ks = {mat->ks.x, mat->ks.y, mat->ks.z};
                omat->kr = {mat->kr.x, mat->kr.y, mat->kr.z};
                omat->kt = {mat->kt.x, mat->kt.y, mat->kt.z};
                omat->ns = (mat->rs) ? 2 / pow(mat->rs, 4.0f) - 2 : 1e6;
                omat->op = mat->op;
                omat->kd_txt = make_texture_info(mat->kd_txt);
                omat->ks_txt = make_texture_info(mat->ks_txt);
                omat->kr_txt = make_texture_info(mat->kr_txt);
                omat->kt_txt = make_texture_info(mat->kt_txt);
                omat->op_txt = make_texture_info(mat->op_txt);
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
                    omat->kd_txt = make_texture_info(mat->kd_txt);
                } else {
                    omat->ks_txt = make_texture_info(mat->ks_txt);
                }
            } break;
            case material_type::specular_glossiness: {
                omat->kd = {mat->kd.x, mat->kd.y, mat->kd.z};
                omat->ks = {mat->ks.x, mat->ks.y, mat->ks.z};
                omat->ns = (mat->rs) ? 2 / pow(1 - mat->rs, 4.0f) - 2 : 1e6;
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
        if (!shp->group_ids.empty()) elem.groupid = shp->group_ids.at(eid);
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

    // convert shapes
    for (auto shp : scn->shapes) {
        auto oobj = new obj_object();
        oobj->name = shp->name;
        oobj->frame = shp->frame;
        if (shp->subdivision)
            oobj->props["subdivision"].push_back(
                std::to_string(shp->subdivision));
        if (shp->catmullclark) oobj->props["catmullclark"].push_back("1");
        auto offset = obj_vertex{(int)obj->pos.size(),
            (int)obj->texcoord.size(), (int)obj->norm.size(),
            (int)obj->color.size(), (int)obj->radius.size()};
        if (shp->frame == identity_frame3f) {
            for (auto& v : shp->pos) obj->pos.push_back(v);
            for (auto& v : shp->norm) obj->norm.push_back(v);
        } else {
            for (auto& v : shp->pos)
                obj->pos.push_back(transform_point(shp->frame, v));
            for (auto& v : shp->norm)
                obj->norm.push_back(transform_direction(shp->frame, v));
        }
        for (auto& v : shp->texcoord) obj->texcoord.push_back(v);
        for (auto& v : shp->color) obj->color.push_back({v.x, v.y, v.z, v.w});
        for (auto& v : shp->radius) obj->radius.push_back(v);
        for (auto grp : shp->groups)
            oobj->groups.push_back(
                {grp.name, (grp.mat) ? grp.mat->name : ""s, grp.faceted});
        for (auto eid = 0; eid < shp->points.size(); eid++) {
            auto vid = shp->points[eid];
            add_elem(shp, oobj, obj_element_type::point, 1, eid);
            add_vert(shp, oobj, offset, vid);
        }
        for (auto eid = 0; eid < shp->lines.size(); eid++) {
            auto line = shp->lines[eid];
            add_elem(shp, oobj, obj_element_type::line, 2, eid);
            for (auto vid : line) add_vert(shp, oobj, offset, vid);
        }
        for (auto eid = 0; eid < shp->triangles.size(); eid++) {
            auto triangle = shp->triangles[eid];
            add_elem(shp, oobj, obj_element_type::face, 3, eid);
            for (auto vid : triangle) add_vert(shp, oobj, offset, vid);
        }
        for (auto eid = 0; eid < shp->quads.size(); eid++) {
            auto quad = shp->quads[eid];
            add_elem(shp, oobj, obj_element_type::face,
                (uint16_t)((quad.z == quad.w) ? 3 : 4), eid);
            if (oobj->elems.back().size == 3) {
                for (auto vid : {quad.x, quad.y, quad.z})
                    add_vert(shp, oobj, offset, vid);
            } else {
                for (auto vid : quad) add_vert(shp, oobj, offset, vid);
            }
        }
        for (auto eid = 0; eid < shp->beziers.size(); eid++) {
            auto bezier = shp->beziers[eid];
            add_elem(shp, oobj, obj_element_type::bezier, 4, eid);
            for (auto vid : bezier) add_vert(shp, oobj, offset, vid);
        }
        for (auto eid = 0; eid < shp->quads_pos.size(); eid++) {
            add_elem(shp, oobj, obj_element_type::face, 4, eid);
            auto last_vid = -1;
            for (auto i = 0; i < 4; i++) {
                if (last_vid == shp->quads_pos[eid][i]) continue;
                auto vert = obj_vertex{-1, -1, -1, -1, -1};
                if (!shp->pos.empty() && !shp->quads_pos.empty())
                    vert.pos = offset.pos + shp->quads_pos[eid][i];
                if (!shp->texcoord.empty() && !shp->quads_texcoord.empty())
                    vert.texcoord =
                        offset.texcoord + shp->quads_texcoord[eid][i];
                if (!shp->norm.empty() && !shp->quads_norm.empty())
                    vert.norm = offset.norm + shp->quads_norm[eid][i];
                oobj->verts.push_back(vert);
                last_vid = shp->quads_pos[eid][i];
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
        auto omat = new obj_material();
        omat->name = env->name + "_mat";
        omat->ke = env->ke;
        omat->ke_txt = make_texture_info(env->ke_txt);
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
            if (nde->shp) onde->objname = nde->shp->name;
            if (nde->env) onde->envname = nde->env->name;
            onde->local = nde->local;
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
    const std::string& filename, const scene* scn, const save_options& opts) {
    auto oscn = scene_to_obj(scn);
    save_obj(filename, oscn, opts.save_textures, opts.skip_missing,
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
            mat->kd = {gmr->baseColorFactor[0], gmr->baseColorFactor[1],
                gmr->baseColorFactor[2]};
            mat->op = gmr->baseColorFactor[3];
            mat->ks = {
                gmr->metallicFactor, gmr->metallicFactor, gmr->metallicFactor};
            mat->rs = gmr->roughnessFactor;
            mat->kd_txt = make_texture_info(gmr->baseColorTexture);
            mat->ks_txt = make_texture_info(gmr->metallicRoughnessTexture);
        }
        if (gmat->pbrSpecularGlossiness) {
            mat->type = material_type::specular_glossiness;
            auto gsg = gmat->pbrSpecularGlossiness;
            mat->kd = {gsg->diffuseFactor[0], gsg->diffuseFactor[1],
                gsg->diffuseFactor[2]};
            mat->op = gsg->diffuseFactor[3];
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
    for (auto gmesh : gltf->meshes) {
        if (gmesh->primitives.empty()) {
            auto shp = new shape();
            shp->name = gmesh->name;
            scn->shapes.push_back(shp);
            continue;
        }
        // primitives
        auto shps = std::vector<shape*>();
        for (auto gprim : gmesh->primitives) {
            auto shp = new shape();
            if (gprim->material) {
                shp->groups.push_back(
                    {"", scn->materials[(int)gprim->material], false});
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
            shps.push_back(shp);
        }
        if (shps.size() == 1) {
            shps[0]->name = gmesh->name;
            scn->shapes.push_back(shps[0]);
        } else {
            auto shp = merge_shapes(shps, true);
            shp->name = gmesh->name;
            scn->shapes.push_back(shp);
            for (auto shp : shps) delete shp;
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
        if (gnde->mesh) nde->shp = scn->shapes[(int)gnde->mesh];
        //            for (auto shp : node->msh->shapes) {
        //                if (node->morph_weights.size() <
        //                shp->morph_targets.size()) {
        //                    node->morph_weights.resize(shp->morph_targets.size());
        //                }
        //            }
        nde->translation = gnde->translation;
        nde->rotation = gnde->rotation;
        nde->scaling = gnde->scale;
        nde->local = mat_to_frame(gnde->matrix);
        //        nde->weights = gnde->weights;
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
        };

    // convert animations
    for (auto ganm : gltf->animations) {
        auto agr = new animation_group();
        agr->name = ganm->name;
        auto sampler_map = std::unordered_map<vec2i, int>();
        for (auto gchannel : ganm->channels) {
            if (sampler_map.find({(int)gchannel->sampler,
                    (int)gchannel->target->path}) == sampler_map.end()) {
                auto gsampler = ganm->get(gchannel->sampler);
                auto anm_ = animation();
                auto anm = &anm_;
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
                sampler_map[{(int)gchannel->sampler,
                    (int)gchannel->target->path}] = (int)agr->animations.size();
                agr->animations.push_back(anm_);
            }
            agr->animations[sampler_map.at({(int)gchannel->sampler,
                                (int)gchannel->target->path})]
                .targets.push_back(scn->nodes[(int)gchannel->target->node]);
        }
        scn->animations.push_back(agr);
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
                gsg->diffuseFactor = {
                    mat->kd[0], mat->kd[1], mat->kd[2], mat->op};
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
                gsg->diffuseFactor = {
                    mat->kd[0], mat->kd[1], mat->kd[2], mat->op};
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
    for (auto shp : scn->shapes) {
        auto gmesh = new glTFMesh();
        gmesh->name = shp->name;
        auto gbuffer = add_opt_buffer(shp->path);
        auto shapes = (is_shape_simple(shp, true)) ? std::vector<shape*>{shp} :
                                                     split_shape(shp, true);
        for (auto shp : shapes) {
            auto gprim = new glTFMeshPrimitive();
            if (!shp->groups.empty())
                gprim->material = glTFid<glTFMaterial>(
                    index(scn->materials, shp->groups.at(0).mat));
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
        for (auto s : shapes)
            if (s != shp) delete s;
        gltf->meshes.push_back(gmesh);
    }

    // hierarchy
    if (scn->nodes.empty()) {
        // shapes
        for (auto shp : scn->shapes) {
            auto gnode = new glTFNode();
            gnode->name = shp->name;
            gnode->mesh = glTFid<glTFMesh>(index(scn->shapes, shp));
            gnode->matrix = frame_to_mat(shp->frame);
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
            if (nde->shp) {
                gnode->mesh = glTFid<glTFMesh>(index(scn->shapes, nde->shp));
            }
            gnode->matrix = frame_to_mat(nde->local);
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
        };

    // animation
    for (auto agr : scn->animations) {
        auto ganm = new glTFAnimation();
        ganm->name = agr->name;
        auto gbuffer = add_opt_buffer(agr->path);
        auto count = 0;
        for (auto& anm_ : agr->animations) {
            auto anm = &anm_;
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

#if YGL_SVG

// Converts an svg scene
scene* svg_to_scene(const svg_scene* sscn, const load_options& opts) {
    log_error("do SVG");
#if 0
    auto scn = new scene();
    auto sid = 0;
    for (auto sshp : sscn->shapes) {
        // TODO: recursive shape
        auto shp = new shape();
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
#else
    return nullptr;
#endif
}

// Load an svg scene
scene* load_svg_scene(const std::string& filename, const load_options& opts) {
    auto sscn = load_svg(filename);
    return svg_to_scene(sscn, opts);
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
    return lerp(ks, fks, 1 - rs);
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
    const shape* shp = nullptr;                      // shape
    const environment* env = nullptr;                // environment
    trace_point_type type = trace_point_type::none;  // type
    vec3f pos = zero3f;                              // pos
    vec3f norm = {0, 0, 1};                          // norm
    vec2f texcoord = zero2f;                         // texcoord
    vec3f ke = zero3f;                               // emission
    vec3f kd = {0, 0, 0};                            // diffuse
    vec3f ks = {0, 0, 0};                            // specular
    vec3f ksg = {1, 1, 1};                           // specular at grazing
    float rs = 0;                                    // specular roughness
    vec3f kr = {0, 0, 0};                            // clear coat
    vec3f krg = {1, 1, 1};                           // clear coat at grazing
    vec3f ktr = {0, 0, 0};                           // thin glass transmission
    vec3f kto = {0, 0, 0};                           // opacity transmission
    bool double_sided = false;                       // double sided
};

// Evaluates the BRDF albedo at normal incidence
vec3f eval_brdf_albedo(const trace_point& pt) {
    return pt.kd + pt.ks + pt.kr + pt.ktr + pt.kto;
}

// Evaluates the weights of each BRDF lobe
std::array<float, 5> eval_brdf_weights(const trace_point& pt) {
    auto w = std::array<float, 5>{{max_element_value(pt.kd),
        max_element_value(pt.ks), max_element_value(pt.kr),
        max_element_value(pt.ktr), max_element_value(pt.kto)}};
    auto s = w[0] + w[1] + w[2] + w[3] + w[4];
    if (!s) return w;
    return {{w[0] / s, w[1] / s, w[2] / s, w[3] / s, w[4] / s}};
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
    return abs(dot(wi, normalize(wn * 2.0f * dot(wo, wn) - wo))) < 0.001f;
}

// Evaluates the BRDF scaled by the cosine of the incoming direction.
// - ggx from [Heitz 2014] and [Walter 2007] and [Lagarde 2014]
// "Understanding the Masking-Shadowing Function in Microfacet-Based
// BRDFs" http://jcgt.org/published/0003/02/03/
// - "Microfacet Models for Refraction through Rough Surfaces" EGSR 07
// https://www.cs.cornell.edu/~srm/publications/EGSR07-btdf.pdf
vec3f eval_surface_brdfcos(const trace_point& pt, const vec3f& wo,
    const vec3f& wi, bool delta = false) {
    auto brdfcos = zero3f;
    auto wh = normalize(wo + wi);
    auto wn = pt.norm;
    if (pt.double_sided && dot(wo, wn) < 0) wn = -wn;

    auto ndo = dot(wn, wo), ndi = dot(wn, wi),
         ndh = clamp(dot(wh, wn), -1.0f, 1.0f);

    if (ndi > 0 && ndo > 0) {
        brdfcos += pt.kd * ndi / pif;
        if (ndh > 0 && pt.rs) {
            auto dg = eval_ggx(pt.rs, ndh, ndi, ndo);
            auto odh = clamp(dot(wo, wh), 0.0f, 1.0f);
            auto ks = fresnel_schlick(pt.ks, odh, pt.rs);
            brdfcos += ks * ndi * dg / (4 * ndi * ndo);
        }
        if (!pt.rs && delta && check_near_mirror(wn, wo, wi)) {
            auto ks = fresnel_schlick(pt.ks, ndo, pt.rs);
            brdfcos += ks;
        }
        if (delta && check_near_mirror(wn, wo, wi)) brdfcos += pt.kr;
    }
    if (wo == -wi && delta) brdfcos += pt.ktr + pt.kto;

    assert(isfinite(brdfcos.x) && isfinite(brdfcos.y) && isfinite(brdfcos.z));
    return brdfcos;
}

// Evaluates the BRDF scaled by the cosine of the incoming direction.
// - uses Kajiya-Kay for hair
vec3f eval_curve_brdfcos(const trace_point& pt, const vec3f& wo,
    const vec3f& wi, bool delta = false) {
    auto brdfcos = zero3f;
    auto wh = normalize(wo + wi);
    auto wn = pt.norm;

    auto ndo = dot(wn, wo), ndi = dot(wn, wi),
         ndh = clamp(dot(wn, wh), 0.0f, 1.0f);
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
    if (wo == -wi && delta) brdfcos += pt.kto;

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
    if (wo == -wi && delta) brdfcos += pt.kto;

    assert(isfinite(brdfcos.x) && isfinite(brdfcos.y) && isfinite(brdfcos.z));
    return brdfcos;
}

// Evaluates the BRDF scaled by the cosine of the incoming direction.
vec3f eval_brdfcos(const trace_point& pt, const vec3f& wo, const vec3f& wi,
    bool delta = false) {
    if (eval_brdf_albedo(pt) == zero3f) return zero3f;
    switch (pt.type) {
        case trace_point_type::none: return zero3f;
        case trace_point_type::surface:
            return eval_surface_brdfcos(pt, wo, wi, delta);
        case trace_point_type::curve:
            return eval_curve_brdfcos(pt, wo, wi, delta);
        case trace_point_type::point:
            return eval_point_brdfcos(pt, wo, wi, delta);
        case trace_point_type::environment: return zero3f;
    }
}

// Compute the weight for sampling the BRDF
float weight_surface_brdfcos(const trace_point& pt, const vec3f& wo,
    const vec3f& wi, bool delta = false) {
    auto weights = eval_brdf_weights(pt);
    auto wh = normalize(wi + wo);
    auto wn = pt.norm;
    if (pt.double_sided && dot(wo, wn) < 0) wn = -wn;

    auto ndo = dot(wn, wo), ndi = dot(wn, wi), ndh = dot(wn, wh);

    auto pdf = 0.0f;
    if (ndo > 0 && ndi > 0) {
        pdf += weights[0] * ndi / pif;
        if (ndh > 0 && pt.rs) {
            auto d = sample_ggx_pdf(pt.rs, ndh);
            auto hdo = dot(wo, wh);
            pdf += weights[1] * d / (4 * hdo);
        }
        if (!pt.rs && delta && check_near_mirror(wn, wo, wi)) pdf += weights[1];
        if (delta && check_near_mirror(wn, wo, wi)) pdf += weights[2];
    }
    if (wi == -wo && delta) pdf += weights[3] + weights[4];

    assert(isfinite(pdf));
    if (!pdf) return 0;
    return 1 / pdf;
}

// Compute the weight for sampling the BRDF
float weight_curve_brdfcos(const trace_point& pt, const vec3f& wo,
    const vec3f& wi, bool delta = false) {
    auto weights = eval_brdf_weights(pt);

    auto pdf = 0.0f;
    pdf += (weights[0] + weights[1] + weights[2]) / (4 * pif);
    if (wi == -wo && delta) pdf += weights[3] + weights[4];

    assert(isfinite(pdf));
    if (!pdf) return 0;
    return 1 / pdf;
}

// Compute the weight for sampling the BRDF
float weight_point_brdfcos(const trace_point& pt, const vec3f& wo,
    const vec3f& wi, bool delta = false) {
    auto weights = eval_brdf_weights(pt);

    auto pdf = 0.0f;
    pdf += (weights[0] + weights[1] + weights[2]) / (4 * pif);
    if (wi == -wo && delta) pdf += weights[3] + weights[4];

    assert(isfinite(pdf));
    if (!pdf) return 0;
    return 1 / pdf;
}

// Compute the weight for sampling the BRDF
float weight_brdfcos(const trace_point& pt, const vec3f& wo, const vec3f& wi,
    bool delta = false) {
    if (eval_brdf_albedo(pt) == zero3f) return 0;
    switch (pt.type) {
        case trace_point_type::none: return 0;
        case trace_point_type::surface:
            return weight_surface_brdfcos(pt, wo, wi, delta);
        case trace_point_type::curve:
            return weight_curve_brdfcos(pt, wo, wi, delta);
        case trace_point_type::point:
            return weight_point_brdfcos(pt, wo, wi, delta);
        case trace_point_type::environment: return 0;
    }
}

// Picks a direction based on the BRDF
std::tuple<vec3f, bool> sample_surface_brdfcos(
    const trace_point& pt, const vec3f& wo, float rnl, const vec2f& rn) {
    auto weights = eval_brdf_weights(pt);
    auto lid = sample_index(weights, rnl);
    auto wn = pt.norm;
    if (pt.double_sided && dot(wo, wn) < 0) wn = -wn;

    auto ndo = dot(wn, wo);
    if (ndo <= 0) return {zero3f, false};

    // sample according to diffuse
    if (lid == 0) {
        auto fp = make_frame_fromz(pt.pos, wn);
        auto rz = sqrtf(rn.y), rr = sqrtf(1 - rz * rz), rphi = 2 * pif * rn.x;
        auto wi_local = vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
        return {transform_direction(fp, wi_local), false};
    }
    // sample according to specular GGX
    else if (lid == 1 && pt.rs) {
        auto fp = make_frame_fromz(pt.pos, wn);
        auto wh_local = sample_ggx(pt.rs, rn);
        auto wh = transform_direction(fp, wh_local);
        return {normalize(wh * 2.0f * dot(wo, wh) - wo), false};
    }
    // sample according to specular mirror
    else if (lid == 1 && !pt.rs) {
        return {normalize(wn * 2.0f * dot(wo, wn) - wo), true};
    }
    // sample according to specular mirror
    else if (lid == 2) {
        return {normalize(wn * 2.0f * dot(wo, wn) - wo), true};
    }
    // transmission hack
    else if (lid == 3 || lid == 4) {
        return {-wo, true};
    } else
        assert(false);

    return {zero3f, false};
}

// Picks a direction based on the BRDF
std::tuple<vec3f, bool> sample_curve_brdfcos(
    const trace_point& pt, const vec3f& wo, float rnl, const vec2f& rn) {
    auto weights = eval_brdf_weights(pt);
    auto lid = sample_index(weights, rnl);
    auto wn = pt.norm;

    // diffuse and specular: samnple a uniform spherical direction
    if (lid == 0 || lid == 1 || lid == 2) {
        auto fp = make_frame_fromz(pt.pos, wn);
        auto rz = 2 * rn.y - 1, rr = sqrtf(1 - rz * rz), rphi = 2 * pif * rn.x;
        auto wi_local = vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
        return {transform_direction(fp, wi_local), false};
    }
    // transmission hack
    else if (lid == 3 || lid == 4) {
        return {-wo, true};
    } else
        assert(false);
    return {zero3f, false};
}

// Picks a direction based on the BRDF
std::tuple<vec3f, bool> sample_point_brdfcos(
    const trace_point& pt, const vec3f& wo, float rnl, const vec2f& rn) {
    auto weights = eval_brdf_weights(pt);
    auto lid = sample_index(weights, rnl);
    auto wn = pt.norm;
    // diffuse and specular: samnple a uniform spherical direction
    if (lid == 0 || lid == 1 || lid == 2) {
        auto fp = make_frame_fromz(pt.pos, wn);
        auto rz = 2 * rn.y - 1, rr = sqrtf(1 - rz * rz), rphi = 2 * pif * rn.x;
        auto wi_local = vec3f{rr * cosf(rphi), rr * sinf(rphi), rz};
        return {transform_direction(fp, wi_local), false};
    }
    // transmission hack
    else if (lid == 3 || lid == 4) {
        // continue ray direction
        return {-wo, true};
    } else
        assert(false);
    return {zero3f, false};
}

// Picks a direction based on the BRDF
std::tuple<vec3f, bool> sample_brdfcos(
    const trace_point& pt, const vec3f& wo, float rnl, const vec2f& rn) {
    if (eval_brdf_albedo(pt) == zero3f) return {zero3f, false};
    switch (pt.type) {
        case trace_point_type::none: return {zero3f, false};
        case trace_point_type::surface:
            return sample_surface_brdfcos(pt, wo, rnl, rn);
        case trace_point_type::curve:
            return sample_curve_brdfcos(pt, wo, rnl, rn);
        case trace_point_type::point:
            return sample_point_brdfcos(pt, wo, rnl, rn);
        case trace_point_type::environment: return {zero3f, false};
    }
}

// Create a point for an environment map. Resolves material with textures.
trace_point eval_environment_point(
    const environment* env, const frame3f& frame, const vec2f& uv) {
    auto pt = trace_point();
    pt.env = env;
    pt.type = trace_point_type::environment;

    pt.pos = eval_pos(env, uv);
    pt.norm = eval_norm(env, uv);
    pt.texcoord = eval_texcoord(env, uv);

    pt.pos = transform_point(frame, pt.pos);
    pt.norm = transform_direction(frame, pt.norm);

    pt.ke = env->ke;
    if (env->ke_txt.txt) {
        auto txt = eval_texture(env->ke_txt, pt.texcoord);
        pt.ke *= {txt.x, txt.y, txt.z};
    }
    return pt;
}

// Create a point for a shape. Resolves geometry and material with textures.
trace_point eval_shape_point(const shape* shp, const frame3f& frame, int eid,
    const vec2f& euv, bool force_double_sided) {
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
    pt.shp = shp;
    switch (get_shape_type(pt.shp)) {
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
    pt.pos = eval_pos(pt.shp, eid, euv);
    pt.norm = eval_norm(pt.shp, eid, euv);
    pt.texcoord = eval_texcoord(pt.shp, eid, euv);

    // shortcuts
    auto mat =
        get_material(pt.shp, eid) ? get_material(pt.shp, eid) : def_material;

    // handle normal map
    if (mat->norm_txt.txt) {
        auto tangsp = eval_tangsp(pt.shp, eid, euv);
        auto txt =
            eval_texture(mat->norm_txt, pt.texcoord, false) * 2.0f - vec4f{1};
        auto ntxt = normalize(vec3f{txt.x, -txt.y, txt.z});
        auto frame = make_frame_fromzx(
            {0, 0, 0}, pt.norm, {tangsp.x, tangsp.y, tangsp.z});
        frame.y *= tangsp.w;
        pt.norm = transform_direction(frame, ntxt);
    }

    // move to world coordinates
    pt.pos = transform_point(frame, pt.pos);
    pt.norm = transform_direction(frame, pt.norm);

    // double-sided
    pt.double_sided = mat->double_sided || force_double_sided;

    // initialized material values
    auto kx = vec3f{1, 1, 1};
    auto op = 1.0f;
    if (!pt.shp->color.empty()) {
        auto col = eval_color(pt.shp, eid, euv);
        kx *= {col.x, col.y, col.z};
        op *= col.w;
    }

    // handle occlusion
    if (mat->occ_txt.txt) {
        auto txt = eval_texture(mat->occ_txt, pt.texcoord);
        kx *= {txt.x, txt.y, txt.z};
    }

    // handle opacity
    op *= mat->op;
    if(mat->op_txt.txt) {
        auto txt = eval_texture(mat->op_txt, pt.texcoord);
        op *= (txt.x + txt.y + txt.z) / 3;
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
                op *= txt.w;
            }
            pt.ks = mat->ks * kx;
            pt.rs = mat->rs;
            if (mat->ks_txt.txt) {
                auto txt = eval_texture(mat->ks_txt, pt.texcoord);
                pt.ks *= {txt.x, txt.y, txt.z};
            }
            pt.kr = mat->kr * kx;
            if (mat->kr_txt.txt) {
                auto txt = eval_texture(mat->kr_txt, pt.texcoord);
                pt.kr *= {txt.x, txt.y, txt.z};
            }
            pt.ktr = mat->kt * kx;
            if (mat->kt_txt.txt) {
                auto txt = eval_texture(mat->kt_txt, pt.texcoord);
                pt.ktr *= {txt.x, txt.y, txt.z};
            }
        } break;
        case material_type::metallic_roughness: {
            auto kb = mat->kd * kx;
            if (mat->kd_txt.txt) {
                auto txt = eval_texture(mat->kd_txt, pt.texcoord);
                kb *= {txt.x, txt.y, txt.z};
                op *= txt.w;
            }
            auto km = mat->ks.x;
            pt.rs = mat->rs;
            if (mat->ks_txt.txt) {
                auto txt = eval_texture(mat->ks_txt, pt.texcoord);
                km *= txt.y;
                pt.rs *= txt.z;
            }
            pt.kd = kb * (1 - km);
            pt.ks = kb * km + vec3f{0.04f} * (1 - km);
        } break;
        case material_type::specular_glossiness: {
            pt.kd = mat->kd * kx;
            if (mat->kd_txt.txt) {
                auto txt = eval_texture(mat->kd_txt, pt.texcoord);
                pt.kd *= {txt.x, txt.y, txt.z};
                op *= txt.w;
            }
            pt.ks = mat->ks * kx;
            pt.rs = mat->rs;
            if (mat->ks_txt.txt) {
                auto txt = eval_texture(mat->ks_txt, pt.texcoord);
                pt.ks *= {txt.x, txt.y, txt.z};
                pt.rs *= txt.w;
            }
            pt.rs = 1 - pt.rs;  // glossiness -> roughnes
            pt.kr = mat->kr * kx;
            if (mat->kr_txt.txt) {
                auto txt = eval_texture(mat->kr_txt, pt.texcoord);
                pt.kr *= {txt.x, txt.y, txt.z};
            }
            pt.ktr = mat->kt * kx;
            if (mat->kt_txt.txt) {
                auto txt = eval_texture(mat->kt_txt, pt.texcoord);
                pt.ktr *= {txt.x, txt.y, txt.z};
            }
        } break;
    }

    // set up final values
    pt.ke *= op;
    pt.kd *= op;
    if (pt.ks != zero3f && pt.rs < 0.9999f) {
        pt.ks *= op;
        pt.ksg *= op;
        pt.rs = pt.rs * pt.rs;
        pt.rs = clamp(pt.rs, 0.02f * 0.02f, 1.0f);
    } else {
        pt.ks = zero3f;
        pt.rs = 0;
    }
    pt.kr *= op;
    pt.krg *= op;
    pt.ktr *= op;
    pt.kto = {1.0f - op, 1.0f - op, 1.0f - op};

    // done
    return pt;
}

// Sample weight for a light point.
float weight_light(
    const trace_lights& lights, const trace_point& lpt, const trace_point& pt) {
    switch (lpt.type) {
        case trace_point_type::none: {
            throw std::runtime_error("should not have gotten here");
        } break;
        case trace_point_type::point: {
            auto area = sample_distribution_weightsum(
                lights.shape_distribs.at(lpt.shp));
            auto dist = length(lpt.pos - pt.pos);
            return area / (dist * dist);
        } break;
        case trace_point_type::curve: {
            throw std::runtime_error("not implemented yet");
        } break;
        case trace_point_type::surface: {
            auto area = sample_distribution_weightsum(
                lights.shape_distribs.at(lpt.shp));
            auto dist = length(lpt.pos - pt.pos);
            return area * abs(dot(lpt.norm, normalize(lpt.pos - pt.pos))) /
                   (dist * dist);
        } break;
        case trace_point_type::environment: {
            return 4 * pif;
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
    if (lgt.shp) {
        auto& dst = lights.shape_distribs.at(lgt.shp);
        auto sample = sample_shape(lgt.shp, dst, rel, ruv);
        return eval_shape_point(lgt.shp, lgt.frame, sample.first, sample.second,
            params.double_sided);
    } else if (lgt.env) {
        // BUG: this is not uniform sampling
        return eval_environment_point(lgt.env, lgt.frame, ruv);
    } else {
        throw std::runtime_error("should not have gotten here");
    }
    return {};
}

// Picks a point on a light.
trace_point sample_lights(const trace_lights& lights, const trace_point& pt,
    float rnl, float rne, const vec2f& ruv, const trace_params& params) {
    auto lidx = sample_distribution_discrete(lights.light_distrib, rnl);
    auto& lgt = lights.lights.at(lidx);
    return sample_light(lights, lgt, pt, rne, ruv, params);
}

// Intersects a ray with the scn and return the point (or env point).
trace_point intersect_scene(const scene* scn, const bvh_tree* bvh,
    const ray3f& ray, const trace_params& params) {
    auto iid = 0, sid = 0, eid = 0;
    auto euv = zero2f;
    auto ray_t = 0.0f;
    if (intersect_bvh(bvh, ray, false, ray_t, iid, sid, eid, euv)) {
        if (scn->nodes.empty()) {
            return eval_shape_point(scn->shapes[iid], scn->shapes[iid]->frame,
                eid, euv, params.double_sided);
        } else {
            return eval_shape_point(scn->nodes[iid]->shp,
                scn->nodes[iid]->frame_, eid, euv, params.double_sided);
        }
    } else if (!scn->environments.empty()) {
        return eval_environment_point(scn->environments[0],
            scn->environments[0]->frame,
            eval_uv(
                scn->environments[0], transform_direction_inverse(
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
            if (!cpt.shp) break;
            weight *= cpt.ktr + cpt.kto;
            if (weight == zero3f) break;
        }
        return weight;
    }
}

// Mis weight.
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
    if (eval_brdf_albedo(pt) == zero3f || lights.empty()) return l;

    // trace path
    auto weight = vec3f{1, 1, 1};
    auto emission = false;
    for (auto bounce = 0; bounce < params.max_depth; bounce++) {
        // emission
        if (emission) l += weight * eval_emission(pt, wo);

        // early exit
        if (eval_brdf_albedo(pt) == zero3f) break;

        // direct – light
        auto rll = sample_next1f(pxl, params.rng, params.nsamples);
        auto rle = sample_next1f(pxl, params.rng, params.nsamples);
        auto rluv = sample_next2f(pxl, params.rng, params.nsamples);
        auto& lgt = lights.lights[sample_distribution_discrete(
            lights.light_distrib, rll)];
        auto lpt = sample_light(lights, lgt, pt, rle, rluv, params);
        auto lw = weight_light(lights, lpt, pt) *
                  sample_distribution_weightsum(lights.light_distrib);
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
        auto bpt = intersect_scene(scn, bvh, make_ray(pt.pos, bwi), params);
        auto bw = weight_brdfcos(pt, wo, bwi, bdelta);
        auto bke = eval_emission(bpt, -bwi);
        auto bbc = eval_brdfcos(pt, wo, bwi, bdelta);
        auto bld = bke * bbc * bw;
        if (bld != zero3f) {
            // TODO: possible BUG; check for light distribution weight here
            l += weight * bld * weight_mis(bw, weight_light(lights, bpt, pt));
        }

        // skip recursion if path ends
        if (bounce == params.max_depth - 1) break;

        // continue path
        weight *= eval_brdfcos(pt, wo, bwi, bdelta) *
                  weight_brdfcos(pt, wo, bwi, bdelta);
        if (weight == zero3f) break;

        // roussian roulette
        if (bounce > 2) {
            auto rrprob =
                1.0f - min(max_element_value(eval_brdf_albedo(pt)), 0.95f);
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
    if (eval_brdf_albedo(pt) == zero3f || lights.empty()) return l;

    // trace path
    auto weight = vec3f{1, 1, 1};
    auto emission = false;
    for (auto bounce = 0; bounce < params.max_depth; bounce++) {
        // emission
        if (emission) l += weight * eval_emission(pt, wo);

        // early exit
        if (eval_brdf_albedo(pt) == zero3f) break;

        // direct
        auto rll = sample_next1f(pxl, params.rng, params.nsamples);
        auto rle = sample_next1f(pxl, params.rng, params.nsamples);
        auto rluv = sample_next2f(pxl, params.rng, params.nsamples);
        auto& lgt = lights.lights[sample_distribution_discrete(
            lights.light_distrib, rll)];
        auto lpt = sample_light(lights, lgt, pt, rle, rluv, params);
        auto lwi = normalize(lpt.pos - pt.pos);
        auto ld = eval_emission(lpt, -lwi) * eval_brdfcos(pt, wo, lwi) *
                  weight_light(lights, lpt, pt) *
                  sample_distribution_weightsum(lights.light_distrib);
        if (ld != zero3f) {
            l += weight * ld * eval_transmission(scn, bvh, pt, lpt, params);
        }

        // skip recursion if path ends
        if (bounce == params.max_depth - 1) break;

        // roussian roulette
        if (bounce > 2) {
            auto rrprob =
                1.0f - min(max_element_value(eval_brdf_albedo(pt)), 0.95f);
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

        auto bpt = intersect_scene(scn, bvh, make_ray(pt.pos, bwi), params);
        emission = false;

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
    if (eval_brdf_albedo(pt) == zero3f || lights.empty()) return l;

    // trace path
    auto weight = vec3f{1, 1, 1};
    for (auto bounce = 0; bounce < params.max_depth; bounce++) {
        // early exit
        if (eval_brdf_albedo(pt) == zero3f) break;

        // direct
        auto rll = sample_next1f(pxl, params.rng, params.nsamples);
        auto rle = sample_next1f(pxl, params.rng, params.nsamples);
        auto rluv = sample_next2f(pxl, params.rng, params.nsamples);
        auto& lgt = lights.lights[sample_distribution_discrete(
            lights.light_distrib, rll)];
        auto lpt = sample_light(lights, lgt, pt, rle, rluv, params);
        auto lwi = normalize(lpt.pos - pt.pos);
        auto ld = eval_emission(lpt, -lwi) * eval_brdfcos(pt, wo, -lwi) *
                  weight_light(lights, lpt, pt) *
                  sample_distribution_weightsum(lights.light_distrib);
        if (ld != zero3f) {
            l += weight * ld * eval_transmission(scn, bvh, pt, lpt, params);
        }

        // skip recursion if path ends
        if (bounce == params.max_depth - 1) break;

        // roussian roulette
        if (bounce > 2) {
            auto rrprob =
                1.0f - min(max_element_value(eval_brdf_albedo(pt)), 0.95f);
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

        auto bpt = intersect_scene(scn, bvh, make_ray(pt.pos, bwi), params);

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
    if (eval_brdf_albedo(pt) == zero3f || lights.empty()) return l;

    // ambient
    l += params.ambient * eval_brdf_albedo(pt);

    // direct
    for (auto& lgt : lights.lights) {
        auto rle = sample_next1f(pxl, params.rng, params.nsamples);
        auto rluv = sample_next2f(pxl, params.rng, params.nsamples);
        auto lpt = sample_light(lights, lgt, pt, rle, rluv, params);
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
             trace_direct(scn, bvh, lights, rpt, -wi, bounce + 1, pxl, params);
    }

    // opacity
    if (pt.ktr + pt.kto != zero3f) {
        auto opt = intersect_scene(scn, bvh, make_ray(pt.pos, -wo), params);
        l += (pt.ktr + pt.kto) *
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
    if (eval_brdf_albedo(pt) == zero3f) return l;

    // brdf*light
    l += eval_brdfcos(pt, wo, wo) * pif;

    // opacity
    if (bounce >= params.max_depth) return l;
    if (pt.ktr + pt.kto != zero3f) {
        auto opt = intersect_scene(scn, bvh, make_ray(pt.pos, -wo), params);
        l += (pt.ktr + pt.kto) *
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
    auto wn = pt.norm;
    if (pt.double_sided && dot(wn, wo) < 0) wn = -wn;
    return wn * 0.5f + vec3f{0.5f, 0.5f, 0.5f};
}

// Debug frontfacing.
vec3f trace_debug_frontfacing(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt, const vec3f& wo,
    trace_pixel& pxl, const trace_params& params) {
    auto wn = pt.norm;
    if (pt.double_sided && dot(wn, wo) < 0) wn = -wn;
    return dot(wn, wo) > 0 ? vec3f{0, 1, 0} : vec3f{1, 0, 0};
}

// Debug previewing.
vec3f trace_debug_albedo(const scene* scn, const bvh_tree* bvh,
    const trace_lights& lights, const trace_point& pt, const vec3f& wo,
    trace_pixel& pxl, const trace_params& params) {
    return eval_brdf_albedo(pt);
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
    {trace_shader_type::debug_frontfacing, trace_debug_frontfacing},
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
    auto pt = intersect_scene(scn, bvh, ray, params);
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
    auto pt = intersect_scene(scn, bvh, ray, params);
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
    const trace_params& params, const std::function<void(int, int)>& callback) {
    pixels = make_trace_pixels(img, params);

    // render preview
    if (params.preview_resolution) {
        auto pparams = params;
        pparams.resolution = params.preview_resolution;
        pparams.nsamples = 1;
        pparams.filter = ygl::trace_filter_type::box;
        auto pimg = image4f((int)std::round(cam->aspect * pparams.resolution),
            pparams.resolution);
        auto ppixels = make_trace_pixels(pimg, pparams);
        trace_samples(scn, cam, bvh, lights, pimg, ppixels, 1, pparams);
        resize_image(pimg, img, ygl::resize_filter::box);
    } else {
        for (auto& p : img) p = zero4f;
    }
    if (callback) callback(0, 0);

    // start rendering
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
                    if (!tid && callback) callback(s, j);
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

    for (auto shp : scn->shapes) {
        if (!has_emission(shp)) continue;
        if (!contains(lights.shape_distribs, shp))
            lights.shape_distribs[shp] = make_shape_distribution(shp);
    }

    if (scn->nodes.empty()) {
        for (auto shp : scn->shapes) {
            if (!has_emission(shp)) continue;
            auto lgt = trace_light();
            lgt.shp = shp;
            lgt.frame = shp->frame;
            lights.lights.push_back(lgt);
        }
    } else {
        for (auto nde : scn->nodes) {
            if (!nde->shp) continue;
            auto shp = nde->shp;
            if (!has_emission(shp)) continue;
            auto lgt = trace_light();
            lgt.shp = shp;
            lgt.frame = nde->frame_;
            lights.lights.push_back(lgt);
        }
    }

    for (auto env : scn->environments) {
        if (env->ke == zero3f) continue;
        auto lgt = trace_light();
        lgt.env = env;
        lgt.frame = env->frame;
        lights.lights.push_back(lgt);
    }

    lights.light_distrib =
        make_distribution(lights.lights.size(), [](auto) { return 1; });

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
            pxl.rng = make_rng(params.seed, (j * img.width() + i) * 2 + 1);
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
        } else if (obj_streq(cmd, "of")) {
            obj_parse(ss, oobj->frame);
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
            obj_parse(ss, env->matname);
            obj_parse(ss, env->frame);
            obj->environments.push_back(env);
        } else if (obj_streq(cmd, "n")) {
            auto nde = new obj_node();
            obj_parse(ss, nde->name);
            obj_parse(ss, nde->parent);
            obj_parse(ss, nde->camname);
            obj_parse(ss, nde->objname);
            obj_parse(ss, nde->envname);
            obj_parse(ss, nde->local);
            obj_parse(ss, nde->translation);
            obj_parse(ss, nde->rotation);
            obj_parse(ss, nde->scaling);
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
// Dumps a value
inline void obj_dump(char*& s, const char* val) {
    while (*val) *s++ = *val++;
}
// Dumps a value
inline void obj_dump(char*& s, const std::string& val) {
    auto val_ = val.c_str();
    while (*val_) *s++ = *val_++;
}

// Dumps a value
inline void obj_dump(char*& s, int val) { s += sprintf(s, "%d", val); }
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
        obj_dump_sp(fs, env->matname);
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
        obj_dump_sp(fs, nde->local);
        obj_dump_sp(fs, nde->translation);
        obj_dump_sp(fs, nde->rotation);
        obj_dump_sp(fs, nde->scaling);
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
        if (oobj->frame != identity_frame3f)
            obj_dump_line(fs, "of", oobj->frame);
        for (auto& kv : oobj->props) {
            obj_dump_line(fs, "op", join({kv.first}, kv.second));
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
void load_images(
    const glTF* gltf, const std::string& dirname, bool skip_missing) {
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

    auto make_quad = [](std::string name, material* mat, vec3f pos,
                         vec3f rot = {0, 0, 0}, vec2f size = {2, 2}) {
        auto shp = new shape();
        shp->name = name;
        shp->frame = translation_frame(pos);
        if (rot != zero3f) {
            shp->frame = shp->frame *
                         rotation_frame(vec3f{0, 0, 1}, rot[2] * pif / 180) *
                         rotation_frame(vec3f{0, 1, 0}, rot[1] * pif / 180) *
                         rotation_frame(vec3f{1, 0, 0}, rot[0] * pif / 180);
        }
        shp->groups.push_back({"", mat, false});
        ygl::make_quad(shp->quads, shp->pos, shp->norm, shp->texcoord, {1, 1},
            size, {1, 1});
        return shp;
    };

    auto make_box = [](std::string name, material* mat, vec3f pos,
                        vec3f rot = {0, 0, 0}, vec3f size = {2, 2, 2}) {
        auto shp = new shape();
        shp->name = name;
        shp->frame = translation_frame(pos);
        if (rot != zero3f) {
            shp->frame = shp->frame *
                         rotation_frame(vec3f{0, 0, 1}, rot[2] * pif / 180) *
                         rotation_frame(vec3f{0, 1, 0}, rot[1] * pif / 180) *
                         rotation_frame(vec3f{1, 0, 0}, rot[0] * pif / 180);
        }
        shp->groups.push_back({"", mat, false});
        shp->name = name;
        make_cube(shp->quads, shp->pos, shp->norm, shp->texcoord, {1, 1, 1},
            size, {1, 1, 1});
        return shp;
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
    scn->shapes.push_back(
        make_quad("cb_floor", scn->materials[0], {0, 0, 0}, {-90, 0, 0}));
    scn->shapes.push_back(
        make_quad("cb_ceiling", scn->materials[0], {0, 2, 0}, {90, 0, 0}));
    scn->shapes.push_back(
        make_quad("cb_back", scn->materials[0], {0, 1, -1}, {0, 0, 0}));
    scn->shapes.push_back(
        make_quad("cb_left", scn->materials[2], {+1, 1, 0}, {0, -90, 0}));
    scn->shapes.push_back(
        make_quad("cb_right", scn->materials[1], {-1, 1, 0}, {0, 90, 0}));
    scn->shapes.push_back(make_box("cb_tallbox", scn->materials[0],
        {-0.33f, 0.6f, -0.29f}, {0, 15, 0}, {0.6f, 1.2f, 0.6f}));
    scn->shapes.push_back(make_box("cb_shortbox", scn->materials[0],
        {0.33f, 0.3f, 0.33f}, {0, -15, 0}, {0.6f, 0.6f, 0.6f}));
    scn->shapes.push_back(make_quad("cb_light", scn->materials[3],
        {0, 1.999f, 0}, {90, 0, 0}, {0.5f, 0.5f}));
    return scn;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// EXAMPLE SCENES
// -----------------------------------------------------------------------------
namespace ygl {

// Makes/updates a test texture
void update_proc_elem(scene* scn, texture* txt, const proc_texture* ptxt) {
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
void update_proc_elem(scene* scn, material* mat, const proc_material* pmat) {
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
    mat->ke_txt = {};
    mat->kd_txt = {};
    mat->ks_txt = {};
    mat->kr_txt = {};
    mat->kt_txt = {};

    switch (pmat->type) {
        case proc_material_type::none: break;
        case proc_material_type::emission: {
            mat->ke = pmat->emission * pmat->color;
            mat->ke_txt.txt = txt;
        } break;
        case proc_material_type::matte: {
            mat->kd = pmat->color;
            mat->kd_txt.txt = txt;
        } break;
        case proc_material_type::plastic: {
            mat->kd = pmat->color;
            mat->ks = {0.04f, 0.04f, 0.04f};
            mat->rs = pmat->roughness;
            mat->kd_txt.txt = txt;
        } break;
        case proc_material_type::metal: {
            mat->ks = pmat->color;
            mat->rs = pmat->roughness;
            mat->ks_txt.txt = txt;
        } break;
        case proc_material_type::transparent: {
            mat->kd = pmat->color;
            mat->op = pmat->opacity;
            mat->kd_txt.txt = txt;
        } break;
        default: throw std::runtime_error("should not have gotten here");
    }

    mat->norm_txt.txt = norm;
}

// Makes/updates a test shape
void update_proc_elem(scene* scn, shape* shp, const proc_shape* pshp) {
    if (pshp->name == "") throw std::runtime_error("cannot use empty name");
    shp->name = pshp->name;
    shp->frame = pshp->frame;
    shp->groups.clear();
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

    shp->groups.push_back(
        {shp->name, find_named_elem(scn->materials, pshp->material), false});

    switch (pshp->type) {
        case proc_shape_type::floor: {
            make_quad(shp->quads, shp->pos, shp->norm, shp->texcoord,
                {pshp->tesselation.x, pshp->tesselation.y},
                {pshp->size.x, pshp->size.y}, {pshp->uvsize.x, pshp->uvsize.y});
            for (auto& p : shp->pos) p = {-p.x, p.z, p.y};
            for (auto& n : shp->norm) n = {n.x, n.z, n.y};
        } break;
        case proc_shape_type::quad: {
            make_quad(shp->quads, shp->pos, shp->norm, shp->texcoord,
                {pshp->tesselation.x, pshp->tesselation.y},
                {pshp->size.x, pshp->size.y}, {pshp->uvsize.x, pshp->uvsize.y});
        } break;
        case proc_shape_type::cube: {
            make_cube(shp->quads, shp->pos, shp->norm, shp->texcoord,
                pshp->tesselation, pshp->size, pshp->uvsize);
        } break;
        case proc_shape_type::cube_rounded: {
            make_cube_rounded(shp->quads, shp->pos, shp->norm, shp->texcoord,
                pshp->tesselation, pshp->size, pshp->uvsize, pshp->rounded);
        } break;
        case proc_shape_type::sphere: {
            make_sphere(shp->quads, shp->pos, shp->norm, shp->texcoord,
                {pshp->tesselation.x, pshp->tesselation.y}, pshp->size.x,
                {pshp->uvsize.x, pshp->uvsize.y});
        } break;
        case proc_shape_type::sphere_cube: {
            make_sphere_cube(shp->quads, shp->pos, shp->norm, shp->texcoord,
                pshp->tesselation.x, pshp->size.x, pshp->uvsize.x);
        } break;
        case proc_shape_type::disk: {
            make_disk(shp->quads, shp->pos, shp->norm, shp->texcoord,
                {pshp->tesselation.x, pshp->tesselation.y}, pshp->size.x,
                {pshp->uvsize.x, pshp->uvsize.y});
        } break;
        case proc_shape_type::disk_quad: {
            make_disk_quad(shp->quads, shp->pos, shp->norm, shp->texcoord,
                pshp->tesselation.x, pshp->size.x, pshp->uvsize.x);
        } break;
        case proc_shape_type::disk_bulged: {
            make_disk_bulged(shp->quads, shp->pos, shp->norm, shp->texcoord,
                pshp->tesselation.x, pshp->size.x, pshp->uvsize.x,
                pshp->rounded);
        } break;
        case proc_shape_type::cylinder: {
            make_cylinder(shp->quads, shp->pos, shp->norm, shp->texcoord,
                pshp->tesselation, {pshp->size.x, pshp->size.y}, pshp->uvsize);
        } break;
        case proc_shape_type::cylinder_rounded: {
            make_cylinder_rounded(shp->quads, shp->pos, shp->norm,
                shp->texcoord, pshp->tesselation, {pshp->size.x, pshp->size.y},
                pshp->uvsize, pshp->rounded);
        } break;
        case proc_shape_type::geodesic_sphere: {
            auto level = 0;
            auto ntris = (pshp->tesselation.x * pshp->tesselation.y) / 4;
            while (ntris > 4) {
                ntris /= 4;
                level += 1;
            }
            make_geodesic_sphere(shp->triangles, shp->pos, level);
            shp->norm = shp->pos;
            shp->texcoord.assign(shp->pos.size(), {0, 0});
        } break;
        case proc_shape_type::sphere_flipcap: {
            make_sphere_flipcap(shp->quads, shp->pos, shp->norm, shp->texcoord,
                {pshp->tesselation.x, pshp->tesselation.y}, pshp->size.x,
                {pshp->uvsize.x, pshp->uvsize.y},
                {-pshp->rounded, pshp->rounded});
        } break;
        case proc_shape_type::suzanne: {
            make_suzanne(shp->quads, shp->pos, 0);
        } break;
        case proc_shape_type::cubep: {
            make_cube(shp->quads, shp->pos, 0);
        } break;
        case proc_shape_type::fvcube: {
            make_fvcube(shp->quads_pos, shp->pos, shp->quads_norm, shp->norm,
                shp->quads_texcoord, shp->texcoord, 0);
        } break;
        case proc_shape_type::fvsphere: {
            make_sphere(shp->quads, shp->pos, shp->norm, shp->texcoord,
                {pshp->tesselation.x, pshp->tesselation.y}, pshp->size.x,
                {pshp->uvsize.x, pshp->uvsize.y});
        } break;
        case proc_shape_type::matball: {
            make_sphere_flipcap(shp->quads, shp->pos, shp->norm, shp->texcoord,
                {pshp->tesselation.x, pshp->tesselation.y}, pshp->size.x,
                {pshp->uvsize.x, pshp->uvsize.y},
                {-pshp->rounded, pshp->rounded});
            auto shp1 = new shape();
            shp1->name = "interior";
            shp1->frame = shp->frame;
            make_sphere(shp1->quads, shp1->pos, shp1->norm, shp1->texcoord,
                {pshp->tesselation.x, pshp->tesselation.y}, pshp->size.x * 0.8f,
                {pshp->uvsize.x, pshp->uvsize.y});
            shp1->groups.push_back({"interior",
                find_named_elem(scn->materials, pshp->interior), false});
            merge_into(shp, shp1, false);
            delete shp1;
        } break;
        case proc_shape_type::point: {
            make_point(
                shp->points, shp->pos, shp->norm, shp->texcoord, shp->radius);
        } break;
        case proc_shape_type::pointscube: {
            make_random_points(shp->points, shp->pos, shp->norm, shp->texcoord,
                shp->radius, pshp->tesselation.x, pshp->size, pshp->uvsize.x,
                pshp->radius);
        } break;
        case proc_shape_type::hairball: {
            auto shp1 = new shape();
            shp1->name = "interior";
            shp1->frame = shp->frame;
            auto steps = pow2(5);
            make_sphere_cube(shp1->quads, shp1->pos, shp1->norm, shp1->texcoord,
                steps, pshp->size.x * 0.8f, 1);
            shp1->radius.assign(shp1->pos.size(), 0);
            shp1->groups.push_back({"interior",
                find_named_elem(scn->materials, pshp->interior), false});
            make_hair(shp->lines, shp->pos, shp->norm, shp->texcoord,
                shp->radius, {pshp->tesselation.x, pshp->tesselation.y}, {},
                shp1->quads, shp1->pos, shp1->norm, shp1->texcoord,
                pshp->hair_params);
            merge_into(shp, shp1, false);
            delete shp1;
        } break;
        case proc_shape_type::beziercircle: {
            make_bezier_circle(shp->beziers, shp->pos);
            shp->subdivision = 2;
        } break;
        default: throw std::runtime_error("should not have gotten here");
    }

    if (pshp->flip_yz) {
        for (auto& p : shp->pos) std::swap(p.y, p.z);
        for (auto& n : shp->norm) std::swap(n.y, n.z);
    }

    for (auto i = 0; i < pshp->subdivision; i++) {
        subdivide_shape_once(shp, pshp->catmull_clark);
    }

    if (pshp->faceted) facet_shape(shp);
}

// Makes/updates a test shape
void update_proc_elem(scene* scn, camera* cam, const proc_camera* pcam) {
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
    scene* scn, environment* env, const proc_environment* penv) {
    if (penv->name == "") throw std::runtime_error("cannot use empty name");
    env->name = penv->name;
    env->frame = penv->frame * rotation_frame(vec3f{0, 1, 0}, penv->rotation);
    env->ke = penv->emission * penv->color;
    env->ke_txt.txt = find_named_elem(scn->textures, penv->texture);
}

// Makes/updates a test shape
void update_proc_elem(scene* scn, node* nde, const proc_node* pnde) {
    if (pnde->name == "") throw std::runtime_error("cannot use empty name");
    nde->name = pnde->name;
    nde->local = pnde->frame;
    nde->translation = pnde->translation;
    nde->rotation = pnde->rotation;
    nde->scaling = pnde->scaling;
    nde->cam = find_named_elem(scn->cameras, pnde->camera);
    nde->shp = find_named_elem(scn->shapes, pnde->shape);
    nde->env = find_named_elem(scn->environments, pnde->environment);
}

// Fake function for template programming.
void update_proc_elem(scene* scn, animation* shp, const proc_animation* pshp) {}

// Makes/updates a test animation
void update_proc_elem(
    scene* scn, animation_group* agr, const proc_animation* panm) {
    if (panm->name == "") throw std::runtime_error("cannot use empty name");
    if (agr->animations.size() != 1) {
        agr->animations.clear();
        agr->animations.push_back(animation());
    }
    agr->name = panm->name;
    auto anm = &agr->animations.front();
    anm->name = panm->name;
    anm->type = (!panm->bezier) ? keyframe_type::linear : keyframe_type::bezier;
    anm->times = panm->times;
    for (auto& v : anm->times) v *= panm->speed;
    anm->translation = panm->translation;
    anm->rotation = panm->rotation;
    anm->scaling = panm->scaling;
    for (auto& v : anm->translation) v *= panm->scale;
    for (auto& v : anm->scaling) v *= panm->scale;
    anm->targets.clear();
    for (auto& nde : panm->nodes)
        anm->targets.push_back(find_named_elem(scn->nodes, nde));
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
    telems = ntelems;
}

// remove duplicates
void remove_duplicates(proc_scene* scn) {
    remove_duplicate_elems(scn->cameras);
    remove_duplicate_elems(scn->textures);
    remove_duplicate_elems(scn->materials);
    remove_duplicate_elems(scn->shapes);
    remove_duplicate_elems(scn->environments);
    remove_duplicate_elems(scn->nodes);
    remove_duplicate_elems(scn->animations);
}

std::vector<proc_texture*>& proc_texture_presets() {
    static auto presets = std::vector<proc_texture*>();
    if (!presets.empty()) return presets;

    auto make_texture = [](const std::string& name, proc_texture_type type) {
        auto params = new proc_texture();
        params->name = name;
        params->type = type;
        return params;
    };

    presets.push_back(make_texture("grid", proc_texture_type::grid));
    presets.push_back(make_texture("checker", proc_texture_type::checker));
    presets.push_back(make_texture("colored", proc_texture_type::colored));
    presets.push_back(make_texture("rcolored", proc_texture_type::rcolored));
    presets.push_back(make_texture("bump", proc_texture_type::bump));
    presets.back()->tile_size = 32;
    presets.push_back(make_texture("tgrid", proc_texture_type::bump));
    presets.back()->tile_size = 32;
    presets.push_back(make_texture("uv", proc_texture_type::uv));
    presets.push_back(make_texture("gamma", proc_texture_type::gamma));
    presets.push_back(make_texture("gridn", proc_texture_type::grid));
    presets.back()->bump_to_normal = true;
    presets.back()->bump_to_normal = true;
    presets.back()->bump_scale = 4;
    presets.push_back(make_texture("tgridn", proc_texture_type::grid));
    presets.back()->tile_size = 32;
    presets.back()->bump_to_normal = true;
    presets.back()->bump_to_normal = true;
    presets.back()->bump_scale = 4;
    presets.push_back(make_texture("bumpn", proc_texture_type::bump));
    presets.back()->tile_size = 32;
    presets.back()->bump_to_normal = true;
    presets.back()->bump_scale = 4;
    presets.push_back(make_texture("noise", proc_texture_type::noise));
    presets.push_back(make_texture("ridge", proc_texture_type::ridge));
    presets.push_back(make_texture("fbm", proc_texture_type::fbm));
    presets.push_back(
        make_texture("turbulence", proc_texture_type::turbulence));

    presets.push_back(make_texture("gammaf", proc_texture_type::gammaf));
    presets.push_back(make_texture("sky1", proc_texture_type::sky));
    presets.back()->sky_sunangle = pif / 4;
    presets.push_back(make_texture("sky2", proc_texture_type::sky));
    presets.back()->sky_sunangle = pif / 2;

    return presets;
}

std::vector<proc_material*>& proc_material_presets() {
    static auto presets = std::vector<proc_material*>();
    if (!presets.empty()) return presets;

    auto make_material = [](const std::string& name, proc_material_type type,
                             const vec3f& color, float roughness = 1) {
        auto params = new proc_material();
        params->name = name;
        params->type = type;
        params->color = color;
        params->roughness = roughness;
        return params;
    };
    auto make_materialt = [](const std::string& name, proc_material_type type,
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
    auto mirror = 0.0f;

    auto params = std::vector<proc_material>();

    presets.push_back(make_materialt("matte_floor", matte, "grid"));

    presets.push_back(make_material("matte_gray", matte, gray));
    presets.push_back(make_material("matte_red", matte, red));
    presets.push_back(make_material("matte_green", matte, green));
    presets.push_back(make_material("matte_blue", matte, blue));
    presets.push_back(make_materialt("matte_grid", matte, "grid"));
    presets.push_back(make_materialt("matte_colored", matte, "colored"));
    presets.push_back(make_materialt("matte_uv", matte, "uv"));

    presets.push_back(make_material("plastic_red", plastic, red, rough));
    presets.push_back(make_material("plastic_green", plastic, green, rough));
    presets.push_back(make_material("plastic_blue", plastic, blue, sharp));
    presets.push_back(
        make_materialt("plastic_colored", plastic, "colored", rough));
    presets.push_back(
        make_material("plastic_blue_bumped", plastic, blue, sharp));
    presets.back()->normal = "bumpn";
    presets.push_back(
        make_materialt("plastic_colored_bumped", plastic, "colored", rough));
    presets.back()->normal = "bumpn";

    presets.push_back(make_material("silver_mirror", metal, lgray, mirror));
    presets.push_back(make_material("silver_sharp", metal, lgray, sharp));
    presets.push_back(make_material("silver_rough", metal, lgray, rough));
    presets.push_back(make_material("gold_mirror", metal, gold, mirror));
    presets.push_back(make_material("gold_sharp", metal, gold, sharp));
    presets.push_back(make_material("gold_rough", metal, gold, rough));

    presets.push_back(make_material("transparent_red", transparent, red));
    presets.back()->opacity = 0.9f;
    presets.push_back(make_material("transparent_green", transparent, green));
    presets.back()->opacity = 0.5f;
    presets.push_back(make_material("transparent_blue", transparent, blue));
    presets.back()->opacity = 0.2f;

    presets.push_back(make_material("pointlight", emission, white));
    presets.back()->emission = 160;
    presets.push_back(make_material("arealight", emission, white));
    presets.back()->emission = 80;

    return presets;
}

std::vector<proc_shape*>& proc_shape_presets() {
    static auto presets = std::vector<proc_shape*>();
    if (!presets.empty()) return presets;

    auto make_shape = [](const std::string& name, proc_shape_type type,
                          const vec3i& steps, const vec3f& size,
                          const vec3f& uvsize, int subdivision = 0,
                          bool catmullclark = false, bool faceted = false) {
        auto params = new proc_shape();
        params->name = name;
        params->type = type;
        params->tesselation = steps;
        params->subdivision = subdivision;
        params->size = size;
        params->uvsize = uvsize;
        params->catmull_clark = catmullclark;
        params->faceted = faceted;
        return params;
    };

    presets.push_back(make_shape("floor", proc_shape_type::floor, {64, 64, 0},
        {40, 40, 40}, {20, 20, 20}));
    presets.push_back(make_shape(
        "quad", proc_shape_type::quad, {1, 1, 1}, {2, 2, 2}, {1, 1, 1}));
    presets.push_back(make_shape(
        "cube", proc_shape_type::cube, {1, 1, 1}, {2, 2, 2}, {1, 1, 1}));
    presets.push_back(make_shape("cube_rounded", proc_shape_type::cube_rounded,
        {32, 32, 32}, {2, 2, 2}, {1, 1, 1}));
    presets.back()->rounded = 0.15f;
    presets.push_back(make_shape(
        "sphere", proc_shape_type::sphere, {128, 64, 0}, {2, 2, 2}, {1, 1, 1}));
    presets.push_back(make_shape("sphere_cube", proc_shape_type::sphere_cube,
        {32, 0, 0}, {2, 2, 2}, {1, 1, 1}));
    presets.push_back(make_shape("sphere_flipcap",
        proc_shape_type::sphere_flipcap, {128, 64, 0}, {2, 2, 2}, {1, 1, 1}));
    presets.back()->rounded = 0.75f;
    presets.push_back(make_shape(
        "disk", proc_shape_type::disk, {128, 32, 0}, {2, 2, 2}, {1, 1, 1}));
    presets.push_back(make_shape("disk_quad", proc_shape_type::disk_quad,
        {32, 0, 0}, {2, 2, 2}, {1, 1, 1}));
    presets.push_back(make_shape("disk_bulged", proc_shape_type::disk_bulged,
        {32, 0, 0}, {2, 2, 2}, {1, 1, 1}));
    presets.back()->rounded = 0.25;
    presets.push_back(make_shape("cylinder", proc_shape_type::cylinder,
        {128, 32, 32}, {2, 2, 2}, {1, 1, 1}));
    presets.back()->flip_yz = true;
    presets.push_back(
        make_shape("cylinder_rounded", proc_shape_type::cylinder_rounded,
            {128, 32, 32}, {2, 2, 2}, {1, 1, 1}));
    presets.back()->flip_yz = true;
    presets.back()->rounded = 0.15f;
    presets.push_back(make_shape("geodesic_sphere",
        proc_shape_type::geodesic_sphere, {0, 0, 0}, {2, 2, 2}, {1, 1, 1}, 5));
    presets.push_back(
        make_shape("geodesic_spheref", proc_shape_type::geodesic_sphere,
            {0, 0, 0}, {2, 2, 2}, {1, 1, 1}, 5, 0, true));
    presets.push_back(
        make_shape("geodesic_spherel", proc_shape_type::geodesic_sphere,
            {0, 0, 0}, {2, 2, 2}, {1, 1, 1}, 4, 0, true));
    presets.push_back(make_shape(
        "cubep", proc_shape_type::cubep, {1, 1, 1}, {2, 2, 2}, {1, 1, 1}));
    presets.push_back(make_shape("cubes", proc_shape_type::cubep, {1, 1, 1},
        {2, 2, 2}, {1, 1, 1}, 4, true));
    presets.push_back(make_shape(
        "suzanne", proc_shape_type::suzanne, {0, 0, 0}, {2, 2, 2}, {1, 1, 1}));
    presets.push_back(make_shape("suzannes", proc_shape_type::suzanne,
        {0, 0, 0}, {2, 2, 2}, {1, 1, 1}, 2, true));
    presets.push_back(make_shape(
        "cubefv", proc_shape_type::fvcube, {1, 1, 1}, {2, 2, 2}, {1, 1, 1}));
    presets.push_back(make_shape("cubefvs", proc_shape_type::fvcube, {1, 1, 1},
        {2, 2, 2}, {1, 1, 1}, 0, 4));
    presets.push_back(make_shape("spherefv", proc_shape_type::fvsphere,
        {128, 64, 0}, {2, 2, 2}, {1, 1, 1}));
    presets.push_back(make_shape("matball", proc_shape_type::matball,
        {128, 64, 0}, {2, 2, 2}, {1, 1, 1}));
    presets.back()->rounded = 0.75f;
    presets.push_back(make_shape("matballi", proc_shape_type::sphere,
        {128, 64, 0}, vec3f{2 * 0.8f}, {1, 1, 1}));
    presets.push_back(make_shape("pointscube", proc_shape_type::pointscube,
        {10000, 0, 0}, {2, 2, 2}, {1, 1, 1}));
    presets.push_back(make_shape("hairball1", proc_shape_type::hairball,
        {8, 65536, 0}, {2, 2, 2}, {1, 1, 1}));
    presets.back()->hair_params.radius = {0.001f, 0.0001f};
    presets.back()->hair_params.length = {0.1f, 0.1f};
    presets.back()->hair_params.noise = {0.5f, 8};
    presets.push_back(make_shape("hairball2", proc_shape_type::hairball,
        {8, 65536, 0}, {2, 2, 2}, {1, 1, 1}));
    presets.back()->hair_params.radius = {0.001f, 0.0001f};
    presets.back()->hair_params.length = {0.1f, 0.1f};
    presets.back()->hair_params.clump = {0.5f, 128};
    presets.push_back(make_shape("hairball3", proc_shape_type::hairball,
        {8, 65536, 0}, {2, 2, 2}, {1, 1, 1}));
    presets.back()->hair_params.radius = {0.001f, 0.0001f};
    presets.back()->hair_params.length = {0.1f, 0.1f};
    presets.push_back(make_shape("hairballi", proc_shape_type::sphere,
        {64, 32, 0}, vec3f{2 * 0.8f}, {1, 1, 1}));
    presets.push_back(make_shape("beziercircle", proc_shape_type::beziercircle,
        {0, 0, 0}, {2, 2, 2}, {1, 1, 1}));
    presets.push_back(make_shape(
        "point", proc_shape_type::point, {0, 0, 0}, {2, 2, 2}, {1, 1, 1}));

    return presets;
}

std::vector<proc_environment*>& proc_environment_presets() {
    static auto presets = std::vector<proc_environment*>();
    if (!presets.empty()) return presets;

    auto make_environment = [](const std::string& name,
                                const std::string& texture) {
        auto params = new proc_environment();
        params->name = name;
        params->color = {1, 1, 1};
        params->texture = texture;
        return params;
    };

    presets.push_back(make_environment("const", ""));
    presets.push_back(make_environment("sky1", "sky1"));
    presets.push_back(make_environment("sky2", "sky2"));

    return presets;
}

std::vector<proc_animation*>& proc_animation_presets() {
    static auto presets = std::vector<proc_animation*>();
    if (!presets.empty()) return presets;

    auto make_animation = [](const std::string& name, bool bezier,
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

    presets.push_back(make_animation(
        "bounce", false, {0, 1, 2}, {{0, 0, 0}, {0, 1, 0}, {0, 0, 0}}, {}, {}));
    presets.push_back(make_animation("scale", false, {0, 1, 2}, {}, {},
        {{1, 1, 1}, {0.1f, 0.1f, 0.1f}, {1, 1, 1}}));
    presets.push_back(make_animation("rotation", false, {0, 1, 2}, {},
        {rotation_quat<float>({0, 1, 0}, 0),
            rotation_quat<float>({0, 1, 0}, pif),
            rotation_quat<float>({0, 1, 0}, 0)},
        {}));

    return presets;
}

std::vector<proc_camera*>& proc_camera_presets() {
    static auto presets = std::vector<proc_camera*>();
    if (!presets.empty()) return presets;

    auto make_camera = [](const std::string& name, const vec3f& from,
                           const vec3f& to, float yfov, float aspect) {
        auto params = new proc_camera();
        params->name = name;
        params->from = from;
        params->to = to;
        params->yfov = yfov;
        params->aspect = aspect;
        return params;
    };

    presets.push_back(
        make_camera("cam1", {0, 4, 10}, {0, 1, 0}, 15 * pif / 180, 1));
    presets.push_back(make_camera(
        "cam2", {0, 4, 10}, {0, 1, 0}, 15 * pif / 180, 16.0f / 9.0f));
    presets.push_back(make_camera(
        "cam3", {0, 6, 24}, {0, 1, 0}, 7.5f * pif / 180, 2.35f / 1.0f));

    return presets;
}

std::vector<proc_split_scene*>& proc_split_scene_presets() {
    static auto presets = std::vector<proc_split_scene*>();
    if (!presets.empty()) return presets;

    auto make_split_scene = [](proc_scene* scn,
                                const std::vector<proc_scene*>& views) {
        auto params = new proc_split_scene();
        params->scn = scn;
        params->views = views;
        return params;
    };

    auto make_scene = [](const std::string& name) {
        auto params = new proc_scene();
        params->name = name;
        return params;
    };

    auto make_node = [](const std::string& name, const std::string& cam,
                         const std::string& shape, const std::string& env,
                         const vec3f& pos) {
        auto params = new proc_node();
        params->name = name;
        params->camera = cam;
        params->shape = shape;
        params->environment = env;
        params->frame.o = pos;
        return params;
    };
    auto make_texture = [](const std::string& name) {
        for (auto p : proc_texture_presets())
            if (p->name == name) return new proc_texture(*p);
        throw std::runtime_error("missing texture preset " + name);
    };
    auto make_shape = [](const std::string& name, const vec3f& pos) {
        for (auto p : proc_shape_presets()) {
            if (p->name != name) continue;
            auto shp = new proc_shape(*p);
            shp->frame.o = pos;
            return shp;
        }
        return (proc_shape*)nullptr;
    };
    auto make_environment = [](const std::string& name) {
        for (auto p : proc_environment_presets())
            if (p->name == name) return new proc_environment(*p);
        return (proc_environment*)nullptr;
    };
    auto make_camera = [](const std::string& name) {
        for (auto p : proc_camera_presets())
            if (p->name == name) return new proc_camera(*p);
        return (proc_camera*)nullptr;
    };
    auto make_material = [](const std::string& name) {
        for (auto p : proc_material_presets())
            if (p->name == name) return new proc_material(*p);
        return (proc_material*)nullptr;
    };
    auto make_animation = [](const std::string& name) {
        for (auto p : proc_animation_presets())
            if (p->name == name) return new proc_animation(*p);
        return (proc_animation*)nullptr;
    };

    // textures
    presets.push_back(make_split_scene(make_scene("textures"), {}));
    for (auto txt : proc_texture_presets())
        presets.back()->scn->textures.push_back(make_texture(txt->name));

    // shapes
    presets.push_back(make_split_scene(make_scene("shapes"), {}));
    for (auto shp : proc_shape_presets())
        presets.back()->scn->shapes.push_back(make_shape(shp->name, zero3f));

    // envmap
    presets.push_back(make_split_scene(make_scene("envmaps"), {}));
    for (auto env : proc_environment_presets())
        presets.back()->scn->environments.push_back(
            make_environment(env->name));

    // simple scenes shared functions
    auto make_simple_scene =
        [&](const std::string& name, const std::vector<std::string>& shapes,
            const std::vector<std::string>& mats, bool nodes = false,
            const std::vector<std::string>& animations = {}) {
            auto pos =
                std::vector<vec3f>{{-2.50f, 1, 0}, {0, 1, 0}, {+2.50f, 1, 0}};
            auto preset = make_split_scene(make_scene(name), {});
            auto scn = preset->scn;
            scn->cameras.push_back(make_camera("cam3"));
            scn->materials.push_back(make_material("matte_floor"));
            if (scn->materials.back()->texture != "")
                scn->textures.push_back(
                    make_texture(scn->materials.back()->texture));
            scn->shapes.push_back(make_shape("floor", zero3f));
            scn->shapes.back()->material = scn->materials.back()->name;
            for (auto i = 0; i < shapes.size(); i++) {
                auto name = "obj" + std::to_string(i + 1);
                scn->materials.push_back(make_material(mats[i]));
                scn->shapes.push_back(make_shape(shapes[i], pos[i]));
                scn->shapes.back()->name = name;
                scn->shapes.back()->material = mats[i];
                if (scn->shapes.back()->type == proc_shape_type::hairball ||
                    scn->shapes.back()->type == proc_shape_type::matball) {
                    scn->materials.push_back(make_material("matte_gray"));
                    scn->shapes.back()->interior = "matte_gray";
                }
                if (!animations.empty()) {
                    scn->animations.push_back(make_animation(animations[i]));
                    scn->animations.back()->nodes.push_back(name);
                }
            }

            auto& views = preset->views;
            for (auto lights :
                {"pointlights"s, "arealights"s, "arealights1"s, "envlights"s}) {
                auto view = make_scene(lights);
                view->cameras.push_back(make_camera("cam3"));
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
                        scale = 16;
                    }
                    if (lights == "arealights1") {
                        emission = 80;
                        shp = "quad";
                        mat = "arealight";
                        pos = {{-4, 5, 8}, {+4, 5, 8}};
                        scale = 4;
                    }
                    for (auto i = 0; i < 2; i++) {
                        auto name = "light" + std::to_string(i + 1);
                        view->materials.push_back(make_material(mat));
                        view->materials.back()->name = name;
                        view->materials.back()->emission = emission;
                        view->shapes.push_back(make_shape(shp, pos[i]));
                        view->shapes.back()->name = name;
                        view->shapes.back()->material = name;
                        view->shapes.back()->size = {scale, scale, scale};
                        if (lights == "arealights" || lights == "arealights1")
                            view->shapes.back()->frame = lookat_frame(
                                pos[i], {0, 1, 0}, {0, 0, 1}, true);
                    }
                }
                if (lights == "envlights") {
                    view->environments.push_back(make_environment("sky1"));
                    view->textures.push_back(make_texture("sky1"));
                }
                views.push_back(view);
            }
            auto scns = std::vector<proc_scene*>{scn};
            append(scns, views);
            for (auto scn : scns) {
                if (!animations.empty() || nodes) {
                    for (auto cam : scn->cameras) {
                        auto nde = new proc_node();
                        nde->name = cam->name;
                        nde->frame =
                            lookat_frame(cam->from, cam->to, vec3f{0, 1, 0});
                        nde->camera = cam->name;
                        scn->nodes.push_back(nde);
                    }
                    for (auto shp : scn->shapes) {
                        auto nde = new proc_node();
                        nde->name = shp->name;
                        nde->frame = shp->frame;
                        nde->shape = shp->name;
                        scn->nodes.push_back(nde);
                    }
                    for (auto env : scn->environments) {
                        auto nde = new proc_node();
                        nde->name = env->name;
                        nde->frame = env->frame;
                        nde->environment = env->name;
                        scn->nodes.push_back(nde);
                    }
                }
            }
            return preset;
        };

    // plane only
    presets.push_back(make_simple_scene("plane", {}, {}));

    // basic shapes
    presets.push_back(make_simple_scene("basic",
        {"sphere_flipcap", "sphere_cube", "cube_rounded"},
        {"plastic_red", "plastic_green", "plastic_blue"}));

    // simple shapes
    presets.push_back(make_simple_scene("simple",
        {"sphere_flipcap", "sphere_cube", "cube_rounded"},
        {"plastic_colored", "plastic_colored", "plastic_colored"}));

    // simple shapes 1
    presets.push_back(make_simple_scene("spheres",
        {"sphere_flipcap", "sphere_cube", "sphere"},
        {"plastic_colored", "plastic_colored", "plastic_colored"}));

    // simple shapes 2
    presets.push_back(
        make_simple_scene("cubes", {"quad", "cube", "cube_rounded"},
            {"plastic_colored", "plastic_colored", "plastic_colored"}));

    // simple shapes 3
    presets.push_back(
        make_simple_scene("cylinders", {"disk", "cylinder", "cylinder_rounded"},
            {"plastic_colored", "plastic_colored", "plastic_colored"}));

    // simple shapes 3
    presets.push_back(
        make_simple_scene("disks", {"disk", "disk_quad", "disk_bulged"},
            {"plastic_colored", "plastic_colored", "plastic_colored"}));

    // transparent shapes
    presets.push_back(make_simple_scene("transparent", {"quad", "quad", "quad"},
        {"transparent_red", "transparent_green", "transparent_blue"}));

    // lines shapes
    presets.push_back(
        make_simple_scene("lines", {"hairball1", "hairball2", "hairball3"},
            {"matte_gray", "matte_gray", "matte_gray"}));

    // subdiv shapes
    presets.push_back(
        make_simple_scene("subdiv", {"cubes", "suzannes", "suzannes"},
            {"plastic_red", "plastic_green", "plastic_blue"}));

    // plastics shapes
    presets.push_back(
        make_simple_scene("plastics", {"matball", "matball", "matball"},
            {"matte_green", "plastic_green", "plastic_colored"}, true));

    // metals shapes
    presets.push_back(
        make_simple_scene("metals", {"matball", "matball", "matball"},
            {"gold_rough", "gold_sharp", "silver_mirror"}));

    // tesselation shapes
    presets.push_back(make_simple_scene("tesselation",
        {"geodesic_spherel", "geodesic_spheref", "geodesic_sphere"},
        {"matte_gray", "matte_gray", "matte_gray"}));

    // textureuv shapes
    presets.push_back(make_simple_scene("textureuv",
        {"sphere_flipcap", "sphere_flipcap", "sphere_flipcap"},
        {"matte_green", "matte_colored", "matte_uv"}));

    // normalmap shapes
    presets.push_back(make_simple_scene("normalmap",
        {"sphere_flipcap", "sphere_flipcap", "sphere_flipcap"},
        {"plastic_blue", "plastic_blue_bumped", "plastic_colored_bumped"}));

    // animated shapes
    presets.push_back(make_simple_scene("animated",
        {"sphere_flipcap", "sphere_cube", "cube_rounded"},
        {"plastic_colored", "plastic_colored", "plastic_colored"}, true,
        {"bounce", "scale", "rotation"}));

    // instances shared functions
    auto make_random_scene = [&](const std::string& name, const vec2i& num,
                                 const bbox2f& bbox, uint64_t seed = 13) {
        auto rscale = 0.9f * 0.25f *
                      min((bbox.max.x - bbox.min.x) / num.x,
                          (bbox.max.x - bbox.min.x) / num.y);

        auto preset = make_split_scene(make_scene(name), {});
        auto scn = preset->scn;
        scn->cameras.push_back(make_camera("cam3"));
        scn->materials.push_back(make_material("matte_floor"));
        scn->shapes.push_back(make_shape("floor", zero3f));
        scn->shapes.back()->material = scn->materials.back()->name;
        auto shapes = std::vector<std::string>();
        for (auto mat : {"plastic_red", "plastic_green", "plastic_blue"})
            scn->materials.push_back(make_material(mat));
        for (auto shp : {"sphere", "sphere_flipcap", "cube"}) {
            for (auto mat : {"plastic_red", "plastic_green", "plastic_blue"}) {
                scn->shapes.push_back(make_shape(shp, zero3f));
                scn->shapes.back()->name += "_"s + mat;
                scn->shapes.back()->material = mat;
                scn->shapes.back()->size *= 2 * rscale;
                shapes.push_back(scn->shapes.back()->name);
            }
        }

        auto rng = make_rng(seed, 7);
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
                scn->nodes.push_back(
                    make_node("instance" + std::to_string(count++), "",
                        shapes[next_rand1i(rng, (int)shapes.size())], "", pos));
            }
        }

        auto& views = preset->views;
        for (auto name : {"pointlights"s}) {
            auto view = make_scene(name);
            view->cameras.push_back(make_camera("cam3"));
            auto pos = std::vector<vec3f>{{-2, 10, 8}, {+2, 10, 8}};
            for (auto i = 0; i < 2; i++) {
                auto name = "light" + std::to_string(i + 1);
                view->materials.push_back(make_material("pointlight"));
                view->materials.back()->name = name;
                view->materials.back()->emission = 80;
                view->shapes.push_back(make_shape("point", zero3f));
                view->shapes.back()->name = name;
                view->shapes.back()->material = name;
                view->nodes.push_back(make_node(name, "", name, "", pos[i]));
            }
            views.push_back(view);
        }

        return preset;
    };

    // instances
    presets.push_back(
        make_random_scene("instances", {10, 10}, {{-3, -3}, {3, 3}}));
    presets.push_back(
        make_random_scene("instancel", {100, 100}, {{-3, -3}, {3, 3}}));

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
    auto add_textures = [&](proc_scene* preset) {
        auto used = std::unordered_set<std::string>();
        for (auto mat : preset->materials) used.insert(mat->texture);
        for (auto mat : preset->materials) used.insert(mat->normal);
        for (auto env : preset->environments) used.insert(env->texture);
        used.erase("");
        for (auto txt : preset->textures) used.erase(txt->name);
        for (auto txt : used) preset->textures.push_back(make_texture(txt));
    };

    // add missing textures
    for (auto& preset : presets) {
        add_textures(preset->scn);
        for (auto& view : preset->views) add_textures(view);
    }

    // remove duplicates
    for (auto& preset : presets) {
        remove_duplicates(preset->scn);
        for (auto& view : preset->views) remove_duplicates(view);
    }

    return presets;
}

std::vector<proc_scene*>& proc_scene_presets() {
    static auto presets = std::vector<proc_scene*>();
    if (!presets.empty()) return presets;

    for (auto sscn : proc_split_scene_presets()) {
        auto scn = sscn->scn;
        for (auto view : sscn->views) {
            auto preset = new proc_scene();
            preset->name = scn->name + "_" + view->name;
            for (auto mscn : {scn, view}) {
                for (auto v : mscn->cameras)
                    preset->cameras.push_back(new proc_camera(*v));
                for (auto v : mscn->textures)
                    preset->textures.push_back(new proc_texture(*v));
                for (auto v : mscn->materials)
                    preset->materials.push_back(new proc_material(*v));
                for (auto v : mscn->shapes)
                    preset->shapes.push_back(new proc_shape(*v));
                for (auto v : mscn->environments)
                    preset->environments.push_back(new proc_environment(*v));
                for (auto v : mscn->nodes)
                    preset->nodes.push_back(new proc_node(*v));
                for (auto v : mscn->animations)
                    preset->animations.push_back(new proc_animation(*v));
            }
            presets.push_back(preset);
        }
    }

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
    serialize_attr(val.catmull_clark, js, "catmull_clark", reading, false,
        def.catmull_clark);
    serialize_attr(val.size, js, "size", reading, false, def.size);
    serialize_attr(val.uvsize, js, "uvsize", reading, false, def.uvsize);
    serialize_attr(val.radius, js, "radius", reading, false, def.radius);
    serialize_attr(val.faceted, js, "faceted", reading, false, def.faceted);
    serialize_attr(val.catmull_clark, js, "catmull_clark", reading, false,
        def.catmull_clark);
    // TODO: hair parameters
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

// Parses a test_node object
void serialize(proc_node& val, json& js, bool reading) {
    static auto def = proc_node();
    serialize_obj(js, reading);
    serialize_attr(val.name, js, "name", reading);
    serialize_attr(val.parent, js, "parent", reading, false, def.parent);
    serialize_attr(val.camera, js, "camera", reading, false, def.camera);
    serialize_attr(val.shape, js, "shape", reading, false, def.shape);
    serialize_attr(
        val.environment, js, "environment", reading, false, def.environment);
    serialize_attr(val.frame, js, "frame", reading, false, def.frame);
    serialize_attr(
        val.translation, js, "translation", reading, false, def.translation);
    serialize_attr(val.rotation, js, "rotation", reading, false, def.rotation);
    serialize_attr(val.scaling, js, "scaling", reading, false, def.scaling);
}

// Parses a test_node object
void serialize(proc_animation& val, json& js, bool reading) {
    static auto def = proc_animation();
    serialize_obj(js, reading);
    serialize_attr(val.name, js, "name", reading);
    serialize_attr(val.bezier, js, "bezier", reading, false, def.bezier);
    serialize_attr(val.speed, js, "speed", reading, false, def.speed);
    serialize_attr(val.scale, js, "scale", reading, false, def.scale);
    serialize_attr(val.times, js, "times", reading, false, def.times);
    serialize_attr(
        val.translation, js, "translation", reading, false, def.translation);
    serialize_attr(val.rotation, js, "rotation", reading, false, def.rotation);
    serialize_attr(val.scaling, js, "scaling", reading, false, def.scaling);
    serialize_attr(val.nodes, js, "nodes", reading, false, def.nodes);
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
    serialize_attr(val.nodes, js, "nodes", reading, false, def.nodes);
    serialize_attr(
        val.animations, js, "animations", reading, false, def.animations);
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
    const gl_program& prog, int pos, const int* val, int ncomp, int count) {
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
    const gl_program& prog, int pos, const float* val, int ncomp, int count) {
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
    const gl_program& prog, int pos, const gl_texture_info& tinfo, uint tunit) {
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

// Set uniform texture id tid and unit tunit for program pid and
// variable var.
bool set_program_uniform_texture(
    const gl_program& prog, int pos, const gl_texture& txt, uint tunit) {
    assert(gl_check_error());
    if (pos < 0) return false;
    if (is_texture_valid(txt)) {
        bind_texture(txt, tunit);
        glUniform1i(pos, tunit);
    } else {
        unbind_texture(txt, tunit);
        glUniform1i(pos, tunit);
    }
    assert(gl_check_error());
    return true;
}

// Set uniform texture id tid and unit tunit for program pid and
// variable var.
bool set_program_uniform_texture(const gl_program& prog,
    const std::string& name, const gl_texture& txt, uint tunit) {
    return set_program_uniform_texture(
        prog, get_program_uniform_location(prog, name), txt, tunit);
}

// Sets a constant value for a vertex attribute for program pid and
// variable var. The attribute has nc components.
bool set_program_vertattr(
    const gl_program& prog, int pos, const float* value, int nc) {
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
bool set_program_vertattr(
    const gl_program& prog, int pos, const int* value, int nc) {
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
bool set_program_vertattr(const gl_program& prog, const std::string& var,
    const gl_vertex_buffer& buf) {
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
bool set_program_vertattr(const gl_program& prog, int pos,
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

// Init shading
void update_gl_texture(const texture* txt, gl_texture& gtxt) {
    if (!txt) {
        clear_texture(gtxt);
    } else {
        if (!txt->hdr.empty()) {
            update_texture(gtxt, txt->hdr, true, true, true);
        } else if (!txt->ldr.empty()) {
            update_texture(gtxt, txt->ldr, true, true, true);
        } else
            assert(false);
    }
}

// Update shading
void update_gl_shape(const shape* shp, gl_shape& gshp) {
    auto update_vert_buffer = [](auto& buf, const auto& vert) {
        if (vert.empty()) {
            clear_vertex_buffer(buf);
        } else {
            update_vertex_buffer(buf, vert);
        }
    };
    auto update_elem_buffer = [](auto& buf, const auto& elem) {
        if (elem.empty()) {
            clear_element_buffer(buf);
        } else {
            update_element_buffer(buf, elem);
        }
    };
    if (!shp) {
        clear_gl_shape(gshp);
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
            update_elem_buffer(gshp.quads, convert_quads_to_triangles(quads));
            update_elem_buffer(gshp.edges, get_edges({}, {}, shp->quads));
            update_vert_buffer(gshp.color, std::vector<vec4f>{});
            update_vert_buffer(gshp.tangsp, std::vector<vec4f>{});
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

// clear gl_shape
void clear_gl_shape(gl_shape& gshp) {
    clear_vertex_buffer(gshp.pos);
    clear_vertex_buffer(gshp.norm);
    clear_vertex_buffer(gshp.texcoord);
    clear_vertex_buffer(gshp.color);
    clear_vertex_buffer(gshp.tangsp);
    clear_element_buffer(gshp.points);
    clear_element_buffer(gshp.lines);
    clear_element_buffer(gshp.triangles);
    clear_element_buffer(gshp.quads);
    clear_element_buffer(gshp.beziers);
    clear_element_buffer(gshp.edges);
}

// Add gl lights
void add_gl_lights(gl_lights& lights, const frame3f& frame, const shape* shp) {
    if (!has_emission(shp)) return;
    if (lights.pos.size() >= 16) return;
    if (!shp->points.empty()) {
        for (auto p : shp->points) {
            if (lights.pos.size() >= 16) break;
            lights.pos.push_back(transform_point(frame, shp->pos[p]));
            lights.ke.push_back(shp->groups.at(0).mat->ke);
            lights.type.push_back(gl_light_type::point);
        }
    } else {
        auto bbox = make_bbox(shp->pos.size(), shp->pos.data());
        auto pos = bbox_center(bbox);
        auto area = 0.0f;
        for (auto l : shp->lines)
            area += line_length(shp->pos[l.x], shp->pos[l.y]);
        for (auto t : shp->triangles)
            area += triangle_area(shp->pos[t.x], shp->pos[t.y], shp->pos[t.z]);
        for (auto t : shp->quads)
            area += quad_area(
                shp->pos[t.x], shp->pos[t.y], shp->pos[t.z], shp->pos[t.w]);
        auto ke = shp->groups.at(0).mat->ke * area;
        if (lights.pos.size() < 16) {
            lights.pos.push_back(transform_point(frame, pos));
            lights.ke.push_back(ke);
            lights.type.push_back(gl_light_type::point);
        }
    }
}

// Initialize gl lights
gl_lights make_gl_lights(const scene* scn) {
    auto lights = gl_lights();
    if (scn->nodes.empty()) {
        for (auto shp : scn->shapes) { add_gl_lights(lights, shp->frame, shp); }
    } else {
        for (auto nde : scn->nodes) {
            if (!nde->shp) continue;
            add_gl_lights(lights, nde->frame_, nde->shp);
        }
    }
    return lights;
}

}  // namespace ygl

// -----------------------------------------------------------------------------
// IMPLEMENTATION FOR OPRNGL IMAGE SHADER FUNCTIONS
// -----------------------------------------------------------------------------
namespace ygl {

// Initialize the program. Call with true only after the GL is initialized.
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
void draw_image(const gl_stdimage_program& prog, const gl_texture& txt,
    const vec2i& win_size, const vec2f& offset, float zoom, float exposure,
    tonemap_type tonemap) {
    assert(is_texture_valid(txt));

    bind_program(prog.prog);

    gl_enable_blending(true);
    gl_set_blend_over();

    bind_texture(txt, 0);
    set_program_uniform(prog.prog, "zoom", zoom);
    set_program_uniform(
        prog.prog, "win_size", vec2f{(float)win_size.x, (float)win_size.y});
    set_program_uniform(prog.prog, "offset", offset);
    set_program_uniform(prog.prog, "tonemap.exposure", exposure);
    set_program_uniform(prog.prog, "tonemap.type", (int)tonemap);
    set_program_uniform_texture(prog.prog, "img", txt, 0);

    set_program_vertattr(prog.prog, "vert_texcoord", prog.vbo, vec2f{0, 0});
    draw_elems(prog.ebo);

    unbind_program(prog.prog);

    gl_enable_blending(false);

    assert(gl_check_error());
}

// Computes the image uv coordinates corresponding to the view parameters.
vec2i get_draw_image_coords(
    const vec2f& mouse_pos, const gl_stdimage_params& params) {
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
void begin_stdsurface_frame(const gl_stdsurface_program& prog,
    bool shade_eyelight, float exposure, tonemap_type tonemap,
    const mat4f& camera_xform, const mat4f& camera_xform_inv,
    const mat4f& camera_proj) {
    static auto eyelight_id =
        get_program_uniform_location(prog.prog, "lighting.eyelight");
    static auto exposure_id =
        get_program_uniform_location(prog.prog, "tonemap.exposure");
    static auto tonemap_id =
        get_program_uniform_location(prog.prog, "tonemap.type");
    static auto xform_id =
        get_program_uniform_location(prog.prog, "camera.xform");
    static auto xform_inv_id =
        get_program_uniform_location(prog.prog, "camera.xform_inv");
    static auto proj_id =
        get_program_uniform_location(prog.prog, "camera.proj");
    assert(gl_check_error());
    bind_program(prog.prog);
    set_program_uniform(prog.prog, eyelight_id, shade_eyelight);
    set_program_uniform(prog.prog, exposure_id, exposure);
    set_program_uniform(prog.prog, tonemap_id, (int)tonemap);
    set_program_uniform(prog.prog, xform_id, camera_xform);
    set_program_uniform(prog.prog, xform_inv_id, camera_xform_inv);
    set_program_uniform(prog.prog, proj_id, camera_proj);
    assert(gl_check_error());
}

// Ends a frame.
void end_stdsurface_frame(const gl_stdsurface_program& prog) {
    assert(gl_check_error());
    unbind_program(prog.prog);
    //    glBindVertexArray(0);
    //    glUseProgram(0);
    assert(gl_check_error());
}

// Set num lights with position pos, color ke, type ltype. Also set the
// ambient illumination amb.
void set_stdsurface_lights(const gl_stdsurface_program& prog, const vec3f& amb,
    const gl_lights& lights) {
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
void begin_stdsurface_shape(const gl_stdsurface_program& prog,
    const mat4f& xform, float normal_offset) {
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
void end_stdsurface_shape(const gl_stdsurface_program& prog) {
    assert(gl_check_error());
    for (int i = 0; i < 16; i++) unbind_vertex_buffer(i);
    assert(gl_check_error());
}

// Sets normal offset.
void set_stdsurface_normaloffset(
    const gl_stdsurface_program& prog, float normal_offset) {
    static auto normal_offset_id =
        get_program_uniform_location(prog.prog, "shape_normal_offset");
    assert(gl_check_error());
    set_program_uniform(prog.prog, normal_offset_id, normal_offset);
    assert(gl_check_error());
}

// Set the object as highlighted.
void set_stdsurface_highlight(
    const gl_stdsurface_program& prog, const vec4f& highlight) {
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
void set_stdsurface_material(const gl_stdsurface_program& prog,
    material_type type, gl_elem_type etype, const vec3f& ke, const vec3f& kd,
    const vec3f& ks, float rs, float op, const gl_texture_info& ke_txt,
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
void set_stdsurface_vert(const gl_stdsurface_program& prog,
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
void set_stdsurface_vert_skinning(const gl_stdsurface_program& prog,
    gl_vertex_buffer& weights, gl_vertex_buffer& joints, int nxforms,
    const mat4f* xforms) {
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
void set_stdsurface_vert_gltf_skinning(const gl_stdsurface_program& prog,
    gl_vertex_buffer& weights, gl_vertex_buffer& joints, int nxforms,
    const mat4f* xforms) {
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
void set_stdsurface_vert_skinning_off(const gl_stdsurface_program& prog) {
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
    static auto default_material = new material();
    default_material->kd = {0.2f, 0.2f, 0.2f};

    begin_stdsurface_shape(prog, xform);

    auto etype = gl_elem_type::triangle;
    if (!shp->lines.empty()) etype = gl_elem_type::line;
    if (!shp->points.empty()) etype = gl_elem_type::point;

    auto txt = [&textures](const texture_info& info) -> gl_texture_info {
        auto ginfo = gl_texture_info();
        if (!info.txt) return ginfo;
        ginfo.txt = textures.at(info.txt);
        return ginfo;
    };

    auto mat = get_material(shp, 0) ? get_material(shp, 0) : default_material;
    set_stdsurface_material(prog, mat->type, etype, mat->ke, mat->kd, mat->ks,
        mat->rs, mat->op, txt(mat->ke_txt), txt(mat->kd_txt), txt(mat->ks_txt),
        txt(mat->rs_txt), txt(mat->norm_txt), txt(mat->occ_txt), false,
        mat->double_sided || params.double_sided, params.cutout);

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
    const vec2i& viewport_size, const void* highlighted,
    const gl_stdsurface_params& params, const tonemap_params& tmparams) {
    // begin frame
    gl_enable_depth_test(true);
    gl_enable_culling(params.cull_backface && !params.double_sided);
    gl_enable_wireframe(params.wireframe);
    gl_set_viewport(viewport_size);

    auto camera_xform = frame_to_mat(cam->frame);
    auto camera_view = frame_to_mat(inverse(cam->frame));
    auto camera_proj = perspective_mat(cam->yfov,
        (float)viewport_size.x / (float)viewport_size.y, cam->near, cam->far);

    begin_stdsurface_frame(prog, params.eyelight, tmparams.exposure,
        tmparams.type, camera_xform, camera_view, camera_proj);

    if (!params.eyelight) {
        set_stdsurface_lights(prog, params.ambient, lights);
    }

    if (scn->nodes.empty()) {
        for (auto shp : scn->shapes) {
            draw_stdsurface_shape(shp, frame_to_mat(shp->frame),
                shp == highlighted, prog, shapes, textures, params);
        }
    } else {
        for (auto nde : scn->nodes) {
            if (!nde->shp) continue;
            draw_stdsurface_shape(nde->shp, frame_to_mat(nde->frame_),
                nde == highlighted || nde->shp == highlighted, prog, shapes,
                textures, params);
        }
    }

    end_stdsurface_frame(prog);
    gl_enable_wireframe(false);
}

gl_window::~gl_window() {
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
    printf("GLFW error: %s\n", description);
}

// Support
void _glfw_text_cb(GLFWwindow* gwin, unsigned key) {
    auto win = (gl_window*)glfwGetWindowUserPointer(gwin);
    if (win->widget_enabled) { ImGui_ImplGlfwGL3_CharCallback(win->gwin, key); }
    if (win->text_cb) win->text_cb(key);
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
    if (win->mouse_cb) win->mouse_cb(button, action == GLFW_PRESS, mods);
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
    if (win->refresh_cb) win->refresh_cb();
}

// Initialize gl_window
gl_window* make_window(
    int width, int height, const std::string& title, bool opengl4) {
    auto win = new gl_window();

    // gl_window
    if (!glfwInit()) throw std::runtime_error("cannot open gl_window");

    // profile creation
    if (opengl4) {
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    } else {
        glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
        glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    }
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

// Check if the alt key is down
bool get_alt_key(gl_window* win) {
    return glfwGetKey(win->gwin, GLFW_KEY_LEFT_ALT) == GLFW_PRESS ||
           glfwGetKey(win->gwin, GLFW_KEY_RIGHT_ALT) == GLFW_PRESS;
}

// Check if the alt key is down
bool get_ctrl_key(gl_window* win) {
    return glfwGetKey(win->gwin, GLFW_KEY_LEFT_CONTROL) == GLFW_PRESS ||
           glfwGetKey(win->gwin, GLFW_KEY_RIGHT_CONTROL) == GLFW_PRESS;
}

// Check if the alt key is down
bool get_shift_key(gl_window* win) {
    return glfwGetKey(win->gwin, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS ||
           glfwGetKey(win->gwin, GLFW_KEY_RIGHT_SHIFT) == GLFW_PRESS;
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
    auto alt_down = get_alt_key(win);

    // updated
    auto updated = false;

    // handle mouse and keyboard for navigation
    if (mouse_button && alt_down && !get_widget_active(win)) {
        if (mouse_button == 1 && get_shift_key(win)) mouse_button = 3;
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

// Handle scene selection
bool handle_scene_selection(gl_window* win, const scene* scn, const camera* cam,
    const bvh_tree* bvh, int res, const gl_stdimage_params& params,
    scene_selection& sel) {
    auto mouse_pos = get_mouse_posf(win);
    auto mouse_button = get_mouse_button(win);

    if (!(mouse_button == 1 && !get_widget_active(win))) return false;
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
bool draw_multilinetext_widget(
    gl_window* win, const std::string& lbl, std::string& str) {
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
    std::cout << str;
    auto ret = ImGui::InputTextMultiline(lbl.c_str(), buf, buf_size);
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
void draw_image_widget(
    gl_window* win, const gl_texture& txt, const vec2i& size) {
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
void draw_groupid_widget_push(gl_window* win, int gid) { ImGui::PushID(gid); }
// Group ids
void draw_groupid_widget_push(gl_window* win, void* gid) { ImGui::PushID(gid); }
// Group ids
void draw_groupid_widget_push(gl_window* win, const void* gid) {
    ImGui::PushID(gid);
}
// Group ids
void draw_groupid_widget_push(gl_window* win, const char* gid) {
    ImGui::PushID(gid);
}
// Group ids
void draw_groupid_widget_pop(gl_window* win) { ImGui::PopID(); }

// Widget style
void draw_style_widget_push(gl_window* win, const vec4f& col) {
    ImGui::PushStyleColor(ImGuiCol_Text, {col.x, col.y, col.z, col.w});
}

// Widget style
void draw_style_widget_pop(gl_window* win) { ImGui::PopStyleColor(); }

// Image inspection widgets.
void draw_imageinspect_widgets(gl_window* win, const std::string& lbl,
    const image4f& hdr, const image4b& ldr, const vec2f& mouse_pos,
    const gl_stdimage_params& params) {
    auto ij = get_draw_image_coords(mouse_pos, params);
    auto v4f = zero4f;
    auto v4b = zero4b;
    if (!hdr.empty() && contains(hdr, ij.x, ij.y)) {
        v4f = hdr.at(ij.x, ij.y);
        v4b = linear_to_srgb(hdr.at(ij.x, ij.y));
    } else if (!ldr.empty() && contains(ldr, ij.x, ij.y)) {
        v4f = srgb_to_linear(ldr.at(ij.x, ij.y));
        v4b = ldr.at(ij.x, ij.y);
    }
    char buf[1024];
    sprintf(buf, "%5d %5d", ij.x, ij.y);
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

// template <typename T>
// const std::vector<T*>& elems(const T* aux) {
//     static auto v = std::vector<T*>();
//     return v;
// }

const std::vector<texture*>& get_scene_elems(
    const scene* scn, const texture* aux) {
    return scn->textures;
}
const std::vector<shape*>& get_scene_elems(const scene* scn, shape* aux) {
    return scn->shapes;
}
const std::vector<camera*> get_scene_elems(const scene* scn, camera* aux) {
    return scn->cameras;
}
const std::vector<material*> get_scene_elems(const scene* scn, material* aux) {
    return scn->materials;
}
const std::vector<environment*> get_scene_elems(
    const scene* scn, environment* aux) {
    return scn->environments;
}
const std::vector<node*> get_scene_elems(const scene* scn, node* aux) {
    return scn->nodes;
}
const std::vector<animation_group*> get_scene_elems(
    const scene* scn, animation_group* aux) {
    return scn->animations;
}

// Implementation of camera selection
bool draw_camera_widgets(gl_window* win, const std::string& lbl, camera* cam) {
    if (!cam) return false;
    auto edited = false;
    auto from = cam->frame.o;
    auto to = cam->frame.o - cam->frame.z * cam->focus;
    edited += draw_value_widget(win, lbl + " name", cam->name);
    edited += draw_value_widget(win, lbl + " from", from, -10, 10);
    edited += draw_value_widget(win, lbl + " to", to, -10, 10);
    edited += draw_value_widget(win, lbl + " yfov", cam->yfov, 0.1, 10);
    edited += draw_value_widget(win, lbl + " aspect", cam->aspect, 1, 3);
    edited += draw_value_widget(win, lbl + " aperture", cam->aperture, 0, 1);
    edited += draw_value_widget(win, lbl + " focus", cam->focus, 0.1, 100);
    edited += draw_value_widget(win, lbl + " near", cam->near, 0.01f, 1);
    edited += draw_value_widget(win, lbl + " far", cam->far, 1, 1000);
    if (edited) cam->frame = lookat_frame(from, to, {0, 1, 0});
    return edited;
}

static const std::unordered_map<std::string, vec4f>
    draw_visitor_highlight_colors = {{"red", {1, 0.5f, 0.5f, 1}},
        {"green", {0.5f, 1, 0.5f, 1}}, {"blue", {0.5f, 0.5f, 1, 1}}};

vec4f get_highlight_color(
    const std::unordered_map<std::string, std::string>& highlights,
    const std::string& name) {
    if (!highlights.empty() && contains(highlights, name)) {
        return draw_visitor_highlight_colors.at(highlights.at(name));
    }
    return zero4f;
}

template <typename T>
void draw_scene_tree_widgets(gl_window* win, const std::string& lbl_, T& val,
    scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {}

template <typename T>
void draw_scene_tree_widgets(gl_window* win, const std::string& lbl_, T* val,
    scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    if (!val) return;
    auto lbl = val->name;
    if (!lbl_.empty()) lbl = lbl_ + ": " + val->name;
    auto selection = sel.get_raw();
    auto color = get_highlight_color(highlights, val->name);
    if (color != zero4f) draw_style_widget_push(win, color);
    auto open = draw_tree_widget_begin(win, lbl, selection, val);
    if (color != zero4f) draw_style_widget_pop(win);
    if (selection == val) sel = val;
    if (open) {
        visit(val, [win, &sel, &highlights](auto& val, const auto& var) {
            draw_scene_tree_widgets(win, var.name, val, sel, highlights);
        });
        draw_tree_widget_end(win);
    }
}

void draw_scene_tree_widgets(gl_window* win, const std::string& lbl,
    texture_info& val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    draw_scene_tree_widgets(win, lbl, val.txt, sel, highlights);
}

void draw_scene_tree_widgets(gl_window* win, const std::string& lbl,
    std::vector<shape_group_props>& val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    for (auto gid = 0; gid < val.size(); gid++) {
        draw_scene_tree_widgets(win, "mat" + std::to_string(gid + 1),
            val[gid].mat, sel, highlights);
    }
}

template <typename T>
void draw_scene_tree_widgets(gl_window* win, const std::string& lbl,
    std::vector<T*>& val, scene_selection& sel,
    const std::unordered_map<std::string, std::string>& highlights) {
    if (draw_tree_widget_begin(win, lbl)) {
        for (auto v : val) draw_scene_tree_widgets(win, "", v, sel, highlights);
        draw_tree_widget_end(win);
    }
}

template <typename T>
bool draw_scene_elem_widgets(gl_window* win, T& val, const visit_var& var,
    scene* scn, const std::unordered_map<std::string, std::string>& highlights,
    const std::string& elem_name) {
    auto color = get_highlight_color(highlights, elem_name + "___" + var.name);
    if (color != zero4f) draw_style_widget_push(win, color);
    auto edited = 0;
    edited += draw_value_widget(win, var.name, val, var.min, var.max,
        var.type == visit_var_type::color);
    if (color != zero4f) draw_style_widget_pop(win);
    return edited;
}

bool draw_scene_elem_widgets(gl_window* win, make_hair_params& val,
    const visit_var& var, scene* scn,
    const std::unordered_map<std::string, std::string>& highlights,
    const std::string& elem_name) {
    // TODO: hair params
    return false;
}

bool draw_scene_elem_widgets(gl_window* win, texture_info& val,
    const visit_var& var, scene* scn,
    const std::unordered_map<std::string, std::string>& highlights,
    const std::string& elem_name) {
    auto color = get_highlight_color(highlights, elem_name + "___" + var.name);
    if (color != zero4f) draw_style_widget_push(win, color);
    auto edited = 0;
    edited += draw_combo_widget(
        win, var.name, val.txt, get_scene_elems(scn, val.txt));
    if (val.txt) {
        edited += draw_value_widget(win, var.name + " wrap s", val.wrap_s,
            var.min, var.max, var.type == visit_var_type::color);
        draw_continue_widget(win);
        edited += draw_value_widget(win, var.name + " wrap t", val.wrap_t,
            var.min, var.max, var.type == visit_var_type::color);
        edited += draw_value_widget(win, var.name + " linear", val.linear,
            var.min, var.max, var.type == visit_var_type::color);
        draw_continue_widget(win);
        edited += draw_value_widget(win, var.name + " mipmap", val.mipmap,
            var.min, var.max, var.type == visit_var_type::color);
    }
    if (color != zero4f) draw_style_widget_pop(win);
    return edited;
}

bool draw_scene_elem_widgets(gl_window* win,
    std::vector<shape_group_props>& val, const visit_var& var, scene* scn,
    const std::unordered_map<std::string, std::string>& highlights,
    const std::string& elem_name) {
    auto color = get_highlight_color(highlights, elem_name + "___" + var.name);
    if (color != zero4f) draw_style_widget_push(win, color);
    auto edited = 0;
    for (auto gid = 0; gid < val.size(); gid++) {
        auto& grp = val[gid];
        auto name = var.name + "[" + std::to_string(gid) + "].";
        visit(grp, [win, name, scn, &highlights, &elem_name, &edited](
                       auto& val, const auto& var) {
            auto var_ = var;
            var_.name = name + var.name;
            edited += draw_scene_elem_widgets(
                win, val, var_, scn, highlights, elem_name);
        });
    }
    if (draw_button_widget(win, "add " + var.name)) {
        val.push_back({});
        edited += 1;
    }
    draw_continue_widget(win);
    if (draw_button_widget(win, "del " + var.name)) {
        if (!val.empty()) val.pop_back();
        edited += 1;
    }
    if (color != zero4f) draw_style_widget_pop(win);
    return edited;
}

bool draw_scene_elem_widgets(gl_window* win, std::vector<animation>& val,
    const visit_var& var, scene* scn,
    const std::unordered_map<std::string, std::string>& highlights,
    const std::string& elem_name) {
    auto color = get_highlight_color(highlights, elem_name + "___" + var.name);
    if (color != zero4f) draw_style_widget_push(win, color);
    auto edited = 0;
    for (auto gid = 0; gid < val.size(); gid++) {
        auto& grp = val[gid];
        auto name = var.name + "[" + std::to_string(gid) + "].";
        visit(grp, [win, name, scn, &highlights, &elem_name, &edited](
                       auto& val, const auto& var) {
            auto var_ = var;
            var_.name = name + var.name;
            edited += draw_scene_elem_widgets(
                win, val, var_, scn, highlights, elem_name);
        });
    }
    if (color != zero4f) draw_style_widget_pop(win);
    return edited;
}

template <typename T>
bool draw_scene_elem_widgets(gl_window* win, image<T>& val,
    const visit_var& var, scene* scn,
    const std::unordered_map<std::string, std::string>& highlights,
    const std::string& elem_name) {
    if (empty(val)) return false;
    auto color = get_highlight_color(highlights, elem_name + "___" + var.name);
    if (color != zero4f) draw_style_widget_push(win, color);
    auto size = format("{} x {}", val.width(), val.height());
    draw_label_widget(win, var.name, size);
    if (color != zero4f) draw_style_widget_pop(win);
    return false;
}

template <typename T>
bool draw_scene_elem_widgets(gl_window* win, std::vector<T*>& val,
    const visit_var& var, scene* scn,
    const std::unordered_map<std::string, std::string>& highlights,
    const std::string& elem_name) {
    auto color = get_highlight_color(highlights, elem_name + "___" + var.name);
    if (color != zero4f) draw_style_widget_push(win, color);
    auto edited = 0;
    for (auto idx = 0; idx < val.size(); idx++) {
        auto var_ = var;
        var_.name = var.name + "[" + std::to_string(idx) + "]";
        edited += draw_scene_elem_widgets(
            win, val[idx], var_, scn, highlights, elem_name);
    }
    if (color != zero4f) draw_style_widget_pop(win);
    return edited;
}

template <typename T>
bool draw_scene_elem_widgets(gl_window* win, std::vector<T>& val,
    const visit_var& var, scene* scn,
    const std::unordered_map<std::string, std::string>& highlights,
    const std::string& elem_name) {
    if (val.empty()) return false;
    auto color = get_highlight_color(highlights, elem_name + "___" + var.name);
    if (color != zero4f) draw_style_widget_push(win, color);
    draw_label_widget(win, var.name, (int)val.size());
    if (color != zero4f) draw_style_widget_pop(win);
    return false;
}

template <typename T>
bool draw_scene_elem_widgets(gl_window* win, T*& val, const visit_var& var,
    scene* scn, const std::unordered_map<std::string, std::string>& highlights,
    const std::string& elem_name) {
    auto edited = false;
    auto color = get_highlight_color(highlights, elem_name + "___" + var.name);
    if (color != zero4f) draw_style_widget_push(win, color);
    edited = draw_combo_widget(win, var.name, val, get_scene_elems(scn, val));
    if (color != zero4f) draw_style_widget_pop(win);
    return edited;
}

template <typename T, typename TT>
bool draw_scene_elem_widgets(gl_window* win, T* elem, std::vector<TT*>& telems,
    scene* scn,
    const std::unordered_map<std::string, std::string>& highlights) {
    if (!elem) return false;
    auto telem = (TT*)nullptr;
    for (auto te : telems)
        if (te->name == elem->name) telem = te;
    auto edited = false;
    if (telem) {
        visit(*telem, [win, scn, &edited](auto& val, const auto& var) {
            auto edited_ = draw_scene_elem_widgets(win, val, var, scn, {}, "");
            edited = edited || edited_;
        });
        if (edited) update_proc_elem(scn, elem, telem);
    }
    visit(*elem,
        [win, scn, elem, &edited, &highlights](auto& val, const auto& var) {
            auto edited_ = draw_scene_elem_widgets(
                win, val, var, scn, highlights, elem->name);
            edited = edited || edited_;
        });
    return edited;
}

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
        sel = elems.back();
        update_list.push_back(sel);
        update_proc_elem(scn, elems.back(), proc_elems.back());
        return true;
    }
    draw_continue_widget(win);
    return false;
}

bool draw_scene_tree_widgets(gl_window* win, const std::string& lbl, scene* scn,
    scene_selection& sel, std::vector<ygl::scene_selection>& update_list,
    proc_scene* test_scn,
    const std::unordered_map<std::string, std::string>& inspector_highlights) {
    static auto test_scn_def = proc_scene();

    if (!scn) return false;
    if (!lbl.empty() && !draw_header_widget(win, lbl)) return false;
    draw_groupid_widget_push(win, scn);
    // draw_scroll_widget_begin(win, "scene #$%^!@", 240, false);
    visit(scn, [win, &sel, &inspector_highlights](auto& val, const auto& var) {
        draw_scene_tree_widgets(win, var.name, val, sel, inspector_highlights);
    });
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
        if (shp_add && !scn->nodes.empty()) {
            scn->nodes.push_back(new node());
            scn->nodes.back()->name = scn->shapes.back()->name;
            scn->nodes.back()->shp = scn->shapes.back();
            test_scn->nodes.push_back(new proc_node());
            test_scn->nodes.back()->name = scn->nodes.back()->name;
            test_scn->nodes.back()->shape = scn->nodes.back()->shp->name;
            update_list.push_back({});
            update_list.back() = scn->nodes.back();
        }
        draw_add_elem_widgets(
            win, scn, "nde", scn->nodes, test_scn->nodes, sel, update_list);
        draw_add_elem_widgets(win, scn, "env", scn->environments,
            test_scn->environments, sel, update_list);
        draw_add_elem_widgets(win, scn, "anim", scn->animations,
            test_scn->animations, sel, update_list);
    }

    draw_groupid_widget_pop(win);
    return update_list.size() != update_len;
}

bool draw_scene_elem_widgets(gl_window* win, const std::string& lbl, scene* scn,
    scene_selection& sel, std::vector<ygl::scene_selection>& update_list,
    proc_scene* test_scn,
    const std::unordered_map<std::string, std::string>& inspector_highlights) {
    static auto test_scn_def = proc_scene();

    if (!scn || !sel.get_raw()) return false;
    if (!lbl.empty() && !draw_header_widget(win, lbl)) return false;
    draw_groupid_widget_push(win, sel.get_raw());

    auto update_len = update_list.size();

    auto test_scn_res = (test_scn) ? test_scn : &test_scn_def;
    auto edited = false;
    if (sel.is<camera>())
        edited = draw_scene_elem_widgets(win, sel.get<camera>(),
            test_scn_res->cameras, scn, inspector_highlights);
    if (sel.is<shape>())
        edited = draw_scene_elem_widgets(win, sel.get<shape>(),
            test_scn_res->shapes, scn, inspector_highlights);
    if (sel.is<texture>())
        edited = draw_scene_elem_widgets(win, sel.get<texture>(),
            test_scn_res->textures, scn, inspector_highlights);
    if (sel.is<material>())
        edited = draw_scene_elem_widgets(win, sel.get<material>(),
            test_scn_res->materials, scn, inspector_highlights);
    if (sel.is<environment>())
        edited = draw_scene_elem_widgets(win, sel.get<environment>(),
            test_scn_res->environments, scn, inspector_highlights);
    if (sel.is<node>())
        edited = draw_scene_elem_widgets(win, sel.get<node>(),
            test_scn_res->nodes, scn, inspector_highlights);
    if (sel.is<animation>())
        edited = draw_scene_elem_widgets(win, sel.get<animation>(),
            test_scn_res->animations, scn, inspector_highlights);
    if (sel.is<animation_group>())
        edited = draw_scene_elem_widgets(win, sel.get<animation_group>(),
            test_scn_res->animations, scn, inspector_highlights);
    if (edited) update_list.push_back(sel);

    draw_groupid_widget_pop(win);
    return update_list.size() != update_len;
}

}  // namespace ygl

#endif
